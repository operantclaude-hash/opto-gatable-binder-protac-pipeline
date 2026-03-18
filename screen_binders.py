#!/usr/bin/env python3
"""
Generate binder design + screening scripts for a new target.

Reads campaign config and generates step-by-step bash scripts that run
RFdiffusion → LigandMPNN → Boltz-2 YAML generation → Boltz-2 prediction → scoring.
Each script activates the correct conda environment.

Usage:
    python3 screen_binders.py campaign_config.json
    python3 screen_binders.py campaign_config.json --n_gpus 8

Generates:
    campaigns/<name>/binder_screen/
        step1_rfdiffusion.sh     (conda: rfdiffusion_clean)
        step2_ligandmpnn.sh      (conda: ligandmpnn)
        step3_boltz2_yamls.sh    (conda: boltz_only)
        step4_boltz2_predict.sh  (conda: boltz_only, multi-GPU)
        step5_score.sh           (conda: boltz_only)
        run_all.sh               (sequential runner for all steps)
"""

import argparse
import json
import sys
from pathlib import Path

from linker_library import LOV2_DARK, LOV2_LIT, UBQ_SEQ

PIPELINE_DIR = Path(__file__).parent


def validate_binder_design_config(cfg):
    """Validate binder_design section of campaign config."""
    errors = []

    target = cfg.get("target", {})
    if not target.get("sequence"):
        errors.append("Missing 'target.sequence'")
    if not target.get("hotspot_residues_1idx"):
        errors.append("Missing 'target.hotspot_residues_1idx'")

    bd = cfg.get("binder_design", {})
    if not bd.get("target_pdb"):
        errors.append("Missing 'binder_design.target_pdb' — PDB file of target for RFdiffusion")
    if not bd.get("rfdiffusion_path"):
        errors.append("Missing 'binder_design.rfdiffusion_path' — path to RFdiffusion installation")
    if not bd.get("ligandmpnn_path"):
        errors.append("Missing 'binder_design.ligandmpnn_path' — path to LigandMPNN installation")
    if not bd.get("hotspot_residues_pdb"):
        errors.append("Missing 'binder_design.hotspot_residues_pdb' — e.g. ['A1436','A1456']")
    if not bd.get("contig"):
        errors.append("Missing 'binder_design.contig' — RFdiffusion contig string, e.g. 'A1287-1519/0'")

    if errors:
        raise ValueError(
            "binder_design config validation failed:\n  " + "\n  ".join(errors)
        )


def generate_rfdiffusion_script(cfg, screen_dir):
    """Generate step 1: RFdiffusion backbone design script."""
    bd = cfg["binder_design"]
    rfd_path = bd["rfdiffusion_path"]
    target_pdb = bd["target_pdb"]
    contig = bd["contig"]
    hotspots = bd["hotspot_residues_pdb"]
    lengths = bd.get("lengths", [20, 30, 40, 50])
    designs_per_length = bd.get("designs_per_length", 25)
    max_retries = bd.get("max_retries", 3)

    hotspot_str = ",".join(hotspots)
    outdir = screen_dir / "rfdiffusion_outputs"

    lines = [
        "#!/bin/bash",
        "# Step 1: RFdiffusion backbone generation",
        f"# Designs {len(lengths)} lengths × {designs_per_length} = {len(lengths) * designs_per_length} backbones",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate rfdiffusion_clean",
        "",
        f'RFDIFFUSION_DIR="{rfd_path}"',
        f'TARGET_PDB="{target_pdb}"',
        f'OUTDIR="{outdir}"',
        f"MAX_RETRIES={max_retries}",
        "",
        'mkdir -p "$OUTDIR"',
        'cp "$TARGET_PDB" "$OUTDIR/target.pdb"',
        "",
        'cd "$RFDIFFUSION_DIR"',
        "",
        "TOTAL_SUCCESS=0",
        "TOTAL_FAIL=0",
        "",
    ]

    for length in lengths:
        lines += [
            f'echo "=== Designing {length}aa binders ({designs_per_length} designs) ==="',
            f"for i in $(seq 0 {designs_per_length - 1}); do",
            f'    PREFIX="binder_{length}aa"',
            f'    OUTPUT_PDB="${{OUTDIR}}/${{PREFIX}}_${{i}}.pdb"',
            "",
            '    if [ -f "$OUTPUT_PDB" ]; then',
            "        ((TOTAL_SUCCESS++))",
            "        continue",
            "    fi",
            "",
            "    success=false",
            "    for attempt in $(seq 1 $MAX_RETRIES); do",
            "        timeout 120s ./scripts/run_inference.py \\",
            '            inference.output_prefix="${OUTDIR}/${PREFIX}" \\',
            '            inference.input_pdb="${OUTDIR}/target.pdb" \\',
            f"            'contigmap.contigs=[\"{contig} {length}-{length}\"]' \\",
            f"            'ppi.hotspot_res=[{hotspot_str}]' \\",
            f"            'potentials.guiding_potentials=[\"type:interface_ncontacts,binderlen:{length},weight:2,r_0:6,d_0:4\"]' \\",
            "            'potentials.guide_scale=2' \\",
            "            inference.num_designs=1 \\",
            "            inference.design_startnum=${i} \\",
            "            denoiser.noise_scale_ca=0 \\",
            "            denoiser.noise_scale_frame=0 \\",
            '            > "${OUTDIR}/rfdiff_${PREFIX}_${i}_attempt${attempt}.log" 2>&1',
            "",
            '        if [ $? -eq 0 ] && [ -f "$OUTPUT_PDB" ]; then',
            "            ((TOTAL_SUCCESS++))",
            "            success=true",
            '            rm -f "${OUTDIR}/rfdiff_${PREFIX}_${i}_attempt"*.log',
            "            break",
            "        fi",
            "    done",
            "",
            '    if [ "$success" = false ]; then',
            "        ((TOTAL_FAIL++))",
            "    fi",
            "",
            "    if [ $((i % 10)) -eq 0 ] && [ $i -gt 0 ]; then",
            f'        echo "  {length}aa: ${{TOTAL_SUCCESS}} ok, ${{TOTAL_FAIL}} failed (${{i}}/{designs_per_length})"',
            "    fi",
            "done",
            "",
        ]

    lines += [
        'echo ""',
        f'echo "RFdiffusion complete: ${{TOTAL_SUCCESS}} successes, ${{TOTAL_FAIL}} failures"',
        f'echo "Outputs: $OUTDIR/"',
    ]

    return "\n".join(lines) + "\n"


def generate_ligandmpnn_script(cfg, screen_dir):
    """Generate step 2: LigandMPNN sequence design script."""
    bd = cfg["binder_design"]
    lmpnn_path = bd["ligandmpnn_path"]
    seqs_per_backbone = bd.get("seqs_per_backbone", 5)

    rfd_outdir = screen_dir / "rfdiffusion_outputs"
    lmpnn_outdir = screen_dir / "ligandmpnn_outputs"

    lines = [
        "#!/bin/bash",
        f"# Step 2: LigandMPNN sequence design ({seqs_per_backbone} seqs per backbone)",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate ligandmpnn",
        "",
        f'LIGANDMPNN_DIR="{lmpnn_path}"',
        f'INDIR="{rfd_outdir}"',
        f'OUTDIR="{lmpnn_outdir}"',
        "",
        'mkdir -p "$OUTDIR"',
        'cd "$LIGANDMPNN_DIR"',
        "",
        'count=0',
        'total=$(ls "${INDIR}"/binder_*_*.pdb 2>/dev/null | wc -l)',
        'echo "Processing ${total} backbone PDBs..."',
        "",
        'for PDB_FILE in "${INDIR}"/binder_*_*.pdb; do',
        '    [ -f "$PDB_FILE" ] || continue',
        '    BASENAME=$(basename "$PDB_FILE" .pdb)',
        '    MPNN_OUT="${OUTDIR}/${BASENAME}"',
        "",
        '    # Skip if already done',
        '    if ls "${MPNN_OUT}/seqs/"*.fa 1>/dev/null 2>&1; then',
        "        ((count++))",
        "        continue",
        "    fi",
        "",
        '    mkdir -p "$MPNN_OUT"',
        "",
        "    python run.py \\",
        '        --model_type "ligand_mpnn" \\',
        '        --checkpoint_ligand_mpnn "${LIGANDMPNN_DIR}/model_params/ligandmpnn_v_32_010_25.pt" \\',
        '        --pdb_path "$PDB_FILE" \\',
        '        --out_folder "$MPNN_OUT" \\',
        f"        --number_of_batches {seqs_per_backbone} \\",
        "        --batch_size 1 \\",
        "        --temperature 0.1 \\",
        '        --chains_to_design "B" \\',
        "        2>&1 | tail -1",
        "",
        "    ((count++))",
        "    if [ $((count % 25)) -eq 0 ]; then",
        '        echo "  LigandMPNN: ${count}/${total}"',
        "    fi",
        "done",
        "",
        "# Consolidate all sequences into a single FASTA",
        'COMBINED_FASTA="${OUTDIR}/all_binder_sequences.fa"',
        '> "$COMBINED_FASTA"',
        "",
        'for MPNN_DIR in "${OUTDIR}"/binder_*; do',
        '    [ -d "$MPNN_DIR" ] || continue',
        '    BASENAME=$(basename "$MPNN_DIR")',
        "",
        '    for FA in "${MPNN_DIR}"/seqs/*.fa; do',
        '        [ -f "$FA" ] || continue',
        "        python3 -c \"",
        "with open('$FA') as f:",
        "    lines = f.readlines()",
        "for i in range(0, len(lines), 2):",
        "    if i+1 >= len(lines): break",
        "    header = lines[i].strip()",
        "    seq = lines[i+1].strip()",
        "    if ':' in seq:",
        "        seq = seq.split(':')[-1]",
        "    # Skip seq0 (polyglycine wild-type backbone from RFdiffusion)",
        "    if i == 0 and set(seq) <= {'G', 'X'}:",
        "        continue",
        "    name = '${BASENAME}_seq' + str(i//2)",
        "    print(f'>{name}')",
        "    print(seq)",
        '\" >> "$COMBINED_FASTA"',
        "    done",
        "done",
        "",
        'TOTAL_SEQS=$(grep -c "^>" "$COMBINED_FASTA" 2>/dev/null || echo 0)',
        'echo "Combined FASTA: ${TOTAL_SEQS} sequences -> ${COMBINED_FASTA}"',
    ]

    return "\n".join(lines) + "\n"


def generate_boltz2_yamls_script(cfg, screen_dir):
    """Generate step 3: Boltz-2 YAML generation for binder+target 2-chain predictions."""
    target_seq = cfg["target"]["sequence"]

    lmpnn_outdir = screen_dir / "ligandmpnn_outputs"
    yaml_dir = screen_dir / "boltz2_yamls"

    lines = [
        "#!/bin/bash",
        "# Step 3: Generate Boltz-2 YAML files for 2-chain binder+target prediction",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate boltz_only",
        "",
        f'FASTA="{lmpnn_outdir}/all_binder_sequences.fa"',
        f'YAML_DIR="{yaml_dir}"',
        f'TARGET_SEQ="{target_seq}"',
        "",
        'mkdir -p "$YAML_DIR"',
        "",
        'if [ ! -f "$FASTA" ]; then',
        '    echo "ERROR: $FASTA not found. Run step 2 first."',
        "    exit 1",
        "fi",
        "",
        "count=0",
        'while IFS= read -r line; do',
        '    if [[ "$line" == ">"* ]]; then',
        '        NAME="${line#>}"',
        '        NAME="${NAME%% *}"',
        "    else",
        '        BINDER_SEQ="$line"',
        '        YAML_FILE="${YAML_DIR}/${NAME}.yaml"',
        "",
        '        cat > "$YAML_FILE" << YAMLEOF',
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        '      sequence: ${TARGET_SEQ}',
        "      msa: empty",
        "  - protein:",
        "      id: B",
        '      sequence: ${BINDER_SEQ}',
        "      msa: empty",
        "YAMLEOF",
        "",
        "        ((count++))",
        "    fi",
        'done < "$FASTA"',
        "",
        'echo "Generated ${count} YAML files in ${YAML_DIR}/"',
    ]

    return "\n".join(lines) + "\n"


def generate_boltz2_predict_script(cfg, screen_dir, n_gpus):
    """Generate step 4: Boltz-2 prediction script with multi-GPU support."""
    bd = cfg.get("binder_design", {})
    pred_cfg = cfg.get("prediction", {})
    n_samples = pred_cfg.get("diffusion_samples", bd.get("boltz2_samples", 5))
    n_recycle = pred_cfg.get("recycling_steps", 3)

    yaml_dir = screen_dir / "boltz2_yamls"
    results_dir = screen_dir / "boltz2_results"

    lines = [
        "#!/bin/bash",
        f"# Step 4: Boltz-2 predictions ({n_gpus} GPU{'s' if n_gpus > 1 else ''})",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate boltz_only",
        "",
        f'YAML_DIR="{yaml_dir}"',
        f'RESULTS_DIR="{results_dir}"',
        f"N_GPUS={n_gpus}",
        f"N_SAMPLES={n_samples}",
        f"N_RECYCLE={n_recycle}",
        "",
        'mkdir -p "$RESULTS_DIR"',
        "",
        'YAMLS=($(ls "$YAML_DIR"/*.yaml 2>/dev/null))',
        'TOTAL=${#YAMLS[@]}',
        "",
        'if [ $TOTAL -eq 0 ]; then',
        '    echo "ERROR: No YAML files in $YAML_DIR. Run step 3 first."',
        "    exit 1",
        "fi",
        "",
        'echo "Running Boltz-2 on ${TOTAL} predictions across ${N_GPUS} GPU(s)..."',
        "",
    ]

    if n_gpus > 1:
        lines += [
            "# Multi-GPU: split predictions across GPUs",
            'PER_GPU=$(( (TOTAL + N_GPUS - 1) / N_GPUS ))',
            "",
            "predict_chunk() {",
            "    local GPU_ID=$1",
            "    shift",
            '    local CHUNK=("$@")',
            '    local DONE=0',
            "",
            '    for YAML in "${CHUNK[@]}"; do',
            '        [ -f "$YAML" ] || continue',
            '        NAME=$(basename "$YAML" .yaml)',
            '        RESULT_DIR="${RESULTS_DIR}/${NAME}"',
            "",
            "        # Skip if already complete",
            '        if find "$RESULT_DIR" -name "*.cif" 2>/dev/null | grep -q cif; then',
            "            ((DONE++))",
            "            continue",
            "        fi",
            "",
            '        CUDA_VISIBLE_DEVICES=$GPU_ID boltz predict "$YAML" \\',
            '            --out_dir "$RESULT_DIR" \\',
            "            --model boltz2 \\",
            "            --diffusion_samples $N_SAMPLES \\",
            "            --recycling_steps $N_RECYCLE \\",
            '            --override 2>&1 | tail -2',
            "",
            "        ((DONE++))",
            '        echo "  [GPU $GPU_ID] ${DONE}/${#CHUNK[@]}"',
            "    done",
            "}",
            "",
            "# Launch parallel GPU workers",
            "for GPU in $(seq 0 $((N_GPUS-1))); do",
            "    START=$((GPU * PER_GPU))",
            '    CHUNK=("${YAMLS[@]:$START:$PER_GPU}")',
            '    if [ ${#CHUNK[@]} -gt 0 ]; then',
            '        echo "GPU $GPU: ${#CHUNK[@]} predictions"',
            '        predict_chunk $GPU "${CHUNK[@]}" &',
            "    fi",
            "done",
            "",
            'echo "Waiting for all GPUs to finish..."',
            "wait",
        ]
    else:
        lines += [
            "# Single GPU: sequential predictions",
            "COUNT=0",
            'for YAML in "${YAMLS[@]}"; do',
            '    NAME=$(basename "$YAML" .yaml)',
            '    RESULT_DIR="${RESULTS_DIR}/${NAME}"',
            "",
            "    # Skip if already complete",
            '    if find "$RESULT_DIR" -name "*.cif" 2>/dev/null | grep -q cif; then',
            "        ((COUNT++))",
            "        continue",
            "    fi",
            "",
            '    python3 -c "import torch; torch.cuda.empty_cache()" 2>/dev/null',
            '    boltz predict "$YAML" \\',
            '        --out_dir "$RESULT_DIR" \\',
            "        --model boltz2 \\",
            "        --diffusion_samples $N_SAMPLES \\",
            "        --recycling_steps $N_RECYCLE \\",
            '        --override 2>&1 | tail -2',
            "",
            "    ((COUNT++))",
            "    if [ $((COUNT % 25)) -eq 0 ]; then",
            '        echo "  Progress: ${COUNT}/${TOTAL}"',
            "    fi",
            "done",
        ]

    lines += [
        "",
        'echo ""',
        'echo "All predictions complete."',
        'CIF_COUNT=$(find "$RESULTS_DIR" -name "*.cif" 2>/dev/null | wc -l)',
        'echo "Total CIF files: ${CIF_COUNT}"',
        'echo "Results: $RESULTS_DIR/"',
    ]

    return "\n".join(lines) + "\n"


def generate_score_script(cfg, screen_dir):
    """Generate step 5: scoring script that runs score_binders.py."""
    config_path = screen_dir.parent / "campaign_config.json"

    lines = [
        "#!/bin/bash",
        "# Step 5: Score binder predictions and select top candidates",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate boltz_only",
        "",
        f'PIPELINE_DIR="{PIPELINE_DIR}"',
        f'SCREEN_DIR="{screen_dir}"',
        f'CONFIG="{config_path}"',
        "",
        'cd "$PIPELINE_DIR"',
        "",
        'python3 score_binders.py "$SCREEN_DIR" \\',
        f'    --top_n {cfg.get("binder_design", {}).get("top_n_binders", 10)} \\',
        '    --config "$CONFIG"',
        "",
        'echo ""',
        'echo "Next steps:"',
        'echo "  1. Review ranked_binders.json and top_binders.fasta"',
        f'echo "  2. python3 setup_campaign.py {screen_dir.parent}/campaign_config_opto.json"',
        'echo "  3. bash the generated run_predictions.sh"',
        f'echo "  4. python3 analyze_campaign.py {screen_dir.parent}/<campaign_name>_opto/"',
    ]

    return "\n".join(lines) + "\n"


def generate_run_all_script(screen_dir):
    """Generate master script that runs all steps sequentially."""
    lines = [
        "#!/bin/bash",
        "# Run all binder screening steps sequentially",
        "# Each step activates its own conda environment",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        "",
        'echo "=============================================="',
        'echo "Binder Screening Pipeline — Full Run"',
        'echo "=============================================="',
        'echo ""',
        "",
        "STEP=${1:-all}",
        "",
        "run_step() {",
        '    local STEP_NUM=$1',
        '    local SCRIPT=$2',
        '    echo "=== Step ${STEP_NUM}: $(head -2 "$SCRIPT" | tail -1 | sed \'s/^# //\') ==="',
        '    bash "$SCRIPT"',
        '    EXIT_CODE=$?',
        '    if [ $EXIT_CODE -ne 0 ]; then',
        '        echo "WARNING: Step ${STEP_NUM} exited with code $EXIT_CODE"',
        "    fi",
        '    echo ""',
        "}",
        "",
        'case "$STEP" in',
        '    1) run_step 1 "$SCRIPT_DIR/step1_rfdiffusion.sh" ;;',
        '    2) run_step 2 "$SCRIPT_DIR/step2_ligandmpnn.sh" ;;',
        '    3) run_step 3 "$SCRIPT_DIR/step3_boltz2_yamls.sh" ;;',
        '    4) run_step 4 "$SCRIPT_DIR/step4_boltz2_predict.sh" ;;',
        '    5) run_step 5 "$SCRIPT_DIR/step5_score.sh" ;;',
        "    all)",
        '        run_step 1 "$SCRIPT_DIR/step1_rfdiffusion.sh"',
        '        run_step 2 "$SCRIPT_DIR/step2_ligandmpnn.sh"',
        '        run_step 3 "$SCRIPT_DIR/step3_boltz2_yamls.sh"',
        '        run_step 4 "$SCRIPT_DIR/step4_boltz2_predict.sh"',
        '        run_step 5 "$SCRIPT_DIR/step5_score.sh"',
        "        ;;",
        "    *)",
        '        echo "Usage: $0 [1|2|3|4|5|all]"',
        "        exit 1",
        "        ;;",
        "esac",
        "",
        'echo "=============================================="',
        'echo "Binder screening pipeline complete!"',
        'echo "=============================================="',
    ]

    return "\n".join(lines) + "\n"


def screen_binders(config_path, n_gpus_override=None):
    """Main entry: read config, generate all screening scripts."""
    config_path = Path(config_path).resolve()
    with open(config_path) as f:
        cfg = json.load(f)

    validate_binder_design_config(cfg)

    name = cfg["campaign_name"]
    campaign_dir = PIPELINE_DIR / "campaigns" / name
    screen_dir = campaign_dir / "binder_screen"
    screen_dir.mkdir(parents=True, exist_ok=True)

    # Copy config
    (campaign_dir / "campaign_config.json").write_text(
        json.dumps(cfg, indent=2) + "\n"
    )

    bd = cfg.get("binder_design", {})
    n_gpus = n_gpus_override or bd.get("n_gpus", 1)
    lengths = bd.get("lengths", [20, 30, 40, 50])
    dpg = bd.get("designs_per_length", 25)
    spb = bd.get("seqs_per_backbone", 5)
    total_backbones = len(lengths) * dpg
    total_seqs = total_backbones * spb

    print(f"Campaign: {name}")
    print(f"  Binder lengths: {lengths}")
    print(f"  Designs per length: {dpg}")
    print(f"  Total backbones: {total_backbones}")
    print(f"  Seqs per backbone: {spb}")
    print(f"  Total sequences: ~{total_seqs}")
    print(f"  GPUs: {n_gpus}")
    print()

    # Generate scripts
    scripts = [
        ("step1_rfdiffusion.sh", generate_rfdiffusion_script(cfg, screen_dir)),
        ("step2_ligandmpnn.sh", generate_ligandmpnn_script(cfg, screen_dir)),
        ("step3_boltz2_yamls.sh", generate_boltz2_yamls_script(cfg, screen_dir)),
        ("step4_boltz2_predict.sh", generate_boltz2_predict_script(cfg, screen_dir, n_gpus)),
        ("step5_score.sh", generate_score_script(cfg, screen_dir)),
        ("run_all.sh", generate_run_all_script(screen_dir)),
    ]

    for filename, content in scripts:
        path = screen_dir / filename
        path.write_text(content)
        path.chmod(0o755)
        print(f"  Generated: {path}")

    # Estimate GPU time
    est_rfd_min = total_backbones * 1.5  # ~90s per design
    est_lmpnn_min = total_backbones * 0.3  # ~20s per backbone
    est_boltz_min = total_seqs * 3  # ~3 min per 2-chain prediction
    est_boltz_parallel = est_boltz_min / n_gpus

    print(f"\nEstimated GPU time:")
    print(f"  RFdiffusion:  ~{est_rfd_min:.0f} min ({total_backbones} designs)")
    print(f"  LigandMPNN:   ~{est_lmpnn_min:.0f} min ({total_backbones} backbones)")
    print(f"  Boltz-2:      ~{est_boltz_parallel:.0f} min ({total_seqs} predictions ÷ {n_gpus} GPUs)")
    print(f"  Total:        ~{(est_rfd_min + est_lmpnn_min + est_boltz_parallel) / 60:.1f} hrs")

    print(f"\nNext: bash {screen_dir / 'run_all.sh'}")
    print(f"  Or run steps individually: bash {screen_dir / 'step1_rfdiffusion.sh'}")

    return screen_dir


def main():
    parser = argparse.ArgumentParser(
        description="Generate binder screening scripts for a new target"
    )
    parser.add_argument("config", help="Campaign config JSON with binder_design section")
    parser.add_argument("--n_gpus", type=int, default=None,
                        help="Override number of GPUs (default: from config)")
    args = parser.parse_args()

    screen_binders(args.config, args.n_gpus)


if __name__ == "__main__":
    main()
