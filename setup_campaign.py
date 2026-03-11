#!/usr/bin/env python3
"""
Setup a new opto-gatable binder campaign.

Reads a campaign config JSON, assembles Architecture C fusion sequences
for each binder x linker x state combination, generates Boltz-2 YAMLs,
and creates a manifest with domain boundaries.

Usage:
    python3 setup_campaign.py campaign_config.json
    python3 setup_campaign.py defaults/default_config.json
"""

import json
import sys
from pathlib import Path

from linker_library import (
    LOV2_DARK,
    LOV2_LIT,
    UBQ_SEQ,
    LINKER_LIBRARY,
)

PIPELINE_DIR = Path(__file__).parent


def validate_config(cfg):
    """Validate required fields in campaign config."""
    errors = []

    if "campaign_name" not in cfg:
        errors.append("Missing 'campaign_name'")

    target = cfg.get("target", {})
    if not target.get("sequence"):
        errors.append("Missing 'target.sequence'")
    if not target.get("hotspot_residues_1idx"):
        errors.append("Missing 'target.hotspot_residues_1idx'")

    binders = cfg.get("binders", [])
    if not binders:
        errors.append("No binders specified")
    for i, b in enumerate(binders):
        if not b.get("sequence"):
            errors.append(f"Binder {i}: missing 'sequence'")
        if not b.get("name"):
            errors.append(f"Binder {i}: missing 'name'")

    scaffold = cfg.get("scaffold", {})
    if not scaffold.get("lov2_dark_seq"):
        errors.append("Missing 'scaffold.lov2_dark_seq'")
    if not scaffold.get("lov2_lit_seq"):
        errors.append("Missing 'scaffold.lov2_lit_seq'")

    if errors:
        raise ValueError("Config validation failed:\n  " + "\n  ".join(errors))


def get_linker_subset(cfg):
    """Return the linker library subset specified in config."""
    linker_spec = cfg.get("scaffold", {}).get("linkers", "all")

    if linker_spec == "all" or linker_spec is None:
        return dict(LINKER_LIBRARY)

    subset = {}
    for name in linker_spec:
        if name not in LINKER_LIBRARY:
            print(f"  WARNING: linker '{name}' not in library, skipping")
            continue
        subset[name] = LINKER_LIBRARY[name]

    if not subset:
        raise ValueError("No valid linkers specified")
    return subset


def build_archC_fusion(binder_seq, lov2_seq, linker_seq, cap_seq, spacer=""):
    """
    Assemble Architecture C fusion: [Binder]-[spacer]-[LOV2]-[linker]-[Cap]

    Returns (fusion_seq, domain_map) where domain_map has 0-indexed
    [start, end) ranges for each domain.
    """
    parts = [binder_seq]
    if spacer:
        parts.append(spacer)
    parts.extend([lov2_seq, linker_seq, cap_seq])
    fusion = "".join(parts)

    # Compute domain boundaries (0-indexed, half-open)
    pos = 0
    domains = {}

    domains["binder"] = [pos, pos + len(binder_seq)]
    pos += len(binder_seq)

    if spacer:
        domains["spacer"] = [pos, pos + len(spacer)]
        pos += len(spacer)
    else:
        domains["spacer"] = [pos, pos]

    domains["lov2"] = [pos, pos + len(lov2_seq)]
    pos += len(lov2_seq)

    domains["linker"] = [pos, pos + len(linker_seq)]
    pos += len(linker_seq)

    domains["cap"] = [pos, pos + len(cap_seq)]
    pos += len(cap_seq)

    assert pos == len(fusion), f"Domain sum {pos} != fusion length {len(fusion)}"
    return fusion, domains


def make_yaml_1chain(seq, name):
    """Generate Boltz-2 YAML for single-chain structure prediction."""
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        f"      id: A",
        f"      sequence: {seq}",
        "      msa: empty",
    ]
    return "\n".join(lines) + "\n"


def make_yaml_2chain(fusion_seq, target_seq, contacts=None):
    """
    Generate Boltz-2 YAML for 2-chain prediction with optional contact constraints.

    contacts: list of dicts with 'fusion_res_1idx' and 'target_res_1idx' keys
    """
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {fusion_seq}",
        "      msa: empty",
        "  - protein:",
        "      id: B",
        f"      sequence: {target_seq}",
        "      msa: empty",
    ]

    if contacts:
        lines.append("constraints:")
        for c in contacts:
            lines.extend([
                "  - contact:",
                f"      token1: [A, {c['fusion_res_1idx']}]",
                f"      token2: [B, {c['target_res_1idx']}]",
                f"      max_distance: {c.get('max_distance', 10.0)}",
            ])

    return "\n".join(lines) + "\n"


def setup_campaign(config_path):
    """Main entry: read config, assemble fusions, generate YAMLs + manifest."""
    config_path = Path(config_path)
    with open(config_path) as f:
        cfg = json.load(f)

    validate_config(cfg)

    name = cfg["campaign_name"]
    campaign_dir = PIPELINE_DIR / "campaigns" / name
    yaml_dir = campaign_dir / "boltz2_yamls"
    results_dir = campaign_dir / "boltz2_results"
    analysis_dir = campaign_dir / "analysis"
    sweep_dir = campaign_dir / "sweep_volumes"
    chimerax_dir = campaign_dir / "chimerax"

    for d in [yaml_dir, results_dir, analysis_dir, sweep_dir, chimerax_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Copy config into campaign directory
    (campaign_dir / "campaign_config.json").write_text(
        json.dumps(cfg, indent=2) + "\n"
    )

    # Resolve sequences
    scaffold = cfg["scaffold"]
    lov2_dark = scaffold.get("lov2_dark_seq", LOV2_DARK)
    lov2_lit = scaffold.get("lov2_lit_seq", LOV2_LIT)
    cap_seq = scaffold.get("cap_seq", UBQ_SEQ)
    cap_name = scaffold.get("cap_name", "UBQ")
    spacer = scaffold.get("binder_lov2_spacer", "")
    linkers = get_linker_subset(cfg)

    pred_cfg = cfg.get("prediction", {})
    n_samples = pred_cfg.get("diffusion_samples", 5)
    n_recycle = pred_cfg.get("recycling_steps", 3)

    manifest = []
    yaml_names = []

    print(f"Campaign: {name}")
    print(f"  Binders: {len(cfg['binders'])}")
    print(f"  Linkers: {len(linkers)}")
    print(f"  Cap: {cap_name} ({len(cap_seq)}aa)")
    print(f"  Spacer: '{spacer}' ({len(spacer)}aa)")
    print(f"  States: dark ({len(lov2_dark)}aa LOV2) + lit ({len(lov2_lit)}aa LOV2)")
    print()

    for binder in cfg["binders"]:
        binder_name = binder["name"]
        binder_seq = binder["sequence"]

        for linker_name, linker_info in linkers.items():
            linker_seq = linker_info["seq"]

            for state, lov2_seq in [("dark", lov2_dark), ("lit", lov2_lit)]:
                construct_name = f"{binder_name}_{linker_name}_{state}"

                fusion, domains = build_archC_fusion(
                    binder_seq, lov2_seq, linker_seq, cap_seq, spacer
                )

                # Generate single-chain YAML
                yaml_content = make_yaml_1chain(fusion, construct_name)
                yaml_path = yaml_dir / f"{construct_name}.yaml"
                yaml_path.write_text(yaml_content)
                yaml_names.append(construct_name)

                # Manifest entry
                entry = {
                    "construct_name": construct_name,
                    "binder_name": binder_name,
                    "binder_seq": binder_seq,
                    "binder_iptm": binder.get("iptm"),
                    "reference_complex_cif": binder.get("reference_complex_cif"),
                    "linker_type": linker_name,
                    "linker_seq": linker_seq,
                    "linker_length": linker_info["length"],
                    "linker_geometry": linker_info["type"],
                    "state": state,
                    "lov2_length": len(lov2_seq),
                    "cap_name": cap_name,
                    "cap_length": len(cap_seq),
                    "spacer": spacer,
                    "spacer_length": len(spacer),
                    "fusion_length": len(fusion),
                    "fusion_seq": fusion,
                    "domains": domains,
                    "yaml_path": str(yaml_path),
                }

                # Add linker-specific geometry
                if linker_info["type"] == "kinked":
                    entry["pre_helix"] = linker_info.get("pre_helix")
                    entry["kink_res"] = linker_info.get("kink_res")
                    entry["post_helix"] = linker_info.get("post_helix")
                    entry["multi_kink"] = linker_info.get("multi_kink", False)

                manifest.append(entry)

        print(f"  {binder_name}: {len(linkers) * 2} constructs")

    # Write manifest
    manifest_path = campaign_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    # Generate run_predictions.sh
    run_script = generate_run_script(
        yaml_dir, results_dir, yaml_names, n_samples, n_recycle
    )
    run_path = campaign_dir / "run_predictions.sh"
    run_path.write_text(run_script)
    run_path.chmod(0o755)

    total = len(manifest)
    print(f"\nGenerated:")
    print(f"  {total} YAMLs in {yaml_dir}/")
    print(f"  Manifest: {manifest_path}")
    print(f"  Runner:   {run_path}")
    print(f"\nNext: bash {run_path}")

    return manifest


def generate_run_script(yaml_dir, results_dir, yaml_names, n_samples, n_recycle):
    """Generate bash script to run all Boltz-2 predictions sequentially."""
    lines = [
        "#!/bin/bash",
        "# Auto-generated Boltz-2 prediction runner",
        "# Run predictions sequentially (GPU memory constraint)",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate boltz_only",
        "",
        f'YAML_DIR="{yaml_dir}"',
        f'RESULTS_DIR="{results_dir}"',
        "",
        f"TOTAL={len(yaml_names)}",
        "COUNT=0",
        "",
    ]

    for name in yaml_names:
        lines.extend([
            'COUNT=$((COUNT + 1))',
            f'echo "[${{COUNT}}/${{TOTAL}}] {name}..."',
            '',
            '# Skip if already has CIF files',
            f'if find "$RESULTS_DIR/{name}" -name "*.cif" 2>/dev/null | grep -q cif; then',
            '    echo "  (already complete, skipping)"',
            'else',
            '    python3 -c "import torch; torch.cuda.empty_cache()" 2>/dev/null',
            f'    boltz predict "$YAML_DIR/{name}.yaml" \\',
            f'        --out_dir "$RESULTS_DIR/{name}" \\',
            '        --model boltz2 \\',
            f'        --diffusion_samples {n_samples} \\',
            f'        --recycling_steps {n_recycle} \\',
            '        --override 2>&1 | tail -3',
            f'    CIF_COUNT=$(find "$RESULTS_DIR/{name}" -name "*.cif" 2>/dev/null | wc -l)',
            '    echo "  $CIF_COUNT models generated"',
            'fi',
            '',
        ])

    lines.extend([
        'echo ""',
        'echo "All predictions complete."',
        f'echo "Results: $RESULTS_DIR/"',
    ])

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 setup_campaign.py <campaign_config.json>")
        print("       python3 setup_campaign.py defaults/default_config.json")
        sys.exit(1)

    setup_campaign(sys.argv[1])
