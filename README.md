# Opto-Gatable Binder-PROTAC Pipeline

Design optogenetically controllable protein binders for any target. Given a target protein and its functional hotspot residues, the pipeline generates de novo binders, integrates them into a light-switchable scaffold, and ranks constructs by predicted gating efficacy.

```
Target PDB + hotspots
        |
        v
  [RFdiffusion] ──> backbone designs
        |
        v
  [LigandMPNN]  ──> sequence designs
        |
        v
  [Boltz-2]     ──> structure prediction + scoring
        |
        v
  Top N binders
        |
        v
  [Opto-scaffold assembly]
  [Binder]-[GGSGGS]-[LOV2]-[linker]-[Cap]
        |
        v
  [Boltz-2 sweep predictions]
        |
        v
  [Sweep volume analysis] ──> ranked opto-gatable constructs
        |
        v
  (Optional) [VVD + E3 ligase] ──> full degrader constructs
```

## Quick Start

```bash
# 1. Clone and install
git clone <repo-url>
cd opto_scaffold_pipeline
conda activate boltz_only  # or any env with Python 3.10+
pip install -r requirements.txt

# 2. Run tests
python -m unittest discover -s tests -v

# 3. Try the demo (opto-scaffold analysis only, no GPU needed for setup)
python3 setup_campaign.py defaults/default_config.json
# This generates Boltz-2 YAMLs + prediction scripts in campaigns/p300_hat_demo/
```

## Full Pipeline

### Prerequisites

| Tool | Conda Env | Python | Purpose |
|------|-----------|--------|---------|
| Boltz-2 v2.2.1 | `boltz_only` | 3.10 | Structure prediction + scoring |
| RFdiffusion v1.1.0 | `rfdiffusion_clean` | 3.9 | De novo backbone design |
| LigandMPNN | `ligandmpnn` | 3.10 | Sequence design |

**Hardware:** 1x 24GB GPU minimum, 8x 24GB GPUs recommended for throughput.

**Environment setup:**
```bash
bash setup_environment.sh all
```

### Step-by-Step Workflow

#### 1. Create Campaign Config

Copy and modify the template:

```bash
cp defaults/default_config.json my_target_config.json
```

Edit `my_target_config.json`:
- `campaign_name`: unique name for this campaign
- `target.sequence`: full amino acid sequence of your target protein
- `target.hotspot_residues_1idx`: functional residues you want the binder to contact (1-indexed)
- `binder_design.target_pdb`: path to target PDB file for RFdiffusion
- `binder_design.contig`: RFdiffusion contig string (e.g., `"A1-300/0"` for chain A residues 1-300)
- `binder_design.hotspot_residues_pdb`: PDB-format hotspot residues (e.g., `["A150", "A160"]`)
- `binder_design.rfdiffusion_path`: path to your RFdiffusion installation
- `binder_design.ligandmpnn_path`: path to your LigandMPNN installation
- `binder_design.n_gpus`: number of GPUs available

#### 2. Design Binders

```bash
conda activate boltz_only
python3 screen_binders.py my_target_config.json --n_gpus 8
```

This generates bash scripts in `campaigns/<name>/binder_screen/`:

| Script | Conda Env | What It Does |
|--------|-----------|-------------|
| `step1_rfdiffusion.sh` | rfdiffusion_clean | Generate backbone designs |
| `step2_ligandmpnn.sh` | ligandmpnn | Design sequences for backbones |
| `step3_boltz2_yamls.sh` | boltz_only | Generate Boltz-2 input files |
| `step4_boltz2_predict.sh` | boltz_only | Run Boltz-2 (multi-GPU) |
| `step5_score.sh` | boltz_only | Score + rank binders |
| `run_all.sh` | (switches) | Run all steps sequentially |

Run all steps:
```bash
bash campaigns/<name>/binder_screen/run_all.sh
```

Or run individual steps:
```bash
bash campaigns/<name>/binder_screen/step1_rfdiffusion.sh
bash campaigns/<name>/binder_screen/step2_ligandmpnn.sh
# ... etc
```

#### 3. Score Binders

```bash
python3 score_binders.py campaigns/<name>/binder_screen/ \
    --top_n 10 --config my_target_config.json
```

Outputs:
- `binder_screen/ranked_binders.json` — all binders ranked by ipTM
- `binder_screen/top_binders.fasta` — top N binder sequences
- `campaign_config_opto.json` — auto-generated config for opto-scaffold setup

#### 4. Setup Opto-Scaffold Campaign

```bash
python3 setup_campaign.py campaigns/<name>/campaign_config_opto.json --n_gpus 8
```

Assembles Architecture C fusions for each binder × linker × state:
```
[Binder]-[GGSGGS]-[LOV2_dark or LOV2_lit]-[linker]-[UBQ]
```

12 linkers × 2 states × N binders = constructs to predict.

#### 5. Run Boltz-2 Sweep Predictions

```bash
bash campaigns/<name>_opto/run_predictions.sh
```

With 8 GPUs, predictions run in parallel across devices.

#### 6. Analyze Results

```bash
python3 analyze_campaign.py campaigns/<name>_opto/
```

Outputs in `campaigns/<name>_opto/analysis/`:
- `sweep_results.json` — full sweep volume analysis per construct
- `ranked_constructs.json` — constructs ranked by gating ratio
- `power_analysis.json` — predicted success rates + how many to test

#### 7. (Optional) Degron Assembly

For full degrader constructs with VVD dimerizer + E3 ligase:

```bash
python3 assemble_degron.py campaigns/<name>/campaign_config_opto.json --pathway proteasomal
```

### Single-Command Alternative

The `run_pipeline.py` orchestrator manages all phases:

```bash
# Phase 1: Generate binder design scripts
python3 run_pipeline.py my_target_config.json

# Phase 2: After screening completes, score + setup opto
python3 run_pipeline.py my_target_config.json --continue

# Phase 3: After opto predictions complete, analyze
python3 run_pipeline.py my_target_config.json --analyze
```

## Multi-GPU Support

Set `n_gpus` in your config or use `--n_gpus` flag. The pipeline distributes Boltz-2 predictions across GPUs using `CUDA_VISIBLE_DEVICES`:

```bash
python3 screen_binders.py config.json --n_gpus 8
python3 setup_campaign.py config.json --n_gpus 8
```

**Estimated runtime per target (100-variant screen, 8 GPUs):**
- RFdiffusion: ~2.5 hrs (100 backbones)
- LigandMPNN: ~0.5 hrs (100 backbones × 5 seqs)
- Boltz-2 binder screen: ~1 hr (500 predictions ÷ 8 GPUs)
- Boltz-2 opto-scaffold: ~0.5 hr (240 predictions ÷ 8 GPUs)
- Total: ~4.5 hrs

## Architecture

The opto-gatable scaffold uses Architecture C, where the linker is placed between LOV2 and the cap domain:

```
Dark state (binder occluded):
[Binder]─[GGSGGS]─[LOV2═══Jα]─[linker]─[Cap]
                         ↕                  ↕
                    Jα helix           Cap blocks
                    folded             binder face

Light state (binder exposed):
[Binder]─[GGSGGS]─[LOV2   Jα]─[linker]─[Cap]
                         ↕                  ↕
                    Jα helix           Cap swings
                    unfolds            away → binder
                                       exposed
```

The Jα helix unfolding at the LOV2 C-terminus acts as a pivot point. A rigid/kinked linker after the pivot amplifies the angular displacement into large cap movement.

**12 linkers** are tested per binder, including:
- Kinked linkers with proline breaks (E2pppE2, E2dpgnE2, etc.)
- Straight alpha-helical linkers (EAAAK×3, EAAAK×5)
- Polyproline II helices (P16)

**Sweep volume analysis** models the cap's reachable volume as a cone from the LOV2 C-terminal pivot:
- Dark: ±25° cone (Jα constrains)
- Lit: ±60° cone (Jα released)
- Gating ratio = dark hotspot blocking / lit hotspot blocking

## Config Reference

```jsonc
{
  "campaign_name": "my_campaign",        // Unique name, becomes folder
  "target": {
    "name": "display name",
    "sequence": "PROTEIN...",            // Full target sequence
    "hotspot_residues_1idx": [91, 110],  // Functional residues (1-indexed)
    "pdb_path": null                     // Optional: PDB for visualization
  },
  "binders": [                           // Pre-validated binders (or auto-filled by score_binders.py)
    {"name": "binder1", "sequence": "AAA...", "iptm": 0.85}
  ],
  "scaffold": {
    "architecture": "C",                 // Only "C" supported currently
    "lov2_dark_seq": "...",              // LOV2 dark state (136aa)
    "lov2_lit_seq": "...",               // LOV2 lit state (121aa, Jα truncated)
    "cap_seq": "...",                    // Steric cap (default: ubiquitin, 76aa)
    "cap_name": "UBQ",
    "binder_lov2_spacer": "GGSGGS",     // Flexible spacer between binder and LOV2
    "linkers": "all",                    // "all" or list of specific linker names
    "dark_cone_half_angle": 25.0,        // Sweep analysis parameter (degrees)
    "lit_cone_half_angle": 60.0
  },
  "prediction": {
    "diffusion_samples": 5,              // Boltz-2 samples per prediction
    "recycling_steps": 3,
    "model": "boltz2"
  },
  "binder_design": {                     // For screen_binders.py
    "lengths": [20, 30, 40, 50],         // Binder lengths to design
    "designs_per_length": 25,            // RFdiffusion backbones per length
    "seqs_per_backbone": 5,              // LigandMPNN sequences per backbone
    "boltz2_samples": 5,
    "top_n_binders": 10,
    "n_gpus": 8,
    "rfdiffusion_path": "/path/to/RFdiffusion",
    "ligandmpnn_path": "/path/to/LigandMPNN",
    "target_pdb": "/path/to/target.pdb",
    "hotspot_residues_pdb": ["A150", "A160"],  // PDB-format for RFdiffusion
    "contig": "A1-300/0"                       // RFdiffusion contig map
  }
}
```

## Output Files

| File | Description |
|------|-------------|
| `binder_screen/ranked_binders.json` | All binders ranked by ipTM |
| `binder_screen/top_binders.fasta` | Top N binder sequences |
| `campaign_config_opto.json` | Auto-generated opto-scaffold config |
| `manifest.json` | All construct metadata + domain boundaries |
| `analysis/sweep_results.json` | Sweep volume analysis per construct |
| `analysis/ranked_constructs.json` | Constructs ranked by gating ratio |
| `analysis/power_analysis.json` | Statistical power analysis |
| `chimerax/*.cxc` | ChimeraX visualization scripts |

## Tests

```bash
conda activate boltz_only
cd opto_scaffold_pipeline
python -m unittest discover -s tests -v
```

6 test modules covering all pipeline components:
- `test_linker_library.py` — sequence constants, linker library, physical parameters
- `test_setup_campaign.py` — config validation, fusion assembly, YAML generation
- `test_analyze_campaign.py` — sweep geometry, ranking, power analysis
- `test_analyze_binding_validation.py` — confidence extraction, PAE analysis
- `test_degron_library.py` — degron component sequences
- `test_assemble_degron.py` — degron fusion assembly

## Troubleshooting

**Boltz-2 OOM:** Reduce `diffusion_samples` to 3, or ensure only 1 prediction per GPU.

**RFdiffusion timeout:** Some designs take >2 min. Increase timeout in step1 script or reduce binder length.

**LigandMPNN missing weights:** Run `bash get_model_params.sh` in the LigandMPNN directory.

**Empty FASTA after LigandMPNN:** Check that RFdiffusion produced PDB files with 2 chains (target=A, binder=B).

**Low ipTM scores (<0.3):** Normal for fusion constructs. Binder validation uses isolated binder+target predictions (ipTM >0.5 is promising, >0.7 is strong).

**`import error: No module named linker_library`:** Run from the pipeline root directory, or add it to PYTHONPATH.
