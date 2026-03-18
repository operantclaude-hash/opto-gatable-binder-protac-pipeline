# Opto-Gatable Binder Pipeline — Claude Code Guide

## Project Structure
- Pipeline root: this directory (`opto_scaffold_pipeline/`)
- Core modules: `setup_campaign.py`, `analyze_campaign.py`, `screen_binders.py`, `score_binders.py`, `run_pipeline.py`
- Data modules: `linker_library.py`, `degron_library.py`
- Degron assembly: `assemble_degron.py`
- Tests: `tests/` (6 files, ~191 unit tests)
- Config template: `defaults/default_config.json`
- Campaign outputs: `campaigns/<name>/`

## Running Tests
```bash
conda activate boltz_only
python -m unittest discover -s tests -v
```

## Key Conventions
- Architecture C = `[Binder]-[spacer]-[LOV2]-[linker]-[Cap]`
- **Never** call binders "nanobodies" — they are de novo computationally designed protein binders
- Boltz-2 pLDDT is on a **0-1 scale** (not 0-100). Threshold: >0.70
- Never run two Boltz-2 predictions on the same GPU simultaneously
- ipTM 0.2-0.3 on fusion constructs is noise, not signal. Validate binders separately.

## Dependencies
- Python packages: `numpy`, `scipy`, `biopython` (see `requirements.txt`)
- External tools (separate conda envs):
  - `boltz_only`: Boltz-2 v2.2.1 (Python 3.10)
  - `rfdiffusion_clean`: RFdiffusion v1.1.0 (Python 3.9)
  - `ligandmpnn`: LigandMPNN (Python 3.10)

## Workflow
1. Create config from `defaults/default_config.json`
2. `python3 screen_binders.py config.json` — generates binder design scripts
3. Run generated bash scripts (step1-5)
4. `python3 score_binders.py campaigns/<name>/binder_screen/` — rank binders
5. `python3 setup_campaign.py campaign_config_opto.json` — assemble opto-scaffold fusions
6. Run Boltz-2 sweep predictions
7. `python3 analyze_campaign.py campaigns/<name>/` — sweep volume analysis + ranking

Or use the orchestrator: `python3 run_pipeline.py config.json`

## Architecture Diagram
```
Gate 2 (target binding):
[Binder]-[GGSGGS]-[LOV2]-[E2pppE2]-[UBQ]
                    |                 |
                    Ja helix         Cap (occludes binder in dark)
                    unfolds in       Released in light → binder exposed
                    blue light

Gate 1 (E3 ligase activation, optional):
[E3]-[GS]-[VVD]-[GS]-[Gate 2 construct]
             |
             VVD dimerizes in blue light → E3 activates
```

## Config Reference
See `defaults/default_config.json` for all fields. Key sections:
- `target`: sequence + hotspot residues (1-indexed)
- `binders`: list of {name, sequence, iptm}
- `scaffold`: LOV2 sequences, cap, linker selection
- `prediction`: Boltz-2 parameters
- `binder_design`: RFdiffusion/LigandMPNN parameters (for `screen_binders.py`)
