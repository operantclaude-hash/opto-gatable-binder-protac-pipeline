#!/usr/bin/env python3
"""
End-to-end pipeline orchestrator.

Chains all pipeline steps with checkpoint tracking so you can resume.

Usage:
    # Phase 1: Generate binder design scripts
    python3 run_pipeline.py campaign_config.json

    # Phase 2: After binder screening completes, score + setup opto-scaffold
    python3 run_pipeline.py campaign_config.json --continue

    # Phase 3: After opto-scaffold Boltz-2 predictions complete, analyze
    python3 run_pipeline.py campaign_config.json --analyze
"""

import argparse
import json
import sys
from pathlib import Path

PIPELINE_DIR = Path(__file__).parent


def get_checkpoint(campaign_dir):
    """Read checkpoint file to determine pipeline state."""
    cp_path = campaign_dir / ".pipeline_checkpoint"
    if cp_path.exists():
        return json.loads(cp_path.read_text())
    return {"phase": 0}


def set_checkpoint(campaign_dir, phase, **kwargs):
    """Write checkpoint file."""
    cp_path = campaign_dir / ".pipeline_checkpoint"
    data = {"phase": phase, **kwargs}
    cp_path.write_text(json.dumps(data, indent=2) + "\n")


def phase1_design(config_path, n_gpus=None):
    """Phase 1: Generate binder screening scripts."""
    from screen_binders import screen_binders

    with open(config_path) as f:
        cfg = json.load(f)

    name = cfg["campaign_name"]
    campaign_dir = PIPELINE_DIR / "campaigns" / name

    print("=" * 60)
    print("PHASE 1: Binder Design")
    print("=" * 60)

    screen_dir = screen_binders(config_path, n_gpus)

    set_checkpoint(campaign_dir, 1, screen_dir=str(screen_dir))

    print()
    print("=" * 60)
    print("ACTION REQUIRED")
    print("=" * 60)
    print(f"  Run the binder screening scripts:")
    print(f"    bash {screen_dir / 'run_all.sh'}")
    print()
    print(f"  When complete, continue with:")
    print(f"    python3 run_pipeline.py {config_path} --continue")


def phase2_opto(config_path, n_gpus=None):
    """Phase 2: Score binders, setup opto-scaffold campaign."""
    from score_binders import score_binder_screen, generate_opto_config
    from setup_campaign import setup_campaign

    with open(config_path) as f:
        cfg = json.load(f)

    name = cfg["campaign_name"]
    campaign_dir = PIPELINE_DIR / "campaigns" / name
    screen_dir = campaign_dir / "binder_screen"

    print("=" * 60)
    print("PHASE 2: Score Binders + Setup Opto-Scaffold")
    print("=" * 60)

    # Check binder screen results exist
    results_dir = screen_dir / "boltz2_results"
    if not results_dir.exists() or not any(results_dir.iterdir()):
        print(f"ERROR: No binder screen results in {results_dir}")
        print("Run the binder screening scripts first (Phase 1).")
        sys.exit(1)

    # Score binders
    hotspots = cfg.get("target", {}).get("hotspot_residues_1idx")
    top_n = cfg.get("binder_design", {}).get("top_n_binders", 10)
    all_results, top_binders = score_binder_screen(screen_dir, hotspots, top_n)

    if not top_binders:
        print("ERROR: No binders passed scoring.")
        sys.exit(1)

    # Save ranked results
    ranked_path = screen_dir / "ranked_binders.json"
    ranked_path.write_text(json.dumps(all_results, indent=2) + "\n")

    # Generate opto config
    opto_config_path = campaign_dir / "campaign_config_opto.json"
    generate_opto_config(top_binders, str(config_path), str(opto_config_path))

    # Setup opto-scaffold campaign
    print()
    setup_campaign(str(opto_config_path))

    opto_name = f"{name}_opto"
    opto_dir = PIPELINE_DIR / "campaigns" / opto_name
    set_checkpoint(campaign_dir, 2, opto_dir=str(opto_dir))

    run_script = opto_dir / "run_predictions.sh"
    print()
    print("=" * 60)
    print("ACTION REQUIRED")
    print("=" * 60)
    print(f"  Run opto-scaffold Boltz-2 predictions:")
    print(f"    bash {run_script}")
    print()
    print(f"  When complete, analyze with:")
    print(f"    python3 run_pipeline.py {config_path} --analyze")


def phase3_analyze(config_path):
    """Phase 3: Analyze opto-scaffold predictions."""
    from analyze_campaign import analyze_campaign

    with open(config_path) as f:
        cfg = json.load(f)

    name = cfg["campaign_name"]
    campaign_dir = PIPELINE_DIR / "campaigns" / name
    checkpoint = get_checkpoint(campaign_dir)

    opto_dir = checkpoint.get("opto_dir")
    if not opto_dir:
        opto_dir = str(PIPELINE_DIR / "campaigns" / f"{name}_opto")

    print("=" * 60)
    print("PHASE 3: Analyze Opto-Scaffold Campaign")
    print("=" * 60)

    if not Path(opto_dir).exists():
        print(f"ERROR: Opto campaign directory not found: {opto_dir}")
        sys.exit(1)

    analyze_campaign(opto_dir)

    set_checkpoint(campaign_dir, 3)

    print()
    print("=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"  Results: {Path(opto_dir) / 'analysis/'}")
    print(f"  Ranked constructs: {Path(opto_dir) / 'analysis' / 'ranked_constructs.json'}")
    print(f"  Power analysis: {Path(opto_dir) / 'analysis' / 'power_analysis.json'}")
    print()
    print("  Optional next steps:")
    print(f"    python3 assemble_degron.py {campaign_dir / 'campaign_config_opto.json'}")


def main():
    parser = argparse.ArgumentParser(
        description="End-to-end opto-gatable binder pipeline"
    )
    parser.add_argument("config", help="Campaign config JSON")
    parser.add_argument("--continue", dest="continue_phase", action="store_true",
                        help="Continue to Phase 2 (score + opto-scaffold setup)")
    parser.add_argument("--analyze", action="store_true",
                        help="Run Phase 3 analysis")
    parser.add_argument("--n_gpus", type=int, default=None,
                        help="Override number of GPUs")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: Config not found: {config_path}")
        sys.exit(1)

    if args.analyze:
        phase3_analyze(str(config_path))
    elif args.continue_phase:
        phase2_opto(str(config_path), args.n_gpus)
    else:
        phase1_design(str(config_path), args.n_gpus)


if __name__ == "__main__":
    main()
