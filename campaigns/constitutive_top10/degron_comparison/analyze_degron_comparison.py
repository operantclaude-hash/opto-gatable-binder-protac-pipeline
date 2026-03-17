#!/usr/bin/env python3
"""
Analyze degron comparison results: N-term vs C-term placement.

Extracts ipTM, pLDDT, and confidence from Boltz-2 predictions to determine
which degron placement preserves the binder-target interface better.
"""

import json
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).parent / "boltz2_results"

CONSTRUCTS = [
    "degron_nterm_dark_target",
    "degron_nterm_lit_target",
    "degron_cterm_dark_target",
    "degron_cterm_lit_target",
    "degron_nterm_dark_fold",
    "degron_nterm_lit_fold",
    "degron_cterm_dark_fold",
    "degron_cterm_lit_fold",
]

# Domain boundaries (1-indexed) for binder residues in chain A
BINDER_RANGES = {
    "nterm_dark": (425, 474),
    "nterm_lit": (425, 474),
    "cterm_dark": (1, 50),
    "cterm_lit": (1, 50),
}


def extract_metrics(name):
    """Extract confidence metrics from a Boltz-2 result directory."""
    rdir = RESULTS_DIR / name
    if not rdir.exists():
        return None

    conf_files = sorted(rdir.rglob("confidence_*_model_0.json"))
    if not conf_files:
        conf_files = sorted(rdir.rglob("confidence_*.json"))
    if not conf_files:
        return None

    with open(conf_files[0]) as f:
        conf = json.load(f)

    metrics = {
        "iptm": conf.get("protein_iptm", conf.get("iptm", 0.0)),
        "ptm": conf.get("ptm", 0.0),
        "plddt": conf.get("complex_plddt", 0.0),
        "confidence": conf.get("confidence_score", 0.0),
    }

    # Extract cross-chain PAE for 2-chain constructs
    pae_files = sorted(rdir.rglob("pae_*_model_0.npz"))
    if not pae_files:
        pae_files = sorted(rdir.rglob("pae_*.npz"))

    if pae_files and "_target" in name:
        pae_data = np.load(pae_files[0])
        pae_key = list(pae_data.keys())[0]
        pae_matrix = pae_data[pae_key]

        # Chain A length varies by construct
        key = name.replace("degron_", "").replace("_target", "")
        if "dark" in key:
            chain_a_len = 715
        else:
            chain_a_len = 700

        if pae_matrix.shape[0] >= chain_a_len:
            pae_cross = np.concatenate([
                pae_matrix[:chain_a_len, chain_a_len:].flatten(),
                pae_matrix[chain_a_len:, :chain_a_len].flatten(),
            ])
            metrics["min_ipae"] = float(np.min(pae_cross))
            metrics["mean_ipae"] = float(np.mean(pae_cross))

            # Binder-specific cross-chain PAE
            binder_start, binder_end = BINDER_RANGES[key]
            binder_to_target = pae_matrix[binder_start-1:binder_end, chain_a_len:]
            target_to_binder = pae_matrix[chain_a_len:, binder_start-1:binder_end]
            binder_cross = np.concatenate([binder_to_target.flatten(),
                                           target_to_binder.flatten()])
            metrics["binder_min_ipae"] = float(np.min(binder_cross))
            metrics["binder_mean_ipae"] = float(np.mean(binder_cross))

    return metrics


def main():
    print("=" * 80)
    print("DEGRON PLACEMENT COMPARISON: N-TERMINAL vs C-TERMINAL")
    print("Binder: pep4s2 (50aa)  |  E3: TRIM21 RBCC  |  VVD: 158aa")
    print("=" * 80)

    results = {}
    for name in CONSTRUCTS:
        metrics = extract_metrics(name)
        results[name] = metrics
        if metrics:
            print(f"  {name}: ipTM={metrics['iptm']:.3f} pLDDT={metrics['plddt']:.3f}", end="")
            if "min_ipae" in metrics:
                print(f" min_iPAE={metrics['min_ipae']:.2f}"
                      f" binder_iPAE={metrics['binder_min_ipae']:.2f}", end="")
            print()
        else:
            print(f"  {name}: NO RESULTS")

    # Comparison table
    print()
    print("=" * 80)
    print("2-CHAIN COMPARISON (construct + p300 HAT target)")
    print("=" * 80)
    print(f"{'Construct':<32} {'ipTM':<8} {'pLDDT':<8} {'min_iPAE':<10} {'binder_iPAE':<12}")
    print("-" * 80)

    for name in CONSTRUCTS:
        if "_target" not in name:
            continue
        m = results.get(name)
        if not m:
            print(f"{name:<32} ---")
            continue
        print(f"{name:<32} {m['iptm']:.3f}    {m['plddt']:.3f}    "
              f"{m.get('min_ipae', 0):.2f}       {m.get('binder_min_ipae', 0):.2f}")

    print()
    print("=" * 80)
    print("FOLD-ONLY COMPARISON (single chain, structural quality)")
    print("=" * 80)
    print(f"{'Construct':<32} {'pTM':<8} {'pLDDT':<8} {'confidence':<12}")
    print("-" * 80)

    for name in CONSTRUCTS:
        if "_fold" not in name:
            continue
        m = results.get(name)
        if not m:
            print(f"{name:<32} ---")
            continue
        print(f"{name:<32} {m['ptm']:.3f}    {m['plddt']:.3f}    {m['confidence']:.3f}")

    # Verdict
    print()
    print("=" * 80)
    print("VERDICT")
    print("=" * 80)

    nterm_dark = results.get("degron_nterm_dark_target")
    cterm_dark = results.get("degron_cterm_dark_target")
    nterm_lit = results.get("degron_nterm_lit_target")
    cterm_lit = results.get("degron_cterm_lit_target")

    if nterm_dark and cterm_dark:
        if nterm_dark["iptm"] > cterm_dark["iptm"]:
            print("  N-terminal placement shows BETTER binding (higher ipTM)")
        elif cterm_dark["iptm"] > nterm_dark["iptm"]:
            print("  C-terminal placement shows BETTER binding (higher ipTM)")
        else:
            print("  Both placements show SIMILAR binding")

        binder_n = nterm_dark.get("binder_min_ipae", float("inf"))
        binder_c = cterm_dark.get("binder_min_ipae", float("inf"))
        if binder_n < binder_c:
            print("  N-terminal has BETTER binder-target interface (lower binder iPAE)")
        elif binder_c < binder_n:
            print("  C-terminal has BETTER binder-target interface (lower binder iPAE)")

    # Save results
    out_path = Path(__file__).parent / "degron_comparison_results.json"
    out_path.write_text(json.dumps(results, indent=2, default=str) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
