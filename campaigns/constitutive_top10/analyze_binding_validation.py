#!/usr/bin/env python3
"""
Analyze binding validation results: 22 binders x 4 conditions.

Extracts ipTM, min_iPAE, ipDE from Boltz-2 confidence outputs and
PAE matrices. Produces:
  1. Per-condition comparison tables (ipTM and min_iPAE)
  2. Holistic ranking table with all properties for final binder selection
  3. Summary statistics on integration and context effects
"""

import json
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

SCRIPT_DIR = Path(__file__).parent
VALIDATION_DIR = SCRIPT_DIR / "binding_validation"
RESULTS_DIR = VALIDATION_DIR / "boltz2_results"
MANIFEST_PATH = VALIDATION_DIR / "manifest.json"
OUTPUT_DIR = VALIDATION_DIR / "analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def extract_confidence(result_dir, construct_name):
    """Extract confidence metrics from Boltz-2 results.

    Returns dict with ipTM, pTM, pLDDT, ipDE, ipLDDT, min_iPAE, mean_iPAE.
    """
    # Find confidence JSON (best model = model_0)
    conf_files = sorted(result_dir.rglob("confidence_*_model_0.json"))
    if not conf_files:
        # Try any model
        conf_files = sorted(result_dir.rglob("confidence_*.json"))
    if not conf_files:
        return None

    with open(conf_files[0]) as f:
        conf = json.load(f)

    metrics = {
        "iptm": conf.get("protein_iptm", conf.get("iptm", 0.0)),
        "ptm": conf.get("ptm", 0.0),
        "plddt": conf.get("complex_plddt", 0.0),
        "ipde": conf.get("complex_ipde", 0.0),
        "iplddt": conf.get("complex_iplddt", 0.0),
        "confidence": conf.get("confidence_score", 0.0),
    }

    # Extract min/mean iPAE from PAE matrix
    pae_files = sorted(result_dir.rglob("pae_*_model_0.npz"))
    if not pae_files:
        pae_files = sorted(result_dir.rglob("pae_*.npz"))

    if pae_files:
        pae_data = np.load(pae_files[0])
        pae_key = list(pae_data.keys())[0]
        pae_matrix = pae_data[pae_key]

        # Get chain lengths from manifest to extract cross-chain block
        manifest_entry = None
        if MANIFEST_PATH.exists():
            with open(MANIFEST_PATH) as f:
                manifest = json.load(f)
            for entry in manifest:
                if entry["construct_name"] == construct_name:
                    manifest_entry = entry
                    break

        if manifest_entry:
            chain_a_len = manifest_entry["chain_a_length"]
            # Cross-chain PAE blocks
            pae_a_to_b = pae_matrix[:chain_a_len, chain_a_len:]
            pae_b_to_a = pae_matrix[chain_a_len:, :chain_a_len]
            cross_chain = np.concatenate([pae_a_to_b.flatten(), pae_b_to_a.flatten()])

            metrics["min_ipae"] = float(np.min(cross_chain))
            metrics["mean_ipae"] = float(np.mean(cross_chain))
        else:
            metrics["min_ipae"] = None
            metrics["mean_ipae"] = None
    else:
        metrics["min_ipae"] = None
        metrics["mean_ipae"] = None

    return metrics


def main():
    if not MANIFEST_PATH.exists():
        print(f"Manifest not found: {MANIFEST_PATH}")
        print("Run generate_binding_validation.py first.")
        return

    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)

    print(f"Analyzing {len(manifest)} predictions...")
    print()

    results = []
    for entry in manifest:
        name = entry["construct_name"]
        result_dir = RESULTS_DIR / name

        if not result_dir.exists():
            print(f"  {name}: no results")
            results.append({**entry, "metrics": None})
            continue

        metrics = extract_confidence(result_dir, name)
        if metrics:
            results.append({**entry, "metrics": metrics})
            print(f"  {name}: ipTM={metrics['iptm']:.3f} min_iPAE={metrics.get('min_ipae', 'N/A')}")
        else:
            results.append({**entry, "metrics": None})
            print(f"  {name}: no confidence data")

    # Build comparison tables
    binder_names = list(dict.fromkeys(e["binder_name"] for e in manifest))
    conditions = ["binder_HAT", "fusion_HAT", "binder_catcore", "fusion_catcore"]

    # ipTM table
    print("\n" + "=" * 100)
    print("ipTM COMPARISON (higher = better binding confidence)")
    print("=" * 100)
    header = f"{'Binder':<20} {'binder+HAT':<12} {'fusion+HAT':<12} {'delta':<8} {'binder+cat':<12} {'fusion+cat':<12} {'delta':<8}"
    print(header)
    print("-" * 100)

    summary_data = []

    for bname in binder_names:
        row = {"binder": bname}
        vals = {}
        for cond in conditions:
            match = [r for r in results if r["binder_name"] == bname and r["condition"] == cond]
            if match and match[0]["metrics"]:
                vals[cond] = match[0]["metrics"]["iptm"]
            else:
                vals[cond] = None

        bh = vals.get("binder_HAT")
        fh = vals.get("fusion_HAT")
        bc = vals.get("binder_catcore")
        fc = vals.get("fusion_catcore")

        d1 = (fh - bh) if (bh is not None and fh is not None) else None
        d2 = (fc - bc) if (bc is not None and fc is not None) else None

        row["binder_HAT"] = bh
        row["fusion_HAT"] = fh
        row["fusion_delta_HAT"] = d1
        row["binder_catcore"] = bc
        row["fusion_catcore"] = fc
        row["fusion_delta_catcore"] = d2
        summary_data.append(row)

        def fmt(v):
            return f"{v:.3f}" if v is not None else "---"
        def fmt_delta(v):
            if v is None:
                return "---"
            sign = "+" if v >= 0 else ""
            return f"{sign}{v:.3f}"

        print(f"{bname:<20} {fmt(bh):<12} {fmt(fh):<12} {fmt_delta(d1):<8} "
              f"{fmt(bc):<12} {fmt(fc):<12} {fmt_delta(d2):<8}")

    # min_iPAE table
    print("\n" + "=" * 100)
    print("min_iPAE COMPARISON (lower = better interface confidence)")
    print("=" * 100)
    print(header)
    print("-" * 100)

    for bname in binder_names:
        vals = {}
        for cond in conditions:
            match = [r for r in results if r["binder_name"] == bname and r["condition"] == cond]
            if match and match[0]["metrics"] and match[0]["metrics"].get("min_ipae") is not None:
                vals[cond] = match[0]["metrics"]["min_ipae"]
            else:
                vals[cond] = None

        bh = vals.get("binder_HAT")
        fh = vals.get("fusion_HAT")
        bc = vals.get("binder_catcore")
        fc = vals.get("fusion_catcore")

        d1 = (fh - bh) if (bh is not None and fh is not None) else None
        d2 = (fc - bc) if (bc is not None and fc is not None) else None

        def fmt(v):
            return f"{v:.3f}" if v is not None else "---"
        def fmt_delta(v):
            if v is None:
                return "---"
            sign = "+" if v >= 0 else ""
            return f"{sign}{v:.3f}"

        print(f"{bname:<20} {fmt(bh):<12} {fmt(fh):<12} {fmt_delta(d1):<8} "
              f"{fmt(bc):<12} {fmt(fc):<12} {fmt_delta(d2):<8}")

    # Summary statistics
    complete = [r for r in results if r["metrics"] is not None]
    if complete:
        print("\n" + "=" * 100)
        print("SUMMARY")
        print("=" * 100)

        for cond in conditions:
            cond_results = [r for r in complete if r["condition"] == cond]
            if cond_results:
                iptms = [r["metrics"]["iptm"] for r in cond_results]
                print(f"{cond:<20} n={len(cond_results):<3} "
                      f"ipTM: mean={np.mean(iptms):.3f} "
                      f"min={np.min(iptms):.3f} max={np.max(iptms):.3f}")

        # Compute average drop from binder→fusion
        hat_drops = []
        cat_drops = []
        for bname in binder_names:
            bh = next((r["metrics"]["iptm"] for r in complete
                       if r["binder_name"] == bname and r["condition"] == "binder_HAT"), None)
            fh = next((r["metrics"]["iptm"] for r in complete
                       if r["binder_name"] == bname and r["condition"] == "fusion_HAT"), None)
            bc = next((r["metrics"]["iptm"] for r in complete
                       if r["binder_name"] == bname and r["condition"] == "binder_catcore"), None)
            fc = next((r["metrics"]["iptm"] for r in complete
                       if r["binder_name"] == bname and r["condition"] == "fusion_catcore"), None)
            if bh is not None and fh is not None:
                hat_drops.append(fh - bh)
            if bc is not None and fc is not None:
                cat_drops.append(fc - bc)

        if hat_drops:
            print(f"\nAvg ipTM change (binder→fusion, HAT):     {np.mean(hat_drops):+.3f}")
        if cat_drops:
            print(f"Avg ipTM change (binder→fusion, catcore):  {np.mean(cat_drops):+.3f}")

    # ===== RAW DATA TABLE — ALL PROPERTIES =====
    print("\n" + "=" * 150)
    print("ALL BINDERS — RAW DATA (sorted by Stage 1 ipTM)")
    print("=" * 150)
    print(f"{'Binder':<12} {'Len':<5} "
          f"{'S1_ipTM':<8} {'S1_iPAE':<8} "
          f"{'b+HAT':<8} {'f+HAT':<8} {'intg':<7} "
          f"{'b+cat':<8} {'f+cat':<8} {'ctx':<7} "
          f"{'iPAE_bH':<8} {'iPAE_fH':<8} {'iPAE_fC':<8}")
    print("-" * 150)

    holistic = []
    for bname in binder_names:
        binder_entry = next((e for e in manifest if e["binder_name"] == bname), {})
        blen = binder_entry.get("binder_length", 0)
        s1_iptm = binder_entry.get("binder_iptm_stage1", 0)
        s1_ipae = binder_entry.get("binder_min_ipae_stage1")

        metrics_by_cond = {}
        for cond in conditions:
            match = [r for r in results if r["binder_name"] == bname and r["condition"] == cond]
            if match and match[0]["metrics"]:
                metrics_by_cond[cond] = match[0]["metrics"]

        bh_iptm = metrics_by_cond.get("binder_HAT", {}).get("iptm")
        fh_iptm = metrics_by_cond.get("fusion_HAT", {}).get("iptm")
        bc_iptm = metrics_by_cond.get("binder_catcore", {}).get("iptm")
        fc_iptm = metrics_by_cond.get("fusion_catcore", {}).get("iptm")

        bh_ipae = metrics_by_cond.get("binder_HAT", {}).get("min_ipae")
        fh_ipae = metrics_by_cond.get("fusion_HAT", {}).get("min_ipae")
        fc_ipae = metrics_by_cond.get("fusion_catcore", {}).get("min_ipae")

        # Integration delta (binder→fusion, HAT context)
        intg_delta = (fh_iptm - bh_iptm) if (fh_iptm and bh_iptm) else None
        # Context delta (HAT→catcore, binder alone)
        ctx_delta = (bc_iptm - bh_iptm) if (bc_iptm and bh_iptm) else None

        holistic.append({
            "binder": bname,
            "length": blen,
            "s1_iptm": s1_iptm,
            "s1_min_ipae": s1_ipae,
            "binder_HAT_iptm": bh_iptm,
            "fusion_HAT_iptm": fh_iptm,
            "integration_delta": intg_delta,
            "binder_catcore_iptm": bc_iptm,
            "fusion_catcore_iptm": fc_iptm,
            "context_delta": ctx_delta,
            "binder_HAT_min_ipae": bh_ipae,
            "fusion_HAT_min_ipae": fh_ipae,
            "fusion_catcore_min_ipae": fc_ipae,
        })

    # Sort by Stage 1 ipTM (the metric with most literature backing)
    holistic.sort(key=lambda x: x["s1_iptm"] or 0, reverse=True)

    def f(v):
        return f"{v:.3f}" if v is not None else "---"
    def fd(v):
        if v is None: return "---"
        sign = "+" if v >= 0 else ""
        return f"{sign}{v:.3f}"

    for h in holistic:
        print(f"{h['binder']:<12} {h['length']:<5} "
              f"{f(h['s1_iptm']):<8} {f(h['s1_min_ipae']):<8} "
              f"{f(h['binder_HAT_iptm']):<8} {f(h['fusion_HAT_iptm']):<8} {fd(h['integration_delta']):<7} "
              f"{f(h['binder_catcore_iptm']):<8} {f(h['fusion_catcore_iptm']):<8} {fd(h['context_delta']):<7} "
              f"{f(h['binder_HAT_min_ipae']):<8} {f(h['fusion_HAT_min_ipae']):<8} {f(h['fusion_catcore_min_ipae']):<8}")

    print()
    print("Column legend:")
    print("  S1_ipTM / S1_iPAE = Original Stage 1 binder-only metrics (binder + HAT, no fusion)")
    print("  b+HAT / f+HAT     = ipTM for binder-alone / fusion vs HAT domain (this run)")
    print("  intg              = integration delta: f+HAT minus b+HAT (negative = fusion hurts)")
    print("  b+cat / f+cat     = ipTM for binder-alone / fusion vs catalytic core (this run)")
    print("  ctx               = context delta: b+cat minus b+HAT (negative = catcore hurts)")
    print("  iPAE_bH/fH/fC     = min_iPAE for binder+HAT / fusion+HAT / fusion+catcore")

    # Save everything
    output = {
        "raw_data": holistic,
        "condition_tables": summary_data,
        "all_results": [{k: v for k, v in r.items()} for r in results],
    }
    output_path = OUTPUT_DIR / "binding_validation_results.json"
    output_path.write_text(json.dumps(output, indent=2, default=str) + "\n")
    print(f"\nResults saved: {output_path}")


if __name__ == "__main__":
    main()
