#!/usr/bin/env python3
"""
Integrated multi-layer ranking for top-10 binder selection.

Combines all validation layers:
  Layer 1: Epitope contact analysis (hotspot overlap, binder-alone)
  Layer 2: Cap-truncated fusion (binder domain contacts WITH scaffold)
  Layer 3: Multi-sample consistency (reproducibility across diffusion samples)
  Layer 4: Interface quality (pLDDT, PAE, ipTM stability)

Outputs a ranked table with per-layer scores and integrated confidence.
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent

# Load all layer results
with open(BASE / "epitope_contact_results.json") as f:
    epitope_data = json.load(f)

with open(BASE / "multisample_consistency_results.json") as f:
    consistency_data = json.load(f)

with open(BASE / "interface_quality_results.json") as f:
    quality_data = json.load(f)

with open(BASE / "binding_validation" / "analysis" / "binding_validation_results.json") as f:
    bv_data = json.load(f)

# Index by binder name
epitope_map = {r["binder"]: r for r in epitope_data}
consistency_map = {r["binder"]: r for r in consistency_data}
quality_map = {r["binder"]: r for r in quality_data}
bv_map = {r["binder"]: r for r in bv_data["raw_data"]}

# Hotspot residues
HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}

# All binders to rank (peptides only — nanobodies all off-target)
PEPTIDE_BINDERS = [
    "pep30s5", "pep25s1", "pep21s5", "pep15s3", "pep48s5",
    "pep16s5", "pep36s1", "pep29s5", "pep21s2", "pep24s3",
    "pep13s3", "pep17s4", "pep4s2", "pep10s5", "pep32s1",
    "pep32s2",
]


def compute_layer_scores(binder):
    """Compute normalized scores for each validation layer."""
    scores = {}

    # --- Layer 1: Epitope contact (0-1 scale) ---
    ep = epitope_map.get(binder, {})
    hotspot_pct = ep.get("hotspot_pct", 0) / 100.0  # normalize to 0-1
    scores["L1_hotspot_pct"] = hotspot_pct
    scores["L1_verdict"] = ep.get("verdict", "UNKNOWN")

    # --- Layer 2: Cap-truncated fusion (binder domain contacts WITH scaffold) ---
    ct_path = BASE / "captrunc_validation" / "captrunc_results.json"
    ct_map = {}
    if ct_path.exists():
        with open(ct_path) as f:
            for entry in json.load(f):
                ct_map[entry["binder"]] = entry

    ct = ct_map.get(binder, {})
    mean_bdr_hs = ct.get("mean_bdr_hotspots", 0)
    driver = ct.get("driver", "UNKNOWN")

    # Score based on whether binder domain specifically drives contacts
    # Baseline: scaffold_only control shows ~3-4 hotspot contacts from scaffold alone
    if driver == "BINDER-DRIVEN":
        scores["L2_captrunc"] = 1.0
    elif driver == "MIXED" and mean_bdr_hs >= 2.0:
        scores["L2_captrunc"] = 0.7
    elif driver == "SCAFFOLD-DOM" and mean_bdr_hs > 0:
        scores["L2_captrunc"] = 0.3
    elif driver == "SCAFFOLD-ONLY":
        scores["L2_captrunc"] = 0.1  # Binder doesn't reach target when on scaffold
    elif scores["L1_verdict"] == "OFF-TARGET":
        scores["L2_captrunc"] = 0.0
    else:
        scores["L2_captrunc"] = 0.2  # No cap-truncated data
    scores["L2_driver"] = driver
    scores["L2_bdr_hotspots"] = mean_bdr_hs

    # --- Layer 3: Multi-sample consistency ---
    con = consistency_map.get(binder, {})
    consistency = con.get("consistency", 0)
    n_consistent = con.get("n_consistent_hotspots", 0)
    mean_hs = con.get("mean_hotspots", 0)

    # Composite: weight consistency, number of consistent hotspots, mean count
    scores["L3_consistency"] = consistency
    scores["L3_n_consistent"] = n_consistent
    scores["L3_mean_hotspots"] = mean_hs
    # Normalized score: consistent hotspots / 14 * consistency
    scores["L3_score"] = (n_consistent / 14.0) * consistency

    # --- Layer 4: Interface quality ---
    q = quality_map.get(binder, {})
    best_iptm = q.get("best_iptm", 0)
    mean_iptm = q.get("mean_iptm", 0)
    std_iptm = q.get("std_iptm", 1.0)
    iface_plddt = q.get("mean_interface_plddt", 0)
    cross_pae = q.get("mean_cross_pae", 20.0)

    # Stability: lower std = better (normalized: 1 - std/0.3, clamped)
    stability = max(0.0, min(1.0, 1.0 - std_iptm / 0.3))

    # PAE quality: lower = better (normalized: 1 - pae/20, clamped)
    pae_quality = max(0.0, min(1.0, 1.0 - (cross_pae if cross_pae else 20.0) / 20.0))

    # pLDDT quality (on 0-100 scale from CIF): normalize to 0-1
    plddt_quality = min(1.0, iface_plddt / 100.0) if iface_plddt > 1 else iface_plddt

    scores["L4_best_iptm"] = best_iptm
    scores["L4_stability"] = stability
    scores["L4_plddt"] = plddt_quality
    scores["L4_pae"] = pae_quality
    scores["L4_score"] = 0.4 * best_iptm + 0.2 * stability + 0.2 * plddt_quality + 0.2 * pae_quality

    # --- Layer 5: CatCore epitope agreement ---
    cc_path = BASE / "catcore_epitope_results.json"
    cc_map = {}
    if cc_path.exists():
        with open(cc_path) as f:
            for entry in json.load(f):
                cc_map[entry["binder"]] = entry

    cc = cc_map.get(binder, {})
    cc_max_hs = cc.get("max_hotspot", 0)
    cc_avg_hs = cc.get("avg_hotspot", 0)
    cc_n_consistent = cc.get("n_consistent", 0)
    cc_best_iptm = cc.get("best_iptm", 0)

    # Catcore agreement: is the binder also on-target against full catcore?
    # Score based on consistent catcore hotspot contacts
    cc_score = min(1.0, cc_n_consistent / 8.0)  # 8+ consistent = perfect
    # Bonus for high catcore ipTM
    if cc_best_iptm >= 0.80:
        cc_score = min(1.0, cc_score + 0.2)
    elif cc_best_iptm < 0.40:
        cc_score *= 0.5  # Penalize catcore failures

    scores["L5_cc_max_hs"] = cc_max_hs
    scores["L5_cc_avg_hs"] = cc_avg_hs
    scores["L5_cc_consistent"] = cc_n_consistent
    scores["L5_cc_iptm"] = cc_best_iptm
    scores["L5_score"] = cc_score

    return scores


def compute_integrated_score(scores):
    """Compute final integrated confidence score (0-1)."""
    # Weights reflect importance for predicting empirical binding success
    # Since we're testing against full p300 (HAT + catcore), catcore agreement matters
    return (
        0.20 * scores["L1_hotspot_pct"] +       # Correct epitope (HAT alone)
        0.10 * scores["L2_captrunc"] +            # Scaffold integration
        0.20 * scores["L3_score"] +               # Reproducibility
        0.25 * scores["L4_score"] +               # Interface quality/stability
        0.25 * scores.get("L5_score", 0)          # CatCore agreement (critical for full p300)
    )


def main():
    results = []

    for binder in PEPTIDE_BINDERS:
        bv = bv_map.get(binder, {})
        scores = compute_layer_scores(binder)
        integrated = compute_integrated_score(scores)

        results.append({
            "binder": binder,
            "length": bv.get("length", 0),
            "screen_iptm": bv.get("s1_iptm", 0),
            "binder_HAT_iptm": bv.get("binder_HAT_iptm", 0),
            "binder_catcore_iptm": bv.get("binder_catcore_iptm", 0),
            "scores": scores,
            "integrated": integrated,
        })

    results.sort(key=lambda x: x["integrated"], reverse=True)

    print("=" * 140)
    print("INTEGRATED MULTI-LAYER RANKING — TOP BINDERS FOR EMPIRICAL TESTING")
    print("Target: p300 HAT domain | Construct: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]-[UBQ]")
    print("=" * 140)
    print()

    print("Layer weights: L1(epitope)=30%, L2(cap-trunc)=10%, L3(consistency)=30%, L4(quality)=30%")
    print()

    hdr = (f"{'Rk':<4} {'Binder':<10} {'Len':<5} "
           f"{'L1:HS%':<8} {'L2:CapTr':<10} "
           f"{'L3:nHS':<7} {'L4:ipTM':<8} "
           f"{'L5:CC_HS':<9} {'L5:CC_ipTM':<11} "
           f"{'SCORE':<8} {'TIER'}")
    print(hdr)
    print("-" * 140)

    for i, r in enumerate(results, 1):
        s = r["scores"]
        tier = ""
        if r["integrated"] >= 0.50:
            tier = "TIER 1"
        elif r["integrated"] >= 0.35:
            tier = "TIER 2"
        elif r["integrated"] >= 0.20:
            tier = "TIER 3"
        else:
            tier = "EXCLUDE"

        cc_hs = f"{s.get('L5_cc_consistent', 0)}/{s.get('L5_cc_max_hs', 0)}"
        cc_iptm = f"{s.get('L5_cc_iptm', 0):.3f}" if s.get('L5_cc_iptm', 0) > 0 else "—"

        print(f"{i:<4} {r['binder']:<10} {r['length']:>3}  "
              f"{s['L1_hotspot_pct']:>5.1%}   {s.get('L2_driver', '?'):<10} "
              f"{s['L3_n_consistent']:>3}    {s['L4_best_iptm']:>6.3f}  "
              f"{cc_hs:<9} {cc_iptm:<11} "
              f"{r['integrated']:>6.3f}  {tier}")

    # Tier summaries
    tier1 = [r for r in results if r["integrated"] >= 0.50]
    tier2 = [r for r in results if 0.35 <= r["integrated"] < 0.50]
    tier3 = [r for r in results if 0.20 <= r["integrated"] < 0.35]
    excluded = [r for r in results if r["integrated"] < 0.20]

    print()
    print("=" * 140)
    print("TIER SUMMARY")
    print("=" * 140)

    print(f"\n  TIER 1 — HIGH CONFIDENCE ({len(tier1)} binders):")
    for r in tier1:
        print(f"    {r['binder']:<10} {r['length']:>3}aa  integrated={r['integrated']:.3f}  "
              f"screen={r['screen_iptm']:.3f}  HAT_val={r['binder_HAT_iptm']:.3f}  "
              f"catcore={r['binder_catcore_iptm']:.3f}")

    print(f"\n  TIER 2 — MODERATE CONFIDENCE ({len(tier2)} binders):")
    for r in tier2:
        print(f"    {r['binder']:<10} {r['length']:>3}aa  integrated={r['integrated']:.3f}  "
              f"screen={r['screen_iptm']:.3f}  HAT_val={r['binder_HAT_iptm']:.3f}  "
              f"catcore={r['binder_catcore_iptm']:.3f}")

    print(f"\n  TIER 3 — LOW CONFIDENCE ({len(tier3)} binders):")
    for r in tier3:
        print(f"    {r['binder']:<10} {r['length']:>3}aa  integrated={r['integrated']:.3f}  "
              f"HAT_val={r['binder_HAT_iptm']:.3f}")

    print(f"\n  EXCLUDED ({len(excluded)} binders):")
    for r in excluded:
        print(f"    {r['binder']:<10} — {r['scores']['L1_verdict']}, integrated={r['integrated']:.3f}")

    # Recommendations
    print()
    print("=" * 140)
    print("RECOMMENDATIONS FOR TOP 10 SELECTION")
    print("=" * 140)

    if len(tier1) >= 10:
        print(f"\n  We have {len(tier1)} Tier 1 binders — enough for a full 10.")
    elif len(tier1) + len(tier2) >= 10:
        print(f"\n  {len(tier1)} Tier 1 + {len(tier2)} Tier 2 = {len(tier1) + len(tier2)} candidates.")
        print("  This is enough for 10 selections, but confidence is mixed.")
    else:
        total = len(tier1) + len(tier2)
        shortfall = 10 - total
        print(f"\n  CRITICAL: Only {total} Tier 1+2 binders ({len(tier1)} + {len(tier2)}).")
        print(f"  SHORT BY {shortfall} for a full top-10 selection.")
        print()
        print("  OPTIONS:")
        print("    1. Run expanded pool validation (16 additional untested binders with ipTM>0.85)")
        print("       These all have >= 10/14 hotspot contacts in the original screen")
        print("       YAMLs ready at: expanded_pool/boltz2_yamls/")
        print()
        print("    2. Run additional RFdiffusion + LigandMPNN screening")
        print("       Could generate 200+ new binder designs per day")
        print()
        print("    3. Accept Tier 3 binders with caveats")

    # Save
    out_path = BASE / "integrated_ranking_results.json"
    save_data = {
        "description": "Integrated multi-layer ranking for binder selection",
        "layers": {
            "L1": "Epitope contact analysis (hotspot overlap)",
            "L2": "Cap-truncated fusion (scaffold integration)",
            "L3": "Multi-sample consistency (reproducibility)",
            "L4": "Interface quality (pLDDT, PAE, stability)",
        },
        "weights": {"L1": 0.30, "L2": 0.10, "L3": 0.30, "L4": 0.30},
        "tier_thresholds": {"tier1": 0.50, "tier2": 0.35, "tier3": 0.20},
        "rankings": [
            {
                "rank": i + 1,
                "binder": r["binder"],
                "length": r["length"],
                "integrated_score": r["integrated"],
                "screen_iptm": r["screen_iptm"],
                "binder_HAT_iptm": r["binder_HAT_iptm"],
                "binder_catcore_iptm": r["binder_catcore_iptm"],
                "layer_scores": r["scores"],
            }
            for i, r in enumerate(results)
        ],
        "summary": {
            "tier1_count": len(tier1),
            "tier2_count": len(tier2),
            "tier3_count": len(tier3),
            "excluded_count": len(excluded),
            "expanded_pool_available": 16,
        },
    }

    # Convert numpy types
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    out_path.write_text(json.dumps(convert(save_data), indent=2) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
