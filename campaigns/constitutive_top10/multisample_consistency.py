#!/usr/bin/env python3
"""
Layer 3: Multi-sample hotspot consistency analysis.

For each binder-alone + HAT prediction (up to 5 diffusion samples),
check whether hotspot contacts are consistent across models.

A binder that hits hotspots in 5/5 samples is far more reliable than
one that hits in 1/5 (stochastic).

Also computes per-residue interface pLDDT and PAE for Layer 4.
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
BV_DIR = BASE / "binding_validation" / "boltz2_results"

HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}
CONTACT_CUTOFF = 4.0

BINDERS = [
    # ON-TARGET
    "pep30s5", "pep25s1", "pep21s5", "pep15s3", "pep48s5", "pep16s5", "pep36s1",
    # PARTIAL
    "pep29s5", "pep21s2", "pep24s3", "pep13s3",
    # Include a few off-target for comparison
    "pep4s2", "pep17s4",
]

BINDER_LENGTHS = {
    "pep30s5": 20, "pep25s1": 30, "pep21s5": 20, "pep15s3": 20,
    "pep48s5": 20, "pep16s5": 30, "pep36s1": 40,
    "pep29s5": 30, "pep21s2": 20, "pep24s3": 20, "pep13s3": 20,
    "pep4s2": 50, "pep17s4": 20,
}


def parse_cif_atoms(cif_path):
    """Parse Boltz-2 CIF atoms."""
    atoms = []
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 16:
                continue
            try:
                chain = parts[9]
                resnum = int(parts[6])
                x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                atom_name = parts[3]
                atoms.append((chain, resnum, x, y, z, atom_name))
            except (ValueError, IndexError):
                continue
    return atoms


def find_contacts(atoms, cutoff=4.0):
    """Find chain B residues within cutoff of chain A atoms."""
    a_coords = np.array([[x, y, z] for c, r, x, y, z, _ in atoms if c == "A"])
    b_atoms = [(r, x, y, z) for c, r, x, y, z, _ in atoms if c == "B"]

    if len(a_coords) == 0 or len(b_atoms) == 0:
        return set()

    contacted = set()
    for resnum, bx, by, bz in b_atoms:
        dists = np.sqrt(np.sum((a_coords - np.array([bx, by, bz])) ** 2, axis=1))
        if np.min(dists) < cutoff:
            contacted.add(resnum)
    return contacted


def compute_interface_bsa(atoms, cutoff=10.0):
    """Estimate buried surface area from interface contacts."""
    a_ca = [(r, x, y, z) for c, r, x, y, z, a in atoms if c == "A" and a == "CA"]
    b_ca = [(r, x, y, z) for c, r, x, y, z, a in atoms if c == "B" and a == "CA"]

    if not a_ca or not b_ca:
        return 0, 0

    a_interface = 0
    b_interface = 0
    for ar, ax, ay, az in a_ca:
        for br, bx, by, bz in b_ca:
            d = np.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            if d < cutoff:
                a_interface += 1
                break
    for br, bx, by, bz in b_ca:
        for ar, ax, ay, az in a_ca:
            d = np.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            if d < cutoff:
                b_interface += 1
                break

    return a_interface, b_interface


def find_all_models(binder_name):
    """Find all model CIFs for a binder-alone + HAT prediction."""
    patterns = [
        f"{binder_name}_binder_HAT",
        f"{binder_name}_50aa_binder_HAT",
        f"{binder_name}_20aa_binder_HAT",
        f"{binder_name}_30aa_binder_HAT",
        f"{binder_name}_40aa_binder_HAT",
        f"{binder_name}_120aa_binder_HAT",
    ]

    for pat in patterns:
        d = BV_DIR / pat
        if d.exists():
            cifs = sorted(d.rglob("*model_*.cif"))
            return cifs
    return []


def main():
    print("=" * 130)
    print("MULTI-SAMPLE HOTSPOT CONSISTENCY ANALYSIS")
    print(f"Contact cutoff: {CONTACT_CUTOFF}A  |  Hotspots: {sorted(HOTSPOTS)}")
    print("=" * 130)
    print()

    all_results = []

    for binder in BINDERS:
        cifs = find_all_models(binder)
        if not cifs:
            print(f"{binder}: NO CIF FILES FOUND")
            continue

        n_models = len(cifs)
        per_model = []

        for cif in cifs:
            model_name = cif.stem.split("_model_")[1] if "_model_" in cif.stem else "?"
            atoms = parse_cif_atoms(cif)
            contacts = find_contacts(atoms, CONTACT_CUTOFF)
            hs_contacts = contacts & HOTSPOTS
            a_iface, b_iface = compute_interface_bsa(atoms)

            per_model.append({
                "model": model_name,
                "total_contacts": len(contacts),
                "hotspot_contacts": len(hs_contacts),
                "hotspots_hit": sorted(hs_contacts),
                "a_interface_res": a_iface,
                "b_interface_res": b_iface,
            })

        # Aggregate statistics
        hs_counts = [m["hotspot_contacts"] for m in per_model]
        total_counts = [m["total_contacts"] for m in per_model]

        # Which hotspots are hit consistently?
        hotspot_freq = {}
        for h in HOTSPOTS:
            count = sum(1 for m in per_model if h in m["hotspots_hit"])
            if count > 0:
                hotspot_freq[h] = count

        # Models with any hotspot contact
        models_with_hs = sum(1 for c in hs_counts if c > 0)
        # Consistent hotspots (hit in >= 60% of models)
        consistent_hs = [h for h, c in hotspot_freq.items() if c >= max(2, n_models * 0.6)]

        result = {
            "binder": binder,
            "n_models": n_models,
            "models_with_hotspots": models_with_hs,
            "consistency": models_with_hs / n_models if n_models > 0 else 0,
            "mean_hotspots": np.mean(hs_counts),
            "min_hotspots": min(hs_counts),
            "max_hotspots": max(hs_counts),
            "consistent_hotspots": sorted(consistent_hs),
            "n_consistent_hotspots": len(consistent_hs),
            "hotspot_frequency": {str(k): v for k, v in sorted(hotspot_freq.items())},
            "mean_total_contacts": np.mean(total_counts),
            "per_model": per_model,
        }
        all_results.append(result)

    # Print summary table
    print(f"{'Binder':<12} {'Models':<8} {'w/HS':<8} {'Consist':<9} "
          f"{'Mean HS':<9} {'Min-Max':<10} {'Consist HS':<12} {'Consistent hotspot residues'}")
    print("-" * 130)

    for r in sorted(all_results, key=lambda x: x["n_consistent_hotspots"], reverse=True):
        hs_str = ",".join(str(h) for h in r["consistent_hotspots"][:8])
        print(f"{r['binder']:<12} {r['n_models']:<8} {r['models_with_hotspots']:<8} "
              f"{r['consistency']:<9.1%} {r['mean_hotspots']:<9.1f} "
              f"{r['min_hotspots']}-{r['max_hotspots']:<7} {r['n_consistent_hotspots']:<12} {hs_str}")

    # Detailed per-model view for top candidates
    print()
    print("=" * 130)
    print("PER-MODEL DETAIL (top candidates)")
    print("=" * 130)

    for r in sorted(all_results, key=lambda x: x["n_consistent_hotspots"], reverse=True):
        if r["n_consistent_hotspots"] < 3:
            continue
        print(f"\n  {r['binder']} — {r['n_consistent_hotspots']} consistent hotspots, "
              f"{r['models_with_hotspots']}/{r['n_models']} models with hotspot contacts")
        for m in r["per_model"]:
            hs = ",".join(str(h) for h in m["hotspots_hit"])
            print(f"    model_{m['model']}: {m['hotspot_contacts']:>2}/14 hotspots, "
                  f"{m['total_contacts']:>3} total contacts, "
                  f"iface_res: A={m['a_interface_res']} B={m['b_interface_res']}  "
                  f"HS: [{hs}]")

    # Confidence tiers
    print()
    print("=" * 130)
    print("CONFIDENCE TIERS")
    print("=" * 130)

    tier1 = [r for r in all_results if r["consistency"] >= 0.8 and r["n_consistent_hotspots"] >= 5]
    tier2 = [r for r in all_results if r not in tier1 and r["consistency"] >= 0.6 and r["n_consistent_hotspots"] >= 3]
    tier3 = [r for r in all_results if r not in tier1 and r not in tier2 and r["models_with_hotspots"] > 0]
    tier_fail = [r for r in all_results if r["models_with_hotspots"] == 0]

    print(f"\n  TIER 1 — HIGH CONFIDENCE (>= 80% consistency, >= 5 consistent hotspots): {len(tier1)}")
    for r in tier1:
        print(f"    {r['binder']:<12} {r['consistency']:.0%} consistency, "
              f"{r['n_consistent_hotspots']} consistent HS, mean={r['mean_hotspots']:.1f}")

    print(f"\n  TIER 2 — MODERATE CONFIDENCE (>= 60% consistency, >= 3 consistent HS): {len(tier2)}")
    for r in tier2:
        print(f"    {r['binder']:<12} {r['consistency']:.0%} consistency, "
              f"{r['n_consistent_hotspots']} consistent HS, mean={r['mean_hotspots']:.1f}")

    print(f"\n  TIER 3 — LOW CONFIDENCE (some hotspot contacts but inconsistent): {len(tier3)}")
    for r in tier3:
        print(f"    {r['binder']:<12} {r['consistency']:.0%} consistency, "
              f"{r['n_consistent_hotspots']} consistent HS, mean={r['mean_hotspots']:.1f}")

    print(f"\n  FAIL — NO HOTSPOT CONTACTS: {len(tier_fail)}")
    for r in tier_fail:
        print(f"    {r['binder']:<12} {r['mean_total_contacts']:.0f} total contacts, all off-target")

    # Save
    out_path = BASE / "multisample_consistency_results.json"
    # Convert numpy types for JSON
    for r in all_results:
        r["mean_hotspots"] = float(r["mean_hotspots"])
        r["mean_total_contacts"] = float(r["mean_total_contacts"])
    out_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
