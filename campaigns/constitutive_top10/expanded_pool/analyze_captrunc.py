#!/usr/bin/env python3
"""
Analyze cap-truncated fusion results for expanded pool.

For each binder, checks whether the BINDER DOMAIN specifically contacts
hotspot residues when attached to the LOV2 scaffold (without UBQ cap).

Tests against both HAT and CatCore targets.
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
CT_DIR = BASE / "captrunc_results"

HAT_HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}
CATCORE_HOTSPOTS = {h + 241 for h in HAT_HOTSPOTS}
HAT_EQUIV = {h + 241: h for h in HAT_HOTSPOTS}

BINDER_LENGTHS = {
    "pep16s4": 30, "pep29s4": 50, "pep11s1": 20, "pep37s3": 20,
    "pep9s5": 20, "pep41s5": 20, "pep19s2": 20,
}


def parse_cif_atoms(cif_path):
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
                atoms.append((chain, resnum, x, y, z))
            except (ValueError, IndexError):
                continue
    return atoms


def domain_contacts(atoms, bdr_len, hotspots, cutoff=4.0):
    """Decompose contacts into binder vs scaffold contributions."""
    target_atoms = [(r, x, y, z) for c, r, x, y, z in atoms if c == "B"]
    if not target_atoms:
        return 0, set(), 0, set()

    target_coords = np.array([[x, y, z] for _, x, y, z in target_atoms])

    bdr_contacts = set()
    scaf_contacts = set()

    for c, r, x, y, z in atoms:
        if c != "A":
            continue
        dists = np.sqrt(np.sum((target_coords - np.array([x, y, z])) ** 2, axis=1))
        if np.min(dists) < cutoff:
            idx = np.argmin(dists)
            tres = target_atoms[idx][0]
            if r <= bdr_len:
                bdr_contacts.add(tres)
            else:
                scaf_contacts.add(tres)

    bdr_hs = bdr_contacts & hotspots
    scaf_hs = scaf_contacts & hotspots
    return len(bdr_contacts), bdr_hs, len(scaf_contacts), scaf_hs


def main():
    print("=" * 140)
    print("EXPANDED POOL — CAP-TRUNCATED FUSION DOMAIN DECOMPOSITION")
    print("Construct: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2] (no cap) + target")
    print("=" * 140)

    all_results = []

    for target_name, hotspots, target_label in [
        ("HAT", HAT_HOTSPOTS, "HAT"),
        ("catcore", CATCORE_HOTSPOTS, "CatCore"),
    ]:
        print(f"\n{'─'*60}")
        print(f"  TARGET: {target_label}")
        print(f"{'─'*60}")
        print(f"{'Binder':<10} {'BdrLen':<7} {'Bdr→tgt':<9} {'Bdr→HS':<9} "
              f"{'Scaf→tgt':<10} {'Scaf→HS':<9} {'Driver':<16} {'Best model bdr HS'}")
        print("-" * 120)

        for binder in sorted(BINDER_LENGTHS.keys()):
            bdr_len = BINDER_LENGTHS[binder]
            d = CT_DIR / f"{binder}_captrunc_{target_name}"
            if not d.exists():
                continue

            cifs = sorted(d.rglob("*model_*.cif"))
            if not cifs:
                continue

            bdr_hs_counts = []
            scaf_hs_counts = []
            best_bdr_hs = set()

            for cif in cifs:
                atoms = parse_cif_atoms(cif)
                bdr_n, bdr_hs, scaf_n, scaf_hs = domain_contacts(
                    atoms, bdr_len, hotspots
                )
                bdr_hs_counts.append(len(bdr_hs))
                scaf_hs_counts.append(len(scaf_hs))
                if len(bdr_hs) > len(best_bdr_hs):
                    best_bdr_hs = bdr_hs

            mean_bdr_hs = np.mean(bdr_hs_counts)
            mean_scaf_hs = np.mean(scaf_hs_counts)
            max_bdr_hs = max(bdr_hs_counts)

            if mean_bdr_hs > mean_scaf_hs and max_bdr_hs >= 3:
                driver = "BINDER-DRIVEN"
            elif mean_bdr_hs > 0 and mean_bdr_hs >= mean_scaf_hs * 0.5:
                driver = "MIXED"
            elif mean_bdr_hs > 0:
                driver = "SCAFFOLD-DOM"
            else:
                driver = "SCAFFOLD-ONLY"

            # Convert catcore hotspots to HAT equivalents for display
            if target_label == "CatCore":
                hs_display = sorted([HAT_EQUIV.get(h, h) for h in best_bdr_hs])
            else:
                hs_display = sorted(best_bdr_hs)

            hs_str = ",".join(str(h) for h in hs_display[:8])

            print(f"{binder:<10} {bdr_len:<7} {np.mean([domain_contacts(parse_cif_atoms(c), bdr_len, hotspots)[0] for c in cifs]):<9.1f} "
                  f"{mean_bdr_hs:<9.1f} "
                  f"{np.mean([domain_contacts(parse_cif_atoms(c), bdr_len, hotspots)[2] for c in cifs]):<10.1f} "
                  f"{mean_scaf_hs:<9.1f} {driver:<16} [{hs_str}]")

            all_results.append({
                "binder": binder,
                "target": target_label,
                "binder_length": bdr_len,
                "mean_bdr_hotspots": float(mean_bdr_hs),
                "max_bdr_hotspots": max_bdr_hs,
                "mean_scaf_hotspots": float(mean_scaf_hs),
                "driver": driver,
                "best_bdr_hotspots_hit": hs_display,
            })

    # Summary
    print()
    print("=" * 140)
    print("SUMMARY — BINDER DOMAIN REACHING TARGET THROUGH SCAFFOLD")
    print("=" * 140)

    # Group by binder
    by_binder = {}
    for r in all_results:
        by_binder.setdefault(r["binder"], {})[r["target"]] = r

    print(f"\n{'Binder':<10} {'HAT driver':<16} {'HAT bdr HS':<12} "
          f"{'CC driver':<16} {'CC bdr HS':<12} {'Overall'}")
    print("-" * 100)

    for binder in sorted(by_binder.keys()):
        hat = by_binder[binder].get("HAT", {})
        cc = by_binder[binder].get("CatCore", {})

        hat_driver = hat.get("driver", "?")
        cc_driver = cc.get("driver", "?")
        hat_hs = hat.get("mean_bdr_hotspots", 0)
        cc_hs = cc.get("mean_bdr_hotspots", 0)

        if "BINDER" in hat_driver or "BINDER" in cc_driver:
            overall = "STRONG"
        elif "MIXED" in hat_driver or "MIXED" in cc_driver:
            overall = "MODERATE"
        elif hat_hs > 0 or cc_hs > 0:
            overall = "WEAK"
        else:
            overall = "POOR"

        print(f"{binder:<10} {hat_driver:<16} {hat_hs:<12.1f} "
              f"{cc_driver:<16} {cc_hs:<12.1f} {overall}")

    # Save
    out_path = BASE / "expanded_captrunc_results.json"
    out_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
