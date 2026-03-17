#!/usr/bin/env python3
"""
CatCore Epitope Contact Analysis
=================================
For each binder+catcore Boltz-2 prediction, parse CIF models 0-4,
find chain B (catcore) residues within 4A of chain A (binder),
and check overlap with remapped HAT hotspot residues.

Compare results with the existing HAT epitope analysis.
"""

import os
import json
import math
from collections import defaultdict

# =============================================================================
# Constants
# =============================================================================

BASE_DIR = "/home/thinkingscopeanalysis/als_therapeutics/nanobody_optogenetic_protac/opto_scaffold_pipeline/campaigns/constitutive_top10"
RESULTS_DIR = os.path.join(BASE_DIR, "binding_validation/boltz2_results")

# HAT hotspots (1-indexed in HAT numbering)
HAT_HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}

# Catcore hotspots: catcore_1idx = 241 + hat_1idx (pre-gap, all hotspots are pre-gap)
CATCORE_HOTSPOTS = {241 + h for h in HAT_HOTSPOTS}
# = {332, 351, 352, 391, 399, 400, 401, 404, 411, 422, 426, 460, 461, 462}

CONTACT_CUTOFF = 4.0  # Angstroms

BINDERS = [
    "pep4s2", "pep13s3", "pep30s5", "pep29s5", "pep25s1", "pep21s5",
    "pep21s2", "pep15s3", "pep17s4", "pep10s5", "pep24s3", "pep32s1",
    "pep32s2", "pep16s5", "pep48s5", "pep36s1"
]

# =============================================================================
# CIF Parsing
# =============================================================================

def parse_cif_atoms(cif_path):
    """Parse ATOM lines from Boltz-2 CIF file.
    Returns dict: chain -> list of (resnum, x, y, z)
    """
    atoms = defaultdict(list)
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            chain = parts[9]       # label_asym_id
            resnum = int(parts[6]) # label_seq_id
            x = float(parts[10])
            y = float(parts[11])
            z = float(parts[12])
            atoms[chain].append((resnum, x, y, z))
    return atoms


def find_contact_residues(atoms, chain_binder="A", chain_target="B", cutoff=4.0):
    """Find target residues with any atom within cutoff of any binder atom.
    Returns set of target residue numbers.
    """
    binder_atoms = atoms.get(chain_binder, [])
    target_atoms = atoms.get(chain_target, [])

    if not binder_atoms or not target_atoms:
        return set()

    cutoff_sq = cutoff * cutoff
    contact_residues = set()

    # Group target atoms by residue for efficiency
    target_by_res = defaultdict(list)
    for resnum, x, y, z in target_atoms:
        target_by_res[resnum].append((x, y, z))

    # For each binder atom, check against target atoms
    for _, bx, by, bz in binder_atoms:
        for resnum, coords_list in target_by_res.items():
            if resnum in contact_residues:
                continue  # already found
            for tx, ty, tz in coords_list:
                dx = bx - tx
                dy = by - ty
                dz = bz - tz
                if dx*dx + dy*dy + dz*dz < cutoff_sq:
                    contact_residues.add(resnum)
                    break

    return contact_residues


def get_iptm(binder_name, model_idx):
    """Get ipTM from confidence JSON."""
    pred_dir = os.path.join(RESULTS_DIR, f"{binder_name}_binder_catcore",
                            f"boltz_results_{binder_name}_binder_catcore",
                            "predictions", f"{binder_name}_binder_catcore")
    conf_path = os.path.join(pred_dir, f"confidence_{binder_name}_binder_catcore_model_{model_idx}.json")
    if not os.path.exists(conf_path):
        return None
    with open(conf_path) as f:
        data = json.load(f)
    return data.get("iptm", data.get("protein_iptm"))


def get_cif_path(binder_name, model_idx):
    """Get CIF path for a binder+catcore prediction."""
    pred_dir = os.path.join(RESULTS_DIR, f"{binder_name}_binder_catcore",
                            f"boltz_results_{binder_name}_binder_catcore",
                            "predictions", f"{binder_name}_binder_catcore")
    return os.path.join(pred_dir, f"{binder_name}_binder_catcore_model_{model_idx}.cif")


def classify_verdict(n_hotspot_contacts, total_hotspots=14):
    """Classify epitope targeting verdict."""
    pct = n_hotspot_contacts / total_hotspots * 100
    if pct >= 50:
        return "ON-TARGET"
    elif pct >= 25:
        return "PARTIAL"
    elif pct > 0:
        return "MARGINAL"
    else:
        return "OFF-TARGET"


# =============================================================================
# Main Analysis
# =============================================================================

def main():
    # Load existing HAT epitope results
    hat_results_path = os.path.join(BASE_DIR, "epitope_contact_results.json")
    with open(hat_results_path) as f:
        hat_results_list = json.load(f)
    hat_results = {r["binder"]: r for r in hat_results_list}

    print("=" * 120)
    print("CATCORE EPITOPE CONTACT ANALYSIS")
    print("=" * 120)
    print()
    print(f"Contact cutoff: {CONTACT_CUTOFF} A")
    print(f"HAT hotspots (1-idx): {sorted(HAT_HOTSPOTS)}")
    print(f"CatCore hotspots (1-idx): {sorted(CATCORE_HOTSPOTS)}")
    print(f"Mapping: catcore_res = hat_res + 241 (pre-gap formula)")
    print()

    # Analyze each binder across all 5 models
    catcore_results = []

    for binder in BINDERS:
        print("-" * 120)
        print(f"BINDER: {binder}")
        print("-" * 120)

        model_data = []
        all_hotspot_contacts_union = set()
        all_contacts_union = set()

        for model_idx in range(5):
            cif_path = get_cif_path(binder, model_idx)
            if not os.path.exists(cif_path):
                print(f"  Model {model_idx}: CIF not found")
                continue

            atoms = parse_cif_atoms(cif_path)
            chains = sorted(atoms.keys())
            chain_sizes = {c: len(set(r for r, _, _, _ in atoms[c])) for c in chains}

            contacts = find_contact_residues(atoms, "A", "B", CONTACT_CUTOFF)
            hotspot_contacts = contacts & CATCORE_HOTSPOTS

            # Convert catcore hotspot contacts back to HAT numbering for display
            hat_equiv = sorted([c - 241 for c in hotspot_contacts])

            iptm = get_iptm(binder, model_idx)

            model_data.append({
                "model": model_idx,
                "iptm": iptm,
                "total_contacts": len(contacts),
                "hotspot_contacts": len(hotspot_contacts),
                "hotspot_set": hotspot_contacts,
                "hat_equiv": hat_equiv,
                "all_contacts": contacts,
            })

            all_hotspot_contacts_union |= hotspot_contacts
            all_contacts_union |= contacts

            verdict = classify_verdict(len(hotspot_contacts))
            print(f"  Model {model_idx}: ipTM={iptm:.4f}  contacts={len(contacts):3d}  "
                  f"hotspot={len(hotspot_contacts):2d}/14  "
                  f"HAT-equiv={hat_equiv}  [{verdict}]")

        # Summary across models
        hotspot_counts = [m["hotspot_contacts"] for m in model_data]
        avg_hotspot = sum(hotspot_counts) / len(hotspot_counts) if hotspot_counts else 0
        min_hotspot = min(hotspot_counts) if hotspot_counts else 0
        max_hotspot = max(hotspot_counts) if hotspot_counts else 0

        # Best model by hotspot contacts
        best_model = max(model_data, key=lambda m: (m["hotspot_contacts"], m["iptm"] or 0)) if model_data else None

        # Consistency: hotspots seen in at least 3/5 models
        from collections import Counter
        hotspot_freq = Counter()
        for m in model_data:
            for h in m["hotspot_set"]:
                hotspot_freq[h] += 1
        consistent_hotspots = {h for h, c in hotspot_freq.items() if c >= 3}

        union_hat_equiv = sorted([c - 241 for c in all_hotspot_contacts_union])
        consistent_hat_equiv = sorted([c - 241 for c in consistent_hotspots])

        print(f"\n  SUMMARY: avg={avg_hotspot:.1f}  min={min_hotspot}  max={max_hotspot}  "
              f"union={len(all_hotspot_contacts_union)}/14  consistent(>=3/5)={len(consistent_hotspots)}/14")
        print(f"  Union hotspots (HAT equiv): {union_hat_equiv}")
        print(f"  Consistent hotspots (HAT equiv): {consistent_hat_equiv}")

        # Store result (using best model for primary comparison)
        catcore_results.append({
            "binder": binder,
            "best_model_idx": best_model["model"] if best_model else None,
            "best_iptm": best_model["iptm"] if best_model else None,
            "best_total_contacts": best_model["total_contacts"] if best_model else 0,
            "best_hotspot_contacts": best_model["hotspot_contacts"] if best_model else 0,
            "best_hotspot_pct": (best_model["hotspot_contacts"] / 14 * 100) if best_model else 0,
            "best_verdict": classify_verdict(best_model["hotspot_contacts"]) if best_model else "N/A",
            "best_hat_equiv": sorted([c - 241 for c in best_model["hotspot_set"]]) if best_model else [],
            "avg_hotspot": avg_hotspot,
            "min_hotspot": min_hotspot,
            "max_hotspot": max_hotspot,
            "union_hotspots_hat": union_hat_equiv,
            "consistent_hotspots_hat": consistent_hat_equiv,
            "n_consistent": len(consistent_hotspots),
            "all_model_data": [
                {
                    "model": m["model"],
                    "iptm": m["iptm"],
                    "total_contacts": m["total_contacts"],
                    "hotspot_contacts": m["hotspot_contacts"],
                    "hat_equiv": m["hat_equiv"],
                }
                for m in model_data
            ],
        })
        print()

    # =========================================================================
    # Comparison Table: HAT vs CatCore
    # =========================================================================
    print()
    print("=" * 140)
    print("COMPARISON TABLE: HAT vs CatCore Epitope Contacts (best model per binder)")
    print("=" * 140)
    print()

    header = (f"{'Binder':<10} | {'HAT ipTM':>8} | {'HAT Hot':>7} | {'HAT Vrd':<10} | "
              f"{'CC ipTM':>8} | {'CC Hot':>6} | {'CC Vrd':<10} | "
              f"{'CC Avg':>6} | {'CC Min':>6} | {'CC Max':>6} | {'CC Con':>6} | "
              f"{'Match?':<8} | HAT Hotspots Contacted (HAT numbering)")
    print(header)
    print("-" * 140)

    both_on_target = []
    hat_only = []
    catcore_only = []
    neither = []

    for cr in catcore_results:
        binder = cr["binder"]
        hr = hat_results.get(binder, {})

        hat_iptm = hr.get("binder_HAT_iptm", 0)
        hat_hot = hr.get("hotspot_contacts", 0)
        hat_verdict = hr.get("verdict", "N/A")
        hat_hotspots_contacted = sorted(hr.get("contacted_hotspots", []))

        cc_iptm = cr["best_iptm"] or 0
        cc_hot = cr["best_hotspot_contacts"]
        cc_verdict = cr["best_verdict"]
        cc_avg = cr["avg_hotspot"]
        cc_min = cr["min_hotspot"]
        cc_max = cr["max_hotspot"]
        cc_con = cr["n_consistent"]
        cc_hat_equiv = cr["best_hat_equiv"]

        # Determine if "on-target" in each context (>= PARTIAL = >= 25%)
        hat_on = hat_verdict in ("ON-TARGET", "PARTIAL")
        cc_on = cc_verdict in ("ON-TARGET", "PARTIAL")

        if hat_on and cc_on:
            match_str = "BOTH"
            both_on_target.append(binder)
        elif hat_on and not cc_on:
            match_str = "HAT-only"
            hat_only.append(binder)
        elif not hat_on and cc_on:
            match_str = "CC-only"
            catcore_only.append(binder)
        else:
            match_str = "NEITHER"
            neither.append(binder)

        # Show HAT hotspots contacted in both HAT and CatCore
        hat_set = set(hat_hotspots_contacted)
        cc_set = set(cc_hat_equiv)
        overlap = sorted(hat_set & cc_set)
        hat_unique = sorted(hat_set - cc_set)
        cc_unique = sorted(cc_set - hat_set)

        hotspot_str = ""
        if overlap:
            hotspot_str += f"BOTH:{overlap} "
        if hat_unique:
            hotspot_str += f"HAT-only:{hat_unique} "
        if cc_unique:
            hotspot_str += f"CC-only:{cc_unique}"
        if not hotspot_str:
            hotspot_str = "(none)"

        print(f"{binder:<10} | {hat_iptm:>8.4f} | {hat_hot:>4d}/14 | {hat_verdict:<10} | "
              f"{cc_iptm:>8.4f} | {cc_hot:>3d}/14 | {cc_verdict:<10} | "
              f"{cc_avg:>6.1f} | {cc_min:>6d} | {cc_max:>6d} | {cc_con:>3d}/14 | "
              f"{match_str:<8} | {hotspot_str}")

    print("-" * 140)
    print()

    # =========================================================================
    # Summary Statistics
    # =========================================================================
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"ON-TARGET in BOTH HAT and CatCore:  {len(both_on_target):2d}/16  {both_on_target}")
    print(f"ON-TARGET in HAT only:              {len(hat_only):2d}/16  {hat_only}")
    print(f"ON-TARGET in CatCore only:          {len(catcore_only):2d}/16  {catcore_only}")
    print(f"OFF-TARGET in BOTH:                 {len(neither):2d}/16  {neither}")
    print()

    # Consistency analysis
    print("=" * 80)
    print("MULTI-MODEL CONSISTENCY (CatCore, 5 diffusion samples)")
    print("=" * 80)
    print()
    print(f"{'Binder':<10} | {'M0':>3} {'M1':>3} {'M2':>3} {'M3':>3} {'M4':>3} | {'Avg':>5} | {'Consistent (>=3/5)':>20} | Verdict")
    print("-" * 80)

    for cr in catcore_results:
        binder = cr["binder"]
        counts = [m["hotspot_contacts"] for m in cr["all_model_data"]]
        while len(counts) < 5:
            counts.append(0)
        avg = cr["avg_hotspot"]
        con = cr["consistent_hotspots_hat"]

        # Verdict based on consistency
        if cr["n_consistent"] >= 7:
            con_verdict = "CONSISTENT ON-TARGET"
        elif cr["n_consistent"] >= 4:
            con_verdict = "CONSISTENT PARTIAL"
        elif cr["n_consistent"] >= 1:
            con_verdict = "VARIABLE"
        else:
            con_verdict = "CONSISTENTLY OFF"

        print(f"{binder:<10} | {counts[0]:>3} {counts[1]:>3} {counts[2]:>3} {counts[3]:>3} {counts[4]:>3} | "
              f"{avg:>5.1f} | {str(con):>20} | {con_verdict}")

    print("-" * 80)
    print()

    # =========================================================================
    # Per-hotspot residue analysis
    # =========================================================================
    print("=" * 100)
    print("PER-HOTSPOT RESIDUE FREQUENCY ACROSS ALL BINDERS (CatCore, best model)")
    print("=" * 100)
    print()

    hotspot_binder_map = defaultdict(list)
    for cr in catcore_results:
        for h in cr["best_hat_equiv"]:
            hotspot_binder_map[h].append(cr["binder"])

    print(f"{'HAT Res':<8} {'CC Res':<8} {'#Binders':>8}  Binders")
    print("-" * 100)
    for h in sorted(HAT_HOTSPOTS):
        cc_res = h + 241
        binders_hit = hotspot_binder_map.get(h, [])
        print(f"{h:<8} {cc_res:<8} {len(binders_hit):>8}  {binders_hit}")

    print()

    # =========================================================================
    # Final ranking: binders ON-TARGET in both HAT and CatCore
    # =========================================================================
    print("=" * 100)
    print("FINAL RANKING: Binders ON-TARGET (>=PARTIAL) in BOTH HAT and CatCore")
    print("=" * 100)
    print()

    if both_on_target:
        ranked = []
        for binder in both_on_target:
            hr = hat_results[binder]
            cr = [c for c in catcore_results if c["binder"] == binder][0]
            # Score: average of HAT and CatCore hotspot counts
            combined = (hr["hotspot_contacts"] + cr["best_hotspot_contacts"]) / 2
            ranked.append((binder, hr, cr, combined))

        ranked.sort(key=lambda x: (-x[3], -(x[1].get("binder_HAT_iptm", 0))))

        print(f"{'Rank':<5} {'Binder':<10} {'HAT ipTM':>8} {'HAT Hot':>7} {'HAT Vrd':<10} "
              f"{'CC ipTM':>8} {'CC Hot':>6} {'CC Vrd':<10} {'CC Consist':>10} {'Combined':>8}")
        print("-" * 100)
        for i, (binder, hr, cr, combined) in enumerate(ranked, 1):
            print(f"{i:<5} {binder:<10} {hr['binder_HAT_iptm']:>8.4f} {hr['hotspot_contacts']:>4d}/14 "
                  f"{hr['verdict']:<10} {cr['best_iptm']:>8.4f} {cr['best_hotspot_contacts']:>3d}/14 "
                  f"{cr['best_verdict']:<10} {cr['n_consistent']:>7d}/14 {combined:>8.1f}")
        print("-" * 100)
    else:
        print("  NO binders are ON-TARGET in both HAT and CatCore predictions!")

    print()

    # Save results
    output_path = os.path.join(BASE_DIR, "catcore_epitope_results.json")
    with open(output_path, "w") as f:
        json.dump(catcore_results, f, indent=2, default=str)
    print(f"Results saved to: {output_path}")

    # Also save comparison summary
    comparison = {
        "both_on_target": both_on_target,
        "hat_only": hat_only,
        "catcore_only": catcore_only,
        "neither": neither,
        "catcore_hotspots_1idx": sorted(CATCORE_HOTSPOTS),
        "hat_hotspots_1idx": sorted(HAT_HOTSPOTS),
        "mapping": "catcore_1idx = hat_1idx + 241 (pre-gap)",
    }
    comp_path = os.path.join(BASE_DIR, "hat_vs_catcore_comparison.json")
    with open(comp_path, "w") as f:
        json.dump(comparison, f, indent=2)
    print(f"Comparison saved to: {comp_path}")


if __name__ == "__main__":
    main()
