#!/usr/bin/env python3
"""
Check where each binder ACTUALLY contacts HAT — hotspot overlap analysis.

For each binder-alone + HAT prediction, compute:
  1. Which HAT residues are within 4A of the binder (contact residues)
  2. How many of the 14 hotspot residues are contacted
  3. Which face of HAT the binder sits on

This catches off-target binders that have high ipTM but bind the wrong face.

Hotspot residues (1-indexed in HAT sequence):
  91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221
These are the substrate-binding groove residues — the functional target.
"""

import json
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
BV_DIR = BASE / "binding_validation" / "boltz2_results"

# HAT hotspot residues (1-indexed in the HAT sequence = chain B)
HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}

CONTACT_CUTOFF = 4.0  # Angstroms

# All binders tested
BINDERS = [
    "pep4s2", "pep13s3", "pep30s5", "pep29s5", "pep25s1",
    "pep21s5", "pep21s2", "pep15s3", "pep17s4", "pep10s5",
    "pep24s3", "pep32s1", "pep32s2", "pep16s5", "pep48s5",
    "pep36s1",
    "nb62s3", "nb134s2", "nb178s5", "nb87s2", "nb39s3", "nb153s2",
]


def parse_cif_atoms(cif_path):
    """Parse atom coordinates from a Boltz-2 CIF file.

    Boltz-2 CIF column layout:
      0=ATOM 1=id 2=element 3=atom_name 4=. 5=resname
      6=label_seq_id 7=auth_seq_id 8=? 9=label_asym_id
      10=x 11=y 12=z 13=occ 14=entity_id 15=auth_asym_id
      16=resname 17=bfactor 18=model

    Returns list of (chain, resnum, x, y, z, atom_name).
    """
    atoms = []
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 16:
                continue
            try:
                chain = parts[9]       # label_asym_id (A or B)
                resnum = int(parts[6]) # label_seq_id
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                atom_name = parts[3]
                atoms.append((chain, resnum, x, y, z, atom_name))
            except (ValueError, IndexError):
                continue
    return atoms


def find_contacts(atoms, chain_a="A", chain_b="B", cutoff=4.0):
    """Find chain B residues within cutoff of any chain A atom."""
    a_coords = []
    b_atoms = []

    for chain, resnum, x, y, z, aname in atoms:
        if chain == chain_a:
            a_coords.append([x, y, z])
        elif chain == chain_b:
            b_atoms.append((resnum, x, y, z))

    if not a_coords or not b_atoms:
        return set()

    a_arr = np.array(a_coords)  # (N, 3)
    contacted_residues = set()

    for resnum, bx, by, bz in b_atoms:
        dists = np.sqrt(np.sum((a_arr - np.array([bx, by, bz])) ** 2, axis=1))
        if np.min(dists) < cutoff:
            contacted_residues.add(resnum)

    return contacted_residues


def find_cif(binder_name, target="HAT"):
    """Find the binder-alone CIF file."""
    # Try different naming patterns
    patterns = [
        f"{binder_name}_binder_{target}",
        f"{binder_name}_50aa_binder_{target}",
        f"{binder_name}_20aa_binder_{target}",
        f"{binder_name}_30aa_binder_{target}",
        f"{binder_name}_40aa_binder_{target}",
        f"{binder_name}_120aa_binder_{target}",
    ]

    for pat in patterns:
        d = BV_DIR / pat
        if d.exists():
            cifs = list(d.rglob("*model_0.cif"))
            if cifs:
                return cifs[0]
    return None


def main():
    print("=" * 120)
    print("EPITOPE CONTACT ANALYSIS — WHERE DOES EACH BINDER ACTUALLY BIND HAT?")
    print(f"Contact cutoff: {CONTACT_CUTOFF}A  |  Hotspots: {sorted(HOTSPOTS)}")
    print("=" * 120)
    print()

    print(f"{'Binder':<12} {'HAT ipTM':<10} {'Contacts':<10} {'Hotspot':<10} "
          f"{'%Hotspot':<10} {'Verdict':<15} {'Contacted hotspots'}")
    print("-" * 120)

    # Load ipTM data
    bv_path = BASE / "binding_validation" / "analysis" / "binding_validation_results.json"
    with open(bv_path) as f:
        bv = json.load(f)
    iptm_map = {r["binder"]: r["binder_HAT_iptm"] for r in bv["raw_data"]}

    results = []

    for binder in BINDERS:
        cif = find_cif(binder, "HAT")
        iptm = iptm_map.get(binder, 0)

        if not cif:
            print(f"{binder:<12} {iptm:<10.4f} {'—':<10} {'—':<10} {'—':<10} NO CIF")
            continue

        atoms = parse_cif_atoms(cif)
        contacts = find_contacts(atoms, "A", "B", CONTACT_CUTOFF)
        hotspot_contacts = contacts & HOTSPOTS
        pct = len(hotspot_contacts) / len(HOTSPOTS) * 100 if HOTSPOTS else 0

        if pct >= 50:
            verdict = "ON-TARGET"
        elif pct >= 25:
            verdict = "PARTIAL"
        elif len(hotspot_contacts) > 0:
            verdict = "MARGINAL"
        else:
            verdict = "OFF-TARGET"

        hs_str = ",".join(str(h) for h in sorted(hotspot_contacts)) if hotspot_contacts else "none"

        print(f"{binder:<12} {iptm:<10.4f} {len(contacts):<10} "
              f"{len(hotspot_contacts):>2}/14     {pct:>5.1f}%    {verdict:<15} {hs_str}")

        results.append({
            "binder": binder,
            "binder_HAT_iptm": iptm,
            "total_contacts": len(contacts),
            "hotspot_contacts": len(hotspot_contacts),
            "hotspot_pct": pct,
            "verdict": verdict,
            "contacted_hotspots": sorted(hotspot_contacts),
            "all_contact_residues": sorted(contacts),
        })

    # Also check fusion constructs
    print()
    print("=" * 120)
    print("FUSION CONSTRUCT CONTACTS (lit state + HAT)")
    print("What's contacting HAT — is it the binder domain, or the scaffold/cap?")
    print("=" * 120)
    print()

    # Binder residue ranges in fusion (chain A)
    BINDER_RANGES = {
        "pep4s2": (1, 50),
        "pep13s3": (1, 20), "pep30s5": (1, 20), "pep29s5": (1, 30),
        "pep25s1": (1, 30), "pep21s5": (1, 20), "pep21s2": (1, 20),
        "pep15s3": (1, 20), "pep17s4": (1, 20), "pep10s5": (1, 20),
        "pep24s3": (1, 20), "pep32s1": (1, 20), "pep32s2": (1, 20),
        "pep16s5": (1, 30), "pep48s5": (1, 20), "pep36s1": (1, 40),
        "nb62s3": (1, 120), "nb134s2": (1, 120), "nb178s5": (1, 120),
        "nb87s2": (1, 120), "nb39s3": (1, 120), "nb153s2": (1, 120),
    }

    print(f"{'Binder':<12} {'Fus ipTM':<10} {'Bdr→HAT':<10} {'Scaffold→HAT':<14} "
          f"{'Cap→HAT':<10} {'Bdr hotspot':<12} {'Verdict'}")
    print("-" * 120)

    for binder in BINDERS:
        # Find fusion CIF
        fus_patterns = [
            f"{binder}_fusion_HAT",
            f"{binder}_50aa_fusion_HAT",
            f"{binder}_20aa_fusion_HAT",
            f"{binder}_30aa_fusion_HAT",
            f"{binder}_120aa_fusion_HAT",
        ]
        fus_cif = None
        for pat in fus_patterns:
            d = BV_DIR / pat
            if d.exists():
                cifs = list(d.rglob("*model_0.cif"))
                if cifs:
                    fus_cif = cifs[0]
                    break

        if not fus_cif:
            continue

        fus_iptm = {r["binder"]: r["fusion_HAT_iptm"] for r in bv["raw_data"]}.get(binder, 0)

        atoms = parse_cif_atoms(fus_cif)
        binder_start, binder_end = BINDER_RANGES[binder]

        # Separate binder atoms, scaffold atoms, and target atoms
        binder_coords = []
        scaffold_coords = []
        target_atoms = []

        # Determine cap range (last 76 residues of chain A)
        chain_a_max = max(r for c, r, x, y, z, a in atoms if c == "A")
        cap_start = chain_a_max - 75  # cap is 76 residues

        binder_atoms_list = []
        scaffold_atoms_list = []
        cap_atoms_list = []

        for chain, resnum, x, y, z, aname in atoms:
            if chain == "A":
                if binder_start <= resnum <= binder_end:
                    binder_coords.append([x, y, z])
                    binder_atoms_list.append((resnum, x, y, z))
                elif resnum >= cap_start:
                    cap_atoms_list.append((resnum, x, y, z))
                else:
                    scaffold_atoms_list.append((resnum, x, y, z))
            elif chain == "B":
                target_atoms.append((resnum, x, y, z))

        # Count contacts from each domain to HAT
        def count_domain_contacts(domain_atoms, target_atoms, cutoff):
            contacted = set()
            for _, dx, dy, dz in domain_atoms:
                for tres, tx, ty, tz in target_atoms:
                    dist = np.sqrt((dx-tx)**2 + (dy-ty)**2 + (dz-tz)**2)
                    if dist < cutoff:
                        contacted.add(tres)
            return contacted

        bdr_contacts = count_domain_contacts(binder_atoms_list, target_atoms, CONTACT_CUTOFF)
        scaf_contacts = count_domain_contacts(scaffold_atoms_list, target_atoms, CONTACT_CUTOFF)
        cap_contacts = count_domain_contacts(cap_atoms_list, target_atoms, CONTACT_CUTOFF)

        bdr_hs = bdr_contacts & HOTSPOTS

        if len(bdr_contacts) > len(scaf_contacts) + len(cap_contacts) and len(bdr_hs) >= 4:
            verdict = "BINDER-DRIVEN"
        elif len(cap_contacts) > len(bdr_contacts):
            verdict = "CAP-DOCKING"
        elif len(scaf_contacts) > len(bdr_contacts):
            verdict = "SCAFFOLD-DOCKING"
        elif len(bdr_contacts) == 0:
            verdict = "NO CONTACT"
        else:
            verdict = "MIXED"

        print(f"{binder:<12} {fus_iptm:<10.4f} {len(bdr_contacts):<10} "
              f"{len(scaf_contacts):<14} {len(cap_contacts):<10} "
              f"{len(bdr_hs):>2}/14       {verdict}")

    # Summary
    print()
    print("=" * 120)
    print("CRITICAL FINDING SUMMARY")
    print("=" * 120)

    on_target = [r for r in results if r["verdict"] == "ON-TARGET"]
    partial = [r for r in results if r["verdict"] == "PARTIAL"]
    marginal = [r for r in results if r["verdict"] == "MARGINAL"]
    off_target = [r for r in results if r["verdict"] == "OFF-TARGET"]

    print(f"\n  ON-TARGET  (>= 50% hotspots): {len(on_target)}")
    for r in on_target:
        print(f"    {r['binder']:<12} ipTM={r['binder_HAT_iptm']:.3f}  "
              f"hotspots={r['hotspot_contacts']}/14  contacts={r['total_contacts']}")

    print(f"\n  PARTIAL    (25-50% hotspots): {len(partial)}")
    for r in partial:
        print(f"    {r['binder']:<12} ipTM={r['binder_HAT_iptm']:.3f}  "
              f"hotspots={r['hotspot_contacts']}/14  contacts={r['total_contacts']}")

    print(f"\n  MARGINAL   (1-24% hotspots): {len(marginal)}")
    for r in marginal:
        print(f"    {r['binder']:<12} ipTM={r['binder_HAT_iptm']:.3f}  "
              f"hotspots={r['hotspot_contacts']}/14  contacts={r['total_contacts']}")

    print(f"\n  OFF-TARGET (0 hotspots):     {len(off_target)}")
    for r in off_target:
        print(f"    {r['binder']:<12} ipTM={r['binder_HAT_iptm']:.3f}  "
              f"contacts={r['total_contacts']} (all on wrong face)")

    # Save
    out_path = BASE / "epitope_contact_results.json"
    out_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
