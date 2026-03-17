#!/usr/bin/env python3
"""
Layer 4: Interface quality metrics from Boltz-2 confidence outputs.

For each binder-alone + HAT prediction, extract:
  1. Interface pLDDT (binder residues at interface)
  2. Cross-chain PAE (binder→HAT and HAT→binder blocks)
  3. Overall confidence score and ipTM
  4. Per-model best metrics across diffusion samples

Higher interface pLDDT and lower cross-chain PAE = more confident interface.
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
BV_DIR = BASE / "binding_validation" / "boltz2_results"

HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}
CONTACT_CUTOFF = 4.0

BINDERS = {
    # ON-TARGET
    "pep30s5": 20, "pep25s1": 30, "pep21s5": 20, "pep15s3": 20,
    "pep48s5": 20, "pep16s5": 30, "pep36s1": 40,
    # PARTIAL
    "pep29s5": 30, "pep21s2": 20, "pep24s3": 20, "pep13s3": 20,
    # OFF-TARGET (for comparison)
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


def find_interface_residues(atoms, cutoff=8.0):
    """Find binder (chain A) residues at the interface with HAT (chain B)."""
    a_ca = [(r, x, y, z) for c, r, x, y, z, a in atoms if c == "A" and a == "CA"]
    b_ca = [(r, x, y, z) for c, r, x, y, z, a in atoms if c == "B" and a == "CA"]

    a_interface = set()
    b_interface = set()
    for ar, ax, ay, az in a_ca:
        for br, bx, by, bz in b_ca:
            d = np.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            if d < cutoff:
                a_interface.add(ar)
                b_interface.add(br)

    return sorted(a_interface), sorted(b_interface)


def find_all_models(binder_name):
    """Find all model CIFs and confidence JSONs."""
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
            return d, cifs
    return None, []


def extract_confidence(result_dir, model_idx):
    """Extract confidence metrics from Boltz-2 JSON output."""
    # Try different confidence file patterns
    patterns = [
        result_dir / "predictions" / f"*_confidence_{model_idx}.json",
    ]

    # Search recursively
    jsons = sorted(result_dir.rglob(f"*confidence*{model_idx}*"))
    if not jsons:
        jsons = sorted(result_dir.rglob(f"*confidence*"))

    if not jsons:
        return None

    # Try to find the right one for this model
    for j in jsons:
        if f"_{model_idx}" in j.stem or f"model_{model_idx}" in j.stem:
            with open(j) as f:
                return json.load(f)

    # Fallback: load first one
    if jsons:
        with open(jsons[0]) as f:
            return json.load(f)
    return None


def extract_pae(result_dir, model_idx):
    """Extract PAE matrix from Boltz-2 .npz output."""
    npzs = sorted(result_dir.rglob(f"*pae*{model_idx}*"))
    if not npzs:
        npzs = sorted(result_dir.rglob("*pae*"))

    for npz_path in npzs:
        if f"_{model_idx}" in npz_path.stem or f"model_{model_idx}" in npz_path.stem:
            data = np.load(npz_path)
            for key in data.files:
                return data[key]

    if npzs:
        data = np.load(npzs[0])
        for key in data.files:
            return data[key]
    return None


def main():
    print("=" * 130)
    print("INTERFACE QUALITY METRICS — BOLTZ-2 CONFIDENCE ANALYSIS")
    print(f"Contact cutoff: {CONTACT_CUTOFF}A")
    print("=" * 130)
    print()

    all_results = []

    for binder, binder_len in BINDERS.items():
        result_dir, cifs = find_all_models(binder)
        if not cifs:
            print(f"{binder}: NO CIF FILES FOUND")
            continue

        best_iptm = 0
        best_model = None

        per_model = []

        for i, cif in enumerate(cifs):
            model_idx = cif.stem.split("_model_")[1] if "_model_" in cif.stem else str(i)

            # Parse structure
            atoms = parse_cif_atoms(cif)
            a_interface, b_interface = find_interface_residues(atoms, cutoff=8.0)

            # Extract confidence
            conf = extract_confidence(result_dir, model_idx)
            iptm = 0
            ptm = 0
            confidence = 0
            if conf:
                iptm = conf.get("protein_iptm", conf.get("iptm", 0))
                ptm = conf.get("ptm", 0)
                confidence = conf.get("confidence_score", 0)

            # Extract PAE
            pae_matrix = extract_pae(result_dir, model_idx)
            cross_pae_ab = None
            cross_pae_ba = None
            if pae_matrix is not None and pae_matrix.ndim == 2:
                total_len = pae_matrix.shape[0]
                if total_len == binder_len + 317:
                    # Cross-chain PAE blocks
                    # A→B: rows 0:binder_len, cols binder_len:total
                    ab_block = pae_matrix[:binder_len, binder_len:]
                    ba_block = pae_matrix[binder_len:, :binder_len]
                    cross_pae_ab = float(np.mean(ab_block))
                    cross_pae_ba = float(np.mean(ba_block))

                    # Interface-specific PAE (only interface residues)
                    if a_interface and b_interface:
                        a_idx = [r - 1 for r in a_interface if r - 1 < binder_len]
                        b_idx = [r - 1 + binder_len for r in b_interface
                                 if r - 1 + binder_len < total_len]
                        if a_idx and b_idx:
                            iface_pae = pae_matrix[np.ix_(a_idx, b_idx)]
                            cross_pae_ab = float(np.mean(iface_pae))

            # Interface pLDDT from CIF B-factors
            binder_bfactors = []
            target_bfactors = []
            for chain, resnum, x, y, z, aname in atoms:
                if aname != "CA":
                    continue
                # B-factor column in Boltz-2 CIF
                # Actually need to re-parse with bfactor
            # Re-parse for bfactors
            bfactors_a = {}
            bfactors_b = {}
            with open(cif) as f:
                for line in f:
                    if not line.startswith("ATOM"):
                        continue
                    parts = line.split()
                    if len(parts) < 18:
                        continue
                    try:
                        chain = parts[9]
                        resnum = int(parts[6])
                        aname = parts[3]
                        bfactor = float(parts[17])
                        if aname == "CA":
                            if chain == "A":
                                bfactors_a[resnum] = bfactor
                            elif chain == "B":
                                bfactors_b[resnum] = bfactor
                    except (ValueError, IndexError):
                        continue

            # Interface pLDDT
            iface_plddt_a = [bfactors_a[r] for r in a_interface if r in bfactors_a]
            iface_plddt_b = [bfactors_b[r] for r in b_interface if r in bfactors_b]
            mean_iface_plddt = (
                np.mean(iface_plddt_a + iface_plddt_b)
                if iface_plddt_a or iface_plddt_b else 0
            )

            # Binder overall pLDDT
            binder_plddt = np.mean(list(bfactors_a.values())) if bfactors_a else 0

            model_data = {
                "model": model_idx,
                "iptm": iptm,
                "ptm": ptm,
                "confidence": confidence,
                "binder_plddt": float(binder_plddt),
                "interface_plddt": float(mean_iface_plddt),
                "n_interface_res_a": len(a_interface),
                "n_interface_res_b": len(b_interface),
                "cross_pae_ab": cross_pae_ab,
                "cross_pae_ba": cross_pae_ba,
            }
            per_model.append(model_data)

            if iptm > best_iptm:
                best_iptm = iptm
                best_model = model_data

        # Aggregate
        iptms = [m["iptm"] for m in per_model]
        iface_plddts = [m["interface_plddt"] for m in per_model if m["interface_plddt"] > 0]
        paes = [m["cross_pae_ab"] for m in per_model if m["cross_pae_ab"] is not None]

        result = {
            "binder": binder,
            "binder_length": binder_len,
            "n_models": len(per_model),
            "best_iptm": max(iptms) if iptms else 0,
            "mean_iptm": float(np.mean(iptms)) if iptms else 0,
            "std_iptm": float(np.std(iptms)) if iptms else 0,
            "best_interface_plddt": max(iface_plddts) if iface_plddts else 0,
            "mean_interface_plddt": float(np.mean(iface_plddts)) if iface_plddts else 0,
            "best_cross_pae": min(paes) if paes else None,
            "mean_cross_pae": float(np.mean(paes)) if paes else None,
            "best_model": best_model,
            "per_model": per_model,
        }
        all_results.append(result)

    # Print summary table
    print(f"{'Binder':<12} {'Best ipTM':<10} {'Mean ipTM':<10} {'Std':<8} "
          f"{'Iface pLDDT':<12} {'Cross PAE':<10} {'Iface Res':<10} {'Quality'}")
    print("-" * 130)

    for r in sorted(all_results, key=lambda x: x["best_iptm"], reverse=True):
        plddt_str = f"{r['mean_interface_plddt']:.3f}" if r['mean_interface_plddt'] > 0 else "—"
        pae_str = f"{r['mean_cross_pae']:.2f}" if r['mean_cross_pae'] is not None else "—"
        best_model = r.get("best_model", {})
        iface_n = f"{best_model.get('n_interface_res_a', 0)}+{best_model.get('n_interface_res_b', 0)}"

        # Quality assessment
        quality = []
        if r["best_iptm"] >= 0.85:
            quality.append("STRONG")
        elif r["best_iptm"] >= 0.70:
            quality.append("good")
        else:
            quality.append("weak")

        if r["mean_interface_plddt"] >= 0.80:
            quality.append("confident")
        elif r["mean_interface_plddt"] >= 0.60:
            quality.append("ok-conf")
        elif r["mean_interface_plddt"] > 0:
            quality.append("low-conf")

        if r["std_iptm"] < 0.03:
            quality.append("stable")
        elif r["std_iptm"] > 0.10:
            quality.append("VARIABLE")

        if r["mean_cross_pae"] is not None and r["mean_cross_pae"] < 5.0:
            quality.append("tight-iface")
        elif r["mean_cross_pae"] is not None and r["mean_cross_pae"] > 15.0:
            quality.append("loose-iface")

        print(f"{r['binder']:<12} {r['best_iptm']:<10.4f} {r['mean_iptm']:<10.4f} "
              f"{r['std_iptm']:<8.4f} {plddt_str:<12} {pae_str:<10} "
              f"{iface_n:<10} {', '.join(quality)}")

    # Per-model detail for top candidates
    print()
    print("=" * 130)
    print("PER-MODEL DETAIL")
    print("=" * 130)

    for r in sorted(all_results, key=lambda x: x["best_iptm"], reverse=True):
        if r["best_iptm"] < 0.70:
            continue
        print(f"\n  {r['binder']} — best ipTM={r['best_iptm']:.4f}, "
              f"iface pLDDT={r['mean_interface_plddt']:.3f}")
        for m in r["per_model"]:
            pae_str = f"{m['cross_pae_ab']:.2f}" if m['cross_pae_ab'] is not None else "—"
            print(f"    model_{m['model']}: ipTM={m['iptm']:.4f}  "
                  f"pLDDT={m['binder_plddt']:.3f}  "
                  f"iface_pLDDT={m['interface_plddt']:.3f}  "
                  f"cross_PAE={pae_str}  "
                  f"iface={m['n_interface_res_a']}+{m['n_interface_res_b']}")

    # Save
    out_path = BASE / "interface_quality_results.json"
    out_path.write_text(json.dumps(all_results, indent=2, default=str) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
