#!/usr/bin/env python3
"""
Full analysis of expanded pool binder-alone + HAT predictions.

Runs all validation layers:
  Layer 1: Epitope contact analysis (hotspot overlap)
  Layer 3: Multi-sample consistency (across 5 diffusion samples)
  Layer 4: Interface quality (ipTM, pLDDT, PAE)

Identifies candidates worth advancing to catcore validation.
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
RESULTS_DIR = BASE / "boltz2_results"

HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}
CONTACT_CUTOFF = 4.0

BINDERS = {
    "pep6s1": 30, "pep11s1": 20, "pep49s5": 30, "pep9s5": 20,
    "pep29s4": 50, "pep4s4": 20, "pep30s4": 30, "pep41s5": 20,
    "pep34s1": 50, "pep46s3": 30, "pep37s3": 20, "pep19s2": 20,
    "pep16s4": 30, "pep40s3": 40, "pep43s1": 20, "pep14s4": 20,
}

SCREEN_IPTM = {
    "pep6s1": 0.9049, "pep11s1": 0.8857, "pep49s5": 0.8848, "pep9s5": 0.8821,
    "pep29s4": 0.8815, "pep4s4": 0.8765, "pep30s4": 0.8711, "pep41s5": 0.8689,
    "pep34s1": 0.8687, "pep46s3": 0.8673, "pep37s3": 0.8666, "pep19s2": 0.8650,
    "pep16s4": 0.8635, "pep40s3": 0.8624, "pep43s1": 0.8591, "pep14s4": 0.8583,
}

SCREEN_HOTSPOTS = {
    "pep6s1": 12, "pep11s1": 11, "pep49s5": 11, "pep9s5": 11,
    "pep29s4": 11, "pep4s4": 11, "pep30s4": 10, "pep41s5": 11,
    "pep34s1": 10, "pep46s3": 12, "pep37s3": 11, "pep19s2": 11,
    "pep16s4": 11, "pep40s3": 11, "pep43s1": 12, "pep14s4": 10,
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
                atom_name = parts[3]
                atoms.append((chain, resnum, x, y, z, atom_name))
            except (ValueError, IndexError):
                continue
    return atoms


def find_contacts(atoms, cutoff=4.0):
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


def find_result_dir(binder_name, binder_len):
    patterns = [
        f"{binder_name}_{binder_len}aa_binder_HAT",
        f"{binder_name}_binder_HAT",
    ]
    for pat in patterns:
        d = RESULTS_DIR / pat
        if d.exists():
            return d
    return None


def main():
    print("=" * 140)
    print("EXPANDED POOL ANALYSIS — 16 NEW BINDERS (BINDER-ALONE + HAT)")
    print(f"Contact cutoff: {CONTACT_CUTOFF}A  |  Hotspots: {sorted(HOTSPOTS)}")
    print("=" * 140)
    print()

    all_results = []

    for binder, binder_len in sorted(BINDERS.items()):
        result_dir = find_result_dir(binder, binder_len)
        if not result_dir:
            print(f"{binder}: NO RESULTS FOUND")
            continue

        cifs = sorted(result_dir.rglob("*model_*.cif"))
        if not cifs:
            continue

        per_model = []
        for cif in cifs:
            model_idx = cif.stem.split("_model_")[1] if "_model_" in cif.stem else "0"

            atoms = parse_cif_atoms(cif)
            contacts = find_contacts(atoms, CONTACT_CUTOFF)
            hs_contacts = contacts & HOTSPOTS

            # Get confidence
            conf_files = sorted(result_dir.rglob(f"*confidence*{model_idx}*"))
            iptm = 0
            plddt = 0
            if conf_files:
                with open(conf_files[0]) as f:
                    conf = json.load(f)
                iptm = conf.get("protein_iptm", conf.get("iptm", 0))
                plddt = conf.get("complex_plddt", 0)

            # Interface pLDDT from B-factors
            bfactors = {}
            with open(cif) as f:
                for line in f:
                    if not line.startswith("ATOM"):
                        continue
                    parts = line.split()
                    if len(parts) < 18:
                        continue
                    try:
                        if parts[3] == "CA" and parts[9] == "A":
                            bfactors[int(parts[6])] = float(parts[17])
                    except (ValueError, IndexError):
                        continue

            binder_plddt = np.mean(list(bfactors.values())) if bfactors else 0

            per_model.append({
                "model": model_idx,
                "iptm": iptm,
                "total_contacts": len(contacts),
                "hotspot_contacts": len(hs_contacts),
                "hotspots_hit": sorted(hs_contacts),
                "binder_plddt": float(binder_plddt),
            })

        # Aggregate
        iptms = [m["iptm"] for m in per_model]
        hs_counts = [m["hotspot_contacts"] for m in per_model]

        # Consistency: which hotspots hit in >= 60% of models
        hotspot_freq = {}
        for h in HOTSPOTS:
            count = sum(1 for m in per_model if h in m["hotspots_hit"])
            if count > 0:
                hotspot_freq[h] = count

        n_models = len(per_model)
        consistent_hs = [h for h, c in hotspot_freq.items()
                         if c >= max(2, n_models * 0.6)]
        models_with_hs = sum(1 for c in hs_counts if c > 0)

        # Epitope verdict
        best_hs = max(hs_counts) if hs_counts else 0
        best_pct = best_hs / len(HOTSPOTS) * 100

        if best_pct >= 50:
            verdict = "ON-TARGET"
        elif best_pct >= 25:
            verdict = "PARTIAL"
        elif best_hs > 0:
            verdict = "MARGINAL"
        else:
            verdict = "OFF-TARGET"

        result = {
            "binder": binder,
            "binder_length": binder_len,
            "screen_iptm": SCREEN_IPTM.get(binder, 0),
            "screen_hotspots": SCREEN_HOTSPOTS.get(binder, 0),
            "best_iptm": max(iptms) if iptms else 0,
            "mean_iptm": float(np.mean(iptms)) if iptms else 0,
            "std_iptm": float(np.std(iptms)) if iptms else 0,
            "best_hotspots": best_hs,
            "mean_hotspots": float(np.mean(hs_counts)) if hs_counts else 0,
            "consistency": models_with_hs / n_models if n_models > 0 else 0,
            "n_consistent_hotspots": len(consistent_hs),
            "consistent_hotspots": sorted(consistent_hs),
            "verdict": verdict,
            "mean_binder_plddt": float(np.mean([m["binder_plddt"] for m in per_model])),
            "per_model": per_model,
        }
        all_results.append(result)

    # Summary table
    print(f"{'Binder':<10} {'Len':<5} {'Scrn ipTM':<10} {'Val ipTM':<10} {'Std':<7} "
          f"{'BestHS':<7} {'MeanHS':<8} {'Consist':<8} {'nConsHS':<8} {'pLDDT':<7} "
          f"{'Verdict':<12} {'Advance?'}")
    print("-" * 140)

    # Sort by a quick composite
    all_results.sort(key=lambda r: (
        r["n_consistent_hotspots"] * 0.4 +
        r["best_iptm"] * 0.3 +
        r["consistency"] * 0.3
    ), reverse=True)

    advance = []
    for r in all_results:
        should_advance = (
            r["verdict"] in ("ON-TARGET", "PARTIAL") and
            r["consistency"] >= 0.60 and
            r["best_iptm"] >= 0.50 and
            r["n_consistent_hotspots"] >= 3
        )

        flag = "YES → catcore" if should_advance else "no"
        if should_advance:
            advance.append(r)

        hs_str = ",".join(str(h) for h in r["consistent_hotspots"][:6])
        print(f"{r['binder']:<10} {r['binder_length']:>3}  "
              f"{r['screen_iptm']:<10.4f} {r['best_iptm']:<10.4f} {r['std_iptm']:<7.3f} "
              f"{r['best_hotspots']:<7} {r['mean_hotspots']:<8.1f} "
              f"{r['consistency']:<8.0%} {r['n_consistent_hotspots']:<8} "
              f"{r['mean_binder_plddt']:<7.1f} "
              f"{r['verdict']:<12} {flag}")

    # Advancement summary
    print()
    print("=" * 140)
    print(f"ADVANCE TO CATCORE VALIDATION: {len(advance)}/{len(all_results)} binders")
    print("=" * 140)

    for r in advance:
        print(f"  {r['binder']:<10} {r['binder_length']:>3}aa  "
              f"ipTM={r['best_iptm']:.3f}  HS={r['best_hotspots']}/14  "
              f"consistent={r['n_consistent_hotspots']}  "
              f"pLDDT={r['mean_binder_plddt']:.1f}")

    not_advanced = [r for r in all_results if r not in advance]
    if not_advanced:
        print(f"\n  NOT ADVANCING ({len(not_advanced)}):")
        for r in not_advanced:
            reason = []
            if r["verdict"] in ("OFF-TARGET", "MARGINAL"):
                reason.append(f"epitope={r['verdict']}")
            if r["consistency"] < 0.60:
                reason.append(f"consistency={r['consistency']:.0%}")
            if r["best_iptm"] < 0.50:
                reason.append(f"ipTM={r['best_iptm']:.3f}")
            if r["n_consistent_hotspots"] < 3:
                reason.append(f"consistHS={r['n_consistent_hotspots']}")
            print(f"    {r['binder']:<10} — {', '.join(reason)}")

    # Save
    out_path = BASE / "expanded_pool_results.json"
    out_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nResults saved: {out_path}")

    # Generate catcore validation YAMLs for advancing binders
    if advance:
        print()
        print("=" * 140)
        print("GENERATING CATCORE VALIDATION YAMLS")
        print("=" * 140)

        catcore_yaml_dir = BASE / "catcore_yamls"
        catcore_yaml_dir.mkdir(exist_ok=True)

        # Read catcore sequence from existing campaign
        catcore_seq = None
        campaign_base = BASE.parent
        bv_yamls = sorted(campaign_base.rglob("*binder_catcore*.yaml"))
        if bv_yamls:
            with open(bv_yamls[0]) as f:
                for line in f:
                    if "sequence:" in line:
                        seq = line.strip().split("sequence:")[1].strip()
                        if len(seq) > 400:  # catcore is ~563aa
                            catcore_seq = seq
                            break

        if not catcore_seq:
            print("WARNING: Could not find catcore sequence from existing YAMLs")
            # Fallback: check for catcore sequence in binding_validation
            bv_yamls2 = sorted(Path(campaign_base / "binding_validation").rglob("*catcore*.yaml"))
            if bv_yamls2:
                with open(bv_yamls2[0]) as f:
                    for line in f:
                        if "sequence:" in line:
                            seq = line.strip().split("sequence:")[1].strip()
                            if len(seq) > 400:
                                catcore_seq = seq
                                break

        if catcore_seq:
            # Get binder sequences from the prediction YAMLs
            binder_seqs = {}
            for r in advance:
                yaml_dir = BASE / "boltz2_yamls"
                yaml_files = sorted(yaml_dir.glob(f"{r['binder']}_*_binder_HAT.yaml"))
                if yaml_files:
                    with open(yaml_files[0]) as f:
                        for line in f:
                            if "sequence:" in line and "KFSAKRL" not in line:
                                binder_seqs[r["binder"]] = line.strip().split("sequence:")[1].strip()
                                break

            count = 0
            run_lines = ["#!/bin/bash", "set -e", f"cd {BASE}", ""]
            for r in advance:
                name = r["binder"]
                blen = r["binder_length"]
                bseq = binder_seqs.get(name, "")
                if not bseq:
                    continue

                yaml_name = f"{name}_{blen}aa_binder_catcore"
                yaml_path = catcore_yaml_dir / f"{yaml_name}.yaml"
                yaml_path.write_text(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {bseq}
      msa: empty
  - protein:
      id: B
      sequence: {catcore_seq}
      msa: empty
""")
                out_dir = BASE / "catcore_results" / yaml_name
                run_lines.append(
                    f"echo 'Running {yaml_name}...' && "
                    f"conda run -n boltz_only boltz predict {yaml_path} "
                    f"--out_dir {out_dir} "
                    f"--recycling_steps 3 --diffusion_samples 5 --model boltz2"
                )
                count += 1

            run_lines.append("")
            run_lines.append("echo 'All catcore validations complete!'")
            run_script = BASE / "run_catcore_predictions.sh"
            run_script.write_text("\n".join(run_lines) + "\n")
            run_script.chmod(0o755)
            print(f"Generated {count} catcore YAMLs in {catcore_yaml_dir}")
            print(f"Run script: {run_script}")
        else:
            print("ERROR: Could not find catcore sequence. Manual setup needed.")


if __name__ == "__main__":
    main()
