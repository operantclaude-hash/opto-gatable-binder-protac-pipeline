#!/usr/bin/env python3
"""
Analyze dark-state vs lit-state Boltz-2 predictions for constitutive top-10 binders.

Compares dark (binder+LOV2-dark+E2pppE2+UBQ fused with target) vs lit (binder-only
with target) to assess whether binder-target binding is preserved when the
opto-gating scaffold is added in the dark conformation.

Binders:
  - pep4s2  (50aa peptide, seq2)
  - pep13s3 (20aa peptide #13, seq3)
  - pep30s5 (20aa peptide #30, seq5)
  - nb62s3  (120aa nanobody-scaffold #62, seq3)

Targets:
  - HAT     (p300 HAT domain, isolated)
  - catcore (p300 catalytic core, 4BHW-based)

Also includes dark-fold (single-chain) constructs for structural quality assessment,
and binder-only controls from binding_validation.

Note: ipTM 0.2-0.3 on fusion constructs is expected due to ipTM dilution from the
large scaffold — it does NOT indicate lack of binding. Binder-only ipTM > 0.85 proves
binding; the dark/lit comparison here checks relative change.
"""

import json
from pathlib import Path


# --- Paths ---
SCRIPT_DIR = Path(__file__).parent
DARK_RESULTS_DIR = SCRIPT_DIR / "boltz2_results"
LIT_RESULTS_DIR = SCRIPT_DIR.parent / "binding_validation" / "boltz2_results"

# --- Binder definitions ---
BINDERS = [
    {"name": "pep4s2",  "size": "50aa",  "type": "peptide"},
    {"name": "pep13s3", "size": "20aa",  "type": "peptide"},
    {"name": "pep30s5", "size": "20aa",  "type": "peptide"},
    {"name": "nb62s3",  "size": "120aa", "type": "nanobody-scaffold"},
]

TARGETS = ["HAT", "catcore"]

# --- Construct name mappings ---
def dark_construct_name(binder, target):
    """Dark state: binder fused with LOV2-dark + E2pppE2 linker + UBQ cap, 2-chain with target."""
    return f"{binder['name']}_{binder['size']}_dark_{target}"

def lit_construct_name(binder, target):
    """Lit state: binder-only (no scaffold), 2-chain with target — serves as binding control."""
    return f"{binder['name']}_fusion_{target}"

def binder_only_name(binder, target):
    """Binder-only control (from binding_validation): binder + target, no scaffold."""
    return f"{binder['name']}_binder_{target}"

def dark_fold_name(binder):
    """Dark fold: single-chain scaffold fold quality (no target)."""
    return f"{binder['name']}_{binder['size']}_E2pppE2_dark"


def extract_best_metrics(results_dir, name):
    """
    Extract confidence metrics from Boltz-2 results, selecting the best model
    by confidence_score across all diffusion samples.

    Returns dict with metrics or None if no results found.
    """
    rdir = results_dir / name
    if not rdir.exists():
        return None

    # Find all confidence JSONs (multiple diffusion samples)
    conf_files = sorted(rdir.rglob("confidence_*_model_*.json"))
    if not conf_files:
        conf_files = sorted(rdir.rglob("confidence_*.json"))
    if not conf_files:
        return None

    best_metrics = None
    best_confidence = -1.0

    for cf in conf_files:
        with open(cf) as f:
            conf = json.load(f)

        metrics = {
            "iptm": conf.get("protein_iptm", conf.get("iptm", 0.0)),
            "ptm": conf.get("ptm", 0.0),
            "plddt": conf.get("complex_plddt", 0.0),
            "iplddt": conf.get("complex_iplddt", 0.0),
            "confidence": conf.get("confidence_score", 0.0),
            "ipde": conf.get("complex_ipde", 0.0),
            "model_file": cf.name,
            "n_models": len(conf_files),
        }

        # Extract per-chain pTM if available
        chains_ptm = conf.get("chains_ptm", {})
        if chains_ptm:
            metrics["chain_ptms"] = {k: float(v) for k, v in chains_ptm.items()}

        # Extract pair-chain ipTM matrix
        pair_iptm = conf.get("pair_chains_iptm", {})
        if pair_iptm and "1" in pair_iptm:
            # For 2-chain: pair_chains_iptm["1"]["0"] = target-to-binder ipTM
            metrics["cross_iptm_1to0"] = float(pair_iptm["1"].get("0", 0.0))
            metrics["cross_iptm_0to1"] = float(pair_iptm["0"].get("1", 0.0))

        cs = metrics["confidence"]
        if cs > best_confidence:
            best_confidence = cs
            best_metrics = metrics

    return best_metrics


def print_separator(char="=", width=110):
    print(char * width)


def main():
    print_separator()
    print("DARK vs LIT STATE COMPARISON: Constitutive Top-10 Binders")
    print("Scaffold: [Binder]-[LOV2-dark]-[E2pppE2]-[UBQ] (dark) vs binder-only (lit)")
    print_separator()

    all_results = {
        "dark_target": {},
        "lit_target": {},
        "binder_only": {},
        "dark_fold": {},
    }

    # --- 1. Collect all metrics ---
    for binder in BINDERS:
        for target in TARGETS:
            # Dark state (scaffold + target)
            dk_name = dark_construct_name(binder, target)
            dk_metrics = extract_best_metrics(DARK_RESULTS_DIR, dk_name)
            all_results["dark_target"][dk_name] = dk_metrics

            # Lit state / fusion (binder fused to target context)
            lt_name = lit_construct_name(binder, target)
            lt_metrics = extract_best_metrics(LIT_RESULTS_DIR, lt_name)
            all_results["lit_target"][lt_name] = lt_metrics

            # Binder-only control
            bo_name = binder_only_name(binder, target)
            bo_metrics = extract_best_metrics(LIT_RESULTS_DIR, bo_name)
            all_results["binder_only"][bo_name] = bo_metrics

        # Dark fold (single chain, no target)
        df_name = dark_fold_name(binder)
        df_metrics = extract_best_metrics(DARK_RESULTS_DIR, df_name)
        all_results["dark_fold"][df_name] = df_metrics

    # --- 2. Print binder-only controls ---
    print()
    print_separator()
    print("BINDER-ONLY CONTROLS (binder + target, no scaffold)")
    print("These establish baseline binding — ipTM > 0.85 = strong binder")
    print_separator()
    print(f"{'Binder':<12} {'Target':<10} {'ipTM':<8} {'pLDDT':<8} {'Conf':<8} {'Models':<8}")
    print_separator("-")

    for binder in BINDERS:
        for target in TARGETS:
            bo_name = binder_only_name(binder, target)
            m = all_results["binder_only"].get(bo_name)
            if m:
                print(f"{binder['name']:<12} {target:<10} {m['iptm']:.4f}  "
                      f"{m['plddt']:.4f}  {m['confidence']:.4f}  {m['n_models']}")
            else:
                print(f"{binder['name']:<12} {target:<10} ---")

    # --- 3. Print dark vs lit comparison table ---
    print()
    print_separator()
    print("DARK vs LIT COMPARISON (2-chain: construct + target)")
    print("Dark = binder+LOV2-dark+E2pppE2+UBQ fused, with target")
    print("Lit  = binder+LOV2-lit+linker fused, with target (binding_validation)")
    print_separator()
    print(f"{'Binder':<12} {'Target':<10} {'Dark_ipTM':<10} {'Lit_ipTM':<10} "
          f"{'Delta':<8} {'Dark_pLDDT':<12} {'Lit_pLDDT':<10} {'Verdict':<14}")
    print_separator("-")

    comparison_rows = []
    for binder in BINDERS:
        for target in TARGETS:
            dk_name = dark_construct_name(binder, target)
            lt_name = lit_construct_name(binder, target)
            dk = all_results["dark_target"].get(dk_name)
            lt = all_results["lit_target"].get(lt_name)

            row = {
                "binder": binder["name"],
                "size": binder["size"],
                "target": target,
                "dark_construct": dk_name,
                "lit_construct": lt_name,
            }

            if dk and lt:
                delta = lt["iptm"] - dk["iptm"]
                if delta > 0.02:
                    verdict = "LIT > DARK"
                elif delta < -0.02:
                    verdict = "DARK > LIT"
                else:
                    verdict = "~EQUAL"

                row.update({
                    "dark_iptm": dk["iptm"],
                    "lit_iptm": lt["iptm"],
                    "delta_iptm": delta,
                    "dark_plddt": dk["plddt"],
                    "lit_plddt": lt["plddt"],
                    "dark_confidence": dk["confidence"],
                    "lit_confidence": lt["confidence"],
                    "verdict": verdict,
                })

                print(f"{binder['name']:<12} {target:<10} {dk['iptm']:.4f}    "
                      f"{lt['iptm']:.4f}    {delta:+.4f}  {dk['plddt']:.4f}      "
                      f"{lt['plddt']:.4f}    {verdict}")
            elif dk and not lt:
                row.update({"dark_iptm": dk["iptm"], "dark_plddt": dk["plddt"],
                            "dark_confidence": dk["confidence"], "verdict": "LIT MISSING"})
                print(f"{binder['name']:<12} {target:<10} {dk['iptm']:.4f}    "
                      f"---       ---     {dk['plddt']:.4f}      ---       LIT MISSING")
            elif lt and not dk:
                row.update({"lit_iptm": lt["iptm"], "lit_plddt": lt["plddt"],
                            "lit_confidence": lt["confidence"], "verdict": "DARK MISSING"})
                print(f"{binder['name']:<12} {target:<10} ---       "
                      f"{lt['iptm']:.4f}    ---     ---           "
                      f"{lt['plddt']:.4f}    DARK MISSING")
            else:
                row["verdict"] = "BOTH MISSING"
                print(f"{binder['name']:<12} {target:<10} ---       ---       "
                      f"---     ---           ---       BOTH MISSING")

            comparison_rows.append(row)

    # --- 4. Print dark fold quality ---
    print()
    print_separator()
    print("DARK FOLD QUALITY (single-chain, no target — structural quality only)")
    print("pTM > 0.5 and pLDDT > 0.6 = reasonable fold")
    print_separator()
    print(f"{'Construct':<35} {'pTM':<8} {'pLDDT':<8} {'Conf':<8} {'Models':<8}")
    print_separator("-")

    for binder in BINDERS:
        df_name = dark_fold_name(binder)
        m = all_results["dark_fold"].get(df_name)
        if m:
            print(f"{df_name:<35} {m['ptm']:.4f}  {m['plddt']:.4f}  "
                  f"{m['confidence']:.4f}  {m['n_models']}")
        else:
            print(f"{df_name:<35} ---")

    # --- 5. Cross-chain ipTM detail ---
    print()
    print_separator()
    print("CROSS-CHAIN ipTM DETAIL (pair_chains_iptm target->scaffold)")
    print("This shows how strongly the target chain interacts with the scaffold chain")
    print_separator()
    print(f"{'Construct':<35} {'State':<6} {'ipTM(0->1)':<12} {'ipTM(1->0)':<12} {'Overall':<10}")
    print_separator("-")

    for binder in BINDERS:
        for target in TARGETS:
            dk_name = dark_construct_name(binder, target)
            lt_name = lit_construct_name(binder, target)

            dk = all_results["dark_target"].get(dk_name)
            lt = all_results["lit_target"].get(lt_name)

            if dk and "cross_iptm_0to1" in dk:
                print(f"{dk_name:<35} {'DARK':<6} {dk['cross_iptm_0to1']:.4f}      "
                      f"{dk['cross_iptm_1to0']:.4f}      {dk['iptm']:.4f}")
            if lt and "cross_iptm_0to1" in lt:
                print(f"{lt_name:<35} {'LIT':<6} {lt['cross_iptm_0to1']:.4f}      "
                      f"{lt['cross_iptm_1to0']:.4f}      {lt['iptm']:.4f}")

    # --- 6. Summary verdict ---
    print()
    print_separator()
    print("SUMMARY VERDICT")
    print_separator()

    valid_rows = [r for r in comparison_rows
                  if "dark_iptm" in r and "lit_iptm" in r]

    if not valid_rows:
        print("  No complete dark/lit pairs found for comparison.")
    else:
        lit_better = [r for r in valid_rows if r.get("verdict") == "LIT > DARK"]
        dark_better = [r for r in valid_rows if r.get("verdict") == "DARK > LIT"]
        equal = [r for r in valid_rows if r.get("verdict") == "~EQUAL"]

        print(f"  Total comparisons: {len(valid_rows)}")
        print(f"  Lit > Dark:  {len(lit_better)}  (scaffold removal improves binding)")
        print(f"  Dark > Lit:  {len(dark_better)}  (scaffold does NOT block binding)")
        print(f"  ~Equal:      {len(equal)}  (within noise threshold +/-0.02)")
        print()

        # Average deltas
        avg_delta = sum(r["delta_iptm"] for r in valid_rows) / len(valid_rows)
        print(f"  Mean delta(lit-dark) ipTM: {avg_delta:+.4f}")

        hat_rows = [r for r in valid_rows if r["target"] == "HAT"]
        cat_rows = [r for r in valid_rows if r["target"] == "catcore"]

        if hat_rows:
            avg_hat = sum(r["delta_iptm"] for r in hat_rows) / len(hat_rows)
            print(f"  Mean delta for HAT target:     {avg_hat:+.4f}")
        if cat_rows:
            avg_cat = sum(r["delta_iptm"] for r in cat_rows) / len(cat_rows)
            print(f"  Mean delta for catcore target:  {avg_cat:+.4f}")

        print()
        print("  INTERPRETATION:")
        print("  - ipTM 0.2-0.3 on fusion constructs is expected (ipTM dilution from scaffold)")
        print("  - Binder-only ipTM > 0.85 confirms strong binding (see controls above)")
        print("  - Dark vs Lit delta indicates whether scaffold addition affects binding signal")
        print("  - Small delta (~0) = scaffold does not disrupt predicted binding interface")
        if avg_delta > 0.02:
            print("  >>> Lit state shows BETTER binding — scaffold may partially occlude target")
        elif avg_delta < -0.02:
            print("  >>> Dark state shows BETTER binding — unexpected, may be artifact")
        else:
            print("  >>> Dark and Lit states show SIMILAR binding — scaffold is not disrupting interface")

    # --- 7. Save results ---
    output = {
        "analysis": "dark_vs_lit_comparison",
        "description": "Compares dark-state (scaffold+target) vs lit-state (binder-only+target) Boltz-2 predictions",
        "binders": [b["name"] for b in BINDERS],
        "targets": TARGETS,
        "comparison": comparison_rows,
        "binder_only_controls": {},
        "dark_fold_quality": {},
        "all_metrics": {},
    }

    # Add binder-only controls
    for binder in BINDERS:
        for target in TARGETS:
            bo_name = binder_only_name(binder, target)
            m = all_results["binder_only"].get(bo_name)
            if m:
                # Remove non-serializable items
                clean = {k: v for k, v in m.items() if k != "model_file"}
                output["binder_only_controls"][bo_name] = clean

    # Add dark fold quality
    for binder in BINDERS:
        df_name = dark_fold_name(binder)
        m = all_results["dark_fold"].get(df_name)
        if m:
            clean = {k: v for k, v in m.items() if k != "model_file"}
            output["dark_fold_quality"][df_name] = clean

    # Add all raw metrics
    for category, entries in all_results.items():
        output["all_metrics"][category] = {}
        for name, m in entries.items():
            if m:
                clean = {k: v for k, v in m.items() if k != "model_file"}
                output["all_metrics"][category][name] = clean

    out_path = SCRIPT_DIR / "dark_vs_lit_results.json"
    out_path.write_text(json.dumps(output, indent=2, default=str) + "\n")
    print(f"\nResults saved: {out_path}")


if __name__ == "__main__":
    main()
