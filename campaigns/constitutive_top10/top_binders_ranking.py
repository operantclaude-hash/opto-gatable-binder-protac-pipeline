#!/usr/bin/env python3
"""
Rank top constitutive binders for empirical testing.

Aggregates data from:
  - Stage 1/2 binder funnel (original ipTM screening)
  - Binding validation (binder-only + fusion vs HAT and catcore)
  - Dark vs lit comparison (gating signal)

Outputs a ranked table of ~22 binders with key metrics to help select
the top 10 for empirical validation.

Ranking criteria (for constitutive binding, no VVD-degron):
  1. Binder-HAT ipTM (direct binding to isolated HAT domain)
  2. Binder-catcore ipTM (binding to full catalytic core context)
  3. Context consistency (HAT vs catcore agreement)
  4. min_iPAE (interface quality — lower is better)
  5. Binder size (smaller = easier to express, better AAV packaging)
"""

import json
from pathlib import Path

BASE = Path(__file__).parent

# Load binding validation results
bv_path = BASE / "binding_validation" / "analysis" / "binding_validation_results.json"
with open(bv_path) as f:
    bv_data = json.load(f)

# Binder sequences (from campaign config + YAMLs)
SEQUENCES = {
    "pep4s2": "SPEEEILEESEAGLKRAEEWGKVAEACRKAAAAGDAAAVAACQAQLDRFG",
    "pep13s3": "AAQSGFSPASEAAAAELLAA",
    "pep30s5": "APPGKTAAEDAARRAAAAAA",
    "pep25s1": "AKGVANSTTDAAANAARIGELAAARLAAGL",
    "pep21s5": "SIGSSSSGASGSLGAVTAED",
    "pep10s5": "FGSVTGSGAADAVDDALAAA",
    "pep24s3": "SLAEAERARQRARREAAHAA",
    "pep32s1": "MGPQVTRENAARAAARVAAR",
    "pep15s3": "LGKGSTPAATRSWEAMAAKL",
    "pep16s5": "PRMPAGGRTEAERRTSENAGRGYAENRARS",
    "pep17s4": "MSELTEVEGPGEASPAAKAH",
    "pep21s2": "SIGSSSSGGSGSLGSVSGSS",
    "pep29s5": "VEEWGRRGEEMADVHPGGRGAARMARAVLE",
    "pep32s2": "MGPAVTRENAARAAARVAAR",
    "pep36s1": "PVLEYDESGSSGEGGLGSVGRGLVEAREGGERGEEEDRRR",
    "pep48s5": "PIETGPGGAGKDAQRLDEKL",
    "nb153s2": "SLRERAAANSKEAAAQTVVGAARLRAAVDAVDPAVADAVAAQYTANRNTLNAANAALAAIAAESLDRIKATYEAEGSVAALLAAEEAEGRRQQAAVLAAREAASAASRAIGAIVAWLRAQ",
    "nb62s3": "DLAAAALEVCEETGAQSRAVAAANQAAYAAQIAALAARIAALRAEAAARAAELAARAEATGSEEVRRERERAEAAAAAELAALEAQRAAAERNRAATAERLAANEEGLAAERARLEAEAA",
    "nb39s3": "MVLLGAIEGGKEGVRNRREREELALAGEDERARELDRRAGEGGGAFGGLATGEDVHEAYKIAEKGYEKREKEIKARMAELEKTGEPAAAAQAEAWRKAMEEEKKAEKEANKDAGAAGVGG",
    "nb134s2": "SAREAWLEEARAALDEITAENLAGAAANRRAYAEQLAAARARKAALEAEDAALLARVRAECAATGDEERLREVEREIAARLAAIDKEIANLERALAHAEQLQVDIVRARERDLAELEAEA",
    "nb178s5": "SAAKQAGIENAKALGEEAVAATEANAAALAADAAARLAELEKEIEEYKKKGEEELKELKKKSKSEEEIEKKKKEYEEKIKEIEEQMKTLERRVAAARENVEKAKEGKVVGDKEDAAFFNA",
    "nb87s2": "SQEIIEAYEEIGRESEEMSRRNCEEAEAYAAEVRARRAEIEARVAEAEAVLAADAAGAAAAEARERLAAAEAELARAEVELRRAEGQLRWCKAAEERLKEMKVTGEKEAAELAAFDEKAK",
}


def compute_composite_score(row):
    """Composite score balancing binding strength, context consistency, and interface quality.

    Components (all normalized 0-1):
      - binder_HAT_iptm: primary binding metric (weight 0.35)
      - binder_catcore_iptm: contextual binding (weight 0.25)
      - context_consistency: 1 - |HAT - catcore| / max(HAT, catcore) (weight 0.15)
      - interface_quality: 1 - min_ipae / 15.0, clamped (weight 0.15)
      - size_bonus: smaller binders get slight bonus (weight 0.10)
    """
    hat = row["binder_HAT_iptm"]
    cat = row["binder_catcore_iptm"]

    # Context consistency: how well HAT and catcore agree
    max_val = max(hat, cat, 0.01)
    context = 1.0 - abs(hat - cat) / max_val

    # Interface quality from min_iPAE (lower = tighter interface)
    ipae = row.get("binder_HAT_min_ipae", 10.0)
    interface = max(0.0, min(1.0, 1.0 - ipae / 15.0))

    # Size bonus: 20aa=1.0, 30aa=0.9, 50aa=0.7, 120aa=0.2
    length = row["length"]
    size = max(0.0, min(1.0, 1.0 - (length - 20) / 125.0))

    score = (0.35 * hat + 0.25 * cat + 0.15 * context +
             0.15 * interface + 0.10 * size)
    return score


def main():
    rows = bv_data["raw_data"]

    # Compute composite scores
    for row in rows:
        row["composite"] = compute_composite_score(row)
        row["sequence"] = SEQUENCES.get(row["binder"], "???")
        row["avg_binding"] = (row["binder_HAT_iptm"] + row["binder_catcore_iptm"]) / 2

    # Sort by composite score
    rows.sort(key=lambda r: r["composite"], reverse=True)

    print("=" * 120)
    print("TOP CONSTITUTIVE BINDERS — RANKED FOR EMPIRICAL TESTING")
    print("Target: p300 HAT domain | No VVD-degron | Architecture C scaffold")
    print("=" * 120)
    print()

    # Main ranking table
    print(f"{'Rank':<5} {'Binder':<12} {'Size':<6} "
          f"{'HAT ipTM':<10} {'CatC ipTM':<10} {'Avg':<8} "
          f"{'min iPAE':<10} {'S1 ipTM':<9} {'Composite':<10} {'Notes'}")
    print("-" * 120)

    for i, row in enumerate(rows, 1):
        hat = row["binder_HAT_iptm"]
        cat = row["binder_catcore_iptm"]
        avg = row["avg_binding"]
        ipae = row.get("binder_HAT_min_ipae", 0)
        s1 = row["s1_iptm"]
        comp = row["composite"]
        length = row["length"]
        name = row["binder"]

        # Generate notes
        notes = []
        if hat > 0.85:
            notes.append("strong HAT")
        if cat > 0.85:
            notes.append("strong catcore")
        if hat < 0.60:
            notes.append("WEAK HAT")
        if cat < 0.40:
            notes.append("WEAK catcore")
        if ipae < 1.5:
            notes.append("tight iface")
        if abs(hat - cat) > 0.30:
            notes.append("context-dep")
        if length <= 20:
            notes.append("compact")
        elif length >= 120:
            notes.append("large")

        print(f"{i:<5} {name:<12} {length:>3}aa  "
              f"{hat:<10.4f} {cat:<10.4f} {avg:<8.4f} "
              f"{ipae:<10.3f} {s1:<9.3f} {comp:<10.4f} {', '.join(notes)}")

    # Detailed view with sequences
    print()
    print("=" * 120)
    print("BINDER SEQUENCES (for synthesis)")
    print("=" * 120)
    for i, row in enumerate(rows, 1):
        name = row["binder"]
        seq = row["sequence"]
        length = row["length"]
        print(f"  {i:>2}. {name:<12} ({length:>3}aa): {seq}")

    # Selection criteria summary
    print()
    print("=" * 120)
    print("SELECTION CRITERIA FOR TOP 10")
    print("=" * 120)
    print("""
  INCLUDE if:
    - Binder-HAT ipTM > 0.80 (strong direct binding to isolated HAT)
    - Binder-catcore ipTM > 0.55 (maintains binding in structural context)
    - min_iPAE < 5.0 (confident interface geometry)

  PREFER:
    - Smaller binders (20-30aa) for better expression and AAV packaging
    - Context-consistent binders (HAT ~ catcore ipTM) for robustness
    - Sequence diversity (avoid near-identical sequences like pep32s1/pep32s2)

  EXCLUDE if:
    - Binder-HAT ipTM < 0.60 (unreliable binding prediction)
    - Both HAT and catcore ipTM < 0.50 (no binding signal)
    - Nanobody with binder-HAT ipTM < 0.50 (poor revalidation)

  NOTE on fusion ipTM:
    - Fusion ipTM 0.2-0.4 is EXPECTED dilution from scaffold mass
    - Do NOT use fusion ipTM to exclude binders
    - Binder-only ipTM is the primary selection criterion
""")

    # Category breakdown
    print("=" * 120)
    print("BINDERS BY CATEGORY")
    print("=" * 120)

    strong = [r for r in rows if r["binder_HAT_iptm"] > 0.80]
    moderate = [r for r in rows if 0.60 <= r["binder_HAT_iptm"] <= 0.80]
    weak = [r for r in rows if r["binder_HAT_iptm"] < 0.60]

    print(f"\n  STRONG binder-HAT (ipTM > 0.80): {len(strong)} binders")
    for r in strong:
        print(f"    {r['binder']:<12} HAT={r['binder_HAT_iptm']:.3f}  "
              f"catcore={r['binder_catcore_iptm']:.3f}  {r['length']}aa")

    print(f"\n  MODERATE binder-HAT (0.60-0.80): {len(moderate)} binders")
    for r in moderate:
        print(f"    {r['binder']:<12} HAT={r['binder_HAT_iptm']:.3f}  "
              f"catcore={r['binder_catcore_iptm']:.3f}  {r['length']}aa")

    print(f"\n  WEAK binder-HAT (ipTM < 0.60): {len(weak)} binders")
    for r in weak:
        print(f"    {r['binder']:<12} HAT={r['binder_HAT_iptm']:.3f}  "
              f"catcore={r['binder_catcore_iptm']:.3f}  {r['length']}aa")

    # Save ranked results
    out = {
        "ranked_binders": [
            {
                "rank": i + 1,
                "binder": r["binder"],
                "length": r["length"],
                "sequence": r["sequence"],
                "binder_HAT_iptm": r["binder_HAT_iptm"],
                "binder_catcore_iptm": r["binder_catcore_iptm"],
                "avg_binding": r["avg_binding"],
                "binder_HAT_min_ipae": r.get("binder_HAT_min_ipae", None),
                "s1_iptm": r["s1_iptm"],
                "composite": r["composite"],
                "fusion_HAT_iptm": r["fusion_HAT_iptm"],
                "fusion_catcore_iptm": r["fusion_catcore_iptm"],
            }
            for i, r in enumerate(rows)
        ],
        "selection_notes": {
            "strong_binders_HAT_gt_0.80": len(strong),
            "moderate_binders_0.60_to_0.80": len(moderate),
            "weak_binders_HAT_lt_0.60": len(weak),
            "total_evaluated": len(rows),
        }
    }

    out_path = BASE / "top_binders_ranked.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nRanked results saved: {out_path}")


if __name__ == "__main__":
    main()
