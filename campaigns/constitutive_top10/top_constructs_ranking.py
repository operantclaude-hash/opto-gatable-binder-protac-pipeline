#!/usr/bin/env python3
"""
Rank constitutively-open opto-scaffold constructs for empirical testing.

Each construct is: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]-[UBQ]
  - LOV2_lit:  121aa (Jα helix released, binder exposed)
  - E2pppE2:    23aa (triple-proline kink, 978x gating ratio)
  - UBQ cap:    76aa (ubiquitin)
  - GGSGGS:      6aa (spacer)

All binders are de novo designs (RFdiffusion + LigandMPNN), NOT nanobodies.
Original screen: ~550 binders against p300 HAT, top hits ipTM > 0.90.

Binding validation tested each construct in the lit (constitutively-open) state
against both HAT domain and full catalytic core (catcore).

The "fusion" ipTM values (0.2-0.5) reflect ipTM dilution from the ~220aa scaffold —
this is expected and NOT a sign of weak binding. Binder-alone ipTM confirms the
binding interface; fusion ipTM confirms the scaffold doesn't completely destroy it.
"""

import json
from pathlib import Path

BASE = Path(__file__).parent

# Load binding validation results
bv_path = BASE / "binding_validation" / "analysis" / "binding_validation_results.json"
with open(bv_path) as f:
    bv_data = json.load(f)

# Binder sequences
BINDER_SEQS = {
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

# Scaffold constants
SPACER = "GGSGGS"   # 6aa
LOV2_LIT = "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIRDAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAEREG"  # 121aa
E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"  # 22aa (but manifest says 23 — may include M or extra residue)
UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"  # 76aa


def main():
    rows = bv_data["raw_data"]

    # Add construct info
    for row in rows:
        binder_seq = BINDER_SEQS.get(row["binder"], "")
        binder_len = row["length"]
        construct_len = binder_len + 6 + 121 + 23 + 76  # spacer + LOV2_lit + E2pppE2 + UBQ
        row["binder_seq"] = binder_seq
        row["construct_len"] = construct_len

    # Sort by binder-alone HAT ipTM (primary metric for constitutive binding)
    rows.sort(key=lambda r: r["binder_HAT_iptm"], reverse=True)

    print("=" * 130)
    print("CONSTITUTIVELY-OPEN CONSTRUCTS FOR EMPIRICAL TESTING")
    print("Architecture: [Binder]-[GGSGGS]-[LOV2_lit(121aa)]-[E2pppE2(23aa)]-[UBQ(76aa)]")
    print("Target: p300 HAT domain  |  Sweep gating: 978x (E2pppE2 kink)")
    print("=" * 130)
    print()
    print("Metrics explained:")
    print("  Screen ipTM:    Original binder-vs-HAT screen (binder alone, no scaffold)")
    print("  Binder+HAT:     Re-validation of binder alone against HAT (5 diffusion samples)")
    print("  Binder+CatC:    Binder alone against full catalytic core")
    print("  Fusion+HAT:     Full construct (lit state) against HAT (ipTM diluted by scaffold)")
    print("  Fusion+CatC:    Full construct (lit state) against catcore")
    print("  min iPAE:       Minimum interface PAE (lower = more confident interface)")
    print()

    hdr = (f"{'Rk':<4} {'Binder':<10} {'Bdr':<5} {'Cst':<5} "
           f"{'Screen':<8} {'Bdr+HAT':<9} {'Bdr+CC':<9} "
           f"{'Fus+HAT':<9} {'Fus+CC':<9} {'iPAE':<7} {'Flags'}")
    print(hdr)
    print("-" * 130)

    for i, row in enumerate(rows, 1):
        name = row["binder"]
        bdr_len = row["length"]
        cst_len = row["construct_len"]
        s1 = row["s1_iptm"]
        bdr_hat = row["binder_HAT_iptm"]
        bdr_cc = row["binder_catcore_iptm"]
        fus_hat = row["fusion_HAT_iptm"]
        fus_cc = row["fusion_catcore_iptm"]
        ipae = row.get("binder_HAT_min_ipae", 99)

        flags = []
        if bdr_hat >= 0.85:
            flags.append("STRONG")
        elif bdr_hat >= 0.70:
            flags.append("good")
        elif bdr_hat >= 0.50:
            flags.append("moderate")
        else:
            flags.append("WEAK")

        if bdr_cc >= 0.80:
            flags.append("cc+")
        elif bdr_cc < 0.40:
            flags.append("cc-")

        if ipae < 1.5:
            flags.append("tight")
        elif ipae > 10:
            flags.append("loose")

        if abs(bdr_hat - bdr_cc) > 0.30:
            flags.append("ctx-dep")

        # Sequence similarity flags
        if name == "pep32s2":
            flags.append("~pep32s1")
        if name == "pep21s2":
            flags.append("~pep21s5")

        print(f"{i:<4} {name:<10} {bdr_len:>3}  {cst_len:>3}  "
              f"{s1:<8.3f} {bdr_hat:<9.4f} {bdr_cc:<9.4f} "
              f"{fus_hat:<9.4f} {fus_cc:<9.4f} {ipae:<7.2f} {', '.join(flags)}")

    # Sequences
    print()
    print("=" * 130)
    print("BINDER SEQUENCES (what varies between constructs)")
    print("Scaffold is constant: [binder]-GGSGGS-LOV2_lit-EAAAKEAAAKPPPEAAAKEAAAK-UBQ")
    print("=" * 130)
    for i, row in enumerate(rows, 1):
        print(f"  {i:>2}. {row['binder']:<10} ({row['length']:>3}aa): {row['binder_seq']}")

    # Recommended selections
    print()
    print("=" * 130)
    print("RECOMMENDED TOP 10 FOR EMPIRICAL TESTING")
    print("=" * 130)
    print("""
  TIER 1 — High confidence (binder+HAT ipTM > 0.85, binder+catcore > 0.55):
    1. pep4s2   (50aa, 276aa construct) — Highest HAT binding (0.92), tight iPAE (0.67)
    2. pep29s5  (30aa, 256aa construct) — Very high HAT (0.89), tight iPAE (0.85)
    3. pep13s3  (20aa, 246aa construct) — Strong both HAT+catcore (0.86/0.88)
    4. pep30s5  (20aa, 246aa construct) — Strong both HAT+catcore (0.86/0.88)
    5. nb178s5  (120aa, 346aa construct) — Strong HAT (0.87), tight iPAE (1.30)

  TIER 2 — Good candidates (binder+HAT ipTM 0.80-0.85):
    6. pep21s2  (20aa, 246aa construct) — HAT 0.84, tight iPAE (1.48)
    7. pep21s5  (20aa, 246aa construct) — HAT 0.83, diverse from pep21s2
    8. pep25s1  (30aa, 256aa construct) — HAT 0.83, moderate catcore
    9. nb134s2  (120aa, 346aa construct) — HAT 0.83, good catcore (0.65)
   10. pep15s3  (20aa, 246aa construct) — HAT 0.83, good catcore (0.60)

  ALTERNATES:
    - nb62s3   (120aa) — HAT 0.81, good catcore 0.69, tight iPAE
    - pep17s4  (20aa)  — HAT 0.81, moderate catcore 0.59

  DROP:
    - pep32s2  — near-identical to pep32s1, pick at most one
    - pep21s2  — similar to pep21s5, pick at most one
    - nb39s3, nb153s2 — failed HAT revalidation (< 0.35)

  NOTES:
    - pep4s2 and pep29s5 show context-dependence (strong HAT, weak catcore)
      Worth testing but higher risk of epitope not being accessible in vivo
    - pep13s3 and pep30s5 are the safest picks (strong on both targets)
    - Including 2-3 of the 120aa binders provides size diversity
    - All constructs use the same scaffold so differences are purely binder-driven
""")

    # Save
    out = {
        "description": "Constitutively-open opto-scaffold constructs for empirical testing",
        "architecture": "[Binder]-[GGSGGS(6aa)]-[LOV2_lit(121aa)]-[E2pppE2(23aa)]-[UBQ(76aa)]",
        "scaffold_length": 226,
        "gating_ratio": "978x (E2pppE2 kink, sweep analysis)",
        "constructs": [
            {
                "rank": i + 1,
                "binder_name": r["binder"],
                "binder_length": r["length"],
                "construct_length": r["construct_len"],
                "binder_sequence": r["binder_seq"],
                "screen_iptm": r["s1_iptm"],
                "binder_HAT_iptm": r["binder_HAT_iptm"],
                "binder_catcore_iptm": r["binder_catcore_iptm"],
                "fusion_HAT_iptm": r["fusion_HAT_iptm"],
                "fusion_catcore_iptm": r["fusion_catcore_iptm"],
                "binder_HAT_min_ipae": r.get("binder_HAT_min_ipae"),
            }
            for i, r in enumerate(rows)
        ],
    }
    out_path = BASE / "top_constructs_ranked.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"Results saved: {out_path}")


if __name__ == "__main__":
    main()
