#!/usr/bin/env python3
"""
Generate GenScript order sheet for top 10 constitutive-open constructs.

Construct: M-[Binder]-[GGSGGS]-[LOV2_lit(121aa)]-[E2pppE2(22aa)]-[UBQ(76aa)]
Purpose: Test constitutive binding to p300 HAT domain via BLI
"""

from pathlib import Path

LOV2_LIT = (
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRA"
    "TVRKIRDAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLD"
    "GTEHVRDAAEREG"
)
E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"
UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
SPACER = "GGSGGS"

# Final 10 in priority order
CONSTRUCTS = [
    ("pep9s5",   "MTPRKPSDAGLVRAAEKAAA",                                          20, 0.635, "BINDER-DRIVEN"),
    ("pep21s5",  "SIGSSSSGASGSLGAVTAED",                                          20, 0.685, "MIXED"),
    ("pep16s4",  "AGMPLGGRTEEELRTSRNAAEGDAENRARS",                                30, 0.649, "MIXED"),
    ("pep11s1",  "SSRPNPPRLDARIEESRAGL",                                          20, 0.640, "MIXED"),
    ("pep37s3",  "ETPADRDAAGGRAAAARLNL",                                          20, 0.629, "MIXED"),
    ("pep41s5",  "DDVERGRAAGAAGTGARGEF",                                          20, 0.601, "MIXED"),
    ("pep19s2",  "ATVADPESGGSAVTAGAPVE",                                          20, 0.587, "MIXED"),
    ("pep21s2",  "SIGSSSSGGSGSLGSVSGSS",                                          20, 0.707, "BINDER-DRIVEN"),
    ("pep29s4",  "MMKVVFGGGCGIGTPISEEDIRRGNIGGKEGGAESKARSAARAAAAAAAL",            50, 0.678, "BINDER-DRIVEN"),
    ("pep36s1",  "PVLEYDESGSSGEGGLGSVGRGLVEAREGGERGEEEDRRR",                      40, 0.511, "MIXED"),
]


def build_construct(binder_seq):
    """Build full constitutive-open construct with N-terminal Met."""
    full = binder_seq + SPACER + LOV2_LIT + E2PPPE2 + UBQ
    # Prepend Met if not already present
    if full[0] != 'M':
        full = 'M' + full
    return full


def main():
    base = Path(__file__).parent

    print("=" * 100)
    print("GENSCRIPT ORDER — 10 CONSTITUTIVE-OPEN p300 BINDER CONSTRUCTS")
    print("=" * 100)
    print()
    print("Architecture: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]-[UBQ]")
    print("Purpose: BLI binding assay against p300 HAT domain")
    print("Expression: Cell-free recommended (no glycosylation, no disulfide requirement)")
    print("Tag: C-terminal FLAG (p300 target has N-terminal His)")
    print()

    fasta_lines = []
    order_lines = []

    for i, (name, binder_seq, binder_len, score, driver) in enumerate(CONSTRUCTS, 1):
        full_seq = build_construct(binder_seq)
        construct_len = len(full_seq)
        has_added_met = binder_seq[0] != 'M'

        construct_name = f"OPR-HAT-{i:02d}_{name}"

        print(f"  {i:>2}. {construct_name}")
        print(f"      Binder: {binder_seq} ({binder_len}aa)")
        print(f"      Total:  {construct_len}aa ({construct_len * 3}bp coding)")
        print(f"      Score:  {score:.3f}  |  Cap-trunc: {driver}")
        if has_added_met:
            print(f"      Note:   Met prepended for expression initiation")
        if name == "pep29s4":
            print(f"      Note:   Contains 1 Cys in binder (pos 11) — flag for GenScript")
        if name == "pep21s2":
            print(f"      Note:   High Gly/Ser binder — disordered, folds on binding")
        print()

        # FASTA
        fasta_lines.append(f">{construct_name} | {construct_len}aa | score={score:.3f} | {driver}")
        for j in range(0, len(full_seq), 80):
            fasta_lines.append(full_seq[j:j+80])

        # Order table
        order_lines.append({
            "name": construct_name,
            "sequence": full_seq,
            "length_aa": construct_len,
            "length_bp": construct_len * 3,
            "binder": name,
            "binder_length": binder_len,
            "score": score,
            "driver": driver,
        })

    # Write FASTA
    fasta_path = base / "genscript_constructs.fasta"
    fasta_path.write_text("\n".join(fasta_lines) + "\n")
    print(f"FASTA saved: {fasta_path}")

    # Write CSV for order
    csv_path = base / "genscript_order.csv"
    csv_lines = ["Name,Sequence,Length_AA,Length_BP,Binder,Binder_Length,Score,CapTrunc_Driver"]
    for o in order_lines:
        csv_lines.append(
            f"{o['name']},{o['sequence']},{o['length_aa']},{o['length_bp']},"
            f"{o['binder']},{o['binder_length']},{o['score']:.3f},{o['driver']}"
        )
    csv_path.write_text("\n".join(csv_lines) + "\n")
    print(f"CSV saved: {csv_path}")

    # Summary
    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"  Constructs:  {len(order_lines)}")
    print(f"  Size range:  {min(o['length_aa'] for o in order_lines)}-{max(o['length_aa'] for o in order_lines)}aa")
    print(f"  BP range:    {min(o['length_bp'] for o in order_lines)}-{max(o['length_bp'] for o in order_lines)}bp")
    print(f"  BINDER-DRIVEN: {sum(1 for o in order_lines if o['driver'] == 'BINDER-DRIVEN')}")
    print(f"  MIXED:         {sum(1 for o in order_lines if o['driver'] == 'MIXED')}")
    print()
    print("  For GenScript:")
    print("    - Expression: cell-free (recommended) or E. coli")
    print("    - Codon optimization: human or E. coli depending on system")
    print("    - Tag: C-terminal FLAG-tag (-DYKDDDDK)")
    print("    - pep29s4 has 1 unpaired Cys — use reducing buffer")
    print("    - No signal peptide needed (intracellular target)")


if __name__ == "__main__":
    main()
