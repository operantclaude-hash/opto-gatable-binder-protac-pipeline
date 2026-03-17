#!/usr/bin/env python3
"""
Expression feasibility filter for top-10 constitutive constructs.

Checks:
  1. Furin cleavage sites (R-X-[K/R]-R) — problematic in CHO cells
  2. Net charge at pH 7.4
  3. Isoelectric point estimate
  4. GRAVY (hydropathicity) — high = aggregation risk
  5. Instability index (Guruprasad et al. 1990)
  6. Low-complexity / disorder-prone regions (poly-G, poly-S, poly-A runs)
  7. Cysteine content (unpaired Cys = misfolding in oxidizing environments)
  8. Methionine at N-terminus (expression initiation)
  9. Rare codon clusters (consecutive rare codons in E. coli)
 10. Aggregation-prone hydrophobic patches (5+ consecutive hydrophobic aa)

Also builds the full construct sequences and checks junctions.
"""

import json
import re
from pathlib import Path

BASE = Path(__file__).parent

# Amino acid properties
AA_CHARGE_PH7 = {
    'D': -1, 'E': -1,  # negative
    'K': 1, 'R': 1,     # positive
    'H': 0.1,           # partially protonated at pH 7.4
}

AA_HYDROPATHY = {  # Kyte-Doolittle
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'D': -3.5,
    'E': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5,
}

# Instability index dipeptide weights (Guruprasad et al. 1990)
# Using simplified version — positive DIWV values indicate instability
INSTABILITY_WEIGHTS = {}  # Would need full 400-entry table; use heuristic instead

HYDROPHOBIC = set('IVLFMWYA')

LOV2_LIT = ("LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRA"
            "TVRKIRDAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLD"
            "GTEHVRDAAEREG")
E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"
UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
SPACER = "GGSGGS"


def find_furin_sites(seq):
    """Find R-X-[KR]-R furin cleavage motifs."""
    sites = []
    for i in range(len(seq) - 3):
        if seq[i] == 'R' and seq[i+2] in ('K', 'R') and seq[i+3] == 'R':
            sites.append((i + 1, seq[i:i+4]))  # 1-indexed
    # Also check broader pattern: [RK]-X-X-R
    for i in range(len(seq) - 3):
        if seq[i] in ('R', 'K') and seq[i+3] == 'R':
            motif = seq[i:i+4]
            if (i + 1, motif) not in sites and seq[i+1] in ('R', 'K') or seq[i+2] in ('R', 'K'):
                pass  # Only report the strict pattern
    return sites


def net_charge(seq, pH=7.4):
    """Estimate net charge at given pH."""
    charge = 0
    for aa in seq:
        charge += AA_CHARGE_PH7.get(aa, 0)
    # N-terminus +1 (if pH < ~8), C-terminus -1 (if pH > ~3)
    charge += 0.5  # approximate at pH 7.4
    return charge


def gravy(seq):
    """Grand Average of Hydropathicity."""
    if not seq:
        return 0
    return sum(AA_HYDROPATHY.get(aa, 0) for aa in seq) / len(seq)


def find_hydrophobic_patches(seq, min_len=5):
    """Find runs of consecutive hydrophobic residues."""
    patches = []
    current_start = None
    current_len = 0
    for i, aa in enumerate(seq):
        if aa in HYDROPHOBIC:
            if current_start is None:
                current_start = i
            current_len += 1
        else:
            if current_len >= min_len:
                patches.append((current_start + 1, current_len, seq[current_start:current_start+current_len]))
            current_start = None
            current_len = 0
    if current_len >= min_len:
        patches.append((current_start + 1, current_len, seq[current_start:current_start+current_len]))
    return patches


def find_low_complexity(seq, min_run=5):
    """Find homopolymeric runs and low-complexity regions."""
    runs = []
    for aa in set(seq):
        for m in re.finditer(f'{aa}{{{min_run},}}', seq):
            runs.append((m.start() + 1, len(m.group()), m.group()))
    # Also check di-residue repeats (e.g., GSGSGS)
    for i in range(len(seq) - 5):
        window = seq[i:i+6]
        if len(set(window)) <= 2 and len(window) >= 6:
            # Check it's actually repetitive
            di = window[:2]
            if di * 3 == window:
                runs.append((i + 1, 6, f"({di})x3"))
    return runs


def cysteine_analysis(seq):
    """Count cysteines and check if they're likely paired."""
    n_cys = seq.count('C')
    paired = n_cys % 2 == 0
    return n_cys, paired


def instability_heuristic(seq):
    """Simplified instability estimate based on composition.
    Proteins with instability index > 40 are considered unstable.
    Uses the Guruprasad approach but simplified."""
    # Count destabilizing dipeptides (simplified)
    destabilizing = {'DG', 'GD', 'DP', 'WD', 'WG'}
    stabilizing = {'WW', 'WF', 'FW', 'YY', 'YF', 'FY'}

    score = 0
    for i in range(len(seq) - 1):
        di = seq[i:i+2]
        if di in destabilizing:
            score += 1
        elif di in stabilizing:
            score -= 0.5

    # Normalize per 100 residues
    return score / len(seq) * 100 if seq else 0


def check_expression_signals(seq):
    """Check for problematic expression signals."""
    issues = []

    # Check N-terminal Met
    if seq[0] != 'M':
        issues.append("No N-terminal Met (add for expression)")

    # Check for internal stop-like sequences (rare, but check)
    # Check for PEST sequences (rapid degradation signals)
    pest_pattern = re.compile(r'[PEST]{12,}', re.IGNORECASE)
    for m in pest_pattern.finditer(seq):
        issues.append(f"Potential PEST motif at {m.start()+1}-{m.end()}: {m.group()[:20]}")

    return issues


def domain_at_position(pos, binder_len):
    """Which domain is position in?"""
    boundaries = [
        (1, binder_len, "Binder"),
        (binder_len + 1, binder_len + 6, "GGSGGS spacer"),
        (binder_len + 7, binder_len + 6 + 121, "LOV2_lit"),
        (binder_len + 6 + 122, binder_len + 6 + 121 + 22, "E2pppE2"),
        (binder_len + 6 + 121 + 23, binder_len + 226, "UBQ"),
    ]
    for start, end, name in boundaries:
        if start <= pos <= end:
            return name
    return "unknown"


def main():
    with open(BASE / "final_top10_ranking.json") as f:
        data = json.load(f)

    top10 = data["top10"]

    print("=" * 130)
    print("EXPRESSION FEASIBILITY FILTER — TOP 10 CONSTITUTIVE CONSTRUCTS")
    print("=" * 130)
    print()

    all_issues = []

    for entry in top10:
        rank = entry["rank"]
        name = entry["binder"]
        binder_seq = entry["sequence"]
        binder_len = entry["length"]

        # Build full construct
        full_seq = binder_seq + SPACER + LOV2_LIT + E2PPPE2 + UBQ
        construct_len = len(full_seq)

        print(f"{'='*90}")
        print(f"  #{rank} {name} ({binder_len}aa binder, {construct_len}aa total)")
        print(f"  Score: {entry['integrated_score']:.3f} | Cap-trunc: {entry['cap_trunc_driver']}")
        print(f"{'='*90}")

        issues = []
        warnings = []

        # 1. Furin sites
        furin = find_furin_sites(full_seq)
        if furin:
            for pos, motif in furin:
                domain = domain_at_position(pos, binder_len)
                issues.append(f"FURIN SITE at pos {pos} ({motif}) in {domain}")

        # 2. Net charge
        charge = net_charge(full_seq)
        binder_charge = net_charge(binder_seq)
        if abs(charge) > 15:
            warnings.append(f"High net charge: {charge:+.1f} (binder: {binder_charge:+.1f})")

        # 3. GRAVY
        g = gravy(full_seq)
        binder_g = gravy(binder_seq)
        if g > 0:
            warnings.append(f"Positive GRAVY ({g:.2f}) — aggregation risk")
        if binder_g > 0.5:
            issues.append(f"Binder hydrophobic (GRAVY={binder_g:.2f}) — likely aggregation-prone")

        # 4. Hydrophobic patches
        patches = find_hydrophobic_patches(full_seq, min_len=5)
        for pos, length, patch in patches:
            domain = domain_at_position(pos, binder_len)
            if domain == "Binder":
                issues.append(f"Hydrophobic patch in BINDER at pos {pos}: {patch} ({length}aa)")
            else:
                warnings.append(f"Hydrophobic patch in {domain} at pos {pos}: {patch} ({length}aa)")

        # 5. Low-complexity regions
        lc = find_low_complexity(binder_seq, min_run=4)
        for pos, length, run in lc:
            if length >= 6:
                issues.append(f"Low-complexity run in binder at pos {pos}: {run}")
            else:
                warnings.append(f"Short repeat in binder at pos {pos}: {run}")

        # Also check for Gly/Ser-rich binders (disorder prone)
        gs_frac = (binder_seq.count('G') + binder_seq.count('S')) / len(binder_seq)
        if gs_frac > 0.5:
            warnings.append(f"Binder is {gs_frac:.0%} Gly+Ser — likely disordered, may not fold independently")

        # Ala-rich (helical but potentially low stability)
        ala_frac = binder_seq.count('A') / len(binder_seq)
        if ala_frac > 0.35:
            warnings.append(f"Binder is {ala_frac:.0%} Ala — may form marginal helix with low stability")

        # 6. Cysteine
        n_cys, paired = cysteine_analysis(full_seq)
        binder_cys = binder_seq.count('C')
        scaffold_cys = full_seq.count('C') - binder_cys
        if binder_cys > 0:
            if binder_cys % 2 != 0:
                issues.append(f"Odd number of Cys in binder ({binder_cys}) — unpaired Cys causes misfolding in oxidizing env")
            else:
                warnings.append(f"Binder has {binder_cys} Cys — ensure disulfide bonds form correctly")
        if scaffold_cys > 0:
            warnings.append(f"Scaffold has {scaffold_cys} Cys (from LOV2/UBQ) — normal for these domains")

        # 7. N-terminal Met
        if full_seq[0] != 'M':
            warnings.append(f"No N-terminal Met (starts with {full_seq[0]}). Add Met for expression.")

        # 8. Instability heuristic
        inst = instability_heuristic(full_seq)
        if inst > 5:
            warnings.append(f"Elevated instability score ({inst:.1f})")

        # 9. Proline at junctions (can disrupt folding)
        spacer_start = binder_len
        if binder_seq[-1] == 'P' or binder_seq[-2:] == 'PP':
            warnings.append(f"Pro at binder C-terminus ({binder_seq[-3:]}...) — may affect spacer flexibility")

        # 10. Consecutive rare codons (for E. coli expression)
        rare_in_ecoli = set('PWCM')  # These have some rare codons
        rare_runs = []
        count = 0
        for i, aa in enumerate(binder_seq):
            if aa in rare_in_ecoli:
                count += 1
            else:
                if count >= 3:
                    rare_runs.append((i - count + 1, count))
                count = 0
        if rare_runs:
            for pos, length in rare_runs:
                warnings.append(f"Potential rare codon cluster at binder pos {pos}: {binder_seq[pos-1:pos-1+length]}")

        # Print results
        if issues:
            print(f"\n  ISSUES ({len(issues)}):")
            for issue in issues:
                print(f"    [!] {issue}")

        if warnings:
            print(f"\n  WARNINGS ({len(warnings)}):")
            for w in warnings:
                print(f"    [~] {w}")

        if not issues and not warnings:
            print("\n  [OK] No issues detected")

        # Biophysical summary
        print(f"\n  Biophysics: charge={charge:+.1f} | GRAVY={g:.2f} | "
              f"binder_GRAVY={binder_g:.2f} | Cys={n_cys} | "
              f"GS%={gs_frac:.0%} | Ala%={ala_frac:.0%}")

        # Verdict
        if any("FURIN" in i for i in issues):
            verdict = "REJECT (furin site)"
        elif any("aggregation-prone" in i for i in issues):
            verdict = "CAUTION (aggregation)"
        elif any("Odd number of Cys" in i for i in issues):
            verdict = "CAUTION (unpaired Cys)"
        elif gs_frac > 0.6:
            verdict = "CAUTION (disordered binder)"
        elif ala_frac > 0.45:
            verdict = "CAUTION (low-complexity binder)"
        elif len(issues) > 0:
            verdict = "REVIEW"
        else:
            verdict = "OK for expression"

        print(f"  VERDICT: {verdict}")
        print()

        all_issues.append({
            "rank": rank,
            "binder": name,
            "binder_length": binder_len,
            "construct_length": construct_len,
            "n_issues": len(issues),
            "n_warnings": len(warnings),
            "issues": issues,
            "warnings": warnings,
            "verdict": verdict,
            "charge": float(charge),
            "gravy": float(g),
            "binder_gravy": float(binder_g),
            "gs_fraction": float(gs_frac),
            "ala_fraction": float(ala_frac),
            "n_cys_total": n_cys,
            "n_cys_binder": binder_cys,
            "has_nterm_met": full_seq[0] == 'M',
            "furin_sites": [(pos, motif) for pos, motif in furin] if furin else [],
        })

    # Summary
    print("=" * 130)
    print("EXPRESSION FILTER SUMMARY")
    print("=" * 130)
    print()
    print(f"{'Rk':<4} {'Binder':<10} {'Len':<5} {'Charge':<8} {'GRAVY':<7} "
          f"{'BdrGRAVY':<9} {'GS%':<6} {'Ala%':<6} {'Cys':<5} {'Met1':<5} "
          f"{'Issues':<7} {'Verdict'}")
    print("-" * 110)

    for r in all_issues:
        met1 = "Y" if r["has_nterm_met"] else "N"
        print(f"{r['rank']:<4} {r['binder']:<10} {r['construct_length']:>3}  "
              f"{r['charge']:>+6.1f}  {r['gravy']:>5.2f}  {r['binder_gravy']:>7.2f}  "
              f"{r['gs_fraction']:>4.0%}  {r['ala_fraction']:>4.0%}  "
              f"{r['n_cys_total']:>3}   {met1:<5} "
              f"{r['n_issues']:<7} {r['verdict']}")

    # Recommendations
    print()
    print("=" * 130)
    print("EXPRESSION RECOMMENDATIONS")
    print("=" * 130)

    rejects = [r for r in all_issues if "REJECT" in r["verdict"]]
    cautions = [r for r in all_issues if "CAUTION" in r["verdict"]]
    ok = [r for r in all_issues if r["verdict"] == "OK for expression"]

    if rejects:
        print(f"\n  REJECT ({len(rejects)}):")
        for r in rejects:
            print(f"    {r['binder']}: {'; '.join(r['issues'][:3])}")

    if cautions:
        print(f"\n  CAUTION ({len(cautions)}):")
        for r in cautions:
            top_issue = r["issues"][0] if r["issues"] else r["warnings"][0] if r["warnings"] else "see details"
            print(f"    {r['binder']}: {top_issue}")

    if ok:
        print(f"\n  OK ({len(ok)}):")
        for r in ok:
            print(f"    {r['binder']}: No issues detected")

    # General advice
    no_met = [r for r in all_issues if not r["has_nterm_met"]]
    if no_met:
        print(f"\n  NOTE: {len(no_met)} constructs lack N-terminal Met — prepend M for expression:")
        for r in no_met:
            print(f"    {r['binder']} (starts with {SEQUENCES_SHORT.get(r['binder'], '?')[0]})")

    disordered = [r for r in all_issues if r["gs_fraction"] > 0.4]
    if disordered:
        print(f"\n  NOTE: {len(disordered)} binders have high Gly+Ser content (>40%) — likely disordered.")
        print("    These may benefit from SUMO or MBP fusion tags for solubility.")
        print("    Alternatively, cell-free expression may give better yields than CHO.")

    print(f"\n  GENERAL:")
    print(f"    - All constructs are 246-276aa — adequate size for CHO or cell-free")
    print(f"    - Scaffold (LOV2+E2pppE2+UBQ) has 2 Cys from LOV2 (C450 for FMN binding)")
    print(f"    - For BLI: N-terminal His-tag or C-terminal FLAG compatible with all")
    print(f"    - For CHO: add signal peptide if secretion required")
    print(f"    - For cell-free: no modifications needed beyond Met initiation")

    # Save
    out_path = BASE / "expression_filter_results.json"
    out_path.write_text(json.dumps(all_issues, indent=2) + "\n")
    print(f"\n  Results saved: {out_path}")


# Quick lookup for N-term residue reporting
SEQUENCES_SHORT = {
    "pep21s2": "SIGSSSSGGSGSLGSVSGSS",
    "pep21s5": "SIGSSSSGASGSLGAVTAED",
    "pep29s4": "MMKVVFGGGCGIGTPISEEDIRRGNIGGKEGGAESKARSAARAAAAAAAL",
    "pep16s4": "AGMPLGGRTEEELRTSRNAAEGDAENRARS",
    "pep11s1": "SSRPNPPRLDARIEESRAGL",
    "pep9s5": "MTPRKPSDAGLVRAAEKAAA",
    "pep37s3": "ETPADRDAAGGRAAAARLNL",
    "pep41s5": "DDVERGRAAGAAGTGARGEF",
    "pep19s2": "ATVADPESGGSAVTAGAPVE",
    "pep30s5": "APPGKTAAEDAARRAAAAAA",
}


if __name__ == "__main__":
    main()
