#!/usr/bin/env python3
"""
FINAL INTEGRATED RANKING — All binders, all layers, definitive top-10 selection.

Merges original binding validation pool (16 peptides) + expanded pool (11 advancing)
through the 5-layer scoring system:

  L1: Epitope contacts (HAT hotspot overlap, binder-alone)
  L2: Cap-truncated fusion (binder domain contacts WITH scaffold) — original pool only
  L3: Multi-sample consistency (reproducibility across 5 diffusion samples)
  L4: Interface quality (best ipTM, stability, pLDDT)
  L5: CatCore epitope agreement (hotspot contacts against full catalytic core)

Weights: L1=20%, L2=10%, L3=20%, L4=25%, L5=25%
"""

import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).parent
HOTSPOTS = {91, 110, 111, 150, 158, 159, 160, 163, 170, 181, 185, 219, 220, 221}

# Binder sequences for final output
SEQUENCES = {
    # Original pool
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
    # Expanded pool
    "pep6s1": "AAAGQASDKPEQAAEVAAALGAARDAAVAA",
    "pep11s1": "SSRPNPPRLDARIEESRAGL",
    "pep49s5": "DSLGSGCESGSTGRERLARCSEILKAQLEA",
    "pep9s5": "MTPRKPSDAGLVRAAEKAAA",
    "pep29s4": "MMKVVFGGGCGIGTPISEEDIRRGNIGGKEGGAESKARSAARAAAAAAAL",
    "pep4s4": "GCMPATSEEAARQYAACLAR",
    "pep30s4": "AVDAKQSKETGAEALAAAADAVSAAAAAAE",
    "pep41s5": "DDVERGRAAGAAGTGARGEF",
    "pep34s1": "SEDREKAERHQAALEEVVEATRAAYREDRELGGRVAREELGRLERAEPDD",
    "pep46s3": "AAGGGGPAKDEASRRNMESARRKARRVAAL",
    "pep37s3": "ETPADRDAAGGRAAAARLNL",
    "pep19s2": "ATVADPESGGSAVTAGAPVE",
    "pep16s4": "AGMPLGGRTEEELRTSRNAAEGDAENRARS",
    "pep40s3": "SMETFSSAEAQANYEKGLAGREAGDARNQARRDAARAAAA",
    "pep43s1": "AVASTPAAAKGGGYAEPGAP",
    "pep14s4": "PETEADKAMRARTAAAIAAA",
}


def load_original_pool():
    """Load original binding validation pool data."""
    binders = {}

    # Epitope contacts
    with open(BASE / "epitope_contact_results.json") as f:
        for r in json.load(f):
            binders.setdefault(r["binder"], {}).update({
                "L1_hotspot_pct": r["hotspot_pct"] / 100.0,
                "L1_hotspot_count": r["hotspot_contacts"],
                "L1_verdict": r["verdict"],
            })

    # Multi-sample consistency
    with open(BASE / "multisample_consistency_results.json") as f:
        for r in json.load(f):
            binders.setdefault(r["binder"], {}).update({
                "L3_consistency": r["consistency"],
                "L3_n_consistent": r["n_consistent_hotspots"],
                "L3_mean_hotspots": r["mean_hotspots"],
            })

    # Interface quality
    with open(BASE / "interface_quality_results.json") as f:
        for r in json.load(f):
            binders.setdefault(r["binder"], {}).update({
                "L4_best_iptm": r["best_iptm"],
                "L4_mean_iptm": r["mean_iptm"],
                "L4_std_iptm": r["std_iptm"],
                "L4_iface_plddt": r["mean_interface_plddt"],
            })

    # Cap-truncated fusion
    ct_path = BASE / "captrunc_validation" / "captrunc_results.json"
    if ct_path.exists():
        with open(ct_path) as f:
            for r in json.load(f):
                binders.setdefault(r["binder"], {}).update({
                    "L2_driver": r["driver"],
                    "L2_bdr_hotspots": r["mean_bdr_hotspots"],
                })

    # CatCore epitope
    cc_path = BASE / "catcore_epitope_results.json"
    if cc_path.exists():
        with open(cc_path) as f:
            for r in json.load(f):
                binders.setdefault(r["binder"], {}).update({
                    "L5_cc_max_hs": r["max_hotspot"],
                    "L5_cc_avg_hs": r["avg_hotspot"],
                    "L5_cc_consistent": r["n_consistent"],
                    "L5_cc_iptm": r["best_iptm"],
                })

    # Binding validation ipTMs
    with open(BASE / "binding_validation" / "analysis" / "binding_validation_results.json") as f:
        for r in json.load(f)["raw_data"]:
            binders.setdefault(r["binder"], {}).update({
                "screen_iptm": r["s1_iptm"],
                "hat_iptm": r["binder_HAT_iptm"],
                "cc_iptm": r["binder_catcore_iptm"],
                "length": r["length"],
            })

    return binders


def load_expanded_pool():
    """Load expanded pool data."""
    binders = {}

    ep_path = BASE / "expanded_pool" / "expanded_pool_results.json"
    if not ep_path.exists():
        return binders

    with open(ep_path) as f:
        for r in json.load(f):
            name = r["binder"]
            binders[name] = {
                "screen_iptm": r["screen_iptm"],
                "length": r["binder_length"],
                "L1_hotspot_pct": r["best_hotspots"] / 14.0,
                "L1_hotspot_count": r["best_hotspots"],
                "L1_verdict": r["verdict"],
                "L3_consistency": r["consistency"],
                "L3_n_consistent": r["n_consistent_hotspots"],
                "L3_mean_hotspots": r["mean_hotspots"],
                "L4_best_iptm": r["best_iptm"],
                "L4_mean_iptm": r["mean_iptm"],
                "L4_std_iptm": r["std_iptm"],
                "L4_iface_plddt": r["mean_binder_plddt"],
                "hat_iptm": r["best_iptm"],
                # No cap-truncated data for expanded pool
                "L2_driver": "NOT_TESTED",
                "L2_bdr_hotspots": 0,
            }

    # Add catcore results
    cc_path = BASE / "expanded_pool" / "expanded_catcore_results.json"
    if cc_path.exists():
        with open(cc_path) as f:
            for r in json.load(f):
                name = r["binder"]
                if name in binders:
                    binders[name].update({
                        "L5_cc_max_hs": r["best_cc_hotspots"],
                        "L5_cc_avg_hs": r["mean_cc_hotspots"],
                        "L5_cc_consistent": r["cc_n_consistent"],
                        "L5_cc_iptm": r["best_cc_iptm"],
                        "cc_iptm": r["best_cc_iptm"],
                    })

    # Add cap-truncated fusion results (HAT + CatCore)
    ct_path = BASE / "expanded_pool" / "expanded_captrunc_results.json"
    if ct_path.exists():
        with open(ct_path) as f:
            ct_data = json.load(f)

        # Group by binder, take best driver across HAT and CatCore
        driver_rank = {"BINDER-DRIVEN": 4, "MIXED": 3, "SCAFFOLD-DOM": 2, "SCAFFOLD-ONLY": 1}
        ct_by_binder = {}
        for r in ct_data:
            name = r["binder"]
            existing = ct_by_binder.get(name)
            if existing is None or driver_rank.get(r["driver"], 0) > driver_rank.get(existing["driver"], 0):
                ct_by_binder[name] = r

        for name, r in ct_by_binder.items():
            if name in binders:
                binders[name]["L2_driver"] = r["driver"]
                binders[name]["L2_bdr_hotspots"] = r["mean_bdr_hotspots"]

    return binders


def score_binder(data):
    """Compute 5-layer integrated score."""

    # L1: Epitope (0-1)
    l1 = data.get("L1_hotspot_pct", 0)

    # L2: Cap-truncated fusion
    driver = data.get("L2_driver", "NOT_TESTED")
    if driver == "BINDER-DRIVEN":
        l2 = 1.0
    elif driver == "MIXED":
        l2 = 0.7
    elif driver == "SCAFFOLD-DOM":
        l2 = 0.3
    elif driver == "SCAFFOLD-ONLY":
        l2 = 0.1
    elif driver == "NOT_TESTED":
        l2 = 0.5  # Neutral — no data yet
    else:
        l2 = 0.2

    # L3: Consistency
    consistency = data.get("L3_consistency", 0)
    n_consistent = data.get("L3_n_consistent", 0)
    l3 = (n_consistent / 14.0) * consistency

    # L4: Interface quality
    best_iptm = data.get("L4_best_iptm", 0)
    std_iptm = data.get("L4_std_iptm", 1.0)
    iface_plddt = data.get("L4_iface_plddt", 0)

    stability = max(0, min(1, 1 - std_iptm / 0.3))
    # Handle pLDDT on either 0-1 or 0-100 scale
    plddt_norm = min(1.0, iface_plddt / 100.0) if iface_plddt > 1 else iface_plddt
    l4 = 0.5 * best_iptm + 0.25 * stability + 0.25 * plddt_norm

    # L5: CatCore agreement
    cc_consistent = data.get("L5_cc_consistent", 0)
    cc_iptm = data.get("L5_cc_iptm", 0)
    l5 = min(1.0, cc_consistent / 8.0)
    if cc_iptm >= 0.80:
        l5 = min(1.0, l5 + 0.2)
    elif cc_iptm < 0.40:
        l5 *= 0.5

    return {
        "L1": l1, "L2": l2, "L3": l3, "L4": l4, "L5": l5,
        "integrated": 0.20 * l1 + 0.10 * l2 + 0.20 * l3 + 0.25 * l4 + 0.25 * l5,
    }


def main():
    original = load_original_pool()
    expanded = load_expanded_pool()

    # Merge — expanded pool binders that passed HAT + catcore
    all_binders = {}
    for name, data in original.items():
        if name.startswith("nb"):
            continue  # Skip nanobodies (all off-target)
        all_binders[name] = data
        all_binders[name]["pool"] = "original"

    for name, data in expanded.items():
        if data.get("L1_verdict") in ("OFF-TARGET",):
            continue
        if data.get("L3_consistency", 0) < 0.60:
            continue
        if data.get("L5_cc_consistent") is None:
            continue  # No catcore data
        all_binders[name] = data
        all_binders[name]["pool"] = "expanded"

    # Score all
    results = []
    for name, data in all_binders.items():
        scores = score_binder(data)
        results.append({
            "binder": name,
            "pool": data.get("pool", "?"),
            "length": data.get("length", len(SEQUENCES.get(name, ""))),
            "sequence": SEQUENCES.get(name, "???"),
            "screen_iptm": data.get("screen_iptm", 0),
            "hat_iptm": data.get("hat_iptm", 0),
            "cc_iptm": data.get("cc_iptm", 0),
            "L1_verdict": data.get("L1_verdict", "?"),
            "L2_driver": data.get("L2_driver", "?"),
            "L3_n_consistent": data.get("L3_n_consistent", 0),
            "L5_cc_consistent": data.get("L5_cc_consistent", 0),
            "L5_cc_max_hs": data.get("L5_cc_max_hs", 0),
            **scores,
        })

    results.sort(key=lambda r: r["integrated"], reverse=True)

    print("=" * 150)
    print("DEFINITIVE RANKING — ALL BINDERS, ALL LAYERS")
    print("Target: p300 (HAT + CatCore) | Construct: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]-[UBQ]")
    print("Weights: L1(epitope)=20% | L2(cap-trunc)=10% | L3(consistency)=20% | L4(quality)=25% | L5(catcore)=25%")
    print("=" * 150)
    print()

    hdr = (f"{'Rk':<4} {'Binder':<10} {'Pool':<6} {'Len':<5} "
           f"{'HAT HS':<8} {'L2:Driver':<14} {'ConsHS':<7} "
           f"{'HAT ipTM':<9} {'CC ipTM':<9} {'CC HS':<8} "
           f"{'SCORE':<8} {'TIER'}")
    print(hdr)
    print("-" * 150)

    for i, r in enumerate(results, 1):
        tier = ""
        if r["integrated"] >= 0.50:
            tier = "*** TIER 1 ***"
        elif r["integrated"] >= 0.35:
            tier = "TIER 2"
        elif r["integrated"] >= 0.20:
            tier = "tier 3"
        else:
            tier = "exclude"

        pool = "orig" if r["pool"] == "original" else "NEW"
        cc_hs = f"{r['L5_cc_consistent']}/{r['L5_cc_max_hs']}"

        print(f"{i:<4} {r['binder']:<10} {pool:<6} {r['length']:>3}  "
              f"{r['L1_verdict']:<8} {r['L2_driver']:<14} {r['L3_n_consistent']:>3}    "
              f"{r['hat_iptm']:>7.3f}  {r['cc_iptm']:>7.3f}  {cc_hs:<8} "
              f"{r['integrated']:>6.3f}  {tier}")

    # Top 10 selection
    tier1 = [r for r in results if r["integrated"] >= 0.50]
    tier2 = [r for r in results if 0.35 <= r["integrated"] < 0.50]
    tier3 = [r for r in results if 0.20 <= r["integrated"] < 0.35]

    print()
    print("=" * 150)
    print(f"TIER 1 — HIGH CONFIDENCE: {len(tier1)} binders")
    print("=" * 150)
    for r in tier1:
        pool_tag = " [NEW]" if r["pool"] == "expanded" else ""
        print(f"  {r['binder']:<10} {r['length']:>3}aa  score={r['integrated']:.3f}  "
              f"HAT={r['hat_iptm']:.3f}  CC={r['cc_iptm']:.3f}  "
              f"HAT_HS={r['L1_verdict']}  CC_consHS={r['L5_cc_consistent']}  "
              f"consistency={r['L3_n_consistent']}{pool_tag}")

    print()
    print(f"TIER 2 — MODERATE CONFIDENCE: {len(tier2)} binders")
    for r in tier2:
        pool_tag = " [NEW]" if r["pool"] == "expanded" else ""
        print(f"  {r['binder']:<10} {r['length']:>3}aa  score={r['integrated']:.3f}  "
              f"HAT={r['hat_iptm']:.3f}  CC={r['cc_iptm']:.3f}{pool_tag}")

    # Final top 10 recommendation
    top10 = results[:10]
    print()
    print("=" * 150)
    print("RECOMMENDED TOP 10 FOR EMPIRICAL TESTING")
    print("=" * 150)
    print()

    scaffold_len = 6 + 121 + 23 + 76  # GGSGGS + LOV2_lit + E2pppE2 + UBQ = 226
    for i, r in enumerate(top10, 1):
        construct_len = r["length"] + scaffold_len
        pool_tag = " [NEW from expanded pool]" if r["pool"] == "expanded" else ""
        print(f"  {i:>2}. {r['binder']:<10} ({r['length']:>2}aa binder, {construct_len}aa construct){pool_tag}")
        print(f"      Score: {r['integrated']:.3f}  |  HAT ipTM: {r['hat_iptm']:.3f}  |  CatCore ipTM: {r['cc_iptm']:.3f}")
        print(f"      HAT hotspots: {r['L1_verdict']}  |  CatCore consistent HS: {r['L5_cc_consistent']}/{r['L5_cc_max_hs']}")
        print(f"      Consistency: {r['L3_n_consistent']} hotspots across 5 models  |  Cap-trunc: {r['L2_driver']}")
        print(f"      Sequence: {r['sequence']}")
        print()

    # Full construct sequences
    print("=" * 150)
    print("FULL CONSTRUCT SEQUENCES (ready for synthesis)")
    print("[Binder]-GGSGGS-LOV2_lit(121aa)-EAAAKEAAAKPPPEAAAKEAAAK-UBQ(76aa)")
    print("=" * 150)

    LOV2_LIT = ("LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
                "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAEREG")
    E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"
    UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    SPACER = "GGSGGS"

    for i, r in enumerate(top10, 1):
        full = r["sequence"] + SPACER + LOV2_LIT + E2PPPE2 + UBQ
        print(f"\n  {i}. {r['binder']} ({len(full)}aa):")
        # Print in 80-char lines
        for j in range(0, len(full), 80):
            print(f"     {full[j:j+80]}")

    # Save
    out = {
        "description": "Definitive top-10 binder selection for empirical testing",
        "target": "p300 HAT domain (tested against both isolated HAT and full catalytic core)",
        "construct": "[Binder]-[GGSGGS(6aa)]-[LOV2_lit(121aa)]-[E2pppE2(23aa)]-[UBQ(76aa)]",
        "scaffold_length": scaffold_len,
        "total_binders_evaluated": len(results),
        "tier1_count": len(tier1),
        "tier2_count": len(tier2),
        "scoring": {
            "L1_weight": 0.20, "L2_weight": 0.10, "L3_weight": 0.20,
            "L4_weight": 0.25, "L5_weight": 0.25,
        },
        "top10": [
            {
                "rank": i + 1,
                "binder": r["binder"],
                "pool": r["pool"],
                "length": r["length"],
                "construct_length": r["length"] + scaffold_len,
                "sequence": r["sequence"],
                "integrated_score": r["integrated"],
                "hat_iptm": r["hat_iptm"],
                "catcore_iptm": r["cc_iptm"],
                "hat_verdict": r["L1_verdict"],
                "catcore_consistent_hs": r["L5_cc_consistent"],
                "consistency_hotspots": r["L3_n_consistent"],
                "cap_trunc_driver": r["L2_driver"],
            }
            for i, r in enumerate(top10)
        ],
        "all_rankings": [
            {
                "rank": i + 1,
                "binder": r["binder"],
                "pool": r["pool"],
                "length": r["length"],
                "integrated_score": r["integrated"],
                "hat_iptm": r["hat_iptm"],
                "catcore_iptm": r["cc_iptm"],
                "L1": r["L1"], "L2": r["L2"], "L3": r["L3"],
                "L4": r["L4"], "L5": r["L5"],
            }
            for i, r in enumerate(results)
        ],
    }
    out_path = BASE / "final_top10_ranking.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n\nFinal ranking saved: {out_path}")


if __name__ == "__main__":
    main()
