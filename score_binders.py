#!/usr/bin/env python3
"""
Score Boltz-2 binder predictions and select top N for opto-scaffold campaign.

Parses Boltz-2 confidence outputs (ipTM, pTM, pLDDT), counts hotspot contacts,
ranks binders, and auto-generates a campaign config for setup_campaign.py.

Usage:
    python3 score_binders.py campaigns/<name>/binder_screen/ --top_n 10
    python3 score_binders.py campaigns/<name>/binder_screen/ --top_n 10 --config campaign_config.json
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

try:
    from Bio.PDB import MMCIFParser, PDBParser
except ImportError:
    print("ERROR: BioPython required. Install: pip install biopython")
    sys.exit(1)

from linker_library import LOV2_DARK, LOV2_LIT, UBQ_SEQ

warnings.filterwarnings("ignore")
PIPELINE_DIR = Path(__file__).parent


def extract_confidence(result_dir):
    """Extract confidence metrics from Boltz-2 prediction results.

    Searches for confidence JSON (best model = model_0) and PAE .npz files.

    Returns dict with ipTM, pTM, pLDDT or None if no results found.
    """
    result_dir = Path(result_dir)

    # Find confidence JSON (best model = model_0)
    conf_files = sorted(result_dir.rglob("confidence_*_model_0.json"))
    if not conf_files:
        conf_files = sorted(result_dir.rglob("confidence_*.json"))
    if not conf_files:
        return None

    with open(conf_files[0]) as f:
        conf = json.load(f)

    return {
        "iptm": conf.get("protein_iptm", conf.get("iptm", 0.0)),
        "ptm": conf.get("ptm", 0.0),
        "plddt": conf.get("complex_plddt", 0.0),
        "confidence": conf.get("confidence_score", 0.0),
    }


def extract_hotspot_contacts(cif_path, binder_chain, target_chain,
                             hotspot_residues_1idx, contact_dist=5.0):
    """Count binder atom contacts to target hotspot residues.

    Args:
        cif_path: path to prediction CIF file
        binder_chain: chain ID of binder (usually 'B')
        target_chain: chain ID of target (usually 'A')
        hotspot_residues_1idx: list of 1-indexed residue numbers on target
        contact_dist: distance threshold in Angstroms

    Returns dict with n_hotspot_contacts, hotspot_residues_contacted.
    """
    cif_path = Path(cif_path)
    if not cif_path.exists():
        return None

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("pred", str(cif_path))
    model = structure[0]

    chains = {c.id: c for c in model.get_chains()}
    if binder_chain not in chains or target_chain not in chains:
        return None

    hotspot_set = set(hotspot_residues_1idx)

    # Get target hotspot atoms
    target_hotspot_atoms = []
    for res in chains[target_chain].get_residues():
        res_num = res.get_id()[1]
        if res_num in hotspot_set:
            for atom in res.get_atoms():
                target_hotspot_atoms.append((res_num, atom.get_vector().get_array()))

    if not target_hotspot_atoms:
        return {"n_hotspot_contacts": 0, "hotspot_residues_contacted": []}

    # Get binder atoms
    binder_atoms = []
    for res in chains[binder_chain].get_residues():
        for atom in res.get_atoms():
            binder_atoms.append(atom.get_vector().get_array())

    if not binder_atoms:
        return {"n_hotspot_contacts": 0, "hotspot_residues_contacted": []}

    binder_coords = np.array(binder_atoms)
    contacted_residues = set()
    n_contacts = 0

    for res_num, target_coord in target_hotspot_atoms:
        dists = np.linalg.norm(binder_coords - target_coord, axis=1)
        if np.min(dists) < contact_dist:
            contacted_residues.add(res_num)
            n_contacts += int(np.sum(dists < contact_dist))

    return {
        "n_hotspot_contacts": n_contacts,
        "hotspot_residues_contacted": sorted(contacted_residues),
    }


def parse_binder_fasta(fasta_path):
    """Parse FASTA file to extract binder name → sequence mapping."""
    binders = {}
    current_name = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    binders[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            elif line:
                # Handle LigandMPNN format: chain_A:chain_B
                if ":" in line:
                    line = line.split(":")[-1]
                current_seq.append(line)
        if current_name:
            binders[current_name] = "".join(current_seq)

    return binders


def score_binder_screen(screen_dir, hotspot_residues_1idx=None, top_n=10):
    """Score all Boltz-2 binder predictions in a screen directory.

    Args:
        screen_dir: path to binder_screen/ directory containing boltz2_results/
        hotspot_residues_1idx: target hotspot residues for contact analysis
        top_n: number of top binders to select

    Returns (ranked_results, top_binders) where top_binders is list of dicts
    with name, sequence, iptm, etc.
    """
    screen_dir = Path(screen_dir)
    results_dir = screen_dir / "boltz2_results"

    if not results_dir.exists():
        print(f"ERROR: {results_dir} not found")
        return [], []

    # Find all prediction directories
    pred_dirs = sorted(d for d in results_dir.iterdir() if d.is_dir())
    if not pred_dirs:
        print(f"ERROR: No prediction results in {results_dir}")
        return [], []

    # Try to load binder sequences from FASTA
    binder_seqs = {}
    fasta_files = list(screen_dir.glob("*.fa")) + list(screen_dir.glob("*.fasta"))
    for fa in fasta_files:
        binder_seqs.update(parse_binder_fasta(fa))

    print(f"Scoring {len(pred_dirs)} binder predictions...")
    results = []

    for pred_dir in pred_dirs:
        name = pred_dir.name
        metrics = extract_confidence(pred_dir)

        if metrics is None:
            print(f"  {name}: no confidence data, skipping")
            continue

        result = {
            "name": name,
            "iptm": metrics["iptm"],
            "ptm": metrics["ptm"],
            "plddt": metrics["plddt"],
            "confidence": metrics["confidence"],
            "sequence": binder_seqs.get(name, ""),
        }

        # Hotspot contact analysis if available
        if hotspot_residues_1idx:
            cif_files = sorted(pred_dir.rglob("*_model_0.cif"))
            if not cif_files:
                cif_files = sorted(pred_dir.rglob("*.cif"))
            if cif_files:
                contacts = extract_hotspot_contacts(
                    cif_files[0], "B", "A", hotspot_residues_1idx
                )
                if contacts:
                    result.update(contacts)

        results.append(result)
        iptm_str = f"ipTM={metrics['iptm']:.3f}"
        hs_str = ""
        if "n_hotspot_contacts" in result:
            hs_str = f" HS={len(result.get('hotspot_residues_contacted', []))}"
        print(f"  {name}: {iptm_str}{hs_str}")

    # Rank by ipTM
    results.sort(key=lambda x: x["iptm"], reverse=True)

    # Select top N
    top = results[:top_n]

    print(f"\nTop {len(top)} binders by ipTM:")
    for i, r in enumerate(top, 1):
        seq_len = len(r["sequence"]) if r["sequence"] else "?"
        hs = len(r.get("hotspot_residues_contacted", []))
        print(f"  {i:>2}. {r['name']:<40} ipTM={r['iptm']:.3f}  "
              f"{seq_len}aa  {hs} hotspots")

    return results, top


def generate_opto_config(top_binders, base_config_path, output_path):
    """Auto-generate campaign config for setup_campaign.py with top binders.

    Args:
        top_binders: list of dicts with name, sequence, iptm
        base_config_path: path to base campaign config (for target, scaffold settings)
        output_path: where to write the new config
    """
    with open(base_config_path) as f:
        base_cfg = json.load(f)

    # Update campaign name
    base_name = base_cfg.get("campaign_name", "campaign")
    base_cfg["campaign_name"] = f"{base_name}_opto"

    # Replace binders with top scored binders
    base_cfg["binders"] = []
    for b in top_binders:
        if not b.get("sequence"):
            print(f"  WARNING: {b['name']} has no sequence, skipping")
            continue
        base_cfg["binders"].append({
            "name": b["name"],
            "sequence": b["sequence"],
            "length": len(b["sequence"]),
            "iptm": round(b["iptm"], 4),
        })

    # Ensure scaffold defaults are present
    scaffold = base_cfg.setdefault("scaffold", {})
    scaffold.setdefault("architecture", "C")
    scaffold.setdefault("lov2_dark_seq", LOV2_DARK)
    scaffold.setdefault("lov2_lit_seq", LOV2_LIT)
    scaffold.setdefault("cap_seq", UBQ_SEQ)
    scaffold.setdefault("cap_name", "UBQ")
    scaffold.setdefault("binder_lov2_spacer", "GGSGGS")
    scaffold.setdefault("linkers", "all")

    output_path = Path(output_path)
    output_path.write_text(json.dumps(base_cfg, indent=2) + "\n")
    print(f"\nOpto campaign config saved: {output_path}")
    print(f"  {len(base_cfg['binders'])} binders × linkers → setup_campaign.py")
    return base_cfg


def main():
    parser = argparse.ArgumentParser(
        description="Score Boltz-2 binder predictions and select top N"
    )
    parser.add_argument("screen_dir",
                        help="Path to binder_screen/ directory with boltz2_results/")
    parser.add_argument("--top_n", type=int, default=10,
                        help="Number of top binders to select (default: 10)")
    parser.add_argument("--config", default=None,
                        help="Base campaign config for auto-generating opto config")
    parser.add_argument("--hotspots", default=None,
                        help="Comma-separated 1-indexed hotspot residues (e.g., 91,110,111)")
    args = parser.parse_args()

    screen_dir = Path(args.screen_dir)
    if not screen_dir.exists():
        print(f"ERROR: {screen_dir} not found")
        sys.exit(1)

    # Parse hotspots
    hotspot_residues = None
    if args.hotspots:
        hotspot_residues = [int(x.strip()) for x in args.hotspots.split(",")]

    # If config provided, try to get hotspots from there
    if hotspot_residues is None and args.config:
        with open(args.config) as f:
            cfg = json.load(f)
        hotspot_residues = cfg.get("target", {}).get("hotspot_residues_1idx")

    # Score
    all_results, top_binders = score_binder_screen(
        screen_dir, hotspot_residues, args.top_n
    )

    if not all_results:
        print("No results to score.")
        sys.exit(1)

    # Save ranked results
    ranked_path = screen_dir / "ranked_binders.json"
    ranked_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nAll ranked binders saved: {ranked_path}")

    # Save top binders FASTA
    fasta_lines = []
    for b in top_binders:
        if b.get("sequence"):
            fasta_lines.append(f">{b['name']} | ipTM={b['iptm']:.3f} | {len(b['sequence'])}aa")
            for i in range(0, len(b["sequence"]), 80):
                fasta_lines.append(b["sequence"][i:i+80])
    if fasta_lines:
        fasta_path = screen_dir / "top_binders.fasta"
        fasta_path.write_text("\n".join(fasta_lines) + "\n")
        print(f"Top binders FASTA: {fasta_path}")

    # Auto-generate opto config if base config provided
    if args.config:
        opto_config_path = screen_dir.parent / "campaign_config_opto.json"
        generate_opto_config(top_binders, args.config, opto_config_path)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Total predictions scored: {len(all_results)}")
    print(f"  Top {len(top_binders)} selected (ipTM range: "
          f"{top_binders[-1]['iptm']:.3f} - {top_binders[0]['iptm']:.3f})")
    if hotspot_residues:
        with_hs = [b for b in top_binders if b.get("hotspot_residues_contacted")]
        print(f"  With hotspot contacts: {len(with_hs)}/{len(top_binders)}")

    if args.config:
        print(f"\nNext step:")
        print(f"  python3 setup_campaign.py {screen_dir.parent / 'campaign_config_opto.json'}")


if __name__ == "__main__":
    main()
