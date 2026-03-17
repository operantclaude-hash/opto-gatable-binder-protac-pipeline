#!/usr/bin/env python3
"""
Generate binding validation YAMLs for the expanded binder pool.

These are 16 binder designs from the original Stage 1 screen that had:
  - ipTM > 0.85 in the original screen
  - >= 10/14 hotspot residues contacted (on-target epitope)

but were NEVER put through binding validation (binder-alone + HAT/catcore).

This expands the validation pool from 16 peptides to 32 peptides,
doubling the candidates available for the final top-10 selection.
"""

from pathlib import Path

OUT_DIR = Path(__file__).parent / "boltz2_yamls"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HAT_SEQ = (
    "KFSAKRLPSTRLGTFLENRVNDFLRRQNHPESGEVTVRVVHASDKTVEVKPGMKARFVDS"
    "GEMAESFPYRTKALFAFEEIDGVDLCFFGMHVQEYGSDCPPPNQRRVYISYLDSVHFFR"
    "PKCLRTAVYHEILIGYLEYVKKLGYTTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQE"
    "WYKKMLDKAVSERIVHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKESQKLYA"
    "TMEKHKEVFFVIRLIAGPAANSLPPIVDPDPLIPCDLMDGRDAFLTLARDRHLEFSSLRR"
    "AQWSTGCMLVELHTQSQD"
)

# Expanded pool: 16 untested binders with ipTM > 0.85 and >= 10/14 hotspots
EXPANDED_BINDERS = {
    "pep6s1":  {"seq": "AAAGQASDKPEQAAEVAAALGAARDAAVAA", "screen_iptm": 0.9049, "hotspots": 12},
    "pep11s1": {"seq": "SSRPNPPRLDARIEESRAGL", "screen_iptm": 0.8857, "hotspots": 11},
    "pep49s5": {"seq": "DSLGSGCESGSTGRERLARCSEILKAQLEA", "screen_iptm": 0.8848, "hotspots": 11},
    "pep9s5":  {"seq": "MTPRKPSDAGLVRAAEKAAA", "screen_iptm": 0.8821, "hotspots": 11},
    "pep29s4": {"seq": "MMKVVFGGGCGIGTPISEEDIRRGNIGGKEGGAESKARSAARAAAAAAAL", "screen_iptm": 0.8815, "hotspots": 11},
    "pep4s4":  {"seq": "GCMPATSEEAARQYAACLAR", "screen_iptm": 0.8765, "hotspots": 11},
    "pep30s4": {"seq": "AVDAKQSKETGAEALAAAADAVSAAAAAAE", "screen_iptm": 0.8711, "hotspots": 10},
    "pep41s5": {"seq": "DDVERGRAAGAAGTGARGEF", "screen_iptm": 0.8689, "hotspots": 11},
    "pep34s1": {"seq": "SEDREKAERHQAALEEVVEATRAAYREDRELGGRVAREELGRLERAEPDD", "screen_iptm": 0.8687, "hotspots": 10},
    "pep46s3": {"seq": "AAGGGGPAKDEASRRNMESARRKARRVAAL", "screen_iptm": 0.8673, "hotspots": 12},
    "pep37s3": {"seq": "ETPADRDAAGGRAAAARLNL", "screen_iptm": 0.8666, "hotspots": 11},
    "pep19s2": {"seq": "ATVADPESGGSAVTAGAPVE", "screen_iptm": 0.8650, "hotspots": 11},
    "pep16s4": {"seq": "AGMPLGGRTEEELRTSRNAAEGDAENRARS", "screen_iptm": 0.8635, "hotspots": 11},
    "pep40s3": {"seq": "SMETFSSAEAQANYEKGLAGREAGDARNQARRDAARAAAA", "screen_iptm": 0.8624, "hotspots": 11},
    "pep43s1": {"seq": "AVASTPAAAKGGGYAEPGAP", "screen_iptm": 0.8591, "hotspots": 12},
    "pep14s4": {"seq": "PETEADKAMRARTAAAIAAA", "screen_iptm": 0.8583, "hotspots": 10},
}


def write_yaml(name, chain_a_seq, chain_b_seq):
    """Write a 2-chain Boltz-2 YAML."""
    yaml_path = OUT_DIR / f"{name}.yaml"
    yaml_path.write_text(f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {chain_a_seq}
      msa: empty
  - protein:
      id: B
      sequence: {chain_b_seq}
      msa: empty
""")
    return yaml_path


def main():
    count = 0
    binder_lengths = {}

    for binder_name, info in sorted(EXPANDED_BINDERS.items()):
        binder_seq = info["seq"]
        binder_lengths[binder_name] = len(binder_seq)

        # Binder alone + HAT (primary validation)
        name = f"{binder_name}_{len(binder_seq)}aa_binder_HAT"
        write_yaml(name, binder_seq, HAT_SEQ)
        count += 1

    print(f"Generated {count} binder-alone + HAT YAMLs in {OUT_DIR}")

    # Write run script
    run_script = Path(__file__).parent / "run_expanded_predictions.sh"
    yamls = sorted(OUT_DIR.glob("*.yaml"))
    lines = ["#!/bin/bash", "set -e", f"cd {Path(__file__).parent}", ""]
    for y in yamls:
        out_dir = Path(__file__).parent / "boltz2_results" / y.stem
        lines.append(
            f"echo 'Running {y.stem}...' && "
            f"conda run -n boltz_only boltz predict {y} "
            f"--out_dir {out_dir} "
            f"--recycling_steps 3 --diffusion_samples 5 --model boltz2"
        )
    lines.append("")
    lines.append("echo 'All expanded pool predictions complete!'")
    run_script.write_text("\n".join(lines) + "\n")
    run_script.chmod(0o755)
    print(f"Run script: {run_script}")

    # Summary
    print(f"\nExpanded pool: {len(EXPANDED_BINDERS)} new binders")
    print("Binder sizes:")
    for name, length in sorted(binder_lengths.items()):
        iptm = EXPANDED_BINDERS[name]["screen_iptm"]
        hs = EXPANDED_BINDERS[name]["hotspots"]
        print(f"  {name:<10} {length:>2}aa  screen_ipTM={iptm:.4f}  hotspots={hs}/14")


if __name__ == "__main__":
    main()
