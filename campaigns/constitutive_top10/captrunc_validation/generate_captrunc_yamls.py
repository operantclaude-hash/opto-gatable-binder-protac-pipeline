#!/usr/bin/env python3
"""
Generate cap-truncated fusion YAMLs for Boltz-2 validation.

Construct: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]  (NO UBQ cap)

This removes the cap-docking artifact while testing whether the scaffold
disrupts binder-target binding.
"""

import json
from pathlib import Path

OUT_DIR = Path(__file__).parent / "boltz2_yamls"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Scaffold components
SPACER = "GGSGGS"
LOV2_LIT = (
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
    "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAER"
    "EG"
)
E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"

# Target sequences
HAT_SEQ = (
    "KFSAKRLPSTRLGTFLENRVNDFLRRQNHPESGEVTVRVVHASDKTVEVKPGMKARFVDS"
    "GEMAESFPYRTKALFAFEEIDGVDLCFFGMHVQEYGSDCPPPNQRRVYISYLDSVHFFR"
    "PKCLRTAVYHEILIGYLEYVKKLGYTTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQE"
    "WYKKMLDKAVSERIVHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKESQKLYA"
    "TMEKHKEVFFVIRLIAGPAANSLPPIVDPDPLIPCDLMDGRDAFLTLARDRHLEFSSLRR"
    "AQWSTGCMLVELHTQSQD"
)

# Read catcore sequence from an existing YAML
CATCORE_SEQ = None
catcore_yamls = list(Path(__file__).parent.parent.rglob("*dark_catcore*.yaml"))
if catcore_yamls:
    with open(catcore_yamls[0]) as f:
        for line in f:
            if "sequence:" in line and CATCORE_SEQ is None:
                pass  # skip chain A
            elif "sequence:" in line:
                CATCORE_SEQ = line.strip().split("sequence:")[1].strip()
                break

# On-target and partial binders (from epitope analysis)
BINDERS = {
    # ON-TARGET (>= 50% hotspots)
    "pep30s5": {"seq": "APPGKTAAEDAARRAAAAAA", "hotspots": 8, "iptm": 0.856},
    "pep25s1": {"seq": "AKGVANSTTDAAANAARIGELAAARLAAGL", "hotspots": 8, "iptm": 0.833},
    "pep21s5": {"seq": "SIGSSSSGASGSLGAVTAED", "hotspots": 9, "iptm": 0.833},
    "pep15s3": {"seq": "LGKGSTPAATRSWEAMAAKL", "hotspots": 8, "iptm": 0.825},
    "pep48s5": {"seq": "PIETGPGGAGKDAQRLDEKL", "hotspots": 8, "iptm": 0.776},
    "pep16s5": {"seq": "PRMPAGGRTEAERRTSENAGRGYAENRARS", "hotspots": 7, "iptm": 0.707},
    "pep36s1": {"seq": "PVLEYDESGSSGEGGLGSVGRGLVEAREGGERGEEEDRRR", "hotspots": 7, "iptm": 0.534},
    # PARTIAL (25-49% hotspots, high ipTM)
    "pep29s5": {"seq": "VEEWGRRGEEMADVHPGGRGAARMARAVLE", "hotspots": 6, "iptm": 0.888},
    "pep21s2": {"seq": "SIGSSSSGGSGSLGSVSGSS", "hotspots": 6, "iptm": 0.840},
    "pep24s3": {"seq": "SLAEAERARQRARREAAHAA", "hotspots": 6, "iptm": 0.689},
    "pep13s3": {"seq": "AAQSGFSPASEAAAAELLAA", "hotspots": 4, "iptm": 0.858},
}


def write_yaml(name, chain_a_seq, chain_b_seq, chain_b_name="target"):
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

    for binder_name, info in BINDERS.items():
        binder_seq = info["seq"]
        fusion_seq = binder_seq + SPACER + LOV2_LIT + E2PPPE2

        # Cap-truncated fusion + HAT
        name = f"{binder_name}_captrunc_HAT"
        write_yaml(name, fusion_seq, HAT_SEQ)
        count += 1

        # Cap-truncated fusion + catcore (if we have the sequence)
        if CATCORE_SEQ:
            name = f"{binder_name}_captrunc_catcore"
            write_yaml(name, fusion_seq, CATCORE_SEQ)
            count += 1

    # Also generate a scaffold-only control (no binder, just LOV2+linker)
    scaffold_only = "AAAAAGGSGGS" + SPACER + LOV2_LIT + E2PPPE2  # 5aa dummy + scaffold
    write_yaml("scaffold_only_HAT", scaffold_only, HAT_SEQ)
    count += 1

    print(f"Generated {count} YAMLs in {OUT_DIR}")

    # Write run script
    run_script = Path(__file__).parent / "run_captrunc_predictions.sh"
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
    lines.append("echo 'All cap-truncated predictions complete!'")
    run_script.write_text("\n".join(lines) + "\n")
    run_script.chmod(0o755)
    print(f"Run script: {run_script}")


if __name__ == "__main__":
    main()
