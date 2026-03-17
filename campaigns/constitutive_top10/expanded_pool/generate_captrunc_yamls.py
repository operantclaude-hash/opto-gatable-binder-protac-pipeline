#!/usr/bin/env python3
"""
Generate cap-truncated fusion YAMLs for expanded pool top candidates.

Construct: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2] (NO UBQ cap)
Tests against both HAT and CatCore.
"""

from pathlib import Path

OUT_DIR = Path(__file__).parent / "captrunc_yamls"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SPACER = "GGSGGS"
LOV2_LIT = (
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
    "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAER"
    "EG"
)
E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"

HAT_SEQ = (
    "KFSAKRLPSTRLGTFLENRVNDFLRRQNHPESGEVTVRVVHASDKTVEVKPGMKARFVDS"
    "GEMAESFPYRTKALFAFEEIDGVDLCFFGMHVQEYGSDCPPPNQRRVYISYLDSVHFFR"
    "PKCLRTAVYHEILIGYLEYVKKLGYTTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQE"
    "WYKKMLDKAVSERIVHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKESQKLYA"
    "TMEKHKEVFFVIRLIAGPAANSLPPIVDPDPLIPCDLMDGRDAFLTLARDRHLEFSSLRR"
    "AQWSTGCMLVELHTQSQD"
)

# Get catcore sequence from existing yamls
CATCORE_SEQ = None
for yaml_path in sorted(Path(__file__).parent.rglob("*catcore*.yaml")):
    with open(yaml_path) as f:
        for line in f:
            if "sequence:" in line:
                seq = line.strip().split("sequence:")[1].strip()
                if len(seq) > 400:
                    CATCORE_SEQ = seq
                    break
    if CATCORE_SEQ:
        break

# Top candidates from expanded pool that need cap-truncated testing
BINDERS = {
    "pep16s4": "AGMPLGGRTEEELRTSRNAAEGDAENRARS",
    "pep29s4": "MMKVVFGGGCGIGTPISEEDIRRGNIGGKEGGAESKARSAARAAAAAAAL",
    "pep11s1": "SSRPNPPRLDARIEESRAGL",
    "pep37s3": "ETPADRDAAGGRAAAARLNL",
    "pep9s5":  "MTPRKPSDAGLVRAAEKAAA",
    "pep41s5": "DDVERGRAAGAAGTGARGEF",
    "pep19s2": "ATVADPESGGSAVTAGAPVE",
}


def write_yaml(name, chain_a_seq, chain_b_seq):
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

    for binder_name, binder_seq in sorted(BINDERS.items()):
        fusion_seq = binder_seq + SPACER + LOV2_LIT + E2PPPE2

        # Cap-truncated + HAT
        write_yaml(f"{binder_name}_captrunc_HAT", fusion_seq, HAT_SEQ)
        count += 1

        # Cap-truncated + CatCore
        if CATCORE_SEQ:
            write_yaml(f"{binder_name}_captrunc_catcore", fusion_seq, CATCORE_SEQ)
            count += 1

    print(f"Generated {count} YAMLs in {OUT_DIR}")

    # Write run script
    run_script = Path(__file__).parent / "run_captrunc_predictions.sh"
    yamls = sorted(OUT_DIR.glob("*.yaml"))
    lines = ["#!/bin/bash", "set -e", f"cd {Path(__file__).parent}", ""]
    for y in yamls:
        out_dir = Path(__file__).parent / "captrunc_results" / y.stem
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
    print(f"\nPredictions: {len([y for y in yamls if 'HAT' in y.name])} HAT + "
          f"{len([y for y in yamls if 'catcore' in y.name])} CatCore = {len(yamls)} total")


if __name__ == "__main__":
    main()
