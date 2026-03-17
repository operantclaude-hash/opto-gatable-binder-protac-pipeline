#!/usr/bin/env python3
"""
Generate 2-chain Boltz-2 YAMLs for binding validation of top 22 binders.

Tests a 2x2 matrix:
  - Binder alone vs Fusion (constitutively open)
  - HAT domain (317aa) vs Catalytic core (563aa)

This isolates:
  - Binder→Fusion: does LOV2 integration kill binding?
  - HAT→Catcore: do flanking domains cause steric clashes?

Generates 88 YAMLs: 22 binders × 4 conditions
"""

import json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
YAML_DIR = SCRIPT_DIR / "binding_validation" / "boltz2_yamls"
RESULTS_DIR = SCRIPT_DIR / "binding_validation" / "boltz2_results"
YAML_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ===== Target sequences =====

HAT_SEQ = (
    "KFSAKRLPSTRLGTFLENRVNDFLRRQNHPESGEVTVRVVHASDKTVEVKPGMKARFVDS"
    "GEMAESFPYRTKALFAFEEIDGVDLCFFGMHVQEYGSDCPPPNQRRVYISYLDSVHFFR"
    "PKCLRTAVYHEILIGYLEYVKKLGYTTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQE"
    "WYKKMLDKAVSERIVHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKESQKLYA"
    "TMEKHKEVFFVIRLIAGPAANSLPPIVDPDPLIPCDLMDGRDAFLTLARDRHLEFSSLRR"
    "AQWSTGCMLVELHTQSQD"
)  # 317aa

CATCORE_SEQ = (
    "KKIFKPEELRQALMPTLEALYRQDPESLPFRQPVDPQLLGIPDYFDIVKSPMDLSTIKRK"
    "LDTGQYQEPWQYVDDIWLMFNNAWLYNRKTSRVYKYCSKLSEVFEQEIDPVMQSLGYCC"
    "GRKLEFSPQTLCCYGKQLCTIPRDATYYSYQNRYHFCEKCFNEIQGESVSLGDDPSQPQT"
    "TINKEQFSKRKNDTLDPELFVECTECGRKMHQICVLHHEIIWPAGFVCDGCLKKSARTRKE"
    "NKFSAKRLPSTRLGTFLENRVNDFLRRQNHPESGEVTVRVVHASDKTVEVKPGMKARFVDS"
    "GEMAESFPYRTKALFAFEEIDGVDLCFFGMHVQEYGSDCPPPNQRRVYISYLDSVHFFRP"
    "KCLRTAVYHEILIGYLEYVKKLGYTTGHIWACPPSEGDDYIFHCHPPDQKIPKPKRLQEWF"
    "KKMLDKAVSERIVHDYKDIFKQATEDRLTSAKELPYFEGDFWPNVLEESIKESGGSGSQKL"
    "YATMEKHKEVFFVIRLIAGPAANSLPPIVDPDPLIPCDLMDGRDAFLTLARDKHLEFSSLR"
    "RAQWSTMCMLVELHTQSQD"
)  # 563aa

# ===== Top 22 binders (sorted by ipTM, with min_iPAE from Stage 1/2) =====

BINDERS = [
    # --- Stage 1 peptides (ipTM > 0.90) ---
    {"name": "pep13s3",  "seq": "AAQSGFSPASEAAAAELLAA",   "iptm": 0.942, "min_ipae": 0.822},
    {"name": "pep10s5",  "seq": "FGSVTGSGAADAVDDALAAA",   "iptm": 0.940, "min_ipae": 0.956},
    {"name": "pep32s1",  "seq": "MGPQVTRENAARAAARVAAR",    "iptm": 0.936, "min_ipae": 1.327},
    {"name": "pep25s1",  "seq": "AKGVANSTTDAAANAARIGELAAARLAAGL", "iptm": 0.927, "min_ipae": 1.360},
    {"name": "pep30s5",  "seq": "APPGKTAAEDAARRAAAAAA",    "iptm": 0.916, "min_ipae": 1.507},
    {"name": "pep15s3",  "seq": "LGKGSTPAATRSWEAMAAKL",    "iptm": 0.912, "min_ipae": 1.497},
    {"name": "pep24s3",  "seq": "SLAEAERARQRARREAAHAA",    "iptm": 0.909, "min_ipae": 1.553},
    {"name": "pep16s5",  "seq": "PRMPAGGRTEAERRTSENAGRGYAENRARS", "iptm": 0.908, "min_ipae": 1.410},
    {"name": "pep17s4",  "seq": "MSELTEVEGPGEASPAAKAH",    "iptm": 0.907, "min_ipae": 1.219},
    {"name": "pep21s5",  "seq": "SIGSSSSGASGSLGAVTAED",    "iptm": 0.907, "min_ipae": 1.336},
    # --- Stage 1 peptides (ipTM 0.89-0.91) ---
    {"name": "pep32s2",  "seq": "MGPAVTRENAARAAARVAAR",    "iptm": 0.905, "min_ipae": 1.554},
    {"name": "pep29s5",  "seq": "VEEWGRRGEEMADVHPGGRGAARMARAVLE", "iptm": 0.904, "min_ipae": 1.361},
    {"name": "pep36s1",  "seq": "PVLEYDESGSSGEGGLGSVGRGLVEAREGGERGEEEDRRR", "iptm": 0.904, "min_ipae": 1.435},
    {"name": "pep48s5",  "seq": "PIETGPGGAGKDAQRLDEKL",    "iptm": 0.904, "min_ipae": 1.215},
    {"name": "pep21s2",  "seq": "SIGSSSSGGSGSLGSVSGSS",    "iptm": 0.899, "min_ipae": 1.577},
    {"name": "pep4s2",   "seq": "SPEEEILEESEAGLKRAEEWGKVAEACRKAAAAGDAAAVAACQAQLDRFG", "iptm": 0.898, "min_ipae": 0.903},
    # --- Stage 2 nanobody-scaffold binders ---
    {"name": "nb153s2",  "seq": "SLRERAAANSKEAAAQTVVGAARLRAAVDAVDPAVADAVAAQYTANRNTLNAANAALAAIAAESLDRIKATYEAEGSVAALLAAEEAEGRRQQAAVLAAREAASAASRAIGAIVAWLRAQ", "iptm": 0.880, "min_ipae": 1.213},
    {"name": "nb39s3",   "seq": "MVLLGAIEGGKEGVRNRREREELALAGEDERARELDRRAGEGGGAFGGLATGEDVHEAYKIAEKGYEKREKEIKARMAELEKTGEPAAAAQAEAWRKAMEEEKKAEKEANKDAGAAGVGG", "iptm": 0.861, "min_ipae": 1.580},
    {"name": "nb62s3",   "seq": "DLAAAALEVCEETGAQSRAVAAANQAAYAAQIAALAARIAALRAEAAARAAELAARAEATGSEEVRRERERAEAAAAAELAALEAQRAAAERNRAATAERLAANEEGLAAERARLEAEAA", "iptm": 0.847, "min_ipae": 1.088},
    {"name": "nb134s2",  "seq": "SAREAWLEEARAALDEITAENLAGAAANRRAYAEQLAAARARKAALEAEDAALLARVRAECAATGDEERLREVEREIAARLAAIDKEIANLERALAHAEQLQVDIVRARERDLAELEAEA", "iptm": 0.836, "min_ipae": 2.206},
    {"name": "nb178s5",  "seq": "SAAKQAGIENAKALGEEAVAATEANAAALAADAAARLAELEKEIEEYKKKGEEELKELKKKSKSEEEIEKKKKEYEEKIKEIEEQMKTLERRVAAARENVEKAKEGKVVGDKEDAAFFNA", "iptm": 0.827, "min_ipae": 1.817},
    {"name": "nb87s2",   "seq": "SQEIIEAYEEIGRESEEMSRRNCEEAEAYAAEVRARRAEIEARVAEAEAVLAADAAGAAAAEARERLAAAEAELARAEVELRRAEGQLRWCKAAEERLKEMKVTGEKEAAELAAFDEKAK", "iptm": 0.818, "min_ipae": 1.670},
]

# ===== LOV2 + linker + cap for constitutive fusion =====

LOV2_LIT = (
    "LERIEKNFVITDPRLPDNPIIFASDSFLQLTEYSREEILGRNCRFLQGPETDRATVRKIR"
    "DAIDNQTEVTVQLINYTKSGKKFWNLFHLQPMRDQKGDVQYFIGVQLDGTEHVRDAAER"
    "EG"
)  # 121aa

E2PPPE2 = "EAAAKEAAAKPPPEAAAKEAAAK"  # 23aa (E2pppE2 - top gating linker)

UBQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"  # 76aa

SPACER = "GGSGGS"  # 6aa


def make_2chain_yaml(chain_a_seq, chain_b_seq):
    """Generate Boltz-2 YAML for 2-chain prediction (no constraints)."""
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {chain_a_seq}",
        "      msa: empty",
        "  - protein:",
        "      id: B",
        f"      sequence: {chain_b_seq}",
        "      msa: empty",
    ]
    return "\n".join(lines) + "\n"


def build_fusion(binder_seq):
    """Build constitutively open fusion: [Binder]-[GGSGGS]-[LOV2_lit]-[E2pppE2]-[UBQ]"""
    return binder_seq + SPACER + LOV2_LIT + E2PPPE2 + UBQ


def main():
    conditions = [
        ("binder_HAT",     lambda b: b["seq"],          HAT_SEQ),
        ("fusion_HAT",     lambda b: build_fusion(b["seq"]), HAT_SEQ),
        ("binder_catcore", lambda b: b["seq"],          CATCORE_SEQ),
        ("fusion_catcore", lambda b: build_fusion(b["seq"]), CATCORE_SEQ),
    ]

    manifest = []
    yaml_names = []

    n_binders = len(BINDERS)
    print(f"Binding Validation: {n_binders} binders x 4 conditions = {n_binders * 4} predictions")
    print("=" * 70)
    print(f"{'Condition':<20} {'Chain A':<25} {'Chain B':<15}")
    print("-" * 70)
    print(f"{'binder_HAT':<20} {'binder alone':<25} {'HAT (317aa)':<15}")
    print(f"{'fusion_HAT':<20} {'constitutive fusion':<25} {'HAT (317aa)':<15}")
    print(f"{'binder_catcore':<20} {'binder alone':<25} {'catcore (563aa)':<15}")
    print(f"{'fusion_catcore':<20} {'constitutive fusion':<25} {'catcore (563aa)':<15}")
    print()

    for binder in BINDERS:
        bname = binder["name"]
        for cond_name, chain_a_fn, chain_b_seq in conditions:
            chain_a_seq = chain_a_fn(binder)
            construct_name = f"{bname}_{cond_name}"

            yaml_content = make_2chain_yaml(chain_a_seq, chain_b_seq)
            yaml_path = YAML_DIR / f"{construct_name}.yaml"
            yaml_path.write_text(yaml_content)
            yaml_names.append(construct_name)

            manifest.append({
                "construct_name": construct_name,
                "binder_name": bname,
                "binder_iptm_stage1": binder["iptm"],
                "binder_min_ipae_stage1": binder.get("min_ipae"),
                "condition": cond_name,
                "chain_a": "fusion" if "fusion" in cond_name else "binder",
                "chain_b": "catcore" if "catcore" in cond_name else "HAT",
                "chain_a_length": len(chain_a_seq),
                "chain_b_length": len(chain_b_seq),
                "binder_length": len(binder["seq"]),
            })

        print(f"  {bname}: 4 YAMLs ({len(binder['seq'])}aa binder, "
              f"{len(build_fusion(binder['seq']))}aa fusion)")

    # Write manifest
    manifest_path = SCRIPT_DIR / "binding_validation" / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    # Generate run script
    run_lines = [
        "#!/bin/bash",
        "# Binding validation: 10 binders x 4 conditions = 40 predictions",
        "# Tests: binder vs fusion, HAT vs catcore",
        "",
        "set +e",
        "source ~/miniconda3/etc/profile.d/conda.sh",
        "conda activate boltz_only",
        "",
        f'YAML_DIR="{YAML_DIR}"',
        f'RESULTS_DIR="{RESULTS_DIR}"',
        "",
        f"TOTAL={len(yaml_names)}",
        "COUNT=0",
        "",
    ]

    for name in yaml_names:
        run_lines.extend([
            'COUNT=$((COUNT + 1))',
            f'echo "[${{COUNT}}/${{TOTAL}}] {name}..."',
            '',
            f'if find "$RESULTS_DIR/{name}" -name "*.cif" 2>/dev/null | grep -q cif; then',
            '    echo "  (already complete, skipping)"',
            'else',
            '    python3 -c "import torch; torch.cuda.empty_cache()" 2>/dev/null',
            f'    boltz predict "$YAML_DIR/{name}.yaml" \\',
            f'        --out_dir "$RESULTS_DIR/{name}" \\',
            '        --model boltz2 \\',
            '        --diffusion_samples 5 \\',
            '        --recycling_steps 3 \\',
            '        --override 2>&1 | tail -3',
            f'    CIF_COUNT=$(find "$RESULTS_DIR/{name}" -name "*.cif" 2>/dev/null | wc -l)',
            '    echo "  $CIF_COUNT models generated"',
            'fi',
            '',
        ])

    run_lines.extend([
        'echo ""',
        f'echo "All {len(yaml_names)} predictions complete."',
        f'echo "Results: $RESULTS_DIR/"',
        'echo "Run: python3 analyze_binding_validation.py"',
    ])

    run_path = SCRIPT_DIR / "binding_validation" / "run_predictions.sh"
    run_path.write_text("\n".join(run_lines) + "\n")
    run_path.chmod(0o755)

    print(f"\nGenerated:")
    print(f"  {len(yaml_names)} YAMLs in {YAML_DIR}/")
    print(f"  Manifest: {manifest_path}")
    print(f"  Runner:   {run_path}")
    print(f"\nNext: bash {run_path}")


if __name__ == "__main__":
    main()
