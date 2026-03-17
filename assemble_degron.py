#!/usr/bin/env python3
"""
Assemble dual-gated VVD-degron constructs from validated opto-gated fusions.

Takes Architecture C fusions ([Binder]-[spacer]-[LOV2]-[linker]-[Cap])
and prepends or appends the VVD + E3 ubiquitin ligase degron system.

Supports two placement modes:
  N-terminal (default): [E3]-[GS]-[VVD]-[GS]-[Binder]-[spacer]-[LOV2]-[linker]-[Cap]
    - E3 near target for efficient ubiquitination
    - Standard in bioPROTAC literature (Clift 2023, Columbia DiVa 2025)

  C-terminal:           [Binder]-[spacer]-[LOV2]-[linker]-[Cap]-[GS]-[VVD]-[GS]-[E3]
    - Binder interface unobstructed
    - E3 distal from target (may reduce ubiquitination efficiency)

Supports two degradation pathways:
  proteasomal: TRIM21_RBCC (236aa, K48 chains) — cytosolic targets (e.g., p300)
  lysosomal:   NEDD4L_HECT (377aa, K63 chains) — membrane targets

Usage:
  python assemble_degron.py campaigns/constitutive_top10/campaign_config.json
  python assemble_degron.py campaigns/constitutive_top10/campaign_config.json --pathway proteasomal
  python assemble_degron.py campaigns/constitutive_top10/campaign_config.json --placement c_term
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from degron_library import (
    VVD_SEQ, VVD_LENGTH,
    E3_LIGASE_LIBRARY,
    PATHWAY_DEFAULTS,
    GS_LINKER_LIBRARY,
    AAV_MAX_BP, AAV_MAX_AA,
)
from linker_library import LOV2_DARK, LOV2_LIT
from setup_campaign import make_yaml_1chain, make_yaml_2chain


# ============================================================================
# Configuration validation
# ============================================================================

def validate_degron_config(cfg):
    """Validate the degron section of a campaign config.

    Returns resolved config dict with defaults applied.
    """
    degron = cfg.get("degron", {})
    errors = []

    pathway = degron.get("pathway", "proteasomal")
    if pathway not in PATHWAY_DEFAULTS:
        errors.append(f"Invalid pathway '{pathway}'. Choose: {list(PATHWAY_DEFAULTS.keys())}")

    e3_name = degron.get("e3_ligase", PATHWAY_DEFAULTS.get(pathway, "TRIM21_RBCC"))
    if e3_name not in E3_LIGASE_LIBRARY:
        errors.append(f"Invalid e3_ligase '{e3_name}'. Choose: {list(E3_LIGASE_LIBRARY.keys())}")

    placement = degron.get("placement", "n_term")
    if placement not in ("n_term", "c_term"):
        errors.append(f"Invalid placement '{placement}'. Choose: n_term, c_term")

    for key in ("e3_vvd_linker", "vvd_opto_linker"):
        linker_name = degron.get(key, "GS15")
        if linker_name not in GS_LINKER_LIBRARY:
            errors.append(f"Invalid {key} '{linker_name}'. Choose: {list(GS_LINKER_LIBRARY.keys())}")

    if errors:
        raise ValueError("Degron config errors:\n  " + "\n  ".join(errors))

    return {
        "pathway": pathway,
        "e3_ligase": e3_name,
        "placement": placement,
        "e3_vvd_linker": degron.get("e3_vvd_linker", "GS15"),
        "vvd_opto_linker": degron.get("vvd_opto_linker", "GS15"),
    }


# ============================================================================
# Component resolution
# ============================================================================

def resolve_degron_components(degron_cfg):
    """Resolve all degron component sequences from config.

    Returns dict with 'e3_seq', 'vvd_seq', 'gs1_seq', 'gs2_seq',
    and their names/lengths.
    """
    e3_name = degron_cfg["e3_ligase"]
    e3_info = E3_LIGASE_LIBRARY[e3_name]
    gs1_name = degron_cfg["e3_vvd_linker"]
    gs2_name = degron_cfg["vvd_opto_linker"]

    return {
        "e3_name": e3_name,
        "e3_seq": e3_info["seq"],
        "e3_length": e3_info["length"],
        "e3_pathway": e3_info["pathway"],
        "e3_ub_chain": e3_info["ubiquitin_chain"],
        "vvd_seq": VVD_SEQ,
        "vvd_length": VVD_LENGTH,
        "gs1_name": gs1_name,
        "gs1_seq": GS_LINKER_LIBRARY[gs1_name]["seq"],
        "gs1_length": GS_LINKER_LIBRARY[gs1_name]["length"],
        "gs2_name": gs2_name,
        "gs2_seq": GS_LINKER_LIBRARY[gs2_name]["seq"],
        "gs2_length": GS_LINKER_LIBRARY[gs2_name]["length"],
        "placement": degron_cfg["placement"],
    }


# ============================================================================
# Fusion assembly
# ============================================================================

def build_degron_fusion(opto_fusion_seq, opto_domains, degron_components):
    """Build full degron + opto-gated fusion construct.

    N-terminal: [E3]-[GS1]-[VVD]-[GS2]-[opto_fusion]
    C-terminal: [opto_fusion]-[GS2]-[VVD]-[GS1]-[E3]

    Args:
        opto_fusion_seq: The opto-gated fusion sequence string
        opto_domains: Domain map from build_archC_fusion {name: [start, end)}
        degron_components: Dict from resolve_degron_components

    Returns:
        (full_seq, full_domains) where full_domains includes both
        degron and shifted opto domains.
    """
    dc = degron_components
    placement = dc["placement"]

    if placement == "n_term":
        # [E3]-[GS1]-[VVD]-[GS2]-[opto_fusion]
        prefix = dc["e3_seq"] + dc["gs1_seq"] + dc["vvd_seq"] + dc["gs2_seq"]
        full_seq = prefix + opto_fusion_seq
        prefix_len = len(prefix)

        pos = 0
        domains = {}
        domains["e3_ligase"] = [pos, pos + dc["e3_length"]]
        pos += dc["e3_length"]
        domains["gs1"] = [pos, pos + dc["gs1_length"]]
        pos += dc["gs1_length"]
        domains["vvd"] = [pos, pos + dc["vvd_length"]]
        pos += dc["vvd_length"]
        domains["gs2"] = [pos, pos + dc["gs2_length"]]
        pos += dc["gs2_length"]

    else:
        # [opto_fusion]-[GS2]-[VVD]-[GS1]-[E3]
        suffix = dc["gs2_seq"] + dc["vvd_seq"] + dc["gs1_seq"] + dc["e3_seq"]
        full_seq = opto_fusion_seq + suffix
        prefix_len = 0

        domains = {}

    # Shift opto domains by prefix length
    for name, (start, end) in opto_domains.items():
        domains[name] = [start + prefix_len, end + prefix_len]

    if placement == "c_term":
        pos = len(opto_fusion_seq)
        domains["gs2"] = [pos, pos + dc["gs2_length"]]
        pos += dc["gs2_length"]
        domains["vvd"] = [pos, pos + dc["vvd_length"]]
        pos += dc["vvd_length"]
        domains["gs1"] = [pos, pos + dc["gs1_length"]]
        pos += dc["gs1_length"]
        domains["e3_ligase"] = [pos, pos + dc["e3_length"]]
        pos += dc["e3_length"]

    assert len(full_seq) == max(e for _, (_, e) in domains.items()), \
        "Domain map end does not match sequence length"

    return full_seq, domains


# ============================================================================
# AAV compatibility
# ============================================================================

def check_aav_compatibility(total_length_aa):
    """Check if construct fits within AAV packaging limit.

    Returns dict with 'compatible', 'total_bp', 'remaining_bp', 'warning'.
    """
    total_bp = total_length_aa * 3
    remaining_bp = AAV_MAX_BP - total_bp
    warning = None

    if total_bp > AAV_MAX_BP:
        compatible = False
        warning = (f"Exceeds AAV limit: {total_bp}bp > {AAV_MAX_BP}bp "
                   f"(over by {total_bp - AAV_MAX_BP}bp)")
    elif remaining_bp < 500:
        compatible = True
        warning = f"Near AAV limit: {total_bp}bp, only {remaining_bp}bp remaining"
    else:
        compatible = True

    return {
        "compatible": compatible,
        "total_aa": total_length_aa,
        "total_bp": total_bp,
        "remaining_bp": remaining_bp,
        "aav_max_bp": AAV_MAX_BP,
        "warning": warning,
    }


# ============================================================================
# Campaign assembly
# ============================================================================

def assemble_campaign_degrons(config_path, pathway_override=None,
                              placement_override=None):
    """Load campaign config, build degron fusions for each validated construct.

    Produces:
      - degron_assembly/constructs.fasta
      - degron_assembly/manifest.json
      - degron_assembly/boltz2_yamls/*.yaml (fold predictions)

    Args:
        config_path: Path to campaign_config.json
        pathway_override: Override pathway from config
        placement_override: Override placement from config
    """
    config_path = Path(config_path)
    with open(config_path) as f:
        cfg = json.load(f)

    campaign_dir = config_path.parent
    campaign_name = cfg.get("campaign_name", campaign_dir.name)

    # Apply overrides
    if "degron" not in cfg:
        cfg["degron"] = {}
    if pathway_override:
        cfg["degron"]["pathway"] = pathway_override
    if placement_override:
        cfg["degron"]["placement"] = placement_override

    # Validate and resolve
    degron_cfg = validate_degron_config(cfg)
    degron_components = resolve_degron_components(degron_cfg)

    # Output dirs
    out_dir = campaign_dir / "degron_assembly"
    out_dir.mkdir(exist_ok=True)
    yaml_dir = out_dir / "boltz2_yamls"
    yaml_dir.mkdir(exist_ok=True)

    # Load constructs from campaign manifest
    manifest_path = campaign_dir / "manifests" / "campaign_manifest.json"
    if not manifest_path.exists():
        # Try alternate location
        manifest_path = campaign_dir / "campaign_manifest.json"
    if not manifest_path.exists():
        print(f"No campaign manifest found at {manifest_path}")
        print("Run setup_campaign.py first to generate constructs.")
        return

    with open(manifest_path) as f:
        constructs = json.load(f)

    target_seq = cfg["target"]["sequence"]

    print("=" * 70)
    print(f"DEGRON ASSEMBLY: {campaign_name}")
    print(f"  Pathway:   {degron_cfg['pathway']}")
    print(f"  E3 ligase: {degron_cfg['e3_ligase']} ({degron_components['e3_length']}aa)")
    print(f"  VVD:       {VVD_LENGTH}aa")
    print(f"  Placement: {degron_cfg['placement']}")
    print(f"  Linkers:   {degron_cfg['e3_vvd_linker']}/{degron_cfg['vvd_opto_linker']}")
    print(f"  Constructs: {len(constructs)}")
    print("=" * 70)

    fasta_lines = []
    degron_manifest = []
    yaml_names = []

    for construct in constructs:
        name = construct["name"]
        opto_seq_dark = construct["fusion_seq_dark"]
        opto_seq_lit = construct["fusion_seq_lit"]
        opto_domains = construct["domain_map"]

        # Build dark and lit degron fusions
        for state, opto_seq in [("dark", opto_seq_dark), ("lit", opto_seq_lit)]:
            full_seq, full_domains = build_degron_fusion(
                opto_seq, opto_domains, degron_components
            )
            degron_name = f"{name}_degron_{state}"

            aav = check_aav_compatibility(len(full_seq))

            fasta_lines.append(f">{degron_name} len={len(full_seq)} {aav['total_bp']}bp")
            fasta_lines.append(full_seq)

            # Fold YAML (single chain)
            yaml_content = make_yaml_1chain(full_seq, degron_name)
            yaml_path = yaml_dir / f"{degron_name}.yaml"
            yaml_path.write_text(yaml_content)
            yaml_names.append(degron_name)

            # 2-chain YAML with target
            tch_name = f"{degron_name}_target"
            tch_yaml = make_yaml_2chain(full_seq, target_seq)
            (yaml_dir / f"{tch_name}.yaml").write_text(tch_yaml)
            yaml_names.append(tch_name)

            entry = {
                "name": degron_name,
                "source_construct": name,
                "state": state,
                "total_length_aa": len(full_seq),
                "total_bp": aav["total_bp"],
                "aav_compatible": aav["compatible"],
                "aav_warning": aav["warning"],
                "domain_map": full_domains,
                "degron": {
                    "pathway": degron_cfg["pathway"],
                    "e3_ligase": degron_cfg["e3_ligase"],
                    "placement": degron_cfg["placement"],
                    "e3_vvd_linker": degron_cfg["e3_vvd_linker"],
                    "vvd_opto_linker": degron_cfg["vvd_opto_linker"],
                },
            }
            degron_manifest.append(entry)

            status = "OK" if aav["compatible"] else "OVER"
            warn = f" [{aav['warning']}]" if aav["warning"] else ""
            print(f"  {degron_name}: {len(full_seq)}aa ({aav['total_bp']}bp) [{status}]{warn}")

    # Write outputs
    fasta_path = out_dir / "constructs.fasta"
    fasta_path.write_text("\n".join(fasta_lines) + "\n")

    manifest_out = out_dir / "degron_manifest.json"
    manifest_out.write_text(json.dumps(degron_manifest, indent=2) + "\n")

    print(f"\nOutputs:")
    print(f"  FASTA:    {fasta_path}")
    print(f"  Manifest: {manifest_out}")
    print(f"  YAMLs:    {yaml_dir}/ ({len(yaml_names)} files)")

    # AAV summary
    lengths = [e["total_length_aa"] for e in degron_manifest]
    if lengths:
        print(f"\nSize range: {min(lengths)}-{max(lengths)}aa "
              f"({min(lengths)*3}-{max(lengths)*3}bp)")
        n_compat = sum(1 for e in degron_manifest if e["aav_compatible"])
        print(f"AAV compatible: {n_compat}/{len(degron_manifest)}")


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Assemble VVD-degron fusions from opto-gated constructs",
    )
    parser.add_argument("config", type=str, help="Path to campaign_config.json")
    parser.add_argument("--pathway", choices=["proteasomal", "lysosomal"],
                        help="Degradation pathway (overrides config)")
    parser.add_argument("--placement", choices=["n_term", "c_term"],
                        help="Degron placement (overrides config)")
    args = parser.parse_args()

    assemble_campaign_degrons(
        args.config,
        pathway_override=args.pathway,
        placement_override=args.placement,
    )


if __name__ == "__main__":
    main()
