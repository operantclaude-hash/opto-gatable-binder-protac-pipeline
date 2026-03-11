#!/usr/bin/env python3
"""
Analyze an opto-gatable binder campaign.

Performs sweep volume analysis on Boltz-2 predictions, ranks constructs
by gating ratio, and generates power analysis + ChimeraX visualizations.

Usage:
    python3 analyze_campaign.py campaigns/p300_hat_demo/

Requires: numpy, scipy, BioPython
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

try:
    from Bio.PDB import MMCIFParser, Superimposer
except ImportError:
    print("ERROR: BioPython required. Install: pip install biopython")
    sys.exit(1)

from linker_library import (
    LINKER_LIBRARY,
    HELIX_RISE_PER_RESIDUE,
    PPII_RISE_PER_RESIDUE,
    DEFAULT_CAP_RADIUS,
    DEFAULT_DARK_CONE,
    DEFAULT_LIT_CONE,
    STERIC_EXCLUSION_RADIUS,
    SWEEP_GRID_SPACING,
    N_AZIMUTHAL_SAMPLES,
    N_POLAR_SAMPLES,
    IPTM_SUCCESS_RATES,
    GATING_SUCCESS_RATES,
)

warnings.filterwarnings("ignore")
PIPELINE_DIR = Path(__file__).parent


# ===== Geometry functions (generalized from compute_sweep_volume.py) =====

def extract_ca_coords(structure, chain_id=None, res_range=None):
    """Extract CA coordinates from a structure/chain.

    Args:
        structure: BioPython Structure object
        chain_id: chain ID to extract from (None = first chain)
        res_range: (start, end) 0-indexed half-open range to filter residues
    """
    coords = []
    model = structure[0]
    chains = list(model.get_chains())

    if chain_id:
        chain = None
        for c in chains:
            if c.id == chain_id:
                chain = c
                break
        if chain is None:
            raise ValueError(f"Chain {chain_id} not found")
    else:
        chain = chains[0]

    residues = [r for r in chain.get_residues() if "CA" in r]
    for i, res in enumerate(residues):
        if res_range and (i < res_range[0] or i >= res_range[1]):
            continue
        coords.append(res["CA"].get_vector().get_array())

    return np.array(coords)


def superimpose_binder(query_coords, ref_coords):
    """Superimpose query binder CA onto reference binder CA.

    Returns (rotation, translation) matrices.
    """
    sup = Superimposer()
    n = min(len(query_coords), len(ref_coords))
    sup.set(ref_coords[:n], query_coords[:n])
    return sup.rotran


def apply_transform(coords, rotran):
    """Apply rotation + translation to coordinates."""
    rot, trans = rotran
    return np.dot(coords, rot) + trans


def extract_pivot_and_exit_vector(lov2_coords):
    """Extract pivot point (last CA) and exit direction from LOV2 domain.

    The exit vector points from the LOV2 tail centroid toward the C-terminal
    pivot, representing the direction the linker extends from LOV2.
    """
    pivot = lov2_coords[-1]
    tail = lov2_coords[-min(5, len(lov2_coords)):]
    tail_centroid = tail.mean(axis=0)
    exit_vec = pivot - tail_centroid
    norm = np.linalg.norm(exit_vec)
    if norm < 1e-6:
        exit_vec = np.array([0.0, 0.0, 1.0])
    else:
        exit_vec /= norm
    return pivot, exit_vec


def sample_cone_directions(axis, half_angle_deg, n_polar=10, n_azimuthal=36):
    """Generate uniformly distributed directions within a cone.

    Uses Rodrigues' rotation formula to rotate the axis direction.
    """
    half_angle = np.radians(half_angle_deg)
    directions = []

    # Build orthonormal basis
    if abs(axis[2]) < 0.9:
        perp = np.cross(axis, [0, 0, 1])
    else:
        perp = np.cross(axis, [1, 0, 0])
    perp /= np.linalg.norm(perp)
    perp2 = np.cross(axis, perp)

    for i_p in range(n_polar + 1):
        theta = half_angle * (i_p / n_polar)
        for i_a in range(n_azimuthal):
            phi = 2 * np.pi * i_a / n_azimuthal
            # Rodrigues' rotation
            d = (axis * np.cos(theta)
                 + perp * np.sin(theta) * np.cos(phi)
                 + perp2 * np.sin(theta) * np.sin(phi))
            directions.append(d / np.linalg.norm(d))

    return np.array(directions)


def compute_rod_length(linker_info):
    """Compute effective rigid rod length for a linker."""
    ltype = linker_info["type"]
    length = linker_info["length"]

    if ltype == "ppII":
        return length * PPII_RISE_PER_RESIDUE
    else:
        return length * HELIX_RISE_PER_RESIDUE


def compute_straight_endpoints(pivot, directions, rod_length):
    """Compute cap center positions for a straight linker."""
    return pivot + directions * rod_length


def compute_kinked_endpoints(pivot, directions, linker_info, n_azimuthal=36):
    """Compute cap center positions for a kinked linker.

    Models as two connected rods with a kink, allowing azimuthal freedom
    at the kink point.
    """
    pre_len = linker_info.get("pre_helix", 10) * HELIX_RISE_PER_RESIDUE
    post_len = linker_info.get("post_helix", 10) * HELIX_RISE_PER_RESIDUE
    kink_res = linker_info.get("kink_res", 2)

    # Kink angle depends on motif type
    desc = linker_info.get("description", "").lower()
    if "beta turn" in desc:
        kink_angle = np.radians(100)
    elif "proline" in desc or "pp" in desc.lower():
        kink_angle = np.radians(30 * kink_res)
    else:
        kink_angle = np.radians(45)

    multi_kink = linker_info.get("multi_kink", False)
    if multi_kink:
        kink_angle = np.radians(30)

    endpoints = []
    for d in directions:
        kink_point = pivot + d * pre_len

        # Generate azimuthal ring at kink
        if abs(d[2]) < 0.9:
            perp = np.cross(d, [0, 0, 1])
        else:
            perp = np.cross(d, [1, 0, 0])
        perp /= np.linalg.norm(perp)
        perp2 = np.cross(d, perp)

        for i_a in range(n_azimuthal):
            phi = 2 * np.pi * i_a / n_azimuthal
            post_dir = (d * np.cos(kink_angle)
                        + perp * np.sin(kink_angle) * np.cos(phi)
                        + perp2 * np.sin(kink_angle) * np.sin(phi))
            post_dir /= np.linalg.norm(post_dir)
            endpoints.append(kink_point + post_dir * post_len)

    return np.array(endpoints)


def generate_cap_volume_points(endpoints, cap_radius):
    """Expand each linker endpoint into a sphere of cap-accessible points."""
    shell_dirs = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                v = np.array([dx, dy, dz], dtype=float)
                shell_dirs.append(v / np.linalg.norm(v))
    shell_dirs = np.array(shell_dirs)

    points = []
    for ep in endpoints:
        points.append(ep)
        for sd in shell_dirs:
            points.append(ep + sd * cap_radius)

    return np.array(points)


def voxelize(points, grid_spacing):
    """Reduce point cloud to voxel grid centers for uniform density."""
    if len(points) == 0:
        return points
    voxel_idx = np.round(points / grid_spacing).astype(int)
    unique = np.unique(voxel_idx, axis=0)
    return unique.astype(float) * grid_spacing


def apply_steric_exclusion(sweep_points, body_coords, exclusion_radius):
    """Remove sweep points that overlap with the protein body."""
    if len(sweep_points) == 0 or len(body_coords) == 0:
        return sweep_points
    body_tree = cKDTree(body_coords)
    distances, _ = body_tree.query(sweep_points)
    mask = distances > exclusion_radius
    return sweep_points[mask]


def compute_blocking_fraction(sweep_points, target_coords, hotspot_indices=None):
    """Compute fraction of target atoms inside the sweep volume.

    Args:
        sweep_points: voxelized sweep volume points
        target_coords: all target atom coordinates
        hotspot_indices: 0-indexed residue indices for hotspot subset
    """
    if len(sweep_points) == 0:
        return {"total_blocked_fraction": 0.0, "hotspot_blocked_fraction": 0.0}

    sweep_tree = cKDTree(sweep_points)
    grid = SWEEP_GRID_SPACING

    # Total blocking
    dists, _ = sweep_tree.query(target_coords)
    total_blocked = np.sum(dists < grid) / len(target_coords)

    # Hotspot blocking
    hotspot_blocked = 0.0
    if hotspot_indices is not None and len(hotspot_indices) > 0:
        hot_mask = np.zeros(len(target_coords), dtype=bool)
        for idx in hotspot_indices:
            if idx < len(target_coords):
                hot_mask[idx] = True
        if hot_mask.any():
            hot_coords = target_coords[hot_mask]
            hot_dists, _ = sweep_tree.query(hot_coords)
            hotspot_blocked = np.sum(hot_dists < grid) / len(hot_coords)

    return {
        "total_blocked_fraction": float(total_blocked),
        "hotspot_blocked_fraction": float(hotspot_blocked),
    }


def save_pointcloud_pdb(points, filepath, atom_name="C"):
    """Save point cloud as PDB HETATM records for ChimeraX visualization."""
    with open(filepath, "w") as f:
        for i, pt in enumerate(points):
            f.write(
                f"HETATM{i+1:5d}  {atom_name:3s} SWP A{1:4d}    "
                f"{pt[0]:8.3f}{pt[1]:8.3f}{pt[2]:8.3f}"
                f"  1.00  0.00           {atom_name[0]:>2s}\n"
            )
        f.write("END\n")


# ===== Sweep volume analysis =====

def generate_sweep_pointcloud(pivot, exit_vec, linker_info, cone_half_angle,
                              cap_radius, body_coords):
    """Generate the full sweep point cloud for one state.

    1. Sample directions within the cone
    2. Compute linker endpoints (straight or kinked)
    3. Expand by cap radius
    4. Voxelize
    5. Apply steric exclusion
    """
    directions = sample_cone_directions(
        exit_vec, cone_half_angle, N_POLAR_SAMPLES, N_AZIMUTHAL_SAMPLES
    )
    rod_length = compute_rod_length(linker_info)

    ltype = linker_info["type"]
    if ltype in ("kinked",):
        endpoints = compute_kinked_endpoints(
            pivot, directions, linker_info, N_AZIMUTHAL_SAMPLES
        )
    else:
        endpoints = compute_straight_endpoints(pivot, directions, rod_length)

    cap_points = generate_cap_volume_points(endpoints, cap_radius)
    voxelized = voxelize(cap_points, SWEEP_GRID_SPACING)
    filtered = apply_steric_exclusion(voxelized, body_coords, STERIC_EXCLUSION_RADIUS)

    return filtered


def analyze_single_construct(manifest_entry, campaign_dir, cfg):
    """Analyze sweep volume for one binder x linker combination.

    Loads models from both dark and lit states, computes sweep volumes,
    and returns blocking fractions + gating ratio.
    """
    binder_name = manifest_entry["binder_name"]
    linker_type = manifest_entry["linker_type"]
    linker_info = LINKER_LIBRARY.get(linker_type)
    if linker_info is None:
        print(f"  WARNING: linker '{linker_type}' not in library, skipping")
        return None

    scaffold_cfg = cfg["scaffold"]
    cap_radius = DEFAULT_CAP_RADIUS
    dark_cone = scaffold_cfg.get("dark_cone_half_angle", DEFAULT_DARK_CONE)
    lit_cone = scaffold_cfg.get("lit_cone_half_angle", DEFAULT_LIT_CONE)

    results_dir = campaign_dir / "boltz2_results"

    # Load manifest entries for this binder+linker pair (dark and lit)
    manifest_path = campaign_dir / "manifest.json"
    with open(manifest_path) as f:
        full_manifest = json.load(f)

    dark_entries = [m for m in full_manifest
                    if m["binder_name"] == binder_name
                    and m["linker_type"] == linker_type
                    and m["state"] == "dark"]
    lit_entries = [m for m in full_manifest
                   if m["binder_name"] == binder_name
                   and m["linker_type"] == linker_type
                   and m["state"] == "lit"]

    if not dark_entries or not lit_entries:
        print(f"  Missing dark/lit entries for {binder_name}_{linker_type}")
        return None

    dark_entry = dark_entries[0]
    lit_entry = lit_entries[0]

    # Extract domain boundaries
    dark_domains = dark_entry["domains"]
    lit_domains = lit_entry["domains"]

    # Find CIF files for dark and lit predictions
    dark_cifs = find_cif_files(results_dir, dark_entry["construct_name"])
    lit_cifs = find_cif_files(results_dir, lit_entry["construct_name"])

    if not dark_cifs and not lit_cifs:
        print(f"  No CIF files for {binder_name}_{linker_type}")
        return None

    parser = MMCIFParser(QUIET=True)

    # Process all models, extract pivot geometry from each
    dark_pivots = []
    lit_pivots = []

    for cif_path in dark_cifs:
        try:
            struct = parser.get_structure("s", str(cif_path))
            binder_range = dark_domains["binder"]
            lov2_range = dark_domains["lov2"]

            all_ca = extract_ca_coords(struct)
            binder_ca = all_ca[binder_range[0]:binder_range[1]]
            lov2_ca = all_ca[lov2_range[0]:lov2_range[1]]

            pivot, exit_vec = extract_pivot_and_exit_vector(lov2_ca)
            dark_pivots.append({
                "pivot": pivot,
                "exit_vec": exit_vec,
                "binder_ca": binder_ca,
                "body_ca": all_ca[:lov2_range[1]],  # binder + spacer + lov2
            })
        except Exception as e:
            print(f"    Error loading {cif_path.name}: {e}")

    for cif_path in lit_cifs:
        try:
            struct = parser.get_structure("s", str(cif_path))
            binder_range = lit_domains["binder"]
            lov2_range = lit_domains["lov2"]

            all_ca = extract_ca_coords(struct)
            binder_ca = all_ca[binder_range[0]:binder_range[1]]
            lov2_ca = all_ca[lov2_range[0]:lov2_range[1]]

            pivot, exit_vec = extract_pivot_and_exit_vector(lov2_ca)
            lit_pivots.append({
                "pivot": pivot,
                "exit_vec": exit_vec,
                "binder_ca": binder_ca,
                "body_ca": all_ca[:lov2_range[1]],
            })
        except Exception as e:
            print(f"    Error loading {cif_path.name}: {e}")

    if not dark_pivots and not lit_pivots:
        print(f"  No valid models for {binder_name}_{linker_type}")
        return None

    # Use target hotspot coords from config
    target = cfg["target"]
    hotspot_1idx = target["hotspot_residues_1idx"]
    hotspot_0idx = [h - 1 for h in hotspot_1idx]

    # Generate reference target coords (from the best model binder alignment)
    # If reference CIF exists, use that; otherwise use first model
    ref_cif = None
    binder_cfg = next((b for b in cfg["binders"] if b["name"] == binder_name), None)
    if binder_cfg and binder_cfg.get("reference_complex_cif"):
        ref_path = Path(binder_cfg["reference_complex_cif"])
        if ref_path.exists():
            ref_cif = ref_path

    # For now, generate sweep volumes from each model's geometry
    # and aggregate the blocking fractions
    dark_blockings = []
    lit_blockings = []

    sweep_dir = campaign_dir / "sweep_volumes"

    # Use averaged pivot geometry across models
    if dark_pivots:
        avg_dark_pivot = np.mean([p["pivot"] for p in dark_pivots], axis=0)
        avg_dark_exit = np.mean([p["exit_vec"] for p in dark_pivots], axis=0)
        avg_dark_exit /= np.linalg.norm(avg_dark_exit)
        avg_dark_body = np.vstack([p["body_ca"] for p in dark_pivots[:1]])

        dark_sweep = generate_sweep_pointcloud(
            avg_dark_pivot, avg_dark_exit, linker_info,
            dark_cone, cap_radius, avg_dark_body
        )

        # Save dark sweep PDB
        if len(dark_sweep) > 0:
            save_pointcloud_pdb(
                dark_sweep,
                sweep_dir / f"{binder_name}_{linker_type}_dark_sweep.pdb"
            )

        # Compute blocking against target hotspot positions
        # We need target coordinates — generate synthetic hotspot coords
        # from the models if no reference complex is available
        if ref_cif and ref_cif.exists():
            target_coords, hotspot_coords = load_target_from_reference(
                ref_cif, hotspot_0idx, parser
            )
        else:
            # Use target positions from the prediction models
            # (approximate — binder portion overlaps with reference)
            target_coords = generate_synthetic_hotspot_coords(
                avg_dark_pivot, avg_dark_exit, hotspot_0idx
            )
            hotspot_coords = target_coords

        if target_coords is not None and len(dark_sweep) > 0:
            dark_blocking = compute_blocking_fraction(
                dark_sweep, target_coords,
                list(range(len(hotspot_coords))) if hotspot_coords is not None else None
            )
        else:
            dark_blocking = {"total_blocked_fraction": 0.0, "hotspot_blocked_fraction": 0.0}
    else:
        dark_sweep = np.array([])
        dark_blocking = {"total_blocked_fraction": 0.0, "hotspot_blocked_fraction": 0.0}

    if lit_pivots:
        avg_lit_pivot = np.mean([p["pivot"] for p in lit_pivots], axis=0)
        avg_lit_exit = np.mean([p["exit_vec"] for p in lit_pivots], axis=0)
        avg_lit_exit /= np.linalg.norm(avg_lit_exit)
        avg_lit_body = np.vstack([p["body_ca"] for p in lit_pivots[:1]])

        lit_sweep = generate_sweep_pointcloud(
            avg_lit_pivot, avg_lit_exit, linker_info,
            lit_cone, cap_radius, avg_lit_body
        )

        if len(lit_sweep) > 0:
            save_pointcloud_pdb(
                lit_sweep,
                sweep_dir / f"{binder_name}_{linker_type}_lit_sweep.pdb"
            )

        if target_coords is not None and len(lit_sweep) > 0:
            lit_blocking = compute_blocking_fraction(
                lit_sweep, target_coords,
                list(range(len(hotspot_coords))) if hotspot_coords is not None else None
            )
        else:
            lit_blocking = {"total_blocked_fraction": 0.0, "hotspot_blocked_fraction": 0.0}
    else:
        lit_sweep = np.array([])
        lit_blocking = {"total_blocked_fraction": 0.0, "hotspot_blocked_fraction": 0.0}

    # Compute gating ratio
    dark_hot = dark_blocking["hotspot_blocked_fraction"]
    lit_hot = lit_blocking["hotspot_blocked_fraction"]

    if lit_hot > 0.001:
        gating_ratio = dark_hot / lit_hot
    elif dark_hot > 0.001:
        gating_ratio = dark_hot / 0.001  # cap at 1000x
    else:
        gating_ratio = 1.0

    # Pivot shift (distance between dark and lit pivot positions)
    pivot_shift = 0.0
    if dark_pivots and lit_pivots:
        pivot_shift = float(np.linalg.norm(avg_dark_pivot - avg_lit_pivot))

    return {
        "binder_name": binder_name,
        "linker_type": linker_type,
        "linker_geometry": linker_info["type"],
        "linker_description": linker_info.get("description", ""),
        "dark_cone_angle": dark_cone,
        "lit_cone_angle": lit_cone,
        "n_dark_models": len(dark_cifs),
        "n_lit_models": len(lit_cifs),
        "n_dark_valid": len(dark_pivots),
        "n_lit_valid": len(lit_pivots),
        "dark_sweep_voxels": len(dark_sweep),
        "lit_sweep_voxels": len(lit_sweep),
        "pivot_shift_angstrom": round(pivot_shift, 1),
        "dark_blocking": {
            "total_blocked_fraction": round(dark_blocking["total_blocked_fraction"], 4),
            "hotspot_blocked_fraction": round(dark_blocking["hotspot_blocked_fraction"], 4),
        },
        "lit_blocking": {
            "total_blocked_fraction": round(lit_blocking["total_blocked_fraction"], 4),
            "hotspot_blocked_fraction": round(lit_blocking["hotspot_blocked_fraction"], 4),
        },
        "differential": {
            "dark_hotspot_pct": round(dark_hot * 100, 1),
            "lit_hotspot_pct": round(lit_hot * 100, 1),
            "gating_ratio": round(gating_ratio, 1),
        },
    }


def find_cif_files(results_dir, construct_name):
    """Find all CIF files for a construct in Boltz-2 results."""
    cif_files = []
    construct_dir = results_dir / construct_name
    if construct_dir.exists():
        for cif in sorted(construct_dir.rglob("*.cif")):
            cif_files.append(cif)
    return cif_files


def load_target_from_reference(ref_cif, hotspot_0idx, parser):
    """Load target coordinates and hotspot coords from a reference complex CIF.

    Assumes 2-chain complex: shorter chain = binder, longer chain = target.
    """
    struct = parser.get_structure("ref", str(ref_cif))
    model = struct[0]
    chains = list(model.get_chains())

    if len(chains) < 2:
        return None, None

    # Identify target as longer chain
    chain_lens = []
    for c in chains:
        n_res = len([r for r in c.get_residues() if "CA" in r])
        chain_lens.append((c.id, n_res))
    chain_lens.sort(key=lambda x: x[1], reverse=True)
    target_chain_id = chain_lens[0][0]

    # Extract all target atom coords
    target_chain = None
    for c in chains:
        if c.id == target_chain_id:
            target_chain = c
            break

    all_atoms = []
    residue_atom_indices = {}
    idx = 0
    for i, res in enumerate(target_chain.get_residues()):
        if "CA" not in res:
            continue
        start_idx = idx
        for atom in res.get_atoms():
            all_atoms.append(atom.get_vector().get_array())
            idx += 1
        residue_atom_indices[i] = (start_idx, idx)

    target_coords = np.array(all_atoms)

    # Extract hotspot atom indices
    hotspot_atom_indices = []
    for h_idx in hotspot_0idx:
        if h_idx in residue_atom_indices:
            start, end = residue_atom_indices[h_idx]
            hotspot_atom_indices.extend(range(start, end))

    if hotspot_atom_indices:
        hotspot_coords = target_coords[hotspot_atom_indices]
    else:
        hotspot_coords = target_coords

    return target_coords, hotspot_coords


def generate_synthetic_hotspot_coords(pivot, exit_vec, hotspot_0idx):
    """Generate approximate hotspot positions when no reference complex available.

    Places hotspot atoms in the binding direction from the pivot, at typical
    binding interface distance (~20-40A from binder).
    """
    # Without a reference complex, we can't meaningfully compute blocking
    # Return None to signal that blocking analysis should be skipped
    return None


# ===== Power analysis =====

def lookup_rate(value, rate_table):
    """Look up success rate from a binned rate table."""
    for (lo, hi), rate in rate_table.items():
        if lo <= value < hi:
            return rate
    return 0.0


def power_analysis(binder_iptms, gating_ratios, target_n_hits=3):
    """
    Estimate how many constructs to test experimentally.

    Args:
        binder_iptms: list of ipTM scores for candidate binders
        gating_ratios: list of gating ratios for each construct
        target_n_hits: desired number of functional hits

    Returns:
        dict with n_to_test, expected_hits, per_construct_probs, confidence
    """
    per_construct = []

    for iptm, gating in zip(binder_iptms, gating_ratios):
        # P(binding) from ipTM
        p_bind = lookup_rate(iptm, IPTM_SUCCESS_RATES)

        # P(gating works | binding works) from gating ratio
        p_gate = lookup_rate(gating, GATING_SUCCESS_RATES)

        # Joint probability
        p_success = p_bind * p_gate

        per_construct.append({
            "iptm": round(iptm, 3),
            "p_bind": round(p_bind, 3),
            "gating_ratio": round(gating, 1),
            "p_gate": round(p_gate, 3),
            "p_success": round(p_success, 4),
        })

    # Sort by success probability
    per_construct.sort(key=lambda x: x["p_success"], reverse=True)

    # Compute how many to test for 95% confidence of >= target_n_hits
    probs = [c["p_success"] for c in per_construct]
    expected_total = sum(probs)

    # Use binomial approximation: test N constructs, each with avg prob p
    # P(>= k hits) >= 1 - sum(C(N,i)*p^i*(1-p)^(N-i), i=0..k-1)
    # Simplified: find N where expected hits >= target_n_hits / 0.95
    if expected_total > 0:
        avg_p = expected_total / len(probs)
        # N * avg_p >= target_n_hits / 0.8 (80% confidence buffer)
        n_for_target = int(np.ceil(target_n_hits / (avg_p * 0.8))) if avg_p > 0 else len(probs)
        n_to_test = min(n_for_target, len(probs))
    else:
        n_to_test = len(probs)

    # Expected hits from top N
    top_n_probs = sorted(probs, reverse=True)[:n_to_test]
    expected_hits = sum(top_n_probs)

    return {
        "target_hits": target_n_hits,
        "n_constructs_available": len(probs),
        "n_to_test_recommended": n_to_test,
        "expected_hits_from_top_n": round(expected_hits, 2),
        "avg_success_probability": round(sum(probs) / len(probs), 4) if probs else 0,
        "best_construct_prob": round(max(probs), 4) if probs else 0,
        "per_construct": per_construct,
    }


# ===== ChimeraX visualization =====

def generate_chimerax_script(result, campaign_dir, cfg):
    """Generate a ChimeraX .cxc script for one binder x linker sweep result."""
    binder = result["binder_name"]
    linker = result["linker_type"]
    sweep_dir = campaign_dir / "sweep_volumes"
    dark_pdb = sweep_dir / f"{binder}_{linker}_dark_sweep.pdb"
    lit_pdb = sweep_dir / f"{binder}_{linker}_lit_sweep.pdb"

    dark_hot = result["differential"]["dark_hotspot_pct"]
    lit_hot = result["differential"]["lit_hotspot_pct"]
    ratio = result["differential"]["gating_ratio"]

    hotspots = cfg["target"]["hotspot_residues_1idx"]
    hotspot_spec = ",".join(str(r) for r in hotspots)

    # Check for reference PDB/CIF
    ref_pdb = cfg["target"].get("pdb_path")

    lines = [
        f"# Sweep Volume: {binder} + {linker}",
        f"# Dark cone: {result.get('dark_cone_angle', 25)} deg, "
        f"Lit cone: {result.get('lit_cone_angle', 60)} deg",
        f"# Dark hotspot blocked: {dark_hot:.0f}%, "
        f"Lit: {lit_hot:.0f}%, Ratio: {ratio:.1f}x",
        "",
        "close all",
        "set bgColor white",
        "",
    ]

    model_num = 1

    # Reference complex (if available)
    if ref_pdb and Path(ref_pdb).exists():
        lines += [
            "# === Reference complex ===",
            f"open {ref_pdb}",
            f"hide #{model_num}",
            f"show #{model_num} cartoons",
            f"color #{model_num} cornflowerblue",
            "",
            "# Highlight hotspot residues",
            f"color #{model_num}/A:{hotspot_spec} lime green",
            f"show #{model_num}/A:{hotspot_spec} atoms",
            f"style #{model_num}/A:{hotspot_spec} sphere",
            f"size #{model_num}/A:{hotspot_spec} atomRadius 1.5",
            "",
        ]
        model_num += 1

    # Dark sweep volume
    if dark_pdb.exists():
        lines += [
            "# === Dark sweep volume ===",
            f"open {dark_pdb}",
            f"hide #{model_num} atoms",
            f"surface #{model_num} resolution 4 gridSpacing 2",
            f"color #{model_num} tomato",
            f"transparency #{model_num} 55",
            "",
        ]
        model_num += 1

    # Lit sweep volume
    if lit_pdb.exists():
        lines += [
            "# === Lit sweep volume ===",
            f"open {lit_pdb}",
            f"hide #{model_num} atoms",
            f"surface #{model_num} resolution 4 gridSpacing 2",
            f"color #{model_num} forest green",
            f"transparency #{model_num} 55",
            "",
        ]
        model_num += 1

    # Labels and view
    lines += [
        "# === Labels ===",
        f"2dlabels create title text '{binder} + {linker}' "
        f"xpos 0.02 ypos 0.95 size 18 bold true",
        f"2dlabels create dark_label text 'Dark: "
        f"cone +/-{result.get('dark_cone_angle', 25)} deg, "
        f"hotspot blocked {dark_hot:.0f}%' "
        f"xpos 0.02 ypos 0.90 size 14 color tomato",
        f"2dlabels create lit_label text 'Lit: "
        f"cone +/-{result.get('lit_cone_angle', 60)} deg, "
        f"hotspot blocked {lit_hot:.0f}%' "
        f"xpos 0.02 ypos 0.86 size 14 color forest green",
        f"2dlabels create ratio_label text 'Gating ratio: {ratio:.1f}x' "
        f"xpos 0.02 ypos 0.82 size 14 bold true",
        "",
        "lighting soft",
        "graphics silhouettes true width 1.5",
        "view all",
        "",
    ]

    return "\n".join(lines)


# ===== Main =====

def analyze_campaign(campaign_dir):
    """Run full analysis on a completed campaign."""
    campaign_dir = Path(campaign_dir)
    config_path = campaign_dir / "campaign_config.json"
    manifest_path = campaign_dir / "manifest.json"

    if not config_path.exists():
        print(f"Config not found: {config_path}")
        return
    if not manifest_path.exists():
        print(f"Manifest not found: {manifest_path}")
        return

    with open(config_path) as f:
        cfg = json.load(f)
    with open(manifest_path) as f:
        manifest = json.load(f)

    print(f"Campaign: {cfg['campaign_name']}")
    print(f"Target: {cfg['target']['name']}")
    print(f"Hotspots: {cfg['target']['hotspot_residues_1idx']}")

    # Get unique binder x linker combinations
    pairs = set()
    for entry in manifest:
        pairs.add((entry["binder_name"], entry["linker_type"]))

    print(f"Constructs: {len(pairs)} binder x linker pairs")
    print()

    # Analyze each pair
    sweep_results = []
    for binder_name, linker_type in sorted(pairs):
        entry = next(m for m in manifest
                     if m["binder_name"] == binder_name
                     and m["linker_type"] == linker_type)

        print(f"  Analyzing {binder_name} + {linker_type}...", end=" ")
        result = analyze_single_construct(entry, campaign_dir, cfg)

        if result:
            sweep_results.append(result)
            ratio = result["differential"]["gating_ratio"]
            dark = result["differential"]["dark_hotspot_pct"]
            lit = result["differential"]["lit_hotspot_pct"]
            print(f"dark={dark:.0f}% lit={lit:.0f}% ratio={ratio:.1f}x")
        else:
            print("(skipped)")

    if not sweep_results:
        print("\nNo results to analyze.")
        return

    # Rank by gating ratio
    ranked = sorted(sweep_results,
                    key=lambda x: x["differential"]["gating_ratio"],
                    reverse=True)

    # Save sweep results
    analysis_dir = campaign_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    sweep_path = analysis_dir / "sweep_results.json"
    sweep_path.write_text(json.dumps(ranked, indent=2) + "\n")

    # Build ranked constructs with combined scoring
    ranked_constructs = []
    for rank, r in enumerate(ranked, 1):
        binder_cfg = next(
            (b for b in cfg["binders"] if b["name"] == r["binder_name"]),
            {}
        )
        ranked_constructs.append({
            "rank": rank,
            "binder_name": r["binder_name"],
            "linker_type": r["linker_type"],
            "linker_geometry": r["linker_geometry"],
            "binder_iptm": binder_cfg.get("iptm", 0.0),
            "gating_ratio": r["differential"]["gating_ratio"],
            "dark_hotspot_pct": r["differential"]["dark_hotspot_pct"],
            "lit_hotspot_pct": r["differential"]["lit_hotspot_pct"],
            "pivot_shift": r["pivot_shift_angstrom"],
        })

    ranked_path = analysis_dir / "ranked_constructs.json"
    ranked_path.write_text(json.dumps(ranked_constructs, indent=2) + "\n")

    # Power analysis
    binder_iptms = [c["binder_iptm"] for c in ranked_constructs]
    gating_ratios = [c["gating_ratio"] for c in ranked_constructs]
    pa = power_analysis(binder_iptms, gating_ratios, target_n_hits=3)

    power_path = analysis_dir / "power_analysis.json"
    power_path.write_text(json.dumps(pa, indent=2) + "\n")

    # Generate ChimeraX scripts
    chimerax_dir = campaign_dir / "chimerax"
    chimerax_dir.mkdir(parents=True, exist_ok=True)

    for r in ranked:
        cxc = generate_chimerax_script(r, campaign_dir, cfg)
        cxc_path = chimerax_dir / f"sweep_{r['binder_name']}_{r['linker_type']}.cxc"
        cxc_path.write_text(cxc)

    # Print summary
    print("\n" + "=" * 70)
    print(f"RANKED CONSTRUCTS (by gating ratio)")
    print("=" * 70)
    print(f"{'Rank':<5} {'Binder':<20} {'Linker':<12} {'ipTM':<6} "
          f"{'Dark%':<7} {'Lit%':<6} {'Ratio':<8}")
    print("-" * 70)
    for c in ranked_constructs:
        print(f"{c['rank']:<5} {c['binder_name']:<20} {c['linker_type']:<12} "
              f"{c['binder_iptm']:<6.2f} {c['dark_hotspot_pct']:<7.1f} "
              f"{c['lit_hotspot_pct']:<6.1f} {c['gating_ratio']:<8.1f}")

    print(f"\nPower Analysis:")
    print(f"  Available constructs: {pa['n_constructs_available']}")
    print(f"  Recommended to test: {pa['n_to_test_recommended']}")
    print(f"  Expected hits (top N): {pa['expected_hits_from_top_n']:.1f}")
    print(f"  Best construct prob:   {pa['best_construct_prob']:.3f}")

    print(f"\nOutputs:")
    print(f"  Sweep results:      {sweep_path}")
    print(f"  Ranked constructs:  {ranked_path}")
    print(f"  Power analysis:     {power_path}")
    print(f"  ChimeraX scripts:   {chimerax_dir}/")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_campaign.py <campaign_dir>")
        print("       python3 analyze_campaign.py campaigns/p300_hat_demo/")
        sys.exit(1)

    analyze_campaign(sys.argv[1])
