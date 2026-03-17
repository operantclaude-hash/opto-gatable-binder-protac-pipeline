"""Tests for analyze_campaign.py — geometry, sweep volume, blocking, and power analysis."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from analyze_campaign import (
    extract_pivot_and_exit_vector,
    sample_cone_directions,
    compute_rod_length,
    compute_straight_endpoints,
    compute_kinked_endpoints,
    generate_cap_volume_points,
    voxelize,
    apply_steric_exclusion,
    compute_blocking_fraction,
    superimpose_binder,
    apply_transform,
    lookup_rate,
    power_analysis,
)
from linker_library import IPTM_SUCCESS_RATES, GATING_SUCCESS_RATES


class TestExtractPivotAndExitVector(unittest.TestCase):

    def test_line_along_z(self):
        coords = np.array([[0, 0, i] for i in range(10)], dtype=float)
        pivot, exit_vec = extract_pivot_and_exit_vector(coords)
        np.testing.assert_array_equal(pivot, [0, 0, 9])
        self.assertGreater(exit_vec[2], 0.9)

    def test_line_along_x(self):
        coords = np.array([[i, 0, 0] for i in range(10)], dtype=float)
        pivot, exit_vec = extract_pivot_and_exit_vector(coords)
        np.testing.assert_array_equal(pivot, [9, 0, 0])
        self.assertGreater(exit_vec[0], 0.9)

    def test_pivot_is_last_ca(self):
        coords = np.random.RandomState(42).randn(20, 3)
        pivot, _ = extract_pivot_and_exit_vector(coords)
        np.testing.assert_array_equal(pivot, coords[-1])

    def test_exit_vec_is_unit(self):
        coords = np.random.RandomState(42).randn(15, 3)
        _, exit_vec = extract_pivot_and_exit_vector(coords)
        self.assertAlmostEqual(np.linalg.norm(exit_vec), 1.0, places=6)

    def test_degenerate_same_points(self):
        coords = np.array([[5, 5, 5]] * 5, dtype=float)
        _, exit_vec = extract_pivot_and_exit_vector(coords)
        np.testing.assert_array_almost_equal(exit_vec, [0, 0, 1])

    def test_short_coords(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
        pivot, exit_vec = extract_pivot_and_exit_vector(coords)
        np.testing.assert_array_equal(pivot, [2, 0, 0])
        self.assertAlmostEqual(np.linalg.norm(exit_vec), 1.0, places=6)


class TestSampleConeDirections(unittest.TestCase):

    def test_output_count(self):
        dirs = sample_cone_directions(np.array([0, 0, 1.0]), 25.0, n_polar=10, n_azimuthal=36)
        self.assertEqual(len(dirs), 11 * 36)  # (n_polar + 1) * n_azimuthal

    def test_all_unit_vectors(self):
        dirs = sample_cone_directions(np.array([1, 0, 0.0]), 30.0)
        norms = np.linalg.norm(dirs, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_all_within_cone_angle(self):
        axis = np.array([0, 0, 1.0])
        half_angle = 25.0
        dirs = sample_cone_directions(axis, half_angle, n_polar=10, n_azimuthal=36)
        cos_angles = dirs @ axis
        angles_deg = np.degrees(np.arccos(np.clip(cos_angles, -1, 1)))
        self.assertTrue(np.all(angles_deg <= half_angle + 0.5))

    def test_zero_cone_returns_axis_only(self):
        axis = np.array([0, 1, 0.0])
        dirs = sample_cone_directions(axis, 0.0, n_polar=5, n_azimuthal=12)
        for d in dirs:
            np.testing.assert_allclose(np.abs(d), np.abs(axis), atol=1e-6)

    def test_different_axis_orientations(self):
        for axis in [np.array([1, 0, 0.0]), np.array([0, 1, 0.0]),
                     np.array([0, 0, 1.0]), np.array([1, 1, 1.0]) / np.sqrt(3)]:
            dirs = sample_cone_directions(axis, 30.0, n_polar=5, n_azimuthal=12)
            norms = np.linalg.norm(dirs, axis=1)
            np.testing.assert_allclose(norms, 1.0, atol=1e-6)


class TestComputeRodLength(unittest.TestCase):

    def test_straight_helix(self):
        info = {"type": "straight", "length": 15}
        self.assertAlmostEqual(compute_rod_length(info), 15 * 1.5)

    def test_ppII(self):
        info = {"type": "ppII", "length": 16}
        self.assertAlmostEqual(compute_rod_length(info), 16 * 3.12)

    def test_kinked_uses_helix_rise(self):
        info = {"type": "kinked", "length": 22}
        self.assertAlmostEqual(compute_rod_length(info), 22 * 1.5)


class TestComputeStraightEndpoints(unittest.TestCase):

    def test_single_direction(self):
        pivot = np.array([0, 0, 0.0])
        dirs = np.array([[0, 0, 1.0]])
        endpoints = compute_straight_endpoints(pivot, dirs, 10.0)
        np.testing.assert_allclose(endpoints[0], [0, 0, 10])

    def test_multiple_directions(self):
        pivot = np.array([5, 5, 5.0])
        dirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1.0]])
        endpoints = compute_straight_endpoints(pivot, dirs, 20.0)
        self.assertEqual(len(endpoints), 3)
        np.testing.assert_allclose(endpoints[0], [25, 5, 5])
        np.testing.assert_allclose(endpoints[1], [5, 25, 5])
        np.testing.assert_allclose(endpoints[2], [5, 5, 25])

    def test_offset_pivot(self):
        pivot = np.array([10, 20, 30.0])
        dirs = np.array([[0, 0, 1.0]])
        endpoints = compute_straight_endpoints(pivot, dirs, 5.0)
        np.testing.assert_allclose(endpoints[0], [10, 20, 35])


class TestComputeKinkedEndpoints(unittest.TestCase):

    def test_output_count(self):
        pivot = np.zeros(3)
        dirs = np.array([[0, 0, 1.0], [1, 0, 0.0]])
        info = {"pre_helix": 10, "post_helix": 10, "kink_res": 2,
                "description": "proline kink", "multi_kink": False}
        endpoints = compute_kinked_endpoints(pivot, dirs, info, n_azimuthal=12)
        self.assertEqual(len(endpoints), 2 * 12)

    def test_endpoints_beyond_pre_length(self):
        pivot = np.zeros(3)
        dirs = np.array([[0, 0, 1.0]])
        info = {"pre_helix": 10, "post_helix": 10, "kink_res": 2,
                "description": "proline kink", "multi_kink": False}
        endpoints = compute_kinked_endpoints(pivot, dirs, info, n_azimuthal=36)
        pre_len = 10 * 1.5
        dists = np.linalg.norm(endpoints - pivot, axis=1)
        self.assertTrue(np.all(dists >= pre_len * 0.3))  # tolerance for kink deflection

    def test_beta_turn_kink_angle(self):
        """Beta turn should produce larger deflection than proline."""
        pivot = np.zeros(3)
        dirs = np.array([[0, 0, 1.0]])
        pro_info = {"pre_helix": 10, "post_helix": 10, "kink_res": 2,
                    "description": "proline kink", "multi_kink": False}
        beta_info = {"pre_helix": 10, "post_helix": 10, "kink_res": 4,
                     "description": "beta turn", "multi_kink": False}
        pro_endpoints = compute_kinked_endpoints(pivot, dirs, pro_info, n_azimuthal=36)
        beta_endpoints = compute_kinked_endpoints(pivot, dirs, beta_info, n_azimuthal=36)
        # Beta turn endpoints should spread more in XY plane
        pro_xy = np.mean(np.sqrt(pro_endpoints[:, 0]**2 + pro_endpoints[:, 1]**2))
        beta_xy = np.mean(np.sqrt(beta_endpoints[:, 0]**2 + beta_endpoints[:, 1]**2))
        self.assertGreater(beta_xy, pro_xy)


class TestGenerateCapVolumePoints(unittest.TestCase):

    def test_single_endpoint_27_points(self):
        endpoints = np.array([[0, 0, 0.0]])
        pts = generate_cap_volume_points(endpoints, 14.0)
        self.assertEqual(len(pts), 27)

    def test_multiple_endpoints(self):
        endpoints = np.array([[0, 0, 0], [100, 0, 0.0]])
        pts = generate_cap_volume_points(endpoints, 14.0)
        self.assertEqual(len(pts), 54)

    def test_center_is_endpoint(self):
        endpoints = np.array([[5, 5, 5.0]])
        pts = generate_cap_volume_points(endpoints, 14.0)
        np.testing.assert_allclose(pts[0], [5, 5, 5])

    def test_shell_at_cap_radius(self):
        endpoints = np.array([[0, 0, 0.0]])
        cap_r = 14.0
        pts = generate_cap_volume_points(endpoints, cap_r)
        dists = np.linalg.norm(pts[1:] - pts[0], axis=1)
        np.testing.assert_allclose(dists, cap_r, atol=0.01)


class TestVoxelize(unittest.TestCase):

    def test_deduplication(self):
        pts = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [1.0, 1.0, 1.0]], dtype=float)
        result = voxelize(pts, grid_spacing=2.0)
        self.assertEqual(len(result), 1)

    def test_separate_voxels(self):
        pts = np.array([[0, 0, 0], [10, 10, 10.0]])
        result = voxelize(pts, grid_spacing=2.0)
        self.assertEqual(len(result), 2)

    def test_empty_input(self):
        pts = np.array([]).reshape(0, 3)
        result = voxelize(pts, 2.0)
        self.assertEqual(len(result), 0)

    def test_grid_alignment(self):
        pts = np.array([[3.0, 0, 0], [3.9, 0, 0]], dtype=float)
        result = voxelize(pts, 2.0)
        self.assertEqual(len(result), 1)
        np.testing.assert_allclose(result[0], [4.0, 0, 0])


class TestApplyStericExclusion(unittest.TestCase):

    def test_removes_overlapping(self):
        sweep = np.array([[0, 0, 0], [10, 0, 0], [20, 0, 0.0]])
        body = np.array([[0, 0, 0.0]])
        result = apply_steric_exclusion(sweep, body, 3.0)
        self.assertEqual(len(result), 2)

    def test_no_overlap_returns_all(self):
        sweep = np.array([[100, 100, 100.0]])
        body = np.array([[0, 0, 0.0]])
        result = apply_steric_exclusion(sweep, body, 3.0)
        self.assertEqual(len(result), 1)

    def test_empty_sweep(self):
        sweep = np.array([]).reshape(0, 3)
        body = np.array([[0, 0, 0.0]])
        result = apply_steric_exclusion(sweep, body, 3.0)
        self.assertEqual(len(result), 0)

    def test_empty_body(self):
        sweep = np.array([[0, 0, 0], [1, 1, 1.0]])
        body = np.array([]).reshape(0, 3)
        result = apply_steric_exclusion(sweep, body, 3.0)
        self.assertEqual(len(result), 2)


class TestComputeBlockingFraction(unittest.TestCase):

    def test_full_overlap(self):
        sweep = np.array([[0, 0, 0], [2, 0, 0], [4, 0, 0.0]])
        target = np.array([[0, 0, 0], [2, 0, 0.0]])
        result = compute_blocking_fraction(sweep, target)
        self.assertAlmostEqual(result["total_blocked_fraction"], 1.0)

    def test_no_overlap(self):
        sweep = np.array([[0, 0, 0.0]])
        target = np.array([[100, 100, 100.0]])
        result = compute_blocking_fraction(sweep, target)
        self.assertAlmostEqual(result["total_blocked_fraction"], 0.0)

    def test_empty_sweep(self):
        sweep = np.array([]).reshape(0, 3)
        target = np.array([[0, 0, 0.0]])
        result = compute_blocking_fraction(sweep, target)
        self.assertEqual(result["total_blocked_fraction"], 0.0)
        self.assertEqual(result["hotspot_blocked_fraction"], 0.0)

    def test_hotspot_subset(self):
        sweep = np.array([[0, 0, 0], [2, 0, 0.0]])
        target = np.array([[0, 0, 0], [2, 0, 0], [100, 0, 0.0]])
        result = compute_blocking_fraction(sweep, target, hotspot_indices=[0, 1])
        self.assertAlmostEqual(result["hotspot_blocked_fraction"], 1.0)
        self.assertAlmostEqual(result["total_blocked_fraction"], 2.0 / 3.0, places=1)

    def test_hotspot_out_of_range(self):
        sweep = np.array([[0, 0, 0.0]])
        target = np.array([[0, 0, 0.0]])
        result = compute_blocking_fraction(sweep, target, hotspot_indices=[5])
        self.assertAlmostEqual(result["hotspot_blocked_fraction"], 0.0)


class TestSuperimposeBinder(unittest.TestCase):

    def test_identity_transform(self):
        coords = np.random.RandomState(42).randn(20, 3) * 10
        rot, trans = superimpose_binder(coords.copy(), coords.copy())
        np.testing.assert_allclose(rot, np.eye(3), atol=1e-3)
        np.testing.assert_allclose(trans, np.zeros(3), atol=1e-3)

    def test_known_translation(self):
        ref = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0.0]])
        query = ref + np.array([5, 0, 0])
        rotran = superimpose_binder(query, ref)
        transformed = apply_transform(query, rotran)
        np.testing.assert_allclose(transformed, ref, atol=1e-3)


class TestApplyTransform(unittest.TestCase):

    def test_identity(self):
        coords = np.random.RandomState(42).randn(10, 3)
        result = apply_transform(coords, (np.eye(3), np.zeros(3)))
        np.testing.assert_allclose(result, coords)

    def test_translation_only(self):
        coords = np.array([[1, 2, 3.0]])
        result = apply_transform(coords, (np.eye(3), np.array([10, 20, 30.0])))
        np.testing.assert_allclose(result, [[11, 22, 33]])


class TestLookupRate(unittest.TestCase):

    def test_high_iptm(self):
        result = lookup_rate(0.95, IPTM_SUCCESS_RATES)
        self.assertEqual(result, 0.70)

    def test_mid_iptm(self):
        result = lookup_rate(0.85, IPTM_SUCCESS_RATES)
        self.assertEqual(result, 0.40)

    def test_low_iptm(self):
        result = lookup_rate(0.3, IPTM_SUCCESS_RATES)
        self.assertEqual(result, 0.03)

    def test_boundary_value(self):
        result = lookup_rate(0.9, IPTM_SUCCESS_RATES)
        self.assertEqual(result, 0.70)

    def test_gating_high(self):
        result = lookup_rate(15.0, GATING_SUCCESS_RATES)
        self.assertEqual(result, 0.80)

    def test_gating_low(self):
        result = lookup_rate(0.5, GATING_SUCCESS_RATES)
        self.assertEqual(result, 0.05)

    def test_below_all_bins(self):
        result = lookup_rate(-1.0, IPTM_SUCCESS_RATES)
        self.assertEqual(result, 0.0)


class TestPowerAnalysis(unittest.TestCase):

    def test_output_structure(self):
        result = power_analysis([0.9, 0.85], [10.0, 5.0], target_n_hits=1)
        self.assertIn("target_hits", result)
        self.assertIn("n_constructs_available", result)
        self.assertIn("n_to_test_recommended", result)
        self.assertIn("expected_hits_from_top_n", result)
        self.assertIn("per_construct", result)

    def test_high_iptm_high_gating(self):
        result = power_analysis([0.95], [50.0], target_n_hits=1)
        self.assertEqual(result["n_constructs_available"], 1)
        self.assertGreater(result["best_construct_prob"], 0.3)

    def test_low_everything(self):
        result = power_analysis([0.3], [0.5], target_n_hits=3)
        self.assertLess(result["best_construct_prob"], 0.05)

    def test_per_construct_sorted(self):
        result = power_analysis([0.5, 0.9, 0.7], [5.0, 10.0, 3.0], target_n_hits=1)
        probs = [c["p_success"] for c in result["per_construct"]]
        self.assertEqual(probs, sorted(probs, reverse=True))

    def test_n_to_test_does_not_exceed_available(self):
        result = power_analysis([0.9, 0.9], [10.0, 10.0], target_n_hits=100)
        self.assertLessEqual(result["n_to_test_recommended"], 2)


if __name__ == "__main__":
    unittest.main()
