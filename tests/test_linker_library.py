"""Tests for linker_library.py — shared constants and linker definitions."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from linker_library import (
    LOV2_DARK, LOV2_LIT, UBQ_SEQ,
    LINKER_LIBRARY, EAAAK,
    HELIX_RISE_PER_RESIDUE, PPII_RISE_PER_RESIDUE,
    DEFAULT_CAP_RADIUS, DEFAULT_DARK_CONE, DEFAULT_LIT_CONE,
    STERIC_EXCLUSION_RADIUS, SWEEP_GRID_SPACING,
    N_AZIMUTHAL_SAMPLES, N_POLAR_SAMPLES,
    IPTM_SUCCESS_RATES, GATING_SUCCESS_RATES,
)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


class TestLinkerSequences(unittest.TestCase):
    """Test core protein sequences."""

    def test_lov2_dark_length(self):
        self.assertEqual(len(LOV2_DARK), 136)

    def test_lov2_lit_length(self):
        self.assertEqual(len(LOV2_LIT), 121)

    def test_ubq_length(self):
        self.assertEqual(len(UBQ_SEQ), 76)

    def test_lov2_lit_is_prefix_of_dark(self):
        self.assertTrue(LOV2_DARK.startswith(LOV2_LIT))

    def test_lov2_jalpha_difference(self):
        """Jalpha helix is 15aa (136 - 121)."""
        self.assertEqual(len(LOV2_DARK) - len(LOV2_LIT), 15)

    def test_lov2_dark_valid_amino_acids(self):
        for i, aa in enumerate(LOV2_DARK):
            self.assertIn(aa, VALID_AA, f"Invalid AA '{aa}' at position {i} in LOV2_DARK")

    def test_lov2_lit_valid_amino_acids(self):
        for i, aa in enumerate(LOV2_LIT):
            self.assertIn(aa, VALID_AA, f"Invalid AA '{aa}' at position {i} in LOV2_LIT")

    def test_ubq_valid_amino_acids(self):
        for i, aa in enumerate(UBQ_SEQ):
            self.assertIn(aa, VALID_AA, f"Invalid AA '{aa}' at position {i} in UBQ_SEQ")


class TestLinkerLibraryIntegrity(unittest.TestCase):
    """Test the 12-linker library for data integrity."""

    def test_library_has_12_entries(self):
        self.assertEqual(len(LINKER_LIBRARY), 12)

    def test_all_linkers_have_required_fields(self):
        required = {"seq", "length", "type", "description"}
        for name, info in LINKER_LIBRARY.items():
            for field in required:
                self.assertIn(field, info, f"Linker '{name}' missing field '{field}'")

    def test_sequence_lengths_match_stated_lengths(self):
        for name, info in LINKER_LIBRARY.items():
            self.assertEqual(
                len(info["seq"]), info["length"],
                f"Linker '{name}': len(seq)={len(info['seq'])} != stated length={info['length']}"
            )

    def test_linker_types_are_valid(self):
        valid_types = {"kinked", "straight", "ppII"}
        for name, info in LINKER_LIBRARY.items():
            self.assertIn(
                info["type"], valid_types,
                f"Linker '{name}' has invalid type '{info['type']}'"
            )

    def test_kinked_linkers_have_geometry_fields(self):
        for name, info in LINKER_LIBRARY.items():
            if info["type"] == "kinked":
                self.assertIn("pre_helix", info, f"Kinked linker '{name}' missing 'pre_helix'")
                self.assertIn("kink_res", info, f"Kinked linker '{name}' missing 'kink_res'")
                self.assertIn("post_helix", info, f"Kinked linker '{name}' missing 'post_helix'")

    def test_kinked_geometry_sums(self):
        """pre_helix + kink_res + post_helix should equal linker length."""
        for name, info in LINKER_LIBRARY.items():
            if info["type"] == "kinked":
                total = info["pre_helix"] + info["kink_res"] + info["post_helix"]
                if info.get("multi_kink"):
                    continue  # EP4 has repeating pattern, different accounting
                self.assertEqual(
                    total, info["length"],
                    f"Linker '{name}': geometry sum {total} != length {info['length']}"
                )

    def test_multi_kink_only_on_EP4(self):
        for name, info in LINKER_LIBRARY.items():
            if name == "EP4":
                self.assertTrue(info.get("multi_kink", False))
            else:
                self.assertFalse(info.get("multi_kink", False),
                                 f"Unexpected multi_kink=True on linker '{name}'")

    def test_eaaak_composition(self):
        """EAAAK3 and EAAAK5 should be exact repeats of EAAAK."""
        self.assertEqual(LINKER_LIBRARY["EAAAK3"]["seq"], EAAAK * 3)
        self.assertEqual(LINKER_LIBRARY["EAAAK5"]["seq"], EAAAK * 5)

    def test_p16_is_pure_proline(self):
        self.assertEqual(LINKER_LIBRARY["P16"]["seq"], "P" * 16)

    def test_all_linker_seqs_valid_amino_acids(self):
        for name, info in LINKER_LIBRARY.items():
            for i, aa in enumerate(info["seq"]):
                self.assertIn(aa, VALID_AA,
                              f"Invalid AA '{aa}' at pos {i} in linker '{name}'")


class TestPhysicalConstants(unittest.TestCase):
    """Test physical constants for sweep analysis."""

    def test_helix_rise_value(self):
        self.assertEqual(HELIX_RISE_PER_RESIDUE, 1.5)

    def test_ppii_rise_value(self):
        self.assertEqual(PPII_RISE_PER_RESIDUE, 3.12)

    def test_default_cap_radius(self):
        self.assertEqual(DEFAULT_CAP_RADIUS, 14.0)

    def test_dark_cone_less_than_lit_cone(self):
        self.assertLess(DEFAULT_DARK_CONE, DEFAULT_LIT_CONE)

    def test_cone_angles_positive_and_reasonable(self):
        self.assertGreater(DEFAULT_DARK_CONE, 0)
        self.assertLess(DEFAULT_DARK_CONE, 180)
        self.assertGreater(DEFAULT_LIT_CONE, 0)
        self.assertLess(DEFAULT_LIT_CONE, 180)

    def test_steric_exclusion_positive(self):
        self.assertGreater(STERIC_EXCLUSION_RADIUS, 0)

    def test_grid_spacing_positive(self):
        self.assertGreater(SWEEP_GRID_SPACING, 0)

    def test_sampling_counts_positive(self):
        self.assertGreater(N_AZIMUTHAL_SAMPLES, 0)
        self.assertGreater(N_POLAR_SAMPLES, 0)


class TestRateTables(unittest.TestCase):
    """Test ipTM and gating success rate lookup tables."""

    def test_iptm_bins_cover_0_to_1(self):
        boundaries = sorted(set(lo for (lo, hi) in IPTM_SUCCESS_RATES))
        self.assertEqual(min(boundaries), 0.0)
        hi_boundaries = sorted(set(hi for (lo, hi) in IPTM_SUCCESS_RATES))
        self.assertEqual(max(hi_boundaries), 1.0)

    def test_iptm_rates_monotonically_increase(self):
        """Higher ipTM bins should have higher success rates."""
        sorted_bins = sorted(IPTM_SUCCESS_RATES.items(), key=lambda x: x[0][0])
        rates = [r for _, r in sorted_bins]
        for i in range(1, len(rates)):
            self.assertGreaterEqual(rates[i], rates[i - 1],
                                    f"ipTM rates not monotonic at bin {i}")

    def test_gating_bins_cover_zero(self):
        lo_bounds = [lo for (lo, hi) in GATING_SUCCESS_RATES]
        self.assertIn(0.0, lo_bounds)

    def test_gating_bins_cover_infinity(self):
        hi_bounds = [hi for (lo, hi) in GATING_SUCCESS_RATES]
        self.assertIn(float("inf"), hi_bounds)

    def test_gating_rates_monotonically_increase(self):
        sorted_bins = sorted(GATING_SUCCESS_RATES.items(), key=lambda x: x[0][0])
        rates = [r for _, r in sorted_bins]
        for i in range(1, len(rates)):
            self.assertGreaterEqual(rates[i], rates[i - 1])

    def test_all_rates_between_0_and_1(self):
        for rate in IPTM_SUCCESS_RATES.values():
            self.assertGreaterEqual(rate, 0.0)
            self.assertLessEqual(rate, 1.0)
        for rate in GATING_SUCCESS_RATES.values():
            self.assertGreaterEqual(rate, 0.0)
            self.assertLessEqual(rate, 1.0)


if __name__ == "__main__":
    unittest.main()
