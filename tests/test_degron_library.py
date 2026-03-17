"""Tests for degron_library.py — VVD, E3 ligase, and GS linker definitions."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from degron_library import (
    VVD_SEQ, VVD_LENGTH,
    E3_LIGASE_LIBRARY,
    PATHWAY_DEFAULTS,
    GS_LINKER_LIBRARY,
    AAV_MAX_BP, AAV_MAX_AA,
)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


class TestVVDSequence(unittest.TestCase):
    """Test VVD dimerizer sequence."""

    def test_vvd_length(self):
        self.assertEqual(len(VVD_SEQ), 158)

    def test_vvd_length_constant(self):
        self.assertEqual(VVD_LENGTH, len(VVD_SEQ))

    def test_vvd_valid_amino_acids(self):
        for i, aa in enumerate(VVD_SEQ):
            self.assertIn(aa, VALID_AA, f"Invalid AA '{aa}' at position {i} in VVD_SEQ")

    def test_vvd_starts_with_met(self):
        self.assertEqual(VVD_SEQ[0], "M")


class TestE3LigaseLibrary(unittest.TestCase):
    """Test E3 ubiquitin ligase entries."""

    def test_has_two_ligases(self):
        self.assertEqual(len(E3_LIGASE_LIBRARY), 2)

    def test_trim21_exists(self):
        self.assertIn("TRIM21_RBCC", E3_LIGASE_LIBRARY)

    def test_nedd4l_exists(self):
        self.assertIn("NEDD4L_HECT", E3_LIGASE_LIBRARY)

    def test_required_fields(self):
        required = {"seq", "length", "pathway", "ubiquitin_chain", "source", "description"}
        for name, info in E3_LIGASE_LIBRARY.items():
            for field in required:
                self.assertIn(field, info, f"E3 '{name}' missing field '{field}'")

    def test_sequence_lengths_match(self):
        for name, info in E3_LIGASE_LIBRARY.items():
            self.assertEqual(
                len(info["seq"]), info["length"],
                f"E3 '{name}': len(seq)={len(info['seq'])} != stated length={info['length']}"
            )

    def test_trim21_length(self):
        self.assertEqual(E3_LIGASE_LIBRARY["TRIM21_RBCC"]["length"], 236)

    def test_nedd4l_length(self):
        self.assertEqual(E3_LIGASE_LIBRARY["NEDD4L_HECT"]["length"], 377)

    def test_trim21_proteasomal(self):
        self.assertEqual(E3_LIGASE_LIBRARY["TRIM21_RBCC"]["pathway"], "proteasomal")
        self.assertEqual(E3_LIGASE_LIBRARY["TRIM21_RBCC"]["ubiquitin_chain"], "K48")

    def test_nedd4l_lysosomal(self):
        self.assertEqual(E3_LIGASE_LIBRARY["NEDD4L_HECT"]["pathway"], "lysosomal")
        self.assertEqual(E3_LIGASE_LIBRARY["NEDD4L_HECT"]["ubiquitin_chain"], "K63")

    def test_all_valid_amino_acids(self):
        for name, info in E3_LIGASE_LIBRARY.items():
            for i, aa in enumerate(info["seq"]):
                self.assertIn(aa, VALID_AA,
                              f"Invalid AA '{aa}' at pos {i} in E3 '{name}'")

    def test_trim21_starts_with_met(self):
        self.assertEqual(E3_LIGASE_LIBRARY["TRIM21_RBCC"]["seq"][0], "M")

    def test_nedd4l_starts_with_met(self):
        self.assertEqual(E3_LIGASE_LIBRARY["NEDD4L_HECT"]["seq"][0], "M")


class TestPathwayDefaults(unittest.TestCase):
    """Test pathway -> E3 ligase mapping."""

    def test_proteasomal_default(self):
        self.assertEqual(PATHWAY_DEFAULTS["proteasomal"], "TRIM21_RBCC")

    def test_lysosomal_default(self):
        self.assertEqual(PATHWAY_DEFAULTS["lysosomal"], "NEDD4L_HECT")

    def test_defaults_reference_valid_ligases(self):
        for pathway, e3_name in PATHWAY_DEFAULTS.items():
            self.assertIn(e3_name, E3_LIGASE_LIBRARY,
                          f"Pathway '{pathway}' default '{e3_name}' not in library")

    def test_two_pathways(self):
        self.assertEqual(len(PATHWAY_DEFAULTS), 2)


class TestGSLinkerLibrary(unittest.TestCase):
    """Test GS flexible linker entries."""

    def test_four_linkers(self):
        self.assertEqual(len(GS_LINKER_LIBRARY), 4)

    def test_required_fields(self):
        for name, info in GS_LINKER_LIBRARY.items():
            self.assertIn("seq", info, f"GS linker '{name}' missing 'seq'")
            self.assertIn("length", info, f"GS linker '{name}' missing 'length'")

    def test_lengths_match(self):
        for name, info in GS_LINKER_LIBRARY.items():
            self.assertEqual(
                len(info["seq"]), info["length"],
                f"GS '{name}': len(seq)={len(info['seq'])} != length={info['length']}"
            )

    def test_gs_composition(self):
        """All GS linkers should contain only G and S."""
        for name, info in GS_LINKER_LIBRARY.items():
            for aa in info["seq"]:
                self.assertIn(aa, {"G", "S"},
                              f"Non-GS residue '{aa}' in linker '{name}'")

    def test_gs_pattern(self):
        """GS linkers should be repeats of GGGGS."""
        unit = "GGGGS"
        for name, info in GS_LINKER_LIBRARY.items():
            repeats = info["length"] // 5
            self.assertEqual(info["seq"], unit * repeats,
                             f"GS linker '{name}' is not clean GGGGS repeats")

    def test_gs5_length(self):
        self.assertEqual(GS_LINKER_LIBRARY["GS5"]["length"], 5)

    def test_gs20_length(self):
        self.assertEqual(GS_LINKER_LIBRARY["GS20"]["length"], 20)


class TestAAVConstants(unittest.TestCase):
    """Test AAV packaging constraints."""

    def test_aav_max_bp(self):
        self.assertEqual(AAV_MAX_BP, 4700)

    def test_aav_max_aa(self):
        self.assertEqual(AAV_MAX_AA, 4700 // 3)

    def test_aav_max_aa_consistency(self):
        self.assertEqual(AAV_MAX_AA, AAV_MAX_BP // 3)

    def test_typical_construct_fits(self):
        """Proteasomal construct should fit in AAV (TRIM21+VVD+binder+LOV2+linker+cap)."""
        # TRIM21(236) + GS15(15) + VVD(158) + GS15(15) + binder(~100) + LOV2(136) + linker(~25) + cap(76) ≈ 761aa
        typical = 236 + 15 + 158 + 15 + 100 + 136 + 25 + 76
        self.assertLess(typical, AAV_MAX_AA)


if __name__ == "__main__":
    unittest.main()
