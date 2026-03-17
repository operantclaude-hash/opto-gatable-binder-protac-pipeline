"""Tests for assemble_degron.py — VVD-degron fusion assembly and AAV checks."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from assemble_degron import (
    validate_degron_config,
    resolve_degron_components,
    build_degron_fusion,
    check_aav_compatibility,
)
from degron_library import (
    VVD_SEQ, VVD_LENGTH,
    E3_LIGASE_LIBRARY,
    GS_LINKER_LIBRARY,
    AAV_MAX_BP,
)


def _make_mock_opto(binder_len=20, lov2_len=10, linker_len=5, cap_len=8, spacer_len=4):
    """Create mock opto-fusion sequence and domain map."""
    total = binder_len + spacer_len + lov2_len + linker_len + cap_len
    seq = "A" * binder_len + "S" * spacer_len + "E" * lov2_len + "G" * linker_len + "P" * cap_len
    pos = 0
    domains = {}
    domains["binder"] = [pos, pos + binder_len]; pos += binder_len
    domains["spacer"] = [pos, pos + spacer_len]; pos += spacer_len
    domains["lov2"] = [pos, pos + lov2_len]; pos += lov2_len
    domains["linker"] = [pos, pos + linker_len]; pos += linker_len
    domains["cap"] = [pos, pos + cap_len]; pos += cap_len
    return seq, domains


class TestValidateDegronConfig(unittest.TestCase):

    def test_valid_defaults(self):
        cfg = {}
        result = validate_degron_config(cfg)
        self.assertEqual(result["pathway"], "proteasomal")
        self.assertEqual(result["e3_ligase"], "TRIM21_RBCC")
        self.assertEqual(result["placement"], "n_term")
        self.assertEqual(result["e3_vvd_linker"], "GS15")
        self.assertEqual(result["vvd_opto_linker"], "GS15")

    def test_valid_lysosomal(self):
        cfg = {"degron": {"pathway": "lysosomal", "e3_ligase": "NEDD4L_HECT"}}
        result = validate_degron_config(cfg)
        self.assertEqual(result["pathway"], "lysosomal")
        self.assertEqual(result["e3_ligase"], "NEDD4L_HECT")

    def test_valid_c_term(self):
        cfg = {"degron": {"placement": "c_term"}}
        result = validate_degron_config(cfg)
        self.assertEqual(result["placement"], "c_term")

    def test_invalid_pathway_raises(self):
        cfg = {"degron": {"pathway": "mitochondrial"}}
        with self.assertRaises(ValueError) as ctx:
            validate_degron_config(cfg)
        self.assertIn("pathway", str(ctx.exception))

    def test_invalid_e3_raises(self):
        cfg = {"degron": {"e3_ligase": "FAKE_LIGASE"}}
        with self.assertRaises(ValueError) as ctx:
            validate_degron_config(cfg)
        self.assertIn("e3_ligase", str(ctx.exception))

    def test_invalid_placement_raises(self):
        cfg = {"degron": {"placement": "middle"}}
        with self.assertRaises(ValueError) as ctx:
            validate_degron_config(cfg)
        self.assertIn("placement", str(ctx.exception))

    def test_invalid_linker_raises(self):
        cfg = {"degron": {"e3_vvd_linker": "GS99"}}
        with self.assertRaises(ValueError) as ctx:
            validate_degron_config(cfg)
        self.assertIn("e3_vvd_linker", str(ctx.exception))

    def test_custom_linkers(self):
        cfg = {"degron": {"e3_vvd_linker": "GS5", "vvd_opto_linker": "GS20"}}
        result = validate_degron_config(cfg)
        self.assertEqual(result["e3_vvd_linker"], "GS5")
        self.assertEqual(result["vvd_opto_linker"], "GS20")


class TestResolveDegronComponents(unittest.TestCase):

    def test_proteasomal_components(self):
        cfg = validate_degron_config({})
        dc = resolve_degron_components(cfg)
        self.assertEqual(dc["e3_name"], "TRIM21_RBCC")
        self.assertEqual(dc["e3_seq"], E3_LIGASE_LIBRARY["TRIM21_RBCC"]["seq"])
        self.assertEqual(dc["vvd_seq"], VVD_SEQ)
        self.assertEqual(dc["vvd_length"], VVD_LENGTH)
        self.assertEqual(dc["e3_pathway"], "proteasomal")
        self.assertEqual(dc["e3_ub_chain"], "K48")

    def test_lysosomal_components(self):
        cfg = validate_degron_config({"degron": {"pathway": "lysosomal"}})
        dc = resolve_degron_components(cfg)
        self.assertEqual(dc["e3_name"], "NEDD4L_HECT")
        self.assertEqual(dc["e3_pathway"], "lysosomal")
        self.assertEqual(dc["e3_ub_chain"], "K63")

    def test_gs_linker_seqs(self):
        cfg = validate_degron_config({"degron": {"e3_vvd_linker": "GS10", "vvd_opto_linker": "GS5"}})
        dc = resolve_degron_components(cfg)
        self.assertEqual(dc["gs1_seq"], GS_LINKER_LIBRARY["GS10"]["seq"])
        self.assertEqual(dc["gs2_seq"], GS_LINKER_LIBRARY["GS5"]["seq"])
        self.assertEqual(dc["gs1_length"], 10)
        self.assertEqual(dc["gs2_length"], 5)


class TestBuildDegronFusionNTerm(unittest.TestCase):
    """Test N-terminal degron placement."""

    def setUp(self):
        self.opto_seq, self.opto_domains = _make_mock_opto()
        cfg = validate_degron_config({"degron": {
            "e3_vvd_linker": "GS5", "vvd_opto_linker": "GS5"
        }})
        self.dc = resolve_degron_components(cfg)

    def test_n_term_sequence_order(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        # Should start with E3 sequence
        self.assertTrue(full.startswith(self.dc["e3_seq"]))
        # Should end with opto cap
        self.assertTrue(full.endswith(self.opto_seq[-8:]))  # cap is 8aa 'P's

    def test_n_term_total_length(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        expected = (self.dc["e3_length"] + self.dc["gs1_length"] +
                    self.dc["vvd_length"] + self.dc["gs2_length"] +
                    len(self.opto_seq))
        self.assertEqual(len(full), expected)

    def test_n_term_domains_contiguous(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        all_ranges = sorted(domains.values(), key=lambda r: r[0])
        for i in range(len(all_ranges) - 1):
            self.assertEqual(all_ranges[i][1], all_ranges[i + 1][0],
                             f"Gap between domains at position {all_ranges[i][1]}")

    def test_n_term_opto_domains_shifted(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        prefix_len = (self.dc["e3_length"] + self.dc["gs1_length"] +
                      self.dc["vvd_length"] + self.dc["gs2_length"])
        for name in self.opto_domains:
            self.assertEqual(domains[name][0], self.opto_domains[name][0] + prefix_len)
            self.assertEqual(domains[name][1], self.opto_domains[name][1] + prefix_len)

    def test_n_term_domain_extraction(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        self.assertEqual(full[domains["e3_ligase"][0]:domains["e3_ligase"][1]],
                         self.dc["e3_seq"])
        self.assertEqual(full[domains["vvd"][0]:domains["vvd"][1]], self.dc["vvd_seq"])
        self.assertEqual(full[domains["gs1"][0]:domains["gs1"][1]], self.dc["gs1_seq"])
        self.assertEqual(full[domains["gs2"][0]:domains["gs2"][1]], self.dc["gs2_seq"])

    def test_n_term_binder_extractable(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        binder = full[domains["binder"][0]:domains["binder"][1]]
        self.assertEqual(binder, "A" * 20)

    def test_n_term_e3_starts_at_zero(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        self.assertEqual(domains["e3_ligase"][0], 0)


class TestBuildDegronFusionCTerm(unittest.TestCase):
    """Test C-terminal degron placement."""

    def setUp(self):
        self.opto_seq, self.opto_domains = _make_mock_opto()
        cfg = validate_degron_config({"degron": {
            "placement": "c_term",
            "e3_vvd_linker": "GS5", "vvd_opto_linker": "GS5",
        }})
        self.dc = resolve_degron_components(cfg)

    def test_c_term_sequence_order(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        # Should start with binder (opto-fusion first)
        self.assertTrue(full.startswith("A" * 20))
        # Should end with E3
        self.assertTrue(full.endswith(self.dc["e3_seq"]))

    def test_c_term_total_length(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        expected = (len(self.opto_seq) + self.dc["gs2_length"] +
                    self.dc["vvd_length"] + self.dc["gs1_length"] +
                    self.dc["e3_length"])
        self.assertEqual(len(full), expected)

    def test_c_term_opto_domains_unshifted(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        for name in self.opto_domains:
            self.assertEqual(domains[name], self.opto_domains[name])

    def test_c_term_binder_starts_at_zero(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        self.assertEqual(domains["binder"][0], 0)

    def test_c_term_e3_at_end(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        self.assertEqual(domains["e3_ligase"][1], len(full))

    def test_c_term_domain_extraction(self):
        full, domains = build_degron_fusion(self.opto_seq, self.opto_domains, self.dc)
        self.assertEqual(full[domains["e3_ligase"][0]:domains["e3_ligase"][1]],
                         self.dc["e3_seq"])
        self.assertEqual(full[domains["vvd"][0]:domains["vvd"][1]], self.dc["vvd_seq"])


class TestBuildDegronFusionWithRealSequences(unittest.TestCase):
    """Test with actual E3/VVD sequences for size verification."""

    def test_proteasomal_nterm_size(self):
        opto_seq, opto_domains = _make_mock_opto(binder_len=100, lov2_len=136,
                                                  linker_len=23, cap_len=76, spacer_len=6)
        cfg = validate_degron_config({})
        dc = resolve_degron_components(cfg)
        full, domains = build_degron_fusion(opto_seq, opto_domains, dc)
        # TRIM21(236) + GS15(15) + VVD(158) + GS15(15) + opto(341) = 765
        self.assertEqual(len(full), 236 + 15 + 158 + 15 + 341)

    def test_lysosomal_nterm_size(self):
        opto_seq, opto_domains = _make_mock_opto(binder_len=100, lov2_len=136,
                                                  linker_len=23, cap_len=76, spacer_len=6)
        cfg = validate_degron_config({"degron": {"pathway": "lysosomal"}})
        dc = resolve_degron_components(cfg)
        full, domains = build_degron_fusion(opto_seq, opto_domains, dc)
        # NEDD4L(377) + GS15(15) + VVD(158) + GS15(15) + opto(341) = 906
        self.assertEqual(len(full), 377 + 15 + 158 + 15 + 341)


class TestCheckAAVCompatibility(unittest.TestCase):

    def test_small_construct_compatible(self):
        result = check_aav_compatibility(500)  # 1500bp
        self.assertTrue(result["compatible"])
        self.assertIsNone(result["warning"])

    def test_oversized_construct(self):
        result = check_aav_compatibility(2000)  # 6000bp
        self.assertFalse(result["compatible"])
        self.assertIn("Exceeds", result["warning"])

    def test_near_limit_warning(self):
        # 4700bp / 3 ≈ 1566, so 1500aa = 4500bp → 200bp remaining < 500
        result = check_aav_compatibility(1500)
        self.assertTrue(result["compatible"])
        self.assertIsNotNone(result["warning"])
        self.assertIn("Near", result["warning"])

    def test_typical_proteasomal_fits(self):
        # TRIM21(236) + GS15(15) + VVD(158) + GS15(15) + opto(~350) ≈ 774aa
        result = check_aav_compatibility(774)
        self.assertTrue(result["compatible"])
        self.assertIsNone(result["warning"])

    def test_bp_calculation(self):
        result = check_aav_compatibility(100)
        self.assertEqual(result["total_bp"], 300)
        self.assertEqual(result["total_aa"], 100)
        self.assertEqual(result["remaining_bp"], AAV_MAX_BP - 300)

    def test_exact_limit(self):
        # 4700/3 = 1566.67, so 1567aa = 4701bp → over
        result = check_aav_compatibility(1567)
        self.assertFalse(result["compatible"])


if __name__ == "__main__":
    unittest.main()
