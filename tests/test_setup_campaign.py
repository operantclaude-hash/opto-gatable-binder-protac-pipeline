"""Tests for setup_campaign.py — config validation, fusion assembly, YAML generation."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from setup_campaign import (
    validate_config,
    get_linker_subset,
    build_archC_fusion,
    make_yaml_1chain,
    make_yaml_2chain,
    generate_run_script,
)
from linker_library import LOV2_DARK, LOV2_LIT, UBQ_SEQ, LINKER_LIBRARY


def _make_valid_config():
    """Return a minimal valid campaign config dict."""
    return {
        "campaign_name": "test_campaign",
        "target": {
            "name": "test_target",
            "sequence": "ACDEFGHIKLMNPQRSTVWY",
            "hotspot_residues_1idx": [1, 5, 10],
        },
        "binders": [
            {
                "name": "binder_1",
                "sequence": "AAAAAEEEEE",
                "length": 10,
                "iptm": 0.9,
            }
        ],
        "scaffold": {
            "lov2_dark_seq": LOV2_DARK,
            "lov2_lit_seq": LOV2_LIT,
            "cap_seq": UBQ_SEQ,
        },
    }


class TestValidateConfig(unittest.TestCase):

    def test_valid_config_passes(self):
        cfg = _make_valid_config()
        validate_config(cfg)  # should not raise

    def test_missing_campaign_name_raises(self):
        cfg = _make_valid_config()
        del cfg["campaign_name"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("campaign_name", str(ctx.exception))

    def test_missing_target_sequence_raises(self):
        cfg = _make_valid_config()
        cfg["target"]["sequence"] = ""
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("target.sequence", str(ctx.exception))

    def test_missing_hotspots_raises(self):
        cfg = _make_valid_config()
        del cfg["target"]["hotspot_residues_1idx"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("hotspot", str(ctx.exception))

    def test_no_binders_raises(self):
        cfg = _make_valid_config()
        cfg["binders"] = []
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("binder", str(ctx.exception).lower())

    def test_binder_missing_sequence_raises(self):
        cfg = _make_valid_config()
        del cfg["binders"][0]["sequence"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("sequence", str(ctx.exception))

    def test_binder_missing_name_raises(self):
        cfg = _make_valid_config()
        del cfg["binders"][0]["name"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("name", str(ctx.exception))

    def test_missing_lov2_dark_raises(self):
        cfg = _make_valid_config()
        del cfg["scaffold"]["lov2_dark_seq"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("lov2_dark", str(ctx.exception))

    def test_missing_lov2_lit_raises(self):
        cfg = _make_valid_config()
        del cfg["scaffold"]["lov2_lit_seq"]
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        self.assertIn("lov2_lit", str(ctx.exception))

    def test_multiple_errors_reported(self):
        cfg = {"target": {}, "binders": [], "scaffold": {}}
        with self.assertRaises(ValueError) as ctx:
            validate_config(cfg)
        msg = str(ctx.exception)
        # Should report multiple errors
        self.assertGreater(msg.count("\n"), 1)


class TestGetLinkerSubset(unittest.TestCase):

    def test_all_linkers_returns_full_library(self):
        cfg = {"scaffold": {"linkers": "all"}}
        result = get_linker_subset(cfg)
        self.assertEqual(len(result), len(LINKER_LIBRARY))

    def test_none_linkers_returns_full_library(self):
        cfg = {"scaffold": {"linkers": None}}
        result = get_linker_subset(cfg)
        self.assertEqual(len(result), len(LINKER_LIBRARY))

    def test_specific_linkers_filtered(self):
        cfg = {"scaffold": {"linkers": ["E2pppE2", "EAAAK5"]}}
        result = get_linker_subset(cfg)
        self.assertEqual(len(result), 2)
        self.assertIn("E2pppE2", result)
        self.assertIn("EAAAK5", result)

    def test_all_unknown_linkers_raises(self):
        cfg = {"scaffold": {"linkers": ["NONEXISTENT"]}}
        with self.assertRaises(ValueError):
            get_linker_subset(cfg)

    def test_returns_copy_not_reference(self):
        cfg = {"scaffold": {"linkers": "all"}}
        result = get_linker_subset(cfg)
        result["FAKE"] = {"seq": "AAA"}
        self.assertNotIn("FAKE", LINKER_LIBRARY)


class TestBuildArchCFusion(unittest.TestCase):

    def test_basic_assembly_no_spacer(self):
        binder = "AAAA"
        lov2 = "EEEE"
        linker = "GG"
        cap = "PPPP"
        fusion, domains = build_archC_fusion(binder, lov2, linker, cap)
        self.assertEqual(fusion, "AAAAEEEEGGPPPP")

    def test_assembly_with_spacer(self):
        binder = "AAAA"
        lov2 = "EEEE"
        linker = "GG"
        cap = "PPPP"
        spacer = "SS"
        fusion, domains = build_archC_fusion(binder, lov2, linker, cap, spacer)
        self.assertEqual(fusion, "AAAASSEEEEGGPPPP")

    def test_domain_boundaries_sum_to_total(self):
        binder = "AAAA"
        lov2 = "EEEEE"
        linker = "GGG"
        cap = "PPPP"
        spacer = "SS"
        fusion, domains = build_archC_fusion(binder, lov2, linker, cap, spacer)
        self.assertEqual(domains["cap"][1], len(fusion))

    def test_domain_boundaries_contiguous(self):
        fusion, domains = build_archC_fusion("AAAA", "EEEE", "GG", "PP", "SS")
        keys = ["binder", "spacer", "lov2", "linker", "cap"]
        for i in range(len(keys) - 1):
            self.assertEqual(
                domains[keys[i]][1], domains[keys[i + 1]][0],
                f"Gap between {keys[i]} and {keys[i+1]}"
            )

    def test_empty_spacer_zero_range(self):
        fusion, domains = build_archC_fusion("AAAA", "EEEE", "GG", "PP")
        self.assertEqual(domains["spacer"][0], domains["spacer"][1])

    def test_nonempty_spacer_range(self):
        fusion, domains = build_archC_fusion("AAAA", "EEEE", "GG", "PP", "GGSGGS")
        self.assertEqual(domains["spacer"][1] - domains["spacer"][0], 6)

    def test_binder_starts_at_zero(self):
        fusion, domains = build_archC_fusion("AAAA", "EEEE", "GG", "PP")
        self.assertEqual(domains["binder"][0], 0)

    def test_real_sequences(self):
        binder = "AAQSGFSPASEAAAAELLAA"  # 20aa
        linker_seq = LINKER_LIBRARY["E2pppE2"]["seq"]  # 23aa
        spacer = "GGSGGS"  # 6aa
        fusion, domains = build_archC_fusion(binder, LOV2_DARK, linker_seq, UBQ_SEQ, spacer)
        expected = 20 + 6 + 136 + 23 + 76
        self.assertEqual(len(fusion), expected)

    def test_domain_extraction(self):
        binder = "AAAA"
        lov2 = "EEEE"
        linker = "GG"
        cap = "PPPP"
        spacer = "SS"
        fusion, domains = build_archC_fusion(binder, lov2, linker, cap, spacer)
        self.assertEqual(fusion[domains["binder"][0]:domains["binder"][1]], binder)
        self.assertEqual(fusion[domains["spacer"][0]:domains["spacer"][1]], spacer)
        self.assertEqual(fusion[domains["lov2"][0]:domains["lov2"][1]], lov2)
        self.assertEqual(fusion[domains["linker"][0]:domains["linker"][1]], linker)
        self.assertEqual(fusion[domains["cap"][0]:domains["cap"][1]], cap)


class TestMakeYaml1Chain(unittest.TestCase):

    def test_valid_structure(self):
        yaml = make_yaml_1chain("ACDEF", "test")
        self.assertIn("version: 1", yaml)
        self.assertIn("sequences:", yaml)
        self.assertIn("id: A", yaml)
        self.assertIn("msa: empty", yaml)

    def test_sequence_in_yaml(self):
        yaml = make_yaml_1chain("MKTLQEAAAK", "test")
        self.assertIn("MKTLQEAAAK", yaml)

    def test_ends_with_newline(self):
        yaml = make_yaml_1chain("ACDEF", "test")
        self.assertTrue(yaml.endswith("\n"))


class TestMakeYaml2Chain(unittest.TestCase):

    def test_two_chains_present(self):
        yaml = make_yaml_2chain("AAAA", "EEEE")
        self.assertIn("id: A", yaml)
        self.assertIn("id: B", yaml)

    def test_no_constraints_when_none(self):
        yaml = make_yaml_2chain("AAAA", "EEEE")
        self.assertNotIn("constraints:", yaml)

    def test_constraints_formatted(self):
        contacts = [{"fusion_res_1idx": 5, "target_res_1idx": 10}]
        yaml = make_yaml_2chain("AAAA", "EEEE", contacts=contacts)
        self.assertIn("constraints:", yaml)
        self.assertIn("contact:", yaml)
        self.assertIn("token1:", yaml)
        self.assertIn("token2:", yaml)

    def test_contact_values_in_output(self):
        contacts = [{"fusion_res_1idx": 5, "target_res_1idx": 10, "max_distance": 8.0}]
        yaml = make_yaml_2chain("AAAA", "EEEE", contacts=contacts)
        self.assertIn("[A, 5]", yaml)
        self.assertIn("[B, 10]", yaml)
        self.assertIn("8.0", yaml)

    def test_default_max_distance(self):
        contacts = [{"fusion_res_1idx": 1, "target_res_1idx": 2}]
        yaml = make_yaml_2chain("AAAA", "EEEE", contacts=contacts)
        self.assertIn("10.0", yaml)


class TestGenerateRunScript(unittest.TestCase):

    def test_shebang_line(self):
        script = generate_run_script(Path("/tmp/y"), Path("/tmp/r"), ["test1"], 5, 3)
        self.assertTrue(script.startswith("#!/bin/bash"))

    def test_conda_activation(self):
        script = generate_run_script(Path("/tmp/y"), Path("/tmp/r"), ["test1"], 5, 3)
        self.assertIn("conda activate boltz_only", script)

    def test_all_yaml_names_present(self):
        names = ["construct_a", "construct_b", "construct_c"]
        script = generate_run_script(Path("/tmp/y"), Path("/tmp/r"), names, 5, 3)
        for name in names:
            self.assertIn(name, script)

    def test_total_count_correct(self):
        names = ["a", "b", "c"]
        script = generate_run_script(Path("/tmp/y"), Path("/tmp/r"), names, 5, 3)
        self.assertIn("TOTAL=3", script)

    def test_diffusion_samples_in_command(self):
        script = generate_run_script(Path("/tmp/y"), Path("/tmp/r"), ["test1"], 7, 5)
        self.assertIn("--diffusion_samples 7", script)
        self.assertIn("--recycling_steps 5", script)


if __name__ == "__main__":
    unittest.main()
