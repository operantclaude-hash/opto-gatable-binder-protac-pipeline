"""Tests for analyze_binding_validation.py — confidence extraction from Boltz-2 results."""

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

# We need to patch module-level paths before importing, so import the function
# via importlib after patching.
import importlib


def _make_mock_confidence(result_dir, iptm=0.85, ptm=0.9, plddt=0.75,
                          ipde=0.5, iplddt=0.7, confidence_score=0.8,
                          model_idx=0, use_protein_iptm=True):
    """Create a mock Boltz-2 confidence JSON file."""
    predictions_dir = result_dir / "predictions" / "test"
    predictions_dir.mkdir(parents=True, exist_ok=True)
    conf = {
        "ptm": ptm,
        "complex_plddt": plddt,
        "complex_ipde": ipde,
        "complex_iplddt": iplddt,
        "confidence_score": confidence_score,
    }
    if use_protein_iptm:
        conf["protein_iptm"] = iptm
    else:
        conf["iptm"] = iptm
    path = predictions_dir / f"confidence_test_model_{model_idx}.json"
    path.write_text(json.dumps(conf))
    return path


def _make_mock_pae(result_dir, chain_a_len, chain_b_len, model_idx=0,
                   cross_chain_val=5.0, same_chain_val=2.0):
    """Create a mock PAE .npz file with predictable cross-chain values."""
    total = chain_a_len + chain_b_len
    pae = np.full((total, total), same_chain_val, dtype=np.float32)
    # Cross-chain blocks
    pae[:chain_a_len, chain_a_len:] = cross_chain_val
    pae[chain_a_len:, :chain_a_len] = cross_chain_val
    predictions_dir = result_dir / "predictions" / "test"
    predictions_dir.mkdir(parents=True, exist_ok=True)
    path = predictions_dir / f"pae_test_model_{model_idx}.npz"
    np.savez(path, pae=pae)
    return path


def _make_manifest(manifest_path, entries):
    """Create a mock manifest.json."""
    manifest_path.write_text(json.dumps(entries, indent=2))


class TestExtractConfidence(unittest.TestCase):
    """Test the extract_confidence function with mock Boltz-2 outputs."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        self.results_dir = self.tmp / "binding_validation" / "boltz2_results"
        self.results_dir.mkdir(parents=True)
        self.manifest_path = self.tmp / "binding_validation" / "manifest.json"

    def tearDown(self):
        self.tmpdir.cleanup()

    def _get_extract_fn(self):
        """Import extract_confidence with patched paths."""
        # Read the source and patch constants
        src_file = Path(__file__).parent.parent / "campaigns" / "constitutive_top10" / "analyze_binding_validation.py"
        spec = importlib.util.spec_from_file_location("abv", src_file)
        mod = importlib.util.module_from_spec(spec)
        # Patch paths before exec
        mod.MANIFEST_PATH = self.manifest_path
        mod.RESULTS_DIR = self.results_dir
        spec.loader.exec_module(mod)
        # Re-patch after exec (module-level code sets them)
        mod.MANIFEST_PATH = self.manifest_path
        mod.RESULTS_DIR = self.results_dir
        return mod.extract_confidence

    def test_basic_extraction(self):
        construct = "test_binder_HAT"
        rdir = self.results_dir / construct
        _make_mock_confidence(rdir, iptm=0.85)
        chain_a, chain_b = 20, 30
        _make_mock_pae(rdir, chain_a, chain_b, cross_chain_val=4.5)
        _make_manifest(self.manifest_path, [
            {"construct_name": construct, "chain_a_length": chain_a}
        ])
        extract = self._get_extract_fn()
        metrics = extract(rdir, construct)

        self.assertIsNotNone(metrics)
        self.assertAlmostEqual(metrics["iptm"], 0.85)
        self.assertAlmostEqual(metrics["ptm"], 0.9)
        self.assertAlmostEqual(metrics["min_ipae"], 4.5, places=1)
        self.assertAlmostEqual(metrics["mean_ipae"], 4.5, places=1)

    def test_missing_confidence_returns_none(self):
        rdir = self.results_dir / "missing"
        rdir.mkdir(parents=True)
        extract = self._get_extract_fn()
        result = extract(rdir, "missing")
        self.assertIsNone(result)

    def test_iptm_fallback_key(self):
        """When protein_iptm is absent, should use iptm key."""
        construct = "fallback_test"
        rdir = self.results_dir / construct
        _make_mock_confidence(rdir, iptm=0.77, use_protein_iptm=False)
        _make_manifest(self.manifest_path, [
            {"construct_name": construct, "chain_a_length": 10}
        ])
        extract = self._get_extract_fn()
        metrics = extract(rdir, construct)
        self.assertAlmostEqual(metrics["iptm"], 0.77)

    def test_no_pae_file(self):
        """When PAE file is missing, min_ipae and mean_ipae should be None."""
        construct = "no_pae"
        rdir = self.results_dir / construct
        _make_mock_confidence(rdir, iptm=0.9)
        extract = self._get_extract_fn()
        metrics = extract(rdir, construct)
        self.assertIsNotNone(metrics)
        self.assertAlmostEqual(metrics["iptm"], 0.9)
        self.assertIsNone(metrics["min_ipae"])
        self.assertIsNone(metrics["mean_ipae"])

    def test_cross_chain_pae_block(self):
        """Verify cross-chain PAE extraction with asymmetric values."""
        construct = "asym_pae"
        rdir = self.results_dir / construct
        chain_a, chain_b = 5, 8
        total = chain_a + chain_b

        # Build PAE with distinct cross-chain blocks
        predictions_dir = rdir / "predictions" / "test"
        predictions_dir.mkdir(parents=True)
        pae = np.full((total, total), 10.0, dtype=np.float32)
        pae[:chain_a, chain_a:] = 3.0  # A→B
        pae[chain_a:, :chain_a] = 7.0  # B→A
        np.savez(predictions_dir / "pae_test_model_0.npz", pae=pae)

        _make_mock_confidence(rdir, iptm=0.8)
        _make_manifest(self.manifest_path, [
            {"construct_name": construct, "chain_a_length": chain_a}
        ])
        extract = self._get_extract_fn()
        metrics = extract(rdir, construct)

        # min should be 3.0 (from A→B block)
        self.assertAlmostEqual(metrics["min_ipae"], 3.0, places=1)
        # mean should be (3.0*40 + 7.0*40) / 80 = 5.0
        self.assertAlmostEqual(metrics["mean_ipae"], 5.0, places=1)

    def test_no_manifest_entry(self):
        """When construct not in manifest, iPAE should be None."""
        construct = "orphan"
        rdir = self.results_dir / construct
        _make_mock_confidence(rdir, iptm=0.6)
        chain_a, chain_b = 10, 10
        _make_mock_pae(rdir, chain_a, chain_b, cross_chain_val=5.0)
        _make_manifest(self.manifest_path, [
            {"construct_name": "other_construct", "chain_a_length": 10}
        ])
        extract = self._get_extract_fn()
        metrics = extract(rdir, construct)
        self.assertAlmostEqual(metrics["iptm"], 0.6)
        self.assertIsNone(metrics["min_ipae"])


if __name__ == "__main__":
    unittest.main()
