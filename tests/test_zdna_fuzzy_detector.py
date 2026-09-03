from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "zdna_fuzzy_detector.py"
SPEC = importlib.util.spec_from_file_location("zdna_fuzzy_detector_test_module", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Could not load zdna_fuzzy_detector.py")
zfd = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = zfd
SPEC.loader.exec_module(zfd)


class CoordinateTests(unittest.TestCase):
    def _annotation_file(self, text: str) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
        self.addCleanup(Path(handle.name).unlink, missing_ok=True)
        with handle:
            handle.write(text)
        return Path(handle.name)

    def test_sga_is_detected_as_data_and_auto_converted_from_one_based(self) -> None:
        path = self._annotation_file("chr1\tTSS\t850984\t+\t1\tSAMD11\n")
        self.assertEqual(zfd.load_tss(path, "auto"), {"1": [850983]})

    def test_bed_negative_strand_uses_last_covered_base(self) -> None:
        path = self._annotation_file("chr2\t100\t120\tgene\t0\t-\n")
        self.assertEqual(zfd.load_tss(path, "auto"), {"2": [119]})

    def test_api_inclusive_end_is_normalized_to_half_open(self) -> None:
        row = {"start": "10", "end": "15", "score": "72.5"}
        interval = zfd.result_interval(row, coordinate_base=0, offset=100, chrom="chr1", source="scan")
        self.assertIsNotNone(interval)
        assert interval is not None
        self.assertEqual((interval.start, interval.end, interval.length), (110, 116, 6))


class LocalHunterTests(unittest.TestCase):
    def test_known_alternating_gc_run_has_half_open_coordinates(self) -> None:
        config = zfd.hunter_config("test", 6, 30.0)
        hits = list(zfd.local_hunter_hits("AAAAAGCGCGCNAAAAA", config))
        self.assertEqual(len(hits), 1)
        self.assertEqual((hits[0].start, hits[0].end), (5, 11))
        self.assertAlmostEqual(hits[0].score_percent, 100.0)

    def test_short_run_is_rejected(self) -> None:
        config = zfd.hunter_config("test", 6, 30.0)
        self.assertEqual(list(zfd.local_hunter_hits("AACGCGAA", config)), [])


class ModelSpecificationTests(unittest.TestCase):
    def test_complete_model_specification_is_exposed(self) -> None:
        spec = zfd.model_specification()
        self.assertEqual(spec["name"], "ZDNA-Fuzzy Hunter")
        self.assertEqual(len(spec["rules"]), 12)
        self.assertAlmostEqual(spec["final_score"]["rule_fraction"], 0.149236)
        self.assertAlmostEqual(spec["final_score"]["weighted_component_fraction"], 0.850764)
        self.assertEqual(set(spec["membership_functions"]["main_components"]), set(zfd.TERMS))
        self.assertEqual(spec["rule_evaluation"]["modes"], ["balanced", "moderate", "strict"])


class ShinReproductionTests(unittest.TestCase):
    def test_released_shin_table_reproduces_all_operating_points(self) -> None:
        repo_root = MODULE_PATH.parent
        script = repo_root / "validation" / "reproduce_shin_benchmark.py"
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                [sys.executable, str(script), "--output-dir", tmpdir],
                cwd=repo_root,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            summary = json.loads((Path(tmpdir) / "shin_benchmark_metrics.json").read_text(encoding="utf-8"))
            self.assertTrue(summary["matches_released_expected_values"])
            self.assertTrue(summary["rebuild_performed"])
            self.assertEqual(summary["source_loci"]["total"], 391)
            self.assertEqual(len(summary["excluded_loci"]), 6)
            self.assertEqual(summary["modes"]["balanced"]["tp"], 269)
            self.assertEqual(summary["modes"]["balanced"]["fp"], 7)

    def test_released_shin_supplementary_table_is_generated_output(self) -> None:
        repo_root = MODULE_PATH.parent
        script = repo_root / "validation" / "build_shin_supplementary_table.py"
        released = repo_root / "supplementary" / "Supplementary_Table_S1_Shin_per_locus.csv"
        with tempfile.TemporaryDirectory() as tmpdir:
            rebuilt = Path(tmpdir) / released.name
            completed = subprocess.run(
                [sys.executable, str(script), "--output", str(rebuilt)],
                cwd=repo_root,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertEqual(rebuilt.read_bytes(), released.read_bytes())


if __name__ == "__main__":
    unittest.main()
