"""Small checks of choices that affect the measurement or its denominators."""

import hashlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from ms2 import download
from ms2.protocol import DATA_PATH
from ms2.examples import validate_examples
from ms2.recompute import comparison_sets
from ms2.results import agreement
from ms2.similarity import clean_peaks, match_peaks, maximum_scores, score_pairs
from ms2.tables import peak_list_flags


class MeasurementTests(unittest.TestCase):
    def test_precursor_boundary_and_empty_spectrum(self):
        mz, values = clean_peaks(
            np.array([50.0, 99.5, 99.5001, 110.0]), np.ones(4), 0, 4, 100.0
        )
        np.testing.assert_array_equal(mz, [50.0, 99.5])
        np.testing.assert_array_equal(values, [1.0, 1.0])
        scores = score_pairs(np.array([100.0]), np.ones(1), np.array([0, 1]),
                             np.array([100.0]), np.array([0]), np.array([0]))
        self.assertTrue(np.isnan(scores[0]))

    def test_greedy_closest_peak_is_used_only_once(self):
        first, second = match_peaks(np.array([100.0, 100.006]), np.array([100.004]))
        np.testing.assert_array_equal(first, [1])
        np.testing.assert_array_equal(second, [0])

    def test_absolute_and_relative_tolerances(self):
        self.assertEqual(len(match_peaks(np.array([100.0]), np.array([100.0099]))[0]), 1)
        self.assertEqual(len(match_peaks(np.array([100.0]), np.array([100.0101]))[0]), 0)
        self.assertEqual(len(match_peaks(np.array([2000.0]), np.array([2000.019]))[0]), 1)

    def test_square_root_weights_and_unmatched_peak_denominator(self):
        scores = score_pairs(
            np.array([50.0, 70.0, 50.0]), np.array([1.0, 0.25, 1.0]),
            np.array([0, 2, 3]), np.array([100.0, 100.0]), np.array([0]), np.array([1])
        )
        self.assertAlmostEqual(scores[0], 1 / np.sqrt(1.25), places=14)

    def test_first_reference_attaining_the_maximum(self):
        scores, references = maximum_scores(
            np.array([50.0, 50.0, 50.0]), np.ones(3), np.array([0, 1, 2, 3]),
            np.full(3, 100.0), np.array([0]), [np.array([2, 1])]
        )
        self.assertEqual(scores[0], 1.0)
        self.assertEqual(references[0], 2)

    def test_condition_eligibility_does_not_limit_reference_search(self):
        meta = pd.DataFrame({
            "row": [0, 1, 2, 3, 4], "inchikey": ["A", "A", "B", "B", "C"],
            "formula": ["F"] * 5, "adduct": ["H"] * 5,
            "instrument_type": ["I", "I", "I", None, "I"],
            "collision_energy": [20.0, 20.0, 20.0, None, 40.0],
            "cond": ["H|I|20.0", "H|I|20.0", "H|I|20.0", "H|NA|-1.0", "H|I|40.0"],
        })
        references, repeats, units = comparison_sets(meta)
        np.testing.assert_array_equal(repeats, [0, 1])
        np.testing.assert_array_equal(references[("B", "H")], [2, 3])
        self.assertEqual(set(map(tuple, units[["row", "B"]].to_numpy())), {(0, "B"), (1, "B"), (2, "A")})

    def test_tolerant_and_exact_peak_lists_are_distinct(self):
        tolerant, exact = peak_list_flags(
            np.array([50.0, 70.0, 50.0, 70.0]), np.array([1.0, 0.5, 1.0, 0.500001]),
            np.array([0, 2, 4]), np.array([0]), np.array([1])
        )
        self.assertTrue(tolerant[0])
        self.assertFalse(exact[0])

    def test_ordered_repeat_agreement(self):
        frame = pd.DataFrame({"query_molecule": ["A"] * 3, "reference_molecule": ["B"] * 3,
                              "condition": ["c"] * 3, "within": [True, True, False]})
        result = agreement(frame)
        self.assertEqual(result["given_within"]["n"], 4)
        self.assertEqual(result["given_within"]["fraction"], 0.5)
        self.assertEqual(result["given_outside"]["n"], 2)
        self.assertEqual(result["given_outside"]["fraction"], 0.0)

    def test_included_example_scores(self):
        peaks = pd.read_csv(DATA_PATH / "example_peaks.csv", float_precision="round_trip")
        metadata = pd.read_csv(DATA_PATH / "example_spectra.csv", float_precision="round_trip").set_index("spectrum_id")
        units = pd.read_parquet(DATA_PATH / "isomer_units.parquet")
        pairs = pd.read_csv(DATA_PATH / "example_pairs.csv", float_precision="round_trip")
        for pair in pairs.itertuples(index=False):
            with self.subTest(pair=pair.pair_id):
                ids = [pair.query_id, pair.reference_id]
                parts = [peaks.loc[peaks.spectrum_id == key].sort_values("peak_index") for key in ids]
                masses = np.concatenate([part.mz.to_numpy() for part in parts])
                values = np.concatenate([part.relative_intensity.to_numpy() for part in parts])
                offsets = np.array([0, len(parts[0]), len(parts[0]) + len(parts[1])])
                precursor = metadata.loc[ids, "precursor_mz"].to_numpy()
                score = score_pairs(masses, values, offsets, precursor, np.array([0]), np.array([1]))[0]
                stored = units.loc[(units.query_id == ids[0]) & (units.reference_id == ids[1]), "score"].item()
                self.assertAlmostEqual(score, stored, places=13)
                self.assertEqual(stored, pair.score)

    def test_example_peak_corruption_is_rejected(self):
        peaks = pd.read_csv(DATA_PATH / "example_peaks.csv", float_precision="round_trip")
        metadata = pd.read_csv(DATA_PATH / "example_spectra.csv", float_precision="round_trip")
        pairs = pd.read_csv(DATA_PATH / "example_pairs.csv", float_precision="round_trip")
        units = pd.read_parquet(DATA_PATH / "isomer_units.parquet")
        spectra = pd.read_parquet(DATA_PATH / "spectra.parquet")
        peaks.loc[0, "retained"] = not peaks.loc[0, "retained"]
        with self.assertRaises(ValueError):
            validate_examples(pairs, metadata, peaks, units, spectra)


class DownloadTests(unittest.TestCase):
    def test_verified_download_and_existing_file(self):
        content = b"small public input\n"
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "input.tsv"
            with patch.object(download, "INPUT_BYTES", len(content)), patch.object(
                download, "SHA256", hashlib.sha256(content).hexdigest()
            ), patch.object(download, "urlopen", return_value=io.BytesIO(content)) as fetch:
                download.download(destination)
                download.download(destination)
                self.assertEqual(destination.read_bytes(), content)
                fetch.assert_called_once()

    def test_wrong_hash_leaves_no_published_file(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "input.tsv"
            with patch.object(download, "INPUT_BYTES", 3), patch.object(download, "SHA256", "wrong"), patch.object(
                download, "urlopen", return_value=io.BytesIO(b"abc")
            ):
                with self.assertRaises(ValueError):
                    download.download(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(Path(directory).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
