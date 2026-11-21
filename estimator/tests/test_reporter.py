import io
import json
import os
import tempfile
import unittest

import numpy as np

from estimator import cli, reporter, samcounter
from estimator.matrices import get_slot, read_matrices

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
MAX_READ_LEN = 500
SOURCE = {"tool": "iondude-estimator", "version": "0.4", "date": "2026-11-21", "input": "example5.sam"}
ESTIMATION = {"max_read_length": MAX_READ_LEN, "range": 100, "min_mapping_quality": 20}


class FrequenciesTests(unittest.TestCase):

    def test_row_total_includes_errors(self):
        """
        9 correct calls and 1 error give 0.9 and 0.1, not 1 / 9
        """
        counts = np.diag([9, 10, 10, 10])
        counts[0][1] = 1

        freq, freq_inverted = reporter.frequencies(counts)

        self.assertAlmostEqual(freq[0][0], 0.9)
        self.assertAlmostEqual(freq[0][1], 0.1)
        np.testing.assert_allclose(freq.sum(axis=1), np.ones(4))
        np.testing.assert_allclose(freq @ freq_inverted, np.eye(4), atol=1e-12)

    def test_row_without_observations_is_identity(self):
        """
        A symbol that was never seen gives an identity row, not NaN
        """
        counts = np.diag([5, 0, 5, 5])

        freq, freq_inverted = reporter.frequencies(counts)

        np.testing.assert_array_equal(freq, np.eye(4))
        self.assertTrue(np.all(np.isfinite(freq_inverted)))


class ReporterTests(unittest.TestCase):

    def setUp(self):
        self.sam_counter = samcounter.SamCounter(max_read_len=MAX_READ_LEN, ranges_len=100)
        self.sam_counter.run(os.path.join(DATADIR, "example5.sam"))
        self.reporter = reporter.Reporter(self.sam_counter, SOURCE, ESTIMATION)

    def write(self):
        out = io.StringIO()
        self.reporter.write(out)
        return out.getvalue()

    def read_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "matrices.json")
            self.reporter.save(path)
            return read_matrices(path)

    def test_all_slots_in_the_output(self):
        """
        One slot for each range, no slot for all the positions
        """
        data = json.loads(self.write())

        self.assertEqual(data["slot_type"], "position")
        self.assertEqual(data["range"], 100)
        self.assertEqual(data["max_length"], MAX_READ_LEN)
        self.assertEqual(data["source"], SOURCE)
        self.assertEqual(data["estimation"], ESTIMATION)
        for kind in ("substitutions", "homopolymers"):
            ranges = [(slot["min"], slot["max"]) for slot in data[kind]["slots"]]
            self.assertEqual(ranges, [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500)])

    def test_output_has_no_nan(self):
        """
        Ranges without reads give finite values
        """
        self.assertNotIn("NaN", self.write())

    def test_read_back(self):
        """
        The file loads with read_matrices and keeps the counts
        """
        data = self.read_back()

        subs = sum(slot["counts"] for slot in data["substitutions"]["slots"])
        self.assertEqual(subs.shape, (4, 4))
        # the errors are off the diagonal, the correct calls on the diagonal
        self.assertEqual(subs.sum() - np.trace(subs), self.sam_counter.subs_per_base.sum())
        np.testing.assert_array_equal(np.diag(subs), self.sam_counter.correct_bases.sum(axis=0))
        homs = get_slot(data, "homopolymers", 450)["channel"]
        self.assertEqual(homs.shape, (10, 10))
        np.testing.assert_allclose(homs.sum(axis=1), np.ones(10), atol=1e-12)


class CliTests(unittest.TestCase):

    def test_main_writes_the_matrices(self):
        sam = os.path.join(DATADIR, "example5.sam")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "matrices.json")
            cli.main([sam, "-l", "500", "-o", path])
            data = read_matrices(path)

        self.assertEqual(data["estimation"], {"max_read_length": 500, "range": 20, "min_mapping_quality": 20})
        self.assertEqual(data["source"]["tool"], "iondude-estimator")
        self.assertEqual(data["source"]["input"], sam)
        self.assertEqual(len(data["substitutions"]["slots"]), 25)
        self.assertEqual(get_slot(data, "homopolymers", 499)["min"], 480)
