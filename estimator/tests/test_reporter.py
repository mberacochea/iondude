import io
import os
import tempfile
import unittest

import numpy as np

from estimator import cli, reporter, samcounter
from estimator.read_matrix import read_matrix_file

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
MAX_READ_LEN = 500


class FrequenciesTests(unittest.TestCase):

    def test_row_total_includes_errors(self):
        """
        9 correct calls and 1 error give 0.9 and 0.1, not 1 / 9
        """
        errors = np.zeros((4, 4), dtype=int)
        errors[0][1] = 1
        correct = np.array([9, 10, 10, 10])

        freq, freq_inverted = reporter.frequencies(errors, correct)

        self.assertAlmostEqual(freq[0][0], 0.9)
        self.assertAlmostEqual(freq[0][1], 0.1)
        np.testing.assert_allclose(freq.sum(axis=1), np.ones(4))
        np.testing.assert_allclose(freq @ freq_inverted, np.eye(4), atol=1e-12)

    def test_row_without_observations_is_identity(self):
        """
        A symbol that was never seen gives an identity row, not NaN
        """
        errors = np.zeros((4, 4), dtype=int)
        correct = np.array([5, 0, 5, 5])

        freq, freq_inverted = reporter.frequencies(errors, correct)

        np.testing.assert_array_equal(freq, np.eye(4))
        self.assertTrue(np.all(np.isfinite(freq_inverted)))


class ReporterTests(unittest.TestCase):

    def setUp(self):
        self.sam_counter = samcounter.SamCounter(max_read_len=MAX_READ_LEN, ranges_len=100)
        self.sam_counter.run(os.path.join(DATADIR, "example5.sam"))

    def write(self):
        out = io.StringIO()
        reporter.Reporter(self.sam_counter).write(out)
        return out.getvalue()

    def test_all_matrices_in_the_output(self):
        """
        The matrices for all the positions and for each range are in the same output
        """
        headers = [line for line in self.write().splitlines() if line.startswith("#")]

        self.assertIn("#SF Substitutions frequencies", headers)
        self.assertIn("#SFI Substitutions frequencies inverted", headers)
        self.assertIn("#IF Indels frequencies", headers)
        self.assertIn("#SF 0-100 Substitutions frequencies", headers)
        self.assertIn("#IFI 400-500 Indels frequencies inverted", headers)
        # 2 titles + 3 matrices x (1 + 5 ranges) x 2 types
        self.assertEqual(len(headers), 2 + 3 * 6 * 2)

    def test_output_has_no_nan(self):
        """
        Ranges without reads give finite values
        """
        self.assertNotIn("nan", self.write())

    def test_read_back(self):
        """
        The matrices load with read_matrix_file
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "matrices.txt")
            reporter.Reporter(self.sam_counter).save(path)
            matrices = read_matrix_file(path)

        self.assertEqual(len(matrices), 3 * 6 * 2)
        subs = matrices["#S Substitutions counters"].data
        self.assertEqual(subs.shape, (4, 4))
        self.assertEqual(subs.sum(), self.sam_counter.subs_per_base.sum())
        homs = matrices["#IF 400-500 Indels frequencies"].data
        self.assertEqual(homs.shape, (10, 10))
        np.testing.assert_allclose(homs.sum(axis=1), np.ones((10, 1)), atol=1e-3)


class CliTests(unittest.TestCase):

    def test_main_writes_the_matrices(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "matrices.txt")
            cli.main([os.path.join(DATADIR, "example5.sam"), "-l", "500", "-o", path])
            matrices = read_matrix_file(path)

        self.assertIn("#SF 0-20 Substitutions frequencies", matrices)
        self.assertIn("#IFI 480-500 Indels frequencies inverted", matrices)
