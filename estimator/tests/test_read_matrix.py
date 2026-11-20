import os
import tempfile
import unittest

from estimator.read_matrix import read_matrix_file

MATRICES = """##Substitutions
#SF 0-20 Substitutions frequencies
-\tA\tC\tG\tT
A\t0.9\t0.1\t0.0\t0.0
C\t0.0\t1.0\t0.0\t0.0
G\t0.0\t0.0\t1.0\t0.0
T\t0.0\t0.0\t0.0\t1.0
#S 20-40 Substitutions counters
-\tA\tC\tG\tT\t
A\t0\t1\t0\t0\t
C\t0\t0\t0\t0\t
G\t0\t0\t0\t2\t
T\t0\t0\t0\t0\t
"""


class ReadMatrixTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "matrices.txt")
        with open(self.path, "w") as file:
            file.write(MATRICES)

    def tearDown(self):
        self.tmp.cleanup()

    def test_all_matrices_loaded(self):
        """
        The last matrix of the file is loaded too
        """
        matrices = read_matrix_file(self.path)

        self.assertEqual(list(matrices), ["#SF 0-20 Substitutions frequencies", "#S 20-40 Substitutions counters"])

    def test_values_and_ranges(self):
        """
        Rows with and without a trailing tab are read
        """
        matrices = read_matrix_file(self.path)

        freq = matrices["#SF 0-20 Substitutions frequencies"]
        self.assertEqual((freq.min, freq.max), (0, 20))
        self.assertTrue(freq.is_float)
        self.assertEqual(freq.data.shape, (4, 4))
        self.assertAlmostEqual(freq.data[0, 1], 0.1)

        counters = matrices["#S 20-40 Substitutions counters"]
        self.assertFalse(counters.is_float)
        self.assertEqual(counters.data[2, 3], 2)

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            read_matrix_file(os.path.join(self.tmp.name, "missing.txt"))
