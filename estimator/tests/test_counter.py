import os
import unittest

import numpy as np
from estimator import samcounter

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
MAX_READ_LEN = 500


class EstimatorTests(unittest.TestCase):

    def setUp(self):
        self.sam_counter = samcounter.SamCounter(max_read_len=MAX_READ_LEN, ranges_len=20)

    def test_parser_instantiated(self):
        """
        The SamCounter was instantiated
        """
        self.assertIsInstance(self.sam_counter, samcounter.SamCounter)

    def test_count_substitutions_one(self):
        """
        Count the substitutions based on a .sam file MD TAG.
        """
        self.sam_counter.run(os.path.join(DATADIR, "example1.sam"))

        self.assertEqual(self.sam_counter.subs_per_base[8][3][2], 1)
        self.assertEqual(self.sam_counter.subs_per_base[17][1][3], 1)
        self.assertEqual(self.sam_counter.subs_per_base[116][0][2], 1)

    def test_count_substitutions_with_dels(self):
        """
        Count the substitutions based on a .sam file MD TAG.
        Sam file with no subs but with deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "example2.sam"))

        # no substitutions, only deletions
        self.assertTrue(np.all(self.sam_counter.subs_per_base == np.zeros((MAX_READ_LEN, 4, 4), dtype=int)))

    def test_count_substitutions_with_subs_dels(self):
        """
        Count the substitutions based on a .sam file MD TAG.
        Sam file with subs. and deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "example3.sam"))

        # position 43 - T => A
        self.assertEqual(self.sam_counter.subs_per_base[43][3][0], 1)

    def test_count_substitutions_ex5(self):
        """
        Count the substitutions based on a .sam file MD TAG.
        Sam file with subs. and deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "example4.sam"))

        # on position 7 it has  A => G
        self.assertEqual(self.sam_counter.subs_per_base[7][0][2], 1)

        # on position 35 it has C => G
        self.assertEqual(self.sam_counter.subs_per_base[35][1][2], 1)

        # everything else should be empty
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:6] == np.zeros((6, 4, 4), dtype=int)))
        self.assertTrue(np.all(self.sam_counter.subs_per_base[8:34] == np.zeros((26, 4, 4), dtype=int)))
        self.assertTrue(np.all(self.sam_counter.subs_per_base[36:150] == np.zeros((114, 4, 4), dtype=int)))

        # correct changes per position
        # first position - no changes
        self.assertTrue(np.all(self.sam_counter.correct_bases[0] == np.array([1, 0, 0, 1])))
        # position 11: the base of read 1 is next to a deletion and is not counted
        self.assertTrue(np.all(self.sam_counter.correct_bases[11] == np.array([0, 0, 0, 1])))
        self.assertTrue(np.all(self.sam_counter.correct_bases[56] == np.array([2, 0, 0, 0])))

    def test_count_substitutions_with_subs_dels_two(self):
        """
        Count the substitutions based on a .sam file MD TAG.
        Sam file with subs. and deletions. A bit more complex
        """
        self.sam_counter.run(os.path.join(DATADIR, "example5.sam"))

        # position 47  - C => A is next to a deletion (gap placement), not counted

        # position 50  - A => C
        self.assertEqual(self.sam_counter.subs_per_base[50][0][1], 1)

        # position 52  - G => T
        self.assertEqual(self.sam_counter.subs_per_base[52][2][3], 1)

        # position 54  - T => G
        self.assertEqual(self.sam_counter.subs_per_base[54][3][2], 1)

        # position 375 - T => A
        self.assertEqual(self.sam_counter.subs_per_base[375][3][0], 1)

        # check some empty positions
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:46] == np.zeros((46, 4, 4), dtype=int)))
        self.assertTrue(np.all(self.sam_counter.subs_per_base[55:374] == np.zeros((319, 4, 4), dtype=int)))
        self.assertTrue(np.all(self.sam_counter.subs_per_base[376:378] == np.zeros((2, 4, 4), dtype=int)))

    def test_count_substitutions_reversed_read(self):
        """
        Count the substitutions on reversed reads.
        The counts are in sequencing orientation: the position is flipped
        (cycle = query_length - 1 - read_pos) and the bases are complemented.
        """
        self.sam_counter.run(os.path.join(DATADIR, "example6.sam"))

        # first read (reverse, 120 bases)
        # stored orientation: T => G at 8, C => T at 17, A => G at 116
        # sequencing orientation:
        # A => C at 111
        self.assertEqual(self.sam_counter.subs_per_base[111][0][1], 1)
        # G => A at 102
        self.assertEqual(self.sam_counter.subs_per_base[102][2][0], 1)
        # T => C at 3
        self.assertEqual(self.sam_counter.subs_per_base[3][3][1], 1)

        # second read (reverse, 57 bases)
        # stored orientation: C => G at 2, C => A at 14
        # sequencing orientation:
        # G => C at 54
        self.assertEqual(self.sam_counter.subs_per_base[54][2][1], 1)
        # G => T at 42
        self.assertEqual(self.sam_counter.subs_per_base[42][2][3], 1)

    def test_get_homopolymers(self):
        """
        Get every base run (single bases included) from a DNA sequence
        """
        seq = ""

        self.assertEqual(samcounter.get_homopolymers(seq), [])

        seq = "AAACG"

        self.assertEqual(samcounter.get_homopolymers(seq), [(0, 3, 3), (3, 4, 1), (4, 5, 1)])

        seq = "ACAAAGTGGGGGCTG"

        self.assertEqual(samcounter.get_homopolymers(seq),
                         [(0, 1, 1), (1, 2, 1), (2, 5, 3), (5, 6, 1), (6, 7, 1),
                          (7, 12, 5), (12, 13, 1), (13, 14, 1), (14, 15, 1)])

    def test_subs_ex10(self):
        """
        Count the substitutions on ex10.
        """
        self.sam_counter.run(os.path.join(DATADIR, "example7.sam"))

        # res["A"]["C"] += 1
        # res["A"]["G"] += 1

        self.assertEqual(self.sam_counter.subs_per_base[321][0][1], 1)
        self.assertEqual(self.sam_counter.subs_per_base[325][0][2], 1)

    def test_subs_ex11(self):
        """
        Validation mapping from Loman Ion Torrent reads.
        Both reads are soft-clipped; the aligned part is still counted.
        """
        self.sam_counter.run(os.path.join(DATADIR, "example8.sam"))

        # read 9TBRU:4:92 (reverse, 125 bases): T => G at cycle 74
        self.assertEqual(self.sam_counter.subs_per_base[74][3][2], 1)
        # read 9TBRU:4:774 (forward): T => C at cycle 116
        self.assertEqual(self.sam_counter.subs_per_base[116][3][1], 1)

        # everything else should be empty
        self.assertEqual(int(self.sam_counter.subs_per_base.sum()), 2)

    def test_subs_ex12(self):
        """
        Reverse reads are counted in sequencing orientation:
        the position is flipped and the bases are complemented.
        """
        self.sam_counter.run(os.path.join(DATADIR, "example9_reverted_reads.sam"))

        # stored orientation: T => G at 29 (both reads), A => T at 56 (read2)
        # sequencing orientation (60 bases):
        # A => C at 30
        self.assertEqual(self.sam_counter.subs_per_base[30][0][1], 2)
        # T => A at 3
        self.assertEqual(self.sam_counter.subs_per_base[3][3][0], 1)

    def test_insertions(self):
        """
        Validation of an insertions of one base
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_insertions.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        self.assertEqual(self.sam_counter.indels_per_pos[43][2][3], 1)

    def test_insertions2(self):
        """
        Validation of an insertions of one base
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_insertions2.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        # a run of length 1 extended to 3
        self.assertEqual(self.sam_counter.indels_per_pos[5][0][2], 1)

        self.assertEqual(self.sam_counter.indels_per_pos[55][2][3], 1)

        self.assertEqual(self.sam_counter.indels_per_pos[43][2][3], 1)

    def test_deletions(self):
        """
        Validation of deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_deletions.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        self.assertEqual(self.sam_counter.indels_per_pos[7][1][0], 2)
        self.assertEqual(self.sam_counter.indels_per_pos[42][2][0], 1)

    def test_deletions2(self):
        """
        Validation of an deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_deletions2.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        self.assertEqual(self.sam_counter.indels_per_pos[13][6][4], 1)

        # runs of 10 and more are counted in the 10 class: 11 G -> 9 G
        self.assertEqual(self.sam_counter.indels_per_pos[58][9][8], 1)

        self.assertTrue(np.all(self.sam_counter.indels_per_pos[0:12] == np.zeros((12, 10, 10), dtype=int)))

        zeros_mid = np.zeros((57 - 14, 10, 10), dtype=int)
        self.assertTrue(np.all(self.sam_counter.indels_per_pos[14:57] == zeros_mid))

        zeros_left = np.zeros((MAX_READ_LEN - 59, 10, 10), dtype=int)
        self.assertTrue(np.all(self.sam_counter.indels_per_pos[59:MAX_READ_LEN] == zeros_left))

    def test_indels(self):
        """
        Validation of an deletions
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_indels.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        self.assertEqual(self.sam_counter.indels_per_pos[17][2][3], 1)
        self.assertEqual(self.sam_counter.indels_per_pos[44][2][1], 1)

        self.assertEqual(self.sam_counter.indels_per_pos[50][1][3], 1)

        self.assertEqual(self.sam_counter.indels_per_pos[25][1][0], 1)

    def test_indels_reversed(self):
        """
        Validation of indels on reversed reads: the cycle is flipped
        (cycle = query_length - 1 - position of the last base of the run)
        """
        self.sam_counter.run(os.path.join(DATADIR, "indels/example_indels_reverted.sam"))

        zeros = np.zeros((MAX_READ_LEN, 4, 4), dtype=int)
        self.assertTrue(np.all(self.sam_counter.subs_per_base[0:MAX_READ_LEN] == zeros))

        # read1 reversed (60 bases)
        self.assertEqual(self.sam_counter.indels_per_pos[39][2][3], 1)
        self.assertEqual(self.sam_counter.indels_per_pos[14][2][1], 1)

        self.assertEqual(self.sam_counter.indels_per_pos[50][1][3], 1)

        # read3 reversed (59 bases)
        self.assertEqual(self.sam_counter.indels_per_pos[33][1][0], 1)

if __name__ == '__main__':
    unittest.main()
