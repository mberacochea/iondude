import unittest

from estimator import seq


class SeqHelperTests(unittest.TestCase):

    def test_reverse(self):
        """
        Reverse complement a DNA sequence
        """
        sequence = 'ACGTTGC'
        seq_rev = 'GCAACGT'

        rev = seq.reverse_complement(sequence)

        self.assertEqual(rev, seq_rev)

    def test_reverse_withN(self):
        """
        Reverse complement a DNA sequence, with N bases
        """
        sequence = 'ACGNTTGC'

        with self.assertRaises(ValueError):
            seq.reverse_complement(sequence)
