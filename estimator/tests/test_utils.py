import os
import unittest

from estimator.utils import kmer_for_base, get_context, get_kmers
from estimator.readfq import Seq

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")


class UtilsTests(unittest.TestCase):

    def test_get_count(self):
        self.kmers = {
            'ACGTACGTT': 4,
            'ACGTCCGTT': 8,
            'ACGTGCGTT': 12,
            'ACGTTCGTT': 16,
        }
        self.k = 4
        self.assertEqual(4, kmer_for_base(self.kmers, 'ACGTACGTT', 'A', self.k))
        self.assertEqual(8, kmer_for_base(self.kmers, 'ACGTCCGTT', 'C', self.k))
        self.assertEqual(12, kmer_for_base(self.kmers, 'ACGTGCGTT', 'G', self.k))
        self.assertEqual(16, kmer_for_base(self.kmers, 'ACGTTCGTT', 'T', self.k))

    def test_get_context(self):
        seq = 'ACGTGGTACGTAGATA'
        self.assertEqual('CGTGGTA', get_context(seq, pos=4, k=3))
        self.assertEqual('CGTGGTACG', get_context(seq, pos=5, k=4))
        self.assertEqual('TGGTACG', get_context(seq, pos=6, k=3))

    def test_get_context_short_seq(self):
        seq = 'ACGTGGTACGTAGATAGTAGC'
        # correct
        self.assertEqual('CGTGGTA', get_context(seq, pos=4, k=3))
        self.assertEqual('CGTGGTACG', get_context(seq, pos=5, k=4))
        self.assertEqual('ACGTGGTACGTAGATAGTAGC', get_context(seq, pos=10, k=10))
        # incorrect
        self.assertEqual(None, get_context(seq, pos=10, k=11))
        self.assertEqual(None, get_context(seq, pos=11, k=20))

    def test_get_kmers(self):
        file = open(os.path.join(DATADIR, 'kmers_fasta.fasta'), 'r')
        seq = Seq.readfq(file)

        kmers = [
            'ACG',
             'CGT',
              'GTG',
               'TGG',
                'GGT',
                 'GTA',
                  'TAC',
                   'ACG',
                    'CGT',
                     'GTA',
                      'TAG',
                       'AGA',
                        'GAT',
                         'ATA'
        ]
        k = 3

        kmers_res = []

        for read in seq:
            for ret in get_kmers(read.seq, k):
                kmers_res.append(ret)

        file.close()

        self.assertSequenceEqual(kmers, kmers_res)

    def test_get_kmers_incorrect(self):
        seq = 'ACGTGGT'

        k = 8

        kmers = []
        kmers_res = []

        for ret in get_kmers(seq, k):
            if ret:
                kmers_res.append(ret)

        self.assertSequenceEqual(kmers, kmers_res)

    def test_border_detection(self):
        seq = 'ACGTGGTACGTAGATA'
        seq_len = len(seq)
        k = 3

        kmers = []

        for _i in range(seq_len):
            if _i < k or _i > seq_len - k - 1:
                continue
            kmers.append(get_context(seq, _i, k))

        kmers_expected = [
            'ACGTGGT',
             'CGTGGTA',
              'GTGGTAC',
               'TGGTACG',
                'GGTACGT',
                 'GTACGTA',
                  'TACGTAG',
                   'ACGTAGA',
                    'CGTAGAT',
                     'GTAGATA'
        ]

        self.assertSequenceEqual(kmers, kmers_expected)
