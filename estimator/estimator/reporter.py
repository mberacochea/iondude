import sys
from itertools import pairwise

import numpy as np

from estimator.samcounter import MAX_HOMOPOLYMER

BASES = ["A", "C", "G", "T"]
HOMOPOLYMER_LENGTHS = list(range(1, MAX_HOMOPOLYMER + 1))


def frequencies(errors, correct):
    """
    Get the channel matrix: the frequency of each true symbol (row) read as each symbol (column).

    The total of a row is the number of times the true symbol was seen: the correct calls plus
    the errors. A row without observations has no error data, so it is the identity row.
    :param errors: (n, n) counts, [true][observed], the diagonal is 0
    :param correct: (n,) counts of the correct calls
    :return: the (n, n) channel matrix and its inverse
    """
    counts = errors + np.diag(correct)
    totals = counts.sum(axis=1)
    freq = np.eye(len(correct))
    observed = totals > 0
    freq[observed] = counts[observed] / totals[observed, None]
    return freq, np.linalg.inv(freq)


class Reporter:
    def __init__(self, sam_counter):
        self.sam_counter = sam_counter

    def save(self, out=None):
        """
        Save the matrices
        :param out: output file path, None for stdout
        """
        if out is None:
            self.write(sys.stdout)
        else:
            with open(out, "w") as file:
                self.write(file)

    def write(self, file):
        """
        Write the substitutions and the homopolymers matrices
        :param file: an open file
        """
        counter = self.sam_counter

        print("##Substitutions", file=file)
        self.write_matrices(file, "S", "Substitutions", BASES,
                            counter.subs_per_base, counter.correct_bases)

        print("##Homopolymers", file=file)
        self.write_matrices(file, "I", "Indels", HOMOPOLYMER_LENGTHS,
                            counter.indels_per_pos, counter.correct_homopolymers)

    def write_matrices(self, file, code, label, names, errors_per_pos, correct_per_pos):
        """
        Write the counters, the frequencies and the inverted frequencies for all the positions,
        then for each range of positions.
        :param file: an open file
        :param code: matrix code in the header, S or I
        :param label: matrix name in the header
        :param names: row and column names
        :param errors_per_pos: (positions, n, n) error counts
        :param correct_per_pos: (positions, n) correct counts
        """
        ranges = self.sam_counter.ranges
        blocks = [("", slice(None))]
        blocks += [(f"{low}-{high} ", slice(low, high)) for low, high in pairwise(ranges)]

        for range_name, positions in blocks:
            errors = errors_per_pos[positions].sum(axis=0)
            correct = correct_per_pos[positions].sum(axis=0)
            freq, freq_inverted = frequencies(errors, correct)

            print(f"#{code} {range_name}{label} counters", file=file)
            self.write_matrix(file, names, errors)

            print(f"#{code}F {range_name}{label} frequencies", file=file)
            self.write_matrix(file, names, freq)

            print(f"#{code}FI {range_name}{label} frequencies inverted", file=file)
            self.write_matrix(file, names, freq_inverted)

    @staticmethod
    def write_matrix(file, names, matrix):
        """
        Write a matrix, tab separated, with the row and column names.
        :param file: Output file
        :param names: Row and column names
        :param matrix: 2D numpy array
        """
        print("\t".join(["-"] + [str(name) for name in names]), file=file)
        for name, row in zip(names, matrix, strict=True):
            print("\t".join([str(name)] + [str(value) for value in row]), file=file)
