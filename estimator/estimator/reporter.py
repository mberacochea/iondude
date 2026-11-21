import sys
from itertools import pairwise

import numpy as np

from estimator.matrices import FORMAT, FORMAT_VERSION, SYMBOLS, write_matrices


def frequencies(counts):
    """
    Get the channel matrix: the frequency of each true symbol (row) read as each symbol (column).

    The total of a row is the number of times the true symbol was seen: the correct calls plus
    the errors. A row without observations has no error data, so it is the identity row.
    :param counts: (n, n) counts, [true][observed], the correct calls on the diagonal
    :return: the (n, n) channel matrix and its inverse
    """
    totals = counts.sum(axis=1)
    freq = np.eye(len(counts))
    observed = totals > 0
    freq[observed] = counts[observed] / totals[observed, None]
    return freq, np.linalg.inv(freq)


class Reporter:
    def __init__(self, sam_counter, source, estimation):
        """
        :param sam_counter: a SamCounter that has counted the reads
        :param source: the source fields of the file: tool, version, date and input
        :param estimation: the settings of the estimation, by name
        """
        self.sam_counter = sam_counter
        self.source = source
        self.estimation = estimation

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
        Write the substitutions and the homopolymers matrices as JSON
        :param file: an open file
        """
        write_matrices(file, self.matrices())

    def matrices(self):
        """
        Get the matrices in the in-memory form of estimator.matrices
        :return: a dict with the same shape as the JSON file
        """
        counter = self.sam_counter
        return {
            "format": FORMAT,
            "format_version": FORMAT_VERSION,
            "source": self.source,
            "estimation": self.estimation,
            "slot_type": "position",
            "range": counter.ranges_len,
            "max_length": counter.max_read_len,
            "substitutions": self.section("substitutions", counter.subs_per_base, counter.correct_bases),
            "homopolymers": self.section("homopolymers", counter.indels_per_pos, counter.correct_homopolymers),
        }

    def section(self, kind, errors_per_pos, correct_per_pos):
        """
        Get the counts, the channel matrix and its inverse for each range of positions.
        :param kind: substitutions or homopolymers
        :param errors_per_pos: (positions, n, n) error counts
        :param correct_per_pos: (positions, n) correct counts
        :return: the section, a dict with the symbols and the slots
        """
        slots = []
        for low, high in pairwise(self.sam_counter.ranges):
            counts = errors_per_pos[low:high].sum(axis=0) + np.diag(correct_per_pos[low:high].sum(axis=0))
            channel, inverse = frequencies(counts)
            slots.append({"min": low, "max": high, "counts": counts, "channel": channel, "inverse": inverse})
        return {"symbols": SYMBOLS[kind], "slots": slots}
