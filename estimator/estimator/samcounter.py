import re
import sys

import numpy as np
import pysam


MINIMUM_MAPPING_QUALITY = 20

# iondude models run lengths 1..10, 10 is the "10 or more" class
MAX_HOMOPOLYMER = 10

_BASES = "ACGT"


def get_base_index(base):
    """
    Get the base index, ie. A index is 0. -1 for any other base (N, IUPAC codes, ...).
    :param base: Nucleotide
    :return: index
    """
    return _BASES.find(base.upper())


def get_homopolymers(seq):
    """
    Get the positions and lengths of every base run, single bases included.
    :param seq: A DNA sequence
    :return: A list of tuples (start, end, length). [(1, 3, 2), (4, 7, 3)]
    """
    if not seq:
        return []

    # run-length encode the sequence: every run, also the runs of length 1
    return [(m.start(), m.end(), m.end() - m.start())
            for m in re.finditer(r"([ACGT])\1*", seq.upper())]


def _cycle(read, query_pos):
    """
    Position in the read in sequencing orientation.
    :param read: AlignedSegment
    :param query_pos: position in the stored read (as pysam reports it)
    :return: the cycle, ie. the position as the sequencer read it
    """
    if read.is_reverse:
        return read.query_length - 1 - query_pos
    return query_pos


class SamCounter:

    def __init__(self, max_read_len, ranges_len, min_mapq=MINIMUM_MAPPING_QUALITY):
        """
        This class is in charge of counting the substitutions and homopolymers changes on a sam/bam file.

        All the counts are in the sequencing orientation: the position (cycle) is the position
        in the read as it was sequenced, and the bases are the sequenced bases. Reverse-strand
        reads are flipped, see _cycle.
        :param max_read_len: reads longer than this are not counted
        :param ranges_len: width of each read position range, in bases
        :param min_mapq: reads with a lower mapping quality are not counted
        """
        self.max_read_len = max_read_len
        self.min_mapq = min_mapq
        self.ranges = list(range(0, max_read_len + ranges_len, ranges_len))

        self.reads = 0
        self.skipped = 0

        # Substitutions #
        # keep a counter of changes per position for each possible
        # structure: np array [[[], [], [], []], ..]
        # each inner array represent a base and the a substitution
        # ie. [5][0][1] = 4, this represents 4 substitutions on the position 5 of an A (0) for a C (1)
        self.subs_per_base = np.zeros((max_read_len, 4, 4), dtype=np.int64)

        self.total_bases = 0
        # this count the total bases that didn't suffer a substitution, per position
        self.correct_bases = np.zeros((max_read_len, 4), dtype=np.int64)

        # Homopolymers #

        # total number of homopolymers that didn't suffer a modification
        self.correct_homopolymers = np.zeros((max_read_len, MAX_HOMOPOLYMER), dtype=np.int64)
        self.total_homopolymers = 0

        # keep a counter indels related to homopolymers per position
        # structure: np array [[[], [], [], []], ..]
        # each inner array represent a the change of an homopolymer of length index to value
        # ie. [3][3][4] = 11, this represents the change of an homopolymer of length 3 to 4
        #                     in the position 3
        #
        # Max homopolymer considered is 10
        self.indels_per_pos = np.zeros((max_read_len, MAX_HOMOPOLYMER, MAX_HOMOPOLYMER), dtype=np.int64)

    def count_homopolymers(self, read):
        """
        Count the homopolymer length changes on a read.

        For each run [start, end) of one base in the reference segment of the alignment
        (runs of length 1 included), the observed length is the number of run bases in the
        read between the bases aligned to the flanking reference bases start - 1 and end.
        This does not depend on where the aligner places the gap of an indel.

        A run is not counted when:
        - a flanking reference base is not aligned (deletion or end of the aligned part),
        - a base other than the run base is between the flanks (a substitution or a
          different indel, not a run length change),
        - the observed length is 0 (a fully deleted run has no observed symbol to model).
        :param read: AlignedSegment
        """
        try:
            ref_seq = read.get_reference_sequence()
        except ValueError as error:
            print("homopolymers not counted for {}: {}".format(read.query_name, error), file=sys.stderr)
            self.skipped += 1
            return

        if ref_seq is None or len(ref_seq) == 0:
            return

        ref_hom_positions = get_homopolymers(ref_seq)

        if len(ref_hom_positions) == 0:
            return

        try:
            aligned_positions = read.get_aligned_pairs()
        except ValueError as error:
            print("homopolymers not counted for {}: {}".format(read.query_name, error), file=sys.stderr)
            self.skipped += 1
            return

        read_seq = read.query_sequence
        # reference position of the first base of ref_seq, and a map ref pos -> read pos
        ref_start = read.reference_start
        ref_to_query = {ref_pos: query_pos for query_pos, ref_pos in aligned_positions
                        if query_pos is not None and ref_pos is not None}

        for start, end, length in ref_hom_positions:
            # flanking reference bases of the run
            left_ref = ref_start + start - 1
            right_ref = ref_start + end

            query_left = ref_to_query.get(left_ref)
            query_right = ref_to_query.get(right_ref)

            if query_left is None or query_right is None:
                # a flank is not aligned: deletion of the flank or end of the read
                continue

            # observed run bases: the read bases between the two flanks
            run_base = ref_seq[start]
            observed_length = 0
            other_base = False
            for query_pos in range(query_left + 1, query_right):
                base = read_seq[query_pos].upper()
                if base == run_base:
                    observed_length += 1
                else:
                    # a substitution or a different indel, not a run length change
                    other_base = True
                    break

            if other_base or observed_length == 0:
                continue

            # the run is at the first observed base; class 10 holds the lengths >= 10
            true_class = min(length, MAX_HOMOPOLYMER)
            observed_class = min(observed_length, MAX_HOMOPOLYMER)
            if read.is_reverse:
                # first base of the run in sequencing orientation is the last one stored
                cycle = _cycle(read, query_right - 1)
            else:
                cycle = _cycle(read, query_left + 1)

            self.total_homopolymers += 1
            if observed_class == true_class:
                self.correct_homopolymers[cycle][true_class - 1] += 1
            else:
                self.indels_per_pos[cycle][true_class - 1][observed_class - 1] += 1

    def count_substitutions(self, read):
        """
        Count the substitutions, and the substitutions base => base change.

        A position is only counted when the read base and the reference base are both in an
        M block and not next to an indel: the positions next to a homopolymer indel often
        show false mismatches from the gap placement.
        Positions with a base other than A, C, G or T (in the read or the reference) are skipped.
        :param read: AlignedSegment
        """
        try:
            aligned_positions = read.get_aligned_pairs(with_seq=True)
        except ValueError as error:
            print("substitutions not counted for {}: {}".format(read.query_name, error), file=sys.stderr)
            self.skipped += 1
            return

        read_seq = read.query_sequence

        # read positions next to an indel column are not counted
        near_indel = set()
        for idx, (read_pos, reference_pos, _) in enumerate(aligned_positions):
            if read_pos is None or reference_pos is None:
                for other in (idx - 1, idx + 1):
                    if 0 <= other < len(aligned_positions):
                        other_read_pos = aligned_positions[other][0]
                        if other_read_pos is not None:
                            near_indel.add(other_read_pos)

        for read_pos, reference_pos, base in aligned_positions:
            if read_pos is None or reference_pos is None or base is None:
                # indel column
                continue

            ref_base = base.upper()
            ref_idx = get_base_index(ref_base)
            if ref_idx == -1:
                # non-ACGT base in the reference
                continue

            read_base = read_seq[read_pos].upper()
            read_idx = get_base_index(read_base)
            if read_idx == -1:
                # uncalled or ambiguous base in the read
                continue

            if read_pos in near_indel:
                continue

            cycle = _cycle(read, read_pos)
            if read.is_reverse:
                # complement: A <-> T, C <-> G (the indexes are in ACGT order)
                ref_idx = 3 - ref_idx
                read_idx = 3 - read_idx

            if read_idx == ref_idx:
                # increase the correct base counter, ie base A was read as base A
                self.correct_bases[cycle][ref_idx] += 1
            else:
                # add a substitution on this position
                self.subs_per_base[cycle][ref_idx][read_idx] += 1
            self.total_bases += 1

    def run(self, sam_file):
        """
        Process each mapped read on the supplied SAM/BAM file path.
        :param sam_file: File Path
        """
        with pysam.AlignmentFile(sam_file) as sam:
            for read in sam.fetch(until_eof=True):
                self.count_read(read)

    def count_read(self, read):
        """
        Count the errors of one read, if it passes the filters.
        :param read: AlignedSegment
        """
        if read.is_unmapped or read.is_secondary or read.is_supplementary \
                or read.is_duplicate or read.is_qcfail:
            return

        # mapping quality is lower than the threshold
        if read.mapping_quality < self.min_mapq:
            return

        # reads longer than the max size aren't considered
        # (soft-clipped bases included: the positions are read cycles)
        if not 0 < read.query_length <= self.max_read_len:
            return

        self.count_substitutions(read)
        self.count_homopolymers(read)
        self.reads += 1
