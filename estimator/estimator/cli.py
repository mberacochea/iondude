import argparse
import sys

from estimator.reporter import Reporter
from estimator.samcounter import MINIMUM_MAPPING_QUALITY, SamCounter


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Estimate the substitutions and homopolymers error rates.")
    parser.add_argument("sam", type=str, help="SAM/BAM file with the reads mapped against a reference")
    parser.add_argument("-l", type=int, default=200,
                        help="Max read length, default 200. Longer reads are not counted. You can get the longest "
                             "read length using awk: NR%%4==2{if (length($(0)) > l) l=length($(0))}END{ print l } "
                             "reads.fq")
    parser.add_argument("-r", type=int, default=20,
                        help="Width of each read position range, default 20. The estimator writes one matrix for "
                             "each range, ie. -r 20 gives the ranges 0-20, 20-40 and so on. Use the same value "
                             "for the -r option of iondude.")
    parser.add_argument("-q", type=int, default=MINIMUM_MAPPING_QUALITY,
                        help=f"Minimum mapping quality, default {MINIMUM_MAPPING_QUALITY}")
    parser.add_argument("-o", type=str,
                        help="Output file, otherwise the matrices are printed on stdout")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    counter = SamCounter(args.l, args.r, min_mapq=args.q)
    counter.run(args.sam)

    Reporter(counter).save(args.o)

    print(f"reads counted: {counter.reads}, reads with errors in the alignment: {counter.skipped}, "
          f"bases: {counter.total_bases}, homopolymers: {counter.total_homopolymers}", file=sys.stderr)


if __name__ == "__main__":
    main()
