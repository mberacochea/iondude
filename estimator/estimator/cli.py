import argparse
import sys
from datetime import date
from importlib.metadata import version

from estimator.reporter import Reporter
from estimator.samcounter import MINIMUM_MAPPING_QUALITY, SamCounter

TOOL = "iondude-estimator"

# the arguments that are not settings of the estimation, all the others go to the "estimation" field of the file
NOT_SETTINGS = ("sam", "output")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Estimate the substitutions and homopolymers error rates.")
    parser.add_argument("sam", type=str, help="SAM/BAM file with the reads mapped against a reference")
    parser.add_argument("-l", dest="max_read_length", type=int, default=200,
                        help="Max read length, default 200. Longer reads are not counted. You can get the longest "
                             "read length using awk: NR%%4==2{if (length($(0)) > l) l=length($(0))}END{ print l } "
                             "reads.fq")
    parser.add_argument("-r", dest="range", type=int, default=20,
                        help="Width of each read position range, default 20. The estimator writes one matrix for "
                             "each range, ie. -r 20 gives the ranges 0-20, 20-40 and so on. iondude reads the "
                             "ranges from the file.")
    parser.add_argument("-q", dest="min_mapping_quality", type=int, default=MINIMUM_MAPPING_QUALITY,
                        help=f"Minimum mapping quality, default {MINIMUM_MAPPING_QUALITY}")
    parser.add_argument("-o", dest="output", type=str,
                        help="Output JSON file, otherwise the matrices are printed on stdout")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    counter = SamCounter(args.max_read_length, args.range, min_mapq=args.min_mapping_quality)
    counter.run(args.sam)

    source = {"tool": TOOL, "version": version(TOOL), "date": date.today().isoformat(), "input": args.sam}
    estimation = {name: value for name, value in vars(args).items() if name not in NOT_SETTINGS}
    Reporter(counter, source, estimation).save(args.output)

    print(f"reads counted: {counter.reads}, reads with errors in the alignment: {counter.skipped}, "
          f"bases: {counter.total_bases}, homopolymers: {counter.total_homopolymers}", file=sys.stderr)


if __name__ == "__main__":
    main()
