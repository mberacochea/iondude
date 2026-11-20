#!/usr/bin/env python
"""
Plot the substitutions per read position and the mean base quality per read position.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from estimator.readfq import Seq
from estimator.samcounter import SamCounter


def main():
    parser = argparse.ArgumentParser(description="Plot the errors and the base quality per read position.")
    parser.add_argument("sam", type=str, help="SAM/BAM file with the reads mapped against a reference")
    parser.add_argument("-f", type=str, required=True, help="Fastq file with the reads")
    parser.add_argument("-l", type=int, default=200, help="Max read length, default 200")
    parser.add_argument("-o", type=str, default="quality_errors.png", help="Output image")
    args = parser.parse_args()

    counter = SamCounter(args.l, args.l)
    counter.run(args.sam)

    # mean quality per position
    quality = np.zeros(args.l)
    bases_per_position = np.zeros(args.l)
    with open(args.f) as fastq_file:
        for read in Seq.readfq(fastq_file):
            if len(read) > args.l:
                continue
            for idx, q in enumerate(read.qual_values):
                quality[idx] += q
                bases_per_position[idx] += 1
    quality = np.divide(quality, bases_per_position, out=np.zeros(args.l), where=bases_per_position > 0)

    # substitutions per position
    substitutions = counter.subs_per_base.sum(axis=(1, 2))
    positions = np.arange(args.l)

    fig, axis_error = plt.subplots()
    axis_error.plot(positions, substitutions, "b-")
    axis_error.set_xlabel("position on read")
    axis_error.set_ylabel("errors", color="b")
    axis_error.tick_params(axis="y", colors="b")

    axis_quality = axis_error.twinx()
    axis_quality.plot(positions, quality, "r-")
    axis_quality.set_ylabel("base quality", color="r")
    axis_quality.tick_params(axis="y", colors="r")

    fig.savefig(args.o)


if __name__ == "__main__":
    main()
