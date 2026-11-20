#! /usr/bin/env python

import argparse

import numpy as np
from tqdm import tqdm
from estimator.readfq import Seq
from estimator.utils import print_t, get_kmers


arg_parser = argparse.ArgumentParser(description='Save the kmers histogramas.')

arg_parser.add_argument('-f', '--fa', type=str, help='Reference fasta', required=True)
arg_parser.add_argument('-q', '--fq', type=str, help='Reads', required=True)
arg_parser.add_argument('-k', type=int, help='kmer size', required=True)
arg_parser.add_argument('-o', '--output', type=str, help='file name prefix', required=True)
arg_parser.add_argument('-v', '--verbose', action='store_true', help='verbose')
arg_parser.add_argument('-p', '--print', action='store_true', help='print the dist files?')

args = arg_parser.parse_args()

#############################################
#        count the kmers of size k          #
#############################################

print('\nkmers of the reference starting')
kmers_ref = {}
reference = open(args.fa, 'r')
for read in tqdm(Seq.readfq(reference)):
    for _kref in get_kmers(read.seq, args.k):
        if _kref is not None:
            kmers_ref[_kref] = kmers_ref.get(_kref, 0) + 1
reference.close()

# TODO reverse complementary
# join reverse complementaries
# for k, v in kmers_ref.items():
#     rck = Seq.revcomp(k)
#     if rck in kmers_ref:
#         kmers_ref[rck] = kmers_ref[rck] + v
#         p.append(k)
# for km in p:
#     del kmers_ref[km]

print('\nkmers of the reads starting')
kmers_reads = {}
reads = open(args.fq, 'r')
for read in tqdm(Seq.readfq(reads)):
    for _kread in get_kmers(read.seq, args.k):
        if _kread is not None:
            kmers_reads[_kread] = kmers_reads.get(_kread, 0) + 1
reads.close()

print('\nkmers done')

###########################################################
# print the distribution

# count - number of kmers that appear {count} times
# first 2 cols { col n kmers } for the reads, same for the genome
kf_unique, kf_counts = np.unique(np.array(list(kmers_ref.values())), return_counts=True)
freq_kfg = dict(zip(kf_unique, kf_counts, strict=True))

kr_unique, kr_counts = np.unique(np.array(list(kmers_reads.values())), return_counts=True)
freq_reads = dict(zip(kr_unique, kr_counts, strict=True))

# TODO implement
if args.verbose:
    pass
    # temp_ref = open(args.output + '_raw_counters_ref.txt', 'w')
    # temp_reads = open(args.output + '_raw_counters_reads.txt', 'w')

# for k, v in kmers_reads.items():
#     if args.verbose:
#         print_t(temp_reads, k, str(v))
#     freq_kfr[v] += 1

if args.verbose:
    pass
    # temp_ref.close()
    # temp_reads.close()

if args.print:
    genome_dist = open(args.output + '_genome_dist.txt', 'w')
    reads_dist = open(args.output + '_reads_dist.txt', 'w')

    # print abundance \t count
    for _k, _v in freq_kfg.items():
        print_t(genome_dist, _k, _v)

    for _k, _v in freq_reads.items():
        print_t(reads_dist, _k, _v)

    genome_dist.close()
    reads_dist.close()

print('kmers distribution files printed\n')

print('saving the files in npy format\n')

np.save('{}_kmers_ref.npy'.format(args.output), kmers_ref)
np.save('{}_kmers_reads.npy'.format(args.output), kmers_reads)

print("done")

exit(0)
