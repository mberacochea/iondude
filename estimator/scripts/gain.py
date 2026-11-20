#! /usr/bin/env python

import argparse
import sys

import numpy as np
import pysam
from tqdm import tqdm

from estimator.seq import complement, reverse_complement
from estimator.utils import print_t, kmer_for_base, get_context, get_reference_context

arg_parser = argparse.ArgumentParser(description='Compare the correction between 2 sam files.')

arg_parser.add_argument('-p', '--pre', type=str, help='SAM/BAM file. Must be sorted by queryname', required=True)
arg_parser.add_argument('-a', '--fna', type=str, help='Reference fasta', required=True)
arg_parser.add_argument('-t', '--post', type=str, help='SAM/BAM file of the reads post-correction', required=True)
arg_parser.add_argument('-f', '--kf', type=str, help='Reference kmers', required=True)
arg_parser.add_argument('-q', '--kr', type=str, help='Reads kmers', required=True)
arg_parser.add_argument('-k', type=int, help='kmer size', required=True)
arg_parser.add_argument('-d', '--denoiser', type=str, help='denoiser software', required=True)
arg_parser.add_argument('-c', '--chart', type=str, help='Chart title', required=False)
arg_parser.add_argument('-v', '--verbose', action='store_true', help='Verbose file output')
arg_parser.add_argument('-r', '--run', action='store_true', help='Run and save the values', default=True)

args = arg_parser.parse_args()


def save_log_files(name, inserted, corrected, not_corrected):

    log_broken = open(name + '_broken.txt', 'w')
    log_fixed = open(name + '_fixed.txt', 'w')
    log_not_corrected = open(name + '_not_corrected.txt', 'w')

    # print the inserted errors per context freq
    for _i in range(len(inserted)):
        print_t(log_broken, _i, inserted[_i])

    # print the corrected errors per context freq
    for _c in range(len(corrected)):
        print_t(log_fixed, _c, corrected[_c])

    # print the non-corrected errors per context freq
    for _nc in range(len(not_corrected)):
        print_t(log_not_corrected, _nc, not_corrected[_nc])

    log_broken.close()
    log_fixed.close()
    log_not_corrected.close()

# just print the data
if False: # not args.run:
    inserted_errors_rs_k = np.load('inserted_errors_rs_k' + args.denoiser + '.npy')
    corrected_errors_rs_k = np.load('corrected_errors_rs_k' + args.denoiser + '.npy')
    not_corrected_errors_rs_k = np.load('not_corrected_errors_rs_k' + args.denoiser + '.npy')

    inserted_errors_rf_k = np.load('inserted_errors_rf_k' + args.denoiser + '.npy')
    corrected_errors_rf_k = np.load('corrected_errors_rf_k' + args.denoiser + '.npy')
    not_corrected_errors_rf_k = np.load('not_corrected_errors_rf_k' + args.denoiser + '.npy')

    # save_log_files(args.denoiser + '_rs', inserted_errors_rs_k, corrected_errors_rs_k, not_corrected_errors_rs_k)
    # save_log_files(args.denoiser + '_rf', inserted_errors_rf_k, corrected_errors_rf_k, not_corrected_errors_rf_k)


def get_errors(sam):
    """
    Count the errors of a sam/bam file, just the substitutions
    :param sam           : mapping file
    """
    errors = {}
    not_mapped = {}
    read_subs_zi = {}

    # t = 0
    for pre_aligned_read in tqdm(sam.fetch(until_eof=True)):

        if pre_aligned_read.is_unmapped:
            not_mapped[pre_aligned_read.query_name] = 1
            continue

        if pre_aligned_read.is_supplementary or pre_aligned_read.mapping_quality == 0:
            # BWA puts a mq 0 if mapped elsewhere
            continue

        # TEST skip INDELS
        cigar = pre_aligned_read.cigarstring
        if 'I' in cigar or 'D' in cigar or 'N' in cigar:
            continue

        try:
            aligned_positions = pre_aligned_read.get_aligned_pairs(with_seq=True)
        except ValueError as error:
            print(error, file=sys.stderr)
            print(pre_aligned_read, file=sys.stderr)
            sys.exit(1)

        read_seq = pre_aligned_read.query_sequence
        read_subs_zi[pre_aligned_read.query_name] = ['' for _ in range(pre_aligned_read.query_length)]
        read_errors = np.zeros(pre_aligned_read.query_length, dtype=int)

        for read_pos, _, base in aligned_positions:
            # lower bases represent a substitution
            if base in ['a', 'c', 'g', 't']:
                # this is a substitution
                # the base is the reference base
                read_base = read_seq[read_pos]
                ref_base = base.upper()
                pos = read_pos
                if pre_aligned_read.is_reverse:
                    read_base = complement(read_base)
                    ref_base = complement(ref_base)
                    # pos = ((pre_aligned_read.query_length - 1) - read_pos) - 1
                if read_base == ref_base:
                    print(pre_aligned_read, file=sys.stderr)
                    print('Base substitution incorrect, base substitution with same base', file=sys.stderr)
                    exit(0)
                # add a substitution on this position
                read_subs_zi[pre_aligned_read.query_name][pos] = read_base
                read_errors[pos] += 1

        errors[pre_aligned_read.query_name] = read_errors

    return errors, not_mapped, read_subs_zi

k_len = args.k * 2 + 1

#####################################
# load the kmers dist files (npy files)

kmers_ref = np.load(args.kf, allow_pickle=True).item()
kmers_reads = np.load(args.kr, allow_pickle=True).item()

# values
kmers_ref_items = list(kmers_ref.values())
kmers_reads_items = list(kmers_reads.values())

# get the errors of the post-sam
post_sam = pysam.AlignmentFile(filename=args.post)

# inserted, corrected errors and uncorrected errors per quality [0 to 60]
inserted_errors_q = np.zeros(60, dtype=int)
not_corrected_errors_q = np.zeros(60, dtype=int)
corrected_errors_q = np.zeros(60, dtype=int)

# inserted, corrected errors and uncorrected errors per kmer frequency
inserted_errors_rf_k = np.zeros(max(kmers_ref_items) + 1, dtype=int)
not_corrected_errors_rf_k = np.zeros(max(kmers_ref_items) + 1, dtype=int)
corrected_errors_rf_k = np.zeros(max(kmers_ref_items) + 1, dtype=int)

inserted_errors_rs_k = np.zeros(max(kmers_reads_items) + 1, dtype=int)
not_corrected_errors_rs_k = np.zeros(max(kmers_reads_items) + 1, dtype=int)
corrected_errors_rs_k = np.zeros(max(kmers_reads_items) + 1, dtype=int)

borders = 0
now_maps = 0

#############################################################
#                      histograms data                      #
#############################################################
print('\nbuilding histogram')

# log file
# count the errors per position
print('\nbuilding the error list from the uncorrected seqs')
pre_sam = pysam.AlignmentFile(filename=args.pre)
pre_errors, pre_not_mapped, pre_bases = get_errors(pre_sam)
pre_sam.close()

log_broken_verbose = open(args.denoiser + '_broken_verbose.txt', 'w')
log_fixed_verbose = open(args.denoiser + '_fixed_verbose.txt', 'w')
log_not_corrected_verbose = open(args.denoiser + '_not_corrected_verbose.txt', 'w')

# this apply just for ref sequences
# will print the ref kmers that weren't corrected (all the kmers of the read)
# and will print the times that kmers appear on thr reference
log_not_mapped = open(args.denoiser + '_not_mapped.txt', 'w')

print_t(log_not_mapped, 'read_name')

##############################################################
#                       Reference FASTA                      #
#############################################################
reference_fasta = pysam.Fastafile(args.fna)

# TODO consider if a read is not mapped anymore
for post_aligned_read in tqdm(post_sam.fetch(until_eof=True)):

    if post_aligned_read.is_unmapped:
        # check if this was not mapped before
        if post_aligned_read.query_name not in pre_not_mapped:
            # the denoiser made this file not map anymore
            print_t(log_not_mapped, post_aligned_read.query_name)
        continue

    read_name = post_aligned_read.query_name
    read_seq = post_aligned_read.query_sequence
    ref_name = post_aligned_read.reference_name

    if post_aligned_read.is_supplementary or post_aligned_read.mapping_quality == 0:
        # BWA puts mq = 0 if the reads map elsewhere
        continue

    # TEST skip INDELS
    cigar = post_aligned_read.cigarstring
    if 'I' in cigar or 'D' in cigar or 'N' in cigar:
        continue

    try:
        aligned_positions = post_aligned_read.get_aligned_pairs(with_seq=True)
    except ValueError:
        continue

    # read_qual = reads_qualities[read_name]
    read_len = post_aligned_read.query_length

    if read_len < (args.k * 2) + 1:
        continue

    len_errors = 0
    try:
        base_errors = pre_errors[post_aligned_read.query_name]
        len_errors = len(base_errors)
    except KeyError:
        now_maps += 1
        continue

    ins_errors = []
    del_errors = []

    # reversed then flip the qualities
    #if post_aligned_read.is_reverse:
    #     read_qual = read_qual[::-1]
    # if post_aligned_read.is_reverse:
    #    read_seq = reverse_complement(read_seq)

    position = 0
    ref_base = 'A'
    query_quality = post_aligned_read.query_qualities

    for read_pos, reference_pos, base in aligned_positions:

        if read_pos is None or reference_pos is None or base is None:
            continue

        position = read_pos
        ref_base = base.upper()
        if post_aligned_read.is_reverse:
            ref_base = complement(ref_base)

        if read_pos < args.k or read_pos > read_len - args.k - 1:
            borders += 1
            continue

        # compare against the PRE corrected SAM errors
        if len_errors != read_len or read_pos > len_errors:
            # indel correction or a base was moved, ignore for now
            continue
        else:
            was_incorrect = base_errors[read_pos]

        if base in ['a', 'c', 'g', 't'] or was_incorrect > 0:
            # noisy context - the whole context
            context = get_context(read_seq, position, args.k)
            read_base = read_seq[position]

            # reference context
            ctx_ref = get_reference_context(reference_fasta, ref_name, reference_pos, args.k)
            # this is used to print the value if reversed as the reference has only the forward strand
            # use use that to count the number of times that a context appears
            ctx_ref_rev = ctx_ref

            if post_aligned_read.is_reverse:
                try:
                    ctx_ref = reverse_complement(ctx_ref)
                except Exception as e:
                    print('non conventional base on reference - skipped', e, post_aligned_read)
                    continue

                context = reverse_complement(context)
                read_base = complement(read_base)

            # read context
            context_rs_freq = kmers_reads.get(context, 0)
            context_rf_freq = kmers_ref.get(ctx_ref_rev, 0)

        # lower bases represent a substitution
        if base in ['a', 'c', 'g', 't']:
            # this is a substitution
            if was_incorrect == 0:
                # inserted_errors_q[read_qual] += 1
                # introduced
                if args.verbose:
                    print_t(log_broken_verbose,
                            read_name,
                            context,
                            ctx_ref,
                            context_rf_freq,
                            context_rs_freq,
                            ref_base,
                            read_base,
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'A', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'C', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'G', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'T', args.k),
                            kmer_for_base(kmers_reads, context, 'A', args.k),
                            kmer_for_base(kmers_reads, context, 'C', args.k),
                            kmer_for_base(kmers_reads, context, 'G', args.k),
                            kmer_for_base(kmers_reads, context, 'T', args.k),
                            int(post_aligned_read.is_reverse),
                            query_quality[read_pos],
                            pre_bases[post_aligned_read.query_name][read_pos],
                            position)

                if ref_base == read_base:
                    print(post_aligned_read, context, position, read_pos, aligned_positions, sep='\n')
                    exit(0)

                inserted_errors_rs_k[context_rs_freq] += 1
                inserted_errors_rf_k[context_rf_freq] += 1
            elif was_incorrect > 0:
                if args.verbose:
                    print_t(log_not_corrected_verbose,
                            read_name,
                            context,
                            ctx_ref,
                            context_rf_freq,
                            context_rs_freq,
                            ref_base,
                            read_base,
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'A', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'C', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'G', args.k),
                            kmer_for_base(kmers_ref, ctx_ref_rev, 'T', args.k),
                            kmer_for_base(kmers_reads, context, 'A', args.k),
                            kmer_for_base(kmers_reads, context, 'C', args.k),
                            kmer_for_base(kmers_reads, context, 'G', args.k),
                            kmer_for_base(kmers_reads, context, 'T', args.k),
                            int(post_aligned_read.is_reverse),
                            query_quality[read_pos],
                            pre_bases[post_aligned_read.query_name][read_pos],
                            position)

                if ref_base == read_base:
                    print(post_aligned_read, context, position, read_pos, aligned_positions, sep='\n')
                    exit(0)

                # not corrected
                not_corrected_errors_rs_k[context_rs_freq] += 1
                not_corrected_errors_rf_k[context_rf_freq] += 1
        elif base is not None and was_incorrect > 0:
            # corrected_errors_q[read_qual] += 1
            if args.verbose:
                print_t(log_fixed_verbose,
                        read_name,
                        context,
                        ctx_ref,
                        context_rf_freq,
                        context_rs_freq,
                        ref_base,
                        read_base,
                        kmer_for_base(kmers_ref, ctx_ref_rev, 'A', args.k),
                        kmer_for_base(kmers_ref, ctx_ref_rev, 'C', args.k),
                        kmer_for_base(kmers_ref, ctx_ref_rev, 'G', args.k),
                        kmer_for_base(kmers_ref, ctx_ref_rev, 'T', args.k),
                        kmer_for_base(kmers_reads, context, 'A', args.k),
                        kmer_for_base(kmers_reads, context, 'C', args.k),
                        kmer_for_base(kmers_reads, context, 'G', args.k),
                        kmer_for_base(kmers_reads, context, 'T', args.k),
                        int(post_aligned_read.is_reverse),
                        query_quality[read_pos],
                        pre_bases[post_aligned_read.query_name][read_pos],
                        position)

            corrected_errors_rs_k[context_rs_freq] += 1
            corrected_errors_rf_k[context_rf_freq] += 1

reference_fasta.close()

log_fixed_verbose.close()
log_broken_verbose.close()
log_not_corrected_verbose.close()

# save for future use
if args.run:
    try:
        # np.save('inserted_errors_Q_' + args.denoiser, inserted_errors_q)
        # np.save('not_corrected_errors_Q_' + args.denoiser, not_corrected_errors_q)
        # np.save('corrected_errors_Q_' + args.denoiser, corrected_errors_q)

        np.save('inserted_errors_rs_k' + args.denoiser, inserted_errors_rs_k)
        np.save('corrected_errors_rs_k' + args.denoiser, corrected_errors_rs_k)
        np.save('not_corrected_errors_rs_k' + args.denoiser, not_corrected_errors_rs_k)

        np.save('inserted_errors_rf_k' + args.denoiser, inserted_errors_rf_k)
        np.save('corrected_errors_rf_k' + args.denoiser, corrected_errors_rf_k)
        np.save('not_corrected_errors_rf_k' + args.denoiser, not_corrected_errors_rf_k)
    except Exception as e:
        print(e)

save_log_files(args.denoiser + '_rs', inserted_errors_rs_k, corrected_errors_rs_k, not_corrected_errors_rs_k)
save_log_files(args.denoiser + '_rf', inserted_errors_rf_k, corrected_errors_rf_k, not_corrected_errors_rf_k)

post_sam.close()
