#! /usr/bin/env python

import argparse

arg_parser = argparse.ArgumentParser(description='Group the verbose output by context.')

arg_parser.add_argument('-f', '--file', type=str, help='TSV file with the verbose output', required=True)

args = arg_parser.parse_args()

# context (key: ATAGAT, value: counter)
context_counters = {}
data = open(args.file, 'r')

for line in data:
    # read, context, count reads, count ref
    values = line.split('\t')
    if values[1] not in context_counters:
        context_counters[values[1]] = [0, 0]
    if context_counters[values[1]][0] < int(values[2]):
        context_counters[values[1]][0] = int(values[2])
    if context_counters[values[1]][1] < int(values[3]):
        context_counters[values[1]][1] += int(values[3])

data.close()

for key, value in context_counters.items():
    print('{}\t{}\t{}\t'.format(key, str(value[0]), str(value[1])))
