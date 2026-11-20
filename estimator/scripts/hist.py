#! /usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', type=str)
parser.add_argument('-i', '--image_output', type=str)
parser.add_argument('-m', '--data_modified', type=str, required=False)
parser.add_argument('-s', '--separator', type=str, default='\t')

args = parser.parse_args()

fig, ax1 = plt.subplots()

ddist = np.loadtxt(args.data)

# ax1.semilogy(ddist[:, 0], ddist[:, 1], 'b-')
ax1.plot(ddist[:10, 0], ddist[:10, 1], 'b-')

ax1.set_xlabel('kmer frequency')

# ax1.set_ylabel('exp', color='b')
# ax1.tick_params('y', colors='b')

if args.data_modified:
    mdist = np.loadtxt(args.data_modified)
    ax2 = ax1.twinx()
    ax2.plot(mdist[5:, 0], mdist[5:, 1], 'r.')
    ax2.set_ylabel('modifications', color='r')
    ax2.tick_params('y', colors='r')

fig.tight_layout()

plt.savefig(args.image_output)
