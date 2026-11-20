#! /usr/bin/env python

import argparse

import matplotlib.pyplot as plt
import numpy as np
import numpy


def log_values(values):
    if not args.log:
        return values
    return np.where(values != 0, np.log(values), 0)

plt.style.use('ggplot')

parser = argparse.ArgumentParser()
parser.add_argument('--title', type=str, default=None)
parser.add_argument('--xtitle', type=str, default='k-mer abundance')
parser.add_argument('--y1title', type=str, default='# of k-mers with that abundance')
parser.add_argument('--y2title', type=str, default='# of changes')
parser.add_argument('--log', type=bool, default=False)

parser.add_argument('--xmin', type=int, default=None)
parser.add_argument('--xmax', type=int, default=None)
parser.add_argument('--ymin', type=int, default=None)
parser.add_argument('--ymax', type=int, default=None)
parser.add_argument('--delimiter', type=str, default='\t')

parser.add_argument('dist_file')
# post process
parser.add_argument('dist_file_post')
parser.add_argument('image_output')

args = parser.parse_args()

fig = plt.figure()
ax = fig.add_subplot(111)

data = numpy.loadtxt(args.dist_file, delimiter=args.delimiter)
data_post = numpy.loadtxt(args.dist_file_post, delimiter=args.delimiter)

# axis
data_values = log_values(data[:, 1])
ax.plot(data[:, 0], data_values, 'b-')
ax.set_ylabel(args.y1title if not args.log else 'log(' + args.y1title + ')')
ax.tick_params('y', colors='b')

ax_post = ax.twinx()

data_post_values = data_post[:, 1] # log_values(data_post[:, 1])

ax_post.bar(data[:, 0], data_post_values, 0.30, color='r')
ax_post.set_ylabel(args.y2title)
ax_post.tick_params('y', colors='r')

# title an footer note
if args.title:
    plt.title(args.title)

ax.annotate('File: {}'.format(args.dist_file_post),
            xy=(0.1, 0),
            xytext=(0, 10),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=6, ha='center', va='bottom')

ax_post.axis(xmin=args.xmin, ymin=args.ymin, xmax=args.xmax, ymax=args.ymax)

plt.savefig(args.image_output)
