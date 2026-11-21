#! /usr/bin/env python

import argparse
import sys

import numpy as np
from scipy import stats
from sympy import Matrix, Symbol, solve, symbols

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from estimator.matrices import get_slot, read_matrices
from estimator.seq import complement

BASE_INDEX = {
    'A': 0,
    'a': 0,
    'C': 1,
    'c': 1,
    'G': 2,
    'g': 2,
    'T': 3,
    't': 3
}

BASE = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
}

BASE_INDEX_READS = {
    'A': 11,
    'a': 11,
    'C': 12,
    'c': 12,
    'G': 13,
    'g': 13,
    'T': 14,
    't': 14
}

LOSS = np.matrix([
    [0, 1, 1, 1],
    [1, 0, 1, 1],
    [1, 1, 0, 1],
    [1, 1, 1, 0]
])

# max Q considered, past that value the DUDE has no possibility to correct
# or position, if the matrices are per position
MAX_VALUE = 30


def load_matrices(path):
    """
    Read a matrices file with one slot for each quality score.
    :param path: the JSON file
    :return: the matrices, see estimator.matrices
    """
    matrices = read_matrices(path)
    if matrices['slot_type'] != 'quality':
        sys.exit(f"{path}: the matrices must be by quality score (slot_type quality), "
                 f"not by {matrices['slot_type']}")
    return matrices


def subs_matrix(matrices, name, quality):
    """
    Get a substitutions matrix of a quality score.
    :param matrices: the matrices, see estimator.matrices
    :param name: channel or inverse
    :param quality: the quality score
    :return: a 4 x 4 numpy matrix
    """
    return np.matrix(get_slot(matrices, 'substitutions', quality)[name])


def get_base_index(base: str) -> int:
    return BASE_INDEX[base.decode('UTF-8')]


def main():
    """
    columns
    read name - context_reads - context_ref -  ctx freq ref - ctx freq reads - ref base - read base - map reversed =>

    context freq on reference { A - C - G - T } =>
    context freq on reads     { A - C - G - T }

    on not_corrected.txt
    ref base is the correct one
    read base is the wrong one

    on fixed.txt
    ref base and read base are the same (duh)

    on broken.txt
    ref base is the correct one and read base is the broken

    not corrected - get the difference between the correct context and the incorrect one
    fixed         - get the difference between the correct base and the other values
    broken        - get the difference between the incorrect and the correct

    for each row calculate (the list value that depends on the case)
    stats
    get the average
    get the median
    get min-max value (and the count for each case)
    """

    arg_parser = argparse.ArgumentParser(description='Print a report summary.')

    arg_parser.add_argument('-c', '--corrected', type=str, help='Corrected data from gain_report', required=True)
    arg_parser.add_argument('-n', '--not_corrected', type=str, help='Not corrected data from gain_report',
                            required=True)
    arg_parser.add_argument('-m', '--matrix', type=str, help='Channel Matrix', required=True)
    arg_parser.add_argument('-k', type=int, help='kmer size', required=True)
    # chart data
    # arg_parser.add_argument('-t', '--title', type=str, help='Chart title', required=True)
    # arg_parser.add_argument('-c', '--chart', type=str, help='Chart filename', required=True)

    args = arg_parser.parse_args()

    # read the matrices
    matrices = load_matrices(args.matrix)

    # pA threshold
    threshold = calculate_pA_threshold(matrices)

    # data
    corrected = load_data(args.corrected, args.k)
    not_corrected = load_data(args.not_corrected, args.k)

    yc, xc, fc = process_data(corrected, threshold)

    ync, xnc, fnc = process_data(not_corrected, threshold)

    params = {
        'axes.labelsize': 8,
        'text.fontsize': 8,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'text.usetex': False
    }

    mpl.rcParams.update(params)

    plt.figure(1)

    plt.grid()
    plt.plot(xc, yc, '.', color='#B22400')
    plt.plot(xnc, ync, '.', color='#006BB2')
    plt.xlabel('Base position')
    plt.ylabel(r'$\widehat{pA}$')

    lg = plt.legend(['Corrected', 'Not corrected'], loc=2)
    frame = lg.get_frame()
    frame.set_facecolor('0.9')
    frame.set_edgecolor('0.9')

    plt.ylim([0.4, 1.1])

    # threshold line
    plt.plot([_ for _ in range(1, MAX_VALUE, 1)], threshold, 'k-')

    plt.savefig('err_pos_merged.png')

    plt.figure(2)
    plt.subplot(121)
    plt.hist(fc, color='#B22400')

    plt.xlabel('Context count')
    plt.ylabel('Frequency')
    plt.title('Corrected contexts')

    plt.subplot(122)
    plt.hist(fnc, color='#006BB2')

    plt.xlabel('Context count')
    plt.ylabel('Frequency')
    plt.title('Not corrected contexts')

    plt.savefig('err_pos_merged_hist.png')

    exit(0)


def process_data(data, idx=16):
    """
    Process a gain_report verbose report
    :param data: Data from a report file
    :param idx: index of the column { Quality or Position }
    :return: y, y^1, and qualities|positions on the third parameter as x for the chart
    """
    x = []
    pAs = []
    freq_ctx = []

    for row in data:
        # analyze the context that appear just once on the reference
        ctx_ref_sum = sum([row[7], row[8], row[9], row[10]])
        if ctx_ref_sum != 1:
            continue

        freq = [row[11], row[12], row[13], row[14]]

        if sum(freq) == 0:
            continue

        # check that there are just 2 possibilities for this context
        count_check = 0
        for _i in freq:
            if _i > 0:
                count_check += 1

        if count_check != 2:
            continue

        ref_base = complement(row[5].decode('UTF-8')) if row[15] is True else row[5].decode('UTF-8')
        read_base = complement(row[17].decode('UTF-8')) if row[15] is True else row[17].decode('UTF-8')

        # I am considering only the changes A => C
        # that is, I read C but it is A, so I have to find
        if not (ref_base == 'A' and read_base == 'C') or row[idx] >= MAX_VALUE:
            continue

        # if the correct base doesn't appear on the read - then do not count (it's not possible to correct this)
        if row[BASE_INDEX_READS[ref_base]] == 0:
            continue

        # fA / (fA + fC + fG + fT)
        pA_empirical = row[11] / sum([row[11], row[12], row[13], row[14]])

        pAs.append(pA_empirical)
        x.append(row[idx])
        freq_ctx.append(sum([row[11], row[12], row[13], row[14]]))

    return pAs, x, freq_ctx


def load_data(data_file, k_size):
    """
    Numpy load data
    :param data_file: the file with the data
    :param k_size: kmer size
    :return: the data
    """
    return np.loadtxt(data_file, dtype={'names': ('read_name',
                                                  'context',
                                                  'context_ref',
                                                  'freq_ref',
                                                  'freq_reads',
                                                  'base_ref',        # 5 correct base { on not corrected }
                                                  'base_read',       # 6 incorrect base { on not corrected }
                                                  'freq_ctx_A_ref',   # 7
                                                  'freq_ctx_C_ref',   # 8
                                                  'freq_ctx_G_ref',   # 9
                                                  'freq_ctx_T_ref',   # 10
                                                  'freq_ctx_A_reads',   # 11
                                                  'freq_ctx_C_reads',   # 12
                                                  'freq_ctx_G_reads',   # 13
                                                  'freq_ctx_T_reads',   # 14
                                                  'reverse',            # 15
                                                  'read_quality',       # 16
                                                  'base_read_not_corrected', # 17
                                                  'position'),  # 18
                                        'formats': ('S20',
                                                    'S' + str(k_size),
                                                    'S' + str(k_size),
                                                    np.int_,
                                                    np.int_,
                                                    'S1',
                                                    'S1',
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.int_,
                                                    np.bool_,
                                                    np.int_,
                                                    'S1',
                                                    np.int_)})


def dude_rule(correct_base, read_base, counter, quality, matrices):
    """
    Dude RULE math
    :param correct_base:
    :param read_base:
    :param counter:
    :param quality:
    :param matrices:
    """
    # counters for the modified contexts
    #                   A  C  G  T
    #                     counter
    read_base = read_base.decode('UTF-8')
    ridx = BASE_INDEX[read_base]
    cidx = BASE_INDEX[correct_base]

    # copy
    counter_cp = counter[:]

    for j in range(0, 1000):

        # increment
        counter_cp[cidx] += j

        counter_arr = np.array(counter_cp)

        # INVERTED
        first_op = counter_arr.dot(subs_matrix(matrices, 'inverse', quality))

        # second op counters_inv * shur column product
        # the 2 depends on the observed base

        # loss over putative - channel over the observed
        min_value = sys.float_info.max

        # for each possible base apply the dude rule to get the base that minimize
        base = 'A'
        t = []
        for i in range(4):
            channel = subs_matrix(matrices, 'channel', quality)
            pairwaise_prod = np.array(LOSS[:, i]) * np.array(channel[:, ridx])
            _m = first_op.dot(pairwaise_prod)[0, 0]
            t.append(_m)
            if _m < min_value:
                min_value = _m
                base = BASE[i]

        if correct_base == base:
            # difference after being adjusted
            return (counter[cidx] + j) / (sum(counter) - counter[cidx] + j)

    return None


def dude_rule_qual(correct_base, read_base, counter, quality, matrices):
    """
    FIXME this is a copy of dude_rule but iterates over the quality
    """
    # counters for the modified contexts
    #                   A  C  G  T
    #                     counter
    read_base = read_base.decode('UTF-8')
    ridx = BASE_INDEX[read_base]
    counter_arr = np.array(counter)

    for q in range(49, 1, -1):

        # INVERTED
        first_op = counter_arr.dot(subs_matrix(matrices, 'inverse', q))

        # second op counters_inv * shur column product
        # the 2 depends on the observed base

        # loss over putative - channel over the observed
        min_value = sys.float_info.max

        # for each possible base apply the dude rule to get the base that minimize
        base = 'A'
        t = []
        for i in range(4):
            pairwaise_prod = np.array(LOSS[:, i]) * np.array(subs_matrix(matrices, 'channel', q)[:, ridx])
            _m = first_op.dot(pairwaise_prod)[0, 0]
            t.append(_m)
            if _m < min_value:
                min_value = _m
                base = BASE[i]

        if correct_base == base:
            return q

    return None


def describe_data(title, data):
    """
    Print basic stats data
    :param title: data tile
    :param data: array or np.array
    """
    if type(data) is not np.array:
        np_data = np.array(data)
    else:
        np_data = data

    if len(np_data) == 0:
        return

    print('Stats: ' + title)
    print('\tDescriptive stats: ', stats.describe(np_data))
    print('\tStd              : ', np.std(np_data))
    print('\tMode             : ', stats.mode(np_data))
    print('\tMedian           : ', np.median(np_data))
    print('=' * 180)


def calculate_pA_threshold(matrices):
    """
    For each quality, the min pA values to keep A given A
    :return: A list with the pA min required to fix for each quality value
    """

    ret_values = []

    # loss function
    loss_fn = Matrix(4, 4, lambda i, j: 0 if i == j else 1)

    for Q in range(1, MAX_VALUE, 1):
        # noisy symbol - that is, the one I read in the read
        # currently I am checking which pA I need for bases that are A but are read as C
        z_i = 1
        n_A, n_C, n_G, n_T = symbols('n_A n_C n_G n_T')

        index = {
            n_A: 0,
            n_C: 1,
            n_G: 2,
            n_T: 3
        }

        # counters
        mT = Matrix([[n_A, n_C, n_G, n_T]])

        # pprint(mT)

        pi = Matrix(subs_matrix(matrices, 'channel', Q))
        pi_1 = Matrix(subs_matrix(matrices, 'inverse', Q))

        values = []

        # pprint(pi)
        # pprint(pi_1)

        # ctx = {n_A: 945, n_C: 59, n_G: 0, n_T: 0}
        # #
        # # Keep A given that C was seen
        # print(((mT * pi_1) * (loss_fn.col(0).multiply_elementwise(pi.col(1)))).subs(ctx))
        #
        # # Keep C given that C was seen
        # print(((mT * pi_1) * (loss_fn.col(1).multiply_elementwise(pi.col(1)))).subs(ctx))
        # print(((mT * pi_1) * (loss_fn.col(2).multiply_elementwise(pi.col(1)))).subs(ctx))
        # print(((mT * pi_1) * (loss_fn.col(3).multiply_elementwise(pi.col(1)))).subs(ctx))
        #
        # exit(0)

        for x in range(4):
            e1 = (mT * pi_1)
            e2 = (loss_fn.col(x).multiply_elementwise(pi.col(z_i)))
            # pprint(e1)
            # pprint(e2)
            rule = e1 * e2
            inner_values = [0, 0, 0, 0]
            for symb in rule[0].args:
                try:
                    v = symb.args[0]
                    s = symb.args[1]
                except IndexError:
                    print('IndexError')

                if type(s) is not Symbol:
                    raise ValueError('Not a Symbol!')

                inner_values[index[s]] = v

            values.append(inner_values)

        ###################################################################################

        # min calc

        # x = A and z = z_i
        pA = symbols('pA')

        a_A = values[0][0]
        a_C = values[1][0]

        c_A = values[0][1]
        c_C = values[1][1]

        res = ((c_A - c_C) / (a_C - a_A))
        val_1 = solve((pA / (1 - pA)) - res, pA)[0]

        ################################################

        # x = C and z = z_i
        pA = None
        pA = symbols('pA')

        a_A = values[0][0]
        a_G = values[2][0]

        g_A = values[0][2]
        g_G = values[2][2]

        res = ((g_A - g_G) / (a_G - a_A))
        val_2 = solve((pA / (1 - pA)) - res, pA)[0]  # noqa: F841 - see the min() below

        #################################################

        # x = G and z = z_i
        pA = None
        pA = symbols('pA')

        a_A = values[0][0]
        a_T = values[3][0]

        t_A = values[0][3]
        t_T = values[3][3]

        res = ((t_A - t_T) / (a_T - a_A))
        val_3 = solve((pA / (1 - pA)) - res, pA)[0]  # noqa: F841 - see the min() below

        # ret_values.append(round(min([val_1, val_2, val_3]), 3))
        ret_values.append(round(val_1, 3))

    return ret_values


if __name__ == "__main__":
    main()
