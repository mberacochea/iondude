"""
Load the channel matrices produced by the estimator.
"""

import os

import numpy as np
from parse import search


class Matrix:

    def __init__(self, name):
        """
        Build a new matrix
        :param name: matrix name, the header line without the line break
        """
        self.name = name
        self.data = []

        self.min = 0
        self.max = 0
        s = search('{:d}-{:d}', name)
        if s:
            self.min = s[0]
            self.max = s[1]

        self.inverted = '#SFI' in name or '#IFI' in name
        self.is_float = '#SF' in name or '#IF' in name
        self.dim = 10 if '#I' in name else 4

    def add_row(self, data_line):
        """
        Add a new row to the matrix. The column names line is skipped.
        :param data_line: matrix line, the row name and the values, tab separated
        """
        fields = data_line.split()
        if len(fields) != self.dim + 1 or fields[0] == '-':
            return
        value_type = float if self.is_float else int
        self.data.append([value_type(value) for value in fields[1:]])

    def finish(self):
        """
        Convert the rows to a numpy matrix, rounded to 4 decimals.
        """
        self.data = np.matrix(np.round(np.array(self.data), decimals=4))


def read_matrix_file(file_name):
    """
    Load a set of matrices from a file
    :param file_name: file with the matrices
    :return: a dict { matrix name: matrix }
    """
    if not file_name:
        return None

    if not os.path.exists(file_name):
        raise FileNotFoundError(file_name)

    matrices = {}
    current_matrix = None

    with open(file_name, 'r') as matrix_file:
        for line in matrix_file:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                current_matrix = Matrix(line.strip())
                matrices[current_matrix.name] = current_matrix
                continue
            if current_matrix:
                current_matrix.add_row(line)

    for matrix in matrices.values():
        matrix.finish()

    return matrices
