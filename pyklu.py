# Copyright (C) 2015 kagami-c
# Licensed under LGPL-2.1

'''
Please use this script as a module.

Detailed usage:

1. Import solver in PyKLU

  >>> from pyklu import solve_linear_system

2. Use solver by inputting NumPy ndarray

  >>> solver_linear_system(n, Ap, Ai, Ax, b)

n, Ap, Ai, Ax, b are numpy.ndarray object in this statement.

More details about these inputs are in the Manual of KLU Library
in SuiteSparse(http://faculty.cse.tamu.edu/davis/suitesparse.html)
'''

import ctypes
import numpy
import scipy.sparse as sparse

__all__ = ['solve_linear_system']

if __name__ == '__main__':
    print(__doc__)

libklu = ctypes.cdll.LoadLibrary('libpyklu.so')


def solve_linear_system(A, b):
    check_inumpyut_matrix(A, b)
    compress_A = sparse.csc_matrix(A)
    Ap = compress_A.indptr
    Ai = compress_A.indices
    Ax = compress_A.data
    check_compress_arrays(Ap, Ai, Ax)
    n = Ap.size - 1
    c_Ap = numpy.ctypeslib.as_ctypes(Ap)
    c_Ai = numpy.ctypeslib.as_ctypes(Ai)
    c_Ax = numpy.ctypeslib.as_ctypes(Ax)
    c_b = numpy.ctypeslib.as_ctypes(b)
    libklu.solve_linear_system(n, c_Ap, c_Ai, c_Ax, c_b)
    return b  # c_b shares the same memory space with b


def check_inumpyut_matrix(A, b):
    check_type_is_ndarray(A)
    check_type_is_ndarray(b)
    check_dimension(A, 2)  # left matrix is 2-D
    check_dimension(b, 1)  # right hand is 1-D
    check_dtype_is_double(A)
    check_dtype_is_double(b)


def check_type_is_ndarray(matrix):
    if not isinstance(matrix, numpy.ndarray):
        raise TypeError('Matrix should be numpy.ndarray type, but now {}'
                        .format(matrix.__class__))


def check_dimension(ndarray, dimension):
    if ndarray.ndim != dimension:
        raise TypeError('Matrix dimension should be {}, but now {}'
                        .format(dimension, ndarray.ndim))


def check_dtype_is_double(ndarray_in):
    for element in ndarray_in.flat:
        if not isinstance(element, numpy.double):
            raise TypeError('Element in ndarray should be {}, but now {}'
                            .format(numpy.double, element.__class__))


def check_compress_arrays(Ap, Ai, Ax):
    '''validate the result of scipy.sparse.csc_matrix'''
    assert isinstance(Ap, numpy.ndarray)
    assert isinstance(Ai, numpy.ndarray)
    assert isinstance(Ax, numpy.ndarray)
    assert Ap[-1] == Ai.size == Ax.size
