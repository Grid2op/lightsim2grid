"""
This module provide a function that can serve as a base replacement to the function newtonpf of
`pandapower/pypower/newtonpf.py`


"""
import numpy as np
from scipy import sparse
from pyklu_cpp import KLUSolver


def newtonpf(Ybus, V, Sbus, pv, pq, ppci, options):
    """
    Perform the Newton scheme to compute the AC powerflow of the system provided as input.
    It supports only one single slack bus.

    It is main as being integrated into pandapower as a replacement of the pypower implementation of "newtonpf"

    Parameters
    ----------
    Ybus: ``numpy.ndarray``, ``numpy.sparmatrix``, dtype:complex
        The admittance matrix. If not in a sparse CSC format, it will be converted to it.

    V: ``numpy.ndarray``, dtype:complex
        The initial voltage

    Sbus: ``numpy.ndarray``, dtype:complex
        The origin injection

    pv: ``numpy.ndarray``, dtype:np.int
        Index of the pv buses (slack bus must NOT be on this list)

    pq: ``numpy.ndarray``, dtype:np.int
        Index of the pq buses (slack bus must NOT be on this list)

    ppci: ``dict``
        pandapower internal "ppc", ignored.

    options: ``dict``
        Dictionnary of various pandapower option. Only "max_iteration" and "tolerance_mva" are used at the moment.

    Returns
    -------
    V: ``numpy.ndarray``, dtype:complex
        The final complex voltage vector

    converged: ``bool``
        Whether the powerflow converged or not

    iterations: ``int``
        The number of iterations the solver performed

    J: `numpy.sparmatrix``, dtype:float
        The csc scipy matrix of the jacobian matrix of the system.

    """
    max_it = options["max_iteration"]
    tol = options['tolerance_mva']
    # initialize the solver
    solver = KLUSolver()
    Ybus = sparse.csc_matrix(Ybus)

    # do the newton raphson algorithm
    solver.solve(Ybus, V, Sbus, pv, pq, max_it, tol)

    # extract the results
    Va = solver.get_Va()
    Vm = solver.get_Vm()
    V = Vm * np.exp(1j * Va)
    J = solver.get_J()
    converged = solver.converged()
    iterations = solver.get_nb_iter()

    Vm_it = None
    Va_it = None
    return V, converged, iterations, J, Vm_it, Va_it
