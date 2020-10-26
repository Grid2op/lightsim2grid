# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
This module provide a function that can serve as a base replacement to the function newtonpf of
`pandapower/pypower/newtonpf.py`


"""
import numpy as np
from scipy import sparse

from lightsim2grid_cpp import SparseLUSolver
KLU_solver_available = False
try:
    from lightsim2grid_cpp import KLUSolver
    KLU_solver_available = True
except ImportError:
    pass


def newtonpf(Ybus, Sbus, V0, pv, pq, ppci, options):
    """
    Perform the Newton scheme to compute the AC powerflow of the system provided as input.
    It supports only one single slack bus.

    It is main as being integrated into pandapower as a replacement of the pypower implementation of "newtonpf"

    Parameters
    ----------
    Ybus: ``numpy.ndarray``, ``numpy.sparmatrix``, dtype:complex
        The admittance matrix. If not in a sparse CSC format, it will be converted to it.

    Sbus: ``numpy.ndarray``, dtype:complex
        The power injected at each bus.

    V0: ``numpy.ndarray``, dtype:complex
        The initial voltage

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
    # TODO have that in options maybe (can use GaussSeidel, and NR with KLU -faster- or SparseLU)
    if KLU_solver_available:
        solver = KLUSolver()
    else:
        solver = SparseLUSolver()
    Ybus = sparse.csc_matrix(Ybus)

    # do the newton raphson algorithm
    solver.solve(Ybus, V0, Sbus, pv, pq, max_it, tol)

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
