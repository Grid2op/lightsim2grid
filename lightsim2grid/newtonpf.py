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
import warnings
import numpy as np
from scipy import sparse

from lightsim2grid_cpp import SparseLUSolver
KLU_solver_available = False
try:
    from lightsim2grid_cpp import KLUSolver
    KLU_solver_available = True
except ImportError:
    pass

def _isolate_slack_ids(Sbus, pv, pq):
    # extract the slack bus
    ref = set(np.arange(Sbus.shape[0])) - set(pv) - set(pq)
    ref = np.array(list(ref))
    # build the slack weights
    slack_weights = np.zeros(Sbus.shape[0])
    slack_weights[ref] = 1.0 / ref.shape[0]
    return ref, slack_weights

def newtonpf_old(Ybus, Sbus, V0, pv, pq, ppci, options):
    """
    Perform the Newton scheme to compute the AC powerflow of the system provided as input.
    It supports only one single slack bus.

    It is main as being integrated into pandapower as a replacement of the pypower implementation of "newtonpf"

    .. versionadded:: 0.5.6

    .. note::
        This is a legacy code mainly present for compatibility with older pandapower versions.

    .. warning::
        It considers that all nodes non pv, non pq are slack nodes and assign the same weights
        for all of them even though it was not possible to perform such computation in earlier versions.

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

    # extract the slack bus
    ref, slack_weights = _isolate_slack_ids(Sbus, pv, pq)

    # do the newton raphson algorithm
    solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)

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


def newtonpf(Ybus, Sbus, V0, ref, pv, pq, ppci, options):
    """
    Perform the Newton scheme to compute the AC powerflow of the system provided as input.
    It supports only one single slack bus.

    It is main as being integrated into pandapower as a replacement of the pypower implementation of "newtonpf"

    .. versionchanged:: 0.5.6

    .. note::
        It has been updated in version 0.5.6 to match pandapower new signature (addition of the `ref`
        parameter)

        If you want the old behaviour, please use the `newtonpf_old`function.
    
    Parameters
    ----------
    Ybus: ``numpy.ndarray``, ``numpy.sparmatrix``, dtype:complex
        The admittance matrix. If not in a sparse CSC format, it will be converted to it.

    Sbus: ``numpy.ndarray``, dtype:complex
        The power injected at each bus.

    V0: ``numpy.ndarray``, dtype:complex
        The initial voltage

    ref:
        Id of the slack buses (added in version 0.5.6 to match pandapower changes)

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
    bus = ppci['bus']

    try:
        from pandapower.pypower.idx_bus import SL_FAC  # lazy import for earlier pandapower version (without distributed slack)
        slack_weights = bus[:, SL_FAC].astype(float)  ## contribution factors for distributed slack
        slack_weights /= slack_weights.sum()
    except ImportError:
        # earlier version of pandapower
        warnings.warn("You are using a pandapower version that does not support distributed slack. We will attempt to "
                      "replicate this with lightsim2grid.")
        ref, slack_weights = _isolate_slack_ids(Sbus, pv, pq)

    # initialize the solver
    # TODO have that in options maybe (can use GaussSeidel, and NR with KLU -faster- or SparseLU)
    if KLU_solver_available:
        solver = KLUSolver()
    else:
        solver = SparseLUSolver()
    Ybus = sparse.csc_matrix(Ybus)

    # do the newton raphson algorithm
    solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)

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
