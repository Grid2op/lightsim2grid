# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
from scipy import sparse

from lightsim2grid.solver import DCSolver

try:
    from lightsim2grid.solver import KLUDCSolver
    KLU_solver_available = True
except ImportError:
    KLU_solver_available = False


def _isolate_slack_ids(Sbus, pv, pq):
    # extract the slack bus
    ref = set(np.arange(Sbus.shape[0])) - set(pv) - set(pq)
    ref = np.array(list(ref))
    # build the slack weights
    slack_weights = np.zeros(Sbus.shape[0])
    slack_weights[ref] = 1.0 / ref.shape[0]
    return ref, slack_weights


def _get_valid_solver(options, Bbus):
    # initialize the solver
    # TODO have that in options maybe (can use GaussSeidel, and NR with KLU -faster- or SparseLU)
    solver = KLUDCSolver() if KLU_solver_available else DCSolver()

    if not sparse.isspmatrix_csc(Bbus):
        Bbus = sparse.csc_matrix(Bbus)

    if not Bbus.has_canonical_format:
        Bbus.sum_duplicates()
        if not Bbus.has_canonical_format:
            raise RuntimeError("Your matrix should be in a canonical format. See "
                            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.has_canonical_format.html"
                            " for more information.")
    
    return solver


def dcpf(B, Pbus, Va0, ref, pv, pq):
    """
    Implementation (using lightsim2grid) of the `dcpf` function of pandapower
    that allows to perform DC powerflow.
    
    It is meant to be used from inside pandapower and not directly.
    """
    # initialize the solver and perform some sanity checks
    solver = _get_valid_solver(None, B)
    
    # do the newton raphson algorithm
    slack_weights = np.ones(len(ref), dtype=float)
    slack_weights /= slack_weights.sum()
    solver.solve(B, np.cos(Va0) + 1j * np.sin(Va0), Pbus, ref, slack_weights, pv, pq, 1, 1e-8)
    # extract the results
    Va = solver.get_Va()
    return Va
