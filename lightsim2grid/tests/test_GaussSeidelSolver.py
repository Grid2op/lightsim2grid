# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import unittest
import numpy as np
import pdb
import zipfile
from scipy import sparse
import copy
import time

GaussSeidelSolver_AVAILBLE = False

try:
    from lightsim2grid_cpp import GaussSeidelSolver
    GaussSeidelSolver_AVAILBLE = True
except ImportError:
    # KLU solver is not available, these tests cannot be carried out
    pass


import sys

from numpy import linalg, conj, r_, Inf


def gausspf(Ybus, Sbus, V0, ref, pv, pq, ppopt=None):
    """Solves the power flow using a Gauss-Seidel method.
    Solves for bus voltages given the full system admittance matrix (for
    all buses), the complex bus power injection vector (for all buses),
    the initial vector of complex bus voltages, and column vectors with
    the lists of bus indices for the swing bus, PV buses, and PQ buses,
    respectively. The bus voltage vector contains the set point for
    generator (including ref bus) buses, and the reference angle of the
    swing bus, as well as an initial guess for remaining magnitudes and
    angles. C{ppopt} is a PYPOWER options vector which can be used to
    set the termination tolerance, maximum number of iterations, and
    output options (see C{ppoption} for details). Uses default options
    if this parameter is not given. Returns the final complex voltages,
    a flag which indicates whether it converged or not, and the number
    of iterations performed.
    @see: L{runpf}
    @author: Ray Zimmerman (PSERC Cornell)
    @author: Alberto Borghetti (University of Bologna, Italy)
    """

    ## default arguments
    if ppopt is None:
        ppopt = {"PF_TOL": 1e-8, "PF_MAX_IT_GS": 10000, "VERBOSE": False}

    ## options
    tol     = ppopt['PF_TOL']
    max_it  = ppopt['PF_MAX_IT_GS']
    verbose = ppopt['VERBOSE']

    ## initialize
    converged = 0
    i = 0
    V = V0.copy()
    #Va = angle(V)
    Vm = abs(V)

    ## set up indexing for updating V
    npv = len(pv)
    npq = len(pq)
    pvpq = r_[pv, pq]

    ## evaluate F(x0)
    mis = V * conj(Ybus * V) - Sbus
    F = r_[  mis[pvpq].real,
             mis[pq].imag   ]

    ## check tolerance
    normF = linalg.norm(F, Inf)
    if verbose > 1:
        sys.stdout.write('\n it    max P & Q mismatch (p.u.)')
        sys.stdout.write('\n----  ---------------------------')
        sys.stdout.write('\n%3d        %10.3e' % (i, normF))
    if normF < tol:
        converged = 1
        if verbose > 1:
            sys.stdout.write('\nConverged!\n')

    ## do Gauss-Seidel iterations
    while (not converged and i < max_it):
        ## update iteration counter
        i = i + 1

        ## update voltage
        ## at PQ buses
        for k in pq[list(range(npq))]:
            tmp = (conj(Sbus[k] / V[k]) - Ybus[k, :] * V) / Ybus[k, k]
            V[k] = V[k] + tmp[0]

        ## at PV buses
        if npv:
            for k in pv[list(range(npv))]:
                tmp = (V[k] * conj(Ybus[k,:] * V)).imag
                Sbus[k] = Sbus[k].real + 1j * tmp[0]
                tmp = (conj(Sbus[k] / V[k]) - Ybus[k, :] * V) / Ybus[k, k]
                V[k] = V[k] + tmp[0]
#               V[k] = Vm[k] * V[k] / abs(V[k])
            V[pv] = Vm[pv] * V[pv] / abs(V[pv])

        ## evalute F(x)
        mis = V * conj(Ybus * V) - Sbus
        F = r_[  mis[pv].real,
                 mis[pq].real,
                 mis[pq].imag  ]

        ## check for convergence
        normF = linalg.norm(F, Inf)
        if verbose > 1:
            sys.stdout.write('\n%3d        %10.3e' % (i, normF))
        if normF < tol:
            converged = 1
            if verbose:
                sys.stdout.write('\nGauss-Seidel power flow converged in '
                                 '%d iterations.\n' % i)

    if verbose:
        if not converged:
            sys.stdout.write('Gauss-Seidel power did not converge in %d '
                             'iterations.' % i)

    return V, converged, i


class MakeTests(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        unittest.TestCase.__init__(self, methodName=methodName)
        self.methodName = methodName

        self.max_it = 10000
        self.tol = 1e-5  # tolerance for the solver (good compromise between test speed and solver accuracy)
        self.tol_test = 1e-4  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)
        if not GaussSeidelSolver_AVAILBLE:
            return

        self.solver = GaussSeidelSolver()

        self.path = None
        self.V_init = None
        self.pq = None
        self.pv = None
        self.Sbus = None
        self.Ybus = None

    def load_array(self, myzip, nm):
        arr = myzip.extract(nm)
        res = np.load(arr)
        os.remove(arr)
        return res

    def load_path(self, path):
        try:
            self.path = path
            with zipfile.ZipFile(self.path) as myzip:
                self.V_init = self.load_array(myzip, "V0.npy")
                self.pq = self.load_array(myzip, "pq.npy")
                self.pv = self.load_array(myzip, "pv.npy")
                self.Sbus = self.load_array(myzip, "Sbus.npy")
                self.Ybus = self.load_array(myzip, "Ybus.npy")
                self.Ybus = sparse.csc_matrix(self.Ybus)
            return True
        except:
            return False

    def solver_aux(self):
        self.solver.reset()
        if self.Ybus.shape[0] != 14 and self.Ybus.shape[0] != 30 and self.Ybus.shape[0] != 118:
            # others take too long
            return

        ppopt = {"PF_TOL": self.tol, "PF_MAX_IT_GS": self.max_it, "VERBOSE": False}

        # start = time.perf_counter()
        ref = set(np.arange(self.Sbus.shape[0])) - set(self.pv) - set(self.pq)
        ref = np.array(list(ref))
        # build the slack weights
        slack_weights = np.zeros(self.Sbus.shape[0])
        slack_weights[ref] = 1.0 / ref.shape[0]

        has_conv = self.solver.compute_pf(self.Ybus, self.V_init, self.Sbus, 
                                          ref, slack_weights, self.pv,
                                          self.pq, self.max_it, self.tol)

        Vgs, convergedgs, max_itgs = gausspf(copy.deepcopy(self.Ybus),
                                             copy.deepcopy(self.Sbus),
                                             copy.deepcopy(self.V_init),
                                             ref=None, pv=self.pv, pq=self.pq, ppopt=ppopt)
        # end = time.perf_counter()
        # print("end - start: {}".format(end - start))
        # print("have i conv ? {}".format(has_conv))
        # print("print timer solve: {}".format(self.solver.get_timers()[1]))

        assert has_conv, "the load flow has diverged for {}".format(self.path)
        Va = self.solver.get_Va()
        Vm = self.solver.get_Vm()
        iter_max = self.solver.get_nb_iter()

        Va_pp = np.angle(Vgs)
        Vm_pp = np.abs(Vgs)
        assert max_itgs == iter_max, "wrong number of iteration"
        assert np.sum(np.abs(Va - Va_pp)) <= self.tol_test, "voltages angles are not the same"
        assert np.sum(np.abs(Vm - Vm_pp)) <= self.tol_test, "voltages magnitude are not the same"

    def test_dir(self):
        if not GaussSeidelSolver_AVAILBLE:
            self.skipTest("GaussSeidelSolver is not installed")
        nb_tested = 0
        for path in os.listdir("."):
            _, ext = os.path.splitext(path)
            if ext == ".zip":
                path_ok = self.load_path(path)
                if path_ok:
                    self.solver_aux()
                    nb_tested += 1
        assert nb_tested == 5, "incorrect number of test cases found, found {} while there should be 5".format(nb_tested)


if __name__ == "__main__":
    unittest.main()
