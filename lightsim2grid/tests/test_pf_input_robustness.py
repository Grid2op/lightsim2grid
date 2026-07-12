# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Robustness tests for the powerflow entry points: `LSGrid.ac_pf` / `LSGrid.dc_pf`
and the per-solver `compute_pf` / `solve` methods, in the same spirit as the
binary-file corruption tests (test_binary_serialization.py) and the
out-of-bounds binding tests (test_LSGrid_out_of_bounds.py): ill-formed inputs
must raise a clean Python exception, never crash, read out of bounds, or
silently compute garbage.

Conventions (mirroring test_LSGrid_out_of_bounds.py):
  * ``RuntimeError`` for structural problems: mis-sized V / Sbus /
    slack_weights, non-square Ybus, empty slack_ids, a bus listed in more
    than one of slack/pv/pq, negative max_iter (0 is valid: it returns the
    pre-iteration state), non-positive or non-finite tol
    (``BaseAlgo::check_pf_inputs`` / ``check_iter_tol`` throw
    ``std::runtime_error``);
  * ``IndexError`` for out-of-range / negative bus ids in slack_ids / pv / pq
    (``std::out_of_range``);
  * ``TypeError`` for python arguments of the wrong type (pybind11).
"""

import unittest
import warnings
import numpy as np

import pandapower.networks as pn
from lightsim2grid.network import init_from_pandapower
from lightsim2grid import lightsim2grid_cpp

MAX_IT = 30
TOL = 1e-8

# every AC solver class importable in this build (DC solver objects only
# expose the throwing complex entry point: the DC validation is exercised
# through LSGrid.dc_pf and the C++ compute_pf_dc path instead)
AC_SOLVER_CLASSES = [
    nm for nm in (
        "NR_SparseLU", "NRSing_SparseLU",
        "FDPF_XB_SparseLU", "FDPF_BX_SparseLU",
        "NR_KLU", "NRSing_KLU", "FDPF_XB_KLU", "FDPF_BX_KLU",
        "NR_NICSLU", "NRSing_NICSLU", "FDPF_XB_NICSLU", "FDPF_BX_NICSLU",
        "NR_CKTSO", "NRSing_CKTSO", "FDPF_XB_CKTSO", "FDPF_BX_CKTSO",
        "GaussSeidelAlgo", "GaussSeidelSynchAlgo",
    ) if hasattr(lightsim2grid_cpp, nm)
]


def _make_grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        return init_from_pandapower(net)


class TestAcDcPfInputs(unittest.TestCase):
    """LSGrid.ac_pf / LSGrid.dc_pf with corrupted arguments."""

    def setUp(self):
        self.grid = _make_grid()
        self.n_bus = self.grid.total_bus()
        self.v_init = np.ones(self.n_bus, dtype=complex)

    def _pf_methods(self):
        yield "ac_pf", self.grid.ac_pf
        yield "dc_pf", self.grid.dc_pf

    def test_valid_calls_converge(self):
        for name, pf in self._pf_methods():
            with self.subTest(name):
                V = pf(1.0 * self.v_init, MAX_IT, TOL)
                assert V.shape[0] == self.n_bus, f"{name} diverged on a valid call"

    def test_wrong_v_init_size(self):
        for name, pf in self._pf_methods():
            for descr, v in [("too short", self.v_init[:-1]),
                             ("too long", np.ones(self.n_bus + 3, dtype=complex)),
                             ("empty", np.zeros(0, dtype=complex))]:
                with self.subTest(f"{name}, Vinit {descr}"):
                    with self.assertRaises(RuntimeError):
                        pf(v, MAX_IT, TOL)

    def test_bad_max_iter(self):
        # max_iter=0 is legitimate (returns the pre-iteration state); only
        # negative values are rejected
        for name, pf in self._pf_methods():
            for bad_iter in (-1, -100):
                with self.subTest(f"{name}, max_iter={bad_iter}"):
                    with self.assertRaises(RuntimeError):
                        pf(1.0 * self.v_init, bad_iter, TOL)

    def test_zero_max_iter_is_valid(self):
        # 0 is the well-defined "no iterations, just the initial state" call
        for name, pf in self._pf_methods():
            with self.subTest(name):
                pf(1.0 * self.v_init, 0, TOL)  # must not raise

    def test_bad_tol(self):
        for name, pf in self._pf_methods():
            for bad_tol in (0., -1e-8, float("nan"), float("inf"), -float("inf")):
                with self.subTest(f"{name}, tol={bad_tol}"):
                    with self.assertRaises(RuntimeError):
                        pf(1.0 * self.v_init, MAX_IT, bad_tol)

    def test_wrong_types(self):
        for name, pf in self._pf_methods():
            calls = [
                ("V is a string", lambda pf=pf: pf("not a voltage vector", MAX_IT, TOL)),
                ("V is a list of strings", lambda pf=pf: pf(["a"] * self.n_bus, MAX_IT, TOL)),
                ("max_iter is a string", lambda pf=pf: pf(1.0 * self.v_init, "many", TOL)),
                ("max_iter is a float", lambda pf=pf: pf(1.0 * self.v_init, 3.5, TOL)),
                ("tol is a string", lambda pf=pf: pf(1.0 * self.v_init, MAX_IT, "tiny")),
                ("tol is complex", lambda pf=pf: pf(1.0 * self.v_init, MAX_IT, 1e-8 + 1j)),
            ]
            for descr, fun in calls:
                with self.subTest(f"{name}, {descr}"):
                    with self.assertRaises(TypeError):
                        fun()


class TestSolverComputePfInputs(unittest.TestCase):
    """Direct Solver.compute_pf / Solver.solve with corrupted inputs, for
    every AC solver class available in this build.

    A single valid input set is built from the grid (one converged ac_pf,
    then the *solver-side* accessors), then exactly one aspect is corrupted
    per case. Before the validation in BaseAlgo::check_pf_inputs, several of
    these reached raw Eigen indexing (out-of-bounds reads/writes, invisible
    in Release builds where Eigen's asserts are compiled out).
    """

    @classmethod
    def setUpClass(cls):
        grid = _make_grid()
        n_bus = grid.total_bus()
        grid.ac_pf(np.ones(n_bus, dtype=complex), MAX_IT, TOL)
        cls.Ybus = grid.get_Ybus_solver()  # scipy sparse, n x n
        cls.Sbus = 1.0 * grid.get_Sbus_solver()
        cls.V0 = np.ones(cls.Ybus.shape[0], dtype=complex)
        cls.slack_ids = np.array(grid.get_slack_ids_solver(), dtype=np.int32)
        cls.slack_weights = 1.0 * grid.get_slack_weights_solver()
        cls.pv = np.array(grid.get_pv_solver(), dtype=np.int32)
        cls.pq = np.array(grid.get_pq_solver(), dtype=np.int32)
        cls.n = cls.Ybus.shape[0]

    def _solve(self, solver_cls_name, Ybus=None, V=None, Sbus=None,
               slack_ids=None, slack_weights=None, pv=None, pq=None,
               max_iter=MAX_IT, tol=TOL):
        solver = getattr(lightsim2grid_cpp, solver_cls_name)()
        return solver.compute_pf(
            self.Ybus if Ybus is None else Ybus,
            1.0 * (self.V0 if V is None else V),  # compute_pf mutates V: pass a copy
            self.Sbus if Sbus is None else Sbus,
            self.slack_ids if slack_ids is None else slack_ids,
            self.slack_weights if slack_weights is None else slack_weights,
            self.pv if pv is None else pv,
            self.pq if pq is None else pq,
            max_iter,
            tol,
        )

    def _assert_all_solvers_raise(self, exc, descr, **corrupted):
        for solver_nm in AC_SOLVER_CLASSES:
            with self.subTest(f"{solver_nm}, {descr}"):
                with self.assertRaises(exc):
                    self._solve(solver_nm, **corrupted)

    def test_valid_inputs_converge(self):
        # sanity: the validation must not reject a legitimate system.
        # FDPF solvers cannot run standalone (their B'/B'' matrices are built
        # from the grid, not from Ybus): they must raise a clean RuntimeError
        # instead -- writing this very test found a null-pointer segfault there.
        for solver_nm in AC_SOLVER_CLASSES:
            with self.subTest(solver_nm):
                solver = getattr(lightsim2grid_cpp, solver_nm)()
                if solver_nm.startswith("FDPF"):
                    with self.assertRaises(RuntimeError) as cm:
                        self._solve(solver_nm, max_iter=1000)
                    assert "not attached to an LSGrid" in str(cm.exception)
                    continue
                conv = solver.compute_pf(self.Ybus, 1.0 * self.V0, self.Sbus,
                                         self.slack_ids, self.slack_weights,
                                         self.pv, self.pq, 1000, TOL)
                assert conv, f"{solver_nm} did not converge on a valid system"

    def test_non_square_ybus(self):
        self._assert_all_solvers_raise(
            RuntimeError, "non-square Ybus", Ybus=self.Ybus[:-1, :])

    def test_sbus_size_mismatch(self):
        self._assert_all_solvers_raise(
            RuntimeError, "Sbus too short", Sbus=self.Sbus[:-1])
        self._assert_all_solvers_raise(
            RuntimeError, "Sbus too long",
            Sbus=np.concatenate([self.Sbus, [0. + 0j]]))

    def test_v_size_mismatch(self):
        self._assert_all_solvers_raise(
            RuntimeError, "V too short", V=self.V0[:-1])
        self._assert_all_solvers_raise(
            RuntimeError, "V too long",
            V=np.ones(self.n + 2, dtype=complex))

    def test_slack_weights_size_mismatch(self):
        self._assert_all_solvers_raise(
            RuntimeError, "slack_weights too short",
            slack_weights=self.slack_weights[:-1])

    def test_out_of_bounds_ids(self):
        for arr_name in ("pv", "pq", "slack_ids"):
            valid = getattr(self, arr_name)
            for descr, bad_id in (("id == n", self.n), ("huge id", 10**6), ("negative id", -1)):
                corrupted = np.concatenate([valid, [bad_id]]).astype(np.int32)
                self._assert_all_solvers_raise(
                    IndexError, f"{arr_name} with {descr}", **{arr_name: corrupted})

    def test_empty_slack(self):
        self._assert_all_solvers_raise(
            RuntimeError, "empty slack_ids",
            slack_ids=np.zeros(0, dtype=np.int32))

    def test_overlapping_pv_pq(self):
        # the same bus declared both PV and PQ silently breaks the system of
        # equations: it must be rejected, not solved
        overlap = np.concatenate([self.pq, [self.pv[0]]]).astype(np.int32)
        self._assert_all_solvers_raise(
            RuntimeError, "bus both pv and pq", pq=overlap)
        # a bus listed twice in the same array is just as wrong
        twice = np.concatenate([self.pq, [self.pq[0]]]).astype(np.int32)
        self._assert_all_solvers_raise(
            RuntimeError, "bus twice in pq", pq=twice)

    def test_slack_also_pv(self):
        also_pv = np.concatenate([self.pv, [self.slack_ids[0]]]).astype(np.int32)
        self._assert_all_solvers_raise(
            RuntimeError, "slack bus also in pv", pv=also_pv)

    def test_bad_max_iter(self):
        # max_iter=0 is legitimate (returns the pre-iteration state); only
        # negative values are rejected
        for bad_iter in (-1, -100):
            self._assert_all_solvers_raise(
                RuntimeError, f"max_iter={bad_iter}", max_iter=bad_iter)

    def test_bad_tol(self):
        for bad_tol in (0., -1e-8, float("nan"), float("inf")):
            self._assert_all_solvers_raise(
                RuntimeError, f"tol={bad_tol}", tol=bad_tol)


if __name__ == "__main__":
    unittest.main()
