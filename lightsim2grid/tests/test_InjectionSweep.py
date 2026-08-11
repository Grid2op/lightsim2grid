# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
InjectionSweepCPP computes the very same thing as TimeSeriesCPP -- one powerflow per row of
the injection matrices, on a fixed topology -- but initializes every step with the same
voltage instead of with the result of the step before it.

The tests below pin the two properties that follow from that, and that TimeSeriesCPP does
NOT have: the results do not depend on the order the steps are given in, and they do not
depend on how many threads the sweep was split over.
"""

import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.injectionSweep import InjectionSweep, InjectionSweepCPP
from lightsim2grid.timeSerie import TimeSeriesCPP


class TestInjectionSweepCPP(unittest.TestCase):
    def setUp(self):
        env_name = "l2rpn_case14_sandbox"
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(env_name, backend=LightSimBackend(), param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.Vinit = self.env.backend.V
        data = self.env.chronics_handler.real_data.data
        # a handful of steps is plenty here, and keeps the multi-threaded tests quick
        self.nb_steps = 24
        self.prod_p = 1.0 * data.prod_p[:self.nb_steps]
        self.load_p = 1.0 * data.load_p[:self.nb_steps]
        self.load_q = 1.0 * data.load_q[:self.nb_steps]
        self.sgen_p = np.zeros((self.nb_steps, 0))

    def tearDown(self):
        self.env.close()

    def _compute(self, computer, prod_p=None, load_p=None, load_q=None):
        prod_p = self.prod_p if prod_p is None else prod_p
        load_p = self.load_p if load_p is None else load_p
        load_q = self.load_q if load_q is None else load_q
        status = computer.compute_Vs(prod_p,
                                     np.zeros((prod_p.shape[0], 0)),
                                     load_p,
                                     load_q,
                                     self.Vinit,
                                     self.env.backend.max_it,
                                     self.env.backend.tol)
        assert status == 1, (f"the computation diverged after {computer.nb_solved()} step(s)")
        return 1.0 * computer.get_voltages()

    def test_same_results_as_time_series(self):
        """the init policy changes the starting point, not the solution it converges to"""
        Vs_sweep = self._compute(InjectionSweepCPP(self.grid))
        Vs_ts = self._compute(TimeSeriesCPP(self.grid))
        assert np.max(np.abs(Vs_sweep - Vs_ts)) <= 1e-6, (
            f"max deviation {np.max(np.abs(Vs_sweep - Vs_ts))} between InjectionSweep and TimeSeries")

    def test_matches_grid2op_step_by_step(self):
        """same reference as test_Computers.py: the flows grid2op computes one step at a time"""
        computer = InjectionSweepCPP(self.grid)
        self._compute(computer)
        amps = computer.compute_flows()
        for it_num in range(self.nb_steps - 1):
            obs, *_ = self.env.step(self.env.action_space())
            assert np.max(np.abs(amps[1 + it_num] - obs.a_or * 1e-3)) <= 1e-6, f"error at it {it_num}"

    def test_order_invariance(self):
        """
        The invariant that separates InjectionSweep from TimeSeries: permuting the steps
        permutes the results, and changes nothing else. A chained TimeSeries has no reason
        to satisfy this (each step inherits its starting point from its new neighbour).
        """
        Vs = self._compute(InjectionSweepCPP(self.grid))

        rng = np.random.default_rng(0)
        perm = rng.permutation(self.nb_steps)
        Vs_perm = self._compute(InjectionSweepCPP(self.grid),
                                prod_p=self.prod_p[perm],
                                load_p=self.load_p[perm],
                                load_q=self.load_q[perm])
        assert np.array_equal(Vs_perm, Vs[perm]), (
            f"max deviation {np.max(np.abs(Vs_perm - Vs[perm]))} after permuting the steps")

    def test_thread_independence(self):
        """splitting the sweep over threads must be bit for bit identical"""
        computer1 = InjectionSweepCPP(self.grid)
        assert computer1.nb_thread == 1, "nb_thread should default to 1"
        Vs1 = self._compute(computer1)

        for nb_thread in (2, 4):
            computer = InjectionSweepCPP(self.grid)
            computer.nb_thread = nb_thread
            assert computer.nb_thread == nb_thread
            Vs = self._compute(computer)
            assert np.array_equal(Vs, Vs1), (
                f"results differ with nb_thread={nb_thread} (max deviation "
                f"{np.max(np.abs(Vs - Vs1))})")

    def test_nb_thread_clamped(self):
        computer = InjectionSweepCPP(self.grid)
        computer.nb_thread = 0
        assert computer.nb_thread == 1
        computer.nb_thread = -3
        assert computer.nb_thread == 1

    def test_time_series_rejects_nb_thread(self):
        """
        TimeSeriesCPP chains its steps, so it cannot be split over threads: the attribute
        exists (so the error can explain itself) but refuses anything but 1.
        """
        computer = TimeSeriesCPP(self.grid)
        assert computer.nb_thread == 1
        computer.nb_thread = 1  # legal, no-op
        computer.nb_thread = 0  # clamped to 1, still legal
        assert computer.nb_thread == 1
        with self.assertRaises(RuntimeError) as ctx:
            computer.nb_thread = 2
        assert "InjectionSweep" in str(ctx.exception), (
            f"the error should point at InjectionSweep, got: {ctx.exception}")
        assert computer.nb_thread == 1, "a rejected assignment must not change the value"

    def test_init_from_n_powerflow(self):
        """the seed changes, the converged results do not"""
        computer = InjectionSweepCPP(self.grid)
        assert not computer.init_from_n_powerflow, "should default to False"
        Vs_provided = self._compute(computer)

        computer_n = InjectionSweepCPP(self.grid)
        computer_n.init_from_n_powerflow = True
        assert computer_n.init_from_n_powerflow
        Vs_n = self._compute(computer_n)

        assert np.max(np.abs(Vs_n - Vs_provided)) <= 1e-6, (
            f"max deviation {np.max(np.abs(Vs_n - Vs_provided))} between the two seeds")

    def test_thread_independence_init_from_n(self):
        """the n-powerflow seed is computed once, before the threads are spawned"""
        computer1 = InjectionSweepCPP(self.grid)
        computer1.init_from_n_powerflow = True
        Vs1 = self._compute(computer1)

        computer4 = InjectionSweepCPP(self.grid)
        computer4.init_from_n_powerflow = True
        computer4.nb_thread = 4
        Vs4 = self._compute(computer4)
        assert np.array_equal(Vs4, Vs1)

    def test_dc(self):
        from lightsim2grid.algorithm import AlgorithmType
        dc_algo = AlgorithmType.DC_SparseLU
        computer = InjectionSweepCPP(self.grid)
        computer.change_algorithm(dc_algo)
        Vs_dc = self._compute(computer)

        computer_ts = TimeSeriesCPP(self.grid)
        computer_ts.change_algorithm(dc_algo)
        Vs_ts = self._compute(computer_ts)
        assert np.max(np.abs(Vs_dc - Vs_ts)) <= 1e-6

        computer_th = InjectionSweepCPP(self.grid)
        computer_th.change_algorithm(dc_algo)
        computer_th.nb_thread = 4
        assert np.array_equal(self._compute(computer_th), Vs_dc)

    def test_new_setter_api_matches_legacy_compute_vs(self):
        """the new modify_* + compute() API must give bit-for-bit the same result as the
        legacy bundled compute_Vs() call on the same object -- guards against the SFINAE
        gating silently changing which code path the legacy call takes"""
        Vs_legacy = self._compute(InjectionSweepCPP(self.grid))

        computer = InjectionSweepCPP(self.grid)
        computer.modify_gen_p(self.prod_p)
        computer.modify_sgen_p(self.sgen_p)
        computer.modify_load_p(self.load_p)
        computer.modify_load_q(self.load_q)
        computer.compute(self.Vinit, self.env.backend.max_it, self.env.backend.tol)
        assert computer.get_status() == 1
        Vs_new = 1.0 * computer.get_voltages()

        assert np.array_equal(Vs_new, Vs_legacy)

    def test_timers(self):
        computer = InjectionSweepCPP(self.grid)
        self._compute(computer)
        assert computer.thread_init_time() == 0., "no thread is built when nb_thread is 1"

        computer4 = InjectionSweepCPP(self.grid)
        computer4.nb_thread = 4
        self._compute(computer4)
        assert computer4.nb_solved() == self.nb_steps
        assert computer4.total_time() > 0.


class TestInjectionSweepGrid2op(unittest.TestCase):
    """the grid2op level wrapper, which shares everything but the c++ class with TimeSerie"""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), test=True)

    def tearDown(self):
        self.env.close()

    def test_behaviour(self):
        sweep = InjectionSweep(self.env)
        assert isinstance(sweep.computer, InjectionSweepCPP), (
            "the wrapper should drive an InjectionSweepCPP, not a TimeSeriesCPP")
        sweep.nb_thread = 2
        assert sweep.nb_thread == 2

        sweep.compute_V(scenario_id=0)
        As = sweep.compute_A()
        Ps = sweep.compute_P()

        self.env.set_id(0)
        self.env.reset()
        for it_num in range(100):
            obs, *_ = self.env.step(self.env.action_space())
            assert np.max(np.abs(As[1 + it_num] - obs.a_or)) <= 1e-3, f"error at it {it_num} for A"
            assert np.max(np.abs(Ps[1 + it_num] - obs.p_or)) <= 1e-3, f"error at it {it_num} for P"

    def test_nb_thread_validated(self):
        sweep = InjectionSweep(self.env)
        with self.assertRaises(ValueError):
            sweep.nb_thread = 1.5

    def test_new_setter_api_is_inherited_from_timeserie(self):
        """modify_gen_p / modify_sgen_p / modify_load_p / modify_load_q / compute() are
        defined once on TimeSerie and never overridden by InjectionSweep -- inherited "for
        free", driving the InjectionSweepCPP instance set via _CPP_CLASS. This pins that
        the inherited setter-based API actually works end to end through the grid2op-level
        wrapper (not just on the raw InjectionSweepCPP object, already covered by
        TestInjectionSweepCPP.test_new_setter_api_matches_legacy_compute_vs above)."""
        n_steps = 3
        prod_p = np.tile(self.env.backend.prod_p, (n_steps, 1))
        load_p = np.tile(self.env.backend.load_p, (n_steps, 1))
        load_q = np.tile(self.env.backend.load_q, (n_steps, 1))

        sweep_setter = InjectionSweep(self.env)
        sweep_setter.modify_gen_p(prod_p)
        sweep_setter.modify_load_p(load_p)
        sweep_setter.modify_load_q(load_q)
        Vs_setter = sweep_setter.compute()
        Ps_setter = sweep_setter.compute_P()
        As_setter = sweep_setter.compute_A()

        sweep_legacy = InjectionSweep(self.env)
        Vs_legacy = sweep_legacy.compute_V_from_inj(prod_p, load_p, load_q)
        Ps_legacy = sweep_legacy.compute_P()
        As_legacy = sweep_legacy.compute_A()

        assert np.array_equal(Vs_setter, Vs_legacy)
        assert np.array_equal(Ps_setter, Ps_legacy)
        assert np.array_equal(As_setter, As_legacy)


if __name__ == "__main__":
    unittest.main()
