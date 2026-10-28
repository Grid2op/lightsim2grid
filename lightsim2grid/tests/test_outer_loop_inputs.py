# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The data the outer loops read (tap changers, shunt sections, capability curves, reactive
ranges, the loops themselves) held to docs/security.rst: whatever reaches the grid -- through a
setter, a pickle or a binary file, or a call made after the grid changed -- raises a clean
exception, it never reads or writes outside a table.

The "poisoned state" cases go through ``__setstate__`` exactly as ``pickle.loads`` does; the
binary format restores the same state through the same ``LSGrid::set_state``. The two marked
"segfaulted" crashed the interpreter before their fix."""

import pickle
import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid, ScenarioSweepCPP, RegulationMode
from lightsim2grid.algorithm import ShuntVoltageControl, TransformerVoltageControl, DistributedSlack
from lightsim2grid.tests._exotic_elements_case_ls import build_exotic_elements_case_grid

INT_MAX = 2 ** 31 - 1
INT_MIN = -2 ** 31

# LSGrid::StateRes (LSGrid.hpp), inside the pickled (major, medium, minor, state) tuple
STATE = 3
SHUNT_ID = 8
TRAFO_ID = 9
OUTER_LOOP_NAMES_ID = 24
OUTER_LOOP_INT_PARAMS_ID = 25
GRID_SETTINGS_ID = 29
GEN_REACTIVE_RANGE_MAX_ID = 30
GEN_RAW_TARGET_Q_ID = 33
GEN_CURVE_GEN_ID_ID = 34
# TrafoContainer::StateResIdx, then TapChangers::StateResIdx
RATIO_TAPS = 10
LOW_TAP, POSITION, REGULATED = 0, 1, 12
# ShuntContainer::StateResIdx
SECTION_COUNT, REGULATED_BUS, MIN_SECTION_COUNT = 3, 7, 8


def _rebuild(state):
    """``pickle.loads``' own path: ``__setstate__`` on an object never initialised"""
    obj = LSGrid.__new__(LSGrid)
    obj.__setstate__(state)
    return obj


def _set_at(state, path, value):
    if not path:
        return value
    as_list = list(state)
    as_list[path[0]] = _set_at(as_list[path[0]], path[1:], value)
    return tuple(as_list)


def _get_at(state, path):
    for i in path:
        state = state[i]
    return state


def _flat(grid):
    return np.ones(grid.total_bus(), dtype=complex)


def _trafo_grid(n_steps=2001):
    """bus 0 (138 kV, the slack generator), a transformer 0 -> 1 (20 kV) whose ratio tap changer
    (rho 0.9 .. 1.1, at 1.0) regulates bus 1, a line 1 -> 2 and a load at 2 -- the grid of
    src/tests/test_transformer_voltage_control_loop.cpp"""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.)
    grid.init_bus(3, 1, np.array([138., 20., 20.]), 0, 0)
    grid.init_powerlines(np.full(1, .02), np.full(1, .1), np.zeros(1, complex), np.array([1]), np.array([2]))
    grid.init_trafo(np.array([.005]), np.array([.08]), np.zeros(1, complex), np.ones(1), np.zeros(1),
                    [False], np.array([0]), np.array([1]), False)
    grid.init_loads(np.array([40.]), np.array([20.]), np.array([2]))
    grid.init_generators(np.zeros(1), np.ones(1), np.full(1, -1000.), np.full(1, 1000.), np.array([0]))
    grid.add_gen_slackbus(0, 1.)
    rho = [0.9 + 0.2 * k / (n_steps - 1) for k in range(n_steps)]
    zero = [0.] * n_steps
    grid.set_trafo_ratio_tap_changer(0, 0, (n_steps - 1) // 2, rho, zero, zero, zero, zero)
    grid.set_trafo_ratio_tap_regulation(0, True, 1.0, 0.002, 1)
    grid.change_algorithm("NROuter_SparseLU")
    grid.clear_outer_loops()
    grid.add_outer_loop(TransformerVoltageControl())
    return grid


def _shunt_grid(n_sections=40):
    """buses 0-1-2 in a line, the slack generator at 0, a load and a shunt of `n_sections`
    sections of 5 MVar at 2, regulating bus 2 -- the grid of
    src/tests/test_shunt_voltage_control_loop.cpp"""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.)
    grid.init_bus(3, 1, np.full(3, 138.), 0, 0)
    grid.init_powerlines(np.full(2, .01), np.full(2, .1), np.zeros(2, complex), np.array([0, 1]), np.array([1, 2]))
    grid.init_loads(np.array([60.]), np.array([40.]), np.array([2]))
    grid.init_shunt(np.zeros(1), np.zeros(1), np.array([2]))
    grid.init_generators(np.zeros(1), np.ones(1), np.full(1, -1000.), np.full(1, 1000.), np.array([0]))
    grid.add_gen_slackbus(0, 1.)
    grid.set_shunt_sections(0, 0, [0.] * n_sections, [-5. * k for k in range(1, n_sections + 1)])
    grid.set_shunt_section_regulation(0, True, 0.99, 0., 2)
    grid.change_algorithm("NROuter_SparseLU")
    grid.clear_outer_loops()
    grid.add_outer_loop(ShuntVoltageControl())
    return grid


class TestSettersRejectBadInput(unittest.TestCase):
    """Every setter the outer loops added checks its ids and sizes before writing anything."""

    def setUp(self):
        self.grid = build_exotic_elements_case_grid()
        self.n_gen = len(self.grid.get_generators())
        self.n_svc = len(self.grid.get_svcs())
        self.n_hvdc = len(self.grid.get_dclines())
        self.n_trafo = len(self.grid.get_trafos())
        self.n_shunt = len(self.grid.get_shunts())
        self.n_bus = self.grid.total_bus()

    def test_tap_changer_ids_and_tables(self):
        one = [1.]
        z = [0.]
        for bad in (-1, self.n_trafo):
            with self.assertRaises(IndexError):
                self.grid.set_trafo_ratio_tap_changer(bad, 0, 0, one, z, z, z, z)
            with self.assertRaises(IndexError):
                self.grid.change_trafo_ratio_tap(bad, 0)
            with self.assertRaises(IndexError):
                self.grid.closest_trafo_ratio_tap(bad, 1.)
        with self.assertRaises(RuntimeError):  # no step at all
            self.grid.set_trafo_ratio_tap_changer(0, 0, 0, [], [], [], [], [])
        with self.assertRaises(RuntimeError):  # columns of different lengths
            self.grid.set_trafo_ratio_tap_changer(0, 0, 0, [1., 1.], z, z, z, z)
        with self.assertRaises(RuntimeError):  # the position outside the table
            self.grid.set_trafo_ratio_tap_changer(0, 0, 1, one, z, z, z, z)
        with self.assertRaises(RuntimeError):  # alpha of another length (phase)
            self.grid.set_trafo_phase_tap_changer(0, 0, 0, one, [0., 1.], z, z, z, z)

    def test_tap_table_past_the_largest_position(self):
        # positions INT_MAX and INT_MAX + 1: the last one is not an int
        with self.assertRaises(RuntimeError):
            self.grid.set_trafo_ratio_tap_changer(0, INT_MAX, INT_MAX, [1., 1.], [0.] * 2, [0.] * 2,
                                                  [0.] * 2, [0.] * 2)
        # INT_MAX itself: a loop over the table would step past it
        with self.assertRaises(RuntimeError):
            self.grid.set_trafo_ratio_tap_changer(0, INT_MAX, INT_MAX, [1.], [0.], [0.], [0.], [0.])
        # up to INT_MAX - 1 is fine, and so is the closest tap's loop over it
        self.grid.set_trafo_ratio_tap_changer(0, INT_MAX - 2, INT_MAX - 1, [1., 1.05], [0.] * 2, [0.] * 2,
                                              [0.] * 2, [0.] * 2)
        self.assertEqual(self.grid.get_trafos()[0].ratio_tap_position, INT_MAX - 1)
        self.grid.closest_trafo_ratio_tap(0, 1.)

    def test_a_rejected_table_leaves_the_changer_as_it_was(self):
        z3 = [0.] * 3
        self.grid.set_trafo_ratio_tap_changer(0, 0, 1, [0.95, 1., 1.05], z3, z3, z3, z3)
        before = self.grid.get_trafos()[0].ratio
        with self.assertRaises(RuntimeError):
            self.grid.set_trafo_ratio_tap_changer(0, 0, 7, [1.] * 3, z3, z3, z3, z3)
        info = self.grid.get_trafos()[0]
        self.assertEqual(info.ratio_tap_position, 1)
        self.assertEqual(info.ratio, before)
        with self.assertRaises(RuntimeError):
            self.grid.change_trafo_ratio_tap(0, 3)
        self.assertEqual(self.grid.get_trafos()[0].ratio_tap_position, 1)

    def test_tap_regulation(self):
        for bus in (-2, self.n_bus):
            with self.assertRaises(IndexError):
                self.grid.set_trafo_ratio_tap_regulation(0, True, 1., 0., bus)
        with self.assertRaises(IndexError):
            self.grid.set_trafo_ratio_tap_regulation(self.n_trafo, True, 1., 0., 0)
        for side in (0, 3, -1):
            with self.assertRaises(RuntimeError):
                self.grid.set_trafo_phase_tap_regulation(0, RegulationMode.ACTIVE_POWER, True, 1., 0., side)
        with self.assertRaises(RuntimeError):  # a phase tap changer does not regulate a voltage
            self.grid.set_trafo_phase_tap_regulation(0, RegulationMode.VOLTAGE, True, 1., 0., 1)
        self.grid.set_trafo_ratio_tap_regulation(0, True, 1., 0., -1)  # -1: none

    def test_shunt_sections(self):
        for bad in (-1, self.n_shunt):
            with self.assertRaises(IndexError):
                self.grid.set_shunt_sections(bad, 0, [0.], [1.])
            with self.assertRaises(IndexError):
                self.grid.change_shunt_section_count(bad, 0)
            with self.assertRaises(IndexError):
                self.grid.set_shunt_section_regulation(bad, True, 1., 0., 0)
        with self.assertRaises(RuntimeError):
            self.grid.set_shunt_sections(0, 0, [0., 0.], [1.])
        for count in (-1, 3):
            with self.assertRaises(RuntimeError):
                self.grid.set_shunt_sections(0, count, [0., 0.], [1., 2.])
        self.grid.set_shunt_sections(0, 1, [0., 0.], [1., 2.])
        for count in (-1, 3):
            with self.assertRaises(RuntimeError):
                self.grid.change_shunt_section_count(0, count)
        self.assertEqual(self.grid.get_shunts()[0].section_count, 1)
        for bus in (-2, self.n_bus):
            with self.assertRaises(IndexError):
                self.grid.set_shunt_section_regulation(0, True, 1., 0., bus)

    def test_svc_b0(self):
        with self.assertRaises(RuntimeError):
            self.grid.set_svc_b0(np.zeros(self.n_svc + 1))
        bad = np.zeros(self.n_svc)
        bad[0] = np.nan
        with self.assertRaises(RuntimeError):
            self.grid.set_svc_b0(bad)
        self.grid.set_svc_b0(np.zeros(self.n_svc))

    def test_capability_curves(self):
        p = np.array([0., 100.])
        lo = np.array([-10., -20.])
        hi = np.array([10., 20.])
        for gen in (-1, self.n_gen):
            with self.assertRaises(RuntimeError):
                self.grid.set_gen_capability_curves(np.array([gen, gen], dtype=np.int32), p, lo, hi)
        with self.assertRaises(RuntimeError):  # vectors of different sizes
            self.grid.set_gen_capability_curves(np.array([0], dtype=np.int32), p, lo, hi)
        with self.assertRaises(RuntimeError):
            self.grid.set_gen_capability_curves(np.array([0, 0], dtype=np.int32), np.array([0., np.inf]), lo, hi)
        self.grid.set_gen_capability_curves(np.array([0, 0], dtype=np.int32), p, lo, hi)
        self.assertEqual(self.grid.gen_limits_at(0, 50.), (-15., 15.))
        # an id outside the grid, or a generator without a curve: no limits, no read
        for gen in (-1, 1, self.n_gen, 10 ** 6):
            self.assertEqual(self.grid.gen_limits_at(gen, 50.), (None, None))

    def test_reactive_ranges_and_raw_target_q(self):
        with self.assertRaises(RuntimeError):
            self.grid.set_gen_reactive_range_max(np.ones(self.n_gen + 1))
        with self.assertRaises(RuntimeError):
            self.grid.set_gen_raw_target_q(np.ones(self.n_gen - 1))
        with self.assertRaises(RuntimeError):
            self.grid.set_station_reactive_range_max(np.ones(self.n_hvdc), np.ones(self.n_hvdc + 1))
        with self.assertRaises(RuntimeError):
            self.grid.set_station_reactive_range_max(np.ones(self.n_hvdc + 1), np.ones(self.n_hvdc + 1))
        # "none" is accepted
        self.grid.set_gen_reactive_range_max(np.zeros(0))
        self.grid.set_gen_raw_target_q(np.zeros(0))
        self.grid.set_station_reactive_range_max(np.zeros(0), np.zeros(0))

    def test_outer_loops(self):
        with self.assertRaises(RuntimeError):
            self.grid.add_outer_loop(None)
        self.grid.clear_outer_loops()
        self.grid.add_outer_loop(DistributedSlack())
        with self.assertRaises(RuntimeError):  # one loop of a kind
            self.grid.add_outer_loop(DistributedSlack())

    def test_contingency_masks(self):
        grid = _shunt_grid()
        grid.change_algorithm("NR_SparseLU")
        grid.ac_pf(_flat(grid), 30, 1e-8)
        batch = ScenarioSweepCPP(grid)
        with self.assertRaises(RuntimeError):
            batch.set_contingency_shunts(np.zeros((2, 2), dtype=bool))
        with self.assertRaises(RuntimeError):
            batch.set_contingency_loads(np.zeros((2, 3), dtype=bool))


class TestPoisonedOuterState(unittest.TestCase):
    """A pickle (or binary file) whose new fields disagree with the grid raises on load."""

    def setUp(self):
        grid = build_exotic_elements_case_grid()
        z3 = [0.] * 3
        grid.set_trafo_ratio_tap_changer(0, 0, 1, [0.95, 1., 1.05], z3, z3, z3, z3)
        grid.set_trafo_ratio_tap_regulation(0, True, 1., 0.01, 1)
        grid.set_shunt_sections(0, 1, [0., 0.], [10., 20.])
        grid.set_shunt_section_regulation(0, True, 1., 0.01, 0)
        grid.set_gen_capability_curves(np.array([0, 0], dtype=np.int32), np.array([0., 100.]),
                                       np.array([-10., -20.]), np.array([10., 20.]))
        n_gen = len(grid.get_generators())
        grid.set_gen_reactive_range_max(np.full(n_gen, 50.))
        grid.set_gen_raw_target_q(np.full(n_gen, 1.))
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack())
        self.grid = grid
        self.state = grid.__getstate__()
        self.n_bus = grid.total_bus()

    def _poison(self, path, first_value):
        """the state with the first entry of the list at `path` replaced"""
        path = (STATE,) + tuple(path)
        values = list(_get_at(self.state, path))
        values[0] = first_value
        return _set_at(self.state, path, values)

    def _replace(self, path, value):
        return _set_at(self.state, (STATE,) + tuple(path), value)

    def test_valid_state_round_trips(self):
        grid = pickle.loads(pickle.dumps(self.grid))
        self.assertIsNone(grid.check_grid())
        trafo = grid.get_trafos()[0]
        self.assertEqual(trafo.ratio_tap_position, 1)
        self.assertEqual(trafo.ratio_regulated, 1)
        shunt = grid.get_shunts()[0]
        self.assertEqual(shunt.section_count, 1)
        self.assertEqual(shunt.regulated_bus, 0)
        self.assertEqual(grid.gen_limits_at(0, 50.), (-15., 15.))
        self.assertEqual([loop.name() for loop in grid.get_outer_loops()], ["DistributedSlack"])

    def test_tap_table_overflowing_an_int_is_rejected(self):
        # segfaulted: `position - low_tap` overflowed into a negative number, accepted, then
        # read as an index billions of entries past the table
        state = self._poison((TRAFO_ID, RATIO_TAPS, LOW_TAP), INT_MIN)
        with self.assertRaises(RuntimeError):
            _rebuild(state)
        state = self._poison((TRAFO_ID, RATIO_TAPS, LOW_TAP), INT_MAX)
        state = _set_at(state, (STATE, TRAFO_ID, RATIO_TAPS, POSITION),
                        [INT_MAX] + list(_get_at(state, (STATE, TRAFO_ID, RATIO_TAPS, POSITION)))[1:])
        with self.assertRaises(RuntimeError):  # positions INT_MAX .. INT_MAX + 2
            _rebuild(state)

    def test_tap_position_outside_its_table_is_rejected(self):
        for pos in (-1, 3, INT_MAX, INT_MIN):
            with self.assertRaises(RuntimeError):
                _rebuild(self._poison((TRAFO_ID, RATIO_TAPS, POSITION), pos))

    def test_tap_regulated_bus_out_of_range_is_rejected(self):
        for bus in (-2, self.n_bus, 10 ** 6):
            with self.assertRaises(IndexError):
                _rebuild(self._poison((TRAFO_ID, RATIO_TAPS, REGULATED), bus))

    def test_shunt_min_section_count_outside_its_table_is_rejected(self):
        # segfaulted: the ShuntVoltageControl loop walks the table from that count
        for count in (-1, -10 ** 8, 3, 10 ** 8):
            with self.assertRaises(RuntimeError):
                _rebuild(self._poison((SHUNT_ID, MIN_SECTION_COUNT), count))

    def test_shunt_section_count_outside_its_table_is_rejected(self):
        for count in (-1, 3):
            with self.assertRaises(RuntimeError):
                _rebuild(self._poison((SHUNT_ID, SECTION_COUNT), count))

    def test_shunt_regulated_bus_out_of_range_is_rejected(self):
        for bus in (-2, self.n_bus, 10 ** 6):
            with self.assertRaises(IndexError):
                _rebuild(self._poison((SHUNT_ID, REGULATED_BUS), bus))

    def test_unknown_outer_loop_is_rejected(self):
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((OUTER_LOOP_NAMES_ID,), "NotALoop"))
        with self.assertRaises(RuntimeError):  # names and parameters no longer parallel
            _rebuild(self._replace((OUTER_LOOP_INT_PARAMS_ID,), []))

    def test_grid_settings_of_the_wrong_length_are_rejected(self):
        settings = list(_get_at(self.state, (STATE, GRID_SETTINGS_ID)))
        for bad in (settings[:-1], settings + [0], []):
            with self.assertRaises(RuntimeError):
                _rebuild(self._replace((GRID_SETTINGS_ID,), bad))

    def test_per_generator_vectors_of_the_wrong_size_are_rejected(self):
        for idx in (GEN_REACTIVE_RANGE_MAX_ID, GEN_RAW_TARGET_Q_ID):
            values = list(_get_at(self.state, (STATE, idx)))
            with self.assertRaises(RuntimeError):
                _rebuild(self._replace((idx,), values + [1.]))
            with self.assertRaises(RuntimeError):
                _rebuild(self._replace((idx,), values[:-1]))

    def test_capability_curve_of_an_unknown_generator_is_rejected(self):
        for gen in (-1, len(self.grid.get_generators()), 10 ** 6):
            with self.assertRaises(RuntimeError):
                _rebuild(self._poison((GEN_CURVE_GEN_ID_ID,), gen))
        # the four vectors no longer parallel
        with self.assertRaises(RuntimeError):
            _rebuild(self._replace((GEN_CURVE_GEN_ID_ID,), [0]))


class TestStaleTapPositions(unittest.TestCase):
    """A solve's outer loops leave a transformer at a tap; the table is then replaced by a shorter
    one. What the checks after the solve read at that tap is no longer in the table: they raise
    rather than read past its end (a silent out-of-bounds read before)."""

    def setUp(self):
        self.grid = _trafo_grid()
        self.assertGreater(self.grid.ac_pf(_flat(self.grid), 30, 1e-12).size, 0)
        self.moved = self.grid.get_trafos()[0].res_ratio_tap_position
        self.assertGreater(self.moved, 0)  # the loop moved it, far from position 0

    def test_checks_after_a_shorter_table_raise(self):
        self.grid.set_trafo_ratio_tap_changer(0, 0, 0, [1.], [0.], [0.], [0.], [0.])
        with self.assertRaises(IndexError):
            self.grid.get_violations()
        with self.assertRaises(IndexError):
            self.grid.get_physical_violations()

    def test_a_new_solve_clears_it(self):
        self.grid.set_trafo_ratio_tap_changer(0, 0, 0, [1.], [0.], [0.], [0.], [0.])
        self.assertGreater(self.grid.ac_pf(_flat(self.grid), 30, 1e-12).size, 0)
        self.grid.get_violations()
        self.grid.get_physical_violations()
        self.assertEqual(self.grid.get_trafos()[0].res_ratio_tap_position, 0)

    def test_an_unchanged_table_still_works(self):
        self.grid.get_violations()
        self.grid.get_physical_violations()


class TestPoisonedShuntSectionsSolve(unittest.TestCase):
    """The poisoned minimum section count would have been read by the ShuntVoltageControl
    loop at the next solve: the same grid, untouched, solves."""

    def test_valid_state_solves_with_the_loop(self):
        grid = pickle.loads(pickle.dumps(_shunt_grid()))
        grid.change_algorithm("NROuter_SparseLU")
        self.assertGreater(grid.ac_pf(_flat(grid), 30, 1e-10).size, 0)
        self.assertGreater(grid.get_shunts()[0].res_section_count, 0)

    def test_poisoned_state_is_rejected_before_any_solve(self):
        state = _shunt_grid().__getstate__()
        path = (STATE, SHUNT_ID, MIN_SECTION_COUNT)
        with self.assertRaises(RuntimeError):
            _rebuild(_set_at(state, path, [-10 ** 8]))


if __name__ == "__main__":
    unittest.main()
