# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import tempfile
import unittest
import warnings

import numpy as np

from lightsim2grid.network import init_from_matpower


def _toy_mpc(gen_status_row3=1, bus4_type=1):
    """
    A small 4-bus toy matpower case, on purpose using non-contiguous / non-0-based
    bus numbers (10, 20, 30, 40) to exercise the bus remap, one transformer (branch
    10-40, TAP=1.05), a shunt at bus 40 (GS=5, BS=3) and, most importantly, TWO
    generators connected to the same (reference) bus 10 plus a 3rd, deactivated, one
    at that same bus -- the multi-generator-per-bus scenario this loader exists for.
    """
    bus = np.array([
        [10, 3, 0,  0,  0, 0, 1, 1.0, 0, 138, 1, 1.1, 0.9],
        [20, 1, 50, 20, 0, 0, 1, 1.0, 0, 138, 1, 1.1, 0.9],
        [30, 2, 30, 10, 0, 0, 1, 1.0, 0, 138, 1, 1.1, 0.9],
        [40, bus4_type, 20, 5, 5, 3, 1, 1.0, 0, 138, 1, 1.1, 0.9],
    ], dtype=float)
    gen = np.array([
        [10, 20, 0, 999, -999, 1.0, 100, 1,               999, -999],
        [10, 20, 0, 999, -999, 1.0, 100, 1,               999, -999],
        [10, 20, 0, 999, -999, 1.0, 100, gen_status_row3,  999, -999],
        [30, 40, 0, 999, -999, 1.0, 100, 1,               999, -999],
    ], dtype=float)
    branch = np.array([
        [10, 20, 0.01, 0.05, 0.02, 250, 250, 250, 0,    0, 1, -360, 360],
        [20, 30, 0.01, 0.05, 0.02, 250, 250, 250, 0,    0, 1, -360, 360],
        [30, 40, 0.01, 0.05, 0.02, 250, 250, 250, 0,    0, 1, -360, 360],
        [40, 10, 0.01, 0.08, 0.0,  250, 250, 250, 1.05, 0, 1, -360, 360],
    ], dtype=float)
    return {"baseMVA": 100.0, "bus": bus, "gen": gen, "branch": branch}


def _run_pf(model, n_bus=4):
    return model.ac_pf(np.full(n_bus, 1.0, dtype=complex), 10, 1e-8)


_TOY_M_FILE = """function mpc = case_test
mpc.version = '2';
mpc.baseMVA = 100;

%% bus data
mpc.bus = [
\t10\t3\t0\t0\t0\t0\t1\t1.0\t0\t138\t1\t1.1\t0.9;
\t20\t1\t50\t20\t0\t0\t1\t1.0\t0\t138\t1\t1.1\t0.9;
\t30\t2\t30\t10\t0\t0\t1\t1.0\t0\t138\t1\t1.1\t0.9;
\t40\t1\t20\t5\t5\t3\t1\t1.0\t0\t138\t1\t1.1\t0.9;
];

%% generator data
mpc.gen = [
\t10\t20\t0\t999\t-999\t1.0\t100\t1\t999\t-999;
\t10\t20\t0\t999\t-999\t1.0\t100\t1\t999\t-999;
\t30\t40\t0\t999\t-999\t1.0\t100\t1\t999\t-999;
];

%% branch data
mpc.branch = [
\t10\t20\t0.01\t0.05\t0.02\t250\t250\t250\t0\t0\t1\t-360\t360;
\t20\t30\t0.01\t0.05\t0.02\t250\t250\t250\t0\t0\t1\t-360\t360;
\t30\t40\t0.01\t0.05\t0.02\t250\t250\t250\t0\t0\t1\t-360\t360;
\t40\t10\t0.01\t0.08\t0.0\t250\t250\t250\t1.05\t0\t1\t-360\t360;
];
"""


class TestInitFromMatpowerDict(unittest.TestCase):
    """Tests using an already parsed matpower case (a plain dict), which needs no
    optional dependency (no matpowercaseframes, no scipy)."""

    def test_basic_conversion_converges(self):
        model = init_from_matpower(_toy_mpc())
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0, "powerflow should have converged")

    def test_load_values_no_conversion_needed(self):
        # matpower's PD/QD (bus columns) should be passed straight through, unchanged
        model = init_from_matpower(_toy_mpc())
        _run_pf(model)
        load_p, load_q, _ = model.get_loads_res()
        np.testing.assert_allclose(sorted(load_p), sorted([50., 30., 20.]))
        np.testing.assert_allclose(sorted(load_q), sorted([20., 10., 5.]))

    def test_shunt_bs_sign_is_flipped(self):
        # matpower's BS is "injected" (generator-like sign), lightsim2grid/pandapower's
        # q_mvar is "demanded" (load-like sign) -> BS must be negated, so a positive
        # BS=3 should be reported back as a (roughly, since V != exactly 1pu) negative q
        model = init_from_matpower(_toy_mpc())
        _run_pf(model)
        shunt_p, shunt_q, _ = model.get_shunts_res()
        self.assertEqual(shunt_p.shape[0], 1)
        self.assertGreater(shunt_p[0], 0.)
        self.assertLess(shunt_q[0], 0., "BS=3 (positive, injected) should map to a negative q_mvar (demanded)")

    def test_multiple_generators_same_bus_are_not_aggregated(self):
        """The whole point of this native loader: several `mpc.gen` rows sharing the
        same `GEN_BUS` (here, the reference bus) must survive as independent
        generators, distributed-slack-style, instead of being silently merged like
        pandapower's matpower converter reportedly does."""
        raw = _toy_mpc()
        model = init_from_matpower(raw)
        gen_p, gen_q, gen_v = model.get_gen_res()
        self.assertEqual(len(gen_p), raw["gen"].shape[0], "no generator should be dropped or merged")

        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)
        gen_p, gen_q, gen_v = model.get_gen_res()
        # the two in-service slack generators (rows 0 and 1) share the slack bus and
        # an equal weight -> they should end up with (almost) identical dispatch
        self.assertAlmostEqual(gen_p[0], gen_p[1], places=6)
        self.assertAlmostEqual(gen_q[0], gen_q[1], places=6)
        # the PV generator (row 3) is not slack: its P should stay at its PG setpoint
        self.assertAlmostEqual(gen_p[3], 40., places=4)

    def test_deactivated_generator_status(self):
        raw = _toy_mpc(gen_status_row3=0)  # 3rd row (also at slack bus) is now OFF
        model = init_from_matpower(raw)
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)
        gen_p, gen_q, gen_v = model.get_gen_res()
        self.assertEqual(len(gen_p), raw["gen"].shape[0])
        self.assertEqual(gen_p[2], 0.)
        self.assertEqual(gen_q[2], 0.)

    def test_isolated_bus_is_defensively_deactivated(self):
        """A raw matpower file is not guaranteed (unlike a pandapower net) to keep
        every load/shunt/branch consistent with `BUS_TYPE == 4` (isolated) buses:
        make sure lightsim2grid does not raise and just deactivates what's connected."""
        raw = _toy_mpc(bus4_type=4)  # bus 40 becomes isolated, but still has a load/shunt/branch wired to it
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model = init_from_matpower(raw)
            self.assertTrue(any("isolated" in str(x.message) for x in w))
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)
        load_p, load_q, _ = model.get_loads_res()
        # the load at the now-isolated bus 40 must have been deactivated (reported as 0)
        self.assertIn(0., load_p)

    def test_no_in_service_generator_at_ref_bus_raises(self):
        raw = _toy_mpc()
        raw["gen"][:3, 7] = 0  # deactivate all 3 generators connected to the ref bus (bus 10)
        with self.assertRaises(RuntimeError):
            init_from_matpower(raw)


class TestInitFromMatpowerNBusbarPerSub(unittest.TestCase):
    """There is always exactly one substation per matpower bus (not configurable), but
    `n_busbar_per_sub` (number of busbar sections lightsim2grid allocates per substation,
    for later grid2op-like topology actions) should accept any value >= 1."""

    def test_default_is_one_busbar_per_sub(self):
        model = init_from_matpower(_toy_mpc())
        np.testing.assert_array_equal(model.get_bus_status(), [True, True, True, True])
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)

    def test_explicit_one_busbar_per_sub_same_as_default(self):
        model = init_from_matpower(_toy_mpc(), n_busbar_per_sub=1)
        np.testing.assert_array_equal(model.get_bus_status(), [True, True, True, True])
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)

    def test_several_busbar_per_sub(self):
        model = init_from_matpower(_toy_mpc(), n_busbar_per_sub=3)
        status = model.get_bus_status()
        self.assertEqual(len(status), 4 * 3)
        # only the first busbar section of each substation (the one matpower actually
        # describes) should be active; the 2 extra ones per substation are empty
        np.testing.assert_array_equal(status[:4], [True, True, True, True])
        np.testing.assert_array_equal(status[4:], [False] * 8)
        Vfinal = _run_pf(model, n_bus=12)
        self.assertGreater(Vfinal.shape[0], 0)

    def test_n_busbar_per_sub_must_be_positive(self):
        with self.assertRaises(RuntimeError):
            init_from_matpower(_toy_mpc(), n_busbar_per_sub=0)


class TestInitFromMatpowerDCLine(unittest.TestCase):
    """`mpc.dcline` is optional (most matpower cases have no HVDC line at all) and,
    when present, should be converted with the same loss / sign conventions as
    lightsim2grid's pandapower loader (`model.init_dclines`)."""

    @staticmethod
    def _dcline_row(status=1):
        # F_BUS, T_BUS, BR_STATUS, PF, PT, QF, QT, VF, VT, PMIN, PMAX,
        # QMINF, QMAXF, QMINT, QMAXT, LOSS0, LOSS1
        return [20, 40, status, 10., 9.5, 0., 0., 1.0, 1.0, -100, 100,
                -10, 10, -10, 10, 0.5, 0.05]

    def test_no_dcline_key_at_all(self):
        raw = _toy_mpc()
        self.assertNotIn("dcline", raw)
        model = init_from_matpower(raw)
        self.assertEqual(len(model.get_dclines()), 0)

    def test_dcline_target_p1_is_negated_pf(self):
        raw = _toy_mpc()
        raw["dcline"] = np.array([self._dcline_row()])
        model = init_from_matpower(raw)
        dclines = model.get_dclines()
        self.assertEqual(len(dclines), 1)
        # matpower's PF (10 MW, "from" -> "to") maps to lightsim2grid's "power received
        # at side 1" convention, so it must be negated (same as the pandapower loader)
        self.assertAlmostEqual(dclines[0].target_p1_mw, -10.)

    def test_dcline_losses_match_matpower_formula(self):
        raw = _toy_mpc()
        raw["dcline"] = np.array([self._dcline_row()])
        model = init_from_matpower(raw)
        _run_pf(model)
        dclines = model.get_dclines()
        # PT = PF - (LOSS0 + LOSS1 * PF) = 10 - (0.5 + 0.05 * 10) = 9.0
        self.assertAlmostEqual(dclines[0].res_p1_mw, -10., places=4)
        self.assertAlmostEqual(dclines[0].res_p2_mw, 9.0, places=4)

    def test_deactivated_dcline(self):
        raw = _toy_mpc()
        raw["dcline"] = np.array([self._dcline_row(status=0)])
        model = init_from_matpower(raw)
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)


class TestInitFromMatpowerFile(unittest.TestCase):
    """Tests for the file-path dispatch (".m" via matpowercaseframes, ".mat" via scipy)."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_m_file(self):
        try:
            import matpowercaseframes  # noqa
        except ImportError:
            self.skipTest("matpowercaseframes is not installed")
        path = os.path.join(self.tmpdir.name, "case_test.m")
        with open(path, "w") as f:
            f.write(_TOY_M_FILE)
        model = init_from_matpower(path)
        gen_p, _, _ = model.get_gen_res()
        self.assertEqual(len(gen_p), 3)
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)

    def test_mat_file(self):
        try:
            import scipy.io
        except ImportError:
            self.skipTest("scipy is not installed")
        raw = _toy_mpc()
        path = os.path.join(self.tmpdir.name, "case_test.mat")
        scipy.io.savemat(path, {"mpc": raw})
        model = init_from_matpower(path)
        gen_p, _, _ = model.get_gen_res()
        self.assertEqual(len(gen_p), raw["gen"].shape[0])
        Vfinal = _run_pf(model)
        self.assertGreater(Vfinal.shape[0], 0)

    def test_mat_file_single_row_tables(self):
        # scipy.io.loadmat's squeeze_me=True collapses a (1, ncols) matrix down to a
        # flat 1d array: make sure a case with a single branch/generator still works
        try:
            import scipy.io
        except ImportError:
            self.skipTest("scipy is not installed")
        raw = _toy_mpc()
        raw["bus"] = raw["bus"][:2]
        raw["gen"] = raw["gen"][:1]
        raw["branch"] = raw["branch"][:1]
        path = os.path.join(self.tmpdir.name, "case_test_1row.mat")
        scipy.io.savemat(path, {"mpc": raw})
        model = init_from_matpower(path)
        gen_p, _, _ = model.get_gen_res()
        self.assertEqual(len(gen_p), 1)

    def test_unsupported_extension_raises(self):
        path = os.path.join(self.tmpdir.name, "case_test.json")
        with open(path, "w") as f:
            f.write("{}")
        with self.assertRaises(RuntimeError):
            init_from_matpower(path)


if __name__ == "__main__":
    unittest.main()
