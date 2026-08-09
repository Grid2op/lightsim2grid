# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``TimeSeriesCPP.compute_Vs`` does not use the injection vector the gridmodel
builds: it rebuilds its own, per step, from the four matrices the caller passes --
generator and static-generator ACTIVE power, and load active + reactive power.

``LSGrid::fillSbus_me`` stamps rather more than that: storage units,
REACTIVE_POWER-mode SVCs, the reactive setpoint of a generator whose voltage
regulation is off, and a static generator's reactive setpoint. None of those is a
per-step input, so all are constant across the steps and must be taken from the grid.
Only the HVDC term ever was, so the rest was silently dropped and the powerflow
converged on an injection the caller never asked for.

The invariant checked here: fed the grid's OWN injections, a one-step time series must
reproduce ``ac_pf`` exactly.

Companion C++ suite: ``src/tests/test_timeseries_sbus.cpp``.
"""

import warnings
import unittest
import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import TimeSeriesCPP  # noqa: E402

MAX_IT = 30
TOL = 1e-11


def _model(add_extras):
    """case14, optionally with the element kinds whose injection used to be lost."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        if add_extras:
            # a static generator carrying reactive power
            pp.create_sgen(net, bus=4, p_mw=10.0, q_mvar=25.0)
            # a storage unit (no per-step matrix at all)
            pp.create_storage(net, bus=6, p_mw=15.0, max_e_mwh=100.0, q_mvar=20.0)
            # a generator that does NOT regulate voltage, so its q_mvar is an input
            pp.create_gen(net, bus=9, p_mw=20.0, vm_pu=1.0, q_mvar=30.0,
                          min_q_mvar=-100.0, max_q_mvar=100.0, controllable=True)
        model = init_from_pandapower(net)
    if add_extras:
        # pandapower's converter marks every gen as voltage-regulating; turn the last
        # one off so its reactive setpoint is what Sbus carries
        gens = model.get_generators()
        model.change_v_gen(len(gens) - 1, 1.0)
        model.turnedoff_no_pv()
    model.tell_solver_need_reset()
    return net, model


def _own_injections(model):
    gen_p = np.array(model.get_gen_target_p(), dtype=np.float64)[None, :]
    sgen_p = np.array(model.get_sgen_target_p(), dtype=np.float64)[None, :]
    load_p = np.array(model.get_load_target_p(), dtype=np.float64)[None, :]
    load_q = np.array([ld.target_q_mvar for ld in model.get_loads()],
                      dtype=np.float64)[None, :]
    return gen_p, sgen_p, load_p, load_q


class TestTimeSeriesSbusCompleteness(unittest.TestCase):
    def _check(self, add_extras):
        net, model = _model(add_extras)
        nb_bus = net.bus.shape[0]
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged")

        _, ts_model = _model(add_extras)
        ts = TimeSeriesCPP(ts_model)
        gen_p, sgen_p, load_p, load_q = _own_injections(ts_model)
        status = ts.compute_Vs(gen_p, sgen_p, load_p, load_q,
                               np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
        self.assertEqual(status, 1, "time series diverged")
        Vts = np.array(ts.get_voltages())[0]
        self.assertLess(np.abs(Vts - V).max(), 1e-9,
                        f"TimeSeries disagrees with ac_pf by {np.abs(Vts - V).max():.3e}")

    def test_plain_case14(self):
        # the control: case14 has no sgen, no storage and every generator regulates
        # voltage, so all of its injection IS covered by the four per-step matrices
        self._check(add_extras=False)

    def test_case14_with_sgen_storage_and_a_q_mode_generator(self):
        # the same grid plus the element kinds whose contribution used to vanish
        self._check(add_extras=True)

    def test_extras_actually_change_the_operating_point(self):
        # guard the guard: if the added elements did not move the voltages, the
        # comparison above would pass no matter what TimeSeries did with them
        net, plain = _model(False)
        _, extra = _model(True)
        nb_bus = net.bus.shape[0]
        Vp = plain.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
        Ve = extra.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
        self.assertGreater(np.abs(Ve - Vp).max(), 1e-3,
                           "the added elements do not affect the solution")


if __name__ == "__main__":
    unittest.main()
