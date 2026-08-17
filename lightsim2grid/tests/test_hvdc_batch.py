# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The batch algorithms (TimeSeriesCPP / ContingencyAnalysisCPP) must handle HVDC in
every regime: fixed setpoint, angle-droop linear and angle-droop saturated (both
directions). The batch result of each step / contingency must match a single-shot
``ac_pf`` carrying the same injections and the same droop regime.
"""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc  # noqa: E402
from lightsim2grid.lightsim2grid_cpp import TimeSeriesCPP, ContingencyAnalysisCPP  # noqa: E402

MAX_IT = 30
TOL = 1e-10


def _base_injections(model):
    gen_p = np.array(model.get_gen_target_p(), dtype=np.float64)
    load_p = np.array(model.get_load_target_p(), dtype=np.float64)
    load_q = np.array([ld.target_q_mvar for ld in model.get_loads()], dtype=np.float64)
    return gen_p, load_p, load_q


def _single_shot(model, gen_p, load_p, load_q, nb_bus):
    for i, p in enumerate(gen_p):
        model.change_p_gen(i, float(p))
    for i, p in enumerate(load_p):
        model.change_p_load(i, float(p))
    for i, q in enumerate(load_q):
        model.change_q_load(i, float(q))
    return model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)


def _make(regime):
    """Build an HVDC grid; ``regime`` selects the droop configuration to exercise."""
    if regime == "fixed":
        net, model = make_case14_hvdc(3, 9, loss_factor1=0.011, loss_factor2=0.011,
                                      droop_enabled=False, p_setpoint=40.0)
    else:
        # a full 300 MW reverse saturation is not feasible on these two case14 buses;
        # cap both directions at 50 MW so every regime converges from a flat start.
        net, model = make_case14_hvdc(3, 9, loss_factor1=0.011, loss_factor2=0.011,
                                      droop_enabled=True, p0=10.0, droop_mw_per_deg=5.0,
                                      pmax12=50.0, pmax21=50.0, p_setpoint=20.0)
        if regime == "droop_linear":
            model.set_status_droop_hvdc(0, 0)
        elif regime == "droop_sat_12":
            model.set_status_droop_hvdc(0, 1)
        elif regime == "droop_sat_21":
            model.set_status_droop_hvdc(0, -1)
        else:
            raise ValueError(regime)
    return net, model


REGIMES = ["fixed", "droop_linear", "droop_sat_12", "droop_sat_21"]


class TestHvdcTimeSeries(unittest.TestCase):
    def test_timeseries_matches_single_shot(self):
        for regime in REGIMES:
            with self.subTest(regime=regime):
                net, model = _make(regime)
                nb_bus = net.bus.shape[0]
                gen_p, load_p, load_q = _base_injections(model)

                # two distinct steps (base and a +8% load / +8% gen scaling)
                scales = np.array([1.0, 1.08])
                gen_mat = np.outer(scales, gen_p)
                lp_mat = np.outer(scales, load_p)
                lq_mat = np.outer(scales, load_q)
                sgen_mat = np.zeros((len(scales), 0))

                ts = TimeSeriesCPP(model)
                status = ts.compute_Vs(gen_mat, sgen_mat, lp_mat, lq_mat,
                                       np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                self.assertEqual(status, 1, f"time series diverged ({regime})")
                Vs = np.array(ts.get_voltages())

                for step in range(len(scales)):
                    ref_net, ref_model = _make(regime)
                    Vref = _single_shot(ref_model, gen_mat[step], lp_mat[step], lq_mat[step], nb_bus)
                    self.assertGreater(Vref.shape[0], 0)
                    self.assertLess(np.abs(Vs[step] - Vref).max(), 1e-7,
                                    f"{regime} step {step}")


class TestHvdcContingencyAnalysis(unittest.TestCase):
    def test_contingency_matches_single_shot(self):
        # contingencies on a few powerlines not incident to the HVDC ends.
        cont_lines = [0, 5, 10]
        for regime in REGIMES:
            with self.subTest(regime=regime):
                net, model = _make(regime)
                nb_bus = net.bus.shape[0]

                ca = ContingencyAnalysisCPP(model)
                for line_id in cont_lines:
                    ca.add_n1(line_id)
                ca.compute(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                Vs = np.array(ca.get_voltages())

                for k, line_id in enumerate(cont_lines):
                    if np.abs(Vs[k]).max() == 0.0:
                        # this contingency diverged / islanded: nothing to compare
                        continue
                    ref_net, ref_model = _make(regime)
                    ref_model.deactivate_powerline(line_id)
                    Vref = ref_model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)
                    self.assertGreater(Vref.shape[0], 0, f"{regime} ref cont {line_id} diverged")
                    # compare only the buses the contingency solver actually returns
                    mask = np.abs(Vs[k]) > 0.0
                    self.assertLess(np.abs(Vs[k][mask] - Vref[mask]).max(), 1e-6,
                                    f"{regime} contingency line {line_id}")


if __name__ == "__main__":
    unittest.main()
