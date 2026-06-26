# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tap-changer step r/x/g/b correction when reading a transformer from pypowsybl.

pypowsybl exposes the transformer r/x/g/b at the *neutral* tap and only folds the
tap position into ``rho`` / ``alpha``. The per-step r/x/g/b deltas (in percent)
must be applied on top, otherwise the (phase-shifting) transformer impedance --
and hence its through-flow -- is wrong. On real RTE grids this caused tens of MW
of disagreement with PowSyBl Open Load Flow (see HVDC_OLF_FINDINGS.md).

The ``four_substations`` pypowsybl test network carries a phase-shifting
transformer ``TWT`` whose phase-tap step applies a ~-28.8% r/x correction, which
makes it a compact regression for the fix.
"""

import unittest
import warnings

import numpy as np

try:
    import pypowsybl.network as pn
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import init_from_pypowsybl
    HAS_4SUB = hasattr(pn, "create_four_substations_node_breaker_network")
except ImportError:
    HAS_4SUB = False


class TestPstTapImpedance(unittest.TestCase):
    def setUp(self):
        if not HAS_4SUB:
            self.skipTest("pypowsybl (with the four_substations network) is required")

    def _ls_trafo(self, model, name):
        # the iterator yields a reused proxy: read the attributes during iteration
        for el in model.get_trafos():
            if el.name == name:
                return el.r_pu, el.x_pu
        self.fail(f"transformer {name} not found in the lightsim2grid model")

    def test_phase_tap_step_rx_correction_applied(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            n = pn.create_four_substations_node_breaker_network()
            lf.run_ac(n)

        # base (neutral-tap) per-unit impedance + the current phase-tap step delta
        npu = pn.create_four_substations_node_breaker_network()
        lf.run_ac(npu)
        npu.per_unit = True
        trpu = npu.get_2_windings_transformers(all_attributes=True)
        base_r = float(trpu.loc["TWT", "r"])
        base_x = float(trpu.loc["TWT", "x"])
        ptc = n.get_phase_tap_changers()
        tap = int(ptc.loc["TWT", "tap"])
        step = n.get_phase_tap_changer_steps().loc[("TWT", tap)]
        self.assertGreater(abs(step["x"]), 1.0, "this test needs a non-trivial r/x step correction")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False)

        r_pu, x_pu = self._ls_trafo(model, "TWT")
        exp_r = base_r * (1.0 + step["r"] / 100.0)
        exp_x = base_x * (1.0 + step["x"] / 100.0)
        # the step correction must be applied ...
        self.assertAlmostEqual(r_pu, exp_r, places=9,
                               msg=f"r_pu {r_pu} != corrected {exp_r}")
        self.assertAlmostEqual(x_pu, exp_x, places=9,
                               msg=f"x_pu {x_pu} != corrected {exp_x}")
        # ... and it must actually differ from the uncorrected (neutral) impedance
        self.assertGreater(abs(x_pu - base_x), 1e-4,
                           "the step correction was not applied (x_pu == neutral x)")

    def test_in_place_shift_change_refreshes_rx(self):
        # changing the phase shift in place (no re-import, no "tap") must refresh the
        # series impedance from the alpha -> r/x dependency, NOT leave it stale.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            n = pn.create_four_substations_node_breaker_network()
            lf.run_ac(n)
            npu = pn.create_four_substations_node_breaker_network()
            lf.run_ac(npu)
            npu.per_unit = True
            base_x = float(npu.get_2_windings_transformers(all_attributes=True).loc["TWT", "x"])
            model = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False)

        steps = n.get_phase_tap_changer_steps()
        tap_a = int(n.get_phase_tap_changers().loc["TWT", "tap"])
        # pick another tap whose step x% differs noticeably from the loaded one
        x_a = float(steps.loc[("TWT", tap_a)]["x"])
        cands = [(t, float(steps.loc[("TWT", t)]["x"]))
                 for t in steps.loc["TWT"].index if abs(float(steps.loc[("TWT", t)]["x"]) - x_a) > 5.0]
        tap_b, x_b_pct = cands[0]
        alpha_b_rad = float(np.deg2rad(steps.loc[("TWT", tap_b)]["alpha"]))

        tid = next(i for i, el in enumerate(model.get_trafos()) if el.name == "TWT")
        model.change_shift_trafo(tid, alpha_b_rad)  # in-place, like an OLF tap change
        _, x_pu = self._ls_trafo(model, "TWT")
        self.assertAlmostEqual(x_pu, base_x * (1.0 + x_b_pct / 100.0), places=8,
                               msg="series impedance not refreshed after an in-place shift change")


if __name__ == "__main__":
    unittest.main()
