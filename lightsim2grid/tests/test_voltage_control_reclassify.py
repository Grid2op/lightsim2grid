# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A generator may remotely regulate a bus that another generator already regulates
locally: the two simply form one control group, sharing the reactive duty in
proportion to their ``qmax - qmin`` ranges.

That used to be refused. ``GeneratorContainer::fillpv`` only knew how to keep a
controller's OWN bus out of the PV set; nothing stopped a local regulator from
claiming a bus that somebody else regulated remotely, and the group's bordered
voltage row was then left with no ``Vm`` column to act on. The configuration came
back as "regulates bus X which has no voltage (Vm) unknown" even though it is
perfectly well posed.

What must still be refused is the genuinely unsolvable arrangement: a remote
regulator sharing its own bus with a local regulator. Whatever reactive power it
injects there is absorbed one-for-one by the machine holding that bus' magnitude, so
it has no influence at all on its remote target.

Companion C++ suite: ``src/tests/test_voltage_control_reclassify.cpp``.
"""

import warnings
import unittest
import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower

MAX_IT = 30
TOL = 1e-11
# case14 as lightsim sees it: gens 0..3 on buses 1, 2, 5, 7 and the ext_grid as gen 4
# on bus 0 (the slack).
GEN_AT_BUS1 = 0
GEN_AT_BUS2 = 1
SLACK_BUS = 0


def _case14(extra_gen_bus=None):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        if extra_gen_bus is not None:
            pp.create_gen(net, bus=extra_gen_bus, p_mw=0.0, vm_pu=1.01,
                          min_q_mvar=-100., max_q_mvar=100., controllable=True)
        model = init_from_pandapower(net)
    return net, model


def _solve(model, nb_bus):
    return model.ac_pf(np.ones(nb_bus, dtype=np.complex128), MAX_IT, TOL)


def _gen_q_mvar(model, gen_id):
    return model.get_generators()[gen_id].res_q_mvar


class TestRemoteOntoLocallyRegulatedBus(unittest.TestCase):
    def test_remote_may_target_a_locally_regulated_bus(self):
        # gen 1 sits on bus 2 and regulates it locally; gen 0 (bus 1) aims at bus 2
        # remotely, with the same setpoint.
        net, model = _case14()
        vset = model.get_generators()[GEN_AT_BUS2].target_vm_pu
        model.change_v_gen(GEN_AT_BUS1, vset)
        model.set_gen_regulated_bus(GEN_AT_BUS1, 2)
        model.tell_solver_need_reset()

        V = _solve(model, net.bus.shape[0])
        self.assertGreater(V.shape[0], 0, "ac_pf diverged")
        self.assertLess(abs(abs(V[2]) - vset), 1e-6,
                        f"|V[2]|={abs(V[2]):.8f} != setpoint {vset:.8f}")

        # both machines carry reactive power, split by their qmax - qmin ranges
        gens = model.get_generators()
        w_remote = gens[GEN_AT_BUS1].max_q_mvar - gens[GEN_AT_BUS1].min_q_mvar
        w_local = gens[GEN_AT_BUS2].max_q_mvar - gens[GEN_AT_BUS2].min_q_mvar
        q_remote = _gen_q_mvar(model, GEN_AT_BUS1)
        q_local = _gen_q_mvar(model, GEN_AT_BUS2)
        self.assertGreater(abs(q_remote), 1e-6, "the remote controller carries no Q")
        self.assertLess(abs(q_remote / w_remote - q_local / w_local), 1e-6,
                        f"Q not shared by range: {q_remote}/{w_remote} vs {q_local}/{w_local}")

    def test_conflicting_setpoints_are_named_as_such(self):
        # same shape but disagreeing setpoints: both are in one group now, so the
        # error is about the setpoints rather than about the bus classification
        net, model = _case14()
        vset = model.get_generators()[GEN_AT_BUS2].target_vm_pu
        model.change_v_gen(GEN_AT_BUS1, vset + 0.02)
        model.set_gen_regulated_bus(GEN_AT_BUS1, 2)
        model.tell_solver_need_reset()
        with self.assertRaises(RuntimeError) as ctx:
            _solve(model, net.bus.shape[0])
        self.assertIn("conflicting voltage setpoints", str(ctx.exception))

    def test_remote_may_target_the_slack_bus(self):
        # the slack (bus 0) is pinned by the ext_grid generator, so this used to be
        # refused; that generator now joins the group and MultiSlack grants the
        # free Vm the bordered row needs
        net, model = _case14()
        slack_gen = [i for i, g in enumerate(model.get_generators()) if g.is_slack][0]
        vset = model.get_generators()[slack_gen].target_vm_pu
        model.change_v_gen(GEN_AT_BUS1, vset)
        model.set_gen_regulated_bus(GEN_AT_BUS1, SLACK_BUS)
        model.tell_solver_need_reset()

        V = _solve(model, net.bus.shape[0])
        self.assertGreater(V.shape[0], 0, "ac_pf diverged")
        self.assertLess(abs(abs(V[SLACK_BUS]) - vset), 1e-6)


class TestStillRefused(unittest.TestCase):
    def test_remote_sharing_its_bus_with_a_local_regulator(self):
        # a SECOND generator on bus 2 (which gen 1 pins locally) aiming remotely at
        # bus 9: its reactive injection is swallowed by gen 1, so it has no influence
        # on bus 9 and there is no solution. Must stay refused.
        net, model = _case14(extra_gen_bus=2)
        on_bus2 = [i for i, g in enumerate(model.get_generators()) if g.bus_id == 2]
        self.assertEqual(len(on_bus2), 2, "expected two generators on bus 2")
        model.set_gen_regulated_bus(on_bus2[-1], 9)
        model.tell_solver_need_reset()
        with self.assertRaises(RuntimeError) as ctx:
            _solve(model, net.bus.shape[0])
        self.assertIn("OWN bus has no reactive (Q) equation", str(ctx.exception))


class TestOrdinaryPVUnchanged(unittest.TestCase):
    def test_several_local_generators_on_one_bus(self):
        # the common configuration the reclassification must leave completely alone:
        # two machines on one bus, both regulating it locally, no remote controller
        # anywhere. Voltages must match plain case14 with the extra generator at 0 MW.
        net, model = _case14(extra_gen_bus=2)
        model.tell_solver_need_reset()
        V = _solve(model, net.bus.shape[0])
        self.assertGreater(V.shape[0], 0, "ac_pf diverged")

        _, ref = _case14()
        Vref = _solve(ref, net.bus.shape[0])
        self.assertGreater(Vref.shape[0], 0)
        # the added generator produces 0 MW and regulates the same bus to the same
        # value, so the operating point is unchanged
        self.assertLess(np.abs(V - Vref).max(), 1e-6)


if __name__ == "__main__":
    unittest.main()
