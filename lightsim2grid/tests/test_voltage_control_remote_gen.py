# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Remote voltage control of generators through the bordered VoltageControl
NRSystem extension (no pypowsybl needed: configured directly on a case14 model)."""

import warnings
import unittest
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower


def _model():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case14()
        model = init_from_pandapower(net)
    return net, model


class TestRemoteGenVoltageControl(unittest.TestCase):
    def setUp(self):
        # gen 3 sits on bus 7 (a PV bus locally); bus 9 is a PQ load bus that no
        # generator regulates: a good remote-control probe (7 and 9 are distinct).
        self.gen_id = 3
        self.ctrl_bus = 7
        self.reg_bus = 9
        self.tol = 1e-6

    def _solve(self, model, nb_bus, tol=1e-11):
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, tol)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged with remote voltage control")
        return V

    def test_remote_setpoint_is_met(self):
        net, model = self._model()
        vset = model.get_generators()[self.gen_id].target_vm_pu
        model.set_gen_regulated_bus(self.gen_id, self.reg_bus)
        model.tell_solver_need_reset()
        V = self._solve(model, net.bus.shape[0])
        # the regulated (remote) bus magnitude must equal the generator setpoint
        self.assertLess(abs(abs(V[self.reg_bus]) - vset), self.tol,
                        f"|V[{self.reg_bus}]|={abs(V[self.reg_bus]):.8f} != vset={vset:.8f}")

    def _model(self):
        return _model()

    def test_controller_bus_becomes_pq(self):
        # once it regulates a remote bus, the controller's own bus is no longer PV:
        # it gets a reactive (Q) unknown column in the augmented Jacobian.
        net, model = self._model()
        model.set_gen_regulated_bus(self.gen_id, self.reg_bus)
        model.tell_solver_need_reset()
        self._solve(model, net.bus.shape[0])
        q_cols = model.get_solver().get_q_to_J_col()
        self.assertGreaterEqual(q_cols[self.ctrl_bus], 0,
                                "controller bus has no reactive unknown column")

    def test_kirchhoff_balance(self):
        net, model = self._model()
        model.set_gen_regulated_bus(self.gen_id, self.reg_bus)
        model.tell_solver_need_reset()
        V = self._solve(model, net.bus.shape[0])
        mismatch = model.check_solution(V, False)
        self.assertLess(np.abs(np.real(mismatch)).max(), 1e-6)
        self.assertLess(np.abs(np.imag(mismatch)).max(), 1e-6)

    def test_quadratic_convergence(self):
        net, model = self._model()
        model.set_gen_regulated_bus(self.gen_id, self.reg_bus)
        model.tell_solver_need_reset()
        self._solve(model, net.bus.shape[0])
        nb_iter = model.get_solver().get_nb_iter()
        self.assertLessEqual(nb_iter, 6, f"too many NR iterations ({nb_iter})")

    def test_res_q_written_back(self):
        # the converged reactive output of the remote gen must be reported on res_q
        net, model = self._model()
        model.set_gen_regulated_bus(self.gen_id, self.reg_bus)
        model.tell_solver_need_reset()
        self._solve(model, net.bus.shape[0])
        g = model.get_generators()[self.gen_id]
        self.assertTrue(g.has_res)
        self.assertTrue(np.isfinite(g.res_q_mvar))
        # it actually injects reactive power to hold the remote bus up
        self.assertGreater(abs(g.res_q_mvar), 1e-6)

    def test_two_gens_equal_share(self):
        # two gens (equal reactive range) regulating the SAME remote bus must split
        # the reactive injection equally (the bordered sharing row enforces it).
        net, model = self._model()
        # gens 2 (bus5) and 3 (bus7) both have range 30 MVAr; regulate bus 9 at 1.05
        for gid in (2, 3):
            model.change_v_gen(gid, 1.05)
            model.set_gen_regulated_bus(gid, self.reg_bus)
        model.tell_solver_need_reset()
        V = self._solve(model, net.bus.shape[0])
        self.assertLess(abs(abs(V[self.reg_bus]) - 1.05), self.tol)
        q2 = model.get_generators()[2].res_q_mvar
        q3 = model.get_generators()[3].res_q_mvar
        self.assertLess(abs(q2 - q3), 1e-5, f"equal-range gens should split Q equally: {q2} vs {q3}")

    def test_backward_compatible_when_local(self):
        # with no remote control configured, the result and the Jacobian sparsity
        # must be bit-identical to a model that never touched the new API.
        net, plain = self._model()
        _, same = self._model()
        same.set_gen_regulated_bus(self.gen_id, self.ctrl_bus)  # regulate own bus == local
        same.tell_solver_need_reset()
        nb = net.bus.shape[0]
        Vp = self._solve(plain, nb)
        Vs = self._solve(same, nb)
        self.assertEqual(np.abs(Vp - Vs).max(), 0.0)
        self.assertEqual(plain.get_solver().get_J().nnz, same.get_solver().get_J().nnz)
        self.assertEqual(plain.get_solver().get_nb_iter(), same.get_solver().get_nb_iter())

    def test_regulate_slack_raises(self):
        # the slack bus has no voltage unknown: regulating it is rejected (Phase 0 #3)
        net, model = self._model()
        model.set_gen_regulated_bus(self.gen_id, 0)  # bus 0 is the slack
        model.tell_solver_need_reset()
        with self.assertRaises(Exception):
            model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)


if __name__ == "__main__":
    unittest.main()
