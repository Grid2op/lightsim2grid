# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""BatchCPUPowerFlow: the batch adjoint, as pytorch sees it.

The reference for a gradient here is always a central finite difference of the
same quantity, computed by re-running the batch -- something that knows nothing
about Jacobians, ledgers or adjoint signs. Where the gradient flows through
compute_flows(), the flows themselves are checked against the C++ ones on the
same voltages first, so a failure can be told apart from a failure of the
formulas they go through.
"""

import subprocess
import sys
import unittest
import warnings

import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

try:
    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False


def _two_gens_on_one_bus():
    """case14 with a second generator added on bus 1, regulating it as well. Both
    set-points then aim at the same bus, and only the last one applied reaches the
    solve (VoltageSourceContainer::set_vm is last-writer-wins)."""
    import pandapower as pp
    net = pn.case14()
    bus = int(net.gen.bus.iloc[0])
    pp.create_gen(net, bus=bus, p_mw=float(net.gen.p_mw.iloc[0]) * 0.5,
                  vm_pu=float(net.gen.vm_pu.iloc[0]), min_q_mvar=-500., max_q_mvar=500.,
                  controllable=True)
    return net


@unittest.skipUnless(TORCH_AVAILABLE and PP_AVAILABLE, "needs pytorch and pandapower")
class TestBatchCPUPowerFlow(unittest.TestCase):
    N_SCEN = 4
    TOL = 1e-12

    def setUp(self):
        from lightsim2grid.differentiable import BatchCPUPowerFlow
        self._cls = BatchCPUPowerFlow
        warnings.filterwarnings("ignore")
        self.pf = self._make()
        rng = np.random.default_rng(0)
        base_lp = np.asarray(self.pf._grid.get_load_target_p())
        base_gp = np.asarray(self.pf._grid.get_gen_target_p())
        self.load_p = base_lp[None, :] * (1. + 0.15 * rng.standard_normal((self.N_SCEN, self.pf.n_load)))
        self.load_q = 0.2 * self.load_p
        self.gen_p = base_gp[None, :] * (1. + 0.05 * rng.standard_normal((self.N_SCEN, self.pf.n_gen)))

    def _make(self, **kwargs):
        grid = init_from_pandapower(pn.case14())
        return self._cls.init_from_grid(grid, tol=self.TOL, **kwargs)

    # the loss every gradient below is taken of: linear in Re/Im so that BOTH the
    # angle and the magnitude of every bus carry a cotangent (|V|^2 would leave the
    # angle half of the ledger untested), and bounded away from any kink
    @staticmethod
    def _loss(V):
        n_bus = V.shape[1]
        wr = torch.as_tensor(1. + 0.5 * np.arange(n_bus))
        wi = torch.as_tensor(-0.75 + 0.25 * np.arange(n_bus))
        return (wr * V.real + wi * V.imag).sum()

    def _run(self, pf=None, **inputs):
        pf = pf if pf is not None else self._make()
        tensors = {k: torch.tensor(v) if v is not None else None for k, v in inputs.items()}
        return pf, tensors, pf(**tensors)

    def _fd(self, name, idx, delta=1e-4, loss=None, **inputs):
        """Central finite difference of the loss with respect to one input entry,
        each evaluation on its own batch."""
        loss = loss if loss is not None else (lambda pf, V: self._loss(V))
        out = []
        for sign in (+1., -1.):
            perturbed = {k: (v.copy() if v is not None else None) for k, v in inputs.items()}
            perturbed[name][idx] += sign * delta
            pf, _, V = self._run(**perturbed)
            out.append(float(loss(pf, V)))
        return (out[0] - out[1]) / (2. * delta)

    # ------------------------------------------------------------------ gradient
    def test_gradient_of_every_differentiable_input(self):
        inputs = dict(load_p=self.load_p, load_q=self.load_q, gen_p=self.gen_p)
        pf, tensors, V = self._run(pf=self.pf, **inputs)
        for t in tensors.values():
            t.requires_grad_(True)
        V = pf(**tensors)
        self.assertTrue(pf.converged().all())
        self._loss(V).backward()

        for name in ("load_p", "load_q", "gen_p"):
            grad = tensors[name].grad
            self.assertIsNotNone(grad, name)
            for idx in ((0, 0), (2, 1), (self.N_SCEN - 1, grad.shape[1] - 1)):
                fd = self._fd(name, idx, **inputs)
                self.assertAlmostEqual(
                    grad[idx].item(), fd, places=6,
                    msg=f"{name}{idx}: adjoint {grad[idx].item()} vs finite difference {fd}")

    def test_gradcheck(self):
        load_p = torch.tensor(self.load_p[:2], requires_grad=True)
        pf = self._make()
        self.assertTrue(torch.autograd.gradcheck(
            lambda lp: self._loss(pf(load_p=lp)), (load_p,), eps=1e-6, atol=1e-5))

    def test_gradient_flows_through_the_currents(self):
        load_p = torch.tensor(self.load_p, requires_grad=True)
        pf = self._make()
        loss = pf.compute_flows(pf(load_p=load_p)).i_or_a.sum()
        loss.backward()
        fd = self._fd("load_p", (1, 2), load_p=self.load_p,
                      loss=lambda pf_, V: pf_.compute_flows(V).i_or_a.sum())
        self.assertAlmostEqual(load_p.grad[1, 2].item(), fd, places=4)

    # ------------------------------------------------------------------ the flows
    def test_flows_agree_with_the_cpp_ones(self):
        pf, _, V = self._run(load_p=self.load_p, load_q=self.load_q)
        flows = pf.compute_flows(V)
        np.testing.assert_allclose(flows.p_or_mw.detach().numpy(),
                                   pf.sweep.compute_power_flows(), rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(flows.i_or_a.detach().numpy(),
                                   1000. * np.asarray(pf.sweep.compute_flows()),
                                   rtol=1e-9, atol=1e-9)

    # --------------------------------------------------------------- contingencies
    def test_gradient_with_a_line_contingency(self):
        line_status = np.ones((self.N_SCEN, self.pf.n_line), dtype=bool)
        line_status[:, 3] = False
        line_status[1, 5] = False
        inputs = dict(load_p=self.load_p, line_status=line_status)

        load_p = torch.tensor(self.load_p, requires_grad=True)
        pf = self._make()
        V = pf(load_p=load_p, line_status=torch.as_tensor(line_status))
        self.assertTrue(pf.converged().all())
        self._loss(V).backward()
        for idx in ((0, 0), (1, 4)):
            fd = self._fd("load_p", idx, **inputs)
            self.assertAlmostEqual(load_p.grad[idx].item(), fd, places=6)

    def test_gradient_with_a_generator_contingency(self):
        gen_status = np.ones((self.N_SCEN, self.pf.n_gen), dtype=bool)
        gen_status[2, 3] = False
        inputs = dict(load_p=self.load_p, gen_status=gen_status)

        load_p = torch.tensor(self.load_p, requires_grad=True)
        pf = self._make()
        V = pf(load_p=load_p, gen_status=torch.as_tensor(gen_status))
        self.assertTrue(pf.converged().all())
        self._loss(V).backward()
        for idx in ((2, 0), (0, 1)):
            fd = self._fd("load_p", idx, **inputs)
            self.assertAlmostEqual(load_p.grad[idx].item(), fd, places=6)

    def test_a_call_does_not_depend_on_the_previous_one(self):
        """A mask set by one call must not leak into the next: the same arguments
        always give the same answer, whatever ran before."""
        line_status = np.ones((self.N_SCEN, self.pf.n_line), dtype=bool)
        line_status[:, 3] = False

        pf = self._make()
        pf(load_p=torch.tensor(self.load_p), line_status=torch.as_tensor(line_status))
        after_mask = pf(load_p=torch.tensor(self.load_p))          # no mask this time

        _, _, fresh = self._run(load_p=self.load_p)
        np.testing.assert_allclose(after_mask.detach().numpy(), fresh.detach().numpy(),
                                   rtol=1e-10, atol=1e-10)

        # ... and the other way round: a masked call after an unmasked one
        pf2 = self._make()
        pf2(load_p=torch.tensor(self.load_p))
        after_plain = pf2(load_p=torch.tensor(self.load_p), line_status=torch.as_tensor(line_status))
        _, _, fresh_masked = self._run(load_p=self.load_p, line_status=line_status)
        np.testing.assert_allclose(after_plain.detach().numpy(), fresh_masked.detach().numpy(),
                                   rtol=1e-10, atol=1e-10)

    def test_the_number_of_scenarios_can_change(self):
        pf = self._make()
        pf(load_p=torch.tensor(self.load_p))
        V2 = pf(load_p=torch.tensor(self.load_p[:2]))
        self.assertEqual(V2.shape[0], 2)
        _, _, fresh = self._run(load_p=self.load_p[:2])
        np.testing.assert_allclose(V2.detach().numpy(), fresh.detach().numpy(),
                                   rtol=1e-10, atol=1e-10)

    # -------------------------------------------------------------- failing rows
    def test_a_row_that_does_not_converge_is_nan_and_has_no_gradient(self):
        load_p = self.load_p.copy()
        load_p[1, :] = 1e6      # far outside anything solvable
        tensor = torch.tensor(load_p, requires_grad=True)
        pf = self._make()
        V = pf(load_p=tensor)

        converged = pf.converged()
        self.assertFalse(converged[1])
        self.assertTrue(torch.isnan(V[1]).all())
        self.assertFalse(torch.isnan(V[converged]).any())

        self._loss(torch.where(torch.isnan(V), torch.zeros_like(V), V)).backward()
        np.testing.assert_array_equal(tensor.grad[1].numpy(), np.zeros(pf.n_load))
        self.assertGreater(np.abs(tensor.grad[0].numpy()).max(), 0.)

    # ------------------------------------------------------------------ threading
    def test_threads_do_not_change_the_answer(self):
        big_p = np.repeat(self.load_p, 5, axis=0)
        _, _, v1 = self._run(pf=self._make(nb_thread=1), load_p=big_p)
        pf4 = self._make(nb_thread=4)
        lp = torch.tensor(big_p, requires_grad=True)
        v4 = pf4(load_p=lp)
        np.testing.assert_allclose(v4.detach().numpy(), v1.detach().numpy(), rtol=1e-12, atol=1e-12)

        lp1 = torch.tensor(big_p, requires_grad=True)
        self._loss(self._make(nb_thread=1)(load_p=lp1)).backward()
        self._loss(v4).backward()
        np.testing.assert_allclose(lp.grad.numpy(), lp1.grad.numpy(), rtol=1e-10, atol=1e-10)

    # --------------------------------------------------------------------- errors
    def test_rejects_inputs_that_do_not_line_up(self):
        pf = self._make()
        with self.assertRaises(ValueError):
            pf()                                             # nothing to size the batch by
        with self.assertRaises(ValueError):
            pf(load_p=torch.zeros(self.N_SCEN))              # not 2-D
        with self.assertRaises(ValueError):
            pf(load_p=torch.zeros(self.N_SCEN, self.pf.n_load + 1))
        with self.assertRaises(ValueError):
            pf(load_p=torch.zeros(3, self.pf.n_load), load_q=torch.zeros(4, self.pf.n_load))

    def test_a_gen_v_does_not_leak_into_the_next_call(self):
        """gen_v has no "nothing set" value to send back, so dropping it has to start
        the batch over -- otherwise the previous call's set-points stay applied."""
        pf = self._make()
        gen_v = np.full((self.N_SCEN, self.pf.n_gen), 1.06)
        pf(load_p=torch.tensor(self.load_p), gen_v=torch.tensor(gen_v))
        after = pf(load_p=torch.tensor(self.load_p))

        _, _, fresh = self._run(load_p=self.load_p)
        np.testing.assert_allclose(after.detach().numpy(), fresh.detach().numpy(),
                                   rtol=1e-10, atol=1e-10)

    def test_gen_v_is_usable(self):
        pf = self._make()
        gen_v = np.full((self.N_SCEN, self.pf.n_gen), 1.02)
        V = pf(load_p=torch.tensor(self.load_p), gen_v=torch.tensor(gen_v))
        self.assertTrue(pf.converged().all())
        np.testing.assert_allclose(np.abs(V.detach().numpy()[:, self._regulated_buses(pf)]),
                                   1.02, rtol=1e-8)

    @staticmethod
    def _regulating_groups(pf):
        """generator ids grouped by the bus they regulate, for the buses with more than
        one regulator. Asked of the grid rather than assumed: init_from_pandapower does
        not preserve pandapower's generator order."""
        by_bus = {}
        for g in pf._grid.get_generators():
            if not (g.connected and g.voltage_regulator_on):
                continue
            by_bus.setdefault(g.regulated_bus_id, []).append(g.id)
        return [ids for ids in by_bus.values() if len(ids) > 1]

    @staticmethod
    def _regulated_buses(pf):
        """the buses a gen_v set-point actually pins -- the ones whose magnitude it is
        fair to read back (see get_gen_v_target_bus)"""
        target = np.asarray(pf._sweep.get_gen_v_target_bus())
        return np.unique(target[target >= 0])

    # ------------------------------------------------------- gen_v is differentiable
    def _gen_v_inputs(self):
        """a set-point per row and per generator, all distinct, and away from the
        grid's own so nothing passes by accident"""
        rng = np.random.default_rng(7)
        return 1.02 + 0.02 * rng.random((self.N_SCEN, self.pf.n_gen))

    def test_gradient_of_gen_v(self):
        gen_v = self._gen_v_inputs()
        inputs = dict(load_p=self.load_p, gen_v=gen_v)
        pf = self._make()
        tensors = {k: torch.tensor(v, requires_grad=True) for k, v in inputs.items()}
        V = pf(**tensors)
        self.assertTrue(pf.converged().all())
        self._loss(V).backward()

        grad = tensors["gen_v"].grad
        self.assertIsNotNone(grad)
        target = np.asarray(pf._sweep.get_gen_v_target_bus())
        live = np.flatnonzero(target >= 0)
        self.assertGreater(live.size, 0, "case14 should have at least one live set-point")

        # every entry, live or dead, against a finite difference of the whole sweep
        for g in range(self.pf.n_gen):
            for row in (0, self.N_SCEN - 1):
                fd = self._fd("gen_v", (row, g), delta=1e-6, **inputs)
                self.assertAlmostEqual(
                    grad[row, g].item(), fd, places=5,
                    msg=f"gen_v[{row},{g}]: adjoint {grad[row, g].item()} vs fd {fd}")

    def test_a_set_point_that_never_reaches_the_solve_has_no_gradient(self):
        """Of several generators regulating one bus, only the last one's set-point is
        applied (set_vm is last-writer-wins). The others are dead inputs and their
        gradient must be exactly zero.

        Their set-points must AGREE, which is the constraint the grid now enforces: a
        bus has one magnitude. So they cannot be perturbed one at a time -- the group
        moves together, and what a finite difference then measures is the SUM of its
        gradients, which is the only part of the split the constraint leaves meaningful.

        case14 has one generator per bus and would never exercise any of this, so the
        grid here deliberately doubles one up -- the ordinary case on a real grid,
        where a busbar carries several machines."""
        def make():
            grid = init_from_pandapower(_two_gens_on_one_bus())
            return self._cls.init_from_grid(grid, tol=self.TOL)

        pf = make()
        n_scen, n_gen = 3, pf.n_gen
        rng = np.random.default_rng(3)
        load_p = np.asarray(pf._grid.get_load_target_p())[None, :] * (
            1. + 0.1 * rng.standard_normal((n_scen, pf.n_load)))
        gen_v = 1.01 + 0.02 * rng.random((n_scen, n_gen))
        groups = self._regulating_groups(pf)
        self.assertEqual(len(groups), 1, "the fixture should double up exactly one bus")
        group = groups[0]
        gen_v[:, group] = gen_v[:, group[0]][:, None]   # one bus, one magnitude
        inputs = dict(load_p=load_p, gen_v=gen_v)

        tensors = {k: torch.tensor(v, requires_grad=True) for k, v in inputs.items()}
        self._loss(pf(**tensors)).backward()
        self.assertTrue(pf.converged().all())

        target = np.asarray(pf._sweep.get_gen_v_target_bus())
        dead = np.flatnonzero(target < 0)
        live = np.flatnonzero(target >= 0)
        self.assertGreater(dead.size, 0, "the doubled-up generator should be a dead input")
        self.assertGreater(live.size, 0)

        def fd(name, idx, delta=1e-6):
            out = []
            for sign in (+1., -1.):
                pert = {k: v.copy() for k, v in inputs.items()}
                pert[name][idx] += sign * delta
                out.append(float(self._loss(make()(**{k: torch.tensor(v)
                                                      for k, v in pert.items()}))))
            return (out[0] - out[1]) / (2. * delta)

        for g in dead:
            np.testing.assert_array_equal(tensors["gen_v"].grad[:, g].numpy(),
                                          np.zeros(n_scen))
        # the group moves together: its gradients must sum to what that move costs
        def fd_group(delta=1e-6):
            out = []
            for sign in (+1., -1.):
                pert = {k: v.copy() for k, v in inputs.items()}
                for g in group:
                    pert["gen_v"][0, g] += sign * delta
                out.append(float(self._loss(make()(**{k: torch.tensor(v)
                                                      for k, v in pert.items()}))))
            return (out[0] - out[1]) / (2. * delta)
        self.assertAlmostEqual(sum(tensors["gen_v"].grad[0, g].item() for g in group),
                               fd_group(), places=5)
        # every generator alone on its bus keeps an ordinary, individually measurable one
        for g in live:
            if g in group:
                continue
            self.assertAlmostEqual(tensors["gen_v"].grad[0, g].item(),
                                   fd("gen_v", (0, int(g))), places=5)

    def test_a_row_asking_one_bus_for_two_magnitudes_does_not_converge(self):
        """A bus has one voltage magnitude, so two generators regulating it cannot be
        given two different set-points. Such a row is not solved -- reported like any
        row skipped before the solver -- rather than silently taking whichever
        generator set_vm happens to write last."""
        grid = init_from_pandapower(_two_gens_on_one_bus())
        pf = self._cls.init_from_grid(grid, tol=self.TOL)
        n_scen = 3
        load_p = np.asarray(pf._grid.get_load_target_p())[None, :] * np.ones((n_scen, 1))
        gen_v = np.full((n_scen, pf.n_gen), 1.03)
        group = self._regulating_groups(pf)[0]
        gen_v[1, group[-1]] = 1.05   # row 1 alone contradicts itself

        V = pf(load_p=torch.tensor(load_p), gen_v=torch.tensor(gen_v))
        conv = pf.converged()
        self.assertFalse(conv[1])
        self.assertTrue(conv[0] and conv[2])
        self.assertTrue(torch.isnan(V[1]).all())   # and carries no result at all

    def test_gradcheck_gen_v(self):
        gen_v = torch.tensor(self._gen_v_inputs()[:2], requires_grad=True)
        load_p = torch.tensor(self.load_p[:2])
        pf = self._make()
        self.assertTrue(torch.autograd.gradcheck(
            lambda gv: self._loss(pf(load_p=load_p, gen_v=gv)), (gen_v,),
            eps=1e-6, atol=1e-5))

    def test_gen_v_gradient_with_a_line_contingency(self):
        """the indirect half is taken on the ROW's admittance matrix, so a row that
        drops a line must get that row's gradient, not the base case's"""
        gen_v = self._gen_v_inputs()
        line_status = np.ones((self.N_SCEN, self.pf.n_line), dtype=bool)
        line_status[1, 3] = False
        inputs = dict(load_p=self.load_p, gen_v=gen_v, line_status=line_status)
        pf = self._make()
        tensors = {"load_p": torch.tensor(self.load_p),
                   "gen_v": torch.tensor(gen_v, requires_grad=True),
                   "line_status": torch.as_tensor(line_status)}
        self._loss(pf(**tensors)).backward()

        target = np.asarray(pf._sweep.get_gen_v_target_bus())
        live = np.flatnonzero(target >= 0)
        for g in live[:3]:
            fd = self._fd("gen_v", (1, int(g)), delta=1e-6, **inputs)
            self.assertAlmostEqual(tensors["gen_v"].grad[1, g].item(), fd, places=5,
                                   msg=f"row 1 (line 3 out), gen {g}")


class TestTorchStaysOptional(unittest.TestCase):
    """torch is an optional dependency, and the only thing that may break if it is
    missing is actually building a BatchCPUPowerFlow -- not importing lightsim2grid,
    and not importing the subpackage either."""

    def test_importing_lightsim2grid_does_not_import_torch(self):
        script = ("import sys, lightsim2grid, lightsim2grid.differentiable;"
                  "sys.exit(1 if 'torch' in sys.modules else 0)")
        completed = subprocess.run([sys.executable, "-c", script], capture_output=True)
        self.assertEqual(completed.returncode, 0,
                         "importing lightsim2grid (or its differentiable subpackage) pulled "
                         "torch in; it must stay lazy, see _batch_cpu_power_flow._torch")


if __name__ == "__main__":
    unittest.main()
