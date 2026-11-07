# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""GridPowerFlow: one grid's powerflow, outer loops included, as pytorch sees it.

As for BatchCPUPowerFlow, the reference for a gradient is a central finite difference
of the same loss, the whole powerflow re-run; and the single grid's gradient is also
the one a one-row batch gives for the same solve.
"""

import unittest
import warnings

import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

try:
    import pandapower.networks  # noqa: F401
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False


@unittest.skipUnless(TORCH_AVAILABLE and PP_AVAILABLE, "needs pytorch and pandapower")
class TestGridPowerFlow(unittest.TestCase):
    TOL = 1e-12

    def setUp(self):
        from lightsim2grid.differentiable import GridPowerFlow, BatchCPUPowerFlow
        from test_outer_batch import OUTER
        self._cls = GridPowerFlow
        self._batch_cls = BatchCPUPowerFlow
        self._outer = OUTER
        warnings.filterwarnings("ignore")
        grid = self._grid()
        self.load_p = 1.1 * np.asarray(grid.get_load_target_p()).copy()
        self.load_q = 0.3 * self.load_p
        self.gen_p = np.asarray(grid.get_gen_target_p()).copy()
        self.gen_v = np.asarray([g.target_vm_pu for g in grid.get_generators()])

    def _grid(self):
        """test_outer_batch's case14, short of reactive power: ReactiveLimits switches buses PQ"""
        from test_outer_batch import _grid as _outer_grid
        from lightsim2grid.algorithm import ReactiveLimits
        grid = _outer_grid(distributed=False)
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits())
        return grid

    def _make(self, grid=None):
        return self._cls(grid if grid is not None else self._grid(), tol=self.TOL, max_iter=30)

    @staticmethod
    def _loss(V):
        n_bus = V.shape[-1]
        wr = torch.as_tensor(1. + 0.5 * np.arange(n_bus))
        wi = torch.as_tensor(-0.75 + 0.25 * np.arange(n_bus))
        return (wr * V.real + wi * V.imag).sum()

    def _inputs(self):
        return dict(load_p=self.load_p, load_q=self.load_q, gen_p=self.gen_p, gen_v=self.gen_v)

    def _fd(self, name, idx, delta, **inputs):
        out = []
        for sign in (+1., -1.):
            x = {k: v.copy() for k, v in inputs.items()}
            x[name][idx] += sign * delta
            out.append(float(self._loss(self._make()(**{k: torch.tensor(v) for k, v in x.items()}))))
        return (out[0] - out[1]) / (2. * delta)

    def test_gradient_of_every_input(self):
        inputs = self._inputs()
        pf = self._make()
        tensors = {k: torch.tensor(v, requires_grad=True) for k, v in inputs.items()}
        V = pf(**tensors)
        self.assertTrue(pf.converged())
        self._loss(V).backward()
        # what this is about: the loops switched some bus PQ, whose set-point reaches nothing
        self.assertTrue((np.asarray(pf.grid.get_gen_v_share()) == 0.).any())
        for name, t in tensors.items():
            delta = 1e-6 if name == "gen_v" else 1e-4
            for k in range(t.shape[0]):
                fd = self._fd(name, k, delta, **inputs)
                self.assertAlmostEqual(t.grad[k].item(), fd, places=5 if name == "gen_v" else 6,
                                       msg=f"{name}[{k}]: adjoint {t.grad[k].item()} vs fd {fd}")

    def test_same_gradient_as_a_one_row_batch(self):
        inputs = self._inputs()
        single = {k: torch.tensor(v, requires_grad=True) for k, v in inputs.items()}
        self._loss(self._make()(**single)).backward()
        batch = {k: torch.tensor(v[None, :], requires_grad=True) for k, v in inputs.items()}
        pf = self._batch_cls.init_from_grid(self._grid(), algorithm=self._outer, tol=self.TOL, max_iter=30)
        self._loss(pf(**batch)).backward()
        for name in inputs:
            np.testing.assert_allclose(single[name].grad.numpy(), batch[name].grad[0].numpy(),
                                       rtol=1e-8, atol=1e-10, err_msg=name)

    def test_gradcheck(self):
        pf = self._make()
        load_p = torch.tensor(self.load_p, requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(lambda lp: self._loss(pf(load_p=lp)), (load_p,),
                                                 eps=1e-6, atol=1e-5))

    def test_the_given_grid_is_not_modified(self):
        grid = self._grid()
        before = np.asarray(grid.get_load_target_p()).copy()
        self._make(grid)(load_p=torch.tensor(2. * self.load_p))
        np.testing.assert_array_equal(np.asarray(grid.get_load_target_p()), before)

    def test_a_call_does_not_depend_on_the_previous_one(self):
        pf = self._make()
        V0 = pf(load_q=torch.tensor(self.load_q))
        pf(load_p=torch.tensor(1.2 * self.load_p), load_q=torch.tensor(self.load_q))
        # load_p left out again: the grid's own, not the previous call's
        V2 = pf(load_q=torch.tensor(self.load_q))
        np.testing.assert_allclose(V2.numpy(), V0.numpy(), rtol=0., atol=1e-10)

    def test_backward_of_an_older_forward(self):
        # two graphs, the older one differentiated last: it gets its own gradient
        pf = self._make()
        load_p = torch.tensor(self.load_p, requires_grad=True)
        loss_first = self._loss(pf(load_p=load_p))
        other = torch.tensor(1.01 * self.load_p, requires_grad=True)
        self._loss(pf(load_p=other)).backward()
        loss_first.backward()
        ref = torch.tensor(self.load_p, requires_grad=True)
        self._loss(self._make()(load_p=ref)).backward()
        np.testing.assert_allclose(load_p.grad.numpy(), ref.grad.numpy(), rtol=1e-10, atol=1e-12)
        self.assertGreater(np.abs(other.grad.numpy() - ref.grad.numpy()).max(), 1e-8)

    def test_backward_after_the_grid_changed_is_refused(self):
        pf = self._make()
        load_p = torch.tensor(self.load_p, requires_grad=True)
        loss_first = self._loss(pf(load_p=load_p))
        self._loss(pf(load_p=load_p * 1.01))
        pf.grid.deactivate_powerline(3)   # not an input: the rerun cannot undo it
        with self.assertRaises(RuntimeError):
            loss_first.backward()

    def test_a_powerflow_that_does_not_converge_is_nan_and_has_no_gradient(self):
        pf = self._make()
        load_p = torch.tensor(1e5 * self.load_p, requires_grad=True)
        V = pf(load_p=load_p)
        self.assertFalse(pf.converged())
        self.assertTrue(torch.isnan(V.real).all())
        V.real.nan_to_num(0.).sum().backward()
        self.assertTrue((load_p.grad == 0.).all())

    def test_flows_at_the_taps_the_loops_settled_on(self):
        try:
            from test_outer_batch import _tap_grid
            import pypowsybl  # noqa: F401
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        pf = self._make(_tap_grid())
        V = pf(load_p=torch.tensor(np.asarray(pf.grid.get_load_target_p()).copy()))
        self.assertTrue(pf.converged())
        ids, _ = pf.grid.get_outer_trafo_admittances()
        self.assertGreater(len(ids), 0)        # a tap did move
        flows = pf.compute_flows(V)
        # the grid's own results of that powerflow, at the taps it settled on
        p_or = flows.p_or_mw.detach().numpy()
        np.testing.assert_allclose(p_or[:pf.n_line], np.asarray(pf.grid.get_line_res1()[0]), rtol=1e-9, atol=1e-8)
        np.testing.assert_allclose(p_or[pf.n_line:], np.asarray(pf.grid.get_trafo_res1()[0]), rtol=1e-9, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
