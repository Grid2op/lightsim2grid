# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
GridPowerFlow -- one AC powerflow of one grid, as a differentiable pytorch operation.

    pf = GridPowerFlow(grid)
    V = pf(load_p=load_p, gen_v=gen_v)          # complex (n_bus,)
    loss = pf.compute_flows(V).i_or_a.relu().sum()
    loss.backward()                             # load_p.grad, gen_v.grad

The single-grid counterpart of :class:`BatchCPUPowerFlow`: every input is a
``(n_elements,)`` tensor, or ``None`` for "the grid's own value"; the forward is the
grid's own :func:`LSGrid.ac_pf`, outer loops included (``NROuter_*`` /
``NROuterDist_*``, see :func:`LSGrid.set_outer_loops`). The backward is one transposed
solve of the last Newton's Jacobian (:func:`LSGrid.solve_JT`), mapped onto the inputs
by the same code as the batch's.

Nothing is kept by the forward beyond what the grid's algorithm holds anyway: the
backward reads the system of the grid's last powerflow. A graph whose forward is not
that one any more (another forward ran since, as torch's ``gradcheck`` does) has its
own powerflow run again first, and checked to land on the same voltages.
"""

import numpy as np

from lightsim2grid.algorithm import AlgorithmType

from ._flows import compute_branch_flows, branch_admittances
from . import _adjoint
from ._batch_cpu_power_flow import _torch

__all__ = ["GridPowerFlow"]


class GridPowerFlow:
    """One AC powerflow, differentiable with respect to the grid's injections and
    generator set-points. See the module docstring.

    The grid is copied: the forward writes its inputs into the copy, never into the
    grid it was given.
    """

    def __init__(self, grid, *, algorithm=None, max_iter=10, tol=1e-8, v_init=None):
        torch = _torch()
        self._grid = grid.copy()
        if algorithm is not None:
            self._grid.change_algorithm(algorithm)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.sn_mva = float(self._grid.get_sn_mva())
        self.n_bus = int(self._grid.total_bus())
        self._v_init = (np.ones(self.n_bus, dtype=complex) if v_init is None
                        else np.asarray(v_init, dtype=complex).copy())
        if self._v_init.shape != (self.n_bus,):
            raise ValueError(f"`v_init` must have one entry per grid bus ({self.n_bus}), "
                             f"got {self._v_init.shape}.")

        g = self._grid
        loads, gens, sgens = list(g.get_loads()), list(g.get_generators()), list(g.get_static_generators())
        self._load_bus = np.asarray(g.get_loads().get_bus_id(), dtype=np.int64)
        self._gen_bus = np.asarray(g.get_generators().get_bus_id(), dtype=np.int64)
        self._sgen_bus = np.asarray(g.get_static_generators().get_bus_id(), dtype=np.int64)
        self.n_load, self.n_gen, self.n_sgen = len(loads), len(gens), len(sgens)
        # the grid's own values: what an input left at None means, at every call
        self._base = {
            "load_p": np.array([el.target_p_mw for el in loads], dtype=np.float64),
            "load_q": np.array([el.target_q_mvar for el in loads], dtype=np.float64),
            "gen_p": np.array([el.target_p_mw for el in gens], dtype=np.float64),
            "sgen_p": np.array([el.target_p_mw for el in sgens], dtype=np.float64),
            "gen_v": np.array([el.target_vm_pu for el in gens], dtype=np.float64),
        }
        self._setters = {"load_p": g.change_p_load, "load_q": g.change_q_load,
                         "gen_p": g.change_p_gen, "sgen_p": g.change_p_sgen,
                         "gen_v": g.change_v_gen}

        lines, trafos = g.get_lines(), g.get_trafos()
        self.n_line = np.asarray(lines.get_bus_id_side_1()).shape[0]
        self.n_trafo = np.asarray(trafos.get_bus_id_side_1()).shape[0]
        self.n_branch = self.n_line + self.n_trafo

        def _cat(attr):
            return np.concatenate([np.asarray(getattr(lines, attr)()),
                                   np.asarray(getattr(trafos, attr)())])

        branch_or = _cat("get_bus_id_side_1").astype(np.int64)
        branch_ex = _cat("get_bus_id_side_2").astype(np.int64)
        self._branch_alive = torch.as_tensor((branch_or >= 0) & (branch_ex >= 0))
        self._branch_or = torch.as_tensor(np.clip(branch_or, 0, None))
        self._branch_ex = torch.as_tensor(np.clip(branch_ex, 0, None))
        self._y_11 = torch.as_tensor(_cat("get_yac_eff_11"), dtype=torch.complex128)
        self._y_12 = torch.as_tensor(_cat("get_yac_eff_12"), dtype=torch.complex128)
        self._y_21 = torch.as_tensor(_cat("get_yac_eff_21"), dtype=torch.complex128)
        self._y_22 = torch.as_tensor(_cat("get_yac_eff_22"), dtype=torch.complex128)
        self._bus_vn_kv = torch.as_tensor(np.asarray(g.get_bus_vn_kv()), dtype=torch.float64)

        self._solve_id = 0        # which forward the grid's algorithm holds
        self._converged = None
        self._trafo_moved = None

    @classmethod
    def init_from_grid(cls, grid, *, algorithm=None, **kwargs):
        """Build from an :class:`LSGrid` (copied). ``algorithm`` defaults to the grid's
        own; it must be an AC Newton-Raphson (``NR_*``, ``NROuter_*``,
        ``NROuterDist_*``), since the adjoint is a solve against its Jacobian."""
        return cls(grid, algorithm=algorithm, **kwargs)

    # ----------------------------------------------------------------- forward
    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)

    def forward(self, load_p=None, load_q=None, gen_p=None, sgen_p=None, gen_v=None):
        """Solve; returns complex ``V`` of shape ``(n_bus,)``, per unit, in grid bus
        numbering -- NaN where the powerflow did not converge.

        load_p, load_q : (n_load,) MW / MVAr     differentiable
        gen_p          : (n_gen,)  MW            differentiable
        sgen_p         : (n_sgen,) MW            differentiable
        gen_v          : (n_gen,)  vm_pu         differentiable; generators on one
                                                 bus must agree (see BatchCPUPowerFlow)
        """
        inputs = {"load_p": self._as_input(load_p, self.n_load, "load_p"),
                  "load_q": self._as_input(load_q, self.n_load, "load_q"),
                  "gen_p": self._as_input(gen_p, self.n_gen, "gen_p"),
                  "sgen_p": self._as_input(sgen_p, self.n_sgen, "sgen_p"),
                  "gen_v": self._as_input(gen_v, self.n_gen, "gen_v")}
        return _GridPowerFlowOp.apply(inputs["load_p"], inputs["load_q"], inputs["gen_p"],
                                      inputs["sgen_p"], inputs["gen_v"], self)

    # ----------------------------------------------------------------- results
    def converged(self):
        """Whether the last forward converged (if not, ``V`` is NaN and the gradient 0)."""
        if self._converged is None:
            raise RuntimeError("nothing has been computed yet: call this instance first.")
        return bool(self._converged)

    def compute_flows(self, V):
        """Branch flows from ``V`` (``(n_bus,)``), differentiable -- powerlines then
        transformers, each ``(n_branch,)``. A transformer whose tap the outer loops of
        the last forward moved is taken at that tap."""
        V2 = V.unsqueeze(0)
        ys = (self._y_11, self._y_12, self._y_21, self._y_22)
        if self._trafo_moved is not None:
            ids, blocks = self._trafo_moved
            ys = branch_admittances(*ys, 1, self.n_line, np.zeros(len(ids), dtype=np.int64), ids, blocks)
        flows = compute_branch_flows(V2, *ys, self._branch_or, self._branch_ex, self._bus_vn_kv,
                                     self.sn_mva, self._branch_alive.unsqueeze(0))
        return type(flows)(*(f[0] for f in flows))

    @property
    def grid(self):
        """The grid this object solves (its own copy): its results are the last forward's."""
        return self._grid

    # ----------------------------------------------------------------- helpers
    def _as_input(self, x, n_cols, name):
        torch = _torch()
        if x is None:
            return None
        if not isinstance(x, torch.Tensor):
            x = torch.as_tensor(np.asarray(x, dtype=np.float64))
        if x.device.type != "cpu":
            raise ValueError(f"'{name}' is on device '{x.device}': GridPowerFlow solves on the "
                             f"CPU. Pass `{name}.cpu()`.")
        x = x.to(dtype=torch.float64)
        if tuple(x.shape) != (n_cols,):
            raise ValueError(f"'{name}' must have shape ({n_cols},), got {tuple(x.shape)}.")
        return x

    def _solve(self, values):
        """Write the inputs, run the powerflow; returns its voltages (empty when it did
        not converge). The grid's algorithm then holds this solve's system."""
        self._write_inputs(values)
        res = np.asarray(self._grid.ac_pf(self._v_init, self.max_iter, self.tol))
        self._solve_id += 1
        return res

    def _write_inputs(self, values):
        """Every input onto the grid: the given value, or the grid's own (an input a
        previous call gave must not linger in this one). A no-op change costs nothing."""
        for name, base in self._base.items():
            value = values.get(name)
            target = base if value is None else value
            setter = self._setters[name]
            for k in range(target.shape[0]):
                setter(k, float(target[k]))


class _GridPowerFlowOp:
    """Built lazily, see _BatchCPUPowerFlowOp."""

    _impl = None

    @classmethod
    def apply(cls, load_p, load_q, gen_p, sgen_p, gen_v, pf):
        if cls._impl is None:
            cls._impl = _make_op(_torch())
        return cls._impl.apply(load_p, load_q, gen_p, sgen_p, gen_v, pf)


def _make_op(torch):
    class _Op(torch.autograd.Function):
        @staticmethod
        def forward(ctx, load_p, load_q, gen_p, sgen_p, gen_v, pf):
            names = ("load_p", "load_q", "gen_p", "sgen_p", "gen_v")
            tensors = (load_p, load_q, gen_p, sgen_p, gen_v)
            values = {n: t.detach().numpy().copy() for n, t in zip(names, tensors) if t is not None}
            res = pf._solve(values)
            pf._converged = res.shape[0] == pf.n_bus
            pf._trafo_moved = pf._grid.get_outer_trafo_admittances() if pf._converged else None
            V = (torch.as_tensor(res) if pf._converged
                 else torch.full((pf.n_bus,), float("nan"), dtype=torch.complex128))
            ctx.pf = pf
            ctx.values = values
            ctx.solve_id = pf._solve_id
            ctx.converged = pf._converged
            ctx.needs = tuple(t is not None for t in tensors)
            ctx.save_for_backward(V)
            return V

        @staticmethod
        def backward(ctx, grad_V):
            pf = ctx.pf
            V, = ctx.saved_tensors
            if not ctx.converged:
                return tuple(torch.zeros(n, dtype=torch.float64) if need else None
                             for need, n in zip(ctx.needs, (pf.n_load, pf.n_load, pf.n_gen,
                                                            pf.n_sgen, pf.n_gen))) + (None,)
            if ctx.solve_id != pf._solve_id:
                # another forward ran since: this graph's own powerflow, again, so that the
                # grid holds its system -- the same inputs from the same start give the same
                # solve, which the voltages confirm
                res = pf._solve(ctx.values)
                if res.shape[0] != pf.n_bus or not np.allclose(res, V.numpy(), rtol=0., atol=1e3 * pf.tol):
                    raise RuntimeError(
                        "GridPowerFlow: the powerflow of this graph, run again for its backward, "
                        "does not land where its forward did -- the grid was changed in between "
                        "(through `.grid`?). Call backward before changing it.")
                ctx.solve_id = pf._solve_id
            grads = _adjoint.input_grads(torch, _GridSolver(pf._grid), V.unsqueeze(0), grad_V.unsqueeze(0),
                                         ctx.needs, pf.sn_mva, pf._load_bus, pf._gen_bus, pf._sgen_bus)
            return tuple(None if g is None else g[0] for g in grads) + (None,)

    return _Op


class _GridSolver:
    """What :func:`_adjoint.input_grads` reads of one grid's last solve, as a batch of
    one row."""

    def __init__(self, grid):
        self._grid = grid
        self.dim_J = grid.dim_J()
        self.theta_col = np.asarray(grid.get_theta_col_of_bus())
        self.vm_col = np.asarray(grid.get_vm_col_of_bus())
        self.p_row = np.asarray(grid.get_p_row_of_bus())
        self.q_row = np.asarray(grid.get_q_row_of_bus())
        self.gen_v_target_bus = np.asarray(grid.get_gen_v_target_bus())

    def solve_JT(self, xbar):
        return np.asarray(self._grid.solve_JT(xbar))

    def gen_v_indirect_grad(self, lam):
        return np.asarray(self._grid.gen_v_indirect_grad(lam[0]))[None, :]

    def gen_p_indirect_grad(self, lam):
        return np.asarray(self._grid.gen_p_indirect_grad(lam[0]))[None, :]

    def gen_v_share(self):
        return np.asarray(self._grid.get_gen_v_share())[None, :]
