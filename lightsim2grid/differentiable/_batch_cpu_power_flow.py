# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
BatchCPUPowerFlow -- a batch of AC powerflows as a differentiable pytorch operation.

    pf = BatchCPUPowerFlow.init_from_grid(grid)
    V = pf(load_p=load_p, load_q=load_q, gen_p=gen_p)   # complex (n_scenarios, n_bus)
    loss = pf.compute_flows(V).i_or_a.relu().sum()
    loss.backward()                                     # load_p.grad, load_q.grad, ...

Every input is a ``(n_scenarios, n_elements)`` tensor, or ``None`` for "the grid's
own value in every row". Row ``i`` is one independent scenario: its own injections
and its own set of disconnected elements. The whole batch is one
:class:`ScenarioSweepCPP` call, so it pays one symbolic factorization and
refactorizes per row.

Differentiable: ``load_p``, ``load_q``, ``gen_p``, ``sgen_p``. Discrete and
carrying no gradient: ``line_status``, ``trafo_status``, ``gen_status`` (booleans,
True = connected, grid2op's convention). ``gen_v`` is accepted but is not
differentiable yet -- passing one that requires a gradient raises rather than
silently returning none for it.

How the gradient is obtained
----------------------------
Each row solves ``G(x_i ; u_i) = 0`` and returns ``V_i``. The implicit function
theorem gives the gradient of a scalar loss with respect to every injection from
ONE transposed solve per row (see ``BatchAdjoint`` on the C++ side):

    x̄_i     = the loss's cotangent on V, projected onto the row's unknowns
    λ_i     = J_iᵀ⁻¹ x̄_i                      (keep_jacobian + solve_JT)
    ∂L/∂Sbus_i = λ_i, read at each bus's P and Q equation

and from the per-unit bus injection to an element is an affine map this class
applies directly: a load subtracts from its bus, a generator adds to it, both
scaled by ``sn_mva``. The element inputs are therefore the autograd leaves, and
nothing about the assembly (a generator contingency re-weighting the distributed
slack, a bus released from PV to PQ) has to be mirrored in python -- it stays in
the C++ that the forward already runs.

A row that did not converge -- it diverged, or a contingency islanded it -- comes
back as NaN and gets a zero gradient; see :func:`converged`.
"""

import numpy as np

from lightsim2grid.lightsim2grid_cpp import ScenarioSweepCPP
from lightsim2grid.algorithm import AlgorithmType

from ._flows import compute_branch_flows, BranchFlows

__all__ = ["BatchCPUPowerFlow"]


def _torch():
    """Import torch on first use, with an error that says what to do.

    Deliberately not imported at module scope: torch is an optional dependency of
    lightsim2grid, and `import lightsim2grid` must not drag it in.
    """
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "BatchCPUPowerFlow needs pytorch, which lightsim2grid does not require: "
            "install it with `pip install torch`, or `pip install lightsim2grid[torch]`."
        ) from exc
    return torch


class BatchCPUPowerFlow:
    """A batch of AC powerflows, differentiable with respect to its injections.

    Build it with :func:`init_from_grid`, then call it like a function (or use
    :func:`forward`). See the module docstring.
    """

    def __init__(self, sweep, grid, *, handle_disconnected_grid=False,
                 max_iter=10, tol=1e-8, v_init=None):
        torch = _torch()
        self._sweep = sweep
        self._grid = grid
        self._handle_disconnected_grid = bool(handle_disconnected_grid)
        self.max_iter = int(max_iter)
        self.tol = float(tol)

        self.sn_mva = float(grid.get_sn_mva())
        self.n_bus = int(grid.total_bus())
        self._v_init = (np.ones(self.n_bus, dtype=complex) if v_init is None
                        else np.asarray(v_init, dtype=complex).copy())
        if self._v_init.shape != (self.n_bus,):
            raise ValueError(f"`v_init` must have one entry per grid bus ({self.n_bus}), "
                             f"got {self._v_init.shape}.")

        # --- element -> bus, once. -1 marks an element the grid has disconnected;
        # it contributes to no bus, so it can never carry a gradient.
        self._load_bus = np.asarray(grid.get_loads().get_bus_id(), dtype=np.int64)
        self._gen_bus = np.asarray(grid.get_generators().get_bus_id(), dtype=np.int64)
        self._sgen_bus = np.asarray(grid.get_static_generators().get_bus_id(), dtype=np.int64)
        self.n_load = self._load_bus.shape[0]
        self.n_gen = self._gen_bus.shape[0]
        self.n_sgen = self._sgen_bus.shape[0]

        lines, trafos = grid.get_lines(), grid.get_trafos()
        self.n_line = np.asarray(lines.get_bus_id_side_1()).shape[0]
        self.n_trafo = np.asarray(trafos.get_bus_id_side_1()).shape[0]
        self.n_branch = self.n_line + self.n_trafo

        # --- branch data for compute_flows(), powerlines then transformers
        def _cat(attr):
            return np.concatenate([np.asarray(getattr(lines, attr)()),
                                   np.asarray(getattr(trafos, attr)())])

        branch_or = _cat("get_bus_id_side_1").astype(np.int64)
        branch_ex = _cat("get_bus_id_side_2").astype(np.int64)
        # a branch the grid itself has disconnected has no bus at one end: clamp the
        # index so the gather is in range, and remember never to report a flow for it
        self._branch_alive = (branch_or >= 0) & (branch_ex >= 0)
        self._branch_or = torch.as_tensor(np.clip(branch_or, 0, None))
        self._branch_ex = torch.as_tensor(np.clip(branch_ex, 0, None))
        self._y_11 = torch.as_tensor(_cat("get_yac_eff_11"), dtype=torch.complex128)
        self._y_12 = torch.as_tensor(_cat("get_yac_eff_12"), dtype=torch.complex128)
        self._y_21 = torch.as_tensor(_cat("get_yac_eff_21"), dtype=torch.complex128)
        self._y_22 = torch.as_tensor(_cat("get_yac_eff_22"), dtype=torch.complex128)
        self._bus_vn_kv = torch.as_tensor(np.asarray(grid.get_bus_vn_kv()), dtype=torch.float64)

        # --- call-to-call state (see _apply_status)
        self._n_scen = None
        self._line_off = None
        self._trafo_off = None
        self._gen_off = None
        self._gen_v_set = False
        self._last_connected = None      # (n_scen, n_branch) bool, for compute_flows
        self._converged = None

    # ------------------------------------------------------------------- build
    @classmethod
    def init_from_grid(cls, grid, *, nb_thread=1, handle_disconnected_grid=False,
                       algorithm=None, max_iter=10, tol=1e-8, v_init=None,
                       init_from_n_powerflow=False):
        """Build from an :class:`LSGrid`.

        The grid does not have to be solved first: every batch begins with its own
        "n" warm-up powerflow. ``algorithm`` defaults to the fastest Newton-Raphson
        the build has (KLU when it was compiled in); it must be one of the AC NR
        family, since the adjoint is a solve against that Jacobian.
        """
        sweep = ScenarioSweepCPP(grid)
        if algorithm is None:
            available = sweep.available_default_algorithms()
            algorithm = (AlgorithmType.NR_KLU if AlgorithmType.NR_KLU in available
                         else AlgorithmType.NR_SparseLU)
        sweep.change_algorithm(algorithm)
        sweep.nb_thread = int(nb_thread)
        sweep.init_from_n_powerflow = bool(init_from_n_powerflow)
        if handle_disconnected_grid:
            sweep.handle_disconnected_grid = True
        sweep.keep_jacobian = True
        return cls(sweep, grid, handle_disconnected_grid=handle_disconnected_grid,
                   max_iter=max_iter, tol=tol, v_init=v_init)

    # ----------------------------------------------------------------- forward
    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)

    def forward(self, load_p=None, load_q=None, gen_p=None, sgen_p=None, gen_v=None,
                line_status=None, trafo_status=None, gen_status=None):
        """Solve one scenario per row; returns complex ``V`` of shape
        ``(n_scenarios, n_bus)`` in per unit, in grid bus numbering (the same
        numbering the grid's own accessors use).

        load_p, load_q : (n_scen, n_load) MW / MVAr        differentiable
        gen_p          : (n_scen, n_gen)  MW               differentiable
        sgen_p         : (n_scen, n_sgen) MW               differentiable
        gen_v          : (n_scen, n_gen)  vm_pu            NOT differentiable yet
        line_status    : (n_scen, n_line)  bool, True = connected
        trafo_status   : (n_scen, n_trafo) bool, True = connected
        gen_status     : (n_scen, n_gen)   bool, True = connected

        ``None`` means "whatever the grid itself says", in every row. At least one
        input must be given: it is what fixes the number of scenarios.
        """
        n_scen = self._infer_n_scen(load_p=load_p, load_q=load_q, gen_p=gen_p, sgen_p=sgen_p,
                                    gen_v=gen_v, line_status=line_status,
                                    trafo_status=trafo_status, gen_status=gen_status)

        # The row count is locked by the first setter of a batch, so a different one
        # means starting the batch over. So does dropping gen_v: unlike the three
        # contingency masks, it has no "nothing set" value to send, so the only way to
        # stop applying the previous call's set-points is to start again. clear() also
        # drops settings the constructor made, hence _reconfigure.
        if n_scen != self._n_scen or (gen_v is None and self._gen_v_set):
            self._sweep.clear()
            self._reconfigure()
            self._line_off = self._trafo_off = self._gen_off = None
            self._gen_v_set = False
            self._n_scen = n_scen

        load_p = self._as_input(load_p, self.n_load, n_scen, "load_p")
        load_q = self._as_input(load_q, self.n_load, n_scen, "load_q")
        gen_p = self._as_input(gen_p, self.n_gen, n_scen, "gen_p")
        sgen_p = self._as_input(sgen_p, self.n_sgen, n_scen, "sgen_p")

        if gen_v is not None:
            gen_v = self._as_input(gen_v, self.n_gen, n_scen, "gen_v")
            if gen_v.requires_grad:
                raise NotImplementedError(
                    "`gen_v` is not differentiable yet: its gradient needs the dS/dVm "
                    "column the Jacobian does not store for a voltage-fixed bus. Pass "
                    "`gen_v.detach()` to use it as a (fixed) set-point in the meantime.")

        self._apply_status(line_status, trafo_status, gen_status, n_scen)
        if gen_v is not None:
            self._sweep.modify_gen_v(np.ascontiguousarray(gen_v.detach().numpy()))
            self._gen_v_set = True

        return _BatchCPUPowerFlowOp.apply(load_p, load_q, gen_p, sgen_p, self)

    # ----------------------------------------------------------------- results
    def converged(self):
        """(n_scenarios,) bool: whether each row's powerflow converged. A row that
        did not is NaN in ``V`` and gets a zero gradient."""
        if self._converged is None:
            raise RuntimeError("nothing has been computed yet: call this instance first.")
        return self._converged.copy()

    def compute_flows(self, V, connected=None):
        """Branch flows from ``V``, differentiable -- powerlines then transformers.

        Returns a :class:`BranchFlows`. By default a branch is taken as connected
        where the last :func:`forward` had it connected; pass ``connected``
        ``(n_scen, n_branch)`` to say otherwise.
        """
        torch = _torch()
        if connected is None:
            connected = self._last_connected
            if connected is None:
                connected = torch.ones((V.shape[0], self.n_branch), dtype=torch.bool)
        connected = connected & torch.as_tensor(self._branch_alive).unsqueeze(0)
        return compute_branch_flows(V, self._y_11, self._y_12, self._y_21, self._y_22,
                                    self._branch_or, self._branch_ex, self._bus_vn_kv,
                                    self.sn_mva, connected)

    @property
    def sweep(self):
        """The underlying :class:`ScenarioSweepCPP` (escape hatch: timers, limit
        violations, the adjoint's own linear-solver counters)."""
        return self._sweep

    def adjoint_memory_bytes(self):
        """Bytes the kept Jacobians of the last call occupy. Grows as
        ``n_scenarios * nnz(J)``: worth watching before scaling a batch up."""
        return self._sweep.adjoint_memory_bytes()

    # ----------------------------------------------------------------- helpers
    def _reconfigure(self):
        """Re-apply what ``clear()`` drops. ``keep_jacobian`` and ``nb_thread``
        survive it; ``handle_disconnected_grid`` does not."""
        if self._handle_disconnected_grid:
            self._sweep.handle_disconnected_grid = True
        self._sweep.keep_jacobian = True

    @staticmethod
    def _infer_n_scen(**inputs):
        n = None
        for name, value in inputs.items():
            if value is None:
                continue
            shape = tuple(np.shape(value))
            if len(shape) != 2:
                raise ValueError(f"every input must be 2-D (n_scenarios, n_elements); "
                                 f"'{name}' has shape {shape}.")
            if n is None:
                n = int(shape[0])
            elif int(shape[0]) != n:
                raise ValueError(f"all inputs must have the same number of scenarios; "
                                 f"got {n} and {shape[0]} (for '{name}').")
        if n is None:
            raise ValueError(
                "BatchCPUPowerFlow needs at least one input (load_p, load_q, gen_p, "
                "sgen_p, gen_v, line_status, trafo_status or gen_status): it is what "
                "fixes the number of scenarios.")
        if n <= 0:
            raise ValueError("the number of scenarios must be > 0.")
        return n

    def _as_input(self, x, n_cols, n_scen, name):
        """Bring an input to a CPU float64 tensor of the expected shape, or leave it
        None (which means "do not set this axis at all", so the C++ keeps the grid's
        own value for every row)."""
        torch = _torch()
        if x is None:
            return None
        if not isinstance(x, torch.Tensor):
            x = torch.as_tensor(np.asarray(x, dtype=np.float64))
        if x.device.type != "cpu":
            raise ValueError(f"'{name}' is on device '{x.device}': BatchCPUPowerFlow "
                             "solves on the CPU, and moving a tensor here silently would "
                             "hide that. Use gpusim2grid for a GPU batch, or pass "
                             f"`{name}.cpu()`.")
        x = x.to(dtype=torch.float64)
        if tuple(x.shape) != (n_scen, n_cols):
            raise ValueError(f"'{name}' must have shape ({n_scen}, {n_cols}), "
                             f"got {tuple(x.shape)}.")
        return x

    def _as_status(self, x, n_cols, n_scen, name):
        torch = _torch()
        if x is None:
            return None
        if not isinstance(x, torch.Tensor):
            x = torch.as_tensor(np.asarray(x, dtype=bool))
        x = x.to(dtype=torch.bool, device="cpu")
        if tuple(x.shape) != (n_scen, n_cols):
            raise ValueError(f"'{name}' must have shape ({n_scen}, {n_cols}), "
                             f"got {tuple(x.shape)}.")
        return x

    def _apply_status(self, line_status, trafo_status, gen_status, n_scen):
        """Hand the three contingency masks to the sweep, and remember them.

        A call that drops a mask the previous one had must say so explicitly -- an
        all-False mask -- or the sweep would keep applying the old one. That is what
        makes a call's result depend only on its own arguments.
        """
        torch = _torch()
        line_status = self._as_status(line_status, self.n_line, n_scen, "line_status")
        trafo_status = self._as_status(trafo_status, self.n_trafo, n_scen, "trafo_status")
        gen_status = self._as_status(gen_status, self.n_gen, n_scen, "gen_status")

        for name, status, n_cols, attr, setter in (
                ("line", line_status, self.n_line, "_line_off", self._sweep.set_contingency_lines),
                ("trafo", trafo_status, self.n_trafo, "_trafo_off", self._sweep.set_contingency_trafos),
                ("gen", gen_status, self.n_gen, "_gen_off", self._sweep.set_contingency_gens)):
            if status is not None:
                off = np.ascontiguousarray(~status.numpy())
            elif getattr(self, attr) is not None:
                off = np.zeros((n_scen, n_cols), dtype=bool)   # explicitly: nothing tripped
            else:
                off = None                                      # never set, nothing to undo
            if off is not None:
                setter(off)
            setattr(self, attr, off)

        # what compute_flows() takes as "this branch carries something"
        off_branch = []
        for off, n_cols in ((self._line_off, self.n_line), (self._trafo_off, self.n_trafo)):
            off_branch.append(np.zeros((n_scen, n_cols), dtype=bool) if off is None else off)
        self._last_connected = ~torch.as_tensor(np.concatenate(off_branch, axis=1))


class _BatchCPUPowerFlowOp:
    """Built lazily: torch.autograd.Function must be subclassed, and torch may not
    be importable at import time. See :func:`apply`."""

    _impl = None

    @classmethod
    def apply(cls, load_p, load_q, gen_p, sgen_p, pf):
        if cls._impl is None:
            cls._impl = _make_op(_torch())
        return cls._impl.apply(load_p, load_q, gen_p, sgen_p, pf)


def _make_op(torch):
    class _Op(torch.autograd.Function):
        @staticmethod
        def forward(ctx, load_p, load_q, gen_p, sgen_p, pf):
            sweep = pf._sweep
            for tensor, setter in ((load_p, sweep.modify_load_p),
                                   (load_q, sweep.modify_load_q),
                                   (gen_p, sweep.modify_gen_p),
                                   (sgen_p, sweep.modify_sgen_p)):
                if tensor is not None:
                    setter(np.ascontiguousarray(tensor.detach().numpy()))

            sweep.compute(pf._v_init, pf.max_iter, pf.tol)

            converged = np.asarray(sweep.converged_mask(), dtype=bool)
            pf._converged = converged
            V = torch.as_tensor(np.array(sweep.get_voltages()))
            # a row that was never solved holds zeros, which is indistinguishable from
            # a genuine result: NaN says "there is nothing here" out loud, and
            # propagates into any loss that forgets to mask it
            V[~torch.as_tensor(converged)] = float("nan")

            ctx.pf = pf
            ctx.needs = tuple(t is not None for t in (load_p, load_q, gen_p, sgen_p))
            ctx.save_for_backward(V)
            return V

        @staticmethod
        def backward(ctx, grad_V):
            pf = ctx.pf
            sweep = pf._sweep
            V, = ctx.saved_tensors
            n_scen, n_bus = V.shape
            dim_J = sweep.dim_J()

            theta_col = np.asarray(sweep.get_theta_col_of_bus())
            vm_col = np.asarray(sweep.get_vm_col_of_bus())
            p_row = np.asarray(sweep.get_p_row_of_bus())
            q_row = np.asarray(sweep.get_q_row_of_bus())

            # NaN rows (never solved) carry no cotangent, and must not put a NaN into
            # the system either
            finite = torch.isfinite(V.real) & torch.isfinite(V.imag)
            gV = torch.where(finite, grad_V, torch.zeros_like(grad_V))
            Vs = torch.where(finite, V, torch.ones_like(V))

            # project the loss's sensitivity to V onto the Newton-Raphson unknowns:
            # an angle unknown takes Im(conj(V) . gV), a magnitude one
            # Re(conj(V/|V|) . gV) -- the two halves of dL/dV in polar coordinates
            xbar = torch.zeros(n_scen, dim_J, dtype=torch.float64)
            th_bus = np.flatnonzero(theta_col >= 0)
            vm_bus = np.flatnonzero(vm_col >= 0)
            if th_bus.size:
                xbar[:, torch.as_tensor(theta_col[th_bus])] = \
                    (torch.conj(Vs[:, th_bus]) * gV[:, th_bus]).imag
            if vm_bus.size:
                v_hat = Vs[:, vm_bus] / Vs[:, vm_bus].abs()
                xbar[:, torch.as_tensor(vm_col[vm_bus])] = (torch.conj(v_hat) * gV[:, vm_bus]).real

            lam = torch.as_tensor(sweep.solve_JT(np.ascontiguousarray(xbar.numpy())))

            # lambda is the gradient with respect to the per-unit bus injection; an
            # element's own is that, times how it enters its bus
            def element_grad(bus_of_element, row_of_bus, sign):
                grad = torch.zeros(n_scen, bus_of_element.shape[0], dtype=torch.float64)
                rows = np.where(bus_of_element >= 0, row_of_bus[np.clip(bus_of_element, 0, None)], -1)
                sel = np.flatnonzero(rows >= 0)   # an element off the grid, or on a bus
                if sel.size:                      # with no such equation, keeps its zero
                    grad[:, torch.as_tensor(sel)] = sign * lam[:, torch.as_tensor(rows[sel])] / pf.sn_mva
                return grad

            wants_lp, wants_lq, wants_gp, wants_sp = ctx.needs
            grad_load_p = element_grad(pf._load_bus, p_row, -1.) if wants_lp else None
            grad_load_q = element_grad(pf._load_bus, q_row, -1.) if wants_lq else None
            grad_gen_p = element_grad(pf._gen_bus, p_row, +1.) if wants_gp else None
            grad_sgen_p = element_grad(pf._sgen_bus, p_row, +1.) if wants_sp else None
            return grad_load_p, grad_load_q, grad_gen_p, grad_sgen_p, None

    return _Op
