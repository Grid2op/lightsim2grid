# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Continuation powerflow (CPF): trace the PV curve of a grid up to its voltage-collapse
"nose" point.

Everything numerical happens in C++ (:class:`lightsim2grid.lightsim2grid_cpp.ContinuationSweepCPP`,
see ``src/core/batch_algorithm/ContinuationSweep.hpp``). This module only turns a
loading factor and a pair of steering vectors into the *target injection state* the C++
class continues towards, and presents the traced curve back.
"""

__all__ = ["ContinuationPowerFlow", "CPFResult", "run_cpf", "plot_pv_curve"]

import numpy as np

from .lightsim2grid_cpp import ContinuationSweepCPP


class CPFResult:
    """
    One traced curve. Attributes:

    lam: ``(n_points,)``
        The continuation parameter at each point. ``lam[0]`` is 0 (the base case) and
        ``lam == 1`` is the target state, i.e. the full requested increase. So the load
        of load ``l`` at point ``i`` is ``load_p[l] * (1 + alpha[l] * (k - 1) * lam[i])``
        with ``k`` the loading factor and ``alpha`` the steering vector.
    V: ``(n_points, n_bus)`` complex
        Bus voltages at each point, in the grid's own bus ordering. Buses that are not
        part of the solved system read exactly ``0``.
    Vm, Va_deg: ``(n_points, n_bus)``
        Magnitude (pu) and angle (degrees) of ``V``.
    tangent_lam: ``(n_points,)``
        The tangent's lambda component at each point, in ``(0, 1]``. It tends to zero as
        the curve turns: that is the collapse indicator. The last point has no tangent
        computed after it and reads 0.
    lam_max: float
        The largest lambda reached, i.e. the loading margin along the chosen direction.
    success: bool
        Whether the curve reached its requested end (the nose, or ``stop_at_lam``).
    msg: str
        Why the run stopped, in words.
    nb_retries: int
        How many times a corrector failed and the step had to be halved.
    """

    def __init__(self, lam, V, tangent_lam, success, msg, nb_retries, cpp):
        self.lam = lam
        self.V = V
        self.Vm = np.abs(V)
        self.Va_deg = np.rad2deg(np.angle(V))
        self.tangent_lam = tangent_lam
        self.success = success
        self.msg = msg
        self.nb_retries = nb_retries
        self.lam_max = float(lam[-1]) if lam.size else 0.0
        #: the underlying C++ object, for the timers, the flows and the linear-solver counters
        self.cpp = cpp

    def __repr__(self):
        return (f"CPFResult(success={self.success}, n_points={self.lam.size}, "
                f"lam_max={self.lam_max:.4f}, msg={self.msg!r})")


def _checked_steering(vect, nb_el, name):
    """A steering vector: one coefficient per element, finite, in [0, 1]."""
    arr = np.asarray(vect, dtype=float).ravel()
    if arr.shape[0] != nb_el:
        raise ValueError(f"{name}: expected one coefficient per element ({nb_el}), "
                         f"got {arr.shape[0]}.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name}: every coefficient must be finite.")
    if np.any(arr < 0.0) or np.any(arr > 1.0):
        raise ValueError(f"{name}: every coefficient must be in [0, 1] "
                         f"(got min {arr.min()}, max {arr.max()}). A coefficient of 0 "
                         f"keeps that element at its base value; 1 moves it fully.")
    return arr


class ContinuationPowerFlow:
    """
    Continuation powerflow on a :class:`lightsim2grid.network.LSGrid`.

    The curve runs from the grid's own injection state (``lambda = 0``) to a *target*
    state (``lambda = 1``), following the classical predictor / corrector scheme. This
    is MATPOWER's ``runcpf`` formulation, and the option names and defaults below follow
    its ``cpf.*`` options so that nothing here is a new concept.

    The target is built for you from a loading factor and two steering vectors:

    - ``load_steering`` -- one coefficient per load, in ``[0, 1]``. A load with 0 does
      not move during the continuation; a load with 1 moves fully. Defaults to all ones
      (every load scaled).
    - ``gen_steering`` -- the same, per generator. Defaults to 1 for every non-slack
      generator and 0 for the slack ones, so generation follows the load and the slack
      only picks up the incremental losses -- MATPOWER requires exactly that of a target
      case ("same as base case w.r.t. Qg and slack Pg"). Pass ``0`` to hold generation
      fixed and let the slack supply the whole increase instead; that traces a different
      curve, with a different nose.

    .. warning::
        This is a *natural* parameterisation: the corrector solves at a fixed lambda, so
        the curve stops AT the nose and does not trace the lower branch. That needs an
        arc-length parameterisation (MATPOWER's ``cpf.parameterization`` 2 / 3), which is
        not implemented yet; when it is, it will become the default here as it is there.

    Examples
    ---------

    .. code-block:: python

        import numpy as np
        from lightsim2grid.network import init_from_pandapower
        from lightsim2grid.continuationPowerflow import ContinuationPowerFlow
        import pandapower.networks as pn

        grid = init_from_pandapower(pn.case118())

        # how much margin is there if every load grows together?
        cpf = ContinuationPowerFlow(grid)
        res = cpf.run(loading_factor=3.0)
        print(res.lam_max)

        # ... and if only the loads of one area grow?
        steering = np.zeros(len(grid.get_loads()))
        steering[my_area_load_ids] = 1.0
        res = cpf.run(loading_factor=3.0, load_steering=steering)

    """

    def __init__(self, grid, algorithm=None):
        self._grid = grid
        self._cpp = ContinuationSweepCPP(grid)
        if algorithm is not None:
            self._cpp.change_algorithm(algorithm)

        self._base_load_p = np.array(grid.get_load_target_p(), dtype=float)
        self._base_load_q = np.array([el.target_q_mvar for el in grid.get_loads()], dtype=float)
        self._base_gen_p = np.array(grid.get_gen_target_p(), dtype=float)
        self._base_sgen_p = np.array(grid.get_sgen_target_p(), dtype=float)
        self._gen_is_slack = np.array([el.is_slack for el in grid.get_generators()], dtype=bool)

    @property
    def cpp(self):
        """The underlying C++ object (timers, flows, linear-solver counters)."""
        return self._cpp

    def _build_target(self, loading_factor, load_steering, gen_steering, scale_q, direction):
        """
        Returns the four target axes, as ``{axis_name: values}``; an axis that does not
        move is left out entirely (the C++ side then treats it as "unchanged").
        """
        if direction is not None:
            if load_steering is not None or gen_steering is not None:
                raise ValueError("run(): pass either `direction` (explicit per-element deltas) "
                                 "or `load_steering` / `gen_steering` (a steered loading "
                                 "factor), not both -- they are two ways of saying the same "
                                 "thing and there is no sensible way to combine them.")
            known = {"load_p": self._base_load_p, "load_q": self._base_load_q,
                     "gen_p": self._base_gen_p, "sgen_p": self._base_sgen_p}
            unknown = set(direction) - set(known)
            if unknown:
                raise ValueError(f"run(): unknown key(s) {sorted(unknown)} in `direction`; "
                                 f"expected any of {sorted(known)}.")
            target = {}
            for axis, delta in direction.items():
                base = known[axis]
                delta = np.asarray(delta, dtype=float).ravel()
                if delta.shape[0] != base.shape[0]:
                    raise ValueError(f"run(): direction['{axis}'] has {delta.shape[0]} entries, "
                                     f"expected {base.shape[0]}.")
                target[axis] = base + delta
            return target

        if not np.isfinite(loading_factor) or loading_factor <= 1.0:
            raise ValueError(f"run(): loading_factor must be a finite number strictly greater "
                             f"than 1 (got {loading_factor}); it is the factor the steered "
                             f"elements reach at lambda = 1.")
        growth = loading_factor - 1.0

        if load_steering is None:
            alpha = np.ones(self._base_load_p.shape[0], dtype=float)
        else:
            alpha = _checked_steering(load_steering, self._base_load_p.shape[0], "load_steering")

        if gen_steering is None:
            # MATPOWER-like: generation follows the load, the slack only takes the losses.
            beta = (~self._gen_is_slack).astype(float)
        elif np.isscalar(gen_steering):
            beta = _checked_steering(np.full(self._base_gen_p.shape[0], float(gen_steering)),
                                     self._base_gen_p.shape[0], "gen_steering")
        else:
            beta = _checked_steering(gen_steering, self._base_gen_p.shape[0], "gen_steering")

        target = {"load_p": self._base_load_p * (1.0 + alpha * growth),
                  "gen_p": self._base_gen_p * (1.0 + beta * growth)}
        if scale_q:
            # constant power factor per load: P and Q scaled by the same coefficient
            target["load_q"] = self._base_load_q * (1.0 + alpha * growth)
        return target

    def run(self,
            loading_factor=2.0,
            load_steering=None,
            gen_steering=None,
            scale_q=True,
            direction=None,
            stop_at_lam=None,
            step=0.05,
            step_min=1e-4,
            step_max=0.2,
            adapt_step=False,
            adapt_step_damping=0.7,
            adapt_step_tol=1e-3,
            nose_tol=1e-5,
            max_steps=1000,
            exact_tangent=False,
            v_init=None,
            max_iter=10,
            tol=1e-8):
        """
        Trace the curve and return a :class:`CPFResult`.

        Parameters
        ----------
        loading_factor: ``float``
            What the fully-steered elements reach at ``lambda = 1``: 2.0 means "twice the
            base load". Must be > 1.
        load_steering: ``np.ndarray``, optional
            One coefficient per load, in ``[0, 1]``; see the class docstring. ``None``
            (default) means all ones.
        gen_steering: ``np.ndarray`` or ``float``, optional
            One coefficient per generator, in ``[0, 1]``, or a single number applied to
            every generator. ``None`` (default) means 1 for non-slack generators and 0
            for slack ones, so generation follows the load. Pass ``0`` to let the slack
            supply the whole increase.
        scale_q: ``bool``
            Scale each load's reactive power by the same coefficient as its active power
            (constant power factor, the default). ``False`` holds Q at its base value.
        direction: ``dict``, optional
            Escape hatch: explicit per-element DELTAS from the base state, in MW / MVAr,
            as ``{"load_p": ..., "load_q": ..., "gen_p": ..., "sgen_p": ...}`` (any
            subset). Mutually exclusive with the steering arguments.
        stop_at_lam: ``float``, optional
            Stop once lambda reaches this value (MATPOWER's numeric ``cpf.stop_at``).
            ``None`` (default) traces until the nose.
        step, step_min, step_max, adapt_step, adapt_step_damping, adapt_step_tol, nose_tol:
            The continuation's step control, named and defaulted as MATPOWER's ``cpf.*``.
        max_steps: ``int``
            Hard cap on the number of traced points.
        exact_tangent: ``bool``
            Refactorize the Jacobian at each converged point before taking its tangent.
            Off by default: the tangent then reuses the factorization the corrector left
            standing, which is one Newton iterate behind the converged point.
        v_init: ``np.ndarray``, optional
            Starting voltage, one entry per grid bus. Defaults to a flat 1.04 pu start.
        max_iter, tol:
            Passed to every corrector's Newton-Raphson.
        """
        target = self._build_target(loading_factor, load_steering, gen_steering,
                                    scale_q, direction)

        self._cpp.clear_target()
        if "load_p" in target:
            self._cpp.set_target_load_p(target["load_p"])
        if "load_q" in target:
            self._cpp.set_target_load_q(target["load_q"])
        if "gen_p" in target:
            self._cpp.set_target_gen_p(target["gen_p"])
        if "sgen_p" in target:
            self._cpp.set_target_sgen_p(target["sgen_p"])

        self._cpp.step = step
        self._cpp.step_min = step_min
        self._cpp.step_max = step_max
        self._cpp.adapt_step = adapt_step
        self._cpp.adapt_step_damping = adapt_step_damping
        self._cpp.adapt_step_tol = adapt_step_tol
        self._cpp.nose_tol = nose_tol
        self._cpp.max_steps = max_steps
        self._cpp.exact_tangent = exact_tangent
        self._cpp.stop_at_lam = -1.0 if stop_at_lam is None else float(stop_at_lam)

        if v_init is None:
            v_init = np.full(self._grid.total_bus(), 1.04, dtype=complex)
        self._cpp.compute(v_init, max_iter, tol)

        return CPFResult(lam=np.array(self._cpp.get_lam(), dtype=float, copy=True),
                         V=np.array(self._cpp.get_voltages(), dtype=complex, copy=True),
                         tangent_lam=np.array(self._cpp.get_tangent_lam(), dtype=float, copy=True),
                         success=bool(self._cpp.get_status() == 1),
                         msg=self._cpp.get_msg(),
                         nb_retries=self._cpp.nb_retries(),
                         cpp=self._cpp)


def run_cpf(grid, **kwargs):
    """
    One-shot helper: build a :class:`ContinuationPowerFlow` on ``grid`` and run it.
    Every keyword argument is forwarded to :func:`ContinuationPowerFlow.run`, except
    ``algorithm`` which selects the Newton-Raphson algorithm to use.

    .. code-block:: python

        from lightsim2grid.continuationPowerflow import run_cpf
        res = run_cpf(grid, loading_factor=3.0)
        print(f"collapse at {res.lam_max:.3f} of the requested increase")
    """
    algorithm = kwargs.pop("algorithm", None)
    return ContinuationPowerFlow(grid, algorithm=algorithm).run(**kwargs)


def plot_pv_curve(res, bus_indices=None, figsize=(10, 4)):
    """
    Plot ``|V|`` against lambda for the given buses, and the tangent's lambda component
    against the point index (which falls towards 0 at the nose).

    Needs matplotlib. ``bus_indices`` defaults to the ten buses whose voltage drops most
    along the curve -- the ones that actually say something about the collapse.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("plot_pv_curve needs matplotlib, which is not installed. "
                          "Install it with `pip install matplotlib`, or read "
                          "res.lam / res.Vm yourself.") from exc

    if res.lam.size == 0:
        raise ValueError("plot_pv_curve: this result holds no traced point (the base case "
                         f"did not converge). msg: {res.msg}")

    Vm = res.Vm
    if bus_indices is None:
        # buses that are actually in the solved system (the others read exactly 0)
        alive = np.flatnonzero(Vm[0] > 0.0)
        drop = Vm[0, alive] - Vm[-1, alive]
        bus_indices = alive[np.argsort(drop)[::-1][:10]]

    fig, (ax_pv, ax_t) = plt.subplots(1, 2, figsize=figsize)
    for bus in bus_indices:
        ax_pv.plot(res.lam, Vm[:, bus], marker=".", markersize=3, label=f"bus {bus}")
    ax_pv.set_xlabel(r"$\lambda$")
    ax_pv.set_ylabel(r"$|V|$ [pu]")
    ax_pv.set_title("PV curve")
    if len(bus_indices) <= 10:
        ax_pv.legend(fontsize="small")

    ax_t.plot(res.tangent_lam, marker=".", markersize=3)
    ax_t.axhline(0.0, linestyle="--", linewidth=1, color="k")
    ax_t.set_xlabel("point index")
    ax_t.set_ylabel(r"$t_\lambda$")
    ax_t.set_title(r"tangent $\lambda$ component ($\to 0$ at the nose)")

    fig.tight_layout()
    return fig
