# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["ScenarioSweepCPP", "LimitViolation", "ViolationElementType",
           "LimitViolationType", "PreContingencyResult",
           "ContingencyResult", "SecurityAnalysisResult"]

import numpy as np
import warnings


try:
    from lightsim2grid.lightSimBackend import LightSimBackend
    __all__.append("ScenarioSweep")
    GRID2OP_INSTALLED = True
except ImportError as exc_:  # noqa: F841
    # grid2Op is not installed
    GRID2OP_INSTALLED = False

from lightsim2grid.algorithm import AlgorithmType
from .lightsim2grid_cpp import ScenarioSweepCPP
# reused as-is, not duplicated: both ContingencyAnalysis.run() and ScenarioSweep.run()
# (below) return these same dataclasses / the same pybind11-bound LimitViolation type,
# so a caller handling one can handle the other identically. Defined unconditionally
# (not gated on grid2op being installed) in contingencyAnalysis.py, imported before
# this module in lightsim2grid/__init__.py -- no circular-import risk.
from .contingencyAnalysis import (PreContingencyResult, ContingencyResult, SecurityAnalysisResult,
                                   LimitViolation, ViolationElementType, LimitViolationType)


class ScenarioSweep:
    """
    Batch powerflow that varies **both** the injection and a contingency (line / trafo
    disconnection) per simulation, independently, row-aligned: row `i` of every
    ``modify_*`` input is solved together with row `i` of ``set_contingency_lines`` /
    ``set_contingency_trafos``.

    Unlike :class:`lightsim2grid.timeSerie.TimeSerie` /
    :class:`lightsim2grid.injectionSweep.InjectionSweep` (a single bundled
    ``compute_V_from_inj`` call), this class uses a setter-based API: call as many of
    ``modify_gen_p`` / ``modify_sgen_p`` / ``modify_load_p`` / ``modify_load_q`` /
    ``modify_gen_v`` / ``set_contingency_lines`` / ``set_contingency_trafos`` as are
    relevant, then ``compute()``. Any axis you never set defaults to the grid's own
    state for every row (its own target injection / voltage setpoint, or "nothing
    disconnected"). The first setter you call fixes the number of simulations; every
    later setter is checked against it immediately. Note ``modify_gen_v`` is different
    from the other ``modify_*`` setters: it does not feed the injection (Sbus), it only
    re-seeds the voltage magnitude (in pu, NOT kV) at each voltage-regulating generator's
    regulated bus before that row's solve.

    This is deliberately a *different* API from
    :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis`'s ``add_n1`` /
    ``add_nk`` (a set of distinct scenarios applied to one shared base case): here
    each row pairs its own injection with its own contingency, so a dense,
    row-aligned mask is the natural shape -- the two APIs are not interchangeable.

    Examples
    ---------

    .. code-block:: python

        import numpy as np
        import grid2op
        from lightsim2grid import ScenarioSweep
        from lightsim2grid import LightSimBackend

        env_name = ...
        env = grid2op.make(env_name, backend=LightSimBackend())

        sweep = ScenarioSweep(env)
        sweep.modify_load_p(load_p)  # shape (n_simul, n_load)
        line_mask = np.zeros((load_p.shape[0], env.n_line), dtype=bool)
        line_mask[:, 3] = True  # disconnect line 3 for every simulation
        sweep.set_contingency_lines(line_mask)
        sweep.compute()
        Vs = sweep.get_voltages()
        Ps, amps = sweep.compute_P(), sweep.compute_A()
    """

    #: the c++ class this wrapper drives.
    _CPP_CLASS = ScenarioSweepCPP

    def __init__(self, grid2op_env):
        if not GRID2OP_INSTALLED:
            raise RuntimeError(f"Impossible to use the python wrapper `{type(self).__name__}` "
                               "when grid2op is not installed. Please fall back to the "
                               "c++ version (available in python) with:\n"
                               f"\tfrom lightsim2grid.scenarioSweep import {type(self)._CPP_CLASS.__name__}\n"
                               "and refer to the appropriate documentation.")

        from grid2op.Environment import Environment  # type: ignore # otherwise i got issues...
        if not isinstance(grid2op_env.backend, LightSimBackend):
            raise RuntimeError("This class only works with LightSimBackend")
        if not isinstance(grid2op_env, Environment):
            raise RuntimeError("Please an environment of class \"Environment\", "
                               "and not \"MultimixEnv\" or \"BaseMultiProcessEnv\"")
        self.grid2op_env = grid2op_env.copy()
        self.computer = type(self)._CPP_CLASS(self.grid2op_env.backend._grid)
        self.__computed = False
        # ScenarioSweepCPP.set_contingency_lines/trafos are write-only (no getter) --
        # run() needs each row's disconnected branch ids for
        # ContingencyResult.element_ids/element_names, so this wrapper caches the last
        # validated mask it sent down. Reset by clear() below.
        self._line_mask = None
        self._trafo_mask = None

        self.available_default_algorithms = self.computer.available_default_algorithms()
        if AlgorithmType.NR_KLU in self.available_default_algorithms:
            # use the faster KLU if available
            self.computer.change_algorithm(AlgorithmType.NR_KLU)

    @property
    def nb_thread(self):
        """Number of OS threads used to compute the steps (default: ``1``).

        The steps are split into contiguous ranges, each solved by its own thread with its
        own solver (and its own admittance-matrix copy, since a contingency edits it), writing
        to disjoint rows of the result matrix: the results do **not** depend on the number of
        threads. Values ``< 1`` are clamped to ``1``.

        Must be set before ``compute()`` actually runs; it has no effect on a batch that has
        already been computed.
        """
        return self.computer.nb_thread

    @nb_thread.setter
    def nb_thread(self, val: int):
        if int(val) != val:
            raise ValueError("The `nb_thread` attribute must be an integer.")
        self.computer.nb_thread = int(val)

    @property
    def init_from_n_powerflow(self):
        """Whether to initialize the complex voltages of **each** simulation with the results
        of a "n" powerflow (a powerflow with no injection change and no contingency) instead of
        the vector given to ``compute``. Default: ``False``. Must be set before ``compute()``
        actually runs.
        """
        return self.computer.init_from_n_powerflow

    @init_from_n_powerflow.setter
    def init_from_n_powerflow(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `init_from_n_powerflow` attribute must be a boolean.")
        self.computer.init_from_n_powerflow = bool(val)

    @property
    def handle_disconnected_grid(self):
        """Whether a row whose contingency splits the grid into several connected
        components is simulated on its largest component instead of being skipped.
        Default: ``False``, meaning such a row is not simulated at all (its voltages
        are left at 0.). When ``True``, the buses of the other component(s) are
        "masked" (their voltage reported as 0.) and the largest component is solved
        normally, without triggering any extra matrix re-factorization. Supported by
        the Newton-Raphson family (AC) and by the DC solver; a non Newton-Raphson AC
        algorithm (*eg* Gauss-Seidel or Fast-Decoupled) is rejected. Same
        name/semantics as :attr:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.handle_disconnected_grid`.
        """
        return self.computer.handle_disconnected_grid

    @handle_disconnected_grid.setter
    def handle_disconnected_grid(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `handle_disconnected_grid` attribute must be a boolean.")
        self.computer.handle_disconnected_grid = bool(val)

    @property
    def compute_limit_violations(self):
        """Whether limit violations are computed inline, per row, during ``compute()``
        / ``run()`` (see :func:`get_violations` / :func:`get_violations_n` /
        :func:`run` -- there is no ``converged`` / ``converged_n`` here, unlike
        :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis`: a
        non-converged row's violations already carry a
        ``ViolationElementType.GRID``-typed sentinel that fully encodes that row's
        status by itself). Default: ``False``. Computing violations means an extra
        per-element current / voltage check in every row's solve, so leave this off
        if you only need :func:`compute_flows` / :func:`get_flows`. Changing this
        flag clears any previously-computed results.
        """
        return self.computer.compute_limit_violations

    @compute_limit_violations.setter
    def compute_limit_violations(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `compute_limit_violations` attribute must be a boolean.")
        val = bool(val)
        if val == self.computer.compute_limit_violations:
            return  # no-op, matches the C++ side (which also no-ops and does not clear)
        self.computer.compute_limit_violations = val  # this clear()s the whole C++ object:
        # registered injections, contingency masks, handle_disconnected_grid, the
        # row-count lock, everything -- not just computed results. Drop the Python-side
        # mask cache too, or a stale non-None _line_mask/_trafo_mask survives this reset
        # and the next run() derives element_ids/element_names from masks that no longer
        # match what was actually registered C++-side.
        self.__computed = False
        self._line_mask = None
        self._trafo_mask = None

    @property
    def violation_threshold(self):
        """Threshold (a ``float`` in ``]0., 1.]``, default ``1.0``) applied to every
        limit-violation check performed when :attr:`compute_limit_violations` is
        ``True``. Same meaning as
        :attr:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.violation_threshold`
        -- see that docstring for the full formulas. Lowering it makes every check
        stricter (more violations reported, never fewer) and invalidates any
        already-computed results; raising it back up does not.
        """
        return self.computer.violation_threshold

    @violation_threshold.setter
    def violation_threshold(self, val):
        try:
            val = float(val)
        except (TypeError, ValueError):
            raise ValueError("The `violation_threshold` attribute must be a real number.")
        if not (0. < val <= 1.):
            raise ValueError("The `violation_threshold` attribute must be in the range "
                             f"]0., 1.] (got {val}).")
        if val < self.computer.violation_threshold:
            # mirrors the C++ side, which calls clear_results_only() when the
            # threshold is tightened (results computed under the previous, looser
            # threshold would silently under-report)
            self.__computed = False
        self.computer.violation_threshold = val

    def _check_2d(self, arr, name):
        arr = np.asarray(arr)
        if arr.ndim != 2:
            raise RuntimeError(f"{name} should be a matrix with rows representing simulations "
                               f"and columns representing individual elements.")
        return arr

    def modify_gen_p(self, gen_p):
        """Per-step active generator setpoints, shape ``(n_simul, n_gen)``."""
        gen_p = self._check_2d(gen_p, "gen_p")
        if gen_p.shape[1] != self.grid2op_env.n_gen:
            raise RuntimeError(f"The number of generators on the grid ({self.grid2op_env.n_gen}) "
                               f"differs from the number of columns of gen_p ({gen_p.shape[1]}).")
        self.computer.modify_gen_p(gen_p)
        self.__computed = False

    def modify_sgen_p(self, sgen_p):
        """Per-step active static generator setpoints, shape ``(n_simul, n_sgen)``."""
        sgen_p = self._check_2d(sgen_p, "sgen_p")
        n_sgen = len(self.grid2op_env.backend._grid.get_static_generators())
        if sgen_p.shape[1] != n_sgen:
            raise RuntimeError(f"The number of static generators on the grid ({n_sgen}) "
                               f"differs from the number of columns of sgen_p ({sgen_p.shape[1]}).")
        self.computer.modify_sgen_p(sgen_p)
        self.__computed = False

    def modify_load_p(self, load_p):
        """Per-step active load setpoints, shape ``(n_simul, n_load)``."""
        load_p = self._check_2d(load_p, "load_p")
        if load_p.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid ({self.grid2op_env.n_load}) "
                               f"differs from the number of columns of load_p ({load_p.shape[1]}).")
        self.computer.modify_load_p(load_p)
        self.__computed = False

    def modify_load_q(self, load_q):
        """Per-step reactive load setpoints, shape ``(n_simul, n_load)``."""
        load_q = self._check_2d(load_q, "load_q")
        if load_q.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid ({self.grid2op_env.n_load}) "
                               f"differs from the number of columns of load_q ({load_q.shape[1]}).")
        self.computer.modify_load_q(load_q)
        self.__computed = False

    def modify_gen_v(self, gen_v):
        """Per-step generator target voltage magnitude, shape ``(n_simul, n_gen)``, in pu
        (``vm_pu``), NOT kV. Unlike ``modify_gen_p``/``modify_load_p``/``modify_load_q``,
        this does NOT feed the injection (Sbus) -- it only re-seeds ``|V|`` at each
        voltage-regulating generator's regulated bus before that step's solve."""
        gen_v = self._check_2d(gen_v, "gen_v")
        if gen_v.shape[1] != self.grid2op_env.n_gen:
            raise RuntimeError(f"The number of generators on the grid ({self.grid2op_env.n_gen}) "
                               f"differs from the number of columns of gen_v ({gen_v.shape[1]}).")
        self.computer.modify_gen_v(gen_v)
        self.__computed = False

    def set_contingency_lines(self, mask):
        """Per-step powerline contingency mask, shape ``(n_simul, n_line)``, dtype bool.
        ``True`` means "deactivate this powerline for this simulation"."""
        mask = self._check_2d(mask, "mask").astype(dtype=bool, copy=False)
        n_line = len(self.grid2op_env.backend._grid.get_lines())
        if mask.shape[1] != n_line:
            raise RuntimeError(f"The number of powerlines on the grid ({n_line}) "
                               f"differs from the number of columns of the mask ({mask.shape[1]}).")
        self.computer.set_contingency_lines(mask)
        self._line_mask = mask  # cached: see the note in __init__ (no C++-side getter)
        self.__computed = False

    def set_contingency_trafos(self, mask):
        """Per-step trafo contingency mask, shape ``(n_simul, n_trafo)``, dtype bool.
        See :func:`set_contingency_lines`."""
        mask = self._check_2d(mask, "mask").astype(dtype=bool, copy=False)
        n_trafo = len(self.grid2op_env.backend._grid.get_trafos())
        if mask.shape[1] != n_trafo:
            raise RuntimeError(f"The number of trafos on the grid ({n_trafo}) "
                               f"differs from the number of columns of the mask ({mask.shape[1]}).")
        self.computer.set_contingency_trafos(mask)
        self._trafo_mask = mask  # cached: see the note in __init__ (no C++-side getter)
        self.__computed = False

    def compute(self, v_init=None, max_iter=None, tol=None, ignore_errors=False):
        """
        Run the batch: one powerflow per simulation, using whatever was set by the
        ``modify_*`` / ``set_contingency_*`` setters. Raises if nothing was ever set.

        ``max_iter`` / ``tol`` default to the backend's own values
        (``self.grid2op_env.backend.max_it`` / ``.tol``) when not given.
        """
        if v_init is None:
            v_init_comp = self.grid2op_env.backend.V
        else:
            v_init_comp = 1.0 * v_init  # make a copy !
        if max_iter is None:
            max_iter = self.grid2op_env.backend.max_it
        if tol is None:
            tol = self.grid2op_env.backend.tol
        self.computer.compute(v_init_comp, max_iter, tol)
        status = self.computer.get_status()
        if status != 1 and not ignore_errors:
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        elif status != 1:
            warnings.warn(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        self.__computed = True

    def get_voltages(self):
        """Complex voltage, at each bus, for each simulation. Must be called after ``compute()``."""
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        return self.computer.get_voltages().copy()

    def compute_A(self):
        """
        Current flows (in Amps, A) at the origin (for powerline) / high voltage (for
        transformer) side, per simulation. Does not recompute the voltages; you must
        call ``compute()`` first.
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        ampss = self.computer.compute_flows()
        return 1000. * ampss

    def compute_P(self):
        """
        Active power flows (in MW) at the origin (for powerline) / high voltage (for
        transformer) side, per simulation. Does not recompute the voltages; you must
        call ``compute()`` first.
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        mws = 1. * self.computer.compute_power_flows()  # If I don't copy, lazy eval may break stuff...
        return mws

    def get_flows(self, v_init=None, max_iter=None, tol=None, ignore_errors=False):
        """
        Run ``compute()`` (using whatever was set by the setters) then retrieve the
        active power / current / voltage results for every simulation, in one call.

        Each row of the resulting flow matrix corresponds to a simulation.
        """
        self.compute(v_init, max_iter, tol, ignore_errors)
        amps = self.compute_A()
        Ps = self.compute_P()
        Vs = self.get_voltages()
        return Ps, amps, Vs

    def get_violations(self):
        """Per row (same order as every ``modify_*`` / ``set_contingency_*`` input):
        list of :class:`LimitViolation`. Requires :attr:`compute_limit_violations` to
        be ``True`` (raises otherwise). Prefer :func:`run` for a structured result --
        this is the raw passthrough. See
        :attr:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.compute_limit_violations`'s
        note on the sentinel entry a non-converged row carries (there is no separate
        ``converged`` here, by design -- see :attr:`compute_limit_violations`).
        """
        return self.computer.get_violations()

    def get_violations_n(self):
        """List of :class:`LimitViolation` for the pre-batch ("n") case (no injection
        change, no contingency) shared by every row. Requires
        :attr:`compute_limit_violations` to be ``True`` (raises otherwise)."""
        return self.computer.get_violations_n()

    @staticmethod
    def _row_converged(limit_violations) -> bool:
        """No ``converged()`` / ``converged_n()`` on ScenarioSweep by design (see
        :attr:`compute_limit_violations`) -- recovers the same boolean from
        :func:`get_violations` / :func:`get_violations_n` alone: a non-converged row
        always carries exactly one ``ViolationElementType.GRID``-typed entry."""
        return not any(v.element_type == ViolationElementType.GRID for v in limit_violations)

    def run(self) -> SecurityAnalysisResult:
        """
        Run this batch (calling :func:`compute` if not already done) and report, for
        the pre-batch ("n") case and for each row, the list of limit violations --
        same return type as
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.run`, reusing
        the very same :class:`PreContingencyResult` / :class:`ContingencyResult` /
        :class:`SecurityAnalysisResult` dataclasses so callers can handle either
        result the same way.

        Requires :attr:`compute_limit_violations` to be ``True`` (set via
        ``this_instance.compute_limit_violations = True`` -- there is no constructor
        argument for it, unlike :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis`),
        else a ``RuntimeError`` is raised.

        Unlike :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.run`,
        rows are already in caller-set order here (no dedup / reordering concept,
        since ``ScenarioSweep`` rows are independent scenarios, not a *set* of
        contingencies) -- and every ``ContingencyResult.contingency_name`` is
        ``None`` (no such concept on ``ScenarioSweep``); ``element_ids`` /
        ``element_names`` are instead derived from that row's own
        ``set_contingency_lines`` / ``set_contingency_trafos`` mask.
        """
        if not self.computer.compute_limit_violations:
            raise RuntimeError("`run` requires `compute_limit_violations=True`, set via "
                               "`this_instance.compute_limit_violations = True` before `compute()`.")
        if not self.__computed:
            # a per-row failure (a contingency that islands the grid, a diverging row,
            # ...) is exactly what this method's own sentinel-violation reporting (see
            # get_violations()/_row_converged above) exists to surface -- do not let
            # compute()'s default RuntimeError pre-empt that. Unlike
            # ContingencyAnalysis (whose C++ side throws unconditionally on a
            # diverging pre-batch "n" case, since every contingency is meaningless
            # relative to it), a diverging "n" case here only sets status=0 like any
            # other row failure -- ignore_errors=True lets it through too, and it is
            # still faithfully reported: the n-divergence fix stamps a GRID/DIVERGENCE
            # sentinel into get_violations_n() rather than leaving it an empty list.
            self.compute(ignore_errors=True)

        violations_n = list(self.computer.get_violations_n())
        pre_contingency_result = PreContingencyResult(
            converged=self._row_converged(violations_n),
            limit_violations=violations_n,
        )

        n_line = len(self.grid2op_env.backend._grid.get_lines())
        violations = self.computer.get_violations()
        nb_steps = len(violations)
        line_mask = self._line_mask if self._line_mask is not None else np.zeros((nb_steps, n_line), dtype=bool)
        trafo_mask = self._trafo_mask if self._trafo_mask is not None else np.zeros(
            (nb_steps, len(self.grid2op_env.backend._grid.get_trafos())), dtype=bool)

        post_contingency_results = []
        for row_id in range(nb_steps):
            row_violations = list(violations[row_id])
            element_ids = [int(el_id) for el_id in np.nonzero(line_mask[row_id])[0]]
            element_ids += [n_line + int(el_id) for el_id in np.nonzero(trafo_mask[row_id])[0]]
            post_contingency_results.append(ContingencyResult(
                element_ids=element_ids,
                element_names=[str(self.grid2op_env.name_line[el_id]) for el_id in element_ids],
                contingency_name=None,
                converged=self._row_converged(row_violations),
                limit_violations=row_violations,
            ))
        return SecurityAnalysisResult(pre_contingency_result, post_contingency_results)

    def clear(self):
        """Clear everything, as if nothing had ever been set / computed."""
        self.computer.clear()
        self.__computed = False
        self._line_mask = None
        self._trafo_mask = None

    def close(self):
        """permanently close the object"""
        self.grid2op_env.close()
        self.clear()
        self.computer.close()
