Scenario Sweep
=======================================

Goal
--------------------------

This class varies **both** the injection (generators / static generators / loads) and a
contingency (a line / transformer disconnection) per simulation, independently,
row-aligned: row `i` of every injection input is solved together with row `i` of the
contingency mask.

It is the 4th instantiation of the same underlying C++ template as
:class:`lightsim2grid.timeSerie.TimeSerie` / :class:`lightsim2grid.injectionSweep.InjectionSweep`
(which vary only the injection) and :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis`
(which varies only the topology, over one shared base injection) -- see :doc:`time_series`
and :doc:`security_analysis`.

It can be used as:

.. code-block:: python

    import numpy as np
    import grid2op
    from lightsim2grid import ScenarioSweep
    from lightsim2grid import LightSimBackend

    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())
    obs = env.reset()

    n_simul = 10
    load_p = np.tile(obs.load_p, (n_simul, 1))
    load_p *= 1 + 0.01 * np.arange(n_simul).reshape(-1, 1)  # vary the load a bit per row

    # lightsim2grid's own line / trafo counts -- NOT grid2op's flat `env.n_line`
    # (grid2op merges lines and trafos into one "powerline" list; lightsim2grid keeps
    # them as two separate containers, and the contingency masks below follow that split)
    n_line = len(env.backend._grid.get_lines())
    n_trafo = len(env.backend._grid.get_trafos())
    line_mask = np.zeros((n_simul, n_line), dtype=bool)
    line_mask[5:, 3] = True  # disconnect line 3 for the last 5 simulations

    sweep = ScenarioSweep(env)
    sweep.modify_load_p(load_p)
    sweep.set_contingency_lines(line_mask)
    # modify_gen_p / modify_sgen_p / modify_load_q / set_contingency_trafos were never
    # called: those axes default to the grid's own current state, broadcast across
    # every row (see "Unset axes" below).
    sweep.compute()

    res_v = sweep.get_voltages()
    res_a = sweep.compute_A()
    res_p = sweep.compute_P()

    # res_p[row_id] / res_a[row_id] / res_v[row_id] are the flows / voltages for
    # simulation `row_id`, combining that row's injections and that row's contingency.

The setter-based API
----------------------------

Unlike `TimeSerie` / `InjectionSweep`'s single bundled `compute_V_from_inj` call,
`ScenarioSweep` is built up incrementally: call as many of `modify_gen_p` /
`modify_sgen_p` / `modify_load_p` / `modify_load_q` / `modify_gen_v` /
`set_contingency_lines` / `set_contingency_trafos` as are relevant, then `compute()`.

- **Row-count locking.** The first setter you call fixes the number of simulations for
  the whole object; every later setter is checked against it immediately -- a shape
  mistake raises right at the call that made it, not later inside `compute()`.
- **Unset axes.** Any injection axis you never call (eg `modify_gen_p`) defaults to the
  grid's own current target value, broadcast across every row -- not zero. Any
  contingency mask you never call (`set_contingency_lines` / `set_contingency_trafos`)
  defaults to all-`False` -- nothing disconnected. `compute()` raises if you never call
  *any* setter (a degenerate call -- use `ac_pf()` / `dc_pf()` directly for a single
  powerflow).
- **`modify_gen_v` is different from the other four `modify_*` setters.** It does not
  feed the injection (`Sbus`) at all -- a PV bus's voltage magnitude is never part of
  what Newton-Raphson solves for, only of what it starts from. `modify_gen_v` instead
  re-seeds `|V|` at each voltage-regulating generator's regulated bus immediately before
  that row's solve, shape `(n_simul, n_gen)`, in pu (`vm_pu`), NOT kV. Left unset, every
  row keeps using the grid's own `target_vm_pu`, exactly as if it had never been called.
- `TimeSerie` / `InjectionSweep` also gained this same setter API (`modify_*` +
  `compute()`) alongside their legacy `compute_V_from_inj` call, which remains
  available (deprecated, not removed).

Two contingency APIs, on purpose
------------------------------------

`ScenarioSweep.set_contingency_lines` / `set_contingency_trafos` and
`ContingencyAnalysis.add_n1` / `add_nk` are **deliberately two different APIs**, not a
naming inconsistency:

- `ContingencyAnalysis` applies an arbitrary *set* of distinct contingencies to one
  shared base injection -- `add_n1(3)` registers "disconnect line 3", independent of
  any row ordering or count.
- `ScenarioSweep` pairs each contingency with that same row's own injection, so a
  dense, row-aligned boolean mask (shape `(n_simul, n_line)` / `(n_simul, n_trafo)`,
  `True` = "deactivate this branch for this simulation") is the natural fit instead.

`ScenarioSweep` does not have `add_n1` / `add_nk`; `ContingencyAnalysis` does not have
`set_contingency_lines` / `set_contingency_trafos`.

.. note::

    Dense contingency masks cost real memory at scale (eg 10k simulations x 3k branches
    ~= 30 MB) -- fine for typical use, a lighter-weight representation may be added in a
    future release if that becomes a bottleneck.

Handling disconnected grids and limit violations
------------------------------------------------------

`ScenarioSweep` has the same `handle_disconnected_grid` mode and inline limit-violation
checking (`compute_limit_violations` / `violation_threshold` / `get_violations` /
`get_violations_n`) as `ContingencyAnalysis` -- see :doc:`security_analysis` for the full
description of what each does. Same names, same semantics, and both classes'
`get_violations` / `get_violations_n` return the same `LimitViolation` objects.

.. code-block:: python

    sweep = ScenarioSweep(env)
    sweep.compute_limit_violations = True  # see warning below: set this FIRST
    sweep.handle_disconnected_grid = True
    sweep.modify_load_p(load_p)
    sweep.set_contingency_lines(line_mask)
    sweep.compute()

    row_violations = sweep.get_violations()  # one list of LimitViolation per row
    n_violations = sweep.get_violations_n()  # violations of the shared pre-batch base case

    sa_res = sweep.run()  # PreContingencyResult / ContingencyResult / SecurityAnalysisResult,
                           # reusing lightsim2grid.contingencyAnalysis's dataclasses --
                           # contingency_name is always None here (no such concept on ScenarioSweep)

.. warning::

    Just like on `ContingencyAnalysis`, setting `compute_limit_violations` through the
    property **clears the whole object** -- registered injections, contingency masks,
    `handle_disconnected_grid`, everything, the same reset `clear()` performs. Set it
    *first*, before calling any `modify_*` / `set_contingency_*` setter or
    `handle_disconnected_grid`, or those calls will be silently wiped.

    Unlike `ContingencyAnalysis`, `ScenarioSweep` has no `(grid, compute_limit_violations)`
    constructor shortcut to dodge this ordering -- the property is the only way to set it.

.. note::

    There is deliberately **no** `converged` / `converged_n` on `ScenarioSweep` (unlike
    `ContingencyAnalysis`). A non-converged row's `get_violations()` entry already carries
    a single `LimitViolation` with `element_type == ViolationElementType.GRID` and
    `violation_type` either `NOT_SIMULATED` (skipped before the solver ran, eg it splits the
    grid and `handle_disconnected_grid` is off) or `DIVERGENCE` (the solver ran but did not
    converge) -- a separate convergence flag would be redundant. If the shared pre-batch "n"
    powerflow itself diverges, that same sentinel is stamped into `get_violations_n()` and
    into every row of `get_violations()`, rather than leaving them as empty lists
    indistinguishable from "converged, nothing found".

A future release may add a "combine mode" axis (this row-aligned / zipped behavior vs. a
cartesian "every contingency x every injection profile") -- not yet available.

Detailed usage
--------------------------

.. automodule:: lightsim2grid.scenarioSweep
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
