Contingency Analysis
=======================================

Goal
-----------------

This class aims to make faster (and easier) the computations of a security analysis (which is the results of some 
powerflow after the disconnection of one or more powerlines)

This function is much (much) faster than its pure grid2op counterpart. For example,
on the case 118, to simulate all n-1 contingencies you can expect a **~20x** speed ups 
compared to using the grid2op `obs.simulate(..., time_step=0)` while obtaining the
exact same results (see section `Benchmarks`)

It can be used as:

.. code-block:: python

    import grid2op
    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid import LightSimBackend
    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    security_analysis = ContingencyAnalysis(env)
    security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
    res_p, res_a, res_v = security_analysis.get_flows()

    # in this results, then
    # res_p[row_id] will be the active power flows (origin side), on all powerlines corresponding to the `row_id` contingency.
    # res_a[row_id] will be the current flows, on all powerlines corresponding to step "row_id"
    # res_v[row_id] will be the complex voltage, on all bus of the grid corresponding to the `row_id` contingency.
    # you can retrieve which contingency is id'ed `row_id` with `security_analysis.contingency_order[row_id]`

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources.

.. note::

    A more advanced usage is given in the `examples\\contingency_analysis.py`
    file from the lightsim2grid package.

.. note::

    Set `security_analysis.init_from_n_powerflow = True` **before** the computation actually runs
    (eg before `get_flows` / `compute_V` / `run` are called -- setting it afterwards has no
    effect) to initialize each contingency with the voltage solution of the pre-contingency
    ("n") case instead of a flat start -- this is usually faster. See
    :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.init_from_n_powerflow`.

Handling contingencies that split the grid
-------------------------------------------

By default, a contingency that splits the grid into several connected components (an "islanding")
is **not simulated**: the corresponding row of the results is left at 0. (for the voltages) and
``NaN`` (for the flows). You can list which contingencies split the grid with
:func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.is_grid_connected_after_contingency`.

Starting from lightsim2grid 1.0.0, you can opt in to a mode that *does* simulate these
contingencies, on the largest connected component, by setting the ``handle_disconnected_grid``
attribute to ``True``:

.. code-block:: python

    import grid2op
    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid import LightSimBackend
    env = grid2op.make(..., backend=LightSimBackend())

    security_analysis = ContingencyAnalysis(env)
    security_analysis.add_all_n1_contingencies()

    # opt in: simulate the largest island instead of skipping split contingencies
    security_analysis.handle_disconnected_grid = True

    res_p, res_a, res_v = security_analysis.get_flows()

When this mode is enabled and a contingency splits the grid:

- the **largest** connected component is solved as a regular powerflow;
- the buses of the other component(s) are *masked*: their voltage is reported as ``0.`` (same
  convention as a skipped contingency) and they do not influence the solved component;
- if a slack generator ends up in a masked component, its slack weight is set to 0 and the
  remaining slack weights are rescaled so the slack power is shared only among the live slacks.

This is implemented **without re-triggering the symbolic factorization** of the linear solver
(the Jacobian, resp. the DC matrix, sparsity pattern is left unchanged), so it stays compatible
with the speed of the contingency analysis. To make it possible, the reference slack is chosen
once, before the computation, so as to minimise the number of contingencies that still have to be
skipped (those that would disconnect the chosen reference slack itself).

The mode works both in **AC** (with a Newton-Raphson algorithm) and in **DC**: in DC the masked
buses' rows of the reduced system are forced to the identity (so their angle is 0), the masked
injections are dropped and the slack imbalance is computed on the live component only. In both
cases the masked buses are reported as ``0.``.

.. note::

    This mode **requires a Newton-Raphson algorithm** (AC, the default) or the **DC** solver.
    Selecting an AC *non* Newton-Raphson algorithm (*eg* Gauss-Seidel or Fast-Decoupled) and
    enabling the mode raises an error.

Running the contingencies on multiple threads
---------------------------------------------

Starting from lightsim2grid 1.0.0, the contingencies can be solved on several CPU threads
at once. By default everything runs on a single thread (``nb_thread = 1``), reproducing the
exact same behaviour (and results) as before. Set the ``nb_thread`` attribute to a value
greater than ``1`` to split the work:

.. code-block:: python

    import grid2op
    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid import LightSimBackend
    env = grid2op.make(..., backend=LightSimBackend())

    security_analysis = ContingencyAnalysis(env)
    security_analysis.add_all_n1_contingencies()

    # solve the contingencies on 4 threads
    security_analysis.nb_thread = 4

    res_p, res_a, res_v = security_analysis.get_flows()

Internally the contingency list is split into ``nb_thread`` contiguous ranges, and each range
is solved by its own thread. To stay correct (and lock-free), every thread works on **its own**
solver instance and **its own** copy of the admittance matrix, and writes to a distinct set of
rows of the (shared) result matrix. As a consequence:

- the results do not depend on the number of threads: ``nb_thread`` changes the timing, not the
  numbers. They match the sequential results up to the solver's convergence tolerance (each
  thread keeps its own solver warm-start state, so the converged voltages agree to roughly
  ``1e-13``, far below the powerflow tolerance);
- it works for the AC (Newton-Raphson), DC and ``handle_disconnected_grid`` modes alike;
- there is a small per-thread set-up cost (one extra solver "warm-up" and one admittance-matrix
  copy per additional thread), so the speed-up is sub-linear and most useful when there are many
  contingencies to simulate.

This feature only relies on the C++ standard library (``std::thread``): no additional dependency
(MPI, OpenMP, ...) is required.

.. note::

    ``nb_thread`` is also available on the lower-level ``ContingencyAnalysisCPP`` class (same
    semantics).

Benchmarks (``nb_thread`` scaling)
+++++++++++++++++++++++++++++++++++

The script below scans ``nb_thread`` from 1 to 8 on a single, reasonably large real-topology
grid (``case6515rte``, up to 1001 n-1 contingencies), and is available by running, from the
root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 benchmark_ca_nb_threads.py

Results, made with:

- date: 2026-07-15 15:55 CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.0.2
- pandas version: 2.3.3
- pandapower version: 3.4.0
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.0.0.rc2
- lightsim2grid extra information:

	- klu_solver_available: True
	- nicslu_solver_available: False
	- cktso_solver_available: False
	- compiled_march_native: False
	- compiled_o3_optim: False

===========  ===========  ===========  ========  ======================
  nb_thread    nb solved    time (ms)    pf / s  speed-up vs 1 thread
===========  ===========  ===========  ========  ======================
          1          703      4213.12       167  1.00x
          2          703      2367.70       297  1.78x
          3          703      1791.98       392  2.35x
          4          703      1451.98       484  2.90x
          5          703      1276.28       551  3.30x
          6          703      1250.20       562  3.37x
          7          703      1230.02       572  3.43x
          8          703      1184.67       593  3.56x
===========  ===========  ===========  ========  ======================

As documented above, speed-up is sub-linear (the per-thread set-up cost, plus the
already-fast single-thread baseline, both eat into the theoretical ``nb_thread``-x
speed-up): most of the gain is captured by 4-5 threads, with diminishing returns beyond
that on this grid / contingency count.

Reporting limit violations
---------------------------

Starting from lightsim2grid 1.0.0, if you have set some operating limits on the grid model
(per-bus voltage bounds with :func:`lightsim2grid.network.LSGrid.set_bus_voltage_limits`,
per-side current limits on lines / trafos with
:func:`lightsim2grid.network.LSGrid.set_line_current_limit_side1` /
:func:`lightsim2grid.network.LSGrid.set_line_current_limit_side2` /
:func:`lightsim2grid.network.LSGrid.set_trafo_current_limit_side1` /
:func:`lightsim2grid.network.LSGrid.set_trafo_current_limit_side2`), ``ContingencyAnalysis`` can
report, for the pre-contingency ("n") case and for each simulated contingency, which of these
limits are violated. The API is modeled after pypowsybl's security analysis
(``res.pre_contingency_result.limit_violations`` / ``res.post_contingency_results``), with the
notable difference that ``post_contingency_results`` is a **list** (ordered like the
contingencies were added), not a dictionary:

.. code-block:: python

    import numpy as np
    import grid2op
    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid import LightSimBackend
    env = grid2op.make(..., backend=LightSimBackend())

    # set some limits on the grid model (kV for buses, kA for lines / trafos)
    grid = env.backend._grid
    nb_bus = grid.get_bus_vn_kv().shape[0]
    grid.set_bus_voltage_limits(0.95 * grid.get_bus_vn_kv(), 1.05 * grid.get_bus_vn_kv())
    grid.set_line_current_limit_side1(line_limit_a1_ka)
    grid.set_line_current_limit_side2(line_limit_a2_ka)
    grid.set_trafo_current_limit_side1(trafo_limit_a1_ka)
    grid.set_trafo_current_limit_side2(trafo_limit_a2_ka)

    # the feature must be explicitly requested at construction time
    security_analysis = ContingencyAnalysis(env, compute_limit_violations=True)
    security_analysis.add_single_contingency(0, name="line_0")  # `name` is optional
    security_analysis.add_all_n1_contingencies()

    res = security_analysis.run()  # or run_ac() / run_dc() to also pick the algorithm family

    for violation in res.pre_contingency_result.limit_violations:
        print(violation.element_type, violation.element_id, violation.side,
              violation.violation_type, violation.value, violation.limit)

    for cont in res.post_contingency_results:  # a list, in the order contingencies were added
        print(cont.element_ids, cont.contingency_name, cont.converged, cont.limit_violations)

Each ``LimitViolation`` reports:

- ``element_type``: :class:`lightsim2grid.contingencyAnalysis.ViolationElementType`
  (``BUS``, ``LINE``, ``TRAFO``, or ``GRID`` for a non-converged contingency, see below);
- ``element_id``: the grid-model bus id (for ``BUS``) or the local line / trafo id (for
  ``LINE`` / ``TRAFO``); unused (``-1``) for ``GRID``;
- ``side``: ``1`` or ``2`` for ``LINE`` / ``TRAFO`` (unused, ``0``, for ``BUS`` / ``GRID``);
- ``violation_type``: :class:`lightsim2grid.contingencyAnalysis.LimitViolationType`
  (``LOW_VOLTAGE``, ``HIGH_VOLTAGE``, ``CURRENT``, ``NOT_SIMULATED`` or ``DIVERGENCE``, see below);
- ``value`` / ``limit``: the value reached and the limit that was violated; unused (``NaN``) for
  ``NOT_SIMULATED`` / ``DIVERGENCE``;
- ``name``: for ``LINE`` / ``TRAFO``, the element's own name (see ``LSGrid.set_line_names`` /
  ``set_trafo_names``); for ``BUS``, the name of the *substation* the violating bus belongs to
  (see ``LSGrid.set_substation_names`` -- there is no per-bus name, only per-substation ones);
  empty string if names were never set on the grid for the relevant kind, or for ``GRID``.

Each ``ContingencyResult`` reports ``element_ids`` (the branch ids disconnected by this
contingency -- always present, even without a ``name``) and the corresponding ``element_names``,
the optional user-supplied ``contingency_name`` (set via ``add_single_contingency(..., name=...)``),
whether the contingency ``converged``, and its ``limit_violations``.

.. warning::

    This feature is **opt-in** and must be requested at construction time
    (``ContingencyAnalysis(env, compute_limit_violations=True)``, or by setting
    ``security_analysis.compute_limit_violations = True`` before running). Leaving it to its
    default (``False``) means ``run`` / ``run_ac`` / ``run_dc`` raise a ``RuntimeError``, and
    the usual ``get_flows()`` is completely unaffected -- there is no need to pay for the extra
    per-element voltage / current checks if you only want the flows.

    Setting ``compute_limit_violations`` through the property (rather than at construction time)
    **clears any contingency already registered** (``add_single_contingency`` /
    ``add_all_n1_contingencies`` / ``add_multiple_contingencies``) -- it is the exact same
    reset performed by ``clear()`` / ``change_algorithm``. If you follow this page's own example
    order (add contingencies, *then* flip the flag), you will silently end up running ``run()``
    on zero contingencies. Either set ``compute_limit_violations=True`` at construction time (as
    in the example above), or re-add the contingencies after setting it:

    .. code-block:: python

        security_analysis = ContingencyAnalysis(env)
        security_analysis.add_all_n1_contingencies()
        security_analysis.compute_limit_violations = True  # clears the contingencies above!
        security_analysis.add_all_n1_contingencies()  # so they must be re-added

        res = security_analysis.run()

    Also, if a contingency does not converge, its ``limit_violations`` contains exactly one
    ``LimitViolation`` with ``element_type == ViolationElementType.GRID`` and ``violation_type``
    either ``LimitViolationType.NOT_SIMULATED`` (a pre-check -- eg it splits the grid in multiple
    connected components -- skipped this contingency without ever invoking the solver) or
    ``LimitViolationType.DIVERGENCE`` (the solver ran but did not converge, eg reached
    ``max_iter``). This is unlike a converged case with an empty ``limit_violations`` list, which
    genuinely means no violation was found.

    The pre-contingency ("n") case is different: if it does not converge, `compute` (and thus
    `run` / `run_ac` / `run_dc`) raises a ``RuntimeError`` instead, since every contingency is
    solved relative to that base case -- a diverging base case makes the whole analysis
    meaningless.

Violation checking is fully compatible with the ``handle_disconnected_grid`` and ``nb_thread``
options described above: it is performed **inline, per contingency, inside the thread that
solves it** (not as a separate post-processing pass), so it does not affect the multi-threading
performance characteristics; and masked (disconnected-island) buses are correctly excluded from
the voltage checks.

Adjusting the violation threshold
++++++++++++++++++++++++++++++++++

By default a violation is reported exactly at the configured limit. The
``violation_threshold`` attribute (available both on ``ContingencyAnalysis`` and on the
lower-level ``ContingencyAnalysisCPP``) lets you tighten that margin, so that situations
approaching a limit are reported before they actually breach it. It is a ``float`` in
``]0., 1.]`` and defaults to ``1.0``, which reproduces the behaviour described above:

.. code-block:: python

    security_analysis = ContingencyAnalysis(env, compute_limit_violations=True)
    security_analysis.add_all_n1_contingencies()

    # report anything that reaches 95% of its limit, instead of only actual breaches
    security_analysis.violation_threshold = 0.95

    res = security_analysis.run()

Setting it below ``1.0`` makes **every** check stricter -- more violations are reported,
never fewer. Each of the three checks owns one interval, running from a "healthy" anchor to
the limit that can be violated, and the threshold moves that limit towards its anchor -- a
linear interpolation, identical for all three:

.. code-block:: text

    effective_limit = threshold * limit + (1 - threshold) * anchor

=================  =========  ===========  ================================================
check              anchor     limit        violates when
=================  =========  ===========  ================================================
``CURRENT``        ``0``      ``limit_a``  ``value >= threshold * limit_a``
``LOW_VOLTAGE``    ``vn_kv``  ``vmin_kv``  ``v <= threshold * vmin + (1 - threshold) * vn``
``HIGH_VOLTAGE``   ``vn_kv``  ``vmax_kv``  ``v >= threshold * vmax + (1 - threshold) * vn``
=================  =========  ===========  ================================================

A line's usable range really *is* ``[0, limit_a]``, so its anchor is ``0`` and the rule
reduces to plainly scaling the limit. A voltage bound has no such natural "zero" end -- so
the bus **nominal** voltage is used instead, since that is what operating limits are
conventionally expressed around (:math:`\pm x\%` of ``vn``). Either way, each acceptable
interval keeps a width of exactly ``threshold`` times its original one. For a
:math:`\pm 5\%` band around nominal:

===============  ==========================  ==============
``threshold``    effective band (pu)         band width
===============  ==========================  ==============
1.0              ``[0.9500, 1.0500]``        ``0.100``
0.95             ``[0.9525, 1.0475]``        ``0.095``
0.9              ``[0.9550, 1.0450]``        ``0.090``
0.5              ``[0.9750, 1.0250]``        ``0.050``
===============  ==========================  ==============

.. note::

    The voltage anchor is really ``vn_kv`` *clamped into* ``[vmin_kv, vmax_kv]``. An operating
    band is **not** required to bracket the nominal voltage, and real data does violate it: a
    380 kV level declared with an operating range of ``[390, 450]`` kV is ordinary practice on
    the European 400 kV network. Such a grid loads and analyses normally; the clamp simply pins
    the anchor to the nearer bound, so that bound has no margin left to give while the opposite
    one still tightens as usual.
    Where the band *does* bracket the nominal voltage -- the overwhelmingly common case -- the
    clamp does nothing at all.

.. note::

    Anchoring both voltage checks on the nominal voltage -- rather than on each other -- is
    what keeps them **independent**. The ``LOW_VOLTAGE`` verdict is a function of ``vmin_kv``,
    the anchor and the threshold alone, so setting ``vmax_kv``, changing it, or leaving it at
    ``NaN`` can never alter it; and symmetrically for ``HIGH_VOLTAGE``. It also means the two
    effective bounds each converge towards the anchor from their own side and can therefore
    never cross: no bus is ever reported as both too low and too high, whatever the
    threshold. A useful corollary is that a bus sitting *above* its nominal voltage can never
    be reported as ``LOW_VOLTAGE``, however small the threshold gets.

The reported ``value`` and ``limit`` of each ``LimitViolation`` are unaffected by all this:
they remain the value actually reached and the limit exactly as you configured it. Only the
test deciding whether a violation is reported at all is shifted.

.. note::

    A bus whose limits are inconsistent -- ``vmin_kv > vmax_kv``, the only way the two
    effective bounds can still end up crossed -- raises a ``RuntimeError`` rather than being
    reported as an arbitrary one of the two violation types, which would hide a genuine
    misconfiguration.

.. note::

    Like ``handle_disconnected_grid`` / ``nb_thread``, ``violation_threshold`` is a plain
    runtime knob and only affects the *next* run. There is one asymmetry worth knowing:
    **lowering** it discards any already-computed results (they were computed with a looser
    threshold and would under-report), while **raising** it keeps them (a stricter result
    still contains everything a looser one would find). Unlike ``compute_limit_violations``,
    neither case clears the registered contingencies, so it is enough to run again:

    .. code-block:: python

        res = security_analysis.run()
        security_analysis.violation_threshold = 0.9  # drops the results above ...
        res = security_analysis.run()                # ... no need to re-add anything

.. _sa_benchmarks:

Benchmarks (Contingency Analysis)
----------------------------------

Here are some benchmarks made with:

- date: 2026-04-21 09:05  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.13.5.final.0 (64 bit)
- numpy version: 2.3.5
- pandas version: 2.3.3
- pandapower version: 3.4.0
- pypowsybl version: 1.15.0
- grid2op version: 1.12.4.dev0
- lightsim2grid version: 0.13.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: False 
	- compiled_o3_optim: True 

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 contingency_analysis.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2_small (177 n-1 simulated)
    Total time spent in "computer" to solve everything: 11.1ms (15913 pf / s), 0.06 ms / pf)
        - time to compute the coefficients to simulate line disconnection: 0.28ms
        - time to pre process Ybus: 0.30ms
        - time to perform powerflows: 10.25ms (17276 pf / s, 0.06 ms / pf)
    In addition, it took 0.50 ms to retrieve the current from the complex voltages (in total 15229.8 pf /s, 0.07 ms / pf)

    Comparison with raw grid2op timings
    It took grid2op (with lightsim2grid, using obs.simulate): 0.28s to perform the same computation
        This is a 24.2 speed up from ContingencyAnalysis over raw grid2op (using obs.simulate and lightsim2grid)
    It took grid2op (with pandapower, using obs.simulate): 9.94s to perform the same computation
        This is a 855.2 speed up from ContingencyAnalysis over raw grid2op (using obs.simulate and pandapower)

    In this case then, the `ContingencyAnalysis` module is 24 times faster than raw grid2op (with obs.simulate
    and lightsim2grid) and 855 times faster than raw grid2op (with obs.simulate and pandapower)
    All results match !

.. note::
  The last "In this case then, ..." line above is now printed directly by ``contingency_analysis.py`` (computed
  from ``total_time_glop_ls`` / ``total_time_glop_pp`` and ``full_time_sa``), instead of being a separate,
  hand-rounded restatement kept below the pasted output -- which is how it used to read "more than **22**
  times faster", inconsistent with the 24.2 / 855.2 figures printed just above it.


Detailed usage
--------------------------

.. automodule:: lightsim2grid.contingencyAnalysis
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
