.. currentmodule:: lightsim2grid.algorithm

.. _outer_loops:

OpenLoadFlow's outer loops: the ``NROuter_*`` algorithms
=========================================================

Goal
----

The default algorithms (``NR_*``, ``NRSing_*``, ``FDPF_*``, ...) solve a *single, fixed* problem: the
topology, the injections, the tap positions and the voltage-control set-points you give are exactly
what is solved, once, and nothing is adjusted in response to the result (see
:ref:`available-powerflow-solvers`).

The ``NROuter_*`` algorithms wrap a single-slack Newton-Raphson in the **outer loops** of
`OpenLoadFlow <https://github.com/powsybl/powsybl-open-loadflow>`__ (the default load flow of
pypowsybl): solve, look at the result, change an input if a criterion is not met (a generator is
past its reactive limit, a tap changer is not where its regulation wants it, the slack bus
absorbed too much active power, ...), then solve again, until every loop is satisfied or a
round limit is reached.

They are what you use to

- enforce the reactive limits of the generators (including their capability curves),
- let transformers and phase shifters move their taps, and shunts switch their sections,
- share the active power imbalance on the units that can take it, within their limits,
- reproduce what pypowsybl computes (see :ref:`outer-loops-pypowsybl`).

.. note::

    The default algorithms are untouched and give the same results as before. The ``NROuter_*``
    algorithms are **opt-in**: nothing changes until you call ``change_algorithm`` with one of them.
    The limitations listed in the :doc:`disclaimer <disclaimer>` still apply to every other algorithm.

The loops work **on a fixed Jacobian pattern**: whatever they switch (a bus from PV to PQ, a
transformer tap, a shunt section) is done by rewriting *values* in a sparsity pattern
reserved before the first solve. That pattern is computed **by the loops themselves, from the data of the
grid**: each loop of the list declares the rows, columns and switchable buses it may need in any state it can
reach, from the units that can switch, participate in the slack or regulate. A solve therefore keeps a single
symbolic analysis of the Jacobian, and so does every row of the batch classes (see :ref:`batch_outer_loops`).
The same text, in the list of algorithms, is in :ref:`available-powerflow-solvers`.

Quick start
-----------

.. code-block:: python

    import numpy as np
    from lightsim2grid.network import init_from_pypowsybl

    grid = init_from_pypowsybl(net, olf_rules=True)   # see "From pypowsybl" below
    grid.change_algorithm("NROuter_KLU")               # or "NROuter_SparseLU", ...

    V0 = np.full(grid.total_bus(), 1.0 + 0j)
    V = grid.ac_pf(V0, 30, 1e-8)                        # the grid's outer loops are run

    print(grid.get_outer_loops())                       # the loops that run, in order
    stats = grid.get_algo().get_outer_loop_stats()
    print(stats.status, stats.nb_outer_iterations, stats.nr_iterations)

The algorithms available are ``NROuter_SparseLU`` (always), ``NROuter_KLU`` and
``NROuter_KLUStaticPivot`` (KLU), and ``NROuter_NICSLU`` / ``NROuter_CKTSO`` when compiled in. They
are registered **by name only**: they are not members of
:class:`lightsim2grid.algorithm.AlgorithmType` (see :ref:`algorithm_names`).
``NROuter_KLUStaticPivot`` and :ref:`constrained-ordering` keep KLU's pivots valid when the loops
switch rows.

What a plain grid does: without any call to the functions below, the loops are
OpenLoadFlow's **default list** (see below). A grid read by ``init_from_pypowsybl`` runs
what the installed OpenLoadFlow would run (see :ref:`outer-loops-pypowsybl`).

The loops
---------

Seven loops are available, all of them in :mod:`lightsim2grid.algorithm`. The names are
OpenLoadFlow's.

.. list-table::
   :header-rows: 1
   :widths: 22 12 36 30

   * - Loop
     - In the default list
     - What it does
     - What it reads on the grid
   * - :class:`DistributedSlack`
     - yes (1st)
     - Shares the active power the slack bus absorbed on the units that take part in the slack,
       proportionally to their weight, each one kept within its ``[min_p, max_p]``.
     - ``set_gen_can_participate_slack``, ``set_storage_can_participate_slack``,
       ``set_gen_p_limits``, ``set_storage_p_limits``
   * - :class:`HvdcAcEmulationLimits`
     - yes (2nd)
     - Saturates an angle-droop ("AC emulation") hvdc line at its maximum power in each direction,
       and releases it when its droop asks for less.
     - the droop parameters and the power limits of the hvdc lines
   * - :class:`VoltageMonitoring`
     - yes (3rd)
     - Switches on a standby SVC when the voltage it monitors leaves its thresholds.
     - ``set_svc_standby``, ``set_svc_b0``
   * - :class:`ReactiveLimits`
     - yes (4th)
     - Switches a bus from PV to PQ when its generators are past their reactive limits, and back
       when the voltage allows it; in *robust mode*, stops a remote controller whose own bus has an
       unrealistic voltage.
     - reactive limits, ``set_gen_capability_curves``, ``set_gen_raw_target_q``, ...
   * - :class:`PhaseControl`
     - no
     - Regulates the active power (or limits the current) through a phase shifter, then rounds to
       the closest tap.
     - ``set_trafo_phase_tap_changer``, ``set_trafo_phase_tap_regulation``
   * - :class:`TransformerVoltageControl`
     - no
     - Moves the taps of transformers regulating a voltage, then rounds to the closest tap.
     - ``set_trafo_ratio_tap_changer``, ``set_trafo_ratio_tap_regulation``
   * - :class:`ShuntVoltageControl`
     - no
     - Switches the sections of shunts regulating a voltage.
     - ``set_shunt_sections``, ``set_shunt_section_regulation``

The default list is ``DistributedSlack``, ``HvdcAcEmulationLimits``, ``VoltageMonitoring``,
``ReactiveLimits``, in that order. A loop that has nothing to do on the grid (no hvdc line in AC
emulation, no standby SVC, no regulating transformer, ...) is dropped at run time, as OpenLoadFlow's
``isNeeded`` does, so listing it is harmless.

Not implemented: OpenLoadFlow's area interchange control, secondary voltage control,
transformer reactive power control, the automation systems and the three-winding transformers,
and loops written in Python.

``DistributedSlack``
~~~~~~~~~~~~~~~~~~~~

After a solve, the active power mismatch absorbed by the (single) slack bus is compared to
``slack_bus_p_max_mismatch_mw``. Above it, the *cumulative* mismatch is shared again on the units
flagged "can participate in the slack", from their initial targets, proportionally to their
participation weight, each unit clamped to its active limits and never across 0 MW. A unit that
reaches a bound leaves the sharing for this solve. The loop reports ``FAILED`` when a residual remains and
``fail_on_residue`` is ``True``.

.. code-block:: python

    DistributedSlack(slack_bus_p_max_mismatch_mw=1., fail_on_residue=True,
                     p_residue_eps_mw=1e-3, moved_fraction=0.9)

- ``slack_bus_p_max_mismatch_mw``: OpenLoadFlow's ``slackBusPMaxMismatch``. Below it the slack bus
  keeps what it absorbed.
- ``fail_on_residue``: OpenLoadFlow's ``slackDistributionFailureBehavior``: ``True`` is ``FAIL``,
  ``False`` is ``LEAVE_ON_SLACK_BUS``, used when every unit reached a bound before the mismatch was shared.
- ``p_residue_eps_mw``: what is left to share below it counts as nothing.
- ``moved_fraction``: a unit counts as moved when it moved by more than this fraction of
  ``p_residue_eps_mw``.

The loop only exists on a grid where at least one unit is flagged "can participate in the slack".
On a grid with none, the single slack generator absorbs everything, as OpenLoadFlow does with its
``distributedSlack`` parameter off. The moved targets are published in the results (generator and
storage ``P``), and in the batch results (``OuterRowResult.gen_target_p``).

.. note::

    OpenLoadFlow's default ``DC_VALUES`` starting point already shares the imbalance on the
    participating units. To start from the same angles, use
    ``grid.set_dc_distribute_slack_on_can_participate(True)`` before ``dc_pf``, and give the
    ``ac_pf`` the angles of that ``dc_pf``.

``HvdcAcEmulationLimits``
~~~~~~~~~~~~~~~~~~~~~~~~~

For every hvdc line in angle-droop mode, the flow ``p0 + k.(theta1 - theta2)`` is compared with the
maximum power of each direction. Beyond it the line is saturated (its regime becomes "saturated
in this direction"); a saturated line whose droop asks for less than its limit, or for the other
direction, is released. No parameters. The regimes are in ``OuterRowResult.hvdc_regime``.

``VoltageMonitoring``
~~~~~~~~~~~~~~~~~~~~~

An SVC whose standby automaton is in standby is a *voltage monitor*: loaded idle (off, ``Q = 0``),
flagged with ``LSGrid.set_svc_standby`` with the low and high thresholds of the automaton (in pu of the
nominal voltage of its own bus) and, optionally, the set-points it switches on at below and above them.
When the voltage of the monitored bus leaves the thresholds, the SVC regulates at the matching
set-point, for good. Every other algorithm sees such an SVC as off. ``set_svc_b0`` gives the fixed
susceptance of the automaton, stamped in the admittance matrix whether in standby or not.
No parameters.

``ReactiveLimits``
~~~~~~~~~~~~~~~~~~

The loop that enforces the reactive limits.

- **PV to PQ**: the reactive power a bus needed to hold its voltage is beyond what its units can
  produce (their summed limits, plus a tolerance). The bus is switched to PQ and its units are frozen at the limit.
- **PQ to PV**: the voltage of a frozen bus is back across its target on the right side: the units
  regulate again.
- A bus may go back to PV at most ``max_pq_pv_switch`` times, which stops a bus from switching back and forth.
- **Robust mode**: a *remote* controller whose own bus voltage is outside the realistic band stops
  regulating (it is switched to PQ at its target reactive power) and its bus restarts at ``robust_restart_vm_pu``.

The limits are taken at the **current** active power target of each unit: when ``DistributedSlack``
moved it, a capability curve is read at the new value (``set_gen_capability_curves``; interpolated,
extrapolated past the ends).

.. code-block:: python

    ReactiveLimits(max_pq_pv_switch=3, max_reactive_power_mismatch=1e-4, robust_mode=None,
                   min_realistic_voltage=nan, max_realistic_voltage=nan,
                   mismatch_base_mva=100., realistic_voltage_margin=1.02, robust_restart_vm_pu=1.)

``robust_mode``, ``min_realistic_voltage`` and ``max_realistic_voltage`` left to ``None`` / NaN are
the driver's (see :ref:`outer-loops-driver`), as one parameter does both in OpenLoadFlow. The
other parameters are OpenLoadFlow's ``reactiveLimitsMaxPqPvSwitch``, ``maxReactivePowerMismatch``
(in pu of ``mismatch_base_mva``) and its own constants for the robust mode. Each is readable as an
attribute of the loop; the loops are immutable once built.

``PhaseControl``
~~~~~~~~~~~~~~~~

Phase shifters whose tap changer regulates (``set_trafo_phase_tap_regulation``, mode
``ACTIVE_POWER`` in MW, ``CURRENT_LIMITER`` in A, or ``FIXED``) are handled as
OpenLoadFlow's ``CONTINUOUS_WITH_DISCRETISATION``:

- an *active power* controller has its shift solved continuously inside the first solve, then rounded to the closest tap;
- a *current limiter* above its limit moves one tap at a time.

The positions it settles on are results (``TrafoInfo.res_phase_tap_position``); the inputs are not modified.

``TransformerVoltageControl``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Transformers whose ratio tap changer regulates a voltage (``set_trafo_ratio_tap_regulation``: the
regulated bus, a target and a deadband), as OpenLoadFlow's ``AFTER_GENERATOR_VOLTAGE_CONTROL``: a group whose bus is
already regulated by a generator's voltage control is hidden by it, and only has its taps rounded. The ratios are solved
continuously, the generators of the low-voltage buses being frozen as PQ meanwhile, then rounded to the closest tap.

.. code-block:: python

    TransformerVoltageControl(use_initial_tap_position=False,
                              max_controlled_nominal_voltage=-1.,
                              min_target_deadband_kv=0.1)

- ``use_initial_tap_position``: OpenLoadFlow's ``transformerVoltageControlUseInitialTapPosition``.
- ``max_controlled_nominal_voltage``: OpenLoadFlow's ``generatorVoltageControlMinNominalVoltage`` (kV); a
  negative value is the automatic one.
- ``min_target_deadband_kv``: the deadband used for a group none of whose transformers has one.

The positions are results (``TrafoInfo.res_ratio_tap_position``).

``ShuntVoltageControl``
~~~~~~~~~~~~~~~~~~~~~~~

The regulating shunts of a bus (``set_shunt_section_regulation``) are one controller: the
susceptance is solved continuously in the first solve, then shared on the shunts (largest first) and
rounded to their closest section, followed by one more solve. A bus already regulated by a generator, an
SVC, a hvdc station or (when it is earlier in the list) a transformer is left alone. No parameters.
The counts are results (``ShuntInfo.res_section_count``).

Choosing the loops, and their order
-----------------------------------

The loops of a grid are a list, kept in order, and **the order matters**: in a pass, the driver
checks each loop in the order of the list, re-solving after every loop that changed something,
and repeats the passes until a full pass changed nothing (see below). Different orders can therefore
settle on different states, and OpenLoadFlow's order is the reference to reproduce.

.. code-block:: python

    from lightsim2grid.algorithm import (DistributedSlack, ReactiveLimits,
                                         TransformerVoltageControl, ShuntVoltageControl)

    grid.change_algorithm("NROuter_KLU")

    grid.get_outer_loops()        # the loops that will run, in order (OpenLoadFlow's default list)

    grid.clear_outer_loops()      # start a list of your own...
    grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=2.))
    grid.add_outer_loop(ReactiveLimits(max_pq_pv_switch=5))
    grid.add_outer_loop(TransformerVoltageControl())
    grid.add_outer_loop(ShuntVoltageControl())   # ... after the transformers, as OpenLoadFlow does

    grid.reset_outer_loops()      # back to OpenLoadFlow's default list

- ``add_outer_loop`` **appends**; there is no insertion. To change the order, ``clear_outer_loops()``
  then add the loops in the order wanted. Calling ``add_outer_loop`` on a grid whose list was never
  edited first loads the default list, then appends to it.
- A list cannot hold the same loop twice (OpenLoadFlow refuses it as well): adding it raises a ``RuntimeError``.
  ``None`` is refused as well.
- ``clear_outer_loops()`` leaves a list with no loop: an ``NROuter_*`` algorithm then is the single-slack
  Newton-Raphson alone.
- The list is copied with the grid, kept by a pickle and by ``save_binary`` (binary format 22), and it is the one
  the batch classes run.

Rules of thumb, all from OpenLoadFlow:

- ``DistributedSlack`` first: it shares the imbalance the others will see.
- ``ReactiveLimits`` before the transformer, phase and shunt loops: a generator holding a bus decides who regulates it.
- ``ShuntVoltageControl`` after ``TransformerVoltageControl``: a transformer regulating the same bus takes precedence.
- ``PhaseControl``, ``TransformerVoltageControl`` and ``ShuntVoltageControl`` are *off* in OpenLoadFlow by default,
  like here: add them to run them.

How the loops are run
~~~~~~~~~~~~~~~~~~~~~

The driver follows OpenLoadFlow's engine:

1. keep the loops that apply to the grid, and initialize them in order;
2. first solve of the Newton-Raphson; the loops only run if it converged;
3. *passes*: each loop, in order, repeats its check, re-solving after each ``UNSTABLE`` result, until it is
   stable, a solve diverges, or the global cap on outer iterations is reached. Passes are repeated while
   a pass re-solved;
4. a loop that returns ``FAILED`` stops everything;
5. the unrealistic-voltage check runs after the last loop able to fix it;
6. the loops are cleaned up in reverse order.

.. _outer-loops-driver:

Driver parameters
-----------------

The parameters of this driver are the grid's, a single object for the whole list:

.. code-block:: python

    from lightsim2grid.algorithm import OuterLoopDriverParams

    grid.set_outer_loop_driver_params(OuterLoopDriverParams(
        max_outer_iterations=20,
        voltage_remote_control_robust_mode=True,
        min_realistic_voltage=0.5, max_realistic_voltage=2.,
        min_nominal_voltage_realistic_check=0.))
    grid.get_outer_loop_driver_params()

- ``max_outer_iterations``: OpenLoadFlow's ``maxOuterLoopIterations``; the cap on the outer iterations of one solve, over every loop.
- ``voltage_remote_control_robust_mode``: OpenLoadFlow's ``voltageRemoteControlRobustMode``.
- ``min_realistic_voltage`` / ``max_realistic_voltage`` (pu), ``min_nominal_voltage_realistic_check`` (kV):
  the realistic voltage band, checked on the buses of at least that nominal voltage. A solve ending outside is reported as
  ``unrealistic_state``.

The defaults are OpenLoadFlow's own source defaults. ``init_from_pypowsybl`` sets those of the
installed OpenLoadFlow. They are kept by a copy, a pickle and ``save_binary``.

``resume_outer_state`` (property of an ``NROuter_*`` algorithm, off by default) is **not** one of
them: with it, a solve starts its loops from where the last converged solve left them (buses
frozen PQ, controllers held, standby SVCs switched on, hvdc regimes, taps, sections), instead of from
the grid. Use it for a sequence of ``ac_pf`` calls following each other (a time series). It is
not OpenLoadFlow's. The batch classes do the same on their own for a ``TimeSeries``.

What the loops read on the grid
-------------------------------

Everything a loop decides on is data of the :class:`~lightsim2grid.network.LSGrid`. A grid loaded with
``init_from_pypowsybl`` has it all; for a grid built another way, these are the setters (see their
docstrings for the exact conventions):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Setter
     - What it provides
   * - ``set_gen_can_participate_slack``, ``set_storage_can_participate_slack``, ``..._overshoot``
     - the units that share the slack, with their weight (the overshoot: how far beyond its limit a
       unit sat in a reference solve)
   * - ``set_gen_p_limits``, ``set_storage_p_limits``
     - active power limits (optional, never enforced by an algorithm that runs no loop)
   * - ``set_gen_capability_curves``, ``gen_limits_at``
     - the reactive capability curve of a generator, read at the current active power target
   * - ``set_gen_raw_target_q``
     - the target Q before it was clamped into the limits, for a non-regulating unit whose target P was moved
   * - ``set_reactive_dispatch_olf``, ``set_gen_reactive_range_max``, ``set_station_reactive_range_max``
     - OpenLoadFlow's split of a bus' reactive power between its units (results only)
   * - ``set_gen_can_be_pv``, ``set_svc_can_be_pv``, ``set_hvdc_can_be_pv``, ``set_hvdc_ac_emulation_frozen``
     - elements a previous loop froze at a limit, which could regulate again (read by the physical checks, see :doc:`physical_violations`)
   * - ``set_svc_standby``, ``set_svc_b0``
     - standby automata of SVCs and their fixed susceptance
   * - ``set_trafo_ratio_tap_changer``, ``set_trafo_phase_tap_changer``
     - the steps of the tap changers (ratio, shift, corrections of r, x, g, b in %) and the current position
   * - ``set_trafo_ratio_tap_regulation``, ``set_trafo_phase_tap_regulation``
     - what a tap changer regulates (bus, target, deadband / mode, side)
   * - ``set_shunt_sections``, ``set_shunt_section_regulation``, ``change_shunt_section_count``
     - the sections of a shunt (cumulative p and q at 1 pu) and what they regulate
   * - ``set_hold_frozen_regulators``, ``set_keep_vinit_at_group_controlled_buses``
     - options on the voltage-control groups: hold a regulator frozen at a limit in its group (it can then be released by
       rewriting values), and keep the starting magnitude at the buses a group regulates instead of setting it to the set-point
   * - ``set_dc_distribute_slack_on_can_participate``
     - makes ``dc_pf`` share the imbalance like the ``DistributedSlack`` loop (for the starting angles)

The converters of :mod:`lightsim2grid.network` fill what their source provides.

Reading what the loops did
--------------------------

After a solve, ``grid.get_algo().get_outer_loop_stats()`` returns an :class:`OuterLoopStats`:

- ``status`` (:class:`OuterLoopStatus`: ``STABLE``, ``UNSTABLE`` or ``FAILED``) and ``failed_loop``;
- ``unrealistic_state``: whether the realistic-voltage check failed;
- ``nb_outer_iterations``, ``nb_passes``, ``loop_iterations`` (how many times each loop changed something, in list order) and
  ``nr_iterations`` (the Newton iterations of every inner solve, the first one included);
- ``decisions``: a list of :class:`OuterDecision`, **every** decision the loops took *or declined*, with the
  loop (``loop``), the outer iteration, the ``action``, whether it was ``taken``, the element, the ``reason`` and
  the ``value`` and ``limit`` that were compared. It is how to find out why a bus was (not) switched, and how
  close a decision was to its threshold.

The results of the grid are those of the state the loops settled on: generator output at the targets that
``DistributedSlack`` moved them to, reactive power frozen at the limits, flows at the taps and sections
the loops chose. The tap positions and the section counts are *results*
(``TrafoInfo.res_ratio_tap_position`` / ``res_phase_tap_position``, ``ShuntInfo.res_section_count``): the
inputs of the grid are never rewritten, as pypowsybl's ``solved_tap_position`` / ``solved_section_count``.

``get_physical_violations`` reports what the grid's loops *would* do after a solve by an algorithm that did
not run them: see :doc:`physical_violations`.

In the batch classes
--------------------

``TimeSeries``, ``InjectionSweep``, ``ContingencyAnalysis`` and ``ScenarioSweep`` run the grid's loops on
each row when their algorithm is an ``NROuter_*`` one, with one symbolic analysis per thread, and
``get_outer_results()`` returns an ``OuterRowResult`` per row. See :ref:`batch_outer_loops` for the details
(a ``TimeSeries`` row resumes from the previous one, a row seeded with ``init_from_n_powerflow`` from the
"n" solve's loops, ``redistribute_slack`` with the loops, contingencies on shunts, loads and generators).

.. _outer-loops-pypowsybl:

From pypowsybl
--------------

``init_from_pypowsybl`` is the way to get everything the loops read:

.. code-block:: python

    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import init_from_pypowsybl

    net = pp.network.load("grid.xiidm", {"iidm.die.with-extensions": "all"})
    params = lf.Parameters(distributed_slack=True, use_reactive_limits=True,
                           transformer_voltage_control_on=True)
    grid = init_from_pypowsybl(net, olf_rules=True, olf_parameters=params)
    grid.change_algorithm("NROuter_KLU")

- **The tap changers, the shunt sections with their regulation, the capability curves, the raw target Q, the
  standby SVC set-points, ``b0`` and OpenLoadFlow's reactive dispatch** are read from the network.
- ``olf_parameters`` (a ``pypowsybl.loadflow.Parameters``, ``None``: the installed defaults, which differ between builds of
  pypowsybl) gives the grid's **list of loops, their order and their parameters, and the driver's**: which loops exist
  follows OpenLoadFlow's switches (``distributed_slack``, ``hvdc_ac_emulation``, ``use_reactive_limits``,
  ``phase_shifter_regulation_on``, ``transformer_voltage_control_on``, ``shunt_compensator_voltage_control_on``,
  the ``svcVoltageMonitoring`` provider parameter), and ``outerLoopNames`` gives the list and its order when set. A loop
  OpenLoadFlow would run that lightsim2grid does not have is left out **with a warning**.
- ``olf_rules`` (``True``, or an ``OlfLoadingParameters``) applies **OpenLoadFlow's network-loading rules**: a generator
  OpenLoadFlow would not let regulate a voltage is added without voltage control, a non-regulating generator injects its target Q
  clamped into its limits, the units OpenLoadFlow's distributed slack shares on get its participation key and active power
  range, a line with both ends on the same bus is disconnected, etc. ``False`` (default) reads the network as written.

The loops and the driver parameters can also be taken from a ``Parameters`` object for a grid not read from
pypowsybl:

.. code-block:: python

    from lightsim2grid.network import olf_outer_loops, set_olf_outer_loops

    loops, driver = olf_outer_loops(params)    # (list of loops in OpenLoadFlow's order, OuterLoopDriverParams)
    set_olf_outer_loops(grid, params)           # clear_outer_loops + add_outer_loop for each + set_outer_loop_driver_params

See :doc:`comparison_with_pypowsybl` for how the results compare with pypowsybl's.

Performance and the linear solver
---------------------------------

- One symbolic analysis (``analyze``) and one numeric factorization per solve, whatever the loops do; the loops
  only trigger ``refactorize`` calls. A pivot order frozen by KLU can turn bad when a loop switches a row: the
  refactorization fallback is therefore **always on** for the ``NROuter_*`` algorithms, and counted in
  :class:`LinearSolverStats` (``nb_refactorize_failed``, ``nb_fallback_factorize``, ...), never hidden.
- Use ``NROuter_KLUStaticPivot`` and/or ``constrained_ordering`` to make that fallback rare: see
  :ref:`static-pivot-refactorization` and :ref:`constrained-ordering`.
- ``NROuter_*`` is a single-slack Newton-Raphson: the distributed slack is the ``DistributedSlack`` loop around it,
  unlike ``NR_*`` whose distributed slack is part of the Jacobian.

Known differences with OpenLoadFlow
-----------------------------------

- The slack of lightsim2grid is a generator; OpenLoadFlow's is a bus. Compare with the same slack (see
  :doc:`comparison_with_pypowsybl`).
- OpenLoadFlow stops its Newton-Raphson on the RMS of the mismatch, lightsim2grid on its infinity norm. Use tight
  tolerances on both.
- A loop acts as soon as a value crosses a threshold, like OpenLoadFlow's. A decision taken at a margin smaller than the
  noise on the value (the Q tests of ``ReactiveLimits`` and the deadband test of ``TransformerVoltageControl``
  are the ones that come close) can differ from one machine or build to another, as on the OpenLoadFlow side.
  ``OuterLoopStats.decisions`` shows the margin.
- Not implemented: see the list under "The loops".

Detailed API
------------

The classes are in :mod:`lightsim2grid.algorithm`: :class:`BaseOuterLoop`, :class:`DistributedSlack`,
:class:`HvdcAcEmulationLimits`, :class:`VoltageMonitoring`, :class:`ReactiveLimits`,
:class:`PhaseControl`, :class:`TransformerVoltageControl`, :class:`ShuntVoltageControl`,
:class:`OuterLoopDriverParams`, :class:`OuterLoopStatus`, :class:`OuterLoopStats` (the full API is
in :doc:`solvers`).
