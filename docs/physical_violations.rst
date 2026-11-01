.. currentmodule:: lightsim2grid

.. _physical_violations:

Physical limits and controls a solution does not apply
=======================================================

Goal
----

A converged powerflow is not necessarily a state the grid can reach. The default algorithms
solve the problem they are given: a generator holds its voltage set-point whatever the reactive power
it takes, a distributed slack asks a unit for the power it has, a standby SVC stays idle, and
so on. Nothing checks whether the equipment could do it, or whether a control would have acted.

``get_physical_violations`` (on a grid) and ``compute_physical_violations`` (on the batch classes)
**report**, after a solve, what the solution assumes that the equipment or its controls would
not allow. Nothing is enforced and nothing is re-solved: to *apply* the controls, run the
``NROuter_*`` algorithms instead (see :ref:`outer_loops`). Each check is the trigger of
an OpenLoadFlow outer loop, and it is the same code the loops use: the report of a plain solve is what the first round of
the loops would act on.

.. note::

    These violations are different from the *operational* ones of ``compute_limit_violations``
    (a bus outside its voltage range, a branch above its thermal rating, see
    :doc:`security_analysis`): those are states the grid reaches and should
    not sit in; these are states it cannot reach, or the control a solution did not apply.

Categories
----------

Every violation has a ``violation_type`` (:class:`~lightsim2grid.contingencyAnalysis.LimitViolationType`) and belongs to exactly
one :class:`~lightsim2grid.contingencyAnalysis.ViolationCategory`, which is the first thing to read: the categories must not be acted on the same way.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Category
     - Meaning
   * - ``OPERATIONAL``
     - A limit chosen by an operator, which the grid can leave: ``LOW_VOLTAGE``, ``HIGH_VOLTAGE``, ``CURRENT``.
       Reported by ``compute_limit_violations``, not here.
   * - ``PHYSICAL``
     - A limit of the equipment itself, which nothing can leave: the solution is **not** realizable
       (``LOW_Q``, ``HIGH_Q``, ``LOW_P``, ``HIGH_P``, ``LOW_VOLTAGE_REMOTE_CONTROL``, ``HIGH_VOLTAGE_REMOTE_CONTROL``).
   * - ``CONTROL``
     - Not a limit: a control or an automaton whose action the solution did not apply. The equipment could do it
       and an outer loop would, so the state is reachable, but it is not the one the controls would settle in.
   * - ``SOLVER``
     - What the solver did: ``NOT_SIMULATED`` and ``DIVERGENCE`` (reported by the batch classes).

``ViolationCategory.CONTROL`` is new with the ``NROuter_*`` algorithms.

What is checked
---------------

The AC checks need an algorithm that publishes its per-bus mismatch (every built-in AC algorithm does).
In DC, only the checks that need no reactive power nor voltage magnitude are made.

.. list-table::
   :header-rows: 1
   :widths: 27 12 13 48

   * - Type
     - Category
     - On
     - Meaning (loop whose trigger it is)
   * - ``LOW_Q`` / ``HIGH_Q``
     - PHYSICAL
     - BUS
     - The reactive power the machines holding a bus' voltage had to produce is beyond what they can (the **sum** of the
       capabilities of the voltage-regulating generators, storage units, hvdc stations and voltage-mode SVCs: how a bus' Q is
       split between its units is a convention, the total is not). ``ReactiveLimits``, PV -> PQ.
   * - ``LOW_VOLTAGE_AT_MIN_Q`` / ``HIGH_VOLTAGE_AT_MAX_Q``
     - CONTROL
     - GENERATOR, SVC, HVDC
     - A unit flagged "pinned at a reactive limit" (``set_gen_can_be_pv``, ``set_svc_can_be_pv``,
       ``set_hvdc_can_be_pv``) whose regulated bus is below its target while it absorbs all it can
       (or above while it produces all it can): it would regulate again. ``ReactiveLimits``, PQ -> PV.
       ``value`` and ``limit`` in kV.
   * - ``LOW_VOLTAGE_REMOTE_CONTROL`` / ``HIGH_VOLTAGE_REMOTE_CONTROL``
     - PHYSICAL
     - GENERATOR
     - A generator regulating a *remote* bus whose own bus is outside the realistic voltage band. ``ReactiveLimits``,
       robust mode. kV.
   * - ``LOW_VOLTAGE_SVC_STANDBY`` / ``HIGH_VOLTAGE_SVC_STANDBY``
     - CONTROL
     - SVC
     - An idle SVC flagged by ``set_svc_standby`` whose monitored bus is outside its automaton's thresholds.
       ``VoltageMonitoring``. kV.
   * - ``HIGH_P`` (hvdc)
     - PHYSICAL
     - HVDC
     - An angle-droop hvdc line whose flow ``p0 + k.(theta1 - theta2)`` is beyond the maximum power of its
       direction (``side``: 1 is 1 -> 2). ``HvdcAcEmulationLimits``. AC and DC.
   * - ``HVDC_AC_EMULATION_RELEASE``
     - CONTROL
     - HVDC
     - An angle-droop line flagged "frozen at its limit" (``set_hvdc_ac_emulation_frozen``) whose droop asks for less
       than that limit. ``HvdcAcEmulationLimits``. MW. AC and DC.
   * - ``LOW_P`` / ``HIGH_P``
     - PHYSICAL
     - GENERATOR, STORAGE
     - A unit of the distributed slack whose share of the imbalance lands it past its ``[min_p, max_p]``
       (``set_gen_p_limits`` / ``set_storage_p_limits``; a grid without them reports nothing). ``DistributedSlack``.
       Generator convention (positive: injected). AC and DC.
   * - ``SLACK_MISMATCH``
     - CONTROL
     - GRID
     - The active power the single slack bus absorbed is above ``slackBusPMaxMismatch``: ``DistributedSlack`` would share it.
       ``value`` is the mismatch in MW (positive: the units must inject more).
   * - ``PHASE_CONTROL_P``
     - CONTROL
     - TRAFO
     - A phase shifter regulating an active power is off its target by more than its deadband. ``PhaseControl``. MW.
   * - ``PHASE_LIMITER_CURRENT``
     - CONTROL
     - TRAFO
     - A phase shifter limiting a current is above the limit. ``PhaseControl``. A.
   * - ``TRANSFORMER_VOLTAGE_DEADBAND``
     - CONTROL
     - BUS
     - A bus whose voltage transformers regulate (and no generator does) is off its target by more than half the
       deadband. ``TransformerVoltageControl``. kV.
   * - ``SHUNT_VOLTAGE_CONTROL``
     - CONTROL
     - BUS
     - A bus whose shunts regulate (and no generator nor transformer does) is off its target by more than half the
       deadband. ``ShuntVoltageControl``. kV.

Some checks need data the caller must provide, because lightsim2grid cannot know which unit an earlier loop froze:
the flags ``set_gen_can_be_pv``, ``set_svc_can_be_pv``, ``set_hvdc_can_be_pv``, ``set_svc_standby`` and
``set_hvdc_ac_emulation_frozen`` (``init_from_pypowsybl`` fills them, from what ``bake_outer_loops`` froze when
``can_be_pv=`` / ``hvdc_ac_emulation_frozen=`` are given). Without the flag, there is nothing to report.

A single powerflow
------------------

.. code-block:: python

    grid.ac_pf(V0, 30, 1e-8)                                  # any AC algorithm
    for v in grid.get_physical_violations(ac=True, tol_mva=1e-4, tol_vm_pu=1e-4):
        print(v.violation_type, v.element_type, v.element_id, v.value, v.limit)

``tol_mva`` (default ``1e-4``) is the absolute slack on every comparison in MW / MVAr, so that an element resting
on its limit is not reported for solver noise; ``tol_vm_pu`` is the same on a voltage comparison.
``ac=False`` checks the last ``dc_pf``. It raises if no powerflow ran.

``grid.get_violations(threshold=1., ac=True, rel_tol=...)`` is the operational counterpart (voltage bands and
thermal limits).

In the batch classes
--------------------

Every batch class (``TimeSeries``, ``InjectionSweep``, ``ContingencyAnalysis``, ``ScenarioSweep``) has the same
three properties, ``False`` / ``1e-4`` / ``1e-4`` by default (``compute_physical_violations``,
``physical_violation_tol_mva``, ``physical_violation_tol_vm_pu``). With
:class:`~lightsim2grid.contingencyAnalysis.ContingencyAnalysis`:

.. code-block:: python

    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid.contingencyAnalysis import ViolationCategory

    ca = ContingencyAnalysis(env)
    ca.add_all_n1_contingencies()
    ca.compute_physical_violations = True       # drops the base case and results, keeps the contingencies
    res = ca.run()

    res.pre_contingency_result.physical_violations      # the "n" case
    for cont in res.post_contingency_results:           # in the order contingencies were added
        for v in cont.physical_violations:
            if v.category == ViolationCategory.CONTROL:
                print(v.violation_type, v.element_type, v.element_id, v.value, v.limit)

``get_physical_violations()`` and ``get_physical_violations_n()`` give the same lists in the C++-side order. The
property can be set at any point before the computation; unlike ``compute_limit_violations``, setting it does
not clear the registered contingencies. A row that did not converge is not checked.

Each ``LimitViolation`` has ``element_type`` (:class:`~lightsim2grid.contingencyAnalysis.ViolationElementType`:
``BUS``, ``LINE``, ``TRAFO``, ``GRID``, ``HVDC``, ``GENERATOR``, ``STORAGE``, ``SVC``), ``element_id``,
``violation_type``, ``value``, ``limit``, ``side`` (where it means something), ``name`` and ``category``.

.. warning::

    ``PHYSICAL`` and ``CONTROL`` violations are *reported on top of* the converged state. A ``DIVERGENCE``
    says nothing about the grid: the solver may have failed on a feasible state.

Changes of behaviour with the outer loops
-----------------------------------------

The reactive range, hvdc droop and standby SVC checks are the outer loops' **own** tests, no longer a separate
copy of them: the tolerance is OpenLoadFlow's, the limits are taken at the current active power target, the load's Q is part of
the bus balance, and the robust mode and the realistic band are those of the grid's ``ReactiveLimits`` loop (or of its
driver parameters, see :ref:`outer-loops-driver`). A bus that was reported before can therefore not be reported
now, or the other way round, close to a limit.
