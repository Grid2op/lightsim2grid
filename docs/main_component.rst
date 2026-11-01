.. currentmodule:: lightsim2grid.network

.. _main_component:

Islands, the main component and the distributed slack
======================================================

Goal
----

A grid that is not fully connected cannot be solved as it is: a powerflow needs one synchronous
component, with a slack bus. :meth:`LSGrid.consider_only_main_component` keeps the main
synchronous component, disconnects everything else, and says what that cost in active power.
It follows OpenLoadFlow, which also keeps one component only.

.. code-block:: python

    report = grid.consider_only_main_component()          # redistribute_slack=True
    print(report.mismatch_mw, report.nb_participants, report.nb_saturated, report.all_saturated)

What it does
------------

1. **The main component** is the connected component of the branch graph (lines, transformers, and any other
   connecting element) that has the **most buses**, whether a slack bus is in it or not. This is OpenLoadFlow's
   rule, and *not* "the component of the slack": a grid whose slack sits on a small island keeps its large
   component.
2. Every element with no bus in that component is disconnected.
3. **A hvdc line with a single converter station in the main component** is not fully disconnected: the converter
   in the component stays active (still injecting / regulating its scheduled power), and only the one outside
   is opened.
4. **The generators and storage units of the distributed slack that sat on a disconnected bus** stay in the slack
   but take no share of it while disconnected: the slack is only distributed on the main component. They take their
   share again once reconnected (``reactivate_gen``, ``update_topo``, ...).
5. **The reference slack bus** forced with :meth:`LSGrid.set_reference_slack_bus` that is outside the component
   is cleared.
6. **If no slack unit is left in the component**, it is given one as :meth:`LSGrid.assign_slack_to_most_connected`
   does, with the buses compared on their nominal voltage first (OpenLoadFlow's most meshed bus is among the buses of
   the highest nominal voltage): the highest-voltage, most connected bus that has a producing unit absorbs what the islanding
   took out. This **replaces the whole slack**, the stranded units included: work on a copy (``grid.copy()``) if the original slack is
   needed afterwards. It raises, leaving the grid unchanged, when no connected generator nor storage unit is in the
   component to take it, and also when no slack was defined at all before the call.

Redistributing the lost power
-----------------------------

Disconnecting an island removes its generation and its consumption from the balance. With ``redistribute_slack=True``
(the default), that **active power imbalance** is shared on the remaining units of the distributed slack, as
OpenLoadFlow's ``DistributedSlack`` outer loop would (see :ref:`outer_loops`):

- the imbalance is the set-points of the stranded generators and static generators, minus the stranded loads, storage
  units and shunts, plus the stranded hvdc converter stations (a rectifier leaving with its island takes its
  consumption out of the balance);
- it is shared proportionally to the slack weights;
- each unit is kept within its ``[min_p, max_p]`` (``set_gen_p_limits`` / ``set_storage_p_limits``; unbounded when
  not set), a unit that reaches a bound leaves the sharing and what it could not take is shared again on the others;
- a unit never changes the sign of its injection: an injecting unit stops at 0 MW on a negative imbalance, a
  drawing unit stops at 0 MW on a positive one;
- the new set-points are written in the grid, and the units that reached a bound are taken out of the distributed slack,
  so the next powerflow only shares what is left (the change in the losses) on the units that can still move.
  If every unit reaches a bound, all stay in the slack, since a powerflow needs one, and the report says how much could not
  be placed.

The same sharing is available alone, for any imbalance, by :meth:`LSGrid.redistribute_active_power`, and the batch
classes use it for their ``redistribute_slack`` option (also with the ``NROuter_*`` algorithms).

Pass ``redistribute_slack=False`` to leave the set-points untouched: the grid2op backend and ``init_from_pypowsybl``
do. The imbalance is then absorbed by the slack of the powerflow.

Without any active power limit set, nothing saturates: the converged state is the same as without the redistribution,
only the set-points move.

The report
----------

Both functions return a :class:`SlackRedistributionReport`:

.. list-table::
   :widths: 25 75

   * - ``mismatch_mw``
     - the imbalance shared, in MW (positive: the units must inject more). It is filled even when
       ``redistribute_slack`` is ``False``
   * - ``nb_participants``
     - how many units took part
   * - ``nb_saturated``
     - how many reached a bound (and left the distributed slack)
   * - ``nb_rounds``
     - how many rounds the sharing took
   * - ``not_distributed_mw``
     - what could not be placed
   * - ``all_saturated``
     - whether every unit of the distributed slack saturated (in which case none left it)

Units that sat at a limit in a reference solve
----------------------------------------------

OpenLoadFlow leaves a unit out of the distributed slack when it sat at an active limit in the reference
solution. These units can be flagged (``set_gen_can_participate_slack``, ``set_storage_can_participate_slack``,
``..._overshoot``; ``init_from_pypowsybl(can_participate_slack=...)``): the Newton solve does not read them, but this
bounded redistribution does, within their ``[min_p, max_p]``, so that they only move away from the limit they sit at. A unit that
leaves the slack this way keeps its flag, and goes back into the slack when its set-point is moved off the limit or when it
is reconnected off it. When the first slack generator leaves the slack, the angle reference moves to the next one (a constant
angle shift).

Related functions
-----------------

- :meth:`LSGrid.assign_slack_to_most_connected`: picks the bus where units able to hold the slack produce power and with the
  most line / transformer ends, makes every producing unit on it a participant (weighted by its share of the bus' production)
  and clears the rest of the slack. Returns ``(bus_id, gen_id)``.
- :meth:`LSGrid.set_reference_slack_bus`: forces which slack bus is the angle reference, without changing the slack
  set or weights; it is kept by a copy, a pickle and ``save_binary``, and used by the batch classes built from the grid
  as the reference of the whole batch (a contingency that strands it is skipped with ``handle_disconnected_grid``).
- :meth:`LSGrid.cap_slack_at_active_limits`: after a solve, takes the units pushed past an active limit out of the slack at
  that limit and solves again (see :doc:`physical_violations` for the ``LOW_P`` / ``HIGH_P`` records it acts on).

.. note::

    The changelog of this version marks a change: ``consider_only_main_component`` keeps the **largest** component
    and not the slack's, and a component that lost every slack unit now gets one. Code that relied on the slack's
    component being the one kept should check which component remains (``get_gen_status`` and friends after the call).
