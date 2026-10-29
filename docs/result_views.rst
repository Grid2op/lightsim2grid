.. _result_views:

Result arrays: views, copies and the lifetime of the grid
==========================================================

Most of the getters of results return **numpy arrays that are views on the memory of the object they
come from**, without copying it. This is what makes reading the results of a step cheap (the grid2op backend reads
all of them after every powerflow), but it has consequences you must know about.

What this means
---------------

- **A view follows the next solve.** ``p, q, v = grid.get_gen_res()`` returns arrays that show the *current*
  results of the grid: after the next ``ac_pf`` they hold the new values. If you need the values of a given solve, copy them
  (``p.copy()``, ``np.array(p)``).
- **The result getters of** :class:`~lightsim2grid.network.LSGrid` (``get_gen_res``, ``get_line_res1``, ``get_load_res``, ...) **do not keep the grid
  alive.** If the grid is garbage-collected (its last reference is dropped), the array points to freed
  memory. Reading it is undefined behaviour: it often still shows the right numbers, and sometimes garbage, and it can
  work on one machine and fail on another.
- **The batch classes** (``TimeSeries``, ``InjectionSweep``, ``ContingencyAnalysis``, ``ScenarioSweep``) return views
  that keep the object they come from alive.
- A few getters always return a **copy**: the reactive ranges (``get_gen_reactive_range_max``,
  ``get_station_reactive_range_max``) and the raw target Q (``get_gen_raw_target_q``). The admittance matrices
  ``get_Ybus_solver()`` / ``get_dcYbus_solver()`` return a read-only view of the grid's matrix, valid until a powerflow rebuilds it.

The pattern to avoid
--------------------

.. code-block:: python

    def reference_production(net):
        grid = init_from_pandapower(net)      # local variable
        grid.ac_pf(V0, 30, 1e-8)
        p, _, _ = grid.get_gen_res()          # a view on the memory of `grid`
        return p                              # `grid` is freed here: `p` dangles

and the fix:

.. code-block:: python

        return p.copy()                       # or keep the grid alive as long as p is used

Rule of thumb: **copy what must outlive the grid, or outlive the next powerflow.**

Status
------

Making every result getter keep its grid alive (and what that costs in memory and time on the hot path of the
grid2op ``step``) is an open item of the ``[TODO]`` list of the changelog. Until it is done, the rule above applies
to the :class:`~lightsim2grid.network.LSGrid` getters.
