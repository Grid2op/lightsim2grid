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
- **A view keeps the grid alive.** The arrays the result getters of :class:`~lightsim2grid.network.LSGrid` return
  (``get_gen_res``, ``get_line_res1``, ``get_load_theta``, ``get_gen_target_p``, ...), and the ``get_bus_id`` /
  ``get_yac_eff_*`` of its element containers, have the grid (or the container) as their ``.base``: the grid is not
  garbage-collected while one of them is in use.
- **A diverged powerflow resets the results to NaN**, in the same arrays: there is no result, and the next powerflow
  that converges writes its own there.
- **A view does not survive the re-initialisation of the elements it shows** (``init_generators`` and the other
  ``init_*``, ``set_state``): the grid replaces their buffers. Call the getter again after that.
- **The batch classes** (``TimeSeries``, ``InjectionSweep``, ``ContingencyAnalysis``, ``ScenarioSweep``) return views
  that keep the object they come from alive. Call their getters again after each ``compute`` (and after ``clear``):
  the next one can replace the buffers an earlier array shows.
- **What comes from the solver side is a copy**, taken when the getter is called: the ``*_solver`` getters of
  :class:`~lightsim2grid.network.LSGrid` (``get_V_solver``, ``get_Sbus_solver``, ``get_pv_solver``,
  ``get_Ybus_solver``, ...) and the voltages of an algorithm (``get_V``, ``get_Va``, ``get_Vm``). Later powerflows
  rebuild those buffers, which no view could survive. Such a copy does not follow the next solve: call the getter again.
- **Each of these has an** ``_internal`` **twin** (``get_Ybus_solver_internal``, ``get_V_solver_internal``,
  ``get_V_internal``, ...) that returns the view instead, without a copy, for code that cannot afford the copy. The view
  is valid only until a later call rebuilds what it shows, and nothing checks it: using it correctly is your own
  responsibility, and misused it can corrupt memory (:ref:`internal_methods`).
- A few other getters return a **copy** as well: the reactive ranges (``get_gen_reactive_range_max``,
  ``get_station_reactive_range_max``) and the raw target Q (``get_gen_raw_target_q``).

Keeping one array keeps the whole grid
--------------------------------------

.. code-block:: python

    def reference_production(net):
        grid = init_from_pandapower(net)      # local variable
        grid.ac_pf(V0, 30, 1e-8)
        p, _, _ = grid.get_gen_res()          # a view on the memory of `grid`
        return p                              # `p` keeps `grid` alive

This is safe, and it keeps the whole grid in memory (its matrices, its solver...) for as long as ``p`` is in use.
Return ``p.copy()`` when only the values are needed.

Rule of thumb: **copy what must outlive the next powerflow, or must not hold the grid.**

.. versionchanged:: 1.3.0
    The views of :class:`~lightsim2grid.network.LSGrid` keep the grid alive, a diverged powerflow resets them in place,
    and what comes from the solver side is copied. Before, a view read freed memory once its grid was gone or after a
    diverged powerflow, and a view of the solver side once a later powerflow had rebuilt it.
