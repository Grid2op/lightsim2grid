LSGrid module
====================================

The main class of the lightsim2grid python package is the `LSGrid` class, that is a python class created
from the c++ `LSGrid` (thanks fo pybind11).

This class basically represents a powergrid (what elements it is made for, their electro technical properties etc.)

.. _network-init-formats:

Supported source formats
--------------------------

An `LSGrid` can be built from several source formats, each with a dedicated ``init_from_*`` function in
``lightsim2grid.network`` (none of them model every element the source format itself supports):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Function
     - Source format
   * - :func:`~lightsim2grid.network.init_from_pandapower`
     - a `pandapower <http://www.pandapower.org/>`_ network (``pandapowerNet``)
   * - :func:`~lightsim2grid.network.init_from_pypowsybl`
     - a `pypowsybl <https://pypowsybl.readthedocs.io/>`_ network (iidm format)
   * - :func:`~lightsim2grid.network.init_from_matpower`
     - a `MATPOWER <https://matpower.org/>`_ case (``.m`` or ``.mat`` file, or an
       already-parsed dict)
   * - :func:`~lightsim2grid.network.init_from_powermodels`
     - a `PowerModels.jl <https://lanl-ansi.github.io/PowerModels.jl/stable/>`_ network
       data dictionary
   * - :func:`~lightsim2grid.network.init_from_pf_delta`
     - a row of the PFΔ dataset (either already parsed into a dict, or a path to its
       ``.json`` file) -- wraps a PowerModels network dict under a ``"network"`` key
       and delegates to ``init_from_powermodels``

See the "Detailed documentation" section below for the full signature and caveats of each.

For example, you can init it from a pandapower grid like (NOT RECOMMENDED, though sometimes needed):

.. code-block:: python

    from lightsim2grid.network import init_from_pandapower
    pp_net = ...  # any pandapower grid eg. pp_net = pn.case118()

    lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

A better initialization is through the :class:`lightsim2grid.lightSimBackend.LightSimBackend` class:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend
    # create a lightsim2grid "LSGrid"
    env_name = ... # eg. "l2rpn_case14_sandbox"
    env = grid2op.make(env_name, backend=LightSimBackend())
    grid_model = env.backend._grid

.. warning::
    We do not recommend to manipulate directly the :class:`lightsim2grid.network.LSGrid` directy, but to use
    it via the backend class. This is much more tested this way.

.. _bus-labelling:

Bus labelling conventions
-------------------------

A recurring source of confusion is that lightsim2grid manipulates bus ids in
**three different conventions**. A given integer (say ``2``) does *not* refer to
the same bus in all of them, so it is important to know which convention a given
method expects or returns.

.. note::
    Internally (c++ side) these conventions are even distinct *types*
    (``LocalBusId``, ``GridModelBusId`` / ``GlobalBusId`` and ``SolverBusId``), so
    that an accidental conversion between them is caught at compile time. Python
    only sees plain integers, hence this section.

1. **Local bus id** — the busbar number *inside a substation*. It is ``-1`` for a
   disconnected element, or between ``1`` and ``n_busbar_per_sub``. This is the
   grid2op convention: the value you put in a ``set_bus`` action and what you read
   in grid2op's ``topo_vect``. It is the convention of
   :func:`lightsim2grid.network.LSGrid.update_topo` (the bulk topology update used
   by :class:`lightsim2grid.lightSimBackend.LightSimBackend`), whose ``new_values``
   array is indexed by the position in the topology vector (``pos_topo_vect``) and
   holds *local* busbar ids. An element's substation is given by its ``sub_id`` and
   its slot in ``topo_vect`` by ``pos_topo_vect``.

2. **GridModel bus id** (a.k.a. *global* bus id) — the index of a bus in the whole
   ``LSGrid``, between ``0`` and ``n_sub * n_busbar_per_sub - 1``. **This is the
   convention of essentially every "by id" public ``LSGrid`` method** (``change_bus_*``
   / ``get_bus_*``, ``deactivate_bus`` / ``reactivate_bus``, ``set_gen_regulated_bus``,
   the ``bus_id`` / ``bus1_id`` / ``bus2_id`` fields of the ``*Info`` objects, …) and
   of the "user facing" matrices/vectors (``get_Ybus``, ``get_Sbus``, ``get_V``,
   ``get_pv``, …).

3. **Solver bus id** — a *compact* index (``0`` … ``nb_connected_bus() - 1``) that depends
   on the current topology: only the buses actually in service get a solver id.
   It is what is passed to the linear/powerflow solver, so it is the convention of
   everything with a ``_solver`` suffix (``get_Ybus_solver``, ``get_V_solver``,
   ``get_pv_solver``, ``get_J_solver``, …) and of the Jacobian-column mappings
   returned by the solver itself (``get_theta_to_J_col`` / ``get_vm_to_J_col`` /
   ``get_q_to_J_col``, see :ref:`use-solver`).

The mapping between conventions 2 and 3 is available (as numpy arrays) through:

================================   ===========================================================
Method                              Meaning
================================   ===========================================================
``lsgrid.id_ac_solver_to_me()``     array indexed by *AC solver* bus id -> *GridModel* bus id
``lsgrid.id_me_to_ac_solver()``     array indexed by *GridModel* bus id -> *AC solver* bus id
``lsgrid.id_dc_solver_to_me()``     array indexed by *DC solver* bus id -> *GridModel* bus id
``lsgrid.id_me_to_dc_solver()``     array indexed by *GridModel* bus id -> *DC solver* bus id
``lsgrid.total_bus()``              total number of buses (``n_sub * n_busbar_per_sub``)
``lsgrid.nb_connected_bus()``       number of buses currently seen by the solver
================================   ===========================================================

Which convention each "by id" method uses:

.. list-table::
   :header-rows: 1
   :widths: 50 18 32

   * - Method(s)
     - get / set
     - Bus id convention
   * - ``update_topo(has_changed, new_values)``
     - set (bulk)
     - Local (per-substation)
   * - ``change_bus_load`` / ``change_bus_gen`` / ``change_bus_sgen`` /
       ``change_bus_shunt`` / ``change_bus_storage`` / ``change_bus_svc``
     - set
     - GridModel (global)
   * - ``change_bus1_powerline`` / ``change_bus2_powerline`` /
       ``change_bus1_trafo`` / ``change_bus2_trafo`` /
       ``change_bus1_dcline`` / ``change_bus2_dcline``
     - set
     - GridModel (global)
   * - ``deactivate_bus`` / ``reactivate_bus``
     - set
     - GridModel (global)
   * - ``set_gen_regulated_bus(gen_id, regulated_bus)``
     - set
     - GridModel (global)
   * - ``get_bus_load`` / ``get_bus_gen`` / ``get_bus_sgen`` /
       ``get_bus_shunt`` / ``get_bus_storage`` / ``get_bus_svc`` /
       ``get_bus1_powerline`` / ``get_bus2_powerline`` / ...
     - get
     - GridModel (global)
   * - ``LoadInfo.bus_id`` (and ``bus1_id`` / ``bus2_id`` / ...,
       ``regulated_bus_id`` of the ``*Info`` objects)
     - get (read-only)
     - GridModel (global)
   * - ``get_Ybus`` / ``get_Sbus`` / ``get_V`` / ``get_Va`` / ``get_Vm`` /
       ``get_pv`` / ``get_pq`` / ``get_slack_ids``
     - get
     - GridModel (global)
   * - ``get_Ybus_solver`` / ``get_Sbus_solver`` / ``get_V_solver`` /
       ``get_pv_solver`` / ``get_pq_solver`` / ``get_slack_ids_solver`` /
       ``get_J_solver``
     - get
     - Solver
   * - ``solver.get_J`` / ``solver.get_theta_to_J_col`` /
       ``solver.get_vm_to_J_col`` / ``solver.get_q_to_J_col``
     - get
     - Solver

.. _elements-modeled:

Elements modeled
------------------

Substations
+++++++++++++++++++++

:func:`~lightsim2grid.network.LSGrid.get_substations` (alias ``get_voltage_levels``) returns a
:class:`lightsim2grid.elements.SubstationContainer`: like every other ``*Container`` on this page it
supports ``len(...)``, indexing and iteration over :class:`lightsim2grid.elements.SubstationInfo`
objects.

.. autoclass:: lightsim2grid.elements.SubstationContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.SubstationInfo
    :members:
    :autosummary:

Generators (standard)
+++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.GeneratorContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.GenInfo
    :members:
    :autosummary:

A generator can also perform **remote voltage control**, *ie* regulate the
voltage of a bus different from the one it is connected to. Use
:func:`lightsim2grid.network.LSGrid.set_gen_regulated_bus` to set the regulated
bus (it defaults to the generator's own bus, which corresponds to local control).
This is read automatically when initializing the grid from pypowsybl. The same
mechanism is used by the :ref:`svc-section` below.

.. warning::
    When the grid is read from pypowsybl, the regulated bus is resolved **once**, at
    import time, and stored by its (fixed) lightsim2grid global bus id. If the
    regulated element is later moved to another bus *inside lightsim2grid* (e.g.
    through a ``change_bus_*`` / topology change), the controller keeps regulating the
    bus resolved at import: the lightsim2grid grid and the original pypowsybl grid then
    desynchronise. Re-import the grid (or call ``set_gen_regulated_bus`` again) if you
    need to follow such a topology change.

Static Generators (more exotic)
++++++++++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.SGenContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.SGenInfo
    :members:
    :autosummary:

Loads and Storage Units
++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.LoadContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.LoadInfo
    :members:
    :autosummary:

Storage units (batteries) are modeled as PQ injections too, but exposed through a
dedicated container. They use the **load convention**: a positive ``target_p`` means
the unit is charging (power drawn from the grid), a negative ``target_p`` means it is
discharging (power injected in the grid). Note that this is the opposite of the
PowSyBl / IIDM (generator) convention; :func:`lightsim2grid.network.init_from_pypowsybl`
negates the battery setpoints accordingly.

.. autoclass:: lightsim2grid.elements.StorageContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.StorageInfo
    :members:
    :autosummary:

.. _svc-section:

Static Var Compensators (SVC)
+++++++++++++++++++++++++++++++

Static Var Compensators (SVC) are shunt-connected devices that can regulate
voltage (or reactive power). Each SVC has a ``regulation_mode``:

- ``0`` (OFF): the device does not regulate anything;
- ``1`` (VOLTAGE): it maintains ``target_vm_pu`` at its regulated bus, possibly
  with a non-zero ``slope_pu`` (droop);
- ``2`` (REACTIVE_POWER): it injects ``target_q_mvar``.

Like generators, an SVC can regulate a remote bus (see ``regulated_bus_id``).
The susceptance limits ``b_min`` / ``b_max`` are stored for information but are
**never enforced** by the powerflow.

.. autoclass:: lightsim2grid.elements.SvcContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.SvcInfo
    :members:
    :autosummary:

Shunts
++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.ShuntContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.ShuntInfo
    :members:
    :autosummary:

Lines
++++++

.. autoclass:: lightsim2grid.elements.LineContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.LineInfo
    :members:
    :autosummary:

Transformers
+++++++++++++

.. autoclass:: lightsim2grid.elements.TrafoContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.TrafoInfo
    :members:
    :autosummary:


HVDC Lines (more exotic)
++++++++++++++++++++++++++++

HVDC links are modeled inside the AC (Newton-Raphson) and DC powerflow. Each
link is made of two converter stations (VSC or LCC, see
:class:`lightsim2grid.elements.ConverterStationInfo`) and can operate either at
a fixed active power setpoint or in angle-droop (AC emulation) mode. The droop
regime can be inspected / forced with
:func:`lightsim2grid.network.LSGrid.set_status_droop_hvdc` and
:func:`lightsim2grid.network.LSGrid.get_status_droop_hvdc`.

.. note::
    The container used to be called ``DCLineContainer`` (and the info object
    ``DCLineInfo``). These names are still importable from
    ``lightsim2grid.elements`` as deprecated aliases of ``HvdcLineContainer`` /
    ``HvdcLineInfo``.

.. autoclass:: lightsim2grid.elements.HvdcLineContainer
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.HvdcLineInfo
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.ConverterStationInfo
    :members:
    :autosummary:

.. _ptdf-lodf-section:

PTDF / LODF
------------

As long as the topology of the grid is not modified, a DC powerflow is a *linear* function of the
bus injections, so it can be replaced by a matrix multiplication -- much faster than solving the
linear system again for every new injection or contingency (see :ref:`benchmark-dc-solvers` for
numbers).

- :func:`~lightsim2grid.network.LSGrid.get_ptdf` (or :func:`~lightsim2grid.network.LSGrid.get_ptdf_solver`
  for the *solver* bus labelling) returns the Power Transfer Distribution Factor matrix: how much the
  flow on each powerline / transformer changes for a 1 MW injection change at each bus.
- :func:`~lightsim2grid.network.LSGrid.get_lodf` returns the Line Outage Distribution Factor matrix:
  how much the flow on each powerline / transformer changes when another one is disconnected --
  the tool of choice for an n-1 contingency analysis restricted to DC (see also
  :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis` for the general AC/DC case).
- :func:`~lightsim2grid.network.LSGrid.get_Bf` returns the sparse "bus to branch" susceptance
  matrix these are built from.

Both ``get_ptdf`` and ``get_lodf`` require a DC powerflow (``dc_pf``) to have been run first, and
are only valid for the topology that was in place when that powerflow was solved -- any topology
change invalidates them. See each function's own documentation below for a full worked example.

.. _cache-reuse-section:

Solver cache reuse
-------------------

Solving a powerflow is not only the linear algebra: the grid must first be turned into what the
solver consumes -- a compact bus labelling, the admittance matrix ``Ybus``, the injection vector
``Sbus``, the PV / PQ split, the slack weights. On a small grid that assembly is worth a good
fifth of the total time, and it is almost entirely redundant between two consecutive powerflows:
changing one load's setpoint does not move a single admittance coefficient.

So ``LSGrid`` keeps what it built and re-stamps only the parts that changed. Every method that
modifies the grid (``change_*``, ``deactivate_*``, ``reactivate_*``, …) records what it
invalidated, and each powerflow rebuilds exactly that much.

**This is on by default and needs nothing from you.** Every powerflow marks its own solver
family "in sync" on the way out.

.. versionchanged:: 1.0.0
    Before 1.0.0 you had to call ``lsgrid.unset_changes()`` yourself after each powerflow, or
    silently pay for a full rebuild every time. That call is now unnecessary (it does nothing
    when cache reuse is enabled, which is the default) and it is kept only for backward
    compatibility.

The AC and the DC solver cache **independently**: each has its own bus labelling, its own matrix
(``Ybus`` / ``Bbus``), its own injections, its own PV / PQ split and slack weights. An AC
powerflow never marks, invalidates or overwrites anything belonging to the DC family, and the
reverse. Hence the per-family accessors :func:`~lightsim2grid.network.LSGrid.get_ac_pv_solver`,
:func:`~lightsim2grid.network.LSGrid.get_dc_pv_solver`, and their ``pq`` / ``slack_weights``
counterparts.

Controlling it
+++++++++++++++

======================================================   =============================================================
Method                                                    Meaning
======================================================   =============================================================
``lsgrid.allow_cache_reuse(bool)``                        turn reuse on (default) or off, for both families
``lsgrid.allow_ac_cache_reuse(bool)``                     ... for the AC family only
``lsgrid.allow_dc_cache_reuse(bool)``                     ... for the DC family only
``lsgrid.get_allow_cache_reuse()``                        ``True`` iff **both** families may reuse
``lsgrid.get_allow_ac_cache_reuse()``                     is the AC family allowed to reuse?
``lsgrid.get_allow_dc_cache_reuse()``                     is the DC family allowed to reuse?
``lsgrid.prevent_cache_reuse()``                          drop what both families cached, once
``lsgrid.prevent_ac_cache_reuse()``                       ... for the AC family only
``lsgrid.prevent_dc_cache_reuse()``                       ... for the DC family only
======================================================   =============================================================

Note the difference: ``allow_*`` is a **mode** (it stays until you change it back), while
``prevent_*`` is a **one-shot invalidation** -- the family throws away what it had and caches
again from the next powerflow on. ``prevent_cache_reuse()`` is the function historically called
``tell_solver_need_reset()``, which still works.

You should not normally need either. The two cases that call for them are:

- **You suspect a caching bug.** ``allow_cache_reuse(False)`` makes every powerflow rebuild
  everything from the containers; the two runs must agree to the last bit.

  .. code-block:: python

      v_cached = lsgrid.ac_pf(v_init, 10, 1e-8)
      lsgrid.allow_cache_reuse(False)
      v_rebuilt = lsgrid.ac_pf(v_init, 10, 1e-8)
      assert (abs(v_cached - v_rebuilt) < 1e-12).all()

- **You modified the grid behind ``LSGrid``'s back**, through something other than its own
  ``change_*`` / ``deactivate_*`` / ``reactivate_*`` methods -- then nothing recorded the
  invalidation, and ``prevent_cache_reuse()`` (or the narrower ``tell_recompute_ybus`` /
  ``tell_recompute_sbus``) is how you say so.

.. note::
    A wrong "nothing changed" claim can cost you a rebuild you were trying to avoid; it can never
    make lightsim2grid read memory it does not own. Every powerflow checks that the data the flags
    describe is actually there before reusing it, and rebuilds from scratch otherwise.

What is never cached across
++++++++++++++++++++++++++++

**Serialization.** Nothing the solvers cache is written to a pickle or a binary file, and nothing
is read back: a grid restored through :func:`~lightsim2grid.network.LSGrid.load_binary` or
``pickle.loads`` always starts cold and rebuilds on its first powerflow. This is a security
property rather than a performance one. A cache is a second copy of state the elements already
determine; read back from a file it becomes a copy that cannot be checked against the elements it
claims to describe. ``check_grid()`` can validate that an index is in range -- it cannot validate
that a matrix really is the admittance matrix of the grid stored next to it, and one that merely
looked well-formed would be solved without complaint. Files are not trusted input, so the cache is
rebuilt, once, from data that is.

**Copying.** :func:`~lightsim2grid.network.LSGrid.copy` does not carry the cache either: the copy
starts cold and rebuilds on its first powerflow. Unlike the serialization case this is not a safety
requirement -- a copy is the same grid, in the same process, so its cache would be perfectly valid
-- and it may change in a future version. The ``allow_*_cache_reuse`` settings *are* copied.

Detailed documentation
--------------------------

.. automodule:: lightsim2grid.network
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
