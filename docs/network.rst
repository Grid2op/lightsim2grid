LSGrid module (doc in progress)
====================================

The main class of the lightsim2grid python package is the `LSGrid` class, that is a python class created
from the the c++ `LSGrid` (thanks fo pybind11).

This class basically represents a powergrid (what elements it is made for, their electro technical properties etc.)

To create such class, for now the only way is to get it from a pandapower grid (and it does not model every elements there !)

For example, you can init it like (NOT RECOMMENDED, though sometimes needed):

.. code-block:: python

    from lightsim2grid.network import init
    pp_net = ...  # any pandapower grid eg. pp_net = pn.case118()

    lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

A better initialization is through the :class:`lightsim2grid.LightSimBackend.LightSimBackend` class:

.. code-block:: python

    from lightsim2grid.network import init
    # create a lightsim2grid "gridmodel"
    env_name = ... # eg. "l2rpn_case14_test"
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
   by :class:`lightsim2grid.LightSimBackend.LightSimBackend`), whose ``new_values``
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

3. **Solver bus id** — a *compact* index (``0`` … ``total_bus() - 1``) that depends
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
``lsgrid.total_bus()``              number of buses currently seen by the solver
``lsgrid.nb_connected_bus()``       number of connected buses
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


Detailed documentation
--------------------------

.. automodule:: lightsim2grid.network
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
