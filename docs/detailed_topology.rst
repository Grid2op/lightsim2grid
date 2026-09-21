.. _detailed-topology:

Detailed topology: the switches inside a substation
=====================================================

A lightsim2grid substation is, by itself, a set of ``n_busbar_per_sub`` buses that elements can
be moved between (see :ref:`bus-labelling`): it says nothing about *how* an element gets from one
bus to another. Real substations do that with switches -- breakers, disconnectors, couplers --
and pypowsybl describes them in its **node-breaker view**: a voltage level is a set of
*connectivity nodes*, some carrying a *busbar section*, some the terminal of an element, joined
by *switches* (and by *internal connections*, links that are always closed).

Since version 1.0.1 an ``LSGrid`` can carry that description -- its **detailed topology** --
and, more importantly, act on it: operating a switch moves the elements the way the switches
say, exactly the way pypowsybl would.

.. contents::
    :local:
    :depth: 1

Reading the switches
----------------------

Only the pypowsybl loader can read them (the pandapower, MATPOWER and PowerModels formats have no
switches). It is off by default:

.. code-block:: python

    import pypowsybl.network as pn
    from lightsim2grid.network import init_from_pypowsybl

    net = pn.create_four_substations_node_breaker_network()
    grid = init_from_pypowsybl(net, detailed_topology=True)

    grid.has_detailed_topology()      # True
    len(grid.get_switches())          # 59, one per pypowsybl switch
    len(grid.get_busbar_sections())   # 6

``detailed_topology`` is ``True``, ``False`` (the default) or ``"auto"`` (on iff the grid has a
switch or a node-breaker voltage level). With ``False`` the grid is exactly what it was before
this feature existed, down to the bus numbering.

A **node-breaker** voltage level is read as is: its nodes, busbar sections, switches and
internal connections. A **bus-breaker** voltage level (the IEEE test cases, most files exported
from bus-branch tools) has no such description, so it gets the one IIDM itself uses to map a
bus-breaker view onto a node-breaker one: one node per bus-breaker bus, each carrying a busbar
section named after the bus, one node per element terminal, joined to its bus by a breaker
that is open iff the terminal is disconnected. The voltage level's own switches, if any, join
bus nodes.

The loader also checks itself: once the switches are declared, it projects them onto the
elements (see below) and requires that **nothing moves** -- the model's reading of the switches
must put every element exactly where pypowsybl's bus view had put it. A grid where the two
disagree is refused with the list of elements concerned.

.. note::
    Not compatible (yet) with ``buses_for_sub=True``, ``fuse_zero_impedance_branches`` or
    ``convert_dangling_lines``; see the changelog's ``[TODO]`` list for what else is deferred.

What the switches decide
--------------------------

Given the switch positions, each substation's nodes fall into connected components (joined by
closed switches and internal connections). A component is an **electrical bus** by pypowsybl's
own rule (powsybl-core's ``Networks.isBusValid``, checked row by row against
``Network.get_buses()``):

.. code-block:: text

    bus  <=>  (busbar sections >= 1 and feeders >= 1)  or  (branch ends >= 1 and feeders >= 2)

where every element terminal is a *feeder* and a *branch end* is the end of a line, of a
transformer, or an HVDC converter station. So:

======================================================  =========
component                                               a bus?
======================================================  =========
a busbar section alone                                  no
a busbar section + one load, or + one line end          yes
two line ends (no busbar section)                       yes
a line end + a load; a transformer end + a shunt        yes
a single line end (its breaker open)                    no
a load + a generator; two loads; an SVC + a load        no
an HVDC converter station + a load                      yes
======================================================  =========

Every element on a component that is a bus is connected to it; every element on a component
that is not is **disconnected** (a lone line end behind an open breaker, a load whose
disconnector is open, ...). Buses are numbered within the substation, 1-based like the local
bus id: the components holding a busbar section first, in busbar-section order, then the others
by the lowest node one of their terminals stands on. The numbering is a function of the switch
positions only, so the same positions always give the same buses.

.. note::
    Because of this numbering, the local bus ids of a voltage level loaded with
    ``detailed_topology=True`` follow the busbar sections rather than the sorted bus names,
    and ``n_busbar_per_sub`` is inferred as the **most buses any switch configuration can make**
    in a voltage level (busbar sections plus branch ends, capped by the feeders), so that no
    switch configuration is ever refused. That can be well above the number of buses in
    service today, and ``total_bus()`` (hence ``get_V()``, ``get_bus_status()``, ...) grows
    accordingly. An explicit ``n_busbar_per_sub`` below that is refused.

Operating a switch
--------------------

.. code-block:: python

    coupler = [sw for sw in grid.get_switches() if sw.name == "S1VL2_COUPLER"][0]
    grid.set_switch_open(coupler.id, True)     # the two sections of S1VL2 are now two buses
    V = grid.ac_pf(V0, 10, 1e-8)

    # the bulk form, one entry per switch (grid-wide ids)
    has_changed = np.zeros(len(grid.get_switches()), dtype=bool)
    new_open = np.zeros(len(grid.get_switches()), dtype=bool)
    has_changed[coupler.id] = True
    new_open[coupler.id] = False               # close it again
    grid.update_switches(has_changed, new_open)

:func:`~lightsim2grid.network.LSGrid.set_switch_open` and
:func:`~lightsim2grid.network.LSGrid.update_switches` move the switches, then **project** the
substations whose switches moved: their components are labelled again and every terminal there
goes where its component says -- through the same mutators a topology-vector action uses
(:func:`~lightsim2grid.network.LSGrid.update_topo`), so the change flags the solvers read, the
per-bus element counts and the powerflow are exactly what that action would have produced. A
terminal already where it should be costs nothing; a substation whose switches did not move is
not touched. An internal connection cannot be operated.

:func:`~lightsim2grid.network.LSGrid.project_switches` projects every substation at once.

The direct mutators (``deactivate_load``, ``change_bus_gen``, ``update_topo``, ...) stay usable
on a grid with a detailed topology. They are the escape hatch: they change the elements without
telling the switches, and the next projection of that substation puts its elements back where
the switches say. The other substations are untouched.

.. warning::
    Two-sided elements keep the rule they already have: by default a line or transformer whose
    one end loses its bus goes off on **both** ends, as with a topology-vector action
    (``synch_status_both_side``). pypowsybl keeps the other end connected (a *half-open*
    branch): load the grid with ``keep_half_open_lines=True`` to get that behaviour, in
    which case the open end is Kron-reduced out and the projection follows pypowsybl end
    by end. HVDC lines are never synched.

Reading the result
--------------------

Every ``*Info`` object carries the node its terminal stands on (``node_id``, ``node1_id`` /
``node2_id``; ``-1`` on a grid without detailed topology), and
:func:`~lightsim2grid.network.LSGrid.get_node_bus` gives the bus of every node (global ids, see
:ref:`bus-labelling`; ``-1`` for a node in no bus). Each
:class:`~lightsim2grid.elements.BusbarSectionInfo` says whether its section is part of a bus
(``connected``, ``bus_id``) and, after a powerflow, that bus' voltage (``res_v_kv``,
``res_theta_deg``). :func:`~lightsim2grid.network.LSGrid.get_substation_topology` exposes one
voltage level's nodes and switches by local id, for debugging.

The detailed topology is part of the grid's state: it survives ``copy()``, pickling and the
:ref:`binary format <binary-serialization>` (format 10), and a restored grid operates its
switches exactly like the original.

Building one by hand
----------------------

A grid built from scratch declares its switches with
:func:`~lightsim2grid.network.LSGrid.init_detailed_topology` (flat arrays: nodes per substation,
busbar sections and switches with their substation and local nodes, sorted by substation), then
tells each container which node its elements stand on (``set_load_to_node_id``,
``set_line_to_node1_id``, ... after the ``set_*_to_subid`` setters: a node is local to a
substation), then calls :func:`~lightsim2grid.network.LSGrid.project_switches`. See
``lightsim2grid/tests/test_detailed_topology.py`` for a complete example.

The grid2op side
------------------

:class:`~lightsim2grid.lightSimBackend.LightSimBackend` passes ``detailed_topology`` through its
``loader_kwargs`` (pypowsybl loader only), so the backend's grid carries the switches, and
requires that the buses in service fit grid2op's ``n_busbar``. Nothing in the backend acts on
them yet: switch actions, and grid2op's own ``DetailedTopoDescription``, are a later step (see
the changelog's ``[TODO]`` list).
