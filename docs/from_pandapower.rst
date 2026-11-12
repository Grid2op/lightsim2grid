.. currentmodule:: lightsim2grid.network

.. _from_pandapower:

Loading a pandapower grid
=========================

:func:`~lightsim2grid.network.init_from_pandapower` builds an :class:`~lightsim2grid.network.LSGrid` from a ``pandapowerNet``. The aim
is the grid ``pp.runpp`` solves: the same buses, admittances and injections, so that both give the same voltages.
This page says what is converted and how, what is refused, and where the two still differ.

The comparison is a test, run element by element on small grids: ``lightsim2grid/tests/test_pp_converter_elements.py``
for the injections and the slack, ``lightsim2grid/tests/test_pp_converter_trafo_vs_runpp.py`` for the
transformers. Both check the solved voltages against ``pp.runpp`` far below any engineering tolerance.

.. note::

    :class:`~lightsim2grid.lightSimBackend.LightSimBackend` uses this converter for its ``"pandapower"``
    loader. ``pp_orig_file`` selects the transformer formulas: ``"pandapower_v3"`` follows pandapower 3,
    ``"pandapower_v2"`` the older ones, which the grids of grid2op environments released before 2026 were built
    with (see :ref:`from_pandapower_trafo`).

What each table becomes
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - pandapower
     - lightsim2grid
   * - ``bus``
     - a bus, deactivated when out of service. The index does not have to be contiguous.
   * - ``line``
     - a pi model: ``r``, ``x`` from ``*_ohm_per_km * length_km``, the shunt admittance from ``g_us_per_km`` and
       ``c_nf_per_km``, half on each side, in pu of the from-bus voltage.
   * - ``trafo``
     - a transformer: see :ref:`from_pandapower_trafo`.
   * - ``load``
     - a load drawing ``p_mw * scaling``, ``q_mvar * scaling``. The voltage dependent part
       (``const_z_*`` / ``const_i_*``) is **not modeled**: such a load draws a constant power, with a warning.
   * - ``sgen``
     - a static generator injecting ``p_mw * scaling``, ``q_mvar * scaling`` (a PQ injection).
   * - ``gen``
     - a generator injecting ``p_mw * scaling`` at the voltage ``vm_pu``. ``min_q_mvar`` / ``max_q_mvar`` are kept
       but not enforced, as with ``pp.runpp``'s default (see :ref:`from_pandapower_options`).
   * - ``storage``
     - a storage unit drawing ``p_mw * scaling``, ``q_mvar * scaling`` (load convention, as pandapower).
   * - ``shunt``
     - a shunt drawing ``step * p_mw`` and ``step * q_mvar`` at its ``vn_kv``, that is times
       ``(vn_bus / vn_kv) ** 2`` at 1 pu, or the row of its ``step`` in ``shunt_characteristic_table``. Its steps
       become sections (:ref:`outer-loops-pandapower`).
   * - ``ext_grid``
     - a slack, when in service: see :ref:`from_pandapower_slack`.
   * - ``dcline``
     - a DC line: active power, losses (``loss_percent``, ``loss_mw``) and the voltage of each end.
   * - ``switch``
     - **ignored**, with a warning: an open switch does not open the line, transformer or bus coupling it is on.
   * - ``controller``
     - only the tap and shunt controllers, as data for the outer loops (:ref:`outer-loops-pandapower`).

Refused, with an error naming the table: three-winding transformers (``trafo3w``), motors, asymmetric loads and
static generators, impedances, wards and extended wards, svc, tcsc, ssc, the vsc tables and the DC grid
(``bus_dc``, ``line_dc``, ``load_dc``, ``source_dc``), and ``parallel != 1`` on anything but a transformer.

.. _from_pandapower_trafo:

Transformers
------------

The converter follows pandapower's ``_calc_tap_from_dataframe`` and ``_calc_r_x_y_from_dataframe``.

**The tap.** It starts from the transformer's rated voltages ``vn_hv_kv``, ``vn_lv_kv`` -- which may differ
from its buses' -- and its ``shift_degree``. Each tap changer (``tap_*``, then ``tap2_*`` if the columns exist)
with a ``tap_side`` of ``"hv"`` or ``"lv"`` then moves the voltage ``vn`` of that side, ``n = tap_pos -
tap_neutral`` steps away from neutral:

- ``"Ratio"`` or ``"Symmetrical"``: ``V = vn + du . exp(j . tap_step_degree)`` with
  ``du = vn . tap_step_percent / 100 . n``; ``vn`` becomes ``|V|`` and the angle of ``V`` is added to the shift;
- ``"Ideal"`` (a phase shifter): the shift grows by ``n . tap_step_degree``, or by
  ``2 . arcsin(n . tap_step_percent / 200)`` without a step angle; the voltages do not move;
- a tabular one (``tap_dependency_table``): ``vn`` is multiplied by the ``voltage_ratio`` of the row of
  ``tap_pos`` in ``trafo_characteristic_table``, its ``angle_deg`` is added to the shift, and ``vk_percent``,
  ``vkr_percent`` come from the same row.

A tap on the low voltage side subtracts its angle instead of adding it. A tap without a side, or with no
``tap_changer_type`` with pandapower 3 data, is ignored with a warning, as pandapower ignores it.

**The ratio** is then ``(vn_trafo_hv / vn_trafo_lv) / (vn_hv_bus / vn_lv_bus)``, which is not 1 for a transformer
rated off its buses' voltages even at neutral.

**The impedance** is pandapower's, in pu of the low voltage bus: ``vk_percent`` / ``vkr_percent`` on the
transformer's ``sn_mva``, scaled by ``(vn_trafo_lv / vn_lv_bus) ** 2 / parallel``, and the magnetising admittance
from ``pfe_kw`` and ``i0_percent``, scaled by the inverse. With the T model (the default) the leakage impedance is
split between the windings by ``leakage_resistance_ratio_hv`` / ``leakage_reactance_ratio_hv`` (0.5 when
absent) and the result converted to a pi model, which gives each side its own shunt admittance. The pi model
is used as is when ``net._options["trafo_model"]`` is ``"pi"``.

``pp_orig_file="pandapower_v2"`` uses pandapower 2's magnetising admittance instead, refuses the pi model, and
splits the leakage impedance half and half (another split is replaced, with a warning).

**The tap is a state.** A transformer whose ``tap_min`` .. ``tap_max`` are integers gets a tap changer holding
the ratio, shift and impedance of each position, all computed as above, so moving the tap gives the transformer
a fresh conversion would give at that position: a ratio tap changer when only the ratio moves, a phase tap
changer when the shift moves too (an angle step, an ideal phase shifter, a table with angles). The second tap
changer stays at its position. See :ref:`outer-loops-pandapower` for what the outer loops do with them.

.. _from_pandapower_slack:

Slack
-----

- The generators with ``slack=True`` are the slack, sharing the imbalance by ``slack_weight``. With
  ``pp_orig_file="pandapower_v3"``, the other generators with a non-zero ``slack_weight`` share it too. The
  ext_grids are then ignored, with a warning.
- Otherwise each ext_grid in service is a slack: the generators on its bus if every ext_grid bus has one, else a
  generator added on its bus at its ``vm_pu``. Several of them share the imbalance by the ext_grids'
  ``slack_weight``.

pandapower's default treats several ext_grids differently: each keeps its own voltage and angle, and none
shares anything. With several of them, the converted grid solves to ``pp.runpp(net, distributed_slack=True)``,
not to ``pp.runpp(net)``.

**The angle reference** is the angle the slack bus has in the initial voltage handed to the algorithm
(``LSGrid.ac_pf``): ``va_degree`` of the ext_grid is not read. Start from it to get pandapower's angles.

.. _from_pandapower_options:

``pp.runpp`` options the converted grid does not follow
-------------------------------------------------------

The converted grid is one fixed problem; the options of a ``pp.runpp`` call are not part of it:

- ``calculate_voltage_angles=False``: the transformer phase shifts are always applied;
- ``enforce_q_lims=True``: the default algorithms do not enforce reactive limits; the ``ReactiveLimits`` outer
  loop of the ``NROuter_*`` algorithms does (:ref:`outer_loops`);
- ``distributed_slack``: see :ref:`from_pandapower_slack`;
- ``trafo_model``: read from ``net._options`` when a previous ``pp.runpp`` left it there, ``"t"`` otherwise;
- ``init`` and the tolerances: those of the algorithm, given when solving.

Not converted either: the thermal limits (``max_i_ka``, ``max_loading_percent``), which grid2op reads on its own
side, and the controllers other than the tap and shunt ones.
