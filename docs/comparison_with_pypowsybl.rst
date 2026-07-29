.. # with overline, for parts
.. * with overline, for chapters
.. = for sections
.. - for subsections
.. ^ for subsubsections
.. " for paragraphs


Comparison with pypowsybl default load-flow
============================================

In this section of the documentation we attempt to compare lightsim2grid 
and the default implementation of pypowsybl (which is OLF - `Open Load Flow <https://github.com/powsybl/powsybl-open-loadflow>`__)

All the tests were conducted on the same laptop and on publically available grid:

- ieee 9 bus
- ieee 14 bus
- ieee 30 bus
- ieee 57 bus
- ieee 118 bus
- ieee 300 bus

In all cases, the lightsim2grid `LSGrid` (lightsim2grid internal
representation of a powergrid) were initialized from the pypowsybl grid.

Disclaimer
-----------

Compared to pypowsybl, lightsim2grid has only a very (very) limited number of possible
behaviour.

Pypowsybl is likely to be more accurate (if parametrized properly) for industrial grid 
as it can emulate the behaviour of much more elements in much more detail than lightsim2grid.

For example, lighsim2grid does not look at the "reactive power limits" of generators, 
pypowsybl (open load-flow) can. LightSim2Grid does not change the tap ratio of 
any transformers during computation, pypowsybl is perfectly able to do that, when 
the slack bus is distributed in lightsim2grid, lightsim2grid does not check whether 
or not the generators can produce / absorb the the active power they are supposed to, 
pypowsybl is able to dynamically meet this criteria etc. etc.

Also, we want to note that the comparison here will limit to computations that both
lightsim2grid and pypowsybl are able to perform. At time of writing, pypowsybl could
do many more things than lightsim2grid.

.. important::
    The overall message of this page is not to show that lightsim2grid should be
    prefered to pypowsybl. 

    Its goal is rather to explain how to get consistent results between pypowsybl 
    and lightsim2grid.


Methodology
------------

The results will show the difference between pypowsybl and lightsim2grid when running
the same simulation (an AC powerflow) on the same grid.

It will expose:

- the parameters used to initialize the lightsim2grid `LSGrid`
- the parameters used to run the powerflow computation with pypowsybl
- the time it takes to perform these powerflows in different settings
- the mismatch of the voltage angle (in radian) 
  and the voltage angle (in per unit) at each bus of the grid (average and max value)

Reproduce the results
************************

You can run the example by running the script:

.. code-block:: bash

    cd benchmarks
    python compare_lightsim2grid_pypowsybl.py --case $CASE_NAME

For example:


.. code-block:: bash

    cd benchmarks
    python compare_lightsim2grid_pypowsybl.py --case ieee9


Load-flow parameters
**********************

The parameters used to compute the powerflow in these examples are the canonical
"every outer loop disabled" parameters shipped with lightsim2grid, exported as
``get_pypowsybl_loopfree_parameters``. lightsim2grid solves a single power-flow
problem (no outer loops), so to get consistent results pypowsybl must be told to
run no outer loop either:

.. code-block:: python

    from lightsim2grid.network import get_pypowsybl_loopfree_parameters

    # pin the slack on the same bus lightsim2grid uses (depends on the case,
    # e.g. "VL69_0" for ieee118); omit slack_bus_ids to read the slack from
    # the network instead.
    params = get_pypowsybl_loopfree_parameters(slack_bus_ids="VL69_0")

Under the hood this builds a :class:`pypowsybl.loadflow.Parameters` that disables
distributed slack, reactive limits, transformer / shunt / phase-shifter voltage
control, area-interchange and secondary-voltage control, automation systems, etc.
The key mechanism is the empty ``outerLoopNames`` allow-list, which registers
*zero* outer loops regardless of which loops a future pypowsybl release adds.

If you need a *raw* (un-baked) network whose outer loops would actually trigger
to also match lightsim2grid, freeze the converged outer-loop state first with
``lightsim2grid.network.bake_outer_loops`` before solving loop-free.

.. important::
    As you notice from these parameters, a lot of the
    simulation capacity of pypowsybl are switched off when using lightsim2grid.

.. note::
    If you are interested in an "abalation study" on the impact of certain parameters
    above, let us know, for example with a github issue or by reaching out on discord.


Comparing the two engines
**************************

To check that lightsim2grid reproduces an OLF operating point you can use the
``compare_baked`` helper. It solves the network *with* outer loops in OLF, bakes
the converged outer-loop state into the inputs (so the problem becomes a plain
power flow), optionally applies the same outages to both engines, then solves
loop-free in OLF and in lightsim2grid and compares the bus voltages:

.. code-block:: python

    from lightsim2grid.network import compare_baked
    import pypowsybl as pp

    res = compare_baked(
        pp.network.create_ieee14,   # a callable returning a fresh network
        slack_gen_id="B1-G",
        line_outages=["L1-2-1"],    # optional, applied to both engines
    )
    print(res)                      # ComparisonResult(max |dV| = ..., ...)
    print(res.max_dvm_pu)           # largest |Vmag| mismatch (pu)
    print(res.table)                # per-bus detail

The call returns a :class:`lightsim2grid.network.ComparisonResult` summarising the
largest voltage-magnitude and voltage-angle mismatches (plus a per-bus table):

.. autoclass:: lightsim2grid.network.ComparisonResult
    :members:
    :no-index:

.. autofunction:: lightsim2grid.network.compare_baked
    :no-index:

Inspecting results in a pypowsybl-like way
********************************************

``compare_baked`` above only compares bus voltages. If you want to inspect (or write
generic analysis code against) the *full* result of a lightsim2grid powerflow -- lines,
transformers, generators, loads, shunts, HVDC lines, ... -- with the exact same API and
DataFrame shape as a solved pypowsybl :class:`pypowsybl.network.Network`, use
:class:`lightsim2grid.network.LightsimResultNetwork`. It wraps a converged ``LSGrid``
(built by ``init_from_pypowsybl``) and the pypowsybl network it was built from, and
exposes ``get_buses`` / ``get_lines`` / ``get_2_windings_transformers`` /
``get_generators`` / ``get_loads`` / ``get_shunt_compensators`` /
``get_static_var_compensators`` / ``get_batteries`` / ``get_hvdc_lines`` /
``get_vsc_converter_stations`` / ``get_lcc_converter_stations``, each returning a pandas
DataFrame indexed by the pypowsybl element id, with pypowsybl's own column names and
sign conventions (post-solve ``p`` / ``q`` in the load convention):

.. code-block:: python

    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl, LightsimResultNetwork

    net = pp.network.create_ieee14()
    grid = init_from_pypowsybl(net, gen_slack_id="B1-G")
    V = grid.ac_pf(..., 10, 1e-7)
    assert len(V) > 0  # converged

    res_net = LightsimResultNetwork(grid, net)
    res_net.get_generators()  # same shape/columns as net.get_generators() after a solve
    res_net.get_lines(attributes=["p1", "q1", "i1"])

.. autoclass:: lightsim2grid.network.LightsimResultNetwork
    :members:
    :no-index:

.. warning::
    Only powerflow *results* (post-solve quantities such as ``p`` / ``q`` / ``i`` / bus
    voltage) are actually read from the solved ``LSGrid`` and mapped back onto
    pypowsybl's DataFrame shape. Everything else in the returned DataFrame -- topology,
    ratings, static per-element metadata, ... -- is read from the original pypowsybl
    ``net``, not recomputed by lightsim2grid.

.. note::
    ``LightsimResultNetwork`` is a convenience wrapper for pypowsybl-shaped analysis
    code, not the fast path: building each DataFrame has real overhead (id alignment,
    column renaming, sign-convention conversion) on top of the powerflow itself. If you
    only need raw arrays for performance-sensitive code, use lightsim2grid's own
    accessors on ``LSGrid`` directly instead.

Results
-----------------------------------

The benchmarks were run on:

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

The results were obtained by launching:

.. code-block:: bash

    python compare_lightsim2grid_pypowsybl.py --case ieee9
    python compare_lightsim2grid_pypowsybl.py --case ieee14
    python compare_lightsim2grid_pypowsybl.py --case ieee30
    python compare_lightsim2grid_pypowsybl.py --case ieee57
    python compare_lightsim2grid_pypowsybl.py --case ieee118
    python compare_lightsim2grid_pypowsybl.py --case ieee300

And formatting the results in the table below.


Precision of lightsim2grid
*****************************

On average (across all buses) the errors were:

========== ============= ===============
case name   angle (rad)  magnitude (pu)
========== ============= ===============
ieee9       1.82e-08        1.15e-08
ieee14      9.70e-10        1.27e-09 
ieee30      1.58e-09        3.55e-09 
ieee57      1.63e-07        2.71e-07
ieee118     1.06e-07        3.15e-09
ieee300     3.10e-07        1.75e-08
========== ============= ===============

Maximum error, for all buses:

========== ============= ===============
case name   angle (rad)  magnitude (pu)
========== ============= ===============
ieee9       3.35e-08        2.65e-08
ieee14      2.35e-09        2.92e-09 
ieee30      3.23e-09        7.96e-09 
ieee57      9.54e-07        1.20e-06
ieee118     2.54e-07        6.92e-08
ieee300     3.80e-07        1.59e-07
========== ============= ===============

As we can notice in the tables above, the results match up to the 
solver precisions (set to 1e-6 for lightsim2grid).

On these grids, lightsim2grid and pypowsybl give the same exact results.

Computation times (1 powerflow)
********************************

In this part, we report the time to compute the initial powerflow, right
after the initialization of the grid for both lightsim2grid and pypowsybl.

The timings reported here are measured from python using "time.perf_counter()"
before and after the computation are performed.

Only the time to perform the powerflow is measured. In particular, the time
to read back the data is excluded.

Times are expressed in ms.

========== =============== ===============
case name   lightsim2grid    pypowsybl
========== =============== ===============
ieee9       1.10e-01         4.02e+00
ieee14      1.48e-01         3.82e+00
ieee30      2.61e-01         5.19e+00 
ieee57      4.24e-01         5.28e+00 
ieee118     6.78e-01         6.11e+00
ieee300     1.73e+00         1.15e+01
========== =============== ===============

For this initial computation, lightsim2grid seems to be between 30 and 5x faster 
than pypowsybl.

.. warning::
    This is not fair for pypowsybl.

    Pypowsybl is not optimized only for speed and can simulate 
    much more complex grids with an higher fidelity, which is not 
    reported here.


Computation times (100 powerflows)
************************************

In this section, we compare the capacity of lightsim2grid and pypowsybl to 
perform successive powerflow computation when only the loads are modified.

This comparison is done when using "raw" lightsim2grid / pypowsybl code, without
trying to achieve the "best performance". Some performance gain could
be achieved with different optimizations, for example by recycling previous
results (avoiding to allocate memory, preventing copy, re use of some matrix
strucure, taking advantage of the linear solver and avoid costly
call when performing some factorization etc.)

The results in the table bellow are given in ms and report the average 
time it took to perform the 100 powerflows.

========== =============== ===============
case name   lightsim2grid    pypowsybl
========== =============== ===============
ieee9       1.64e-02         9.23e-01
ieee14      2.68e-02         1.14e+00 
ieee30      5.81e-02         1.44e+00 
ieee57      1.37e-01         1.56e+00
ieee118     3.04e-01         3.02e+00
ieee300     1.73e+00         5.71e+00
========== =============== ===============


Computation times security analysis
************************************

In this setting, we compare the time it takes to run a "contingency analysis"
by simulating, in turn, the disconnection of every lines or transformer on
the grid.

The table here is obtained by using `contingencyAnalysis` module of
lightsim2grid and the `pypowsybl.security` module from pypowsybl.

The table below provides the average time it takes to simulate the
effect of 1 contingency in ms. We don't measure the time taken to 
compute the flows from the resulting voltages.

========== =============== ===============
case name   lightsim2grid    pypowsybl
========== =============== ===============
ieee9       2.15e-02         3.47e-01
ieee14      2.72e-02         2.01e-01
ieee30      4.73e-02         1.45e-01
ieee57      1.32e-01         1.86e-01
ieee118     2.04e-01         3.84e-01
ieee300     9.91e-01         1.23e+00
========== =============== ===============
