.. # with overline, for parts
.. * with overline, for chapters
.. = for sections
.. - for subsections
.. ^ for subsubsections
.. " for paragraphs


Comparison with pypowsybl default load-flow
============================================

In this section of the documentation we attempt to compare lightsim2grid 
and the default implementation of pypowsybl (which is OLF - `Open Load Flow <https://github.com/powsybl/powsybl-open-loadflow>__`)

All the tests were conducted on the same laptop and on publically available grid:

- ieee 9 bus
- ieee 14 bus
- ieee 30 bus
- ieee 57 bus
- ieee 118 bus
- ieee 300 bus

In all cases, the lightsim2grid `gridmodel` (lightsim2grid internal
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

- the parameters used to initialize the lightsim2grid `gridmodel`
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

The parameters used to compute the powerflow in these examples are:

.. code-block:: python

    import pypowsybl.loadflow as pypow_lf

    params = pypow_lf.Parameters(
    voltage_init_mode=pypow._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
    transformer_voltage_control_on=False,
    use_reactive_limits=False,
    phase_shifter_regulation_on=False,
    twt_split_shunt_admittance=True,
    shunt_compensator_voltage_control_on=False,
    read_slack_bus=False,
    write_slack_bus=True,
    distributed_slack=False,
    dc_use_transformer_ratio=True,
    hvdc_ac_emulation=False,
    dc_power_factor=1.,
    provider_parameters={
        "useActiveLimits": "false",
        "useReactiveLimits": "false",
        "svcVoltageMonitoring": "false",
        "voltageRemoteControl": "false",
        "writeReferenceTerminals": "false",
        "slackBusSelectionMode" : "NAME",
        "slackBusesIds" : "VL69_0",  # DEPENDS ON CASE_NAME: for case 118
        "voltagePerReactivePowerControl": "false",
        "generatorReactivePowerRemoteControl": "false",
        "secondaryVoltageControl": "false",
        }
    )

.. important::
    As you notice from these parameters, a lot of the
    simulation capacity of pypowsybl are switched off when using lightsim2grid.

.. note::
    If you are interested in an "abalation study" on the impact of certain parameters
    above, let us know, for example with a github issue or by reaching out on discord.

Results
-----------------------------------

The benchmarks were run on:

- date: 2026-01-09 10:46  CET
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.0.2
- pandas version: 2.3.3
- pandapower version: 3.2.1
- pypowsybl version: 1.13.0
- grid2op version: 1.12.2
- lightsim2grid version: 0.12.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
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
ieee9       1.29e-01         3.56e+00
ieee14      1.75e-01         3.98e+00 
ieee30      2.96e-01         3.92e+00 
ieee57      4.77e-01         5.44e+00
ieee118     6.74e-01         6.12e+00
ieee300     2.51e+00         1.28e+01
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
ieee9       1.71e-02         7.26e-01
ieee14      2.83e-02         8.95e-01 
ieee30      6.00e-02         1.26e+00 
ieee57      1.41e-01         1.46e+00
ieee118     3.11e-01         2.48e+00
ieee300     1.76e+00         5.78e+00
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
ieee9       2.29e-02         3.03e-01
ieee14      2.92e-02         1.97e-01
ieee30      4.85e-02         1.68e-01
ieee57      1.32e-01         1.81e-01
ieee118     2.05e-01         3.38e-01
ieee300     9.95e-01         1.31e+00
========== =============== ===============
