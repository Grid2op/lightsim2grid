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
    above, let us know, for example on github or on discord.

Results
-----------

