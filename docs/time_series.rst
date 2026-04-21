Time Series (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Goal
--------------------------

This class aims to make faster (and easier) the computations of the current flows (measured in Amps)
at a certain side of a powerline / transformer when the topology is not modified.

It can be used as:

.. code-block:: python

    from lightsim2grid import TimeSerie
    import grid2op
    from lightsim2grid.lightSimBackend import LightSimBackend

    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    time_series = TimeSerie(env)
    res_p, res_a, res_v = time_series.get_flows(scenario_id=..., seed=...)

    # we have:
    # res_p[row_id] will be the active power flows (origin side), on all powerlines corresponding to step "row_id"
    # res_a[row_id] will be the current flows, on all powerlines corresponding to step "row_id"
    # res_v[row_id] will be the complex voltage, on all bus of the grid at step "row_id"

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources (for now please use the more basic :class:`lightsim2grid.timeSerie.Computers` for such purpose)

Importantly, this method is around **13x** faster than simulating "do nothing" (or "one change then nothing") with grid2op
(see section :ref:`timeserie_benchmark` )

.. note:: 

    A more detailed example is given in the 
    `examples\\time_serie.py` file from the lightsim2grid package.

.. warning:: Topology and injections
    
    The topology is taken from the initial provided grid and cannot be changed when evaluating
    a given "time serie".

    Then, the call to `time_series.compute_V(scenario_id=..., seed=...)` will only read the injections
    (productions and loads) from grid2op to compute the voltages.

.. note:: 
    
    As this class calls a long c++ function, it is possible to use the python `Threading`
    module to achieve high efficient parrallelism. An example is provided in the
    `examples\\computers_with_grid2op_multithreading.py` file.

.. _timeserie_benchmark: 

.. _ts_benchmarks:

Benchmarks (Time Series)
-------------------------

Here are some benchmarks made with:

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

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 time_serie.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2
    Total time spent in "computer" to solve everything: 0.03s (21834 pf / s), 0.05 ms / pf)
        - time to pre process the injections: 0.00s
        - time to perform powerflows: 0.02s (23675 pf / s, 0.04 ms / pf)
    In addition, it took 0.00 s to retrieve the current from the complex voltages (in total 20703.1 pf /s, 0.05 ms / pf)

    Comparison with raw grid2op timings
    It took grid2op (with lightsim2grid): 0.31s to perform the same computation
        This is a 11.3 speed up from TimeSerie over raw grid2op (lightsim2grid)
    It took grid2op (with pandapower): 6.47s to perform the same computation
        This is a 232.6 speed up from TimeSerie over raw grid2op (pandapower)
    All results match !



In this case then, the `TimeSerie` module is more than **15** times faster than raw grid2op.


Detailed usage
--------------------------

.. automodule:: lightsim2grid.timeSerie
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
