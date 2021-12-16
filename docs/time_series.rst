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
    Vs = time_series.compute_V(scenario_id=..., seed=...)
    As = time_series.compute_A()  # will contain the flows, in amps at each step (rows) for each powerline (column)

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources (for now please use the more basic :class:`lightsim2grid.timeSerie.Computers` for such purpose)

Importantly, this method is around 15 times faster than simulating "do nothing" (or "one change then nothing") with grid2op
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

Benchmarks
-----------------

Here are some benchmarks made with:

- system: Linux 5.11.0-38-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.18.5
- pandas version: 1.1.4
- pandapower version: 2.6.0
- lightsim2grid version: 0.5.5
- grid2op version: 1.6.4

Where lightsim2grid has been installed from source with all optimization enabled.

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd examples
    python3 time_serie.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2_small
    Total time spent in "computer" to solve everything: 0.53s (15252 pf / s), 0.07 ms / pf)
        - time to pre process the injections: 0.03s
        - time to perform powerflows: 0.50s (16269 pf / s, 0.06 ms / pf)
    In addition, it took 0.06 s to retrieve the current from the complex voltages (in total 13666.3 pf /s, 0.07 ms / pf)
    Comparison with raw grid2op timings
    It took grid2op: 9.26s to perform the same computation
    This is a 15.7 speed up from TimeSerie over raw grid2op


In this case then, the `TimeSerie` module is more than **15** times faster than raw grid2op.


Detailed usage
--------------------------

.. automodule:: lightsim2grid.timeSerie
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
