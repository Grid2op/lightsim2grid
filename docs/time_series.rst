Time Series
=======================================

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
to / from other data sources (for now please use the more basic :class:`lightsim2grid.timeSerie.TimeSeriesCPP` for such purpose)

Importantly, this method is around **11x** faster than simulating "do nothing" (or "one change then nothing") with grid2op
(and lightsim2grid, see section :ref:`timeserie_benchmark` )

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
    `examples\\timeseries_with_grid2op_multithreading.py` file.

.. note::

    Set `time_series.init_from_n_powerflow = True` **before** calling `compute_V` (setting it
    afterwards has no effect) to initialize the first step of the batch with the voltage
    solution of the grid's current ("n") state instead of a flat start -- this is usually
    faster. See :func:`lightsim2grid.timeSerie.TimeSerie.init_from_n_powerflow`.

Independent scenarios: `InjectionSweep`
----------------------------------------

`TimeSerie` initializes each step with the solution of the step before it. That is the right
thing to do for a *time* series -- two consecutive instants are close to one another, so the
previous solution is an excellent starting point -- but it makes the steps a chain: step `i`
cannot be computed before step `i-1`, and its result depends on it.

If your "steps" are unrelated scenarios rather than consecutive instants (a sample of load /
generation patterns, a set of "what if" injections, a Monte-Carlo draw...), use
:class:`lightsim2grid.injectionSweep.InjectionSweep` instead. It computes exactly the same
thing, with exactly the same interface, but starts every step from the same voltage -- the
one you provide, or the "n" powerflow result if `init_from_n_powerflow` is set -- just like
:class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis` does for each contingency.

Two things follow:

- the results do not depend on the order the steps are given in;
- the batch can be spread over several OS threads, in c++, with `sweep.nb_thread = ...`
  (the results do not depend on that either). `TimeSerie` cannot: splitting a chain into
  per-thread ranges would break it, so setting `nb_thread` to anything but 1 raises there.

.. code-block:: python

    import grid2op
    from lightsim2grid import InjectionSweep, LightSimBackend

    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    sweep = InjectionSweep(env)
    sweep.nb_thread = 4
    res_p, res_a, res_v = sweep.get_flows(scenario_id=..., seed=...)

.. note::

    On a genuine time series `TimeSerie` usually needs fewer solver iterations, precisely
    because each step starts from its neighbour's solution. `InjectionSweep` trades that
    away for independence -- and buys back much more than it costs as soon as you use
    several threads.

.. note::

    `TimeSerie` / `InjectionSweep` also expose a newer, setter-based API
    (`modify_gen_p` / `modify_sgen_p` / `modify_load_p` / `modify_load_q` /
    `modify_gen_v` + `compute()`) as an alternative to the single bundled
    `compute_V_from_inj` call -- see :doc:`scenario_sweep` (which uses that same API,
    plus a per-step contingency) for details on how it behaves when an axis is never
    set. `modify_gen_v` (per-step generator target voltage magnitude, in pu -- `vm_pu`,
    NOT kV) is different from the other four: it does not feed the injection (`Sbus`)
    at all, it only re-seeds `|V|` at each voltage-regulating generator's regulated
    bus before that step's solve.

.. _timeserie_benchmark:

.. _ts_benchmarks:

Benchmarks (Time Series)
-------------------------

Here are some benchmarks made with:

- date: 2026-08-28 16:57  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.3.5
- pandas version: 2.3.3
- pandapower version: 3.4.0
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.0.0
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 time_serie.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2
    Total time spent in "computer" to solve everything: 0.03s (20834 pf / s), 0.05 ms / pf)
        - time to pre process the injections: 0.00s
        - time to perform powerflows: 0.03s (22437 pf / s, 0.04 ms / pf)
    In addition, it took 0.00 s to retrieve the current from the complex voltages (in total 19934.2 pf /s, 0.05 ms / pf)

    Comparison with raw grid2op timings
    It took grid2op (with lightsim2grid): 0.32s to perform the same computation
        This is a 11.0 speed up from TimeSerie over raw grid2op (lightsim2grid)
    It took grid2op (with pandapower): 6.56s to perform the same computation
        This is a 227.1 speed up from TimeSerie over raw grid2op (pandapower)

    In this case then, the `TimeSerie` module is 11 times faster than raw grid2op (lightsim2grid) and 227 times faster than raw grid2op (pandapower)
    All results match !


Detailed usage
--------------------------

.. automodule:: lightsim2grid.timeSerie
    :members:
    :autosummary:

.. automodule:: lightsim2grid.injectionSweep
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
