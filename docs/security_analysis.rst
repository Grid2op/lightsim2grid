Contingency Analysis (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

.. warning::
    This function might give wrong result for lightsim2grid version 0.5.5 were they were a bug : when some contingencies made the grid
    non connex, it made all the other contingencies diverge. This bug has been fixed in version 0.6.0 and this is why we **do not recommend**
    to use this feature with lightsim2grid version < 0.6.0 !
    
Goal
-----------------

This class aims to make faster (and easier) the computations of a security analysis (which is the results of some 
powerflow after the disconnection of one or more powerlines)

This function is much (much) faster than its pure grid2op counterpart. For example,
on the case 118, to simulate all n-1 contingencies you can expect a **~20x** speed ups 
compared to using the grid2op `obs.simulate(..., time_step=0)` while obtaining the
exact same results (see section `Benchmarks`)

It can be used as:

.. code-block:: python

    import grid2op
    from lightsim2grid import ContingencyAnalysis
    from lightsim2grid import LightSimBackend
    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    security_analysis = ContingencyAnalysis(env)
    security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
    res_p, res_a, res_v = security_analysis.get_flows()

    # in this results, then
    # res_p[row_id] will be the active power flows (origin side), on all powerlines corresponding to the `row_id` contingency.
    # res_a[row_id] will be the current flows, on all powerlines corresponding to step "row_id"
    # res_v[row_id] will be the complex voltage, on all bus of the grid corresponding to the `row_id` contingency.
    # you can retrieve which contingency is id'ed `row_id` with `security_analysis.contingency_order[row_id]`

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources.

.. note:: 
    
    A more advanced usage is given in the `examples\\security_analysis.py` 
    file from the lightsim2grid package.

.. _sa_benchmarks:

Benchmarks (Contingency Analysis)
----------------------------------

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

    cd examples
    python3 security_analysis.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2_small (177 n-1 simulated)
    Total time spent in "computer" to solve everything: 11.1ms (15913 pf / s), 0.06 ms / pf)
        - time to compute the coefficients to simulate line disconnection: 0.28ms
        - time to pre process Ybus: 0.30ms
        - time to perform powerflows: 10.25ms (17276 pf / s, 0.06 ms / pf)
    In addition, it took 0.50 ms to retrieve the current from the complex voltages (in total 15229.8 pf /s, 0.07 ms / pf)

    Comparison with raw grid2op timings
    It took grid2op (with lightsim2grid, using obs.simulate): 0.28s to perform the same computation
        This is a 24.2 speed up from SecurityAnalysis over raw grid2op (using obs.simulate and lightsim2grid)
    It took grid2op (with pandapower, using obs.simulate): 9.94s to perform the same computation
        This is a 855.2 speed up from SecurityAnalysis over raw grid2op (using obs.simulate and pandapower)
    All results match !


In this case then, the `SecurityAnalysis` module is more than **22** times faster than raw grid2op (
with obs.simulate as a way to compute the outcome of a contingency)


Detailed usage
--------------------------

.. automodule:: lightsim2grid.securityAnalysis
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
