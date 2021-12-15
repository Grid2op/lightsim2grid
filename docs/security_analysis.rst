Security Analysis (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Goal
-----------------

This class aims to make faster (and easier) the computations of a security analysis (which is the results of some 
powerflow after the disconnection of one or more powerlines)

This function is much (much) faster than its pure grid2op counterpart. For example,
on the case 118, to simulate all n-1 contingencies you can expect a **~100x** speed ups 
compared to using the grid2op `obs.simulate(..., time_step=0)` while obtaining the
exact same results (see section `Benchmarks`)

It can be used as:

.. code-block:: python

    import grid2op
    from lightsim2grid import SecurityAnalysis
    from lightsim2grid import LightSimBackend
    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    security_analysis = SecurityAnalysis(env)
    security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
    res_a, res_v = security_analysis.get_flows()

    # in this results, then
    # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
    # you can retrieve it with `security_analysis.contingency_order[row_id]`

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources.

.. note:: 
    
    A more advanced usage is given in the `examples\\security_analysis.py` 
    file from the lightsim2grid package.

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
    python3 security_analysis.py


For this setting the outputs are:

.. code-block:: bash

    For environment: l2rpn_neurips_2020_track2_small (186 n-1 simulated)
    Total time spent in "computer" to solve everything: 2.8ms (66375 pf / s), 0.02 ms / pf)
        - time to compute the coefficients to simulate line disconnection: 0.04ms
        - time to pre process Ybus: 2.00ms
        - time to perform powerflows: 0.73ms (256458 pf / s, 0.00 ms / pf)
    In addition, it took 0.67 ms to retrieve the current from the complex voltages (in total 53578.7 pf /s, 0.02 ms / pf)
    Comparison with raw grid2op timings
    It took grid2op: 0.41s to perform the same computation
    This is a 116.7 speed up from SecurityAnalysis over raw grid2op (using obs.simulate)


In this case then, the `SecurityAnalysis` module is more than **100** times faster than raw grid2op (
with obs.simulate as a way to compute the outcome of a contingency)


Detailed usage
--------------------------

.. automodule:: lightsim2grid.securityAnalysis
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`