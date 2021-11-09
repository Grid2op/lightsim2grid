Time Series (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Goal
##############################

This class aims to make faster (and easier) the computations of a security analysis (which is the results of some 
powerflow after the disconnection of one or more powerlines)

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


Detailed usage
###############
TODO examples on how to import, and documentation of main methods

.. automodule:: lightsim2grid.securityAnalysis
    :members:
    :autosummary:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`