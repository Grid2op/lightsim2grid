Time Series (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Goal
##############################

This class aims to make faster (and easier) the computations of the current flows (measured in Amps)
at a certain side of a powerline.

It can be used as:

.. code-block:: python

    from lightsim2grid import TimeSerie
    import grid2op
    from lightsim2grid.LightSimBackend import LightSimBackend

    env_name = ...
    env = grid2op.make(env_name, backend=LightSimBackend())

    time_series = TimeSerie(env)
    Vs = time_series.compute_V(scenario_id=..., seed=...)
    As = time_series.compute_A()  # will contain the flows, in amps at each step (rows) for each powerline (column)

For now this relies on grid2op, but we could imagine a version of this class that can read
to / from other data sources.

.. warning:: Topology and injections
    
    The topology is taken from the initial provided grid and cannot be changed when evaluating
    a given "time serie".

    Then, the call to `time_series.compute_V(scenario_id=..., seed=...)` will only read the injections
    (productions and loads) from grid2op to compute the voltages.


Detailed usage
###############
TODO examples on how to import, and documentation of main methods


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`