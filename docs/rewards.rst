Custom Rewards
=======================================

This module provides grid2op reward classes that leverage lightsim2grid-specific features. For now
this is a single reward, :class:`~lightsim2grid.rewards.N1ContingencyReward`, which uses
:class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis` to penalize an agent for leaving the
grid in a state where some n-1 contingencies would exceed their current limits.

Detailed usage
--------------------------

.. automodule:: lightsim2grid.rewards
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`