GridModel module (doc in progress)
====================================

The main class of the lightsim2grid python package is the `GridModel` class, that is a python class created
from the the c++ `GridModel` (thanks fo pybind11).

This class basically represents a powergrid (what elements it is made for, their electro technical properties etc.)

To create such class, for now the only way is to get it from a pandapower grid (and it does not model every elements there !)

For example, you can init it like:

.. code-block:: python

    from lightsim2grid.initGridModel import init
    pp_net = ...  # any pandapower grid

    lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

.. _elements-modeled: 

Elements modeled
------------------

- PandaPowerConverter
- GridModel

example on how to import, and documentation of main methods

TODO doc

Detailed documentation
######################

.. automodule:: lightsim2grid.initGridModel
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
