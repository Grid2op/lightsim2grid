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

Generators (standard)
+++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.DataGen
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.GenInfo
    :members:
    :autosummary:

Static Generators (more exotic)
++++++++++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.DataSGen
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.SGenInfo
    :members:
    :autosummary:

Loads and Storage Units
++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.DataLoad
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.LoadInfo
    :members:
    :autosummary:

Shunts
++++++++++++++++++++++++

.. autoclass:: lightsim2grid.elements.DataShunt
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.ShuntInfo
    :members:
    :autosummary:

Lines
++++++

.. autoclass:: lightsim2grid.elements.DataLine
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.LineInfo
    :members:
    :autosummary:

Transformers
+++++++++++++

.. autoclass:: lightsim2grid.elements.DataTrafo
    :members:
    :autosummary:

.. autoclass:: lightsim2grid.elements.TrafoInfo
    :members:
    :autosummary:


Detailed documentation
######################

.. automodule:: lightsim2grid.initGridModel
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
