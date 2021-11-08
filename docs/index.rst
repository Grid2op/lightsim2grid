.. LightSim2Grid documentation master file, created by
   sphinx-quickstart on Tue Oct  6 09:41:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LightSim2Grid's documentation!
=========================================
LightSim2Grid is a simulator that can be used to accelerate the computation (leading to approximately the
same results) when using the `Grid2Op <https://github.com/rte-france/Grid2Op>`_ platform.

It is a port in c++ of some part of the excellent `PandaPower <https://github.com/e2nIEE/pandapower>`_ package
to make the computation faster.

As from version 0.3.0:

- `LightSimBackend` is available on windows (you don't need to compile KLU anymore)
- `LightSimBackend` is compatible with the "pickle" library, so you can use it with multiprocessing
  on windows, and with frameworks that requires it (for example ray rllib).

.. toctree::
   :maxdepth: 2
   :caption: Quickstart

   quickstart
   disclaimer
   use_with_grid2op
   use_solver


Technical Documentation (work in progress)
-------------------------------------------

This is a work in progress at the moment

.. toctree::
   :maxdepth: 2
   :caption: Technical Documentation

   gridmodel
   lightsimbackend
   solvers
   time_series


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
