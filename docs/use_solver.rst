Use as Pandapower Solver
=========================
LightSim2grid can be used as a specific implementation of the pandapower "newtonpf" function.

Suppose you somehow get:

- `Ybus` the admittance matrix of your powersystem, for example given by pandapower
  (will be converted to a scipy `sparse.csc_matrix` )
- `V0` the (complex) voltage vector at each bus, for example  given by pandapower
- `Sbus` the (complex) power absorb at each bus, for example  as given by pandapower
- `ppci` a ppc internal pandapower test case (for compatibility with  pandapower, currently not used)
- `pv` list of PV buses
- `pq` list of PQ buses
- `options` list of pandapower "options" (or dictionary with keys `max_iteration` and `tolerance_mva`)

You can define replace the `newtonpf` function of `pandapower.pandapower.newtonpf` function with the following
piece of code:

.. code-block:: python

    from lighsim2grid.newtonpf import newtonpf
    V, converged, iterations, J = newtonpf(Ybus, V, Sbus, pv, pq, ppci, options)

This function uses the KLU algorithm and a c++ implementation of a Newton solver for speed.


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`