Check if some given (complex) Voltage meets the KCL
====================================================

This module allows to check if a given vector checks the KCL (Kirchhoff Current Law) or not.

It is a wrapper around the :func:`lightsim2grid.initGridModel.GridModel.check_solution` function.

TODO DOC in progress !

Examples
---------
See the section :ref:`use_with_g2op` for more information and more examples.

For standard grid2op environment, you can use it like:

.. code-block:: python

    import grid2op
    import numpy as np
    from lightsim2grid import PhysicalLawChecker

    # create a grid2op environment
    env_name = "l2rpn_case14_sandbox"
    env = grid2op.make(env_name, ...)

    # create the checker
    checker = PhysicalLawChecker(env)

    # get an observation
    obs = env.reset()

    # retrieve somehow a complex voltage
    v = np.zeros(2*env.n_sub, dtype=complex)
    v[...] = ...  # put here the value of the complex voltage you want to get

    # check if it meets the KCL (Kirchhoff's Current Law)
    mismatch = checker.check_solution(v, obs)
    # mistmatch has same size as v and contains the (complex) current mismatch at each bus of the grid.


Detailed usage
--------------------------

.. automodule:: lightsim2grid.physical_law_checker
    :members:
    :autosummary:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
