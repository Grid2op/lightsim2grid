Available "solvers" (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Type of solvers available
##########################

For now, lightsim2grid ships with at most 3 fully working (and tested) solvers:

- **LS+GS** (LightSimBackend+Gauss Seidel): the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows.
- **LS+GS S** (LightSimBackend+Gauss Seidel Synchronous): the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows.
- **LS+SLU** (Newton Raphson+SparseLU): the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform).
- **LS+KLU** (Newton Raphson+KLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "KLU" from the `SuiteSparse` c package implemented.
- **LS+NICSLU** (Newton Raphson+NICSLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "NICSLU". [**NB** NICSLU is a free software but not open source, in order to use
  it with lightsim2grid, you need to check section 
  [It is required to install lightsim2grid from source for such solver and following the 
  Readme for instructions on how to compile with such solver]

Usage
############
In this section we briefly explain how to switch from one solver to another. An example of code using this feature
is given in the
`"benchmark_solvers.py" <https://github.com/BDonnot/lightsim2grid/blob/master/benchmarks/benchmark_solvers.py>`_
script available in the `"benchmarks" <https://github.com/BDonnot/lightsim2grid/tree/master/benchmarks/>`_
directory of the lightsim2grid repository.

A concrete example of how to change the solver in this backend is:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend

    # create an environment
    env_name = "l2rpn_case14_sandbox"
    env_lightsim = grid2op.make(env_name, backend=LightSimBackend())

    # retrieve the available solver types
    available_solvers = env_lightsim.backend.available_solvers

    # change the solver types (for example let's use the Gauss Seidel algorithm)
    env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidel)

    # customize the solver (available for all solvers)
    env_lightsim.backend.set_solver_max_iter(10000)  # all solvers here are iterative, this is the maximum number of iterations
    env_lightsim.backend.set_tol(1e-7)  # change the tolerance (smaller tolerance gives a more accurate results but takes longer to compute)
    # see the documentation of LightSimBackend for more information



Detailed usage
###############
TODO examples on how to import, and documentation of main methods


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`