.. _solvers_doc:

Available "solvers" (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Type of solvers available
--------------------------

For now, lightsim2grid ships with at most TODO DOC fully working (and tested) solvers:

- **LS+GS** [:class:`lightsim2grid.solver.GaussSeidelSolver`]: the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows.
- **LS+GS S** [:class:`lightsim2grid.solver.GaussSeidelSynchSolver`]: the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows.
- **LS+SLU** [:class:`lightsim2grid.solver.SparseLUSolver`]: the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform). This solver supports distributed slack bus.
- **LS+SLU (single)** [:class:`lightsim2grid.solver.SparseLUSolverSingleSlack`]: same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **LS+KLU** [:class:`lightsim2grid.solver.KLUSolver`]: he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "KLU" from the `SuiteSparse` C package. This solver supports distributed slack bus.
- **LS+KLU (single)** [:class:`lightsim2grid.solver.KLUSolverSingleSlack`]: same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **LS+NICSLU** [:class:`lightsim2grid.solver.NICSLUSolver`]: he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "NICSLU". [**NB** NICSLU is a free software but not open source, in order to use
  it with lightsim2grid, you need to install lightsim2grid from source for such solver]
- **LS+NICSLU (single)** [:class:`lightsim2grid.solver.NICSLUSolverSingleSlack`]: same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **LS+CKTSO** [:class:`lightsim2grid.solver.CKTSOSolver`]: he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "CKTSO". [**NB** CKTSO is a free software but not open source, in order to use
  it with lightsim2grid, you need to install lightsim2grid from source for such solver] (only tested on linux)
- **LS+CKTSO (single)** [:class:`lightsim2grid.solver.CKTSOSolverSingleSlack`]: same as above but this solver does not support distributed slack bus and
  can thus be slightly faster. (only tested on linux)

Usage
--------------------------
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

    env_lightsim.reset()  # do not forget to reset

For the list of availbale solvers, you can consult the "enum" :class:`lightsim2grid.solver.SolverType`.


Detailed usage
--------------------------

.. automodule:: lightsim2grid.solver
    :members:
    :autosummary:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`