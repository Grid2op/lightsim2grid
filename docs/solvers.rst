Available "solvers" (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Type of solvers available
##########################

For now, lightsim2grid ships with at most 3 fully working (and tested) solvers:

- KLUSolver (only under some conditions): it uses the Newton-Raphson algorithm to compute the powerflow, and the
  extremely fast KLU routine to solve the linear systems (currently only available on linux / macos if you compiled
  KLU from source, see the section where the installation of lightsim is described)
- SparseLUSolver: it uses the Newton-Raphson algorithm to compute the powerflow, and the
  slower "SparseLU" solver from the c++ Eigen (available on all platform)
- GaussSeidel: it uses a different algorithm to compute the powerflow. This algorithm is called "Gauss Seidel" and is
  most of the time slower than the Newton Raphson algorithm (available on all platform).


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