.. currentmodule:: lightsim2grid.solver

.. _solvers_doc:

Available "solvers" (doc in progress)
=======================================

The documentation of this section is in progress. It is rather incomplete for the moment, and only expose the most
basic features.

If you are interested in collaborating to improve this section, let us know.

Type of solvers available
--------------------------

In lightsim2grid you can have 4 different types of solvers:

- `GaussSeidel` methods: :class:`lightsim2grid.solver.GaussSeidelSolver` and :class:`lightsim2grid.solver.GaussSeidelSynchSolver`
  solves the AC powerflow using the Gauss Seidel method (an example of this algorithm is available in the
  great matpower library here `gausspf <https://matpower.org/docs/ref/matpower5.0/gausspf.html>`_ )
- `DC` methods: solve the DC approximation of the AC powerflow. To solve them it requires manipulating sparse matrices
  and you can use different linear algebra library for that. This is why you have up to 4 different DC solvers:
  :class:`lightsim2grid.solver.DCSolver` (use `Eigen SparseLU <https://eigen.tuxfamily.org/dox/group__SparseLU__Module.html>`_ ),
  :class:`lightsim2grid.solver.KLUDCSolver` (uses `KLU <https://github.com/DrTimothyAldenDavis/SuiteSparse/tree/dev/KLU>`_ )
  :class:`lightsim2grid.solver.NICSLUDCSolver` (uses `NICSLU <https://github.com/chenxm1986/nicslu>`_  and requires and license 
  and to compile lightsim2grid from source)
  :class:`lightsim2grid.solver.CKTSODCSolver` (uses `CKTSO <https://github.com/chenxm1986/cktso>`_  and requires and license 
  and to compile lightsim2grid from source)
- `AC with single slack` methods: solves the AC equations where only one bus is the slack bus (if multiple slack buses are 
  detected, only the first one will be used as slack bus, the others will be treated as "pv" buses). It also exists
  in different "flavours" that uses different linear albrea libraries (same as DC) which are: 
  :class:`lightsim2grid.solver.SparseLUSolverSingleSlack`, :class:`lightsim2grid.solver.KLUSolverSingleSlack`,
  :class:`lightsim2grid.solver.NICSLUSolverSingleSlack` and :class:`lightsim2grid.solver.CKTSOSolverSingleSlack`
- `AC with distributed slack` methods: solves the AC equations with multple slack buses. As for DC and AC with single slack,
  this is avaialble in 4 different flavours (each using internally a different linear albrea solver): 
  :class:`lightsim2grid.solver.SparseLUSolver`, :class:`lightsim2grid.solver.KLUSolver`,
  :class:`lightsim2grid.solver.NICSLUSolver` and :class:`lightsim2grid.solver.CKTSOSolver`

.. warning::
  Solvers based on `NICSLU` and `CKTSO` require a compilation from source. Solvers based on CKTSO are (for now)
  only tested on linux.

By default, when avaialble, lightsim2grid try to use the `KLU` linear solver, so the :class:`lightsim2grid.solver.KLUDCSolver`,
:class:`lightsim2grid.solver.KLUSolverSingleSlack` and :class:`lightsim2grid.solver.KLUSolver`. If not available
(for example if you compiled from source without including the KLU package) it falls back to the "SparseLU" linear solver
so :class:`lightsim2grid.solver.DCSolver`, :class:`lightsim2grid.solver.SparseLUSolverSingleSlack` and 
:class:`lightsim2grid.solver.SparseLUSolver`.

If it detects that the grid is "single slack" it uses the "SingleSlack" version (:class:`lightsim2grid.solver.KLUSolverSingleSlack` or
:class:`lightsim2grid.solver.SparseLUSolverSingleSlack`).

At any moment, you can change the solver used by lightsim2grid with:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend

    # create an environment
    env_name = "l2rpn_case14_sandbox"
    env_lightsim = grid2op.make(env_name, backend=LightSimBackend())

    env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.KLU)  # for KLU solver

Or alternatively, you can change it when you create the backend:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend

    # create an environment
    env_name = "l2rpn_case14_sandbox"
    env_lightsim = grid2op.make(env_name,
                                backend=LightSimBackend(solver_type=lightsim2grid.SolverType.KLU))

The correspondance between the type of solver used (in the above example :class:`lightsim2grid.solver.KLUSolver` )  
and its "name" in the `lightsim2grid.SolverType` (in the above example `lightsim2grid.SolverType.KLU` ) 
module is :

========================================================   =============================================================================
Solver                                                     name in "SolverType"
========================================================   =============================================================================
:class:`lightsim2grid.solver.GaussSeidelSolver`            `GaussSeidel` (SolverType.GaussSeidel)
:class:`lightsim2grid.solver.GaussSeidelSynchSolver`       `GaussSeidelSynch` (SolverType.GaussSeidelSynch)
:class:`lightsim2grid.solver.DCSolver`                     `DC` (SolverType.DC)
:class:`lightsim2grid.solver.KLUDCSolver`                  `KLUDC` (SolverType.KLUDC)
:class:`lightsim2grid.solver.NICSLUDCSolver`               `NICSLUDC` (SolverType.NICSLUDC) 
:class:`lightsim2grid.solver.CKTSODCSolver`                `CKTSODC` (SolverType.CKTSODC)
:class:`lightsim2grid.solver.SparseLUSolverSingleSlack`    `SparseLUSingleSlack` (SolverType.SparseLUSingleSlack)
:class:`lightsim2grid.solver.KLUSolverSingleSlack`         `KLUSingleSlack` (SolverType.KLUSingleSlack)
:class:`lightsim2grid.solver.NICSLUSolverSingleSlack`      `NICSLUSingleSlack` (SolverType.NICSLUSingleSlack)
:class:`lightsim2grid.solver.CKTSOSolverSingleSlack`       `CKTSOSingleSlack` (SolverType.CKTSOSingleSlack)
:class:`lightsim2grid.solver.SparseLUSolver`               `SparseLU` (SolverType.SparseLU)
:class:`lightsim2grid.solver.KLUSolver`                    `KLU` (SolverType.KLU)
:class:`lightsim2grid.solver.NICSLUSolver`                 `NICSLU` (SolverType.NICSLU)
:class:`lightsim2grid.solver.CKTSOSolver`                  `CKTSO` (SolverType.CKTSO)
========================================================   =============================================================================

Usage
--------------------------
In this section we briefly explain how to switch from one solver to another. An example of code using this feature
is given in the
`"benchmark_solvers.py" <https://github.com/BDonnot/lightsim2grid/blob/master/benchmarks/benchmark_solvers.py>`_
script available in the `"benchmarks" <https://github.com/BDonnot/lightsim2grid/tree/master/benchmarks/>`_
directory of the lightsim2grid repository.

To change the solver used by the backend, the preferred solution is to set it once you create it:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend

    # create an environment
    env_name = "l2rpn_case14_sandbox"
    env_lightsim = grid2op.make(env_name,
                                backend=LightSimBackend(solver_type=lightsim2grid.SolverType.KLU)
                               )

.. note::

  For the list of availbale solvers, you can consult the "enum" :class:`lightsim2grid.solver.SolverType`.


You can also (so it's not recommended) change the solver after the backend is created with:

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



Detailed usage
--------------------------

.. automodule:: lightsim2grid.solver
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`