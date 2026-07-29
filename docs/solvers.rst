.. currentmodule:: lightsim2grid.algorithm

.. _solvers_doc:

.. _available-powerflow-solvers:

Available powerflow algorithms
===============================

The documentation of this section is in progress. It is rather incomplete for the moment, and only exposes the most
basic features.

If you are interested in collaborating to improve this section, let us know.

.. seealso::

   :ref:`algorithm_names` — explains the three distinct meanings of "solver" in lightsim2grid,
   how the :class:`~lightsim2grid.algorithm.AlgorithmType` enum values are named, and the migration
   table from old names (``KLU``, ``SparseLU``, ``DC``, …) to the new canonical names.

Types of powerflow algorithms
------------------------------

LightSim2Grid supports five families of powerflow algorithms:

- **Gauss-Seidel**: :class:`lightsim2grid.algorithm.GaussSeidelAlgo` and
  :class:`lightsim2grid.algorithm.GaussSeidelSynchAlgo`.
  Solve the AC powerflow using the iterative Gauss-Seidel method (see
  `gausspf <https://matpower.org/docs/ref/matpower5.0/gausspf.html>`_ in MATPOWER).

- **DC approximation**: solve the linearised (DC) power-flow equations using a direct sparse
  factorisation.  Up to four linear-solver backends are available:

  - :class:`lightsim2grid.algorithm.DC_SparseLU` — Eigen SparseLU (always available)
  - :class:`lightsim2grid.algorithm.DC_KLU` — SuiteSparse KLU (when compiled with KLU)
  - :class:`lightsim2grid.algorithm.DC_NICSLU` — NICSLU (requires license + source build)
  - :class:`lightsim2grid.algorithm.DC_CKTSO` — CKTSO (requires license + source build)

- **Newton-Raphson (single slack)**: solves the full AC equations with a single slack bus.
  If multiple slack buses are present only the first is used; the others are treated as PV buses.
  Available with four linear-solver backends:

  - :class:`lightsim2grid.algorithm.NRSing_SparseLU`
  - :class:`lightsim2grid.algorithm.NRSing_KLU`
  - :class:`lightsim2grid.algorithm.NRSing_NICSLU`
  - :class:`lightsim2grid.algorithm.NRSing_CKTSO`

- **Newton-Raphson (distributed / multi-slack)**: solves the full AC equations with multiple
  slack buses.  Available with four linear-solver backends:

  - :class:`lightsim2grid.algorithm.NR_SparseLU`
  - :class:`lightsim2grid.algorithm.NR_KLU`
  - :class:`lightsim2grid.algorithm.NR_NICSLU`
  - :class:`lightsim2grid.algorithm.NR_CKTSO`

- **Fast-Decoupled Powerflow (FDPF)**: the XB and BX variants of the fast-decoupled
  Newton-Raphson method.  Available with four linear-solver backends each:

  - :class:`lightsim2grid.algorithm.FDPF_XB_SparseLU`, :class:`lightsim2grid.algorithm.FDPF_BX_SparseLU`
  - :class:`lightsim2grid.algorithm.FDPF_XB_KLU`, :class:`lightsim2grid.algorithm.FDPF_BX_KLU`
  - :class:`lightsim2grid.algorithm.FDPF_XB_NICSLU`, :class:`lightsim2grid.algorithm.FDPF_BX_NICSLU`
  - :class:`lightsim2grid.algorithm.FDPF_XB_CKTSO`, :class:`lightsim2grid.algorithm.FDPF_BX_CKTSO`

.. warning::
   Algorithms based on ``NICSLU`` and ``CKTSO`` require a compilation from source.
   CKTSO algorithms are (for now) only tested on Linux.

Linear-solver diagnostics: ``LinearSolverStats``
--------------------------------------------------

Every algorithm above is backed by a linear solver (``SparseLU``, ``KLU``, ``NICSLU`` or
``CKTSO``) whose ``analyze``/``factorize``/``refactorize``/``solve`` calls are counted and
timed. Call ``get_linear_solver_stats()`` on a solver (e.g. ``env.backend._grid.get_solver().get_linear_solver_stats()``)
to get a :class:`lightsim2grid.algorithm.LinearSolverStats` with:

- ``nb_analyze`` / ``nb_factorize`` / ``nb_refactorize`` / ``nb_solve`` / ``nb_reset``: how
  many times each was called. These accumulate over the whole lifetime of the solver
  object (not reset every powerflow), so a fallback or failure that fires occasionally is
  distinguishable from one that fires systematically.
- ``nb_refactorize_failed`` / ``nb_fallback_factorize`` / ``nb_fallback_factorize_failed``:
  see :class:`~lightsim2grid.algorithm.NRRefactorRetry_KLU` below.
- ``timer_initialize`` / ``timer_factor`` / ``timer_refactor`` / ``timer_solve``: matching
  durations, reset every ``compute_pf``/``compute_pf_dc`` call like
  :class:`~lightsim2grid.algorithm.TimerJac` (returned by ``get_timers_jacobian()``), which
  these numbers also feed into.

The two-linear-solver Fast-Decoupled family (``FDPF_XB_*``/``FDPF_BX_*``) exposes this per
solver instead, as ``get_linear_solver_stats_bp()`` / ``get_linear_solver_stats_bpp()``
(for B' and B'' respectively) on the concrete solver object.

Retrying a failed refactor: ``NRRefactorRetry_*``
----------------------------------------------------

:class:`lightsim2grid.algorithm.NRRefactorRetry_KLU`,
:class:`~lightsim2grid.algorithm.NRRefactorRetry_CKTSO` and
:class:`~lightsim2grid.algorithm.NRRefactorRetry_NICSLU` are Newton-Raphson (multi-slack)
variants of :class:`~lightsim2grid.algorithm.NR_KLU` / ``NR_CKTSO`` / ``NR_NICSLU``: if a
Jacobian ``refactorize()`` call fails, they fall back to a full ``factorize()`` (reusing the
existing symbolic factorization) before reporting an error, instead of failing immediately.
This is a defensive measure recommended by SuiteSparse's own documentation for KLU,
generalized here to any linear solver with a real factorize/refactorize distinction.

.. note::
   There is no ``NRRefactorRetry_SparseLU``: Eigen's ``SparseLU`` has no cheaper
   "reuse pivot order" refactor -- its ``factorize()`` and ``refactorize()`` are already
   the same call, so the fallback would be a no-op.

.. note::
   These are registered by name only (not part of the :class:`~lightsim2grid.algorithm.AlgorithmType`
   enum), the same way externally-loaded algorithm plugins are -- select them with
   ``grid.change_algorithm("NRRefactorRetry_KLU")`` rather than via ``AlgorithmType``.

Use :class:`~lightsim2grid.algorithm.LinearSolverStats` (``get_linear_solver_stats()``, see
above) to inspect how often the fallback actually fires: ``nb_refactorize_failed`` and
``nb_fallback_factorize`` stay at ``0`` on a grid where refactor never fails.

Default algorithm selection
------------------------------

By default, when KLU is available, lightsim2grid uses:

- :class:`~lightsim2grid.algorithm.NR_KLU` (AC multi-slack)
- :class:`~lightsim2grid.algorithm.NRSing_KLU` (AC single-slack, when only one slack bus is detected)
- :class:`~lightsim2grid.algorithm.DC_KLU` (DC approximation)

When KLU is not available (e.g. installed from PyPI without a source build), it falls back to:

- :class:`~lightsim2grid.algorithm.NR_SparseLU`
- :class:`~lightsim2grid.algorithm.NRSing_SparseLU`
- :class:`~lightsim2grid.algorithm.DC_SparseLU`

Correspondence between class and ``AlgorithmType`` enum
---------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Python class
     - ``AlgorithmType`` enum value
   * - :class:`~lightsim2grid.algorithm.GaussSeidelAlgo`
     - ``AlgorithmType.GaussSeidel``
   * - :class:`~lightsim2grid.algorithm.GaussSeidelSynchAlgo`
     - ``AlgorithmType.GaussSeidelSynch``
   * - :class:`~lightsim2grid.algorithm.DC_SparseLU`
     - ``AlgorithmType.DC_SparseLU``
   * - :class:`~lightsim2grid.algorithm.DC_KLU`
     - ``AlgorithmType.DC_KLU``
   * - :class:`~lightsim2grid.algorithm.DC_NICSLU`
     - ``AlgorithmType.DC_NICSLU``
   * - :class:`~lightsim2grid.algorithm.DC_CKTSO`
     - ``AlgorithmType.DC_CKTSO``
   * - :class:`~lightsim2grid.algorithm.NRSing_SparseLU`
     - ``AlgorithmType.NRSing_SparseLU``
   * - :class:`~lightsim2grid.algorithm.NRSing_KLU`
     - ``AlgorithmType.NRSing_KLU``
   * - :class:`~lightsim2grid.algorithm.NRSing_NICSLU`
     - ``AlgorithmType.NRSing_NICSLU``
   * - :class:`~lightsim2grid.algorithm.NRSing_CKTSO`
     - ``AlgorithmType.NRSing_CKTSO``
   * - :class:`~lightsim2grid.algorithm.NR_SparseLU`
     - ``AlgorithmType.NR_SparseLU``
   * - :class:`~lightsim2grid.algorithm.NR_KLU`
     - ``AlgorithmType.NR_KLU``
   * - :class:`~lightsim2grid.algorithm.NR_NICSLU`
     - ``AlgorithmType.NR_NICSLU``
   * - :class:`~lightsim2grid.algorithm.NR_CKTSO`
     - ``AlgorithmType.NR_CKTSO``
   * - :class:`~lightsim2grid.algorithm.FDPF_XB_SparseLU`
     - ``AlgorithmType.FDPF_XB_SparseLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_BX_SparseLU`
     - ``AlgorithmType.FDPF_BX_SparseLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_XB_KLU`
     - ``AlgorithmType.FDPF_XB_KLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_BX_KLU`
     - ``AlgorithmType.FDPF_BX_KLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_XB_NICSLU`
     - ``AlgorithmType.FDPF_XB_NICSLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_BX_NICSLU`
     - ``AlgorithmType.FDPF_BX_NICSLU``
   * - :class:`~lightsim2grid.algorithm.FDPF_XB_CKTSO`
     - ``AlgorithmType.FDPF_XB_CKTSO``
   * - :class:`~lightsim2grid.algorithm.FDPF_BX_CKTSO`
     - ``AlgorithmType.FDPF_BX_CKTSO``

Usage
------

The preferred way to select an algorithm is to pass ``algo_type`` when creating the backend:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend
    from lightsim2grid.algorithm import AlgorithmType

    env_name = "l2rpn_case14_sandbox"
    env = grid2op.make(env_name,
                       backend=LightSimBackend(algo_type=AlgorithmType.NR_KLU))

You can also change the algorithm after creation, using :func:`lightsim2grid.lightSimBackend.LightSimBackend.set_algo_type`:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend
    from lightsim2grid.algorithm import AlgorithmType

    env_name = "l2rpn_case14_sandbox"
    env = grid2op.make(env_name, backend=LightSimBackend())

    # switch to Gauss-Seidel
    env.backend.set_algo_type(AlgorithmType.GaussSeidel)

    # inspect which algorithms are available in this build
    print(env.backend._grid.available_algorithm_names())

    # tune solver parameters
    env.backend.set_solver_max_iter(10000)
    env.backend.set_tol(1e-7)

.. warning::

   Do not call ``env.backend._grid.change_algorithm(...)`` directly: ``env.reset()``
   rebuilds ``env.backend._grid`` from scratch and re-applies whatever algorithm was
   last set through :func:`~lightsim2grid.lightSimBackend.LightSimBackend.set_algo_type`
   (or the ``algo_type`` kwarg at creation), so a change made directly on ``_grid`` is
   silently reverted on the next reset. Always go through ``set_algo_type`` (or the
   ``algo_type`` kwarg) so the change survives resets.

.. note::

   For the complete list of available algorithm types, see :class:`lightsim2grid.algorithm.AlgorithmType`.
   For an explanation of the naming convention and the three distinct meanings of "solver", see
   :ref:`algorithm_names`.


Detailed API
-------------

.. automodule:: lightsim2grid.algorithm
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
