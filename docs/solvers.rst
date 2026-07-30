.. currentmodule:: lightsim2grid.algorithm

.. _solvers_doc:

.. _available-powerflow-solvers:

Available powerflow algorithms
===============================

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


Fine-tuning the Newton-Raphson iteration
------------------------------------------

Every Newton-Raphson-based algorithm above (all the ``NR_*`` / ``NRSing_*`` classes) supports two
independent, orthogonal knobs controlling *how* each iteration is taken:

- :class:`~lightsim2grid.algorithm.ScalingPolicyType` -- whether/how the raw Newton step
  ``(dVa, dVm)`` is scaled down before being applied:

  - ``NoScaling`` (default): the full step is applied (fastest per iteration, but can
    overshoot on a badly-conditioned grid).
  - ``MaxVoltageChange``: clamps the step so it never exceeds ``max_dVa`` (rad) /
    ``max_dVm`` (pu).
  - ``LineSearch``: an Armijo backtracking line search (constants ``ls_c`` / ``ls_rho``).
  - ``Iwamoto``: the Iwamoto optimal multiplier (bounded by ``iw_mu_min`` / ``iw_mu_max``).

- :class:`~lightsim2grid.algorithm.RefactorPolicyType` -- whether the Jacobian is fully
  rebuilt and refactorized every iteration:

  - ``AlwaysRefactor`` (default): rebuild and refactorize ``J`` every iteration.
  - ``EveryN``: refactorize only every ``refactor_every_n`` iterations, updating values
    only (cheaper, at the cost of a slightly stale factorization) in between.
  - ``Chord``: build and factorize ``J`` once, on the first iteration, then reuse that
    factorization for every subsequent one (the "chord method").

The per-policy parameters (``max_dVa``, ``max_dVm``, ``ls_c``, ``ls_rho``, ``ls_max_iter``,
``iw_mu_min``, ``iw_mu_max``, ``refactor_every_n``) are only read by their corresponding
policy; changing them has no effect while a different policy is active.

Setting the policy on a raw solver object (see :ref:`use-solver`) is direct:

.. code-block:: python

    from lightsim2grid.algorithm import NR_KLU, ScalingPolicyType, RefactorPolicyType

    solver = NR_KLU()
    solver.set_scaling_policy(ScalingPolicyType.LineSearch)
    solver.set_refactor_policy(RefactorPolicyType.EveryN)
    solver.set_refactor_every_n(5)

When going through a grid2op / :class:`~lightsim2grid.lightSimBackend.LightSimBackend` powerflow,
use :func:`lightsim2grid.network.LSGrid.get_ac_algo_config` /
:func:`~lightsim2grid.network.LSGrid.set_ac_algo_config` (and their ``_dc_`` counterparts for the
DC solver) instead, which read/write a serialisable
:class:`~lightsim2grid.algorithm.AlgoConfig`:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend
    from lightsim2grid.algorithm import ScalingPolicyType

    env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend())
    grid = env.backend._grid

    config = grid.get_ac_algo_config()
    # config.int_params  == [ScalingPolicyType, RefactorPolicyType, ls_max_iter, refactor_every_n]
    # config.real_params == [max_dVa, max_dVm, ls_c, ls_rho, iw_mu_min, iw_mu_max]

    int_params = list(config.int_params)
    int_params[0] = int(ScalingPolicyType.LineSearch)
    config.int_params = int_params  # reassign the whole list, see warning below
    grid.set_ac_algo_config(config)

.. warning::

   ``AlgoConfig.int_params`` / ``real_params`` are plain lists returned **by value**:
   ``config.int_params[0] = ...`` silently does nothing, because it mutates a temporary
   copy, not the object's actual state. You must build the new list and reassign the
   whole attribute (``config.int_params = int_params``, as above) for the change to take
   effect.

.. warning::

   Just like ``env.backend._grid.change_algorithm(...)`` (see the warning above),
   ``grid.set_ac_algo_config(...)`` / ``set_dc_algo_config(...)`` act directly on
   ``env.backend._grid`` and do **not** survive ``env.reset()``: a reset rebuilds
   ``env.backend._grid`` from scratch with a fresh, default ``AlgoConfig``, silently
   discarding any customisation made this way. Unlike the algorithm type itself, there is
   currently no ``LightSimBackend`` kwarg or setter that re-applies a custom
   ``AlgoConfig`` on every reset -- you need to call ``set_ac_algo_config`` /
   ``set_dc_algo_config`` again after every ``env.reset()`` if you want the change to
   persist.

Detailed API
-------------

.. automodule:: lightsim2grid.algorithm
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
