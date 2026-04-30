.. _algorithm_names:

Naming conventions: "solver" vs "algorithm"
============================================

This page explains how the terms *solver* and *algorithm* are used in
LightSim2Grid, and how to read and build the names of the
:class:`~lightsim2grid.algorithm.AlgorithmType` enum values and the corresponding
Python classes.

Three distinct meanings of "solver"
------------------------------------

The word *solver* appears in three different contexts in LightSim2Grid, and
they should not be confused:

1. **The linear solver** — a numerical library that solves a sparse linear
   system :math:`Ax = b` *eg* at each Newton-Raphson iteration.
   Examples: :class:`~lightsim2grid.algorithm.SparseLULinearSolver` (Eigen built-in),
   KLU (SuiteSparse), NICSLU, CKTSO.
   These names end in ``LinearSolver`` in C++.

2. **The powerflow algorithm** — the outer numerical method used to solve the
   power-flow equations. Examples: Newton-Raphson (NR), DC approximation,
   Fast-Decoupled Power Flow (FDPF), Gauss-Seidel.
   Each algorithm can use one linear solver internally (where applicable).
   This is what :class:`~lightsim2grid.algorithm.AlgorithmType` enumerates, and what
   :class:`~lightsim2grid.algorithm.AlgorithmSelector` dispatches.

3. **Solver bus numbering** — the compact internal bus index used by the sparse
   matrices, which only includes connected buses (not disconnected ones).
   Methods like ``get_Ybus_solver()``, ``get_pv_solver()``,
   ``id_me_to_ac_solver()`` etc. use "solver" in this third sense.
   These are low-level diagnostic accessors; the naming is intentionally kept
   for clarity.

How ``AlgorithmType`` enum names are built
------------------------------------------

Each enum value encodes **both** the outer powerflow algorithm and the linear
solver library used for its internal linear algebra.  The naming pattern is::

    {Algorithm}_{LinearSolver}

Where the algorithm prefix is:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Prefix
     - Meaning
   * - ``NR_``
     - Newton-Raphson with **distributed (multi-) slack**
   * - ``NRSing_``
     - Newton-Raphson with **single slack** (slightly faster when only one slack bus exists)
   * - ``DC_``
     - **DC approximation** (linearised power flow)
   * - ``FDPF_XB_``
     - **Fast-Decoupled Power Flow**, XB variant (``fdxb`` in pypower / pandapower)
   * - ``FDPF_BX_``
     - **Fast-Decoupled Power Flow**, BX variant (``fdbx`` in pypower / pandapower)

And the linear-solver suffix is:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Suffix
     - Library
   * - ``_SparseLU``
     - Eigen's built-in sparse LU — always available, no extra dependency
   * - ``_KLU``
     - SuiteSparse KLU — available when lightsim2grid is compiled with KLU support
   * - ``_NICSLU``
     - NICSLU — requires a license and a source build
   * - ``_CKTSO``
     - CKTSO — requires a license and a source build

**Special cases** (no linear solver choice):

* :class:`~lightsim2grid.algorithm.GaussSeidelAlgo` — the Gauss-Seidel iterative
  method.  No sparse factorisation is involved.
* :class:`~lightsim2grid.algorithm.GaussSeidelSynchAlgo` — synchronous (all-buses-at-once)
  Gauss-Seidel variant.

Complete table
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - ``AlgorithmType`` value
     - Python class
     - Description
   * - ``NR_SparseLU``
     - :class:`~lightsim2grid.algorithm.NR_SparseLU`
     - NR multi-slack + SparseLU (default when KLU unavailable)
   * - ``NRSing_SparseLU``
     - :class:`~lightsim2grid.algorithm.NRSing_SparseLU`
     - NR single-slack + SparseLU
   * - ``DC_SparseLU``
     - :class:`~lightsim2grid.algorithm.DC_SparseLU`
     - DC + SparseLU
   * - ``FDPF_XB_SparseLU``
     - :class:`~lightsim2grid.algorithm.FDPF_XB_SparseLU`
     - Fast-Decoupled XB + SparseLU
   * - ``FDPF_BX_SparseLU``
     - :class:`~lightsim2grid.algorithm.FDPF_BX_SparseLU`
     - Fast-Decoupled BX + SparseLU
   * - ``NR_KLU``
     - :class:`~lightsim2grid.algorithm.NR_KLU`
     - NR multi-slack + KLU (default when KLU available)
   * - ``NRSing_KLU``
     - :class:`~lightsim2grid.algorithm.NRSing_KLU`
     - NR single-slack + KLU
   * - ``DC_KLU``
     - :class:`~lightsim2grid.algorithm.DC_KLU`
     - DC + KLU
   * - ``FDPF_XB_KLU``
     - :class:`~lightsim2grid.algorithm.FDPF_XB_KLU`
     - Fast-Decoupled XB + KLU
   * - ``FDPF_BX_KLU``
     - :class:`~lightsim2grid.algorithm.FDPF_BX_KLU`
     - Fast-Decoupled BX + KLU
   * - ``NR_NICSLU``
     - :class:`~lightsim2grid.algorithm.NR_NICSLU`
     - NR multi-slack + NICSLU (requires license)
   * - ``NRSing_NICSLU``
     - :class:`~lightsim2grid.algorithm.NRSing_NICSLU`
     - NR single-slack + NICSLU (requires license)
   * - ``DC_NICSLU``
     - :class:`~lightsim2grid.algorithm.DC_NICSLU`
     - DC + NICSLU (requires license)
   * - ``FDPF_XB_NICSLU``
     - :class:`~lightsim2grid.algorithm.FDPF_XB_NICSLU`
     - Fast-Decoupled XB + NICSLU
   * - ``FDPF_BX_NICSLU``
     - :class:`~lightsim2grid.algorithm.FDPF_BX_NICSLU`
     - Fast-Decoupled BX + NICSLU
   * - ``NR_CKTSO``
     - :class:`~lightsim2grid.algorithm.NR_CKTSO`
     - NR multi-slack + CKTSO (requires license)
   * - ``NRSing_CKTSO``
     - :class:`~lightsim2grid.algorithm.NRSing_CKTSO`
     - NR single-slack + CKTSO (requires license)
   * - ``DC_CKTSO``
     - :class:`~lightsim2grid.algorithm.DC_CKTSO`
     - DC + CKTSO (requires license)
   * - ``FDPF_XB_CKTSO``
     - :class:`~lightsim2grid.algorithm.FDPF_XB_CKTSO`
     - Fast-Decoupled XB + CKTSO
   * - ``FDPF_BX_CKTSO``
     - :class:`~lightsim2grid.algorithm.FDPF_BX_CKTSO`
     - Fast-Decoupled BX + CKTSO
   * - ``GaussSeidel``
     - :class:`~lightsim2grid.algorithm.GaussSeidelAlgo`
     - Gauss-Seidel iterative (no sparse factorisation)
   * - ``GaussSeidelSynch``
     - :class:`~lightsim2grid.algorithm.GaussSeidelSynchAlgo`
     - Synchronous Gauss-Seidel iterative

Usage example
-------------

.. code-block:: python

    import lightsim2grid
    from lightsim2grid.algorithm import AlgorithmType

    env = grid2op.make("l2rpn_case14_sandbox", backend=lightsim2grid.LightSimBackend())

    # switch to NR with KLU (if available)
    env.backend.set_algo_type(AlgorithmType.NR_KLU)

    # inspect what algorithm is currently active
    print(env.backend.get_algo_types())
    # >>> (<AlgorithmType.NRSing_KLU: 7>, <AlgorithmType.DC_KLU: 9>)

    # list all algorithms available in this build
    env.backend._grid.available_algorithm_names()

.. note::
   The method is called ``change_algorithm(type)`` in the C++ LSGrid
   Python binding (``env.backend._grid.change_algorithm(AlgorithmType.NR_KLU)``),
   consistent with pandapower's ``algorithm=`` parameter and MATPOWER's
   ``pf.alg`` option.

Migration from old names
------------------------

Before version 0.14, the enum values and class names used shorter names that
did not make the algorithm component explicit (``KLU``, ``SparseLU``, ``DC``,
etc.).  These old names are kept as **deprecated aliases** in the Python enum
so that existing code continues to work, but they will be removed in a future
release:

.. list-table::
   :header-rows: 1
   :widths: 40 40

   * - Old name (deprecated)
     - New canonical name
   * - ``AlgorithmType.SparseLU``
     - ``AlgorithmType.NR_SparseLU``
   * - ``AlgorithmType.SparseLUSingleSlack``
     - ``AlgorithmType.NRSing_SparseLU``
   * - ``AlgorithmType.DC``
     - ``AlgorithmType.DC_SparseLU``
   * - ``AlgorithmType.KLU``
     - ``AlgorithmType.NR_KLU``
   * - ``AlgorithmType.KLUSingleSlack``
     - ``AlgorithmType.NRSing_KLU``
   * - ``AlgorithmType.KLUDC``
     - ``AlgorithmType.DC_KLU``
   * - ``AlgorithmType.NICSLU``
     - ``AlgorithmType.NR_NICSLU``
   * - ``AlgorithmType.NICSLUSingleSlack``
     - ``AlgorithmType.NRSing_NICSLU``
   * - ``AlgorithmType.NICSLUDC``
     - ``AlgorithmType.DC_NICSLU``
   * - ``AlgorithmType.CKTSO``
     - ``AlgorithmType.NR_CKTSO``
   * - ``AlgorithmType.CKTSOSingleSlack``
     - ``AlgorithmType.NRSing_CKTSO``
   * - ``AlgorithmType.CKTSODC``
     - ``AlgorithmType.DC_CKTSO``

.. seealso::

   :ref:`available-powerflow-solvers` for performance comparisons and
   practical guidance on choosing an algorithm.
