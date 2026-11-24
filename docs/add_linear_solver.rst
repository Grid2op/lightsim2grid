.. _add_linear_solver:

Adding a linear solver
======================

The power-flow algorithms of lightsim2grid that solve a sparse linear system :math:`Ax = b`
(Newton-Raphson, DC, fast-decoupled, and the batch classes built on them) get that linear algebra
from a **linear solver**. Not all of them do: Gauss-Seidel, for example, updates one bus at a time
and uses none.
The ones shipped today are Eigen's ``SparseLU`` (always available), KLU (bundled from SuiteSparse),
and NICSLU and CKTSO (optional, found at build time). They all live in ``src/core/linear_solvers/``.

This page explains how to add another one, for example a vendor library, a GPU solver, or a solver
you wrote yourself, whether it exposes a C API or a C++ one:

- `The linear solver interface`_ is what you write in either case: one C++ class, what it receives,
  what it must do, and when lightsim2grid calls it;
- `Integrating it`_ gives the two ways of making lightsim2grid use that class: **a plugin**, loaded at
  run time by an installed lightsim2grid, or **a fork** of lightsim2grid compiled from source;
- `Testing it`_ applies to both.

The terms used below follow :ref:`algorithm_names`: a *linear solver* is the ``XxxLinearSolver``
class, an *algorithm* is what uses it (``NR_Xxx``, ``DC_Xxx``, ...), and an algorithm is selected
by its *name*.

The running example wraps a C library, ``mylu``, whose header is:

.. code-block:: c

    /* mylu.h */
    typedef struct mylu_handle mylu_handle;
    int  mylu_create(mylu_handle ** h);
    void mylu_destroy(mylu_handle * h);
    int  mylu_analyze(mylu_handle * h, int n, const int * colptr, const int * rowind);
    int  mylu_factor(mylu_handle * h, const double * values);    /* chooses the pivots */
    int  mylu_refactor(mylu_handle * h, const double * values);  /* reuses them */
    int  mylu_solve(mylu_handle * h, double * b);                /* in place */

All the functions return ``0`` on success. Most sparse direct solvers (KLU, UMFPACK, PARDISO,
NICSLU, CKTSO, cuDSS, ...) have this shape: a handle, a symbolic analysis of the sparsity pattern, a
numeric factorization, a cheaper re-factorization that reuses the analysis, and a solve. A C++
library is wrapped the same way, calling its methods instead of free functions.
``examples/linear_solver_plugin/`` holds a stand-in implementation of this header, so that the
example builds and runs on its own.

The linear solver interface
---------------------------

A linear solver is **not** a virtual class: it is a template parameter of the algorithms
(``NRAlgo<LinearSolverPolicy<XxxLinearSolver>, ...>``), so its type is known at compile time and no
call goes through a vtable. The algorithms reach it through ``LinearSolverPolicy`` (see
``linear_solvers/LinearSolverPolicy.hpp``), a thin wrapper that counts and times every call and adds
the refactorization fallback described below.

What to implement
~~~~~~~~~~~~~~~~~

.. code-block:: cpp

    class XxxLinearSolver final {
    public:
        XxxLinearSolver() noexcept;               // default-constructible

        ErrorType reset();
        ErrorType analyze(const EigenRefConstRealSpMat & J);
        ErrorType factorize(const EigenRefConstRealSpMat & J);
        ErrorType refactorize(const EigenRefConstRealSpMat & J);
        ErrorType solve(Eigen::Ref<RealVect> b);
        ErrorType solve_transpose(Eigen::Ref<RealVect> b);

        static constexpr bool CAN_SOLVE_MAT = false;
        static constexpr bool CAN_SOLVE_TRANSPOSE = false;
    };

``EigenRefConstRealSpMat``, ``RealVect`` and ``ErrorType`` come from ``Utils.hpp``.

``reset()``
    Free everything the previous analysis built, so that the next call can be an ``analyze``. It is
    also called before any ``analyze`` was done, and it can be called twice in a row: it must not
    fail in either case.

``analyze(J)``
    The symbolic step: ordering, symbolic factorization, allocation. Only the pattern of ``J`` is
    meaningful at this point; its values are those of the first iteration, but do not rely on them.

``factorize(J)``
    A numeric factorization that **chooses its pivots**.

``refactorize(J)``
    A numeric factorization with the same pattern, that may **reuse the pivot order** of the last
    ``factorize``. This is the call that makes a solver fast in lightsim2grid. A library that has no
    such call (Eigen's own ``SparseLU``, for example) can simply call its ``factorize`` here: the
    results are the same, the computation is only slower. When the reused pivots turn out to be
    bad (a zero or a tiny pivot), return ``ErrorType::SolverReFactor`` rather than a poor
    factorization: see `When it is called`_.

``solve(b)``
    Solve :math:`J x = b` with the last factorization, **in place**: on return ``b`` holds ``x``. If
    the library cannot work in place, solve into a temporary and copy it back.

``solve_transpose(b)`` (optional)
    The method must exist, but it does not have to solve anything. When implemented, it solves :math:`J^T x = b` with the factorization of ``J``, in place. Only the gradient-based batch
    tools (``BatchAdjoint``) use it, and they do not pick up a new linear solver anyway. If you do not
    implement it, return ``ErrorType::NotImplemented`` and set ``CAN_SOLVE_TRANSPOSE = false``:
    never solve :math:`J x = b` instead.

``CAN_SOLVE_MAT``
    Whether ``solve`` accepts a matrix right hand side. ``false`` is always correct.

Return ``ErrorType::NoError`` on success. On failure, return the code of the step that failed
(``SolverAnalyze``, ``SolverFactor``, ``SolverReFactor``, ``SolverSolve``); ``NotInitError`` for a
call made before ``analyze``; ``LicenseError`` if a licensed library refuses to start. The algorithm
stops at the first error and reports it (``get_error()``).

.. danger::
    **Do not throw** from these methods. A failed factorization is an expected outcome (a diverging
    Newton, a split grid), and the algorithms handle it through the return code: report it as an
    ``ErrorType``, not as an exception.

What you receive
~~~~~~~~~~~~~~~~

``J`` is an ``Eigen::SparseMatrix<double>`` seen through an ``Eigen::Ref``: the Jacobian of the
Newton-Raphson, or the :math:`B` matrix of the DC and fast-decoupled algorithms.

- It is **square** and stored in **compressed sparse column** (CSC) format, with **0-based int**
  indices. The three arrays are:

  .. code-block:: cpp

      const int n = static_cast<int>(J.cols());   // == J.rows()
      const int nnz = static_cast<int>(J.nonZeros());
      const int * colptr = J.outerIndexPtr();     // size n + 1, colptr[n] == nnz
      const int * rowind = J.innerIndexPtr();     // size nnz, sorted within a column
      const double * values = J.valuePtr();       // size nnz, in the order of rowind

  so the entries of column ``j`` are ``values[k]``, at row ``rowind[k]``, for ``k`` from
  ``colptr[j]`` to ``colptr[j + 1] - 1``.
- It may hold **explicit zeros**. The batch algorithms keep one sparsity pattern for a whole batch
  and change a bus's role by rewriting values in it: a value can be ``0`` in one call and non-zero in
  the next. A library that drops zero entries during its analysis cannot be used as is: give it the
  pattern, not the values, or turn that pruning off.
- Between an ``analyze`` and the next ``reset`` / ``analyze``, the **pattern does not change**: only
  the values do. The memory address of the arrays is not guaranteed to stay the same, so read the
  pointers again at each call rather than keeping them.
- It is not symmetric in general (a Newton-Raphson Jacobian is not), and it is not guaranteed to be
  diagonally dominant.

``b`` is a dense vector of ``double`` of size ``n`` (``b.data()`` is a contiguous array).

If your library wants something else, convert in the wrapper:

- **64-bit indices** (``int64_t``, ``long long``): copy the two index arrays once in ``analyze``
  into buffers the class owns (``CKTSOSolver.cpp`` does this);
- **compressed sparse row** (CSR): the CSC arrays of ``J`` are the CSR arrays of :math:`J^T`. Either
  ask the library to solve with the transpose of what it factorized (most have that option), or
  transpose the matrix in ``analyze`` and keep the permutation that maps the values;
- **1-based indices** (Fortran-style libraries): shift the indices once, in ``analyze``;
- **single precision, GPU memory**: copy the values in ``factorize`` / ``refactorize`` and the right
  hand side in ``solve``. The algorithms always work in ``double`` on the host (``real_type`` in
  ``Utils.hpp``).

When it is called
~~~~~~~~~~~~~~~~~

``analyze`` is called **when the sparsity pattern changes**, not once per power flow: at the first
solve, and after a change of topology or of the pv / pq split. Each time, the algorithm first calls
``reset`` and then ``analyze`` + ``factorize``. As long as the pattern stays the same, the next
power flows (the next grid2op steps, the next rows of a batch) reuse that analysis.

The **Newton-Raphson** algorithms (``NRAlgo.tpp``) call, for one power flow::

    [ reset(); analyze(J); ]   # only if the pattern changed since the last power flow
    factorize(J)               # first iteration: right after an analyze, or when asked to
      or refactorize(J)        #   otherwise
    solve(b)
    refactorize(J); solve(b)   # every next iteration: same pattern, new values
    refactorize(J); solve(b)
    ...

``factorize`` without an ``analyze`` happens when the algorithm no longer trusts the pivot order of
the last factorization (``request_factorize``, used by the outer loops after they rewrite rows). The
refactor policy of the Newton (``set_refactor_policy``, ``AlwaysRefactor`` by default) can also skip
the refactorization of some iterations and solve with the previous factors.

**When a refactorization fails**, what happens depends on the algorithm:

- by default (``NR_*``, ``NRSing_*``, and ``DC_*``), the error is reported and the power flow fails;
- with the **fallback** on, ``LinearSolverPolicy`` calls ``factorize`` on the same matrix, which
  chooses new pivots, and the solve goes on; the next iterations refactorize again, on those new
  pivots. If that ``factorize`` fails too, the error is reported. The fallback is on for
  ``NRRefactorRetry_*`` and ``NROuter_*``, and the batch classes turn it on when they change a
  bus's role by value: masking the buses a contingency cut off, switching a bus between PV and PQ,
  and the like (``set_refactor_fallback``).

So a ``refactorize`` that returns ``SolverReFactor`` on bad pivots costs one extra ``factorize``
where the fallback is on, and is a clear error everywhere else. A ``refactorize`` that returns a
poor factorization instead gives a Newton that diverges, or converges to the wrong point.

The **DC** algorithms (``BaseDCAlgo.tpp``) call ``analyze`` + ``factorize`` when the pattern
changes, ``refactorize`` when only the values of the matrix changed (a contingency of a batch), and
one ``solve`` per power flow. The **fast-decoupled** algorithms hold two linear solvers, for
:math:`B'` and :math:`B''`. Both matrices stay constant during a solve, so each is analyzed and
factorized once, then only solved, twice per iteration.

The batch classes (``ContingencyAnalysis``, ``TimeSeries``, ``InjectionSweep``, ``ScenarioSweep``)
keep one pattern for every row of the batch: a contingency analysis of thousands of rows calls
``analyze`` **once**, then ``refactorize`` and ``solve`` for every row.

Good practice
~~~~~~~~~~~~~

None of the following is checked by lightsim2grid, and none of it is required to get a working
solver. The linear solver is your code, outside lightsim2grid: a leak in it, or a licence that
restricts how the library may be used, is the solver's concern, not lightsim2grid's. It can still
make some lightsim2grid features unusable with that solver, which is why these are recommended.

- **Copy and move.** Linear solvers are never copied nor moved: the algorithms hold theirs by
  value. A compile-time check is always preferable to a run-time surprise, so it is good practice
  to delete the four copy and move operations, as the existing solvers do.
- **Memory.** Be careful with memory allocation: we strongly encourage RAII (``std::unique_ptr``
  with a custom deleter for a C handle, as in the example below), so that ``reset`` and the
  destructor cannot leak, and a destructor never throws.
- **Several instances.** The batch classes run several algorithms, sometimes in parallel, each with
  its own linear solver. A library with global state (a licence check, a thread pool, a static work
  buffer) that does not tolerate several instances alive at the same time may not be usable with
  these more advanced features of lightsim2grid.
- **C++ standard.** For a later integration in lightsim2grid itself, the core is written in C++14
  at the time of writing (and compiled up to C++26 in CI), so no ``if constexpr``,
  ``std::optional`` or structured bindings there. A plugin, or code of your own in a fork, can use
  any standard, as long as it is compiled with the same toolchain as lightsim2grid.

Optional capabilities
~~~~~~~~~~~~~~~~~~~~~

``LinearSolverPolicy`` and the algorithms detect a few more members at compile time and use them
only when the class has them. You can ignore them all at first.

``set_rgrowth_threshold(real_type)`` and ``last_rgrowth()``
    For a solver that can measure the reciprocal pivot growth of a refactorization: with the
    fallback on, a refactorization below the threshold is rejected and redone by ``factorize``.
    See ``KLULinearSolver``.

``set_pivot_growth_guard(bool)``
    Lets an algorithm ask for the check above on the Jacobians where it matters.

``set_ordering_constraints(OrderingConstraints)``, ``set_fill_reducing_ordering(FillReducingOrdering)`` and ``ORDERS_WITH_CONSTRAINTS_BY_DEFAULT``
    For a solver that can order with constraints (some indices eliminated after others), which keeps
    a usable pivot in the rows the outer loops switch. See ``OrderingConstraints.hpp`` and
    ``KLUStaticPivotLinearSolver``.

If your library decides on its own between a factorization and a refactorization (NICSLU does), the
counters of ``LinearSolverPolicy`` would report the wrong thing: specialize its ``factorize`` and
``refactorize`` for your class, as ``LinearSolverPolicy.hpp`` does for ``NICSLULinearSolver``.

Example: wrapping ``mylu``
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

    #include <memory>

    #include "Utils.hpp"
    #include "Eigen/Core"
    #include "Eigen/SparseCore"

    // the library's C header; wrap it in extern "C" { } if it has no guard of its own
    #include "mylu.h"

    namespace ls2g {

    class MyLULinearSolver final
    {
        public:
            MyLULinearSolver() noexcept = default;
            ~MyLULinearSolver() noexcept = default;  // handle_ frees itself

            ErrorType reset(){
                handle_.reset();  // the next call is an analyze()
                return ErrorType::NoError;
            }

            ErrorType analyze(const EigenRefConstRealSpMat & J){
                mylu_handle * raw = nullptr;
                if(mylu_create(&raw) != 0) return ErrorType::SolverAnalyze;
                handle_.reset(raw);
                const int n = static_cast<int>(J.cols());
                // J is compressed column storage, 0-based, int indices: handed over as is
                if(mylu_analyze(handle_.get(), n, J.outerIndexPtr(), J.innerIndexPtr()) != 0){
                    handle_.reset();
                    return ErrorType::SolverAnalyze;
                }
                return ErrorType::NoError;
            }

            ErrorType factorize(const EigenRefConstRealSpMat & J){
                if(!handle_) return ErrorType::NotInitError;
                if(mylu_factor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverFactor;
                return ErrorType::NoError;
            }

            ErrorType refactorize(const EigenRefConstRealSpMat & J){
                if(!handle_) return ErrorType::NotInitError;
                // a library with no refactorization of its own can simply
                // `return factorize(J);` here
                // SolverReFactor: what the fallback answers to, where it is on
                if(mylu_refactor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverReFactor;
                return ErrorType::NoError;
            }

            ErrorType solve(Eigen::Ref<RealVect> b){
                if(!handle_) return ErrorType::NotInitError;
                if(mylu_solve(handle_.get(), b.data()) != 0) return ErrorType::SolverSolve;
                if(!b.allFinite()) return ErrorType::InifiniteValue;
                return ErrorType::NoError;
            }

            ErrorType solve_transpose(Eigen::Ref<RealVect> /*b*/){
                return ErrorType::NotImplemented;  // see CAN_SOLVE_TRANSPOSE
            }

            static constexpr bool CAN_SOLVE_MAT = false;
            static constexpr bool CAN_SOLVE_TRANSPOSE = false;

            MyLULinearSolver(const MyLULinearSolver &) = delete;
            MyLULinearSolver(MyLULinearSolver &&) = delete;
            MyLULinearSolver & operator=(const MyLULinearSolver &) = delete;
            MyLULinearSolver & operator=(MyLULinearSolver &&) = delete;

        private:
            // RAII (good practice): the unique_ptr frees the handle, in reset() and in the
            // destructor, whatever happens in between
            struct HandleDeleter {
                void operator()(mylu_handle * h) const noexcept { mylu_destroy(h); }
            };
            std::unique_ptr<mylu_handle, HandleDeleter> handle_;
    };

    }  // namespace ls2g

This is ``examples/linear_solver_plugin/MyLULinearSolver.hpp``, header-only. The fork below splits
the same class in a ``.hpp`` and a ``.cpp``, and adds two things to it.

Integrating it
--------------

A linear solver on its own is not selectable: what is selected is an *algorithm* using it, by name.
There are two ways to build and register those algorithms.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - a) plugin
     - b) fork
   * - lightsim2grid
     - installed (``pip install lightsim2grid``)
     - cloned, modified and compiled from source
   * - your code lives
     - in your own repository, as a shared library
     - in ``src/core/`` of your copy of lightsim2grid
   * - available algorithms
     - Newton-Raphson (``NR_``, ``NRSing_``, ``NROuter_``, ...) and DC
     - the same, plus fast-decoupled (``FDPF_``)
   * - a Python class per algorithm, a ``compilation_options`` flag
     - no
     - can be added
   * - must be loaded before use
     - yes, ``load_algorithm_plugin`` in every process
     - no

Start with the plugin unless you need what only the fork has.

a) As a plugin
~~~~~~~~~~~~~~

lightsim2grid can load algorithms from a shared library at run time (:ref:`solver_plugin`
describes that mechanism in full). Because the linear solver is a template parameter, a plugin can
carry its own: it defines ``MyLULinearSolver`` and registers ``NRAlgo<...>`` and the others
instantiated on it, from the headers an installed lightsim2grid ships. lightsim2grid itself is not
modified. ``examples/linear_solver_plugin/`` is a complete one:

.. code-block:: text

    examples/linear_solver_plugin/
    ├── mylu/mylu.h, mylu/mylu.c   the stand-in C library (replace by yours)
    ├── MyLULinearSolver.hpp       the wrapper above
    ├── MyLU_plugin.cpp            the algorithms and their registration
    ├── CMakeLists.txt             builds the plugin against an installed lightsim2grid
    └── test_plugin.py             loads it and compares it with SparseLU

**The registration** (``MyLU_plugin.cpp``):

.. code-block:: cpp

    #include <memory>

    #include <AlgorithmRegistry.hpp>
    // the algorithm templates, without Solvers.hpp: that one also includes the headers of
    // KLU, NICSLU and CKTSO, which an installed lightsim2grid does not ship
    #include <linear_solvers/LinearSolverPolicy.hpp>
    #include <powerflow_algorithm/BaseDCAlgo.hpp>
    #include <powerflow_algorithm/NRAlgo.hpp>
    #include <powerflow_algorithm/NROuterInner.hpp>
    #include <powerflow_algorithm/OuterLoopAlgo.hpp>

    #include "MyLULinearSolver.hpp"

    namespace {
        using MyLU = ls2g::LinearSolverPolicy<ls2g::MyLULinearSolver>;
        using NR_MyLU = ls2g::NRAlgo<MyLU, ls2g::MultiSlackNRSystem>;
        using NRSing_MyLU = ls2g::NRAlgo<MyLU, ls2g::SingleSlackNRSystem>;
        using NROuter_MyLU = ls2g::OuterLoopAlgo<ls2g::NROuterInner<MyLU>>;
        using DC_MyLU = ls2g::BaseDCAlgo<MyLU>;
    }

    static void register_plugin(ls2g::PluginRegistrar& reg) {
        reg.add("NR_MyLU",      []{ return std::unique_ptr<ls2g::BaseAlgo>(new NR_MyLU()); });
        reg.add("NRSing_MyLU",  []{ return std::unique_ptr<ls2g::BaseAlgo>(new NRSing_MyLU()); });
        reg.add("NROuter_MyLU", []{ return std::unique_ptr<ls2g::BaseAlgo>(new NROuter_MyLU()); });
        reg.add("DC_MyLU",      []{ return std::unique_ptr<ls2g::BaseAlgo>(new DC_MyLU()); });
    }
    LS2G_PLUGIN_ENTRY(register_plugin)

Any other combination the fork lists in `Declaring the algorithms`_ can be registered the same way,
except the fast-decoupled ones (see the limits below). The wrapper class carries no ``LS2G_API``:
that macro is for the symbols of ``lightsim2grid_core`` itself, and means *imported from it* in a
plugin on Windows.

**The build** (``CMakeLists.txt``) follows the other example plugins: ``find_package`` of the
installed ``lightsim2grid_core``, a ``MODULE`` library linked to ``lightsim2grid::core`` and to the
linear solver library, and ``ls2g_match_core_build_flags``. The example builds its stand-in library
with it; a real one is found instead:

.. code-block:: cmake

    find_library(MYLU_LIB NAMES mylu HINTS "${MYLU_ROOT}/lib" REQUIRED)
    find_path(MYLU_INCLUDE mylu.h HINTS "${MYLU_ROOT}/include" REQUIRED)
    target_link_libraries(linear_solver_plugin PRIVATE lightsim2grid::core ${MYLU_LIB})
    target_include_directories(linear_solver_plugin PRIVATE ${MYLU_INCLUDE})

A static library, compiled as position-independent code, ends up inside the plugin. A shared one
must be found by the dynamic loader when the plugin is loaded (``LD_LIBRARY_PATH``, the rpath, or
next to the plugin).

**Building and loading it**:

.. code-block:: bash

    LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
    cmake -S examples/linear_solver_plugin -B examples/linear_solver_plugin/build \
          -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE" -DCMAKE_BUILD_TYPE=Release
    cmake --build examples/linear_solver_plugin/build --config Release
    python examples/linear_solver_plugin/test_plugin.py

(``python examples/build_plugins.py linear_solver_plugin --test`` does the same.) Then, in every
process that uses it:

.. code-block:: python

    import numpy as np
    import pandapower.networks as pn
    import lightsim2grid
    from lightsim2grid.network import init_from_pandapower
    from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

    lightsim2grid.load_algorithm_plugin("examples/linear_solver_plugin/build/liblinear_solver_plugin.so")

    grid = init_from_pandapower(pn.case118())
    grid.change_algorithm("NR_MyLU")
    V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 10, 1e-8)
    print(grid.get_algo().get_linear_solver_stats())   # analyze / factorize / refactorize / solve counts

    ca = ContingencyAnalysisCPP(grid)                    # the batch classes take the name too
    ca.change_algorithm("NR_MyLU")

    from lightsim2grid import LightSimBackend            # and so does grid2op
    backend = LightSimBackend(algo_type="NR_MyLU")

**Limits of the plugin route:**

- **No fast-decoupled algorithm.** Part of their code (``BaseFDPFAlgo::fillBp_Bpp``) is only
  compiled inside ``lightsim2grid_core``, in ``SolverInstantiations.cpp``, for the linear solvers it
  knows. A plugin registering one still builds, but loading it fails on that undefined symbol.
- **The same build flags as the core.** The plugin and ``lightsim2grid_core`` exchange Eigen
  objects, so they must agree on ``-march=native``: ``ls2g_match_core_build_flags`` takes care of
  it (see :ref:`solver_plugin`, *Matching build flags*).
- **Rebuild it for each lightsim2grid version.** The plugin compiles the algorithms from the
  headers of the version it was built against. Loading only checks that Eigen is configured the
  same way on both sides, not that those headers match the installed version.
- **Load it before use, everywhere.** A grid saved (pickled, :ref:`binary-serialization`) while
  using ``NR_MyLU`` only loads where the plugin is loaded; so does a ``LightSimBackend`` copied to a
  worker process.
- **Selected by name only.** There is no ``lightsim2grid.algorithm.NR_MyLU`` Python class and no
  ``AlgorithmType`` value: ``get_algo_type()`` reads ``AlgorithmType.Custom``. The algorithm in use
  is still known by its name, through ``get_algo_name()`` (``"NR_MyLU"``).

b) As part of lightsim2grid (fork)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the solver becomes a built-in, like KLU or CKTSO: it is compiled into ``lightsim2grid_core``
when its library is found at build time. That takes six files:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - File
     - What to add
   * - ``src/core/linear_solvers/MyLUSolver.{hpp,cpp}``
     - the wrapper class
   * - ``src/core/CMakeLists.txt``
     - find the library, link it, define ``MYLU_SOLVER_AVAILABLE``
   * - ``src/core/linear_solvers/CMakeLists.txt``
     - compile ``MyLUSolver.cpp`` when the library was found
   * - ``src/core/Solvers.hpp``
     - the algorithm aliases (``NR_MyLU``, ...) and their ``extern template`` declarations
   * - ``src/core/SolverInstantiations.cpp``
     - the explicit instantiations of those aliases
   * - ``src/core/BuiltinSolversRegistration.cpp``
     - register each alias under a name

**The wrapper.** The class of the example, with two changes:

- it is marked ``LS2G_API`` (from ``Utils.hpp``), so that it is exported from the shared library on
  Windows: ``class LS2G_API MyLULinearSolver final``;
- the whole header is guarded by ``#ifdef MYLU_SOLVER_AVAILABLE`` (and the ``.cpp`` with it), which
  the build defines only when the library was found, so that lightsim2grid still compiles without
  it.

To check that it compiles, without a full build (from the repository root, with the submodules
fetched):

.. code-block:: bash

    g++ -fsyntax-only -std=c++14 -DLS2G_BUILDING_CORE -DMYLU_SOLVER_AVAILABLE \
        -I src -I src/core -I eigen -I /path/to/mylu/include \
        src/core/linear_solvers/MyLUSolver.cpp

**Building it.** The existing optional solvers are found through an environment variable pointing
at the library's root folder (``PATH_NICSLU``, ``PATH_CKTSO``): do the same. In
``src/core/CMakeLists.txt``, find the library. This block must come **before**
``add_subdirectory(linear_solvers)``, which reads ``LIB_MYLU``:

.. code-block:: cmake

    # =============================================================================
    # MyLU
    # =============================================================================
    set(PATH_MYLU "$ENV{PATH_MYLU}" CACHE PATH "Path to the MyLU root directory")
    if(PATH_MYLU)
        unset(LIB_MYLU CACHE)
        find_library(LIB_MYLU NAMES mylu HINTS ${PATH_MYLU}/lib NO_DEFAULT_PATH)
        if(LIB_MYLU)
            message(STATUS "MyLU enabled: ${LIB_MYLU}")
        else()
            message(WARNING "MyLU disabled: no library found in ${PATH_MYLU}/lib")
        endif()
    endif()

and, further down, next to the NICSLU and CKTSO blocks (once ``lightsim2grid_core`` exists), link it:

.. code-block:: cmake

    if(LIB_MYLU)
        target_link_libraries(lightsim2grid_core PRIVATE ${LIB_MYLU})
        # PUBLIC, BUILD_INTERFACE only: the Python bindings include MyLUSolver.hpp too
        target_include_directories(lightsim2grid_core PUBLIC $<BUILD_INTERFACE:${PATH_MYLU}/include>)
        target_compile_definitions(lightsim2grid_core PUBLIC MYLU_SOLVER_AVAILABLE)
    endif()

The compile definition is ``PUBLIC`` because the headers that use it (``Solvers.hpp``) are compiled
by the Python bindings and the C++ tests too. In ``src/core/linear_solvers/CMakeLists.txt``,
compile the wrapper only when the library is there:

.. code-block:: cmake

    if(LIB_MYLU)
        list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/MyLUSolver.cpp")
    endif()

If the library comes as a CMake package, ``find_package(MyLU)`` and linking its imported target
works as well. A shared library outside the system paths must be found by the dynamic loader at
run time (``LD_LIBRARY_PATH``, the rpath, or next to ``lightsim2grid_core``).

Declaring the algorithms
++++++++++++++++++++++++

In ``src/core/Solvers.hpp``, include the header with the others and declare the aliases:

.. code-block:: cpp

    #include "linear_solvers/MyLUSolver.hpp"

    // ...

    #ifdef MYLU_SOLVER_AVAILABLE
        /** Newton-Raphson (multi-slack) with the MyLU linear solver **/
        using NR_MyLU = NRAlgo<LinearSolverPolicy<MyLULinearSolver>, MultiSlackNRSystem>;
        /** Newton-Raphson (single-slack) with the MyLU linear solver **/
        using NRSing_MyLU = NRAlgo<LinearSolverPolicy<MyLULinearSolver>, SingleSlackNRSystem>;
        /** DC approximation with the MyLU linear solver **/
        using DC_MyLU = BaseDCAlgo<LinearSolverPolicy<MyLULinearSolver>>;
        /** Fast-Decoupled Power Flow (XB variant) with the MyLU linear solver **/
        using FDPF_XB_MyLU = BaseFDPFAlgo<LinearSolverPolicy<MyLULinearSolver>, FDPFMethod::XB>;
        /** Single-slack Newton-Raphson with outer loops, MyLU linear solver **/
        using NROuter_MyLU = OuterLoopAlgo<NROuterInner<LinearSolverPolicy<MyLULinearSolver>>>;
    #endif

The full menu, with the solvers that have them, is: ``NR_`` and ``NRSing_`` (Newton-Raphson with a
distributed or a single slack), ``DC_``, ``FDPF_XB_`` and ``FDPF_BX_`` (fast-decoupled),
``NROuter_`` and ``NROuterDist_`` (Newton-Raphson with outer loops, see :ref:`outer_loops`), and
``NRRefactorRetry_`` (``NRAlgo<RefactorRetryLinearSolver<...>, MultiSlackNRSystem>``, only useful
for a solver whose ``refactorize`` differs from its ``factorize``). :ref:`available-powerflow-solvers`
describes them.

In the same file, inside the ``#ifndef LS2G_BUILDING_CORE`` block, declare each of them
``extern template``. This tells every other translation unit (the bindings, the tests, a plugin) that
the code is already compiled in ``lightsim2grid_core``, instead of compiling it again:

.. code-block:: cpp

    #ifdef MYLU_SOLVER_AVAILABLE
        extern template class LS2G_API NRAlgo<LinearSolverPolicy<MyLULinearSolver>, MultiSlackNRSystem>;
        extern template class LS2G_API NRAlgo<LinearSolverPolicy<MyLULinearSolver>, SingleSlackNRSystem>;
        extern template class LS2G_API BaseDCAlgo<LinearSolverPolicy<MyLULinearSolver>>;
        extern template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<MyLULinearSolver>, FDPFMethod::XB>;
        extern template class LS2G_API OuterLoopAlgo<NROuterInner<LinearSolverPolicy<MyLULinearSolver>>>;
    #endif

and in ``src/core/SolverInstantiations.cpp``, instantiate them, with the same list:

.. code-block:: cpp

    #ifdef MYLU_SOLVER_AVAILABLE
    template class LS2G_API NRAlgo<LinearSolverPolicy<MyLULinearSolver>, MultiSlackNRSystem>;
    template class LS2G_API NRAlgo<LinearSolverPolicy<MyLULinearSolver>, SingleSlackNRSystem>;
    template class LS2G_API BaseDCAlgo<LinearSolverPolicy<MyLULinearSolver>>;
    template class LS2G_API BaseFDPFAlgo<LinearSolverPolicy<MyLULinearSolver>, FDPFMethod::XB>;
    template class LS2G_API OuterLoopAlgo<NROuterInner<LinearSolverPolicy<MyLULinearSolver>>>;
    #endif

The fast-decoupled algorithms **must** be instantiated there: part of their code
(``BaseFDPFAlgo::fillBp_Bpp``) is only defined in that file.

Registering them
++++++++++++++++

Finally, give each algorithm a name in ``register_builtin_solvers``
(``src/core/BuiltinSolversRegistration.cpp``):

.. code-block:: cpp

    #ifdef MYLU_SOLVER_AVAILABLE
        add_builtin("NR_MyLU",      []{ return std::make_unique<NR_MyLU>(); });
        add_builtin("NRSing_MyLU",  []{ return std::make_unique<NRSing_MyLU>(); });
        add_builtin("DC_MyLU",      []{ return std::make_unique<DC_MyLU>(); });
        add_builtin("FDPF_XB_MyLU", []{ return std::make_unique<FDPF_XB_MyLU>(); });
        add_builtin("NROuter_MyLU", []{ return std::make_unique<NROuter_MyLU>(); });
    #endif

.. note::
    There is no need to add a value to the ``AlgorithmType`` enum. KLU, NICSLU and CKTSO have one
    for historical reasons, and adding one means editing ``AlgorithmSelector.{hpp,cpp}``,
    ``AlgorithmTypeNames.cpp`` and ``binding_enums.cpp`` as well. Recent built-ins
    (``NRRefactorRetry_*``, ``*_KLUStaticPivot``) are registered by name only, as every plugin is:
    ``get_algo_type()`` reports them as ``AlgorithmType.Custom``, ``get_algo_name()`` gives their
    name, and they work everywhere else.

Compiling and using it
++++++++++++++++++++++

Build lightsim2grid from source as described in :ref:`page_install_from_source`, with the variable
set:

.. code-block:: bash

    git submodule update --init --depth 1 eigen SuiteSparse
    export PATH_MYLU=/path/to/mylu   # holds include/mylu.h and lib/libmylu.so
    pip install -e .

The configure step prints ``MyLU enabled: ...``; if it prints the warning instead, the library was
not found and the build carries on without it. The algorithms are then used by name exactly as in
the plugin route, without the ``load_algorithm_plugin`` call.

Exposing it to Python
+++++++++++++++++++++

Selecting by name is all most uses need. Three optional extras, each following what the code
already does for CKTSO:

- **A Python class** (``lightsim2grid.algorithm.NR_MyLU``), to build and drive the algorithm
  directly on a ``Ybus`` / ``Sbus`` without a grid: add a ``py::class_<NR_MyLU>`` block to
  ``src/bindings/python/binding_solvers.cpp`` (``bind_algo_methods``, ``bind_nr_algo_policies``,
  ``bind_linear_solver_stats``, and ``bind_outer_algo`` for ``NROuter_``; the DC and FDPF blocks bind
  fewer), give it a docstring in ``src/core/help_fun_msg.{hpp,cpp}`` (``DocSolver``), and import it in
  ``lightsim2grid/algorithm/__init__.py`` inside a ``try`` / ``except``, since it does not exist in a
  build without the library.
- **A compilation flag** (``lightsim2grid.compilation_options.mylu_solver_available``): see how
  ``cktso_solver_available`` is set in ``src/bindings/python/binding_module.cpp`` and read in
  ``lightsim2grid/compilation_options/__init__.py``.
- **Read the Docs**: the ``_READ_THE_DOCS`` stand-in classes at the end of ``KLUSolver.hpp`` /
  ``CKTSOSolver.hpp`` let Sphinx document the algorithms without the library. Only needed if the
  documentation of your fork must show them.

Choosing the names
~~~~~~~~~~~~~~~~~~

In both routes the name is how an algorithm is selected, and it is **stored** when a grid is saved
(pickle or :ref:`binary-serialization`). Keep to the ``<algorithm>_<LinearSolver>`` convention of
:ref:`algorithm_names`, stick to ASCII letters, digits and ``_``, and do not rename an algorithm
once grids have been saved with it. A plugin distributed to others should prefix its names with
something of its own (see :ref:`solver_plugin`).

Testing it
----------

The cheapest test that says something is to compare each new algorithm with its ``SparseLU``
counterpart on the same grid: the linear solver must not change the result, only the time it takes
to get it. Do it for a single power flow and for a ``ContingencyAnalysisCPP`` with ``add_all_n1()``
on a grid where one contingency splits it: that run calls ``analyze`` once, then ``refactorize`` for
every row, and it puts explicit zeros and masked rows in the matrix, which is where a wrapper gets
things wrong. Check ``get_linear_solver_stats()`` as well: a contingency analysis that reports more
than one ``analyze`` means the pattern is changing, which is a bug.

``examples/linear_solver_plugin/test_plugin.py`` does exactly this, and is a template to copy for
either route. In a fork, the C++ tests (Catch2, ``src/tests/``) build standalone, without Python,
and pick the solver up as long as ``PATH_MYLU`` is set when CMake configures (see
:ref:`cpp_library`); the ``unittest`` files of ``lightsim2grid/tests`` that compare solvers
(``test_KLUSolver.py``, ``test_CKTSOSolver.py``) are another template, skipping the test when
``"NR_MyLU" not in grid.available_algorithm_names()``.

And then
--------

Once the plugin is loaded, or the fork installed, the new algorithms are used like the built-in
ones: load a grid from any of the sources lightsim2grid reads (:ref:`network-init-formats`), select the algorithm
by name, and run a load flow or a batch.

.. code-block:: python

    import json
    import numpy as np
    import lightsim2grid
    from lightsim2grid.network import (init_from_pandapower, init_from_pypowsybl,
                                       init_from_matpower, init_from_powermodels)
    from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

    # plugin route only: in every process, before the algorithm is selected
    lightsim2grid.load_algorithm_plugin("examples/linear_solver_plugin/build/liblinear_solver_plugin.so")

    # 1. load a grid, whichever the source
    import pandapower.networks as pn
    grid = init_from_pandapower(pn.case118())          # a pandapower network

    import pypowsybl as pp
    grid = init_from_pypowsybl(pp.network.create_ieee118())   # a pypowsybl network

    grid = init_from_matpower("case118.m")             # a MATPOWER case file (.m or .mat)

    with open("network.json") as f:                    # a PowerModels.jl network data dictionary
        grid = init_from_powermodels(json.load(f))

    # 2. select the algorithm that uses the new linear solver
    grid.change_algorithm("NR_MyLU")
    print(grid.get_algo_name())                        # "NR_MyLU"

    # 3. a load flow
    V0 = np.ones(grid.total_bus(), dtype=complex)
    V = grid.ac_pf(V0, 10, 1e-8)                       # empty if it did not converge
    p_or, q_or, v_or, a_or = grid.get_line_res1()      # flows at side 1 of the lines

    # 4. a contingency analysis: every line and transformer disconnected in turn
    ca = ContingencyAnalysisCPP(grid)
    ca.change_algorithm("NR_MyLU")
    ca.handle_disconnected_grid = True                 # solve the rows that cut buses off too
    ca.add_all_n1()
    ca.compute(V0, 10, 1e-8)
    voltages = ca.get_voltages()                       # one row per contingency
    amps = ca.compute_flows()                          # current flows, one row per contingency
    print(ca.get_linear_solver_stats())                # one analyze for the whole batch

In grid2op, pass the name to the backend instead: ``LightSimBackend(algo_type="NR_MyLU")``.
:ref:`network-init-formats` describes the loaders and their options, and :doc:`security_analysis`
the contingency analysis.
