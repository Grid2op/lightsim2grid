.. _add_linear_solver:

Adding a linear solver (from source)
====================================

Every power-flow algorithm in lightsim2grid that needs to solve :math:`Ax = b` (Newton-Raphson, DC,
fast-decoupled, and the batch classes built on them) gets its linear algebra from a **linear solver**.
The ones shipped today are Eigen's ``SparseLU`` (always available), KLU (bundled from SuiteSparse),
and NICSLU and CKTSO (optional, found at build time). They all live in ``src/core/linear_solvers/``.

This page explains how to add another one, for example a vendor library, a GPU solver, or a solver
you wrote yourself, by **cloning lightsim2grid, adding a few files and compiling it from source**. It
works whether the library exposes a C API or a C++ one.

.. note::
    This is the *in-tree* route. If you only want a different *algorithm* (not a different linear
    solver), the run-time plugin mechanism described in :ref:`solver_plugin` needs no fork at all. See
    :ref:`add_linear_solver_plugin_route` at the end of this page for how the two relate.

The terms used below follow :ref:`algorithm_names`: a *linear solver* is the ``XxxLinearSolver``
class, an *algorithm* is what wraps it (``NR_Xxx``, ``DC_Xxx``, ...), and an algorithm is selected
by its *name*.

Overview
--------

A linear solver is **not** a virtual class: it is a template parameter of the algorithms, so its
type is known at compile time and no call goes through a vtable. Adding one therefore means:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - File
     - What to add
   * - ``src/core/linear_solvers/XxxSolver.{hpp,cpp}``
     - the wrapper class ``XxxLinearSolver`` (`The linear solver contract`_)
   * - ``src/core/CMakeLists.txt``
     - find the library, link it, define ``XXX_SOLVER_AVAILABLE`` (`Building it`_)
   * - ``src/core/linear_solvers/CMakeLists.txt``
     - compile ``XxxSolver.cpp`` when the library was found
   * - ``src/core/Solvers.hpp``
     - the algorithm aliases (``NR_Xxx``, ...) and their ``extern template`` declarations
   * - ``src/core/SolverInstantiations.cpp``
     - the explicit instantiations of those aliases
   * - ``src/core/BuiltinSolversRegistration.cpp``
     - register each alias under a name

These six files are enough: once they are done, the new algorithms can be selected by name from C++
and from Python (``grid.change_algorithm("NR_Xxx")``, ``LightSimBackend(algo_type="NR_Xxx")``,
``ContingencyAnalysis``, ``TimeSeries``, ...). `Exposing it to Python`_ lists the optional extras.

The running example on this page wraps a fictional C library, ``mylu``, whose header is:

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
NICSLU, CKTSO, cuDSS, ...) expose this shape: a handle, a symbolic analysis of the sparsity pattern,
a numeric factorization, a cheaper re-factorization that reuses the analysis, and a solve. A C++
library is wrapped in exactly the same way: the only difference is that you call its methods instead
of free functions.

The linear solver contract
--------------------------

The algorithms use the class through ``LinearSolverPolicy<XxxLinearSolver>`` (see
``linear_solvers/LinearSolverPolicy.hpp``), a thin wrapper that counts and times every call and adds
the refactorization fallback described below. The class must provide:

.. code-block:: cpp

    class LS2G_API XxxLinearSolver final {
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

The matrix you are given
~~~~~~~~~~~~~~~~~~~~~~~~

``J`` is an ``Eigen::SparseMatrix<double>`` seen through an ``Eigen::Ref``. What that means in
practice:

- it is **square** and stored in **compressed sparse column** (CSC) format, with **0-based int**
  indices: ``J.outerIndexPtr()`` (size ``n + 1``) is the column pointer, ``J.innerIndexPtr()`` the
  row indices and ``J.valuePtr()`` the values, both of size ``J.nonZeros()``. Within a column, the
  row indices are sorted;
- it may hold **explicit zeros**. The batch algorithms keep one sparsity pattern for a whole batch
  and change a bus's role by rewriting values in it: a value can be ``0`` in one call and non-zero in
  the next. A library that drops zero entries during its analysis cannot be used as is; give it the
  pattern, not the values, or turn that pruning off;
- between an ``analyze`` and the next ``reset`` / ``analyze``, the **pattern does not change**: only
  the values do. The memory address of the arrays is not guaranteed to be the same, so read the
  pointers again at each call rather than keeping them;
- it is not symmetric in general (a Newton-Raphson Jacobian is not), and it is not guaranteed to be
  diagonally dominant.

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

The methods
~~~~~~~~~~~

The Newton-Raphson algorithm calls them in this order (``NRAlgo.tpp``); DC and fast-decoupled do the
same with fewer refactorizations::

    reset()                          # new topology, or the algorithm was reset
    analyze(J); factorize(J)         # first iteration: new pattern
    solve(b)
    refactorize(J); solve(b)         # next iterations: same pattern, new values
    refactorize(J); solve(b)
    ...
    reset()                          # and so on, at the next topology change

``reset()``
    Free everything the previous analysis built, so that the next call can be an ``analyze``. It is
    also called before any ``analyze`` was done, and it can be called twice in a row: it must not
    fail in either case.

``analyze(J)``
    The symbolic step: ordering, symbolic factorization, allocation. Only the pattern of ``J`` is
    meaningful at this point (its values are those of the first iteration, but do not rely on them).
    It is the expensive call, and avoiding it is what most of lightsim2grid's performance work is
    about: a contingency analysis of thousands of rows calls it **once**.

``factorize(J)``
    A numeric factorization that **chooses its pivots**. It is always called right after ``analyze``,
    and again when the algorithm has a reason to think the previous pivots are no longer good.

``refactorize(J)``
    A numeric factorization with the same pattern, that may **reuse the pivot order** of the last
    ``factorize``. This is the call that makes a solver fast in lightsim2grid. If the library has no
    such call, implement ``refactorize`` as ``factorize`` (Eigen's ``SparseLU`` does exactly that).

    When the reused pivots turn out to be bad (a zero or a tiny pivot), return
    ``ErrorType::SolverReFactor`` rather than a poor factorization: when the fallback is on
    (``set_refactor_fallback``; the batch algorithms turn it on), ``LinearSolverPolicy`` then calls
    ``factorize`` on the same matrix, which chooses new pivots.

``solve(b)``
    Solve :math:`J x = b` with the last factorization, **in place**: on return ``b`` holds ``x``. If
    the library cannot work in place, solve into a temporary and copy it back.

``solve_transpose(b)``
    Solve :math:`J^T x = b` with the factorization of ``J``, in place. Only used by the
    gradient-based batch tools (``BatchAdjoint``), which do not pick it up from a new linear solver
    anyway. If you do not implement it, return ``ErrorType::NotImplemented`` and set
    ``CAN_SOLVE_TRANSPOSE = false``: never solve :math:`J x = b` instead.

``CAN_SOLVE_MAT``
    Whether ``solve`` accepts a matrix right hand side. ``false`` is always correct.

Return ``ErrorType::NoError`` on success. On failure, return the code of the step that failed
(``SolverAnalyze``, ``SolverFactor``, ``SolverReFactor``, ``SolverSolve``); ``NotInitError`` for a
call made before ``analyze``; ``LicenseError`` if a licensed library refuses to start. The algorithm
stops at the first error and reports it (``get_error()``). **Do not throw** from these methods: a
failed factorization is an expected outcome (a diverging Newton, a split grid), and the algorithms
handle it through the return code.

Other rules:

- **No copy, no move.** Delete the four copy and move operations, as the existing solvers do. The
  algorithms hold the solver by value and never copy it.
- **Own what you allocate** with RAII (``std::unique_ptr`` with a custom deleter for a C handle), so
  that ``reset`` and the destructor cannot leak, and a destructor never throws.
- **One instance per algorithm, one thread per instance.** The batch classes run several algorithms
  in parallel, each with its own linear solver. Global state in the library (a licence check, a
  thread pool, a static work buffer) must tolerate several instances alive at the same time.
- **C++14.** The core is compiled as C++14 (and C++26 in CI): no ``if constexpr``, ``std::optional``
  or structured bindings in the wrapper.
- The class is marked ``LS2G_API`` (from ``Utils.hpp``) so that it is exported from the shared
  library on Windows.

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

Writing the wrapper
-------------------

The whole class is guarded by ``MYLU_SOLVER_AVAILABLE``, which the build defines only when the
library was found, so that lightsim2grid still compiles without it.

``src/core/linear_solvers/MyLUSolver.hpp``:

.. code-block:: cpp

    #ifdef MYLU_SOLVER_AVAILABLE
    #ifndef MYLUSOLVER_H
    #define MYLUSOLVER_H

    #include <memory>

    #include "Utils.hpp"
    #include "Eigen/Core"
    #include "Eigen/SparseCore"

    // the library's C header; wrap it in extern "C" { } if it has no guard of its own
    #include "mylu.h"

    namespace ls2g {

    class LS2G_API MyLULinearSolver final
    {
        public:
            MyLULinearSolver() noexcept = default;
            ~MyLULinearSolver() noexcept = default;  // handle_ frees itself

            ErrorType reset();
            ErrorType analyze(const EigenRefConstRealSpMat & J);
            ErrorType factorize(const EigenRefConstRealSpMat & J);
            ErrorType refactorize(const EigenRefConstRealSpMat & J);
            ErrorType solve(Eigen::Ref<RealVect> b);
            ErrorType solve_transpose(Eigen::Ref<RealVect> b);

            static constexpr bool CAN_SOLVE_MAT = false;
            static constexpr bool CAN_SOLVE_TRANSPOSE = false;

            MyLULinearSolver(const MyLULinearSolver &) = delete;
            MyLULinearSolver(MyLULinearSolver &&) = delete;
            MyLULinearSolver & operator=(const MyLULinearSolver &) = delete;
            MyLULinearSolver & operator=(MyLULinearSolver &&) = delete;

        private:
            struct HandleDeleter {
                void operator()(mylu_handle * h) const noexcept { mylu_destroy(h); }
            };
            std::unique_ptr<mylu_handle, HandleDeleter> handle_;
    };

    }  // namespace ls2g

    #endif  // MYLUSOLVER_H
    #endif  // MYLU_SOLVER_AVAILABLE

``src/core/linear_solvers/MyLUSolver.cpp``:

.. code-block:: cpp

    #include "MyLUSolver.hpp"

    namespace ls2g {

    ErrorType MyLULinearSolver::reset(){
        handle_.reset();  // the next call is an analyze()
        return ErrorType::NoError;
    }

    ErrorType MyLULinearSolver::analyze(const EigenRefConstRealSpMat & J){
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

    ErrorType MyLULinearSolver::factorize(const EigenRefConstRealSpMat & J){
        if(!handle_) return ErrorType::NotInitError;
        if(mylu_factor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverFactor;
        return ErrorType::NoError;
    }

    ErrorType MyLULinearSolver::refactorize(const EigenRefConstRealSpMat & J){
        if(!handle_) return ErrorType::NotInitError;
        if(mylu_refactor(handle_.get(), J.valuePtr()) != 0) return ErrorType::SolverReFactor;
        return ErrorType::NoError;
    }

    ErrorType MyLULinearSolver::solve(Eigen::Ref<RealVect> b){
        if(!handle_) return ErrorType::NotInitError;
        if(mylu_solve(handle_.get(), b.data()) != 0) return ErrorType::SolverSolve;
        if(!b.allFinite()) return ErrorType::InifiniteValue;
        return ErrorType::NoError;
    }

    ErrorType MyLULinearSolver::solve_transpose(Eigen::Ref<RealVect> /*b*/){
        return ErrorType::NotImplemented;  // see CAN_SOLVE_TRANSPOSE
    }

    }  // namespace ls2g

To check that the wrapper compiles, without a full build (from the repository root, with the
submodules fetched):

.. code-block:: bash

    g++ -fsyntax-only -std=c++14 -DLS2G_BUILDING_CORE -DMYLU_SOLVER_AVAILABLE \
        -I src -I src/core -I eigen -I /path/to/mylu/include \
        src/core/linear_solvers/MyLUSolver.cpp

Building it
-----------

The existing optional solvers are found through an environment variable pointing at the library's
root folder (``PATH_NICSLU``, ``PATH_CKTSO``): do the same.

In ``src/core/CMakeLists.txt``, find the library. This block must come **before**
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
by the Python bindings and the C++ tests too.

In ``src/core/linear_solvers/CMakeLists.txt``, compile the wrapper only when the library is there:

.. code-block:: cmake

    if(LIB_MYLU)
        list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/MyLUSolver.cpp")
    endif()

If the library comes as a CMake package, ``find_package(MyLU)`` and linking its imported target
works as well. If it is a shared library outside the system paths, the dynamic loader must find it
at run time (``LD_LIBRARY_PATH``, the rpath, or next to ``lightsim2grid_core``).

Declaring the algorithms
------------------------

A linear solver on its own is not selectable: what is selected is an *algorithm* using it. Give it
the same algorithms as the other solvers, or only those you need. In ``src/core/Solvers.hpp``,
include the header with the others and declare the aliases:

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
----------------

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

The name is how the algorithm is selected, and it is **stored** when a grid is saved (pickle or
:ref:`binary-serialization`): keep to the ``<algorithm>_<LinearSolver>`` convention of
:ref:`algorithm_names`, stick to ASCII letters, digits and ``_``, and do not rename it once grids
have been saved with it.

.. note::
    There is no need to add a value to the ``AlgorithmType`` enum. KLU, NICSLU and CKTSO have one
    for historical reasons, and adding one means editing ``AlgorithmSelector.{hpp,cpp}``,
    ``AlgorithmTypeNames.cpp`` and ``binding_enums.cpp`` as well. Recent built-ins
    (``NRRefactorRetry_*``, ``*_KLUStaticPivot``) are registered by name only, and so is every
    plugin: ``get_algo_type()`` reports them as ``AlgorithmType.Custom`` and they work everywhere
    else.

Compiling and using it
----------------------

Build lightsim2grid from source as described in :ref:`page_install_from_source`, with the variable
set:

.. code-block:: bash

    git submodule update --init --depth 1 eigen SuiteSparse
    export PATH_MYLU=/path/to/mylu   # holds include/mylu.h and lib/libmylu.so
    pip install -e .

The configure step prints ``MyLU enabled: ...``; if it prints the warning instead, the library was
not found and the build carries on without it. Then:

.. code-block:: python

    import numpy as np
    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

    grid = init_from_pandapower(pn.case118())
    print([nm for nm in grid.available_algorithm_names() if "MyLU" in nm])

    grid.change_algorithm("NR_MyLU")
    V0 = np.ones(grid.total_bus(), dtype=complex)
    V = grid.ac_pf(V0, 10, 1e-8)
    print(grid.get_algo().get_linear_solver_stats())   # analyze / factorize / refactorize / solve counts

    # in the batch classes
    ca = ContingencyAnalysisCPP(grid)
    ca.change_algorithm("NR_MyLU")

    # and in grid2op
    from lightsim2grid import LightSimBackend
    backend = LightSimBackend(algo_type="NR_MyLU")

The grid2op-level ``lightsim2grid.ContingencyAnalysis`` takes the same names in its
``change_algorithm``; for the other grid2op-level wrappers, call it on the C++ object they hold
(``ts.computer.change_algorithm("NR_MyLU")``).

Exposing it to Python
---------------------

Selecting by name, as above, is all most uses need. Three optional extras, each following what the
code already does for CKTSO:

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

Testing it
----------

The cheapest test that says something is to compare each new algorithm with its ``SparseLU``
counterpart on the same grid: the linear solver must not change the result, only the time it takes
to get it. Do it for a single power flow and for a ``ContingencyAnalysis`` with ``add_all_n1()`` on a
grid where one contingency splits it: that run calls ``analyze`` once, then ``refactorize`` for every
row, and it puts explicit zeros and masked rows in the matrix, which is where a wrapper gets things
wrong. Check ``get_linear_solver_stats()`` as well: a contingency analysis that reports one
``analyze`` per row means the pattern is changing, which is a bug.

The C++ tests (Catch2, ``src/tests/``) build standalone, without Python, and pick the solver up as
long as ``PATH_MYLU`` is set when CMake configures (see :ref:`cpp_library`). In Python, the
``unittest`` files of ``lightsim2grid/tests`` that compare solvers (``test_KLUSolver.py``,
``test_CKTSOSolver.py``) are a template to copy; skip the test when
``"NR_MyLU" not in grid.available_algorithm_names()``.

.. _add_linear_solver_plugin_route:

Without forking: what exists today
----------------------------------

lightsim2grid has no run-time plugin mechanism for *linear solvers*, only for *algorithms*
(:ref:`solver_plugin`). Because a linear solver is a template parameter of the algorithms, an
algorithm plugin can, in principle, carry its own linear solver: define ``MyLULinearSolver`` in the
plugin, and register ``NRAlgo<LinearSolverPolicy<MyLULinearSolver>, MultiSlackNRSystem>`` in it,
instantiated from the installed headers. ``examples/lm_algorithm`` already instantiates a Newton
algorithm from those headers, with an existing linear solver.

No example covers a new linear solver that way yet, and two limits are known: the fast-decoupled
algorithms cannot be instantiated outside ``lightsim2grid_core`` (see `Declaring the algorithms`_),
and the plugin must be compiled with the same flags as the core (see
:ref:`solver_plugin`, *Matching build flags*). The in-tree route on this page has neither limit.
