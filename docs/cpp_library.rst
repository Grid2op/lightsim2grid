.. _cpp_library:

Using lightsim2grid as a C++ library
=====================================

Since the split between the C++ core and the python bindings, everything that
computes in lightsim2grid lives in a standalone shared library called
``lightsim2grid_core``, with **no dependency on python, pybind11 or grid2op**.
The python package (``lightsim2grid_cpp``) is only a thin pybind11 wrapper
around it.

This page documents how to build that library, install it, and link your own
C++ code against it. If your goal is to plug a custom powerflow algorithm into
the *python* package, see :ref:`solver_plugin` instead -- this page is about
consuming the C++ API directly.

What the library contains
--------------------------

Everything is in the ``ls2g`` namespace, headers at the root of
``src/core/``. The main entry points:

- ``ls2g::LSGrid`` (``LSGrid.hpp``): the grid model itself -- built
  programmatically through the ``init_*`` methods (buses, powerlines, trafos,
  loads, generators, ...), solved with ``ac_pf`` / ``dc_pf``, inspected with
  the ``get_*_res`` methods. This is the C++ class exposed to python under
  the same name, ``LSGrid``.
- the element containers (``element_container/*.hpp``), time series and
  contingency analysis (``batch_algorithm/``), the powerflow algorithms
  (``powerflow_algorithm/``) and linear solvers (``linear_solvers/``).
- ``AlgorithmRegistry.hpp``: the registry used to add custom powerflow
  algorithms (see :ref:`solver_plugin`).
- ``BinaryArchive.hpp``: the fast binary serialization used by
  ``save_binary`` / ``load_binary`` (see :ref:`binary-serialization`).

The public API is C++14 (the library itself compiles with the newest standard
your compiler supports, see below). Eigen is part of the public interface
(vector arguments are ``Eigen::VectorXd`` / ``Eigen::VectorXi`` and friends);
the matching Eigen headers are shipped alongside the library so you do not
need your own copy.

Option 1: use the copy shipped with the python wheel
-----------------------------------------------------

Every ``pip install lightsim2grid`` (>= 1.0) already contains everything a
C++ consumer needs, inside the installed package directory:

- ``liblightsim2grid_core.so`` (or ``.dylib`` / ``.dll``) next to the python
  extension;
- the headers (including the bundled ``Eigen/``) under
  ``lightsim2grid/include/``;
- a CMake package config under
  ``lightsim2grid/share/cmake/lightsim2grid_core/``.

Two python helpers return those paths, so build scripts do not need to
hard-code anything:

.. code-block:: python

    import lightsim2grid
    lightsim2grid.get_include()    # .../site-packages/lightsim2grid/include
    lightsim2grid.get_cmake_dir()  # .../site-packages/lightsim2grid/share/cmake/lightsim2grid_core

Option 2: build and install the core from source
-------------------------------------------------

``src/core`` is a self-contained CMake project (python is never required):

.. code-block:: bash

    git clone https://github.com/Grid2op/lightsim2grid.git
    cd lightsim2grid
    git submodule update --init eigen        # required (header-only)
    git submodule update --init SuiteSparse  # optional, for the KLU solver

    cmake -S src/core -B build_core -DCMAKE_BUILD_TYPE=Release
    cmake --build build_core -j
    cmake --install build_core --prefix /where/you/want/it

The install lays out a standard prefix: ``lib/liblightsim2grid_core.so``,
``include/lightsim2grid/*.hpp`` (with Eigen bundled next to them) and
``lib/cmake/lightsim2grid_core/`` (the package config).

About KLU: the build looks for the KLU sparse linear solver four ways (system
CMake config e.g. ``libsuitesparse-dev``, legacy system libraries, pre-built
static libraries inside the ``SuiteSparse`` submodule, and finally building it
from the submodule -- this last strategy needs CMake >= 3.22). **If none
succeeds the build does not fail**: the library simply falls back to the
always-available Eigen ``SparseLU`` algorithms (which are also the
``LSGrid`` defaults), and only the ``*_KLU`` entries of ``AlgorithmType``
disappear. When KLU is found, ``KLU_SOLVER_AVAILABLE`` is defined on the
target (and propagated to consumers) and the KLU code is baked into the
shared library.

A few environment variables tune the compilation, following the same
convention as the python build: ``__O3_OPTIM`` (``-O3``, on by default),
``__COMPILE_MARCHNATIVE`` (``-march=native``), ``__SANITIZE``
(ASan + UBSan) and ``__DEBUG_ASSERTS`` (``-UNDEBUG -D_GLIBCXX_ASSERTIONS``,
re-enabling the Eigen bounds assertions).

Linking against it with CMake
------------------------------

The package config exports a single imported target,
``lightsim2grid::core`` (with ``ls2g::core`` as a convenience alias). A
minimal consumer:

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.16)
    project(my_grid_tool LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 14)

    find_package(lightsim2grid_core CONFIG REQUIRED)

    add_executable(my_grid_tool main.cpp)
    target_link_libraries(my_grid_tool PRIVATE lightsim2grid::core)

The imported target carries the include directories (lightsim2grid headers
*and* the bundled Eigen), the ``Threads`` dependency, and the
``KLU_SOLVER_AVAILABLE`` definition when the library was built with KLU --
nothing else to configure.

Point CMake at the config with either of:

.. code-block:: bash

    # a source install (option 2): the prefix you installed into
    cmake -S . -B build -DCMAKE_PREFIX_PATH=/where/you/want/it

    # the python wheel (option 1): the exact config directory
    cmake -S . -B build \
        -Dlightsim2grid_core_DIR=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")

Alternatively, a super-project that vendors the lightsim2grid sources can
skip ``find_package`` entirely and pull the target in directly (this is what
the python bindings and the C++ unit tests do):

.. code-block:: cmake

    if(NOT TARGET lightsim2grid_core)
        add_subdirectory(path/to/lightsim2grid/src/core lightsim2grid_core_build)
    endif()
    target_link_libraries(my_grid_tool PRIVATE lightsim2grid_core)

One MSVC-specific note: the headers decorate the API with ``LS2G_API``
(dllimport when consuming, dllexport when building). This is automatic as
long as you *link* the library. Only if you compile core ``.cpp`` files
directly into your own target must you add the ``LS2G_BUILDING_CORE``
compile definition.

Minimal example
----------------

Build a 3-bus grid (slack generator, two lines, one load) and solve an AC
powerflow -- no input file, no python:

.. code-block:: c++

    #include <iostream>
    #include "LSGrid.hpp"

    int main()
    {
        ls2g::LSGrid grid;
        grid.set_sn_mva(100.);

        // three 138 kV buses (one busbar per substation)
        ls2g::RealVect bus_vn_kv(3);
        bus_vn_kv << 138., 138., 138.;
        grid.init_bus(3, 1, bus_vn_kv, 0, 0);

        // two lines 0-1 and 1-2, r/x/h given in pu on the sn_mva base
        ls2g::RealVect r(2), x(2);
        r << 0.01, 0.01;
        x << 0.10, 0.10;
        ls2g::CplxVect h = ls2g::CplxVect::Zero(2);
        Eigen::VectorXi from_id(2), to_id(2);
        from_id << 0, 1;
        to_id << 1, 2;
        grid.init_powerlines(r, x, h, from_id, to_id);

        // a 50 MW / 10 MVar load at bus 2
        ls2g::RealVect load_p(1), load_q(1);
        load_p << 50.;
        load_q << 10.;
        Eigen::VectorXi load_bus(1);
        load_bus << 2;
        grid.init_loads(load_p, load_q, load_bus);

        // a generator at bus 0 (1.02 pu), registered as the slack
        ls2g::RealVect gen_p(1), gen_v(1), min_q(1), max_q(1);
        gen_p << 0.;
        gen_v << 1.02;
        min_q << -1000.;
        max_q << 1000.;
        Eigen::VectorXi gen_bus(1);
        gen_bus << 0;
        grid.init_generators(gen_p, gen_v, min_q, max_q, gen_bus);
        grid.add_gen_slackbus(0, 1.);

        // Newton-Raphson with Eigen SparseLU (the default, KLU not needed)
        const ls2g::CplxVect Vinit =
            ls2g::CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.));
        const ls2g::CplxVect V = grid.ac_pf(Vinit, 20, 1e-8);

        if (V.size() == 0) {  // diverged: the result vector is empty
            std::cerr << "powerflow diverged" << std::endl;
            return 1;
        }
        for (Eigen::Index bus = 0; bus < V.size(); ++bus) {
            std::cout << "bus " << bus << ": " << std::abs(V(bus)) << " pu, "
                      << std::arg(V(bus)) << " rad" << std::endl;
        }
        const auto gen_res = grid.get_gen_res();  // (p_mw, q_mvar, v_kv)
        std::cout << "slack injects " << std::get<0>(gen_res)(0) << " MW" << std::endl;
        return 0;
    }

Conventions worth knowing (they mirror the python API, which shares this
code): ``init_powerlines`` takes r/x/h **in pu** on the ``sn_mva`` base;
loads are in MW / MVar; generator voltage setpoints in pu; ``ac_pf`` /
``dc_pf`` take a per-bus complex initial voltage of size ``total_bus()`` and
return the per-bus complex voltage, or an **empty vector when the powerflow
diverged** (details via ``grid.get_algo().converged()`` /
``get_error()`` / ``get_nb_iter()``). Errors (bad element ids, wrong
``Vinit`` size, invalid slack) are reported as ``std::runtime_error`` /
``std::out_of_range`` exceptions.

The C++ unit tests
-------------------

The library has its own test suite (Catch2, under ``src/tests/``) covering
the binary serialization layer and the ``LSGrid`` API above -- the
``test_lsgrid.cpp`` file doubles as a working example of building and
solving a grid from C++. It builds standalone in seconds, without python:

.. code-block:: bash

    git submodule update --init eigen Catch2   # SuiteSparse optional
    cmake -S src/tests -B build_tests -DCMAKE_BUILD_TYPE=Debug
    cmake --build build_tests -j
    ctest --test-dir build_tests --output-on-failure

    # the whole suite is a small plain binary: valgrind is practical
    valgrind --error-exitcode=1 --leak-check=full build_tests/lightsim2grid_unit_tests

The same suite can be enabled from the top-level (python) build with
``-DBUILD_TESTING=ON`` (off by default: wheels never build it), and runs in
CI on every push -- once through ctest and once under valgrind
(``.github/workflows/cpp_unit_tests.yml``).
