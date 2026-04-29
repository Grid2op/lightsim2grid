.. _solver_plugin:

External Amgorithm Plugins
===========================

LightSim2grid supports dynamically-loaded algorithm plugins.  A plugin is a
shared library (``.so`` / ``.dll``) that registers one or more custom
"powerflow solvers" at load time.  Once loaded, those solvers / algorithms
behave exactly like the built-in ones: they are accessible by name, selectable via
:func:`GridModel.change_algorithm`, and appear in
:func:`GridModel.available_algorithm_names`.

This mechanism lets you add a new solver algorithm — from your own
repository or a third-party library — **without modifying lightsim2grid's
source code**.

How the registry works
-----------------------

All solvers are stored in a process-wide singleton called
``AlgorithmRegistry``.  On startup the built-in algorithm (NR_SparseLU, NR_KLU,
GaussSeidel, DC_KLU, …) are registered.  

A plugin library extends this table at
``dlopen``/``LoadLibrary`` time by placing a static ``SolverRegistrar``
object in one of its translation units.  The registrar's constructor fires
automatically when the library is mapped into the process, calling
``SolverRegistry::instance().register_solver(name, factory)`` before any
Python code can observe the new library.

The lookup flow is:

.. code-block:: text

    Python: lightsim2grid.load_algorithm_plugin("path/to/plugin.so")
      └─ ctypes.CDLL(..., RTLD_GLOBAL)          # dlopen fires static ctors
           └─ SolverRegistrar { "MySolver", factory }
                └─ SolverRegistry::instance().register_solver(...)

    Python: grid.change_algorithm("MySolver")
      └─ AlgorithmSelector::change_algorithm("MySolver")
           └─ SolverRegistry::instance().make("MySolver")
                └─ factory()  →  unique_ptr<BaseAlgo>

The ``AlgorithmRegistry`` C++ API (defined in ``AlgorithmRegistry.hpp``, installed to ``include/lightsim2grid/``):

.. code-block:: cpp

    class AlgorithmRegistry {
    public:
        // Meyers singleton — one instance per process.
        static AlgorithmRegistry& instance();

        // Register a factory under a name.  Called by SolverRegistrar.
        void register_solver(const std::string& name, Factory f);

        // Instantiate a solver by name.  Throws std::invalid_argument if
        // the name is unknown.
        std::unique_ptr<BaseAlgo> make(const std::string& name) const;

        bool is_registered(const std::string& name) const;

        // List of all registered names (built-in + plugins).
        std::vector<std::string> available_algorithms() const;
    };

    // Drop a static instance of this in an anonymous namespace to
    // register your solver when the .so is loaded — no macro needed.
    class AlgorithmRegistrar {
    public:
        AlgorithmRegistrar(const std::string& name, SolverRegistry::Factory f);
    };


The ``BaseAlgo`` C++ interface
-------------------------------

Every solver — built-in or plugin — must publicly inherit from
``ls2g::BaseAlgo`` (defined in ``powerflow_algorithm/BaseAlgo.hpp``,
installed to ``include/lightsim2grid/powerflow_algorithm/``).

Constructor
~~~~~~~~~~~

.. code-block:: cpp

    explicit BaseAlgo(bool is_ac = true);

Pass ``true`` for AC solvers and ``false`` for DC solvers.  This value is
stored in the public member ``IS_AC``, which GridModel uses to route
:func:`change_solver` calls to the right slot (AC or DC).

Methods to override
~~~~~~~~~~~~~~~~~~~~

``compute_pf`` *(pure virtual — you must implement this)*

.. code-block:: cpp

    virtual bool compute_pf(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        CplxVect& V,
        const CplxVect& Sbus,
        Eigen::Ref<const IntVect> slack_ids,
        const RealVect& slack_weights,
        Eigen::Ref<const IntVect> pv,
        Eigen::Ref<const IntVect> pq,
        int max_iter,
        real_type tol);

Solve the power-flow problem ``V·(Ybus·V)* = Sbus``.

+-------------------+------------------------------------------------------+
| Parameter         | Meaning                                              |
+===================+======================================================+
| ``Ybus``          | Complex bus admittance matrix (sparse, n×n).         |
+-------------------+------------------------------------------------------+
| ``V``             | On input: initial voltage phasors.  On output:       |
|                   | solved voltage phasors (in-place update).            |
+-------------------+------------------------------------------------------+
| ``Sbus``          | Complex bus power injections (loads + generators).   |
+-------------------+------------------------------------------------------+
| ``slack_ids``     | Indices of slack (reference) buses.                  |
+-------------------+------------------------------------------------------+
| ``slack_weights`` | Per-slack participation factors (sum to 1).          |
+-------------------+------------------------------------------------------+
| ``pv``            | Indices of PV buses (voltage-controlled generators). |
+-------------------+------------------------------------------------------+
| ``pq``            | Indices of PQ buses (constant power loads).          |
+-------------------+------------------------------------------------------+
| ``max_iter``      | Maximum number of Newton-Raphson iterations.         |
+-------------------+------------------------------------------------------+
| ``tol``           | Convergence tolerance (per-unit mismatch).           |
+-------------------+------------------------------------------------------+

Return ``true`` on convergence; store results in the protected members
listed below and set ``err_ = ErrorType::NoError`` (or an appropriate
error code on failure).

``set_gridmodel`` *(virtual — override if you need grid data)*

.. code-block:: cpp

    virtual void set_gridmodel(const GridModel* gridmodel);

Called by ``ChooseAlgorithm`` after the solver is activated (and again after
every ``change_algorithm`` call).  The default implementation stores the
pointer in the protected member ``gridmodel_ptr_``.  Override only if your
algorithm / solver needs to cache additional data derived from the grid topology.

``reset`` *(virtual — override if you carry extra state)*

.. code-block:: cpp

    virtual void reset();

Called whenever the solver is swapped out or the grid topology changes.
The base implementation clears all result vectors and resets timers.
Call ``BaseAlgo::reset()`` from your override if you want that baseline
behaviour, then clear your own state.

``get_J`` *(virtual — override only for Newton-Raphson algorithms)*

.. code-block:: cpp

    virtual Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const;

Returns the last Jacobian matrix.  The base implementation throws
``std::runtime_error``.  Override if your solver exposes a Jacobian
(needed only when Python code calls ``get_J()`` on the grid model).

Protected result members
~~~~~~~~~~~~~~~~~~~~~~~~~

Your ``compute_pf`` implementation must populate these fields so that
``GridModel`` can extract flows and inject the results back into the grid
state.

+--------------------+------------------------------------------------------+
| Member             | Content                                              |
+====================+======================================================+
| ``V_``             | ``CplxVect`` — solved complex voltages               |
+--------------------+------------------------------------------------------+
| ``Va_``            | ``RealVect`` — voltage angles (radians)              |
+--------------------+------------------------------------------------------+
| ``Vm_``            | ``RealVect`` — voltage magnitudes (p.u.)             |
+--------------------+------------------------------------------------------+
| ``n_``             | ``int`` — number of buses                            |
+--------------------+------------------------------------------------------+
| ``nr_iter_``       | ``int`` — iterations performed                       |
+--------------------+------------------------------------------------------+
| ``err_``           | ``ErrorType`` — ``NoError`` on success               |
+--------------------+------------------------------------------------------+

Read-only accessors (provided by ``BaseAlgo``, no override needed):
``get_Va()``, ``get_Vm()``, ``get_V()``, ``get_error()``,
``get_nb_iter()``, ``converged()``, ``get_timers()``.


Writing a solver plugin
------------------------

**1 — Create the solver class and register it**

Place this in a single ``.cpp`` file.  The anonymous-namespace static
object ensures the registration fires exactly once, at ``dlopen`` time.

.. code-block:: cpp

    // my_solver_plugin.cpp
    #include <AlgorithmRegistry.hpp>
    #include <powerflow_algorithm/BaseAlgo.hpp>

    class MySolver : public ls2g::BaseAlgo {
    public:
        MySolver() : ls2g::BaseAlgo(/*is_ac=*/true) {}

        bool compute_pf(
            const Eigen::SparseMatrix<ls2g::cplx_type>& Ybus,
            ls2g::CplxVect& V,
            const ls2g::CplxVect& Sbus,
            Eigen::Ref<const ls2g::IntVect> slack_ids,
            const ls2g::RealVect& slack_weights,
            Eigen::Ref<const ls2g::IntVect> pv,
            Eigen::Ref<const ls2g::IntVect> pq,
            int max_iter,
            ls2g::real_type tol) override
        {
            // ... your algorithm here ...

            // populate result fields when done:
            V_      = V;             // solved voltages (in-place already updated)
            Va_     = V.array().arg();
            Vm_     = V.array().abs();
            n_      = static_cast<int>(V.size());
            nr_iter_= 1;
            err_    = ls2g::ErrorType::NoError;
            return true;
        }
    };

    // Self-registration — fires when the .so is dlopen'd.
    namespace {
        ls2g::ALgorithmRegistrar _reg(
            "MySolver",
            []{ return std::unique_ptr<ls2g::BaseAlgo>(new MySolver()); }
        );
    }

**2 — Write a CMakeLists.txt**

The recommended approach uses ``find_package(lightsim2grid_core)`` to locate the
installed headers and library.  A source-tree fallback is provided for in-repo
development without a full ``pip install``.

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.15)
    project(my_solver CXX)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # Strategy 1 — installed package (preferred).
    # Get the hint path via:
    #   python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())"
    # then pass it as:
    #   cmake -DLIGHTSIM2GRID_CMAKE_DIR=<path> ...
    find_package(lightsim2grid_core CONFIG QUIET
        HINTS "${LIGHTSIM2GRID_CMAKE_DIR}")

    # Strategy 2 — source tree (in-repo development without pip install).
    if(NOT lightsim2grid_core_FOUND)
        if(NOT DEFINED LIGHTSIM2GRID_SRC)
            set(LIGHTSIM2GRID_SRC "/path/to/lightsim2grid/src/core")
        endif()
        if(NOT DEFINED Eigen3_INCLUDE)
            set(Eigen3_INCLUDE "/path/to/lightsim2grid/eigen")
        endif()
        if(NOT EXISTS "${LIGHTSIM2GRID_SRC}/SolverRegistry.hpp")
            message(FATAL_ERROR
                "lightsim2grid_core not found.\n"
                "Install lightsim2grid and pass:\n"
                "  -DLIGHTSIM2GRID_CMAKE_DIR=<cmake-dir>\n"
                "or pass -DLIGHTSIM2GRID_SRC=<path/to/src/core> for a source build.")
        endif()
        add_library(lightsim2grid_core_iface INTERFACE)
        target_include_directories(lightsim2grid_core_iface INTERFACE
            "${LIGHTSIM2GRID_SRC}" "${Eigen3_INCLUDE}")
        add_library(lightsim2grid::core ALIAS lightsim2grid_core_iface)
    endif()

    # Build as a MODULE (dlopen-able at runtime, not linked at build time).
    add_library(my_solver MODULE my_solver_plugin.cpp)
    target_link_libraries(my_solver PRIVATE lightsim2grid::core)

    if(WIN32)
        # find_package provides the import lib via the IMPORTED target.
    elseif(UNIX AND NOT APPLE)
        if(NOT lightsim2grid_core_FOUND)
            target_link_options(my_solver PRIVATE -Wl,--allow-shlib-undefined)
        endif()
        set_target_properties(my_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
    elseif(APPLE)
        if(NOT lightsim2grid_core_FOUND)
            target_link_options(my_solver PRIVATE -undefined dynamic_lookup)
        endif()
        set_target_properties(my_solver PROPERTIES PREFIX "lib" SUFFIX ".so")
    endif()

**3 — Build**

After ``pip install lightsim2grid``, obtain the CMake directory and pass it:

.. code-block:: bash

    LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")

    mkdir build && cd build
    cmake .. -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
    make

Or without a pip install (source tree, Linux / macOS):

.. code-block:: bash

    mkdir build && cd build
    cmake .. \
        -DLIGHTSIM2GRID_SRC=/path/to/lightsim2grid/src/core \
        -DEigen3_INCLUDE=/path/to/lightsim2grid/eigen
    make

Windows (MSVC, from a Developer Command Prompt):

.. code-block:: bat

    for /f "delims=" %i in ('python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())"') do set LS2G_CMAKE=%i

    mkdir build && cd build
    cmake .. -DLIGHTSIM2GRID_CMAKE_DIR="%LS2G_CMAKE%"
    cmake --build . --config Release

The resulting file is ``Release\my_solver.dll`` on Windows (no ``lib`` prefix)
and ``libmy_solver.so`` on Linux / macOS.
Pass that path to :func:`~lightsim2grid.load_algorithm_plugin`.


Loading and using the plugin from Python
-----------------------------------------

.. code-block:: python

    import lightsim2grid
    from lightsim2grid.lightsim2grid_cpp import GridModel

    # 1. Load the plugin — this fires the C++ static constructors, which
    #    register "MySolver" into the SolverRegistry singleton.
    lightsim2grid.load_algorithm_plugin("build/libmy_solver.so")

    # 2. Confirm the solver is available.
    gm = GridModel()
    print(gm.available_algorithm_names())   # [..., "MySolver", ...]

    # 3. Activate the solver.
    gm.change_algorithm("MySolver")

    # 4. Run a powerflow — lightsim2grid now delegates to MySolver.compute_pf().
    # (set up the grid first via gm.init_from_pandapower() or similar)
    gm.run_pf(...)


Python API reference
---------------------

.. py:function:: lightsim2grid.load_algorithm_plugin(path: str) -> None

    Load a shared library containing one or more lightsim2grid solver / algorithm
    plugins.

    The library must contain at least one static ``AlgorithmRegistrar``
    object in an anonymous namespace (see the example above).  Its
    constructor fires at ``dlopen`` time, registering the new solver name
    into the ``AlgorithmRegistry`` singleton.

    After this call the new solver is usable via
    ``grid.change_algorithm("MyAlgoName")`` and will appear in
    ``grid.available_algorithm_names()``.

    :param path: Absolute or relative path to the ``.so`` / ``.dll`` file.
    :raises OSError: If the library cannot be loaded (missing file,
        ABI mismatch, unresolved symbols, …).

.. py:method:: GridModel.change_algorithm(name: str) -> None

    Select the active solver by name.  The name must be one of the
    strings returned by :py:meth:`GridModel.available_algorithm_names`.

    For built-in solvers the enum overload is also available::

        gm.change_algorithm(AlgorithmType.NR_KLU)

    :param name: Registered solver name (case-sensitive).
    :raises RuntimeError: If *name* is not registered.

.. py:method:: GridModel.available_algorithm_names() -> list[str]

    Return all solver names currently registered, including any that were
    added via :func:`~lightsim2grid.load_solver_plugin`.

    .. code-block:: python

        >>> gm.available_algorithm_names()
        ['NR_SparseLU', 'GaussSeidel', 'NRSing_KLU', 
         # and all other lightsim2grid installed "solver"
         'MySolver']

.. note::

    :attr:`AlgorithmType.Custom` is the enum value assigned to any solver
    loaded via a plugin.  ``gm.get_solver_type()`` returns
    ``AlgorithmType.Custom`` whenever the active solver was registered
    through the plugin mechanism.

.. warning::

    The plugin (``.so`` / ``.dll``) must be compiled against the **same
    version** of lightsim2grid headers that is installed at runtime.  ABI
    mismatches (different ``BaseAlgo`` layout, different Eigen version) will
    cause undefined behaviour or a load-time error.

    On Windows the plugin also links against ``lightsim2grid_cpp.lib``.  This
    import library must match the ``lightsim2grid_cpp.pyd`` that will be loaded
    at runtime — i.e. both must come from the same lightsim2grid build.


Worked example (``examples/external_algorithm/``)
-----------------------------------------------

A minimal but complete example lives in the repository under
``examples/external_algorithm/``.  It implements ``DummyExternalAlgo``:
an AC solver that always "converges" on the first call by returning the
initial voltage vector unchanged — useful as a smoke test for the plugin
mechanism.

Build and run (after ``pip install lightsim2grid``):

.. code-block:: bash

    LS2G_CMAKE=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")
    cd examples/external_algorithm
    cmake -S . -B build -DLIGHTSIM2GRID_CMAKE_DIR="$LS2G_CMAKE"
    cmake --build build
    python test_plugin.py

Or from the source tree without a pip install:

.. code-block:: bash

    cd examples/external_algorithm
    cmake -S . -B build   # falls back to ../../src/core and ../../eigen
    cmake --build build
    python test_plugin.py

Expected output::

    Plugin loaded successfully.
    Registered solvers: ['DC', 'DummyExternal', 'FDPF_BX_SparseLU', ...]
    change_algorithm('DummyExternal') OK — solver type is Custom as expected.
    All checks passed.

The automated regression test is in
``lightsim2grid/tests/test_solver_registry.py``, class
``TestPluginLoading``.  It is skipped automatically when the example
plugin has not been built, and passes once the ``.so`` is present.


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
