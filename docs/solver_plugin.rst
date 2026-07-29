.. _solver_plugin:

External Algorithm Plugins
==========================

LightSim2grid supports dynamically-loaded algorithm plugins.  A plugin is a
shared library (``.so`` / ``.dll``) that registers one or more custom
"powerflow solvers" at load time.  Once loaded, those solvers / algorithms
behave exactly like the built-in ones: they are accessible by name, selectable via
:func:`LSGrid.change_algorithm`, and appear in
:func:`LSGrid.available_algorithm_names`.

This mechanism lets you add a new solver algorithm — from your own
repository or a third-party library — **without modifying lightsim2grid's
source code**.

How the registry works
-----------------------

All solvers are stored in a process-wide singleton called
``AlgorithmRegistry``.  On startup the built-in algorithm (NR_SparseLU, NR_KLU,
GaussSeidel, DC_KLU, …) are registered.  

A plugin library extends this table through an exported entry point named
``ls2g_register_plugin``.  You do not write that function by hand: the
``LS2G_PLUGIN_ENTRY`` macro generates it from a small registration function you
supply.  When ``load_algorithm_plugin()`` loads the library it looks the entry
point up (``dlsym``/``GetProcAddress``) and calls it, staging your solvers and
committing them to the registry atomically.

This deliberately does **not** register from a static constructor firing during
``dlopen``.  Registration can fail — an ABI mismatch, or a solver name that is
already taken — and an exception thrown out of a static constructor while the
dynamic loader is running its init routines cannot be caught: it unwinds into C
loader frames and calls ``std::terminate()``, aborting the interpreter.  Running
registration from an explicitly-called entry point instead means every failure
comes back as a normal, catchable Python exception.

The lookup flow is:

.. code-block:: text

    Python: lightsim2grid.load_algorithm_plugin("path/to/plugin.so")
      └─ C++ _load_algorithm_plugin(path)        # a pybind11 function (try/catch)
           ├─ dlopen(path) / LoadLibrary(path)
           ├─ dlsym("ls2g_register_plugin")      # the LS2G_PLUGIN_ENTRY function
           └─ entry(&AlgorithmRegistry::instance(), errbuf, len)
                └─ AlgorithmRegistry::register_all(...)   # atomic; ABI-checked
                                                          # failure → return code → Python exception

    Python: grid.change_algorithm("MySolver")
      └─ AlgorithmSelector::change_algorithm("MySolver")
           └─ AlgorithmRegistry::instance().make("MySolver")
                └─ factory()  →  unique_ptr<BaseAlgo>

.. important::
    **Choose your solver name carefully: it is an identity, and it is persisted.**
    When a grid is saved (pickle or :ref:`binary-serialization`) it records the
    *name* of the solver it was using, and loading re-selects the solver
    registered under that name.

    Because of that, names are restricted to a small ASCII subset — registering
    anything else is refused, with an error explaining the rule::

        first character  : A-Z a-z _
        other characters : A-Z a-z 0-9 _ .
        length           : 1 to 64 characters

    (``lightsim2grid.lightsim2grid_cpp``'s C++ header exposes this as
    ``ls2g::is_valid_solver_name``.) Non-ASCII characters are excluded so that no
    plugin can register a name that *looks* identical to another solver's — say
    ``NR_SparseLU`` with a Cyrillic ``а`` — and masquerade as it in error
    messages, documentation and ``available_algorithm_names()`` listings while
    being a different registry key. Control characters are excluded because names
    are printed into error messages and logs, where a newline or an ESC byte
    would let a name inject content. The length bound keeps a name from
    overflowing the diagnostic buffer a plugin is given.

    Beyond the character set:

    - a grid saved while using your plugin can only be loaded where that plugin is
      loaded. Otherwise loading raises an error naming the missing solver; the
      reader can fall back to ``LSGrid.load_binary_without_algorithm(path)`` to get
      the grid data with the default solvers;
    - the name is the *only* thing stored — there is no version, checksum or
      identity of the plugin itself. If two plugins register the same name, a grid
      saved with one will silently be loaded with the other, even though they may
      be entirely different algorithms. lightsim2grid cannot detect this, so
      **prefix your solver names** with something specific to your project
      (``"acme_NR_fast"`` rather than ``"NR_fast"``) if your files travel between
      environments.

The ``AlgorithmRegistry`` C++ API (defined in ``AlgorithmRegistry.hpp``, installed to ``include/lightsim2grid/``):

.. code-block:: cpp

    class AlgorithmRegistry {
    public:
        // Meyers singleton — one instance per process.
        static AlgorithmRegistry& instance();

        // Register a single factory under a name (used for the built-ins).
        void register_solver(const std::string& name, Factory f);

        // Register a whole batch atomically: all names are added, or none is
        // and the registry is left unchanged. This is what a plugin goes
        // through. Throws (registry untouched) on an ABI mismatch or a name
        // that is already registered.
        void register_all(FactoryMap batch, Ls2gAbiTag caller_tag = ls2g_current_abi_tag());

        // Instantiate a solver by name.  Throws std::invalid_argument if
        // the name is unknown.
        std::unique_ptr<BaseAlgo> make(const std::string& name) const;

        bool is_registered(const std::string& name) const;

        // List of all registered names (built-in + plugins).
        std::vector<std::string> available_algorithm_names() const;
    };

    // Staging area handed to your registration function; add() one solver per
    // name. The LS2G_PLUGIN_ENTRY macro turns your function into the exported
    // ls2g_register_plugin entry point and commits the batch atomically.
    class PluginRegistrar {
    public:
        void add(const std::string& name, AlgorithmRegistry::Factory f);
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
stored in the public member ``IS_AC``, which LSGrid uses to route
:func:`LSGrid.change_algorithm` calls to the right slot (AC or DC).

Methods to override
~~~~~~~~~~~~~~~~~~~~

``compute_pf`` *(virtual — override to implement your algorithm; the base
implementation throws* ``std::runtime_error`` *if you don't)*

.. code-block:: cpp

    virtual bool compute_pf(
        const Eigen::Ref<const Eigen::SparseMatrix<cplx_type>>& Ybus,
        const Eigen::Ref<const CplxVect>& V,
        const Eigen::Ref<const CplxVect>& Sbus,
        const Eigen::Ref<const IntVect>& slack_ids,
        const Eigen::Ref<const RealVect>& slack_weights,
        const Eigen::Ref<const IntVect>& pv,
        const Eigen::Ref<const IntVect>& pq,
        int max_iter,
        real_type tol);

Solve the power-flow problem ``V·(Ybus·V)* = Sbus``.

+-------------------+------------------------------------------------------+
| Parameter         | Meaning                                              |
+===================+======================================================+
| ``Ybus``          | Complex bus admittance matrix (sparse, n×n).         |
+-------------------+------------------------------------------------------+
| ``V``             | Initial voltage phasors (read-only input -- every    |
|                   | parameter here is a ``const`` reference, so nothing  |
|                   | is updated in place). Store your solved result in    |
|                   | the ``V_`` protected member below instead; callers   |
|                   | retrieve it through ``get_V()``.                     |
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

``set_lsgrid`` *(virtual — override if you need grid data)*

.. code-block:: cpp

    virtual void set_lsgrid(const LSGrid* gridmodel);

Called by ``ChooseAlgorithm`` after the solver is activated (and again after
every ``change_algorithm`` call).  The default implementation stores the
pointer in the protected member ``lsgrid_ptr_``.  Override only if your
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

``get_theta_to_J_col_python`` / ``get_vm_to_J_col_python`` / ``get_q_to_J_col_python`` *(virtual — override only for Newton-Raphson algorithms)*

.. code-block:: cpp

    virtual IntVect get_theta_to_J_col_python() const;
    virtual IntVect get_vm_to_J_col_python()    const;
    virtual IntVect get_q_to_J_col_python()     const;

Expose the mapping between a bus index (solver id) and the position of its
unknown in the Jacobian.  Each returns a vector of size ``n_bus`` where entry
``bus_id`` is the Jacobian **column** holding that bus' unknown, or ``-1`` when
the bus owns no such unknown:

* ``get_theta_to_J_col_python`` — column of the bus' voltage-angle (``theta`` /
  ``ΔVa``) unknown;
* ``get_vm_to_J_col_python`` — column of the bus' voltage-magnitude (``vm`` /
  ``ΔVm``) unknown;
* ``get_q_to_J_col_python`` — column of the bus' reactive (``q``) unknown.

The base implementation throws ``std::runtime_error`` (these only make sense for
solvers that build a Jacobian).  For a Newton-Raphson solver built on
``NRSystem`` you can simply forward to its ``theta_to_J_col()`` /
``vm_to_J_col()`` / ``q_to_J_col()`` accessors.  They surface in Python as
:py:meth:`LSGrid.get_theta_to_J_col`, :py:meth:`LSGrid.get_vm_to_J_col` and
:py:meth:`LSGrid.get_q_to_J_col`.

Protected result members
~~~~~~~~~~~~~~~~~~~~~~~~~

Your ``compute_pf`` implementation must populate these fields so that
``LSGrid`` can extract flows and inject the results back into the grid
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

Place this in a single ``.cpp`` file.  Write a small registration function and
hand it to ``LS2G_PLUGIN_ENTRY``, which generates the exported
``ls2g_register_plugin`` entry point ``load_algorithm_plugin()`` calls.

.. code-block:: cpp

    // my_solver_plugin.cpp
    #include <AlgorithmRegistry.hpp>
    #include <powerflow_algorithm/BaseAlgo.hpp>

    class MySolver : public ls2g::BaseAlgo {
    public:
        MySolver() : ls2g::BaseAlgo(/*is_ac=*/true) {}

        bool compute_pf(
            const Eigen::Ref<const Eigen::SparseMatrix<ls2g::cplx_type>>& Ybus,
            const Eigen::Ref<const ls2g::CplxVect>& V,
            const Eigen::Ref<const ls2g::CplxVect>& Sbus,
            const Eigen::Ref<const ls2g::IntVect>& slack_ids,
            const Eigen::Ref<const ls2g::RealVect>& slack_weights,
            const Eigen::Ref<const ls2g::IntVect>& pv,
            const Eigen::Ref<const ls2g::IntVect>& pq,
            int max_iter,
            ls2g::real_type tol) override
        {
            // ... your algorithm here, starting from the initial guess `V` ...
            ls2g::CplxVect V_solved = V;  // replace with your actual solved voltages

            // populate result fields when done:
            V_      = V_solved;
            Va_     = V_solved.array().arg();
            Vm_     = V_solved.array().abs();
            n_      = static_cast<int>(V_solved.size());
            nr_iter_= 1;
            err_    = ls2g::ErrorType::NoError;
            return true;
        }
    };

    // Registration: load_algorithm_plugin() calls this (via the generated
    // ls2g_register_plugin entry point) after loading the library. Add one
    // solver per name; register a whole batch and they commit atomically.
    static void register_plugin(ls2g::PluginRegistrar& reg) {
        reg.add("MySolver",
                []{ return std::unique_ptr<ls2g::BaseAlgo>(new MySolver()); });
    }
    LS2G_PLUGIN_ENTRY(register_plugin)

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
        if(NOT EXISTS "${LIGHTSIM2GRID_SRC}/AlgorithmRegistry.hpp")
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

    # Match lightsim2grid_core's -march=native / -O3. This is required, not
    # just a performance nicety: lightsim2grid_core and my_solver are two
    # *separate* shared libraries, and -march=native changes which SIMD ISA
    # Eigen sees enabled, which changes EIGEN_MAX_ALIGN_BYTES -- an Eigen
    # object (RealVect, CplxVect, ...) allocated under one alignment
    # assumption and freed under another (possible the moment such an object
    # crosses the BaseAlgo interface) silently corrupts the heap. See the
    # "Matching build flags" section below.
    if(lightsim2grid_core_FOUND)
        # Installed package: it exports the flags it was actually built
        # with -- no guessing needed.
        set(_ls2g_march_native "${lightsim2grid_core_MARCH_NATIVE}")
        set(_ls2g_o3_optim     "${lightsim2grid_core_O3_OPTIM}")
    else()
        # Source tree: no installed package to query -- fall back to the
        # same env vars lightsim2grid itself reads.
        if("$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "1" OR "$ENV{__COMPILE_MARCHNATIVE}" STREQUAL "True")
            set(_ls2g_march_native ON)
        else()
            set(_ls2g_march_native OFF)
        endif()
        # -O3 is ON by default in lightsim2grid_core: only __O3_OPTIM=0/False
        # turns it off, an unset env var means ON.
        if("$ENV{__O3_OPTIM}" STREQUAL "0" OR "$ENV{__O3_OPTIM}" STREQUAL "False")
            set(_ls2g_o3_optim OFF)
        else()
            set(_ls2g_o3_optim ON)
        endif()
    endif()
    if(NOT MSVC)
        if(_ls2g_march_native)
            message(STATUS "my_solver: lightsim2grid_core was built with -march=native -- matching it")
            target_compile_options(my_solver PRIVATE -march=native)
        endif()
        if(_ls2g_o3_optim)
            message(STATUS "my_solver: lightsim2grid_core was built with -O3 -- matching it")
            target_compile_options(my_solver PRIVATE -O3)
        endif()
    endif()

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
    from lightsim2grid.network import init_from_pandapower
    from lightsim2grid.lightsim2grid_cpp import LSGrid

    # 1. Load the plugin — this loads the library, calls its
    #    ls2g_register_plugin entry point, and registers "MySolver" into the
    #    AlgorithmRegistry singleton. A registration failure (ABI mismatch,
    #    duplicate name, ...) raises a catchable exception here.
    lightsim2grid.load_algorithm_plugin("build/libmy_solver.so")

    # 2. Confirm the solver is available (the registry is process-wide, so this
    #    works even on a throwaway, empty grid).
    print(LSGrid().available_algorithm_names())   # [..., "MySolver", ...]

    # 3. Build a grid and activate the plugin solver on it.
    pp_net = ...  # any pandapower grid
    gm = init_from_pandapower(pp_net)
    gm.change_algorithm("MySolver")

    # 4. Run a powerflow — lightsim2grid now delegates to MySolver.compute_pf().
    Vinit = ...  # complex initial voltage guess, size gm.total_bus()
    V = gm.ac_pf(Vinit, 20, 1e-8)


Python API reference
---------------------

.. py:function:: lightsim2grid.load_algorithm_plugin(path: str) -> None

    Load a shared library containing one or more lightsim2grid solver / algorithm
    plugins.

    The library must export a ``ls2g_register_plugin`` entry point, which the
    ``LS2G_PLUGIN_ENTRY`` macro generates from your registration function (see
    the example above).  The library is loaded, the entry point is called, and
    the solver name(s) it declares are registered into the ``AlgorithmRegistry``
    singleton.

    After this call the new solver is usable via
    ``grid.change_algorithm("MyAlgoName")`` and will appear in
    ``grid.available_algorithm_names()``.

    Because registration runs from an explicitly-called entry point rather than
    a static constructor during ``dlopen``, every failure is raised here as a
    normal, catchable exception instead of aborting the interpreter.  Loading a
    plugin whose name collides with an already-registered one — including
    loading the same plugin twice — is rejected the same way, and the registry
    is left unchanged (registration is atomic).

    :param path: Absolute or relative path to the ``.so`` / ``.dll`` file.
    :raises RuntimeError: If the library cannot be loaded (missing file,
        unresolved symbols), does not export the ``ls2g_register_plugin`` entry
        point, or its registration is refused (ABI mismatch, duplicate name).

.. py:method:: LSGrid.change_algorithm(name: str) -> None

    Select the active solver by name.  The name must be one of the
    strings returned by :py:meth:`LSGrid.available_algorithm_names`.

    For built-in solvers the enum overload is also available::

        gm.change_algorithm(AlgorithmType.NR_KLU)

    :param name: Registered solver name (case-sensitive).
    :raises RuntimeError: If *name* is not registered.

.. py:method:: LSGrid.available_algorithm_names() -> list[str]

    Return all solver names currently registered, including any that were
    added via :func:`~lightsim2grid.load_algorithm_plugin`.

    .. code-block:: python

        >>> gm.available_algorithm_names()
        ['NR_SparseLU', 'GaussSeidel', 'NRSing_KLU',
         # and all other lightsim2grid installed "solver"
         'MySolver']

.. py:method:: LSGrid.get_theta_to_J_col() -> numpy.ndarray
.. py:method:: LSGrid.get_vm_to_J_col() -> numpy.ndarray
.. py:method:: LSGrid.get_q_to_J_col() -> numpy.ndarray

    Return the mapping from a bus index (solver id) to the column of the
    Jacobian holding that bus' unknown.  Each call returns an integer array of
    size ``n_bus`` whose entry ``bus_id`` is the Jacobian column of the bus'
    ``theta`` (voltage angle), ``vm`` (voltage magnitude) or ``q`` (reactive)
    unknown respectively, or ``-1`` when the bus owns no such unknown.

    These are only meaningful for Newton-Raphson solvers (those that build a
    Jacobian, see :func:`lightsim2grid.network.LSGrid.get_J_solver`).  They mirror the ``get_J()`` accessor and let
    you locate a given bus' rows/columns inside the augmented Jacobian.

    :raises RuntimeError: If the active solver does not build a Jacobian.

.. note::

    :attr:`AlgorithmType.Custom` is the enum value assigned to any solver
    loaded via a plugin.  ``gm.get_algo_type()`` returns
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

.. _solver_plugin_march_native:

Matching build flags (``-march=native`` / ``-O3``)
----------------------------------------------------

``lightsim2grid_core`` supports two build flags, set as environment
variables *before* building lightsim2grid (see ``env_compile_all.sh`` at the
repository root):

* ``__COMPILE_MARCHNATIVE=1`` — adds ``-march=native`` (opt-in; off by default)
* ``__O3_OPTIM=0`` — removes ``-O3`` (opt-out; ``-O3`` is on by default)

A plugin **must** be compiled with the identical ``-march=native`` setting as
the ``lightsim2grid_core`` it links against, and this is a correctness
requirement, not a performance one. ``lightsim2grid_core`` and a plugin
module are two *separate* shared libraries. ``-march=native`` changes which
SIMD instruction sets Eigen sees enabled (``__AVX__``, ``__AVX2__``, ...),
which changes ``EIGEN_MAX_ALIGN_BYTES`` and thus how Eigen aligns / allocates
/ frees its dynamic-size matrices (``RealVect``, ``CplxVect``, ...). If the
two libraries disagree, an Eigen object allocated under one alignment
assumption and freed under another silently corrupts the heap — surfacing
later, and far away, as ``free(): invalid size`` or ``double free or
corruption``. This is a real, previously-hit bug (see the CHANGELOG entry
for the ``lightsim2grid_core`` / ``lightsim2grid_cpp`` split) — not a
hypothetical one. (``-O3`` alone does not affect ABI/alignment, so a mismatch
there is only a lost optimization, not a correctness risk — but the two
builds should still track each other to avoid surprising perf differences.)

All three example plugins in this repository (``examples/dist_slack_algorithm/``,
``examples/external_algorithm/``, ``examples/lm_algorithm/``) handle this
automatically, and log an
``-- my_solver: lightsim2grid_core was built with -march=native -- matching
it`` status message when they do:

* When linking against an **installed** ``lightsim2grid_core`` (``pip
  install lightsim2grid``, ``find_package(lightsim2grid_core CONFIG)``), the
  installed CMake package exports ``lightsim2grid_core_MARCH_NATIVE`` /
  ``lightsim2grid_core_O3_OPTIM`` — the flags it was *actually* built with —
  so the plugin build queries them directly, no guessing or separate Python
  call needed.
* When falling back to a **source-tree** build (no installed package), there
  is nothing to query, so the plugin build instead reads the same
  ``__COMPILE_MARCHNATIVE`` / ``__O3_OPTIM`` env vars lightsim2grid itself
  reads — set them identically for both builds.

The shared logic lives in
``examples/cmake/MatchLightsim2gridBuildFlags.cmake`` as the
``ls2g_match_core_build_flags(<target>)`` function; ``include()`` it and call
it on your plugin's target right after ``target_link_libraries(<target>
PRIVATE lightsim2grid::core)`` if you are building inside this repository.
The "Writing a solver plugin" CMakeLists.txt template above inlines the same
logic for third-party projects that do not have access to that file.

You can also check what a given lightsim2grid install was built with
directly from Python::

    >>> import lightsim2grid
    >>> lightsim2grid.compilation_options.compiled_march_native
    True
    >>> lightsim2grid.compilation_options.compiled_o3_optim
    True

As defense-in-depth on top of the CMake-level matching above, this is also
checked **at runtime**. Every plugin solver registration (``load_algorithm_plugin()``
/ ``LS2G_PLUGIN_ENTRY``) computes an "ABI tag" — ``EIGEN_MAX_ALIGN_BYTES``, the
resolved ``EIGEN_VECTORIZE_SSE``/``AVX``/``AVX2``/``AVX512``/``FMA``/``NEON``
flags, and the full Eigen version (major.minor.patch) — in the plugin's own
translation unit, and compares it against the tag ``lightsim2grid_core`` was
actually compiled with. If they differ, registration is refused with an error
naming exactly what mismatched, instead of silently registering a plugin that
would corrupt the heap the first time one of its Eigen objects crosses the
``BaseAlgo`` interface. The same comparison runs once between
``lightsim2grid_core`` and ``lightsim2grid_cpp`` themselves at import time —
this is the check that would have caught the ``lightsim2grid_core`` /
``lightsim2grid_cpp`` split bug directly, as a clean Python ``ImportError``
instead of a crash. This does not replace matching the build flags above (a
rejected plugin is still a plugin you can't use), but it turns a silent,
far-away heap corruption into an immediate, actionable error at load time.
See ``src/core/Ls2gAbiTag.hpp``.

The Eigen version is part of the tag, compared exactly (not just the major
version): Eigen is header-only, so its semver promise
(https://libeigen.gitlab.io/news/eigen_5.0.0_released/) is about API
compatibility for code *using* Eigen, not about the internal aligned-malloc
offset arithmetic this check actually cares about staying identical across
two independently-compiled binaries — that isn't something Eigen tests or
promises release-to-release, since there is no shared-object ABI to keep
stable in the first place. Practically, this means **a plugin should be
compiled against the same Eigen headers lightsim2grid itself uses**, not a
separately-installed system Eigen:

* Linking against an **installed** ``lightsim2grid_core`` (the normal case,
  ``find_package(lightsim2grid_core CONFIG)``): the package bundles the exact
  Eigen it was built with alongside its own headers, and
  ``lightsim2grid_core_INCLUDE_DIRS`` already points at it — as long as your
  plugin doesn't add a different Eigen include path ahead of that one (e.g. a
  system ``/usr/include/eigen3``), it picks up the matching Eigen
  automatically and this check is a non-issue.
* Building against the **source tree**: use the ``eigen/`` submodule checked
  out in this repository (the same one the examples' ``Eigen3_INCLUDE``
  fallback points to), not a separately-installed copy.

If you do need a different Eigen version for your plugin for some other
reason, that's a deliberate, advanced choice this check exists specifically
to flag — the error message tells you exactly which field(s) mismatched.


Worked example (``examples/external_algorithm/``)
-------------------------------------------------

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
    Registered solvers: ['DC_KLU', 'DC_SparseLU', 'DummyExternal', 'FDPF_BX_KLU', ...]
    change_algorithm('DummyExternal') OK -- solver type is Custom as expected.
    All checks passed.

The automated regression test is in
``lightsim2grid/tests/test_solver_registry.py``, class
``TestPluginLoading``.  It is skipped automatically when the example
plugin has not been built, and passes once the ``.so`` is present.

Other example plugins
-----------------------

Two more complete, less minimal plugins live in the repository, both following the same build /
``load_algorithm_plugin`` / ``change_algorithm`` pattern as ``external_algorithm`` above:

* ``examples/dist_slack_algorithm/`` registers ``NRDistSlack_SparseLU`` / ``NRDistSlack_KLU``: a
  distributed-slack Newton-Raphson variant (``NRAlgoDistSlack``), with both a ``SparseLU`` and a
  ``KLU`` linear-solver backend.
* ``examples/lm_algorithm/`` registers ``NR_LM_SparseLU`` / ``NR_LM_KLU``: a
  Levenberg-Marquardt-damped Newton-Raphson variant (``LMNRAlgo``), also with ``SparseLU`` and
  ``KLU`` backends. Both names are registered atomically -- if either is already taken, neither is
  added and ``load_algorithm_plugin`` raises instead of half-registering.


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
