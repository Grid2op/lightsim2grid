Change Log
===========

[TODO]
--------
- ``SubstationContainer::sub_vn_kv_`` is dead state: its only writers (``init_sub()``
  and the two-argument constructor) are called from nowhere, and nothing reads it
  back, so it is empty on every grid every loader produces and in every binary file
  saved so far. It is also redundant -- a substation's nominal voltage is the common
  vn_kv of its buses, i.e. ``sub_vn_kv_[I] == bus_vn_kv_[I]`` for ``I < n_sub``, which
  ``check_valid()`` now enforces when it is present. Decide: either drop it from
  ``StateRes`` (needs a ``BINARY_FORMAT_VERSION`` bump) or have ``init_bus()`` fill it
  and make it mandatory (changes what newly-saved files contain, and can only become
  mandatory once no already-saved file needs to load). See the long note on the member
  itself in ``SubstationContainer.hpp``.
- [refacto] have a structure in cpp for the buses
- [refacto] have the id_grid_to_solver and id_solver_to_grid etc. directly in the solver and NOT in the gridmodel.
- [refacto] put some method in the DataGeneric as well as some attribute (_status for example)
- support 3w trafo (as modeled in pandapower)
- improve speed by not performing internal checks 
  (keep check for boundaries and all for python API instead) [see `TODO DEBUG MODE` in c++ code]
- improve speed
- code parrallelism directly in the `Computer` and `SecurityAnalysisCPP` classes
- use the "multi slack hack" (see issue #50) for SecurityAnalysis or Computer for example
- code `helm` powerflow method
- interface with gridpack (to enforce q limits for example)
- maybe have a look at suitesparse "sliplu" tools ?
- pybind11's Eigen support has no ``type_caster`` for ``Eigen::Ref<const
  Eigen::SparseMatrix<...>>`` (only for a concrete, owning ``Eigen::SparseMatrix``):
  its generic sparse caster default-constructs a ``Type value;`` member internally,
  which is impossible for ``Eigen::Ref`` (dense or sparse -- a ``Ref`` must be bound
  to something at construction). pybind11 works around this for *dense* ``Eigen::Ref``
  with a hand-written caster (delayed construction via ``std::unique_ptr``), but ships
  no equivalent for the sparse case. Because of this, every pybind11-exported function
  taking/returning a sparse matrix (``BaseAlgo::compute_pf_with_input_validation``'s
  ``Ybus``/``Bbus``, ``BaseFDPFAlgo::debug_get_Bp_python``/``debug_get_Bpp_python``,
  ``get_J_python``, etc.) uses a bare ``Eigen::SparseMatrix<...>`` instead of the
  ``EigenRefConstCplxSpMat``/``EigenRefConstRealSpMat`` aliases used everywhere else in
  the C++-internal API, which costs one matrix copy per python call. Eigen's
  ``SparseRef.h`` shows the ``Ref`` *can* bind zero-copy onto an
  ``Eigen::Map<SparseMatrix<...>>`` built straight from a scipy CSR/CSC matrix's
  ``data``/``indices``/``indptr`` buffers (our aliases use ``Options=0``, so the
  ``StandardCompressedFormat``-triggered copy path in ``Ref<const SparseMatrix>``'s
  constructor never fires). A custom pybind11 ``type_caster`` specialization for
  ``Eigen::Ref<const Eigen::SparseMatrix<Scalar,...>>`` (mirroring pybind11's own
  dense-``Ref`` ``eigen_map_caster`` pattern: hold the three numpy buffers + a
  ``unique_ptr<Map>``/``unique_ptr<Ref>`` built in ``load()``) would remove that copy
  on the load (python -> C++) direction, which is what matters for ``compute_pf``'s
  hot path. The return (C++ -> python) direction would still copy either way -- python
  needs its own owned object there regardless.
- `GeneratorContainer::is_pseudo_off` (used when ``turnedoff_no_pv`` /
  ``LightSimBackend(turned_off_pv=False)`` is active) is not implemented as
  originally intended: it currently only checks ``target_p == 0`` to decide a
  generator is "pseudo off" (and so should not pin its bus voltage). It was meant
  to emulate PowSyBl OpenLoadFlow's own rule
  (``generatorsWithZeroMwTargetAreNotStarted``,
  ``AbstractLfGenerator.checkIfGeneratorStartedForVoltageControl``), which
  additionally requires ``min_p > 0`` (a generator with ``min_p <= 0`` is allowed
  to legitimately sit at 0 MW and stays a normal voltage-controlling PV bus even
  at ``p == 0``). ``GeneratorContainer`` does not currently store ``min_p`` at
  all. There should be a convention (and a stored ``min_p`` per generator) so
  ``is_pseudo_off`` can check ``target_p == 0 AND min_p > 0``, matching OLF. For
  now this OLF-specific behavior is instead reproduced Python-side, on the
  pypowsybl-loading path only, by ``lightsim2grid.network.from_pypowsybl.
  bake_outer_loops`` (see ``_bake_generator_not_started`` in ``_olf_bake.py``),
  which rewrites ``voltage_regulator_on=False`` in the IIDM network for such
  generators before conversion, rather than relying on this C++ mechanism.
- ``fuse_zero_impedance_branches`` (``from_pypowsybl.init()``) deactivates a fused
  branch and leaves its "loser" bus with no elements, so the *reduced* network
  lightsim2grid actually solves no longer knows that bus's voltage or that
  branch's flow. This is invisible to the C++ solve itself (the reduced problem
  is solved correctly), but it means ``LSGrid.get_V()``/``get_lines()`` alone
  cannot answer "what is bus X's voltage" / "what flows through branch Y" for a
  fused-away bus or branch -- only ``LightsimResultNetwork``
  (``_result_network.py``) currently reconstructs these (fused-bus voltage via
  the new ``LSGrid._bus_fusion_rep``; a leaf-endpoint branch's flow via
  Kirchhoff's current law, see ``_reconstruct_fused_branches``), so anyone
  reading the grid directly through the C++ API (or another Python wrapper)
  still sees ``v_mag=0`` / 0 flow. What's implemented is fine for now, but at
  some point this should be pushed down into the C++ core itself (eg have
  ``LSGrid`` report a fused-away bus's voltage as its representative's
  directly, and/or not deactivate a fused branch attached to only one other
  fused branch) so every consumer benefits, not just callers that go through
  ``LightsimResultNetwork``.

TODO: speed directly update the pv, pq, Sbus and Ybus part when updating the elements
      (less error prone and faster to recompute). Then what is passed to the solver 
      is "only" a "mask" of these data when null
TODO: https://github.com/haranjackson/NewtonKrylov for another type of algorithm ?
TODO: in ContingencyAnalysisCpp: add back the `if(!ac_solver_used)` inside the  `remove_from_Ybus`
      in order to perform the "invertibility" check
TODO: in `main.cpp` check the returned policy of pybind11 and also the `py::call_guard<py::gil_scoped_release>()` stuff
TODO: a cpp class that is able to compute (DC powerflow) ContingencyAnalysis and TimeSeries using PTDF and LODF
TODO: integration test with pandapower (see `pandapower/contingency/contingency.py` and import `lightsim2grid_installed` and check it's True)
TODO: speed: `BaseBatchSolverSynch::compute_amps_flows` / `compute_active_power_flows` build
      `Efrom`/`Eto` via `CplxVect(_voltages.col(...))` for every line/trafo element, once per
      call -- a full nb_steps-sized copy each time. `_voltages` is RowMajor, so a column isn't
      contiguous and can't bind to `Eigen::Ref`, and the surrounding ternary also needs a
      common type with `CplxVect::Zero(nb_steps)` for the open-side case. Removing the copy
      would need a control-flow restructure (separate open-side / closed-side code paths),
      not just a reference-type change.
TODO: building ``examples/dist_slack_algorithm`` against a NICSLU/CKTSO-enabled
      ``lightsim2grid_core`` fails: ``NICSLUSolver.hpp`` (pulled in transitively via
      ``Solvers.hpp``) needs NICSLU's own SDK header ``nicslu_cpp.inl`` (and, presumably,
      CKTSO's own headers for ``CKTSOSolver.hpp``), which are a private implementation
      detail not exposed by an installed ``lightsim2grid_core`` CMake package. The plugin's
      ``CMakeLists.txt`` only special-cases SuiteSparse/KLU's headers this way (see the
      comment there), not NICSLU/CKTSO's. Found while verifying the ``-march=native``
      matching fix below against a full ``env_compile_all.sh`` build; not fixed.
TODO: add a CI job that builds one of the example C++ algorithm plugins
      (``examples/dist_slack_algorithm`` or ``examples/external_algorithm``) against a
      ``lightsim2grid_core`` built with ``__COMPILE_MARCHNATIVE=1``, and actually loads +
      runs it (not just compiles it), to catch a regression of the ``-march=native``
      matching mechanism (``lightsim2grid_core_MARCH_NATIVE`` export +
      ``examples/cmake/MatchLightsim2gridBuildFlags.cmake``, see ``docs/solver_plugin.rst``).
      ``.github/workflows/main.yml``'s ``test_march_native`` job and
      ``test_plugin_against_installed`` job each cover half of this (the former builds
      lightsim2grid itself with the flag but does not touch a plugin; the latter builds a
      plugin but never with the flag) -- neither exercises the combination.
TODO: Levenberg-Marquardt damping (a.k.a. Tikhonov-regularized Newton) : adding small decreasing
      lambba coefficients to the diagonal of J to improve its conditionning.
TODO: a "combine mode" axis for ``ScenarioSweepCPP`` choosing between the current
      row-aligned / zipped semantics (row `i`'s injection paired with row `i`'s
      contingency) and a cartesian one ("every registered contingency x every
      injection profile").

[1.0.1] 2026-xx-yy
--------------------
- [IMPROVED] the last throw sites inside the ``_apply_and_track_buses`` bracket are gone, which
  takes the branch ``fillYbus`` back to **exactly** what it cost before #188 landed: 7,222,254
  instructions, the same figure to the digit. Together with the previous entry that is
  12,085,107 -> 7,222,254 on that function (**-40.2%**), ``pre_process_solver``
  33,503,851 -> 28,640,988, and **-4,862,863 on a whole case9241pegase rebuild solve (-0.56%)**.
  Validating the bus id before the bracket removed the biggest unwind edge but not all of them:
  ``deactivate_no_bus_tracking`` / ``reactivate_no_bus_tracking`` /
  ``change_bus_no_bus_tracking`` each re-checked the element id, and so did
  ``_generic_deactivate`` / ``_generic_reactivate`` / ``_generic_change_bus`` underneath them --
  three layers of the same check, the outermost of which is the user-facing one. The inner two are
  now ``_check_in_range_internal``, the debug-only form: every one of the 18 bracket call sites is
  reached through a public entry point that has already raised for a bad id, and a throw from
  inside the bracket is the thing this whole layer exists to avoid.
  **No user-facing check was lost**, verified from a release build rather than by reading: a
  release-built driver still raises for ``change_bus_load(0, -1)``, ``change_bus_load(0, nb+1000)``,
  the generator equivalents, ``change_bus_load(999999, 0)`` and ``deactivate_load(999999)`` -- six
  out of six. The Debug library carries all ten internal messages, the Release library four (the
  user-facing ones).
  Two stale comments went with it: ``GeneratorContainer::_change_bus`` and
  ``SvcContainer::_change_bus`` both claimed their IndexError came from ``_generic_change_bus``
  "which the caller runs *after* this function". It has not for some time -- 
  ``change_bus_no_bus_tracking`` raises first.
  227 test cases / 590,769 assertions pass in the C++17, C++14 and Debug builds -- including #188's
  refused-``change_bus`` sweep, which is what makes this safe to do -- and results are
  bit-identical across 16 AC and 8 DC configurations on case118 and the three PEGASE cases.
- [IMPROVED] a ``change_bus`` the grid is going to refuse is now rejected **before** the per-bus
  element counts are touched, instead of being compensated for afterwards. -3,899,955 instructions
  on a case9241pegase rebuild solve (**-0.45% of everything**), of which the whole of the branch
  ``fillYbus``' share: 12,085,107 -> 8,185,194, **-32.3%**.
  ``_apply_and_track_buses`` brackets a mutation between "take this element's contribution away"
  and "put it back". The bus-id validation lived inside that bracket, in ``_generic_change_bus``,
  so a refused call was rejected with the contribution already removed -- which is why it needed a
  ``try`` / ``catch(...)`` restoring it on the way out. That catch was not free: an unwind edge
  through ``GenericContainer.hpp`` made GCC keep every ``std::vector<bool>`` access in ``fillYbus``
  live across it, in a function that never calls any of this. Bisecting the nine commits of #188
  against that single number found it exactly -- eight of them at 7,463,013 and the ninth at
  8,859,264 -- and deleting only the catch, keeping everything else in that commit, restored
  7,463,013 to the instruction.
  The new ``GenericContainer::_check_new_bus_id`` is called by the four ``change_bus`` entry points
  before they enter the bracket, so a refused call never touches the counts at all rather than
  touching them and undoing it. It is always active: the id comes from the caller. An exception
  from deeper inside a mutation can still leave the counts short, and that is deliberate -- such a
  grid must be rebuilt and its caches dropped, not carried on with.
  #188's own coverage is what caught the first attempt at this: it found the HVDC and two-sided
  ``change_bus`` paths where the check had landed inside the lambda instead of before the bracket
  (36 failing assertions). All 227 test cases / 590,769 assertions pass in the C++17, C++14 and
  Debug builds, and results are bit-identical on case118 and the three PEGASE cases across 16 AC
  and 8 DC configurations.
- [IMPROVED] the same debug-only treatment for the "connected to a disconnected bus" checks on the
  **results path**, where it is worth considerably more than on the assembly one:
  ``LSGrid::compute_results`` drops **-1,591,221 instructions (-7.65%)** and **-7.4% of its wall
  time** (min of 11 batches of 2000 calls: 954 -> 884 microseconds, four runs, the two ranges not
  overlapping). Nine blocks, in ``TwoSidesContainer_rxh_A::compute_results_tsc_rxha_no_amps`` (four,
  inside a loop over all 16049 branches of case9241pegase),
  ``GenericContainer::v_kv_theta_from_vpu`` (two, over every load, static gen, storage, shunt,
  generator and SVC), ``ShuntContainer::_compute_results`` (two) and
  ``LSGrid::_get_results_back_to_orig_nodes`` (one, per bus). Each was an ``std::ostringstream``
  built inline in the loop.
  Same argument as the assembly path: these run only from ``LSGrid::process_results``, on a bus the
  container just read out of its own ``bus_id_`` for an element already established as connected,
  and on the ``id_me_to_solver`` this very solve built. Not gated, and deliberately: the guard in
  ``HvdcLineContainer::compute_results`` that rejects a half-open droop line, which is a state
  invariant rather than an id check and is the documented alternative to indexing with the open
  side's -1.
  Verified from the binaries: the Release and C++14 libraries carry none of the nine messages, the
  Debug library carries them all. Results bit-identical on case118 and the three PEGASE cases, 16
  AC configurations and 8 DC ones, every element result at 17 digits.
- [IMPROVED] the "connected to a disconnected bus" checks in the **matrix and injection assembly**
  are now debug-only, the same treatment ``_check_in_range`` got and for the same reason: the ids
  they test are ones this library just produced. ``fillYbus`` / ``fillSbus`` / ``fillBdc`` /
  ``hack_Sbus_for_dc_phase_shifter`` are reached only from ``LSGrid::_build_into_cache``; the bus
  comes out of the container's own ``bus_id_`` for an element the loop has already established is
  connected, and the solver id out of the ``id_grid_to_solver`` LSGrid built two steps earlier. A
  failure there is not a caller error but an inconsistency between an element's status and its bus
  -- which ``check_grid()`` validates, and which no public method can produce on its own. Each was
  an ``std::ostringstream`` built inline in the assembly loop, so they also cost the surrounding
  code its registers. 18 blocks across 7 containers.
  Worth -240,759 instructions on the branch ``fillYbus`` (**-3.2%**), -394,671 on ``LSGrid::fillYbus``,
  -359,922 on a whole case9241pegase AC rebuild solve and -169,974 on a DC one. Too small to
  separate from run-to-run noise in wall clock; the instruction counts are exact and reproducible.
  **Nothing a user can reach lost a check**, and the assertion builds keep every one of them: the
  release libraries carry none of the message strings, the Debug library carries them all. Results
  bit-identical on case118 and the three PEGASE cases -- 16 AC configurations (NR / NRSing / FDPF,
  SparseLU and KLU) and 8 DC ones -- across every element result at 17 digits.
- [IMPROVED] ``TwoSidesContainer_rxh_A``'s branch flow results are computed on real and imaginary
  parts instead of through ``std::complex``. -963k instructions of the 21.7M
  ``LSGrid::compute_results`` spends on a case9241pegase solve (-4.4%) and -2.9% on the function's
  wall time (min of 11 batches of 2000 calls, four runs, the two ranges not overlapping). Each of
  the grid's 16049 branches costs six complex-times-complex products -- ``y11.Ehv``, ``y12.Elv``,
  ``y22.Elv``, ``y21.Ehv``, then ``Ehv.conj(I_hvlv)`` and ``Elv.conj(I_lvhv)`` -- and
  ``std::complex`` follows every one of them with a branch that re-derives the result if it came
  out NaN. That was a third (4.24M) of everything the branch flow loop did. Results are
  bit-identical: the products are grouped exactly as ``std::complex`` groups them, ``conj`` is an
  exact sign flip, and the recovery path only fires on a non-finite product. Verified on every
  element result (P/Q/V/theta/amps of lines, trafos, loads, gens, shunts, sgens, storages and
  hvdc, at 17 digits) over case118 and the three PEGASE cases under NR / NRSing / FDPF with
  SparseLU and KLU.
- [IMPROVED] ``GenericContainer::_get_amps`` computes the current in one pass instead of four.
  **-14% on ``LSGrid::compute_results``' wall time** (min of 11 batches of 2000 calls: 1164 -> 1001
  microseconds, four runs) and -1.55M instructions of 23.3M (-6.7%). It read
  ``a = sqrt(p^2 + q^2) / (sqrt(3) . v)`` as four separate expressions -- the sum of squares into a
  vector, a square root over that vector, a full copy of v, then a scan of the copy replacing the
  zeros -- which is two full-length heap allocations and three extra passes over memory, on every
  one of the four calls a solve makes (both ends of the powerlines and of the transformers). The
  guard that stops a disconnected element's zero voltage dividing by zero is now a ternary inside
  the one loop. Same arithmetic in the same order, so results are bit-identical -- verified on
  every element result (P/Q/V/theta/amps of lines, trafos, loads, gens, shunts, sgens, storages and
  hvdc, at 17 digits) over case118 and the three PEGASE cases under NR / NRSing / FDPF with
  SparseLU and KLU. The wall-clock gain is much larger than the instruction count suggests, which
  is the allocations: they cost time, not instructions.
- [IMPROVED] ``GenericContainer::v_kv_from_vpu`` and ``v_deg_from_va`` are one function,
  ``v_kv_theta_from_vpu``. -862k instructions of the 24.1M ``LSGrid::compute_results`` spends on a
  case9241pegase solve (-3.6%), and -2.9% on the function's wall time (min of 11 batches of 2000
  calls, four runs, no overlap between the two). They were called back to back on the same elements
  by ``OneSideContainer::compute_results`` -- for loads, static gens, storages, shunts, generators
  and SVCs -- and everything before the last line of each was the same work: read the element's
  bus, map it through ``id_grid_to_solver``, check that neither is the deactivated bus. Only the
  final assignment differed (Vm times the bus' nominal kV, against Va times 180/pi). The walk, the
  two gathers and the two checks now happen once and both results are written from them. Results
  are bit-identical, verified on every element result (P/Q/V/theta/amps of lines, trafos, loads,
  gens, shunts, sgens, storages and hvdc, at 17 digits) over case118 and the three PEGASE cases
  under NR / NRSing / FDPF with SparseLU and KLU.
- [IMPROVED] ``ShuntContainer::_compute_results`` computes the shunt's power on real and imaginary
  parts instead of through ``std::complex``. -264k instructions of the 24.4M ``LSGrid::compute_results``
  spends on a case9241pegase solve (-1.1% of the post-processing). The expression is
  ``s = E * conj(y * E)`` with ``y = -(p + i.q) / sn_mva``: two complex-times-complex products, each
  followed by ``std::complex``'s NaN-recovery branch, on each of the grid's 7327 shunts on every
  single solve. Written out on parts the branches are gone and the results are bit-identical -- the
  recovery path only fires on a non-finite product, and ``(-1 * x) / s`` and ``-(x / s)`` agree
  exactly in IEEE 754. Verified against the previous build on every element result
  (P/Q/V/theta/amps of lines, trafos, loads, gens, shunts, sgens, storages and hvdc, at 17 digits)
  over case118 and the three PEGASE cases under NR / NRSing / FDPF with SparseLU and KLU.
- [IMPROVED] ``BaseDCAlgo::remove_slack_buses`` writes the slack-free DC matrix' compressed arrays
  itself instead of collecting the surviving coefficients into a vector of ``Eigen::Triplet`` and
  handing that to ``setFromTriplets``. **3.09x on the function** (11.4M -> 3.7M instructions for
  three rebuilds of case9241pegase) and **-4.9% on a whole DC solve** (dc_pf 158.5M -> 150.8M);
  wall clock per rebuild solve 12.97 -> 12.53 ms on case9241pegase, 3.39 -> 3.31 on case2869pegase,
  1.53 -> 1.45 on case1354pegase.
  Nothing here ever needed sorting or merging. ``mat_bus_id_`` is a monotone compaction --
  ``fill_mat_bus_id`` hands out consecutive ids in bus order, skipping the slack buses -- and the
  inner iterator walks each column of the source in increasing row order, so the coefficients that
  survive the row/column deletion come out already in the order a compressed column-major matrix
  stores them, one per coefficient. ``setFromTriplets`` cannot know that: it re-derived that exact
  order the hard way, bucketing by row and transposing back into columns, both through a temporary
  ``SparseMatrix`` carrying a double per entry, for 5% of a DC solve. The matrix is now built in
  place in two passes over the source, one to size each column (straight into the outer array,
  which a prefix sum turns into the column starts) and one to fill it. The triplet vector -- half a
  megabyte on case9241pegase -- is gone with it, and ``res_mat`` is resized rather than reassigned,
  so a rebuild at constant size re-uses what the previous one allocated.
  The sortedness this rests on was measured, not assumed: an assertion over every consecutive pair
  found the old triplet list strictly ordered column-major with no duplicates in all 136 builds of
  the C++ test suite and on the four benchmark grids under both DC solvers. The new build was then
  cross-checked against the old one coefficient by coefficient -- same nnz, same outer array, same
  inner array, same values bit for bit -- over the same 136 builds and the same grids. What remains
  in the code is the cheaper permanent form of that check, in debug builds: the matrix is
  compressed, each column is filled exactly to the size the first pass computed, and its row
  indices are strictly increasing (which is also what would catch ``mat_bus_id_`` losing the
  monotonicity the whole function is built on).
  Results are bit-identical -- voltages, line flows and generator P/Q, to 17 digits -- on case118
  and the three PEGASE cases under ``DC_KLU`` and ``DC_SparseLU``, and the AC algorithms are
  untouched.
  ``LSGrid::fillBdc``, the other ``setFromTriplets`` in the DC path (6.6% of a DC solve), is
  deliberately left alone: its entries are unsorted and heavily duplicated -- every branch writes
  both its end buses' diagonals -- which is the one shape Eigen's assembler is already good at.
- [IMPROVED] ``NRSystem::build_J_sparsity`` writes J's compressed-column arrays itself instead of
  going through ``setFromTriplets`` and then looking every contribution back up in the result.
  **1.5-1.6x on the phase that contains it** -- the ``pre_proc`` timer, which also covers
  ``update_state`` and ``init_topology``, goes 5.06-5.21 -> 3.36-3.39 ms per solve on
  case9241pegase, 1.39-1.46 -> 0.88-0.90 on case2869pegase, 0.60-0.63 -> 0.38-0.39 on
  case1354pegase, and a rebuild solve on case9241pegase drops 2% of its total instruction count.
  The function itself goes 58.7M -> 40.4M instructions for three rebuilds (**1.45x**, and 1.77x
  against 1.0.0).
  ``setFromTriplets`` runs the same two counting sorts this function needs -- bucket the entries by
  row, then transpose into columns, the transpose being what sorts each column -- but it carries a
  double per entry through a temporary ``SparseMatrix``, collapses duplicates in a pass of its own,
  and materialises the transpose as a second matrix; and since it hands back only a matrix, finding
  where each contribution landed took ``4 * nnz(Ybus)`` binary searches afterwards (~170k per
  rebuild on a 9241-bus grid, 12.1M instructions of the 58.7M, with Eigen's own machinery
  accounting for 23.5M more). Doing the two sorts here carries a 4-byte entry id instead of the
  value, folds the duplicate collapse into the transpose, and -- the point of the exercise -- knows
  each coefficient's position at the moment it writes its row index, so the dS maps and the feature
  positions are filled straight from that walk and the binary searches are gone entirely.
  Duplicate entries stay supported, which is what rules out the cheaper sorted-triplet path: a
  feature entry may legitimately land on a dS coefficient (an hvdc droop slope adds to the
  dP/dtheta of its end buses), and 856 of them turn up over the C++ test suite.
  The invariants Eigen used to provide are now asserted in the debug build, where they are
  verified over the whole suite: J compressed, every column filled exactly to its declared end,
  row indices strictly increasing. Beyond that, the build was cross-checked entry by entry against
  ``setFromTriplets`` + ``find_J_pos`` over 409 sparsity builds of the C++ test suite (which is
  what exercises multi-slack, voltage control and hvdc) and on the four benchmark grids: same nnz,
  same outer array, same inner array, same position for every single contribution. Results are
  bit-identical on case118 and the three PEGASE cases across NR / NRSing / NRRefactorRetry /
  FDPF XB / FDPF BX / DC with both SparseLU and KLU (32 configurations).
  Two things were tried and dropped because they measured worse: keeping the scratch buffers as
  members to skip the per-rebuild allocation (the extra indirection in the two hot walks cost more
  than the allocations saved -- ``pre_proc`` 3.65-3.87 ms against 3.33-3.40 for the local
  buffers, five runs out of five), and leaving J's value array uninitialised (``fill_J`` does write
  every coefficient, but a freshly built J being all zeros is what ``get_J_python`` and the
  assertions assume).
  ``Base::find_J_pos`` has no caller left in the solve path. It stays part of the component
  protocol -- a component or an external consumer still needs a way to locate a coefficient in a
  built J -- with its comment corrected to say so.
- [IMPROVED] ``NRSystem::build_J_sparsity`` hands its contributions straight to
  ``setFromTriplets`` instead of copying them into a second vector first. The generic dS pass
  records one 16-byte ``Contrib`` per Jacobian coefficient a dS matrix feeds (its J row and
  column, the Ybus nonzero it reads, and which of the four dS parts it takes), and that vector has
  to survive the build because the value maps are resolved from it afterwards. The code then
  allocated an equally large ``std::vector<Eigen::Triplet<real_type> >`` and refilled it with the
  same row/column pairs plus a literal ``0.`` -- 2.07 MB of it on case9241pegase, duplicating the
  2.07 MB already held -- purely because ``setFromTriplets`` wants ``row() / col() / value()`` and
  ``Contrib`` spelled them ``jrow() / jcol()`` and had no value at all. This pass builds the
  sparsity pattern only, so every value written was zero anyway; ``fill_J`` writes the numbers.
  ``Contrib`` now also answers to Eigen's triplet protocol, the feature entries the components
  declare are appended to the same vector so one range covers the whole pattern, and the maps are
  resolved over the leading dS entries. Worth 12.9M instructions of 908M on a case9241pegase
  rebuild solve (**-1.4% of everything, C++ side**): the 8.5M spent constructing triplets is gone,
  and ``set_from_triplets`` itself drops 12.5M -> 6.6M because the zero is now a constant the
  compiler folds into the store rather than a value loaded per entry. Results are bit-identical --
  the pattern handed to Eigen is the same, in the same order -- verified on case118 and the three
  PEGASE cases across NR / NRSing / NRRefactorRetry / FDPF XB / FDPF BX with both SparseLU and KLU.
- [IMPROVED] the bounds check on an element id is now compiled out with NDEBUG **when the id came
  from this library rather than from a caller**. Profiling a powerflow on case9241pegase found
  ``GenericContainer::_check_in_range`` and the ``_get_bus`` it guards reached 128k times per solve
  -- four passes over both ends of every branch, in ``fillYbus``, ``compute_results`` and
  ``reconnect_connected_buses`` -- all of them from loops shaped ``for(el_id = 0; el_id < nb();
  ++el_id)``, where the bound is a property of the loop and the check can only fail if there is a
  bug here rather than in the caller. Both were also out-of-line calls, not inlined compares: the
  error paths build an ``ostringstream``. Worth ~2.9% of everything lightsim2grid itself does in a
  rebuild solve (14.7M instructions of 919M).
  **Nothing a user can reach lost its check.** ``_check_in_range`` is unchanged and still throws
  for every id that crosses the python boundary -- ``get_bus_load`` / ``get_bus_gen`` /
  ``get_bus1_powerline`` and friends, ``deactivate`` / ``reactivate`` / ``change_bus``,
  ``change_ratio`` / ``change_shift``, ``set_regulated_bus``, ``set_status_droop``,
  ``update_slack_weights_by_id``. The new ``_check_in_range_internal`` is the debug-only twin, used
  by the new ``_get_bus_internal`` and by the ``get_bus_side_1_internal`` /
  ``get_bus_side_2_internal`` accessors the branch containers' own loops now call. The assertion and
  sanitizer CI builds keep every one of them, since ``USE_DEBUG_ASSERTS`` clears ``NDEBUG``.
  Verified from both sides: ``lightsim2grid/tests/test_LSGrid_out_of_bounds.py`` (10 tests, 82
  subtests -- the suite written for exactly this contract) still passes, and a release build
  (``-DNDEBUG``) still raises ``out_of_range`` for a bad id on each user-facing accessor. Results
  are bit-identical on case118 and the three PEGASE cases.
- [IMPROVED] ``NRSystem::fill_internal_variables`` -- the dS assembly, and the most expensive thing
  lightsim2grid itself does in a Newton solve (13.7% of one on case9241pegase, against 63% for KLU)
  -- computes both derivatives from ONE product instead of four. 1.26-1.39x on the phase, on every
  grid from 118 to 9241 buses, rebuild or cache-reuse; -36% of its instruction count under
  callgrind. Writing Y for Ybus(i, j), the formulas derived from pandapower are
  ``dS_dVm(i,j) = conj(Y.Vnorm_j).V_i`` and ``dS_dVa(i,j) = -conj(Y.V_j).i.V_i`` -- and Vnorm_j is
  just V_j / |V_j|, so the first is the second's product divided by a real. With
  ``B = conj(Y.V_j).V_i`` (two complex products), dS_dVm is ``B / |V_j|`` (a real scaling) and
  dS_dVa is ``-i.B`` (a swap and a sign flip). The diagonal's two extra terms share a product the
  same way: ``conj(Y.V_i - Ibus_i).V_i`` is ``B - conj(Ibus_i).V_i``, and its dS_dVm term
  ``conj(Ibus_i).Vnorm_i`` is that same value over |V_i|. The arithmetic is now written on real and
  imaginary parts rather than left to ``std::complex``, whose ``operator*`` carries a NaN-recovery
  branch after every product -- on its own, that part was worth 1.19-1.31x and was bit-identical.
  Two of the four nb_bus scratch vectors are gone with it: the unit phasors and the two
  ``conj(Ibus) . *`` products are each used once per bus, on the diagonal, so they are built there
  instead of materialised and carried through the solve.
  Values move by **less than one ulp relative** (measured 1.9e-16 on case9241pegase): the one
  rounding that differs is dividing by |V_j| after the product rather than before. This is the
  Jacobian, not the residual -- it steers the Newton step, it does not define the answer -- and
  across case118 and the three PEGASE cases, on NR / NRSing / SparseLU, iteration counts are
  unchanged everywhere and converged voltages agree to 7.3e-13, five orders inside the 1e-8
  tolerance. FDPF is untouched (it builds no Jacobian) and stays bit-identical.
- [IMPROVED] ``NRSystem::fill_J`` no longer zeroes the Jacobian before filling it: ~1.2x on the
  fill itself for a large grid (case9241pegase: 1.96 -> 1.63 ms per solve, and 2.28 -> 1.86 ms on
  the rebuild path). It zeroed because every write accumulated, several contributions being able
  to land on one coefficient -- so a value array a megabyte wide was rewritten at every
  factorisation, and every coefficient was then read before being written. Only the FEATURE
  entries actually need that: the four dS families live at (p_row, theta_col), (p_row, vm_col),
  (q_row, theta_col) and (q_row, vm_col), and since a row is a P equation or a Q equation and a
  column a theta unknown or a vm unknown -- never both -- they are pairwise disjoint, while within
  a family the ledger hands every bus its own row and column, so no two dS entries ever claim the
  same coefficient. They can therefore assign. What still accumulates is the feature entries,
  because a component may legitimately add to a dS coefficient (an hvdc droop slope adds to the
  dP/dtheta of its end buses), so those positions -- a handful, against every nonzero of J -- are
  the only ones zeroed. Values are unchanged to the bit (``0 + x`` is ``x``), checked on case118
  and the three PEGASE cases across NR / NRSing / SparseLU / FDPF.
- [ADDED] ``build_J_sparsity`` asserts, in debug builds, the two properties the above rests on: no
  Jacobian coefficient is claimed by two dS entries, and every stored coefficient has a writer (one
  that is never written would keep a stale value now that the array is not cleared). Checked where
  the layout is decided rather than assumed in the hot loop, and covered by the assertion / sanitizer
  CI builds. Two tests pin the same thing behaviourally: filling J twice without changing anything
  must not move it, and filling at one voltage state then another must match a system that only ever
  saw the second -- both with and without a droop hvdc line, which is the case whose feature entries
  land on dS coefficients.
- [IMPROVED] ``Base::find_J_pos`` takes the index arrays rather than the matrix, and
  ``build_J_sparsity`` reads them once instead of once per coefficient: ~1.1x on the pre-processing
  of a rebuild (case9241pegase: 7.97 -> 7.19 ms). It was binding an ``Eigen::Ref`` to J on every
  call, and it is called four times the Ybus nonzeros -- 170k times per solve on a 9241-bus grid,
  which callgrind showed as the single most-called function in a powerflow. The Ref itself is free
  to construct and does not copy (it aliases a compressed, same-Options matrix); paying for one per
  coefficient to read two pointers that never change across the loop is what cost. The
  matrix-taking overload stays, and now delegates.
- [IMPROVED] **Fast-Decoupled powerflow is ~1.25-1.35x faster** (whole solve, both the XB and
  the BX flavour, 49 to 1600 buses), for the same iteration count and the same answer.
  ``BaseFDPFAlgo::has_converged`` rebuilds V from (Vm, Va) and then has to put that pair back in
  canonical form -- magnitude >= 0, angle wrapped -- because the Q iteration (``Vm_(pq) -= q_``)
  can overshoot a magnitude past zero and the P iteration accumulates into ``Va_`` without ever
  wrapping it. It did that with ``Vm_ = V_.array().abs(); Va_ = V_.array().arg();``: a hypot and an
  atan2 per bus, twice per iteration, to rediscover two numbers the solver already holds. V is
  built as ``Vm_ * exp(i.Va_)``, so its modulus is just ``|Vm_|`` and its argument just ``Va_``,
  plus half a turn where the magnitude went negative -- a sign flip and an addition. On a 121-bus
  grid that pair cost 2.5 us per call out of a ~150 us solve. The complex voltage ``V_`` the solve
  actually consumes is built exactly as before and is unchanged bit for bit; ``Vm_`` / ``Va_`` now
  differ in their last bits (in their favour: an exact ``|Vm_|`` instead of a rounded
  ``hypot(Vm.cos, Vm.sin)``), which moves a converged solution by ~1e-13 -- four orders of
  magnitude inside the 1e-9 convergence tolerance, with iteration counts unchanged on all 80
  grid / loading / flavour combinations checked, and no change in the distance to a Newton-Raphson
  reference solution.
- [ADDED] ``BaseFDPFAlgo::canonicalise_vm_va``, the (magnitude, angle) canonicalisation above, as a
  named static so the identity it rests on can be tested head-on. It has to be: its two interesting
  branches -- a negative magnitude, an angle several turns out of range -- turn out to be
  unreachable from any converging solve (a trajectory that overshoots a magnitude ends up
  diverging, and a diverged solve clears ``Vm_`` / ``Va_`` / ``V_``), so no solve-level test can
  reach them.
- [ADDED] ``src/tests/test_fdpf_algorithm.cpp``. The Fast-Decoupled family had no answer-level
  coverage at all -- it appeared in the suite only through ``test_cache_reuse.cpp`` (caching, not
  numbers) and ``test_plugin_registration.cpp`` (names) -- so the change above landed on untested
  code. It now pins that both flavours converge to the Newton-Raphson solution of the same grid,
  that the reported ``Vm`` / ``Va`` stay consistent with the reported ``V``, and that
  ``canonicalise_vm_va`` agrees with the ``abs()`` / ``arg()`` pair it replaced over Vm in [-3, 3]
  x Va in [-10, 10] rad. The one place they deliberately differ is Vm == 0 exactly, where the phase
  of a zero phasor is undefined (``arg()`` returned whichever of 0 / +pi / -pi atan2's signed-zero
  rules produced) and where the solve cannot converge anyway -- the mismatch is divided by that
  zero on the next line and fails its ``allFinite()`` check.
- [IMPROVED] **Gauss-Seidel is ~12x faster** (121-bus mesh, 5568 sweeps: 470 ms -> 37 ms; the
  synchronous variant, 1.9x: 422 ms -> 228 ms), for the same sweeps and bit-identical voltages.
  ``Ybus`` reaches a solver column-major, but a Gauss-Seidel sweep needs one ROW of it at a time --
  and a row of a column-major sparse matrix is not stored anywhere: reading it means scanning every
  column. ``GaussSeidelAlgo::one_iter`` did exactly that (``Ybus.row(k) * V_``, once per pq bus and
  twice per pv bus), which made one sweep cost O(nb_bus * nnz) instead of O(nnz); and it re-found
  ``Ybus(k, k)`` by binary search inside column k on every lookup. Both now come from caches
  rebuilt once per ``compute_pf`` (an O(nnz) transpose, against the hundreds or thousands of sweeps
  that follow): a row-major copy of ``Ybus`` and its diagonal. The row-major inner iterator visits a
  row by increasing column index, which is the order Eigen's sparse * dense product summed it in,
  so the sums are unchanged to the last bit.
- [IMPROVED] ``BaseAlgo::_evaluate_Fx`` (the Gauss-Seidel convergence check) computed its
  sparse-times-dense ``Ybus * V`` product **three times** per call, and the multi-slack overload
  twice. The mismatch was kept as a lazy Eigen expression (``auto mis = ...``), and a lazy
  expression holding a sparse * dense product re-evaluates that product every time the expression
  is evaluated -- once per assignment of the three result segments. It is a concrete vector now.
- [IMPROVED] the Newton-Raphson iteration allocates less per iteration: ~7-14% off the
  algorithm-side cost of a solve (everything but the linear solver: dS assembly, Jacobian fill,
  mismatch, voltage update), which on a 121-bus grid is about a fifth of an NR solve with KLU and
  most of one without a fast linear solver. ``NRSystem``'s residual, trial-voltage and dS assembly
  write into buffers that live as long as the system instead of heap-allocating and freeing a fresh
  one on every call -- which matters most for the line-search / Iwamoto step policies, whose
  backtracking evaluates the residual up to ``ls_max_iter`` times per iteration.
  In particular every ``Ybus * V`` product now lands in a named buffer with ``.noalias()`` rather
  than inside a coefficient-wise expression Eigen has to allocate a temporary to evaluate, and
  ``NRAlgo`` scales its step in place instead of handing ``apply_step`` a ``coeff * F`` expression
  that its ``Eigen::Ref`` parameter had to materialise. Values are unchanged bit for bit.
- [IMPROVED] ``NRSystem`` no longer keeps ``dS_dVm`` / ``dS_dVa`` as two full ``Eigen::SparseMatrix``
  copies of ``Ybus``. Nothing ever indexed them by (row, column): ``fill_internal_variables`` writes
  them and ``fill_J`` reads them purely by nonzero position, so they are plain value arrays, one
  coefficient per ``Ybus`` nonzero. Each topology change is two structure copies (outer + inner
  index arrays) lighter -- the "TODO speed: copy only the sparsity pattern and not the values" that
  used to sit in ``init_topology``.
- [ADDED] ``NRSystem::mismatch_into(res)`` and ``NRSystem::mismatch_sq_norm_at_current()``: the
  allocation-free forms of ``mismatch()`` and ``mismatch_sq_norm_at(<zero step>)``, which the NR
  loop and the step-scaling policies now use. The value-returning forms stay, unchanged, for
  out-of-tree algorithms.
- [BREAKING] everything a solver family caches now lives in one object, ``SolverSideCache``
  (``src/core/SolverSideCache.hpp``), instead of eighteen separate ``LSGrid`` members. The bus
  labelling, ``Ybus`` / ``Sbus``, the slack, the PV / PQ split and each family's connectivity
  snapshot are one picture of the grid taken at one instant, expressed in ONE bus labelling; they
  are now built, handed over and retired as one unit. ``LSGrid`` holds ``ac_cache_`` / ``dc_cache_``;
  ``pre_process_solver`` / ``pre_process_dc_solver`` (C++ only, never exposed to python) no longer
  take a cache at all: they are defined by the fact that they build ``ac_cache_`` / ``dc_cache_`` for
  this grid's own algorithm. Building a solver input into a cache the CALLER owns -- what the batch
  algorithms do -- is now a separate pair, ``build_solver_input`` / ``build_dc_solver_input``, rather
  than the same function distinguishing the two by comparing the cache's address. The two really are
  different operations: the own-cache pair applies the reuse policy and resets / re-configures
  ``_algo``; the foreign pair never re-stamps (the change flags describe THIS grid and say nothing
  about someone else's cache), leaves this grid's algorithms alone, and publishes the labelling and
  PV / PQ split into this grid's cache for the NR extensions before retiring it. The assembly itself
  is shared, unchanged, in one private ``_build_into_cache``.
  **This changes ``LSGrid``'s member layout**, so anything that casts an ``LSGrid`` across a module boundary (``gpusim2grid``) must be rebuilt
  against these headers -- which the plugin ABI policy in ``docs/solver_plugin.rst`` already
  requires. Nothing changes for python.
- [FIXED] a powerflow that throws part-way through no longer leaves the grid claiming its solver-side
  cache is up to date. A solve rebuilds that cache in place -- the bus labelling, Ybus / Sbus, the
  PV / PQ split, the slack weights, the algorithm's factorization -- so a throw anywhere in
  ``ac_pf`` / ``dc_pf`` (a plugin solver raising, the wrong-sized-voltage rejection, an internal
  consistency error) left a mixture of the old grid and the new one behind, with the change flags
  still saying whatever they said before the call. ``ac_pf`` / ``dc_pf`` now run the whole solve
  against a *copy* of the change tracking and hold the grid itself at "everything changed, both
  families" for the duration; the copy becomes the grid's change tracking only at a publication
  statement placed after ``compute_results``, which a throw never reaches. So an interrupted
  powerflow leaves ``need_reset_solver()`` true on both families and the next one rebuilds from
  scratch -- by construction, whatever a future step throws and wherever. Deliberately not a
  ``try`` / ``catch``: an unwind edge through this code is not free (the one that used to guard the
  bus-counting bracket cost 4.9M instructions per solve without ever running), while the copy is 24
  bools -- measured at +21 instructions per solve on ``case9241pegase``, against 206M. Note the
  input validation at the top of ``ac_pf`` / ``dc_pf`` (Vinit size, ``max_iter`` / ``tol``, a droop
  grid handed to a solver that cannot do droop) stays *outside* that bracket: it runs before
  anything is touched, so a rejected call still costs the caller's cache nothing.
  ``build_solver_input`` / ``build_dc_solver_input`` (the batch-algorithm entry points) got the same
  treatment from the other side: they now retire this grid's cache *before* copying the caller's
  labelling into it rather than after, so a throw part way through that copy cannot leave a
  half-published mixture behind a control that still says the cache is up to date.
- [FIXED] a batch algorithm no longer writes half of what it builds into the grid's cache.
  ``slack_weights``, ``bus_pv`` and ``bus_pq`` were not parameters of ``pre_process_solver``: they
  were taken from the grid's own members whatever the caller passed for the other six containers.
  So a ``TimeSeries`` / ``ContingencyAnalysis`` build wrote its PV / PQ split and its slack weights
  into the grid's cache, next to a ``Ybus`` built for a different labelling, and left that grid
  still claiming the result was reusable; and the reuse guard compared one owner's container sizes
  against the other's, so a caller that reused its own containers with a "nothing changed" control
  would have skipped ``fillpv_pq`` and stamped the grid's PV / PQ split onto its own system --
  converged, plausible, wrong. Nothing reached that (every batch asks for a full rebuild, and works
  on a private copy of the grid), which is what made it worth closing before something did. A cache
  being one object removes the whole class of mistake: there is no way to hand a function half of
  it. What a foreign build still publishes into the grid -- the labelling and the split, which the
  NR extensions read back through ``lsgrid_ptr`` rather than from what the solver was handed -- now
  retires that grid's cache instead of leaving it a mixture.
- [FIXED] a build into a caller's cache no longer resets or reconfigures the grid's own solver.
  ``_algo.reset()`` / ``tell_solver_control()`` fired for a solve the grid would never perform,
  discarding a factorization its next powerflow would have reused.
- [FIXED] ``LSGrid.unset_changes()`` no longer papers over a change that has not been solved.
  Its contract is "a powerflow has run, past changes can be forgotten" -- but nothing checked that
  a powerflow *had* run, so calling it with a modification outstanding cleared the flags that
  described it and the next powerflow re-solved the OLD system: converged, plausible, and wrong by
  0.038 pu on a 14-bus grid with one line deactivated. It now refuses instead: each family is
  marked only if its change control reports nothing outstanding AND its cache is structurally
  intact and still matches the grid's connectivity; a family that cannot back the claim is retired
  and rebuilds on its next powerflow. It still applies to BOTH families, as it always has.
  Note a cache cannot detect this by self-inspection -- changing a line's impedance, a transformer
  tap or a load target leaves every vector size and every bus status exactly as they were -- which
  is why the new ``AlgoControl.nothing_changed()`` is the load-bearing test. The element containers
  remain the sole authority on what changed: they declare it through ``AlgoControl`` in their own
  modifiers, and ``unset_changes()`` only reads that back. The structural check beside it answers a
  different question (is the data there at all, and the right shape) and says nothing about change.
- [ADDED] ``AlgoControl.nothing_changed()``: has every change this control tracks already been
  consumed by a solve? The exact negation of ``tell_all_changed()``.
- [BREAKING] ``BINARY_FORMAT_VERSION`` 4 -> 6: nothing about bus connectivity is serialized any
  more. ``LSGrid``'s state loses the AC family's bus-connectivity photograph (5) and
  ``SubstationContainer``'s loses ``bus_status_`` (6). Neither is independent state: a bus is part
  of the powerflow if and only if at least one active element sits on it, and the elements' own
  status *is* in the file, so both were caches of something already there -- with a way for a
  crafted file to make the two disagree. A restored grid counts its buses from its elements. That
  photograph was also the last piece of solver-cache metadata a saved grid carried, so there is
  nothing left in the layout to poison.
- [ADDED] ``SubstationContainer`` keeps, per bus, how many elements hold it alive, maintained
  incrementally by every mutator that can change bus membership. A bus is in the solved system iff
  its count is non-zero, so the two transitions that matter (0 -> 1 and 1 -> 0) are exactly what
  raises ``tell_dimension_changed``. Every other change flag stays where the element containers
  decide it. Which buses an element holds is stated once and only once, in
  ``GenericContainer::contribute_to_buses``: the from-scratch recount and every mutator are built
  from it, so the incremental and the rebuilt answer cannot disagree. The second statement of that
  same rule, ``reconnect_connected_buses``, is deleted from all five containers -- it was the one
  that drifted (see the SVC fix below).
- [FIXED] a bus whose only element is a static VAR compensator now counts as connected.
  ``init_bus_status()`` never called ``SvcContainer::reconnect_connected_buses`` -- the container
  inherits a perfectly good one from ``OneSideContainer_PQ``, it was simply left out of the list --
  so such a bus was dropped from the solved system. An SVC injects reactive power; its bus belongs
  there.
- [FIXED] ``LSGrid.update_topo()`` -- the batch entry point grid2op drives every step -- now counts
  the buses it takes out of the system. It is not a thin wrapper over the individual mutators: it
  also calls ``resolve_status``, which sets a branch's ``status_global_`` (the gate
  ``contribute_to_buses`` reads *first*) and opens the opposite side through the
  ``*_no_bus_tracking`` entry points. Both change which buses the branch holds, and both ran
  OUTSIDE the counting bracket, so stranding a line end alone on a busbar and then disconnecting
  the line left that bus counted with nothing on it: the system silently lost a bus, nothing raised
  ``tell_dimension_changed``, and the next powerflow failed to initialise
  (``ErrorType.NotInitError``). The whole per-element update -- both sides *and* ``resolve_status``
  -- is now one bracket, and the sides use the ``*_no_bus_tracking`` entry points for the same
  reason ``deactivate()`` does: a line end does not own its contribution. As a side effect a
  one-sided element being reactivated and moved in the same entry is now one bracket rather than
  two.
- [FIXED] a ``change_bus`` the grid refuses no longer corrupts the per-bus element counts.
  ``GenericContainer::_change_bus`` rejects a bus id below 0 or past the last bus -- and ``-1`` is
  the "no bus" marker, so it is exactly the id a caller reaches for meaning "disconnect this".
  That rejection is raised from *inside* the counting bracket, after the element's contribution has
  been taken away and before it is put back, so the ``-1`` stood alone: every bus the element held
  stayed one short, silently and for the rest of the grid's life, on a call that changed nothing
  else. ``_apply_and_track_buses`` now restores the contribution on the way out of an exception as
  well as on the normal path. Found by review; covered for all twelve ``change_bus`` entry points,
  against both a negative and a past-the-end bus id.
- [FIXED] ``LSGrid.get_bus_status()`` and ``LSGrid.nb_connected_bus()`` no longer report every bus
  disconnected on a grid that has been built but not yet solved. The per-bus element counts are
  disarmed by everything the incremental ``+1`` / ``-1`` cannot see -- a freshly built grid,
  ``set_state`` / ``load_binary``, an ``init_*`` that replaces a whole element container -- and an
  all-zero count is exactly what "never counted" looks like. Only the powerflow re-established them
  (through ``init_bus_status()``), so ``init_from_matpower(...)`` followed by ``get_bus_status()``
  answered from the disarmed state. Reading connectivity now establishes the counts, the same way
  solving does.
- [FIXED] ``LSGrid.consider_only_main_component()`` now tells the solver its dimension changed.
  ``disconnect_if_not_in_main_component`` deactivates the elements it strands, so the buses whose
  last element it takes away leave the solved system -- which renumbers every bus after them. It
  handed those deactivations a local, throwaway ``DualAlgoControl``, so the 0-crossings were
  detected and dropped on the floor; the flag used to reach the real controller only because
  ``init_bus_status()`` rebuilt the status and compared it against each family's photograph. With
  the crossing *being* the notification, it now goes to the controller the solver reads. The
  two-sided containers additionally do their counting through the branch's own rule, in one bracket
  around both sides and the ``status_global_`` flip, instead of letting each side count for itself
  under the one-sided rule. This is what broke grid2op's ``automatically_disconnect=True``.
- [BREAKING] [DEPRECATED] ``LSGrid.deactivate_bus(bus_id)`` and ``LSGrid.reactivate_bus(bus_id)`` are
  now **no-ops** (kept so existing code keeps importing; they will be removed in a later version).
  A bus is part of the powerflow if and only if at least one active element sits on it, which is now
  the only statement of bus connectivity there is, so there is no separate switch left for them to
  flip. They were already all but inert: whatever they set, the next powerflow rebuilt the bus status
  from the elements and threw it away -- a bus with elements on it came straight back, a bus without
  them was already out -- which is why removing the effect changes no powerflow result. To take a bus
  out of the powerflow, disconnect the **elements** on it. The pandapower and powermodels loaders call
  these for out-of-service buses whose elements they also disconnect, and that is what did the work.
  (The C++ ``SubstationContainer::disconnect_all_buses`` / ``reconnect_bus`` / ``disconnect_bus`` /
  ``reset_bus_status`` are no-ops for the same reason.)
- [BREAKING] ``LSGrid.get_bus_status()`` now returns a **new** list of bool built from the per-bus
  element counts, instead of a reference into an internal vector. The values are unchanged (they are
  what the powerflow always used); what changes is that they are always in step with the elements,
  and that the result is a copy -- python callers who kept it around and expected it to update by
  itself now need to call it again.
- ``init_bus_status()`` no longer walks every element of eight containers on every powerflow that
  rebuilds, and no longer maintains a bus-status vector at all: a bus is connected iff its element
  count is non-zero, and ``nb_connected_bus()`` is one integer moved when a count crosses 0. In the
  steady state it does nothing. It recounts from the elements only where there is nothing to be up
  to date with -- a grid that has never been counted (freshly built, or restored by ``set_state`` /
  ``load_binary``). ``_mark_cache_valid`` no longer copies a ``std::vector<bool>`` per powerflow
  either, and ``nb_connected_bus()`` is a constant-time read rather than a walk over every bus.
  What that is worth is ``nb_bus`` times a constant, independent of the number of elements:
  callgrind measures ~13 instructions saved per bus per topology-changing powerflow, on all four
  of 1000x1, 1000x12, 5000x1 and 5000x12 (substations x busbars) -- 0.08% of a solve on the
  1000-bus grid, 0.93% on the 12 000-bus one. It is the sparsely-filled grids that gain: 5000
  substations at 12 busbars each is a 60 000-bus vector walked to solve a 5 000-bus system. Below
  the wall-clock noise floor in every case; the solve dominates.
- ``LSGrid``'s ``init_*`` methods (the ones that replace a whole element container) now disarm the
  per-bus element counts, so the next ``init_bus_status()`` recounts from the elements. Replacing a
  container wholesale is not something the incremental +1 / -1 bookkeeping can see, and counts that
  describe elements which no longer exist must not survive it.
- the powerflow path no longer re-verifies that the data behind a "nothing changed" claim exists.
  The two entry points that can make such a claim without having built anything --
  ``unset_changes()`` and ``check_solution()`` -- now verify it themselves, at an altitude where
  the check is free. What remains on the powerflow path is the ``allow_*_cache_reuse`` switch plus
  a debug-only assertion, so a future third claimant is caught by the C++ suite (run under ASan,
  UBSan and valgrind in CI) at no cost in release wheels.

[1.0.0] 2026-08-28
--------------------
- [BREAKING] **solver cache reuse is now automatic and on by default**, per solver family.
  A powerflow reuses what the previous one of the same family built -- the compact bus labelling,
  ``Ybus`` / ``Sbus``, the PV / PQ split, the slack weights -- and only re-stamps what the grid
  reports as modified since. Every powerflow marks its own family "in sync" on the way out, so
  ``LSGrid.unset_changes()`` is no longer needed: it is a no-op whenever reuse is enabled, and is
  kept purely for backward compatibility. Worth ~24% of the time of a powerflow on a 14-bus grid,
  and nothing has to be called to get it.
  This is breaking for code that modifies the grid through a path that bypasses ``LSGrid``'s own
  ``change_*`` / ``deactivate_*`` / ``reactivate_*`` methods (direct manipulation of the C++
  containers, for instance). Such code used to get away with it as long as it never called
  ``unset_changes()``; now it must say so with ``prevent_cache_reuse()`` (or the narrower
  ``tell_recompute_ybus`` / ``tell_recompute_sbus``), or turn reuse off altogether with
  ``allow_cache_reuse(False)``. Everything that goes through the public API is unaffected: a
  completeness sweep over every mutating ``LSGrid`` method (``src/tests/test_cache_reuse.cpp``)
  checks that the cached and the fully-rebuilt answer agree.
- [ADDED] ``LSGrid.allow_cache_reuse(bool)`` / ``allow_ac_cache_reuse`` / ``allow_dc_cache_reuse``
  and their ``get_*`` counterparts: turn cache reuse off (durably) for one family or both. The
  answer is identical either way, so this is a debugging switch ("is this wrong number a caching
  artefact?") and an escape hatch for code that mutates the containers directly.
- [ADDED] ``LSGrid.prevent_cache_reuse()`` / ``prevent_ac_cache_reuse()`` / ``prevent_dc_cache_reuse()``:
  a one-shot invalidation of what a family cached (as opposed to the ``allow_*`` mode above).
  ``prevent_cache_reuse()`` is the new name of ``tell_solver_need_reset()``, which still exists and
  behaves identically. Both are now documented rather than tagged "internal, do not use".
- [BREAKING] the AC and the DC solver families no longer share ANY solver-side data.
  ``slack_weights``, the PV split and the PQ split were single members overwritten by whichever
  family solved last -- unlike ``Ybus`` / ``Sbus`` / the id maps, which were already per family.
  That is what made a "nothing changed" claim unsound across families, and it prevented per-family
  cache reuse. Each family now owns its own, exposed as ``get_ac_pv_solver`` / ``get_dc_pv_solver``
  and the ``pq`` / ``slack_weights`` equivalents. The family-less ``get_pv_solver`` /
  ``get_pq_solver`` / ``get_slack_weights_solver`` are unchanged in spirit and answer for the AC
  family once an AC powerflow has run, falling back to DC otherwise -- which also fixes them:
  they used to relabel whatever the last solve left behind with the AC bus mapping.
- [BREAKING] ``LightSimBackend`` no longer carries the ``_last_dc`` workaround (a
  ``tell_solver_need_reset()`` before the first AC powerflow following a DC one, comment:
  "otherwise might segfault"). It existed only because the two solver families shared their
  cached solver-side data; they no longer do, so a DC powerflow before an AC one on the same
  backend is simply safe, and the AC family keeps its cache across it. The private
  ``_last_dc`` attribute is removed with it. Covered by ``TestDcThenAcSameBackend`` in
  ``lightsim2grid/tests/test_DoNothingACDC.py``.
- [FIXED] ``LSGrid.set_sn_mva()`` did not invalidate anything. ``sn_mva`` is the base power of the
  whole per-unit system: it is passed to ``fillYbus`` / ``fillBdc`` and divides ``fillSbus_me``, so
  a powerflow run after it, on a grid whose cache was considered in sync, solved the previous base
  power's system. Found by the completeness sweep described above.
- [FIXED] ``dc_pf()`` returned the voltage MAGNITUDES of a previous call. DC solves for angles
  only, so the magnitudes it reports are the ones it was handed -- but ``BaseDCAlgo`` skipped
  ``Vm_ = V.array().abs()`` unless ``_solver_control.has_v_changed()``, which asks whether the
  GRID's voltage setpoints changed and says nothing about the vector just passed in. Two
  ``dc_pf()`` calls with different ``Vinit`` therefore returned the first call's magnitudes for
  every bus not pinned by a controller, whenever the solver control said "nothing changed" -- ie
  after every ``unset_changes()``, which ``LightSimBackend`` calls after each step. The guard is
  removed: it saved one ``abs()`` over ``nb_bus`` doubles next to a sparse triangular solve of the
  same dimension.
- [FIXED] a diverging powerflow left the algorithm's internal state (half-converged iterate,
  factorization of a system it gave up on) in place for the next call, and threw away nothing.
  The algorithm is now reset on divergence, and the next powerflow of that family rebuilds its
  internals -- while the DATA built for the diverged solve, which is a correct picture of the grid,
  is kept and stays reusable.
- [ADDED] the rule that **nothing the solvers cache ever crosses a serialization boundary** is now
  pinned by tests (``[serialization]`` in ``src/tests/test_cache_reuse.cpp``,
  ``TestDeserializedGridStartsCold`` in ``lightsim2grid/tests/test_cache_reuse.py``). ``StateRes``
  carries the elements; the bus labelling, ``Ybus`` / ``Sbus``, the pv-pq split, the slack weights
  and the change-tracking flags are all derived from them and are rebuilt from them, so a grid
  restored from a pickle or a binary file always starts with a cold cache. That was already true
  (``set_state`` resets the solver state), but nothing checked it and the reset carries a "see if
  it's worth the trouble NOT to do it" TODO -- which is exactly the change these tests now block.
  The reasoning is a security one, not a performance one: a cache read from a file is a second copy
  of state that ``check_grid()`` cannot verify against the elements it claims to describe (it can
  validate an index; it cannot validate "this really is the admittance matrix of this grid"), and a
  well-formed-looking matrix would be solved without complaint.
- [FIXED] ``LSGrid.set_state()`` read the serialized bus-connectivity snapshot (``BUS_STATUS_ID``)
  back into the two solver-cache snapshots. It was immediately overwritten by the reset that follows,
  so nothing was ever wrong -- but it made a file the nominal source of cache metadata. The field is
  still written by ``get_state()`` (the layout is unchanged); it is simply no longer read back.
- [FIXED] ``LSGrid.unset_changes()`` could make the next powerflow **segfault**. It told the grid
  "the cached solver-side data matches me", for the AC *and* the DC family at once, and nothing
  checked that claim. Three sequences -- all of them documented usage -- therefore entered the
  "nothing to rebuild" path with a default-constructed (empty) id map / Ybus / Sbus and indexed
  them with bus ids in the hundreds::

      unset_changes(); ac_pf();           # nothing was ever built
      dc_pf(); unset_changes(); ac_pf();  # what was built belongs to the other family
      ac_pf(); unset_changes(); dc_pf();  # idem

  Release wheels are built ``-O3 -DNDEBUG``, so this was an out-of-bounds read, not an error.
  ``LightSimBackend.runpf`` already carried a python-side workaround for the second one (the
  ``self._last_dc`` ``tell_solver_need_reset()``, comment: "otherwise might segfault"), and
  ``check_solution`` guarded against the first one locally with ``id_me_to_ac_solver_.size() > 0``.
  The check now lives where the reuse decision is made instead: ``_pre_process_solver_impl``
  verifies, in a handful of integer comparisons, that the data the flags describe is actually
  there, and falls back to a full rebuild otherwise. A wrong "nothing changed" can now cost time,
  never memory safety.
- [ADDED] ``InjectionSweepCPP`` (python wrapper ``lightsim2grid.injectionSweep.InjectionSweep``), a
  third batch algorithm. It computes exactly what ``TimeSeriesCPP`` computes -- one powerflow per row
  of the injection matrices, on a fixed grid topology -- with exactly the same interface, and differs
  only in how each computation is initialized: ``TimeSeriesCPP`` warm starts a step with the solution
  of the step before it (they are consecutive instants of a *time* series), while ``InjectionSweepCPP``
  starts every step from the same voltage -- the one given to ``compute_Vs``, or the "n" powerflow
  result when ``init_from_n_powerflow`` is set -- exactly like ``ContingencyAnalysisCPP`` does for
  each of its contingencies.
  Use it when the "steps" are independent scenarios rather than consecutive instants: their results
  then do not depend on the steps computed before them, nor on the order they are given in, and the
  batch can be spread over several threads (``nb_thread``).
  Both classes are two instantiations of one C++ template (``ls2g::BaseInjectionSweep<BatchInitKind>``
  in ``batch_algorithm/BaseInjectionSweep.hpp``), so they cannot drift apart: the init policy is a
  single compile-time parameter, and everything else -- inputs, results, bindings -- is shared.
- [BREAKING] (c++ API only) ``batch_algorithm/TimeSeries.hpp`` / ``.cpp`` are renamed to
  ``batch_algorithm/BaseInjectionSweep.hpp`` / ``.cpp``: the file now holds both ``ls2g::TimeSeries``
  and ``ls2g::InjectionSweep`` (two aliases of the same template), so naming it after either one of
  them would be misleading. C++ code including the old path must update the include; nothing changes
  on the python side. Both class names are unchanged.
- [ADDED] ``nb_thread`` moved up from ``ContingencyAnalysisCPP`` to the common batch base class, so
  every batch algorithm whose computations are independent of one another gets it. It is now
  available on ``InjectionSweepCPP`` as well, where each thread solves a contiguous range of steps
  with its own solver (the admittance matrix is shared read-only, since the topology is fixed --
  unlike ``ContingencyAnalysisCPP``, which needs a copy per thread to emulate the disconnections).
  Results are bit-for-bit independent of ``nb_thread``. ``thread_init_time()`` moved up with it.
- [BREAKING] ``TimeSeriesCPP.nb_thread = n`` now raises a ``RuntimeError`` for ``n > 1`` (the
  attribute did not exist before, so no working code can break). Its steps are chained, so splitting
  them into per-thread ranges would break the chain and make the voltages depend on the number of
  threads; the error points at ``InjectionSweepCPP``, which computes the same injections in parallel.
  ``nb_thread = 1`` (and any value below, which is clamped to 1) stays legal.
- [FIXED] ``TimeSeriesCPP`` solved each step against an **incomplete injection vector**. It does not
  use the ``Sbus`` the gridmodel builds: it rebuilds its own, per step, out of the four matrices the
  caller passes to ``compute_Vs`` -- generator and static-generator ACTIVE power, and load active +
  reactive power -- plus the hvdc term taken from the grid. But ``LSGrid::fillSbus_me`` stamps rather
  more than that, and everything else was silently dropped:
  storage units (active AND reactive: absent altogether), REACTIVE_POWER-mode SVCs (absent
  altogether), the reactive setpoint of a generator whose ``voltage_regulator_on`` is false, a static
  generator's reactive setpoint, and -- in DC -- the phase-shifter term
  ``TrafoContainer::hack_Sbus_for_dc_phase_shifter``. The powerflow converged happily on an injection
  the caller never asked for; on a 4-bus feeder the individual omissions were worth 0.032, 0.052,
  0.066 and 0.081 pu of voltage error respectively. DC was affected too, since
  ``compute_one_powerflow`` falls back to the gridmodel's ``Pbus_`` only when NO complex ``Sbus`` is
  supplied, which is not the ``TimeSeries`` case.
  The remainder is now **derived rather than enumerated**: the complete per-unit injection the
  gridmodel just built (``Sbus_`` in ac, ``Pbus_`` in dc) minus the same reconstruction ``compute_Vs``
  performs, evaluated at the gridmodel's own target values. What is left over is by construction
  everything the per-step matrices do not carry, so an element type added to ``fillSbus_me`` later is
  picked up for free -- an enumeration falling behind is precisely what caused this bug.
  ``ContingencyAnalysisCPP`` (and ``SecurityAnalysis``, which wraps it) was never affected: it passes
  the gridmodel's own complete ``Sbus_`` straight to the solver.
- [ADDED] a stale-solver-state case to ``src/tests/test_batch_voltage_control.cpp``: run ``ac_pf``,
  then change the grid so a bus leaves the solver (which makes the AC maps stale AND too long, worse
  than merely empty), then ``dc_pf`` (which rebuilds only the DC maps), then a batch -- which inherits
  the grid's AC algorithm, since ``LSGrid::change_algorithm`` routes a DC type to the separate
  ``_dc_algo`` slot and so asking the GRID for DC never switches a batch over. The batch must still
  reproduce a clean single-shot ``ac_pf``, which it does because it works on a copy of the grid whose
  copy constructor ``reset()`` the solver state and whose ``prepare_solver_input_base`` rebuilds it
  under ``tell_all_changed()``. Being a C++ test it runs under the valgrind pass of the C++ suite and
  the ASan/UBSan + Eigen-assertion builds, where a too-long index that stays inside the heap would
  still be caught.
- [BREAKING] (c++ API only) ``batch_algorithm/BaseInjectionSweep.hpp`` / ``.cpp`` and
  ``batch_algorithm/ContingencyAnalysis.hpp`` / ``.cpp`` are removed, replaced by
  ``batch_algorithm/BaseBatchSweep.hpp`` / ``.cpp`` (plus two new small headers,
  ``batch_algorithm/YbusPolicy.hpp`` / ``.cpp`` and ``batch_algorithm/SbusPolicy.hpp`` /
  ``.cpp``): a single class template ``ls2g::BaseBatchSweep<YbusPolicy, SbusPolicy,
  BatchInitKind>``, policy-parameterized on what varies per step in the admittance
  matrix (nothing, or a contingency) and in the injection (nothing, or per-step
  values). ``TimeSeries``, ``InjectionSweep`` and ``ContingencyAnalysis`` are now three
  aliases of it -- same behavior, same C++ and python public API as before (the
  policy-only, contingency-only and injection-only methods are gated with
  ``std::enable_if`` so each alias only exposes what it always exposed); nothing
  changes on the python side for these three. C++ code including the old paths must
  update its ``#include``.
- [ADDED] ``ScenarioSweepCPP`` (python wrapper ``lightsim2grid.scenarioSweep.ScenarioSweep``),
  the 4th instantiation of ``ls2g::BaseBatchSweep``: varies both the injection AND a
  contingency (a line / trafo disconnection) per row, independently and row-aligned --
  row ``i`` of the injection matrices is solved together with row ``i`` of the
  contingency mask. Built up with a new setter-based API (see below) rather than a
  single bundled call: ``modify_gen_p`` / ``modify_sgen_p`` / ``modify_load_p`` /
  ``modify_load_q`` for the injection, ``set_contingency_lines`` /
  ``set_contingency_trafos`` (dense boolean masks, shape ``(n_simul, n_line)`` /
  ``(n_simul, n_trafo)``, ``True`` = "deactivate this branch for this row") for the
  contingency, then ``compute(Vinit, max_iter, tol)``. Deliberately a *different* API
  from ``ContingencyAnalysisCPP``'s ``add_n1`` / ``add_nk`` (a set of distinct
  contingencies applied to one shared base injection) -- the two usages do not unify
  cleanly, so ``ScenarioSweepCPP`` does not have ``add_n1`` / ``add_nk`` and
  ``ContingencyAnalysisCPP`` does not have ``set_contingency_lines`` /
  ``set_contingency_trafos``. See ``docs/scenario_sweep.rst``.
- [ADDED] ``ScenarioSweepCPP`` now also has ``ContingencyAnalysisCPP``'s
  ``handle_disconnected_grid`` mode and inline limit-violation checking
  (``compute_limit_violations`` / ``violation_threshold`` / ``get_violations`` /
  ``get_violations_n``), same names and semantics on both classes, and both return the
  same ``LimitViolation`` objects. Deliberately **no** ``converged`` / ``converged_n`` on
  ``ScenarioSweepCPP``: a non-converged row's ``get_violations()`` entry already carries a
  ``GRID`` / ``NOT_SIMULATED``-or-``DIVERGENCE`` sentinel violation, so a separate
  convergence flag would be redundant (``ContingencyAnalysisCPP``'s own ``converged`` /
  ``converged_n`` are unchanged). A diverging pre-batch "n" powerflow (the shared base
  case every row in the batch is solved relative to) now also stamps this same
  ``GRID`` / ``DIVERGENCE`` sentinel into ``get_violations_n()`` and every row of
  ``get_violations()``, instead of leaving them as empty lists indistinguishable from
  "converged, nothing found". As on ``ContingencyAnalysisCPP``, toggling
  ``compute_limit_violations`` clears the whole object (registered contingencies,
  injections, ``handle_disconnected_grid``, ...) -- set it *before* configuring
  anything else. ``lightsim2grid.scenarioSweep.ScenarioSweep`` gains a ``run()`` method
  reusing the existing ``PreContingencyResult`` / ``ContingencyResult`` /
  ``SecurityAnalysisResult`` dataclasses from ``lightsim2grid.contingencyAnalysis``
  (``contingency_name`` is always ``None`` there, unlike ``ContingencyAnalysis.run()``).
- [ADDED] a new setter-based API -- ``modify_gen_p`` / ``modify_sgen_p`` /
  ``modify_load_p`` / ``modify_load_q`` + ``compute(Vinit, max_iter, tol)`` -- is now
  also available on ``TimeSeriesCPP`` / ``InjectionSweepCPP``, alongside their existing
  bundled ``compute_Vs`` call (kept, not removed, deprecated via docstring only). Unlike
  ``compute_Vs``, which requires all four matrices every call, any axis you never set
  defaults to the grid's own current target value, broadcast across every row --
  applies to ``ScenarioSweepCPP`` too. Every setter shares one row-count lock: the
  first one called fixes the number of simulations, every later one (across all
  setters, including ``set_contingency_lines`` / ``set_contingency_trafos`` on
  ``ScenarioSweepCPP``) is checked against it immediately, instead of only at
  ``compute()`` time as ``compute_Vs``'s equivalent check did.
- [FIXED] ``ScenarioSweepCPP`` (non-``handle_disconnected_grid`` path) used to abandon every
  row after the first one that failed to converge (a contingency that islands the grid is a
  routine, expected outcome for this class, not an exceptional one) -- later rows were left as
  all-zero voltages with an empty, misleadingly "converged, nothing found"-looking violation
  list. Each row of ``ScenarioSweepCPP`` is independent (unlike ``TimeSeriesCPP``, whose rows
  chain), so it now keeps going past a failing row, exactly like ``ContingencyAnalysisCPP``
  already did; ``get_status()`` still correctly reports 0 for the batch as a whole.
- [FIXED] ``ScenarioSweep.run()`` (the python wrapper) raised on any per-row failure instead of
  reporting it through the sentinel-violation mechanism its own docstring describes --
  ``compute()`` is now called with ``ignore_errors=True`` internally, matching
  ``ContingencyAnalysis.run()``'s behavior (only a diverging shared pre-batch "n" case, which
  makes the whole result meaningless, still surfaces as an unmistakable ``GRID`` /
  ``DIVERGENCE`` sentinel rather than a silent empty list).
- [FIXED] ``ScenarioSweepCPP::compute_flows`` / ``compute_power_flows`` did not detect a
  ``modify_*`` / ``set_contingency_lines`` / ``set_contingency_trafos`` call made after
  ``compute()`` without a following ``compute()`` -- unlike ``ContingencyAnalysisCPP``'s
  equivalent guard (a registered-contingency-count check), a count comparison cannot catch this
  on ``ScenarioSweepCPP``, since its row-count lock already forbids changing the row count
  itself; only the row *content* changes. Both methods now raise in that situation instead of
  silently mixing a new mask/injection with stale voltages from the previous ``compute()``.
- [FIXED] ``lightsim2grid.scenarioSweep.ScenarioSweep``'s ``compute_limit_violations`` setter
  reset ``self.__computed`` but not the cached ``_line_mask`` / ``_trafo_mask`` (there is no
  C++-side getter for them), even though it clear()s the masks C++-side -- the next ``run()``
  could derive ``element_ids`` / ``element_names`` from a stale cache. Both are now reset too.
- [FIXED] ``modify_sgen_p`` (on both ``lightsim2grid.timeSerie.TimeSerie`` -- and so, by
  inheritance, ``lightsim2grid.injectionSweep.InjectionSweep`` -- and
  ``lightsim2grid.scenarioSweep.ScenarioSweep``) validated ``ndim`` but not the column count,
  unlike ``modify_gen_p`` / ``modify_load_p`` / ``modify_load_q``; the C++ side still caught it,
  just with a less immediately actionable error. Consistency fix only.
- [ADDED] ``src/tests/test_timeseries_sbus.cpp`` and ``lightsim2grid/tests/test_timeseries_sbus.py``:
  a one-step time series fed the grid's OWN injections must reproduce ``ac_pf`` / ``dc_pf`` exactly.
  Covers each dropped element kind separately and all of them together, in ac and in dc, plus a
  multi-step case (the remainder is constant, so it has to reach every row) and a control grid whose
  whole injection the per-step matrices already covered. Eight of the nine C++ cases fail without the
  fix.
- [FIXED] a generator could not remotely regulate a bus that another generator already
  regulated *locally* -- the grid was refused outright with
  ``LSGrid::fill_voltage_control_solver_data: generator N regulates bus X which has no voltage (Vm)
  unknown`` -- even though the configuration is perfectly well posed. Two machines holding one bus
  simply form a control group: one bordered voltage row, one reactive unknown each, and a sharing
  row splitting the reactive duty in proportion to their ``qmax - qmin`` ranges, exactly as two
  *remote* controllers on one bus already did.
  The cause was the order in which buses are classified. ``GeneratorContainer::fillpv`` knew how to
  keep a controller's OWN bus out of the PV set (a remote regulator leaves its own bus PQ), but
  nothing stopped a LOCAL regulator from claiming a bus that somebody else regulated remotely. That
  bus then had no ``Vm`` unknown, so ``VoltageControl::register_in`` got ``vm_col == -1`` for it, the
  group's voltage row was never given its ``+1`` entry, and the row came out structurally EMPTY --
  a singular Jacobian. The check that refused the grid was really reporting the misclassification,
  not a property of the configuration.
  Bus classification now runs the other way round: ``LSGrid::get_group_controlled_buses`` collects
  every bus an active *remote* generator or an active voltage-mode SVC aims at, ``LSGrid::fillpv_pq``
  keeps those out of PV, and ``fill_voltage_control_solver_data`` enrols the local regulators sitting
  on them -- generators and voltage-regulating hvdc converter stations alike -- as members of the
  group. Buses regulated only by things standing on them are deliberately untouched: several machines
  sharing a bus and all regulating it locally remains the ordinary PV case, with the per-bus reactive
  redistribution of ``GeneratorContainer::set_q``.
  Three configurations start working as a result: a remote controller onto a locally regulated bus,
  a remote controller onto a slack bus pinned by its own generator (the slack generator joins the
  group and MultiSlack grants the free ``Vm``), and -- separately, by giving the regulated-bus check
  the ``has_free_q`` escape hatch the controller-bus check already had -- a remote controller onto a
  slack bus that nothing pins. Disagreeing setpoints between a local and a remote regulator of one
  bus are now reported by the group's own check (``conflicting voltage setpoints``) rather than as a
  bus-classification error.
  Two things stay refused, on purpose. A remote regulator SHARING ITS OWN BUS with a local regulator
  has no solution at all: whatever reactive power it injects there is absorbed one-for-one by the
  machine holding that bus' magnitude, so it has no influence on its remote target -- no
  reclassification can help, and the "its OWN bus has no reactive (Q) equation" error still fires.
  Voltage-regulating hvdc converter stations take part on the same footing as generators. A VSC
  station with ``voltage_regulator_on`` pins its own bus through the very same PV path, so it is
  enrolled into the group whenever a controller claims that bus -- two new controller kinds,
  ``VoltageControlSolverData::HVDC_SIDE_1`` / ``HVDC_SIDE_2``, with the hvdc LINE id as ``elem_id``
  and the solved reactive output written back to the station. Its sharing key is a reactive range in
  MVAr, the same currency as a generator's, so a mixed generator/station group shares reactive power
  correctly (a group of a local generator, a remote generator and a station splits ``Q`` three ways
  by range).
  Note for voltage-mode SVCs: their regulated bus now leaves the PV path as well, which is what makes
  the classification uniform, but an SVC still has to be the only controller of its bus in v1 (its
  sharing key is a susceptance range, not a reactive range, so mixed SVC/generator sharing has no
  defined semantics yet). Such a grid is therefore still refused -- now by that restriction, which
  names the real limitation.
- [ADDED] ``src/tests/test_voltage_control_reclassify.cpp`` and
  ``lightsim2grid/tests/test_voltage_control_reclassify.py``: the classification rule itself, the
  newly-supported configurations (including the reactive split, and hvdc stations sharing a bus with
  generators on either end of the line), the ones that must stay refused, and regression guards that
  several purely local regulators on one bus -- and a lone voltage-regulating station -- still take
  the classical PV path and produce bit-identical voltages.
- [FIXED] remote voltage control (remote-regulating generators and voltage-mode SVCs) and the free
  ``Vm`` unknown of a distributed-slack participant were **silently ignored** by every batch
  algorithm -- ``TimeSeriesCPP``, ``ContingencyAnalysisCPP`` and ``SecurityAnalysis``, ie everything
  going through ``BaseBatchSolverSynch``. The powerflow still converged and still returned plausible
  voltages, they were just the voltages of a grid with no such regulation: a silent wrong answer, not
  an error. ``LSGrid::ac_pf`` was unaffected.
  The NR extensions do not receive their data as solver arguments; they reach back into the
  ``LSGrid`` through the algorithm's ``lsgrid_ptr`` and read its *member* solver-side labelling
  (``LSGrid::get_free_vm_slack_solver_buses`` for ``Base``,
  ``LSGrid::fill_hvdc_droop_solver_data`` for ``Hvdc``,
  ``LSGrid::fill_voltage_control_solver_data`` for ``VoltageControl``).
  ``ac_pf``/``dc_pf``/``check_solution`` pass those members straight to ``pre_process_solver``, so
  they are always current there; the batch algorithms own local mapping vectors and pass those
  instead, against a private *copy* of the grid whose copy constructor ``reset()`` all of it to
  empty. ``_pre_process_solver_impl`` wrote only the FORWARD map (``id_me_to_ac_solver_``) back onto
  the grid, so ``id_ac_solver_to_me_`` and ``slack_bus_id_ac_solver_`` stayed empty -- and empty is
  indistinguishable from "this grid has no controller": ``fill_voltage_control_solver_data`` returns
  on its ``nb_bus_solver == 0`` guard and ``get_free_vm_slack_solver_buses`` on its ``slack.empty()``
  one, both without a word. The sync now covers the reverse map and both slack vectors, AC and DC.
  The ``Hvdc`` angle-droop extension was never affected: it only ever read the forward map, the one
  that was already synced (there is now a regression test pinning that too).
  Two consequences worth knowing about. A grid whose voltage-control configuration is not supported
  in v1 now raises from ``fill_voltage_control_solver_data`` in the batch path as it always did in
  the single-shot one, where before it was quietly accepted and solved without the regulation. The
  restriction is on buses whose voltage magnitude is *already determined* -- a controller may not
  regulate the slack, nor a PV bus pinned by a local regulator, since neither owns a ``Vm`` unknown
  for the bordered row to act on; several controllers sharing one regulated bus is fine and always
  was (one voltage row, one ``Q`` column each, plus a sharing row splitting ``Q`` by
  ``qmax - qmin``). And under
  ``ContingencyAnalysis``'s ``handle_disconnected_grid`` mode, a contingency that strands a regulated
  bus is now reported as a ``DIVERGENCE`` instead of converging with the regulation dropped: bus
  masking is deliberately a value-only edit that must not touch the ``J`` pattern, so the bordered
  voltage row survives against a bus the mask pins, and the system has no solution. An honest
  non-result rather than a confident wrong one.
- [ADDED] ``src/tests/test_batch_voltage_control.cpp``: solves the same grid single-shot through
  ``LSGrid::ac_pf`` and through ``TimeSeries`` / ``ContingencyAnalysis``, and requires the voltages
  to agree -- for a remote-regulating generator, two generators sharing one regulated bus, a
  voltage-mode SVC (flat and sloped), a free-``Vm`` distributed-slack participant, and an hvdc angle
  droop.
- [ADDED] ``TimeSeriesCPP`` / ``ContingencyAnalysisCPP`` (and so ``SecurityAnalysis``, which wraps the
  latter) gained a string-based ``change_algorithm(name)`` overload and a ``get_algo_name()`` accessor,
  matching what ``LSGrid`` already had (``AlgorithmSelector::change_algorithm(const std::string&)``) --
  needed to select a plugin solver, or a built-in with no dedicated ``AlgorithmType`` enum member (eg
  ``NRRefactorRetry_KLU``), for a batch computation, which was previously impossible: only the
  enum-based overload was exposed, and ``AlgorithmType::Custom`` (what such solvers collapse onto) is
  rejected by ``change_algorithm(AlgorithmType)``. Also bound ``available_algorithm_names()`` on both
  classes for the same reason -- previously only ``available_default_algorithms()`` (the enum-only,
  plugin-blind list) was exposed.
- [FIXED] two bugs this exposed, both the same "propagate the algorithm via ``AlgorithmType`` instead of
  by registry name" mistake: (1) ``BaseBatchSolverSynch``'s constructor inherited the source ``LSGrid``'s
  AC algorithm via ``get_algo_type()`` (the enum), so constructing a ``TimeSeries`` /
  ``ContingencyAnalysis`` / ``SecurityAnalysis`` from a grid running a plugin or no-enum-member solver
  raised immediately; (2) ``ContingencyAnalysis``'s multi-threaded path (``nb_thread > 1``) rebuilt each
  per-thread solver the same lossy way. Both now propagate by name (``get_algo_name()`` /
  ``get_algo().get_name()``). See ``test_batch_algorithm_solver_selection.py`` for the regression tests.
- [ADDED] closed the last two gaps in the benchmark narrative-generation work above:

  - ``benchmarks/benchmark_grid_size.py`` gained a ``generate_narrative`` (recycling vs no-recycling,
    ``TimeSerie`` vs a regular grid2op step, ``ContingencyAnalysis`` vs a regular grid2op step, each as a
    "between Nx and Mx" range across the grid-size scan) and prints it after its tables; the corresponding
    "Comments" section of ``docs/benchmarks_grid_sizes.rst`` (previously absent) is filled in with the
    generated text. The page's TL;DR table is (and always was) printed directly by the same script run
    that produces its 4 detailed tables, so it could not actually drift against them the way it first
    looked -- a clarifying note was added since the visual duplication invited that reading.
  - ``benchmarks/benchmark_binary_serialization.py``'s single-grid comparison table gained a "speedup vs
    pickle" column (computed the same way the grid-size-scan table below it already does it), replacing
    the hand-written "``save_binary`` is **1.5x** faster..." sentence that used to sit below the table as a
    second, independently-typed copy of the same ratio; ``docs/binary_serialization.rst`` updated to match.
- [FIXED] the orchestration layer under ``benchmarks/`` had drifted from the scripts it is supposed to
  run: ``benchmarks/security_analysis.py`` was renamed to ``benchmarks/contingency_analysis.py`` (its
  ``examples/`` counterpart was already renamed in a previous fix, this one was missed --
  ``ContingencyAnalysis`` has been the class name for a while, the benchmark script's filename was the
  last place still carrying the old ``SecurityAnalysis`` name); ``benchmarks/benchmark_ts_ca.sh`` now
  calls that new name, and also runs ``benchmark_ca_nb_threads.py`` (the ``nb_thread``-scaling benchmark
  used in ``docs/security_analysis.rst``), which it previously skipped entirely;
  ``benchmarks/run_all_benchmarks.sh`` had a stale ``./benchmarks_grid_size.sh`` (that file has never
  existed -- the real one is ``benchmark_grid_size.sh``, singular) and never ran
  ``benchmark_binary_serialization.py`` at all; both are now fixed / added.
- [ADDED] ``benchmarks/benchmark_solvers.py`` now generates, from the numbers measured during the run
  itself, the descriptive text that comments on the "computation time" and "differences" tables
  (``generate_narrative`` function). Previously this text lived only in ``docs/benchmarks.rst``, hand
  written and hand updated after each run, which routinely drifted out of sync with the tables above it
  (wrong speed ups, stale percentages, and a unit mistake -- some durations were labelled "ns" while
  actually being microseconds). The script now prints this text (and saves it next to the other
  ``--save_results`` outputs, as ``description.rst``) so updating the docs after a benchmark run is a
  copy / paste instead of manual float formatting.
- [ADDED] the same treatment as above extended to the rest of the ``benchmarks`` folder: (1)
  ``benchmarks/benchmark_dc_solvers.py`` gained an equivalent ``generate_narrative`` and now prints /
  saves the descriptive text for ``docs/benchmarks_dc.rst`` -- whose hand written "TL;DR" numbers had
  drifted out of sync with (and in some cases no longer resembled) its own tables; (2)
  ``benchmarks/compare_lightsim2grid_pypowsybl.py`` can now run several (by default all 6) ieee cases in
  one process and prints, at the end, the 5 tables and narrative text found in
  ``docs/comparison_with_pypowsybl.rst``, instead of requiring one invocation per case and manually
  cherry-picking numbers into the right table cell (``benchmarks/benchmark_pypowysbl.sh`` updated
  accordingly); (3) ``benchmarks/security_analysis.py`` and ``benchmarks/time_serie.py`` now print the
  final "N times faster than raw grid2op" summary line themselves (computed from the same numbers printed
  just above it), instead of leaving it to be hand-restated in ``docs/security_analysis.rst`` /
  ``docs/time_series.rst`` -- where it had drifted to a value inconsistent with the printed figures.
- [ADDED] a "C++ standards" GitHub Actions workflow
  (``.github/workflows/cpp-standards.yml``) that compiles both the standalone
  C++ unit test suite and the python bindings twice: once pinned to C++14
  (the oldest standard the CMake auto-detection cascade claims to support)
  and once pinned to C++26 (the latest). Every other C++ workflow lets the
  compiler auto-pick its newest supported standard, which never actually
  exercises the oldest end of the claimed range and only exercises "the
  latest" by accident of whichever standard the runner's default compiler
  happens to support that day. The pin is a new ``LS2G_CXX_STANDARD`` CMake
  cache variable (root ``CMakeLists.txt``, ``src/core/CMakeLists.txt``,
  ``src/tests/CMakeLists.txt``), settable via
  ``-DLS2G_CXX_STANDARD=14`` / ``pip install -C "cmake.define.LS2G_CXX_STANDARD=14"``;
  left empty (the default), behavior is unchanged. The C++26 job installs g++-14
  and upgrades CMake (>= 3.30 is required for CMake to recognize C++26 support at
  all), since ubuntu-latest's default toolchain does not reliably provide either.
- [FIXED] the three C++14-C++2x auto-detection cascades (root ``CMakeLists.txt``,
  ``src/core/CMakeLists.txt``, ``src/tests/CMakeLists.txt``) were inconsistent with
  each other: only ``src/tests/CMakeLists.txt`` tried C++26 first, and
  ``src/core/CMakeLists.txt`` skipped C++20 entirely. All three now try the same
  sequence, C++26 down to C++14.
- [FIXED] ``LightSimBackend.set_algo_type`` used to reject anything that was not an
  ``AlgorithmType`` enum value, so a plugin or string-only built-in algorithm (eg
  ``NRRefactorRetry_KLU``, which has no ``AlgorithmType`` enum value at all) could only
  be selected through the private ``env.backend._grid.change_algorithm(...)`` -- which
  does not survive ``env.reset()`` / ``backend.copy()`` (it silently reverts to whichever
  ``AlgorithmType`` was last set through ``set_algo_type``). ``set_algo_type`` now also
  accepts a plain ``str`` (any name returned by
  ``env.backend._grid.available_algorithm_names()``, including string-only built-ins and
  registered plugins), validated the same way as before, and this choice now survives both
  ``env.reset()`` and ``backend.copy()`` -- exactly like an ``AlgorithmType`` value already
  did. The ``algo_type`` constructor kwarg accepts a string too, for the same reason.
- [ADDED] ``LightSimBackend.set_ac_algo_config`` / ``set_dc_algo_config`` (and their
  ``get_*`` counterparts): the persistent, supported way to customize the AC/DC
  ``AlgoConfig`` (scaling / refactor policy parameters). Unlike calling
  ``env.backend._grid.set_ac_algo_config(...)`` directly, which is silently reverted on the
  next ``env.reset()``, the customization now survives ``env.reset()`` and is preserved by
  ``backend.copy()``. Both found and fixed while auditing ``docs/algorithm_names.rst`` /
  ``docs/solvers.rst`` for PR review comments, and confirmed with a dedicated reproduction
  (a backend using ``NRRefactorRetry_KLU`` with a custom ``ScalingPolicyType`` used to lose
  both across a reset). See ``test_algo_reset_copy_persistence.py`` for the regression
  tests (algo type / algo config, each across ``reset()`` and ``copy()``).
- [FIXED] a regression from the previous entry: ``AlgoConfig`` (pybind11) supports
  neither pickling nor ``copy.deepcopy``, so the first implementation of
  ``set_ac_algo_config`` / ``set_dc_algo_config``, which stored the ``AlgoConfig``
  object itself as a backend attribute, broke ``pickle.dump(env.backend, ...)`` (as
  used by ``test_save_load`` in ``test_pickleable.py``) the moment either had ever
  been called, with ``TypeError: cannot pickle 'lightsim2grid.lightsim2grid_cpp.AlgoConfig'
  object``. The backend now stores the config's ``int_params`` / ``real_params`` as a
  plain (picklable, deepcopy-able) tuple of lists instead, and rebuilds a real
  ``AlgoConfig`` from it whenever one needs to be re-applied (``reset()``, ``copy()``).
  Added ``test_backend_still_picklable_with_custom_algo_config`` (changes both the algo
  type and its ``AlgoConfig``, then pickles and reloads ``env.backend``) to
  ``test_algo_reset_copy_persistence.py``.
- [DEPRECATED] ``LightSimBackend``'s ``solver_type`` constructor kwarg and its
  ``set_solver_type`` method, in favour of ``algo_type`` / ``set_algo_type``: "solver"
  now refers specifically to the *linear* solver (KLU, SparseLU, NICSLU, CKTSO), not the
  powerflow algorithm nor the combination of both that this kwarg actually selects (see
  ``docs/algorithm_names.rst``). ``solver_type`` still works and is mapped to
  ``algo_type``; passing both with different values raises a ``BackendError``.
- [FIXED] ``-Wsuggest-override`` warnings on ``GeneratorContainer``,
  ``SvcContainer`` and ``ConverterStationContainer``'s ``_change_p`` /
  ``_deactivate`` / ``_reactivate`` / ``_change_bus`` overrides: they were
  marked ``final`` but not ``override``, unlike every other class overriding
  these same ``OneSideContainer_PQ`` virtuals (``ShuntContainer``,
  ``TrafoContainer``, ...), which use ``override final`` or ``override``.
  Harmless (the signatures already matched), but a regression against the
  ``lightsim2grid_core`` target's "these two warnings are clean, so any
  future one fails CI" contract (see ``src/core/CMakeLists.txt``). Verified
  with a clean rebuild (warnings gone) and the SVC / converter-station /
  HVDC unit test suites (22 tests, all green).
- [FIXED] a plain ``import lightsim2grid`` no longer prints "lightsim2grid.solver is
  deprecated...". ``LightSimBackend`` needed ``SolverType`` (for the ``solver_type`` /
  ``SolverType`` back-compat bridging above) and imported it from the deprecated
  ``lightsim2grid.solver`` shim, whose module-level ``DeprecationWarning`` therefore fired
  on every import of ``lightsim2grid`` itself, regardless of whether the deprecated
  ``solver_type`` / ``SolverType`` were ever used. ``SolverType`` now lives in a private
  module (``lightsim2grid/_solver_type.py``); ``lightsim2grid.solver`` still
  re-exports it and still warns when *it* is imported directly.
- [FIXED] a regression from the previous entry: ``SolverType`` was first placed inside
  ``lightsim2grid._utils`` rather than as a standalone module. Importing any submodule of
  a package always runs that package's ``__init__.py`` first, and
  ``lightsim2grid._utils/__init__.py`` does ``from grid2op.Backend import
  PandaPowerBackend`` -- safe only *lazily*, well after grid2op/pandapower have finished
  importing (see the comment at its top: "grid2op -> pandapower -> lightsim2grid ->
  grid2op"). Importing ``SolverType`` eagerly from ``LightSimBackend``'s own top level
  reentered that ``__init__.py`` *while* ``grid2op.Backend.pandaPowerBackend`` was itself
  still mid-import (pandapower's own import chain probes for lightsim2grid), raising
  ``ImportError: cannot import name 'PandaPowerBackend' from partially initialized module
  'grid2op.Backend'`` -- silently swallowed by that file's ``except ImportError: pass``,
  so ``_DoNotUseAnywherePandaPowerBackend`` was never defined and every
  ``grid2op.make(..., backend=LightSimBackend())`` failed outright. Moved to the
  standalone ``lightsim2grid/_solver_type.py`` module (see above), which has no package
  ``__init__.py`` of its own to reenter.
- [FIXED] the documentation now builds with zero Sphinx warnings (was 40, all coming
  from environments without a compiled ``lightsim2grid_cpp`` or without ``numba``, plus
  three real issues): switched ``sphinx.ext.imgmath`` (needs a local LaTeX install) for
  ``sphinx.ext.mathjax`` (client-side, no system dependency); added ``numba`` to the
  ``docs`` extra (``grid2op.Backend.PandaPowerBackend`` warns at import time without it,
  regardless of whether it is used, and autodoc imports it transitively through
  ``LightSimBackend``); fixed a short section-title underline in ``docs/security.rst``;
  and fixed ``LSGrid.change_algorithm``'s C++-side docstring, whose "Examples" RST
  section title broke autodoc once pybind11 concatenated it with the second overload's
  docstring (replaced with ``.. rubric:: Examples``, which does not participate in
  section nesting) -- this docstring also still referenced the pre-1.0 ``DC``/``KLUDC``/
  ``NICSLUDC`` names and the deprecated ``lightsim2grid.solver`` / ``lightsim2grid.gridmodel``
  modules, both now corrected.
- [FIXED] several documentation rendering bugs found by inspecting the built HTML rather
  than just the source: 7 places (``benchmarks.rst``, ``benchmarks_dc.rst``,
  ``benchmarks_grid_sizes.rst``, ``install_from_source.rst``) used Markdown
  ``[text](url)`` links, which RST has no syntax for, so they rendered as literal
  bracketed text -- replaced with plain prose / an actual ``:ref:``;
  ``docs/quickstart.rst``'s docker section had a ``code-block:: bash`` directive glued
  directly to the preceding text line with no blank line, so it never fired and
  "code-block:: bash" rendered verbatim; ``docs/comparison_with_pypowsybl.rst``'s "Open
  Load Flow" external link had its closing backtick before the trailing underscores
  instead of after, so it rendered as inert text instead of a hyperlink. Also removed two
  pieces of dead Sphinx config: ``html_experimental_html5_writer`` (removed from Sphinx
  since version 4) and the ``recommonmark`` extension (no ``.md`` source anywhere in
  ``docs/``, no ``source_suffix`` mapping to make it apply if there were; also dropped
  from the ``docs`` extra in ``pyproject.toml``).
- [FIXED] several dangling Python cross-references, found by temporarily enabling
  Sphinx's ``nitpicky`` mode (off by default, so these silently rendered as plain text
  instead of erroring): ``lightsim2grid.algorithm.SparseLULinearSolver`` doesn't exist in
  Python (it's a C++-only type) so the reference is now plain text;
  ``lightsim2grid.algorithm.TimerJac`` genuinely was never exported from
  ``lightsim2grid.algorithm`` despite being documented as the return type of
  ``get_timers_jacobian()`` -- now exported, mirroring ``LinearSolverStats``; several bare
  (unqualified) references in ``docs/solver_plugin.rst`` and ``docs/security.rst``
  (``change_solver``, ``get_J``, ``LSGrid.check_grid``, ``LSGrid.update_topo``, one more
  found the same way in ``docs/binary_serialization.rst``) now use their fully-qualified
  path, which Sphinx's python domain resolves reliably regardless of which document
  references them; and ``docs/security_analysis.rst``'s ``automodule`` documented the
  deprecated ``lightsim2grid.securityAnalysis`` shim instead of the real
  ``lightsim2grid.contingencyAnalysis`` module, which also cascaded into fixing three
  more dangling references to classes only the real module actually defines
  (``ViolationElementType``, ``LimitViolationType``,
  ``ContingencyAnalysisCPP.is_grid_connected_after_contingency``).
- [FIXED] packaging metadata left over from before the 1.0 release:
  ``requires-python`` said ``>=3.7`` while the README's own compatibility table and
  ``docs/quickstart.rst`` both already said 3.8 is the actual floor; ``build.verbose = true``
  (a debug-only setting, per its own comment) was left on for every user's ``pip install``.
  ``Development Status :: 4 - Beta`` is kept as-is: even at 1.0.0 the package is not
  considered stable yet.
- [FIXED] about two dozen typos across the docs, README and DISCLAIMER, found while
  auditing this branch's documentation. Two were more than cosmetic: `docs/quickstart.rst`,
  `docs/use_with_grid2op.rst` and `docs/lightsimbackend.rst` each had an example that
  would raise ``NameError`` if actually run (looping over the misspelled ``nb_episde``
  instead of the ``nb_episode`` defined two lines above, and calling ``agent.act(...)``
  when the variable actually created was ``my_agent``).
- [FIXED] the remaining "Wrong API" findings from the documentation audit
  (A11-A24), on top of A1-A10 fixed earlier:

  - ``docs/network.rst``: ``total_bus()`` / ``nb_connected_bus()`` were swapped in
    both the prose describing the "solver bus id" convention and the table
    documenting the two methods; ``total_bus()`` is the *total* number of buses,
    ``nb_connected_bus()`` is the number currently seen by the solver.
  - ``docs/use_solver.rst``: ``get_Va()`` / ``get_Vm()`` example comments were
    swapped (``Va`` is the angle, ``Vm`` is the magnitude).
  - ``docs/benchmarks.rst``: the distributed-slack support of ``NR single (SLU)``
    and ``NR (SLU)`` was inverted (the "single" variant is the one that does
    *not* support distributed slack).
  - ``README.md`` / ``docs/install_from_source.rst``: corrected the stale claim
    that ``-O2`` is used by default and ``__O3_OPTIM`` must be set to *enable*
    O3; ``-O3`` has been the default for a while now, and the env var only
    matters to *disable* it (``__O3_OPTIM=0``).
  - ``README.md``: fixed dead source links (``./src/BaseNRSolver.h`` ->
    ``./src/core/powerflow_algorithm/NRAlgo.hpp``, ``KLUSolver.h`` /
    ``SparLUSolver.h`` / ``NICSLU.h`` -> the real ``.hpp`` filenames under
    ``src/core/linear_solvers``); the "Enable CKTSO" section pointed at
    ``PATH_NICSLU`` instead of ``PATH_CKTSO`` in both its prose and its example
    (copy-paste from the NICSLU section above it); the ``compute_pf`` signature
    shown for a custom powerflow solver was the pre-refactor one (same issue as
    A6, just not previously caught here) and the custom-linear-solver interface
    was documented as ``initialize`` / ``solve(J, b, bool)`` / ``reset`` instead
    of the real ``analyze`` / ``factorize`` / ``refactorize`` / ``solve(b)`` /
    ``reset``.
  - ``docs/security_analysis.rst``: the "advanced usage" example pointed at
    ``examples/security_analysis.py`` (renamed to ``examples/contingency_analysis.py``);
    the benchmark reproduction instructions said ``cd examples`` for a script
    that only exists under ``benchmarks/``.
  - ``docs/time_series.rst``: pointed at the old
    ``computers_with_grid2op_multithreading.py`` name, renamed to
    ``timeseries_with_grid2op_multithreading.py``.
  - ``docs/benchmarks_dive.rst``: ``gridmodel.update_topology`` -> the bound
    method is ``update_topo``.
  - ``DISCLAIMER.md`` / ``docs/disclaimer.rst``: dropped the "it does not model
    AC/DC converters" limitation, which is no longer true now that HVDC lines
    with VSC/LCC converter stations are modeled (``elements.HvdcLineContainer``,
    ``elements.ConverterStationInfo``).
  - ``lightsim2grid/solver/__init__.py``: ``FDPF_BX_CKTSOSolver`` (deprecated
    alias) was mapped to ``FDPF_XB_CKTSO`` instead of ``FDPF_BX_CKTSO`` -- a
    copy-paste bug that silently handed anyone still using this deprecated name
    the wrong algorithm variant.
  - ``lightsim2grid/lightSimBackend.py``: ``set_solver_max_iter``'s docstring
    recommended ``AlgorithmType.SparseKLU``, which does not exist (``NR_KLU``);
    ``get_algo_types``'s docstring had a typo (``from ligthsim2grid import
    LightSimBackend``) and its return type annotation was
    ``Union[AlgorithmType, AlgorithmType]`` where it should be
    ``Tuple[AlgorithmType, AlgorithmType]``.
- [FIXED] deprecated names still presented as the current API in prose, examples
  and benchmark scripts (section E of the doc audit):

  - ``SecurityAnalysis`` (renamed to ``ContingencyAnalysis``): dropped from
    ``README.md`` and ``docs/security_analysis.rst`` prose, and from the printed
    output of ``examples/contingency_analysis.py`` and ``benchmarks/security_analysis.py``
    (both of which also had a local variable literally named ``security_analysis``,
    renamed to ``contingency_analysis``).
  - ``gridmodel`` / ``GridModel`` (renamed to ``network`` / ``LSGrid``): fixed in
    ``docs/network.rst``, ``docs/benchmarks.rst``, ``docs/benchmarks_dive.rst`` and
    ``docs/comparison_with_pypowsybl.rst``.
  - ``Computers`` (alias of ``TimeSeriesCPP``): ``docs/time_series.rst`` presented
    it as the low-level class to use; pointed at ``TimeSeriesCPP`` instead.
  - ``gm.get_solver_type()`` (deprecated alias of ``get_algo_type()``):
    ``docs/solver_plugin.rst`` used it in a code example.
  - ``benchmarks/benchmark_ca.py``, ``benchmarks/benchmark_dc_solvers.py``,
    ``benchmarks/benchmark_grid_size.py``: removed the dead
    ``except ImportError: from lightsim2grid import SecurityAnalysis as
    ContingencyAnalysis`` fallback (``SecurityAnalysis`` is not exported from
    the top-level package on this branch, so the fallback would itself raise
    ``ImportError`` if it were ever reached) and the unused/deprecated
    ``from lightsim2grid.solver import SolverType`` import.
  - ``benchmarks/benchmark_grid_size.py``: ``SolverType.KLU`` does not exist,
    not even as a deprecated alias -- this line would raise ``AttributeError``
    the moment the script's ``__main__`` block ran. Replaced with
    ``AlgorithmType.NR_KLU``.
  - ``set_solver_type`` -> ``set_algo_type``, ``change_solver`` -> ``change_algorithm``,
    and ``.available_solvers`` -> ``.available_default_algorithms`` in
    ``benchmarks/benchmark_gauss_seidel.py``, ``benchmarks/benchmark_solvers.py``,
    ``benchmarks/benchmark_dc_solvers.py``, ``benchmarks/benchmark_grid_size.py``,
    ``benchmarks/test_profile.py`` and ``benchmarks/compare_lightsim2grid_pypowsybl.py``
    -- these are the reference scripts people copy, and all of the above still
    *ran* despite using deprecated names.
  - ``benchmarks/test_profile.py`` imported ``AlgorithmType`` from the deprecated
    ``lightsim2grid.solver`` shim instead of ``lightsim2grid.algorithm``.
- [FIXED] the ``compile_gcc_earliest`` CircleCI job (``gcc:8``, a Debian buster
  based image) was failing: it built a venv against that container's default
  ``python3``, which is 3.7, while ``requires-python`` is now ``>=3.8``. Switched
  it to the same ``uv``-based pattern already used by ``compile_clang_earliest``
  and ``test_legacy_grid2op`` in this file (``uv venv venv_test --python 3.8``),
  which pins a supported interpreter regardless of what the container ships by
  default. The container's ``pip``/``python3`` (3.7) turned out to also be too
  old to resolve `uv`'s own PyPI wheels ("Could not find a version that
  satisfies the requirement uv"), so ``uv`` itself is now installed with its
  standalone shell installer (``curl ... | sh``) instead of ``pip install uv``.
  Once that got the build running, it failed a third time, further along:
  ``cc1plus`` was OOM-killed compiling ``binding_lsgrid.cpp`` (the largest
  pybind11 binding translation unit). The compile command showed
  ``-flto=auto -fno-fat-lto-objects``: pybind11's CMake helpers auto-enable
  LTO on any target where the parent project leaves
  ``CMAKE_INTERPROCEDURAL_OPTIMIZATION`` undefined (see
  ``pybind11NewTools.cmake``), and lightsim2grid never sets it. LTO is
  memory-hungry and this job runs on a "medium" (RAM-constrained) resource
  class. Passed ``-C cmake.define.CMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF``
  for this job only -- it exists to check the package still compiles with
  the oldest supported gcc, not to produce an optimized artifact. Verified
  locally that this removes every ``-flto`` flag from the build and that the
  resulting extension still imports and passes tests. Disabling LTO turned out
  not to be enough on its own: the real CircleCI run still OOM-killed
  ``cc1plus`` on the same file, at plain ``-O3`` with no LTO flags at all.
  ``CMAKE_BUILD_TYPE=Release`` (set in ``pyproject.toml``) makes CMake apply
  its own default ``CMAKE_CXX_FLAGS_RELEASE`` (``-O3 -DNDEBUG``) to every
  target project-wide, independent of this project's own ``USE_O3_OPTIM``
  toggle (which only ever applies to the ``lightsim2grid_core`` target, not
  ``lightsim2grid_cpp`` -- the target ``binding_lsgrid.cpp`` belongs to).
  Rather than chase which CMake mechanism to override to strip ``-O3`` here,
  bumped this job's ``resource_class`` from ``medium`` to ``large``.

- [FIXED] stale statements across the docs that were no longer true (section H
  of the doc audit):

  - ``docs/use_with_grid2op.rst``: dropped the "for now it is not easy to
    change [the solver]" note -- ``algo_type`` / ``set_algo_type`` /
    ``change_algorithm`` exist precisely for this, and are now documented.
  - ``docs/use_solver.rst``: solver counts were wrong and CKTSO/FDPF were
    missing from the page entirely. "11 available solvers" -> 22 (across 4
    categories, not 3); the Newton-Raphson list went from 6 to 8 entries
    (added ``NR_CKTSO`` / ``NRSing_CKTSO``); the DC list went from 3 to 4
    (added ``DC_CKTSO``); added a whole missing "Fast Decoupled solvers"
    subsection documenting the 8 ``FDPF_XB_*`` / ``FDPF_BX_*`` variants.
  - ``docs/solvers.rst``: "four families of powerflow algorithms" followed by
    five bullets -> "five families". The example that called
    ``env.backend._grid.change_algorithm(...)`` directly then
    ``env.reset()  # apply the change`` was actively wrong: ``reset()``
    rebuilds ``_grid`` from scratch and re-applies whatever was last set
    through ``set_algo_type``, so a direct ``_grid.change_algorithm()`` call
    is silently reverted on the next reset (verified empirically). Rewrote the
    example to use ``set_algo_type`` and added a warning explaining why.
  - ``docs/algorithm_names.rst``: the "Complete table" omitted
    ``NRRefactorRetry_KLU/NICSLU/CKTSO`` (built-in solvers selectable only by
    string name, with no ``AlgorithmType`` enum value of their own) and
    ``AlgorithmType.Custom``; added both. The usage example was missing
    ``import grid2op`` and, after ``set_algo_type(AlgorithmType.NR_KLU)``,
    showed the output as if ``NRSing_KLU`` had been set instead -- fixed to
    match what the call actually returns (verified empirically).
  - ``docs/security_analysis.rst``: dropped a ``.. warning::`` about a bug in
    lightsim2grid 0.5.5 fixed in 0.6.0 -- irrelevant noise on a 1.0 release.
  - ``docs/quickstart.rst``: dropped the "python 3.6 at time of writing" claim
    about the docker image (impossible to keep accurate, and moot now that
    the package requires >=3.8 anyway); fixed the "Clean-up" section number
    (was "1.", should be "4.", the fourth top-level step after "Install
    docker" / "Get the lightsim2grid image" / "Run a code on this container").
  - ``docs/use_with_grid2op.rst`` and ``docs/lightsimbackend.rst``: replaced
    the long-superseded ``env_name = "rte_case14_realistic"`` with
    ``l2rpn_case14_sandbox``, already used by ``docs/quickstart.rst``.
  - ``docs/lightsimbackend.rst``: fixed a code example that assigned
    ``env_with_iidm_as_the_grid_description = ...`` then called
    ``grid2op.make(env_name, ...)`` -- ``env_name`` was never defined in that
    snippet; now passes the variable that was actually assigned.
  - ``README.md``: "requires grid2op... at least version 0.7.0" immediately
    followed by ``pip install grid2op>=1.6.4`` -- updated the stated minimum
    to match; fixed ``python3 setup.py build`` references (there is no
    ``setup.py``, the build is scikit-build-core driven) and windows-cmd-only
    env var syntax used without a linux/macos equivalent alongside it; fixed
    ``git clone .../grid2op.git`` followed by ``cd Grid2Op`` (wrong case, the
    cloned directory is ``grid2op``).
  - ``docs/install_from_source.rst``: removed ~70 lines of legacy "for
    lightsim2grid < 0.13" SuiteSparse compilation instructions (make / cmake
    option A/B) that are pure noise for a 1.0 release -- SuiteSparse is
    compiled automatically now; removed a stale note about needing to
    ``pip install pybind11`` manually for old versions (build isolation
    handles it via ``pyproject.toml``'s ``[build-system]`` requires); added a
    pointer to :ref:`cpp_library` (which already documents
    ``-DBUILD_TESTING=ON``, scikit-build-core, and the ``lightsim2grid_core``
    shared library) that this page previously never linked to.
  - ``docs/disclaimer.rst``: had drifted from ``DISCLAIMER.md`` -- was missing
    the ``power-grid-model`` entry from the list of alternative tools, linked
    ``github.com/rte-france/grid2op`` instead of ``Grid2Op/grid2op``, and said
    "Eigen (optionally KLU)" instead of "Eigen and KLU". Synced to match.

- [FIXED] structural issues in ``examples/`` and ``benchmarks/`` (section I of
  the doc audit):

  - ``examples/time_serie.py``: ``time_serie.init_from_n_powerflow = True`` was
    set *after* ``time_serie.compute_V(...)`` had already run, so it had no
    effect at all -- the flag is only read inside ``compute_V``'s underlying
    C++ call. Verified empirically: moving it before ``compute_V`` cuts the
    reported ``preprocessing_time()`` roughly 3x on the same scenario. Moved
    the line and documented why the order matters.
  - ``examples/contingency_analysis.py`` and ``benchmarks/security_analysis.py``
    are near-duplicates that measure different things (only the example sets
    ``init_from_n_powerflow``); both had "ignore the protection... by the
    **TimeSerie** module" comments copy-pasted from the time-series example
    (should say ``ContingencyAnalysis``), and the example's own "the 3 lines
    above are the only lines you need" comment had gone stale (4 lines, one of
    them optional) once ``init_from_n_powerflow`` was added. Also removed two
    unused imports (``BaseAction``, ``ChangeNothing``) from the benchmark
    script (confirmed dead with ``pyflakes``).
  - ``docs/security_analysis.rst`` sends readers to the contingency-analysis
    example for "a more advanced usage", but it demonstrated none of the
    capabilities the docs advertise. Added a section exercising
    ``nb_thread``, ``handle_disconnected_grid``, ``compute_limit_violations``
    and ``run_ac``/``run_dc`` -- including working around the fact that
    setting ``compute_limit_violations`` clears any already-registered
    contingency, so the n-1 list has to be re-added afterwards (verified this
    gotcha empirically: without the re-add, the demo silently runs on 0
    contingencies).
  - ``examples/timeseries_with_grid2op.py``: "consult the documentation of
    **TimeSeries**" -- the easier-to-use Python class is ``TimeSerie`` (no
    trailing s).
  - ``examples/Readme.md`` listed the 4 top-level scripts but not the three
    solver-plugin examples (``external_algorithm``, ``dist_slack_algorithm``,
    ``lm_algorithm``) or the shared ``cmake/`` helper; added all four.
  - ``docs/solver_plugin.rst``: fixed the path of ``env_compile_all.sh``
    (repository root, not ``benchmarks/``); fixed a self-contradicting bullet
    that called ``__O3_OPTIM`` one of "two opt-in build flags" and then noted
    "(on by default, actually)" in the same sentence -- it is opt-*out*, only
    ``__COMPILE_MARCHNATIVE`` is opt-in; fixed the expected plugin-example
    output showing a registered solver named ``'DC'``, which does not exist
    (real names are e.g. ``'DC_SparseLU'`` / ``'DC_KLU'``).
  - **Real bug**, not just docs: both the inline CMake template in
    ``docs/solver_plugin.rst`` and the shipped
    ``examples/cmake/MatchLightsim2gridBuildFlags.cmake`` treated an *unset*
    ``__O3_OPTIM`` as OFF when falling back to env-var detection (source-tree
    builds only), while ``lightsim2grid_core`` itself treats unset as ON (only
    ``__O3_OPTIM=0``/``False`` disables it). A plugin built from a source tree
    would silently end up without ``-O3`` while the core it links against has
    it. Fixed both copies to match the core's actual default; verified with a
    standalone CMake project that the "matching -O3" status message now
    prints with the env var unset, and stays silent with ``__O3_OPTIM=0``.

- [FIXED] a **real bug**, found while documenting these classes: ``ContingencyAnalysis``
  and ``TimeSerie`` were defined as ``class __ContingencyAnalysis`` /
  ``class ___TimeSerie`` and then aliased (``ContingencyAnalysis = __ContingencyAnalysis``)
  to their public names. Sphinx autodoc treats a name bound this way as a plain
  alias -- it rendered a one-line "``class ContingencyAnalysis``: alias of
  ``__ContingencyAnalysis``" stub and silently dropped every method's docstring,
  while the real class was itself skipped as "private" (leading underscores).
  In effect, the two main user-facing classes of the whole package had *no*
  API reference at all, on any page, ever -- this is the actual root cause
  behind several "documented nowhere" gaps below. Renamed both classes to their
  public names directly and dropped the alias; verified the full methods /
  properties of both classes (``add_all_n1_contingencies``, ``nb_thread``,
  ``run_ac``, ``compute_V``, etc.) now render on their respective pages, and
  re-ran the full ``ContingencyAnalysis`` / ``TimeSerie`` test suites (32
  tests, all green) plus a direct import/usage smoke test.
- [FIXED] ``ContingencyAnalysis.init_from_n_powerflow`` / ``handle_disconnected_grid``
  and ``TimeSerie.init_from_n_powerflow`` had **no docstring at all** on the
  Python wrapper class (only the underlying C++ ``*CPP`` computer object did),
  so even after the rename above they still would not have appeared under
  ``:members:`` (which skips undocumented members by default). Added real
  docstrings to all three, and a usage note to ``docs/security_analysis.rst`` /
  ``docs/time_series.rst`` explaining that ``init_from_n_powerflow`` must be
  set *before* the computation runs.
- [FIXED] a **real, currently-broken bug** in ``PhysicalLawChecker`` (found while
  verifying its documented example still works): it imported
  ``from grid2op.Environment.Environment import Environment`` and
  ``from grid2op.Environment.MultiMixEnv import MultiMixEnvironment`` -- neither
  submodule path exists on grid2op 1.12 (the actual files are ``environment.py``
  / ``multiMixEnv.py``, lowercase), so construction failed outright with
  ``ModuleNotFoundError`` on any case-sensitive filesystem. Fixed to import both
  names directly from the ``grid2op.Environment`` package (matching the pattern
  already used elsewhere in this codebase, e.g. ``contingencyAnalysis.py`` /
  ``timeSerie.py``). Fixing the import surfaced a second, independent bug:
  ``check_solution`` calls ``backend.update_from_obs(...)``, which (as of
  grid2op 1.11.0) requires ``_load_bus_target`` / ``_gen_bus_target`` / etc to
  have been initialized by ``Backend.load_grid_public(...)`` -- this class was
  still calling the older, lower-level ``Backend.load_grid(...)`` directly,
  leaving those arrays ``None``. Now calls ``load_grid_public`` when available
  (falling back to ``load_grid`` on older grid2op) and ``assert_grid_correct()``
  afterwards. Verified end-to-end: the documented example now runs, and
  checking the actual converged voltage from a real powerflow gives a KCL
  mismatch on the order of ``1e-12`` (as expected); all 8 pre-existing
  ``test_Checker.py`` tests still pass.
- [FIXED / DOCUMENTED] the remaining "missing documentation" findings from the
  doc audit (section G):

  - MATPOWER / PowerModels / PfΔ grid loaders (``init_from_matpower``,
    ``init_from_powermodels``, ``init_from_pf_delta``) were already picked up
    by ``docs/network.rst``'s ``automodule`` (so they did have real API docs),
    but the page's own intro still said "for now the only way is to get it
    from a pandapower grid". Replaced with a table of all five ``init_from_*``
    loaders and what source format each expects.
  - PTDF / LODF (``LSGrid.get_ptdf`` / ``get_ptdf_solver`` / ``get_lodf`` /
    ``get_Bf``): ``docs/benchmarks_dc.rst`` recommends them prominently but
    never explained what they are or linked anywhere. Added a "PTDF / LODF"
    section to ``docs/network.rst`` and cross-linked it from both places.
  - ``ScalingPolicyType`` / ``RefactorPolicyType`` / ``FDPFMethod`` /
    ``AlgoConfig``: bound in C++, never re-exported from
    ``lightsim2grid.algorithm``. Added them to its ``__all__``, and added a
    "Fine-tuning the Newton-Raphson iteration" section to ``docs/solvers.rst``
    covering both the raw-solver setters and ``LSGrid.get_ac_algo_config`` /
    ``set_ac_algo_config`` -- including a real gotcha found while writing the
    example: ``AlgoConfig.int_params`` / ``real_params`` are returned **by
    value**, so ``config.int_params[0] = ...`` silently does nothing; the
    whole list must be reassigned. Verified both code paths directly.
    (``AlgoControl`` / ``PandaPowerConverter`` were left undocumented: the
    former's own C++ bindings mark every method ``"TODO"`` and it is an
    internal implementation detail, and the latter is only ever used
    internally by ``init_from_pandapower``.)
  - ``SubstationContainer`` / ``SubstationInfo``: bound but absent from
    ``lightsim2grid.elements.__all__`` and from ``docs/network.rst``'s
    "Elements modeled" list. Added the export and a "Substations" subsection
    (written by hand rather than via ``automodule``, since two of
    ``SubstationInfo``'s four fields currently share a copy-pasted, wrong C++
    docstring -- a C++-side issue left alone per the standing "C++ docs are
    out of scope for now" agreement).
  - ``examples/lm_algorithm/`` (a Levenberg-Marquardt-damped Newton-Raphson
    solver plugin) was not mentioned anywhere: ``docs/solver_plugin.rst`` said
    "both example plugins" listing only two. Fixed to three, and added an
    "Other example plugins" section describing ``dist_slack_algorithm`` and
    ``lm_algorithm``.
  - ``lightsim2grid.pandapower_compat`` (the actual home of ``newtonpf`` /
    ``dcpf``) was never named; docs and README only ever mention the
    ``lightsim2grid.newtonpf`` re-export shim, and its DC counterpart
    (``pandapower_compat.dcpf``) was undocumented entirely. Added notes to
    ``docs/use_solver.rst`` and ``README.md``.
  - New ``LightSimBackend`` kwargs ``stop_if_storage_disco``,
    ``automatically_disconnect``, ``gen_slack_id`` were absent from
    ``docs/lightsimbackend.rst``; the existing ``stop_if_load_disco`` /
    ``stop_if_gen_disco`` entries also documented the wrong default
    (``Optional[bool] = True``; the real default is ``None``, which now
    defers to grid2op's own ``allow_detachment``) and called the whole
    section "not yet supported by grid2op so not really usable", which is no
    longer true for grid2op >= 1.11.0. Rewrote the section to match the
    actual current behaviour.
  - README.md said nothing about the headline 1.0 features (the
    ``lightsim2grid_core`` C++ library, the algorithm-plugin mechanism,
    binary serialization, multi-threaded contingency analysis, PTDF/LODF, the
    ``gridmodel``/``GridModel`` -> ``network``/``LSGrid`` rename). Added a
    "Key features" section up top, with links to the relevant documentation
    pages, and replaced the "Using a custom powerflow solver" section (a
    stale predecessor of the plugin mechanism, describing embedding a solver
    in lightsim2grid's own source rather than the actual, current
    out-of-tree plugin API) with a short, accurate pointer to
    ``docs/solver_plugin.rst``.
  - Removed now-inaccurate "doc in progress" / "rather incomplete" banners
    from ``docs/network.rst``, ``docs/security_analysis.rst``,
    ``docs/solvers.rst``, ``docs/time_series.rst`` and ``docs/rewards.rst``
    (each has since accumulated a full API reference and worked examples),
    the "``TODO DOC in progress``" line from ``docs/benchmarks.rst`` /
    ``benchmarks_dc.rst`` / ``benchmarks_grid_sizes.rst``, and the two inline
    ``(TODO DOC)`` markers on ``use_solver.rst``'s ``get_timers`` /
    ``get_error`` (replaced with their actual return values).
  - ``bake_outer_loops`` (freezes OLF's converged outer-loop state -- tap
    positions, reactive-limit switches, slack participation, ... -- into a
    pypowsybl network's inputs so it becomes a plain, loop-free power-flow
    problem) already had a full API entry via ``docs/network.rst``'s
    ``automodule``, but on ``docs/comparison_with_pypowsybl.rst`` -- the page
    that actually motivates it -- it was only a passing prose mention, unlike
    its siblings ``compare_baked`` / ``ComparisonResult`` which get a full
    ``autofunction`` / ``autoclass`` treatment there. Gave it the same
    treatment, with a worked example (verified to run) showing the full
    "solve with outer loops -> bake -> solve loop-free" flow. Also documented
    ``get_pypowsybl_loopfree_distributed_slack_parameters`` (loop-free except
    OLF's own ``DistributedSlack``, matching lightsim2grid's default
    distributed slack), which wasn't mentioned on that page at all.

- [DOCUMENTED] ``docs/comparison_with_pypowsybl.rst``'s "Disclaimer" section listed
  several individual gaps (reactive limits, tap ratio, ...) without naming the single
  architectural difference behind all of them: OpenLoadFlow wraps its Newton-Raphson
  solve in "outer loops" (solve, check a criterion, adjust an input, solve again --
  distributed slack, PV<->PQ reactive-limit switching, discrete tap/shunt control, area
  interchange, secondary voltage control, ...), while lightsim2grid's algorithms solve a
  single, fixed problem with no outer-loop mechanism at all. Added a section explaining
  this, and that the two architectures aren't simply "more features vs. fewer": distributed
  slack is the one outer loop lightsim2grid folds directly into the same Newton-Raphson
  Jacobian instead (``MultiSlackNRSystem``, see ``src/core/powerflow_algorithm/NRSystem.hpp``),
  while the others (reactive limits, discrete tap changing, area interchange, secondary
  voltage control) have neither an outer loop nor an in-Jacobian equivalent in lightsim2grid
  today -- which is exactly the gap ``bake_outer_loops`` papers over for comparison purposes.
  Also cross-referenced ``examples/dist_slack_algorithm/``, a solver plugin that
  reimplements distributed slack the *other* way -- as an explicit OLF-style outer loop
  around a single-slack inner solve -- as a concrete demonstration that lightsim2grid's
  plugin mechanism can express an outer-loop-style algorithm at all.

- [ADDED] ``GridModel.check_grid()`` (C++ ``LSGrid::check_grid()``): a whole-grid
  consistency check that verifies every index the grid carries (element bus ids,
  substation ids, position in the topology vector, generator slack and
  remote-regulated bus references) is in range. It raises ``IndexError`` /
  ``RuntimeError`` on an inconsistency.
- [ADDED] the grid is now validated automatically with ``check_grid()`` when it is
  loaded (from a pickle or the binary format, via ``set_state``) and by every grid
  loader (``init_from_pandapower`` / ``init_from_pypowsybl`` / ``init_from_matpower``
  / ``init_from_powermodels``). A well-formed but inconsistent state (e.g. an
  out-of-range bus id in a crafted binary file) now raises a clean exception instead
  of causing an out-of-bounds access during the next powerflow.
- [FIXED] ``sn_mva`` (the base power of the whole per-unit system) and ``init_vm_pu``
  (the flat-start voltage magnitude) were validated **nowhere**: not by their python
  setters, not by ``check_grid()``, and the loaders take them verbatim from their source
  file (``init_from_powermodels`` / ``init_from_matpower`` do
  ``set_sn_mva(float(network["baseMVA"]))``). A degenerate value does not make the
  powerflow fail, it makes it *quietly wrong*: with ``sn_mva = NaN`` the **DC powerflow
  reported convergence and returned NaN branch flows**, and with a negative one it
  returned a plausible-looking but sign-inverted per-unit system. The built-in solvers'
  own finiteness guards cannot catch this -- they check ``Va``, which is finite, since
  ``Bbus`` does not involve ``sn_mva``; the NaN only appears afterwards, when the results
  are scaled back to MW / MVar -- and the size/finiteness check on the solver output is
  deliberately reserved for external solvers. Both scalars must now be finite and
  strictly positive, checked by the setters and re-checked by ``check_grid()`` (which is
  what covers pickle, the binary format and every loader, since ``set_state`` bypasses
  the setters).
- [FIXED] ``PandaPowerConverter::_check_init`` tested ``sn_mva_ <= 0.`` / ``f_hz_ <= 0.``.
  Every comparison with NaN is false, so a NaN sailed through both -- and they are
  divisors in every method of that class, so the whole converted grid came back NaN with
  no error raised anywhere. Both are now checked as finite and strictly positive.
- [FIXED] ``PandaPowerConverter`` (``get_trafo_param_pp2`` / ``get_trafo_param_pp3`` /
  ``get_line_param`` / ``get_line_param_legacy``, all exposed to python) never checked
  that the arrays it is given have the same length -- there was a ``TODO check all
  vectors have the same size`` where the check belonged. Their element count is taken
  from the *first* argument, and the others are then read with unchecked accessors
  (``vect.coeff(i)``, ``is_tap_hv_side[i]``) and combined in coefficient-wise Eigen
  expressions whose result is sized from one of them. Eigen's own size assertions are
  compiled out of the release wheels (``-O3 -DNDEBUG``), so a caller passing a longer
  first array both read *and wrote* past the end of the shorter ones -- reproducible as
  ``free(): invalid next size (fast)``, i.e. heap corruption, from a handful of lines of
  plain python. All the input lengths are now checked up front.
- [FIXED] ``TrafoContainer::set_state`` copied the two per-transformer
  ``alpha -> r/x correction`` tables (the phase-shifter impedance dependency, see
  ``set_shift_dependent_rx``) out of the state with no check at all, and then called
  ``_update_model_coeffs()``, which indexes ``rx_corr_alpha_[el_id]`` /
  ``rx_corr_pct_[el_id]`` for every transformer with an unchecked
  ``std::vector::operator[]``. A pickle or binary file declaring fewer tables than
  transformers therefore read a ``std::vector`` object out of heap memory past the end
  of the array and dereferenced its pointers -- a segfault, or worse. Two tables of
  *different* lengths for the same transformer had the same effect one level down
  (``_shift_rx_corr_pct`` interpolates ``ys`` with indices derived from ``xs.size()``).
  This ran inside ``set_state``, i.e. **before** ``check_grid()``; both shapes are now
  validated there, exactly as ``init()`` / ``set_shift_dependent_rx()`` maintain them.
- [FIXED] ``TwoSidesContainer::set_tsc_state`` (powerlines, transformers, hvdc lines)
  never checked the length of ``status_global_``: ``nb()`` is ``side_1_.nb()``, so a
  pickle or binary file could declare a shorter -- in particular empty --
  ``status_global`` while both sides carried the real element count. That vector is then
  indexed with element ids bounded by ``nb()`` all over the class and the batch solvers
  (``resolve_status``, ``_deactivate``, ``fillYbus``, ``ContingencyAnalysis``...) with an
  unchecked ``operator[]``. Its length must now match exactly.
- [FIXED] ``check_grid()`` did not range-check ``SvcContainer``'s ``regulated_bus_id_``,
  although it is a gridmodel bus id used *directly as an index* by the powerflow --
  ``id_grid_to_solver[regulated_bus_id_(svc_id)]`` in ``SvcContainer::set_vm`` (followed
  by a **write** into ``V``) and in ``LSGrid::fill_voltage_control_solver_data``. Nothing
  bounded it: ``init_svcs`` only checks the vector's length, so an SVC regulating an
  out-of-range bus -- from a crafted pickle / binary file, or straight from a grid file
  through ``init_from_pypowsybl`` -- passed ``check_grid()`` and then read and wrote out
  of bounds on the next ``ac_pf``. It is now validated like the generator field of the
  same name (``-1``, meaning "regulates no bus", stays legal).
- [FIXED] ``check_grid()`` collected the (optional) ``pos_topo_vect`` of *every*
  container -- shunts, static generators, svcs and hvdc lines included -- into the set it
  proves is a permutation of ``[0, dim_topo)``. But ``dim_topo`` there was just how many
  positions it had collected, while ``update_topo()`` sizes its caller arrays as
  ``nb loads + nb gens + nb storages + 2 * nb lines + 2 * nb trafos``, which excludes
  those containers. A state putting positions on a shunt therefore inflated the bound, so
  a *validated* load position could still be past the end of the array ``update_topo()``
  indexes with it. Only the containers ``update_topo()`` actually drives contribute now,
  and a shunt / sgen / svc / hvdc carrying a topology-vector position (there is no setter
  for one: such a state can only come from a crafted file) is rejected.
- [FIXED] ``SubstationContainer::set_state`` performed **no** validation at all, and it
  restores the root of the grid's index space: ``nb_bus()`` (the bound ``check_grid()``
  validates every element bus id against), ``bus_status_`` (the vector those same ids are
  used to index) and ``n_sub_`` / ``nmax_busbar_per_sub_`` were all read from the file
  independently of one another. A pickle or a binary file declaring, say, 4000 buses but a
  1-entry ``bus_status`` passed ``check_grid()`` and then corrupted the heap on the next
  ``init_bus_status()`` (``disconnect_all_buses()`` writes ``nb_bus()`` entries into it).
  All of these are now cross-checked on load, and re-checked by ``check_grid()``.
- [FIXED] ``n_sub_ == 0`` restored from a state reached ``sub_id_of_bus()``, which does
  ``gridmodel_bus_id % n_sub_`` -- an integer division by zero (SIGFPE that kills the
  process), the same class of bug as ``refactor_every_n == 0``. A grid with buses must now
  declare a strictly positive substation count, and ``n_sub_ * nmax_busbar_per_sub_`` is
  rejected if it overflows the ``int`` it is stored in.
- [FIXED] ``LSGrid::set_ls_to_orig`` accepted any value: ``set_ls_to_orig_internal`` sizes
  the reverse mapping from ``lpNorm<Infinity>()`` (the maximum **absolute** value) and then
  indexes it with the values themselves, so an entry of ``-5`` sized the vector from 5 and
  wrote at index ``-5`` -- an out-of-bounds heap write, reachable from the ``_ls_to_orig``
  python property as well as from a pickle / binary file. Entries must now be ``-1`` or a
  sane non-negative original-grid bus id. ``_orig_to_ls`` is range-checked the same way.
- [FIXED] ``_init_kwargs`` is serialized as two parallel vectors whose lengths are stored
  independently; ``set_state`` walked the keys and indexed the values with the same
  counter, so a file declaring more keys than values read past the end of the values
  vector -- constructing ``std::string`` objects from arbitrary heap contents. The two
  lengths must now match.
- [FIXED] iterating ``gridmodel.get_substations()`` crashed on any grid whose substation
  names were never set (``set_substation_names`` is optional and the pandapower / matpower
  / powermodels loaders never call it): ``SubstationInfo`` read ``sub_names_[id]``
  unconditionally, bounded only against the *substation count*. This needed no crafted
  input at all.
- [FIXED] ``update_topo`` did not check the length of the ``has_changed`` / ``new_values``
  arrays it is given. They are indexed **by position in the topology vector** with an
  unchecked Eigen ``operator()``: the positions are validated (``check_grid()`` proves they
  form a permutation of ``[0, dim_topo)``), but a caller-supplied array shorter than
  ``dim_topo`` was simply read past its end. Both must now have exactly ``dim_topo``
  entries.
- [FIXED] ``update_slack_weights`` did not check that its ``could_be_slack`` array had one
  entry per generator, although it indexes it by generator id.
- [FIXED] ``check_grid()`` now also validates the substation container's own internal
  consistency and the bus-id mapping vectors (``_ls_to_orig`` / ``_orig_to_ls`` /
  ``_bus_fusion_rep``), which it previously ignored entirely.
- [FIXED] whether a solver is one of lightsim2grid's own or comes from a plugin was
  decided by ``AlgorithmType``, which cannot answer that question: it is a fixed enum of
  *serialized* solver identities, and a built-in only appears in it if a member was added
  for it. The ``NRRefactorRetry_KLU`` / ``_NICSLU`` / ``_CKTSO`` family never got one, so
  ``name_to_algo_type()`` reported those three **built-in** solvers as
  ``AlgorithmType::Custom`` exactly like a plugin -- and they were consequently paying for
  the external-solver output check (voltage size + finiteness) on *every single solve*,
  which built-in solvers are explicitly supposed to cost nothing. The registry now records
  a ``SolverOrigin`` (``Builtin`` / ``External``) when a solver is registered, which is
  where the answer is actually known; ``AlgorithmSelector::is_builtin_algo()`` caches it so
  the hot path stays a bool read. ``SolverOrigin::External`` is the default, so a solver
  registered without saying anything is treated as untrusted.
- [FIXED] ``SubstationContainer`` did not check that a nominal voltage is a finite,
  strictly positive number; 0, a negative value or NaN silently produced nonsense per-unit
  conversions. NB this deliberately covers ``bus_vn_kv`` / ``sub_vn_kv`` only, never
  ``bus_vmin_kv`` / ``bus_vmax_kv``, which use NaN as the documented "no limit set" sentinel.
- [FIXED] ``n_sub * nmax_busbar_per_sub`` (the total bus count, stored in an ``int`` and
  used in ``int`` index arithmetic) was computed without an overflow check in ``init_bus``
  and ``init_sub``. All three entry points (those two plus ``set_state``) now go through a
  single ``checked_nb_bus`` helper doing the multiplication in 64 bits and refusing a
  product that does not fit, along with the positivity checks.
- [FIXED] ``LSGrid::set_ls_to_orig_internal`` was declared ``noexcept`` although it assigns
  one Eigen vector and allocates another, either of which throws ``std::bad_alloc`` on
  failure -- which in a ``noexcept`` function is an immediate ``std::terminate()`` rather
  than something the caller can handle. Nothing needed the guarantee, so it is gone.
- [FIXED] the "every bus of a substation must have the same nominal voltage" check in
  ``init_bus`` never fired: a gridmodel bus id is ``sub_id + (local_bus_id - 1) * n_sub``
  and ``LocalBusId`` runs ``1..nmax_busbar_per_sub``, but the loop ran local ids
  ``[1, nmax_busbar_per_sub)`` -- so it started by comparing the reference bus with
  *itself* and stopped one busbar short. With the usual ``nmax_busbar_per_sub == 2`` it
  checked nothing at all, and the invariant was silently unenforced on every path. It is
  now checked over local ids ``2..nmax_busbar_per_sub``, and re-checked by
  ``check_grid()`` -- ``set_state`` bypasses ``init_bus``, so a pickle / binary file could
  always carry a grid violating it.
- [FIXED] ``SubstationContainer::init_sub`` wrote ``n_sub`` values into ``sub_vn_kv_``
  (and into ``bus_vn_kv_``) with an unchecked ``operator[]`` without sizing them first.
  Only the two-argument constructor sizes ``sub_vn_kv_``, and nothing uses it: the live
  path is the default constructor followed by ``init_bus()``, which leaves ``sub_vn_kv_``
  empty. Calling ``init_sub()`` on such a container was an out-of-bounds heap write. It is
  currently uncalled and unbound, so this was latent rather than reachable, but it sized
  its destinations wrongly for whoever wired it up next. ``sub_vn_kv`` is also validated
  against ``bus_vn_kv`` by ``check_grid()`` when it is present (it stays optional: it is
  empty on every grid produced by any loader today, and on every already-saved binary
  file including this repo's own format-4 compatibility fixture).
- [FIXED] an ``AlgoConfig`` (part of the serialized grid state) carrying an out-of-range
  ``RefactorPolicyType`` / ``ScalingPolicyType`` is now rejected instead of being cast to
  the enum, silently falling into a ``default:`` branch and round-tripping back out of
  ``get_config()``. ``set_config`` also validates both policies *before* touching any
  member, so a rejected config no longer leaves a half-applied one behind.
- [FIXED] ``LSGrid::set_orig_to_ls`` did not actually build the inverse of the mapping it
  was given: it walked only the first *n* entries (*n* = the number of non-``-1`` ones,
  which are not necessarily at the front) and stored the lightsim bus id where the
  original one belonged. ``_ls_to_orig`` and ``_orig_to_ls`` are now true inverses.
- [BREAKING] solver names are now restricted to ``[A-Za-z_][A-Za-z0-9_.]{0,63}`` (start with
  an ASCII letter or ``_``; then ASCII letters, digits, ``_`` or ``.``; at most 64
  characters). Registering any other name is refused. A solver name is written into every
  serialized grid and is what re-selects the solver on load, so it has to be an identity
  that cannot be spoofed (a non-ASCII homoglyph of a built-in name), cannot inject content
  into an error message or a log (control characters), and is bounded in length. Every
  built-in and every name used by the shipped example plugins already complies; a plugin
  using an exotic name must rename its solver.
- [BREAKING] the binary format version is bumped to 4: the AC / DC algorithm is now stored
  as its registry **name** instead of an ``AlgorithmType`` enum. Files written by an older
  lightsim2grid are rejected with a clear error (they were already only readable by the
  same format version).
- [FIXED] a grid using a solver from a plugin could be saved but never loaded back, and
  could not even be copied: only the ``AlgorithmType`` was stored, and every external
  solver collapses onto ``AlgorithmType::Custom``, which is not a concrete solver. The
  solver name is now stored, so such a grid round-trips (and ``copy()`` works) as long as
  the plugin is loaded.
- [ADDED] when a grid is loaded whose solver is not registered here (a plugin that has not
  been loaded, or a built-in needing an optional KLU / NICSLU / CKTSO backend this build
  lacks), the error now names the solver, says how to obtain it, and lists the solvers
  that *are* available -- instead of an internal message about ``AlgorithmType::Custom``.
- [ADDED] ``LSGrid.load_binary_without_algorithm(path)``: loads the grid data without
  re-selecting the solver it was saved with (keeping the default solvers), so a grid saved
  with an unavailable solver can still be loaded. All other checks are unchanged.
- [ADDED] a sanity check on the voltages returned by an **external** (plugin) solver
  before they are consumed: a wrong-sized result raises, and a non-finite one is
  reported as a non-converged solve instead of propagating ``NaN`` / ``Inf``. Built-in
  solvers skip this check entirely, so they pay nothing for it.
- [ADDED] a "Security" documentation page describing the trust boundaries (pickle,
  solver plugins and grid2op environments are trusted-input-only; the binary format is
  the least dangerous channel, and is validated with ``check_grid`` on load).
- [BREAKING] For plugin developers (C++ side): solver plugins no longer self-register
  from a static ``AlgorithmRegistrar`` constructor firing during ``dlopen``. Instead a
  plugin writes a ``void(ls2g::PluginRegistrar&)`` registration function and exposes it
  with the new ``LS2G_PLUGIN_ENTRY(...)`` macro, which generates the exported
  ``ls2g_register_plugin`` entry point that ``load_algorithm_plugin()`` calls explicitly
  after loading the library. See ``examples/external_algorithm/`` and ``docs/solver_plugin.rst``.
- [FIXED] ``load_algorithm_plugin()`` no longer risks aborting the interpreter: previously
  a registration failure (ABI mismatch, or a solver name already registered) threw a C++
  exception out of a static constructor running inside ``dlopen``, which is uncatchable and
  calls ``std::terminate()``. Registration now runs from an explicitly-called entry point
  wrapped in ``try/catch`` on the C++ side, so every failure surfaces as a catchable Python
  ``RuntimeError``. Plugin registration is also atomic (a plugin exposing several solvers
  either registers all of them or none) and the loader unloads a rejected plugin.
- [BREAKING] ``lightsim2grid.load_algorithm_plugin()`` now raises ``RuntimeError`` (instead
  of ``OSError``) on failure, and rejects loading a plugin whose solver name is already
  registered — which includes loading the same plugin twice.
- [BREAKING] drop python 3.8 support (end of life end back in 2024)
- [BREAKING] `lightsim2grid.SolverType` enum is no more accessible. You can migrate 
  (for now, deprecation pending) to `from lightsim2grid.solver import SolverType` 
  without any other changes. Or you can see the documentation for the (possible)
  new enum names of the `lightsim2grid.algorithm.AlgorithmType` enum.
- [BREAKING] For plugin developers (C++ side): the virtual method ``set_gridmodel`` in
  ``BaseAlgo`` is renamed ``set_lsgrid``, and the protected member ``gridmodel_ptr_``
  is renamed ``lsgrid_ptr_``.
- [BREAKING] For plugin developers (C++ side): ``BaseAlgo::compute_pf`` /
  ``compute_pf_dc`` (and their ``*_with_input_validation`` wrappers) now take
  ``Sbus`` / ``Pbus`` / ``slack_weights`` as ``Eigen::Ref<const CplxVect>`` /
  ``Eigen::Ref<const RealVect>`` instead of ``const CplxVect&`` / ``const RealVect&``,
  matching the ``Eigen::Ref`` convention already used for ``slack_ids`` / ``pv`` / ``pq``.
  This avoids a copy at every solve when the caller already holds a ``Eigen::Ref`` or a
  matrix row/block (e.g. a batch solver's per-timestep ``Sbus``). Any external algorithm
  plugin overriding ``compute_pf`` / ``compute_pf_dc`` (see ``examples/external_algorithm/``)
  needs to update its signature to match.
- [BREAKING] for Newton based solvers, the Jacobian row / columns does not have 
  the same ordering as before, this is because some modularity is being implemented
  at this level to allow for other types of "extensions" (similar to distributed slack)
  such as the handling of HVDC or remote voltage control. Jacobian indices no more
  follow pandapower convention.
- [BREAKING] (cpp only) the linear solvers now supports methods analyze, factorize 
  (was merged into "analyze" before) and the function "refactor" is called 
  "refactorize" for consistency.
- [BREAKING] The "get_timers_jacobian" method of an algorithm now returns 13 elements.
  To avoid any further breaking changes, we advise user to use the `TimerJac` structure
  which is exported to python and contains all the named element instead.
- [BREAKING] The "get_solver_control" method of LSGrid is now removed. Use
  get_ac_algo_controler() for retrieving the AC algo controller or 
  get_dc_algo_controler() for retrieving the DC algo controller
- [DEPRECATION PENDING] some previous solver names (now called algorithm) were 
  rather ambiguous and not very clear. For example, it was not clear that "KLU"
  referenced the Newton Raphson method, with distributed slack variant, that use_buses_for_sub
  the KLU linear solver. For clarity, it is renamed "NR_KLU" now. Old names 
  are still accessible within the `lightsim2grid.solver.SolverType` enum but 
  will be deprecated in the future. Please see the documentation for the full
  migration guide.
- [DEPRECATION PENDING] the reference to "solver" are now replaced with equivalent "algorithm"
  which is clearer and avoid making mistake between linear solver (*eg* SparseLU from 
  Eigen or KLU) and "method to solve powerflow", such as Newton Raphson. 
  This includes the "lightsim2grid.solver" which is not "lightsim2grid.algorithm"
  and most of cpp side names. Old python names should still be usable but will be 
  deprecated in future releases. See documentation for an exhaustive migration guide.
- [DEPRECATION PENDING] the "solver_control" (cpp side and python side) is 
  renamed "algo_controler" and should be accessed with "get_algo_controler()". 
  The old "gridmodel.get_solver_control()" is still accessible for backward 
  compatibility.
- [DEPRECATION PENDING] lightsim_backend.available_solvers() should be replaced with
  lightsim_backend.available_algorithms(). "available_solvers" is still present for 
  backward compatibility.
- [DEPRECATION PENDING] The "GridModel" class has been renamed to "LSGrid" (clearer name)
  The same happened for "lightsim2grid.gridmodel" module which is now "lightsim2grid.network".
  For now, the "lightsim2grid.gridmodel" is still present, with a deprecation warnings
  throw, to ensure backward compatibility.
- [DEPRECATION PENDING] the `DCLineContainer` / `DCLineInfo` classes are renamed
  `HvdcLineContainer` / `HvdcLineInfo` (clearer name, they model HVDC links). The old
  names are still importable from `lightsim2grid.elements` as deprecated aliases.
- [BREAKING] (pickle) the storage units are now stored in a dedicated `StorageContainer`
  (instead of reusing the `LoadContainer`). Grids pickled with a previous version of
  lightsim2grid can no longer be unpickled.
- [BREAKING] (pickle) the HVDC / DC lines are now stored in a dedicated, richer
  `HvdcLineContainer` (with converter stations and droop control). Grids pickled with a
  previous version of lightsim2grid can no longer be unpickled.
- [ADDED] `LSGrid._bus_fusion_rep`: a property (get/set, persisted through pickle and
  ``save_binary``/``load_binary``) set by ``init_from_pypowsybl(...,
  fuse_zero_impedance_branches=True)``, giving for each lightsim2grid bus id the id of
  the "representative" bus it was fused into (identity for a bus not involved in any
  fusion). Used by ``LightsimResultNetwork`` to fix the bug below.
- [FIXED] `LightsimResultNetwork` (`lightsim2grid.network.LightsimResultNetwork`, built
  on a grid loaded with `fuse_zero_impedance_branches=True`) reported `v_mag=0` for a
  bus fused away by a (near-)zero-impedance branch: that bus's own lightsim2grid bus
  ends up with no elements after fusion (they were all repointed to the fusion
  representative), hence disconnected with no solved voltage of its own. It is now
  reported as its representative's voltage (the two are electrically the same node).
  The fused branch itself similarly read 0 flow / disconnected; its flow is now
  reconstructed via Kirchhoff's current law wherever an endpoint bus has exactly one
  fused branch attached (see `_reconstruct_fused_branches` in `_result_network.py`);
  a bus where 2+ fused branches meet is left as-is, the split being genuinely
  indeterminate from the solved state alone. Found on a real RTE grid
  (`PtFige-20241018-0355`, see the `[TODO]` entry at the top of this file for the
  follow-up C++ work this points at).
- [FIXED] `LSGrid.id_me_to_ac_solver()` returned the AC *solver -> gridmodel* mapping
  (a duplicate of `id_ac_solver_to_me()`) instead of the *gridmodel -> AC solver*
  mapping. It now returns the correct direction (the DC counterparts were already fine).
- [FIXED] an issue preventing to use the `init_from_n_powerflow` attribute of `TimeSeries`
- [FIXED] a warning when compiling lightsim2grid_core (redefinition of some classes)
- [FIXED] storage units (batteries) read from a pypowsybl grid in `init_from_pypowsybl`
  now use the correct sign convention: IIDM batteries declare their `target_p` / `target_q`
  in the generator convention (positive = produced) whereas lightsim2grid stores storage
  in the load convention (positive = charging), so the setpoints are now negated. The
  previous (untested) behaviour modeled a producing battery as a consuming load.
- [FIXED] `ContingencyAnalysisCPP.is_grid_connected_after_contingency()` could segfault: it
  built its working matrix from the grid model's own (never populated, hence empty) Ybus,
  causing out-of-bounds writes. It now uses the correctly indexed internal `Ybus_` and builds
  the required inputs on demand, so it works both before and after `compute()`.
- [FIXED] several Python-exposed methods could read or write out of bounds (potential
  segfault / memory corruption) when given an out-of-range or negative element id or a
  mis-shaped array, because the id / array reached unchecked ``operator[]`` /
  Eigen ``operator()`` / ``.col()`` (whose bounds asserts are compiled out in release
  builds) with no prior validation. Out-of-range ids now raise ``IndexError`` and
  mis-shaped arrays raise ``RuntimeError``, instead of corrupting memory. Affected:
  ``deactivate_* / reactivate_* / change_bus_*`` (loads, gens, sgens, storages, shunts,
  powerlines, trafos, dclines — the bounds check ran *after* the access),
  ``update_slack_weights_by_id`` and ``set_gen_regulated_bus`` (out-of-bounds writes, the
  latter also stored an unvalidated bus id), ``change_ratio_trafo`` /
  ``change_shift_trafo`` / ``change_shift_trafo_deg`` (out-of-bounds writes),
  ``get_status_droop_hvdc``, the grid2op fast-update path
  (``update_gens_p`` / ``update_loads_p`` / ... when ``new_values`` is shorter than
  ``has_changed``), ``set_*_pos_topo_vect`` / ``set_*_to_subid``, and
  ``TimeSeriesCPP.compute_Vs`` (injection-matrix column / row counts are now validated
  against the grid element counts).
- [FIXED] pickling (or `save_binary`/`load_binary`) an `LSGrid` whose `_ls_to_orig`
  bus-mapping was set (which `init_from_pypowsybl` and `init_from_pandapower` both do)
  always raised ``"Impossible to set the converter ls_to_orig: the provided vector has
  not the same size as the number of bus on the grid."`` on restore. `LSGrid::set_state`
  validated `_ls_to_orig` against `substations_.nb_bus()` before `substations_` itself
  had been restored on the fresh instance pickle/binary loading constructs, so that size
  was always 0. Fixed by restoring `substations_` first.
- [FIXED] a member-initialization-order bug in `ContingencyAnalysis`'s constructor
  (compiler ``-Wreorder`` warning): ``_compute_limit_violations_`` was listed before
  ``_li_defaults`` in the initializer list although declared after it in the class.
  C++ always initializes members in declaration order regardless of initializer-list
  order, so the list was misleading (harmless today since there is no data dependency
  between the two, but a latent hazard for future refactors). The initializer list now
  matches declaration order.
- [BREAKING] `ContingencyAnalysisCPP.compute()` (and the python `ContingencyAnalysis.run` /
  `run_ac` / `run_dc` / `compute_V`) now raises a ``RuntimeError`` if the pre-contingency ("n",
  no disconnection) powerflow itself does not converge, instead of silently returning with all
  voltages at 0 and no contingency simulated at all. Every contingency is solved relative to
  this base case, so a diverging base case made the whole analysis meaningless; failing loudly
  is safer than silently returning an all-zero result that could be mistaken for "every
  contingency was skipped".
- [ADDED] `ViolationElementType.GRID` and two new `LimitViolationType` values,
  `NOT_SIMULATED` and `DIVERGENCE`, so a non-converged contingency (`converged()` /
  `ContingencyResult.converged` is `False`) now reports exactly one `LimitViolation` in
  `get_violations()` / `limit_violations` instead of an empty list: `NOT_SIMULATED` when a
  pre-check (graph connectivity) skips the contingency without ever invoking the solver (eg it
  splits the grid in multiple connected components), `DIVERGENCE` when the solver is invoked but
  does not converge (eg reaches `max_iter`). This only applies to `get_violations()` /
  `get_violations_n()` -- see the `ContingencyAnalysisCPP.compute()` change above for the
  pre-contingency ("n") case itself.
- [ADDED] `LimitViolation.name`: the name of the LINE / TRAFO element the violation was
  detected on, if names were configured on the grid (`LSGrid.set_line_names` /
  `set_trafo_names`); empty string otherwise, and for `BUS` / `GRID` violations (there is no
  per-bus name in `LSGrid`, only per-substation ones via `get_substation_names()`). Also adds
  `LSGrid.get_line_names()` / `get_trafo_names()` getters (previously write-only) and
  `GenericContainer::get_names()` c++ side.
- [FIXED] `ContingencyAnalysisCPP` / `TimeSeriesCPP` (and so the python `ContingencyAnalysis`,
  `SecurityAnalysis` and `TimeSerie` wrappers) solved every powerflow with a fresh, independent,
  *default* solver (`NR_SparseLU`, no scaling/damping policy), completely ignoring whatever
  algorithm type or `AlgoConfig` (eg a `ScalingPolicyType.MaxVoltageChange` damping) was set on
  the source `LSGrid`. This is the same class of bug as the `LSGrid` copy-constructor fix above,
  but in `BaseBatchSolverSynch`: a grid whose own `ac_pf()` converged fine (thanks to a configured
  damping policy) could see its pre-contingency ("n") powerflow diverge inside
  `ContingencyAnalysisCPP.compute()`, now raising the `RuntimeError` added above instead of being
  silently swallowed. The algorithm type and config are now copied from the grid model once at
  construction time; `get_algo_config()` / `set_algo_config()` were added to
  `ContingencyAnalysisCPP` / `TimeSeriesCPP` to re-apply a config afterwards (eg after
  `change_algorithm()`, which resets to that algorithm's defaults).
- [FIXED] the pre-contingency ("n") powerflow inside `ContingencyAnalysisCPP.compute()` and the
  initial powerflow inside `TimeSeriesCPP.compute_Vs()` (both go through the shared
  `BaseBatchSolverSynch::_finish_preprocessing`) compared the NR mismatch against the raw,
  un-scaled `tol` argument, instead of `tol / sn_mva` like every other powerflow in the codebase
  (`LSGrid::ac_pf`, and even the *per-contingency* / per-step solves in these same two classes):
  `Sbus_` is already per-unit (`sn_mva`-divided, same convention as `LSGrid::ac_pf`'s `acSbus_`),
  so comparing it against a physical-MW `tol` accepted a mismatch up to `sn_mva` times (eg 100x)
  looser than what the caller asked for on this one call.
- [FIXED] every powerflow solved by `ContingencyAnalysisCPP` and `TimeSeriesCPP` (the
  pre-contingency / initial "n" solve in `BaseBatchSolverSynch::_finish_preprocessing` and
  `::warmup_solver`, every per-contingency solve in `ContingencyAnalysis::run_contingency_range`,
  and every per-step solve in `TimeSeries::compute_Vs`) passed `slack_ids_me_` -- the slack bus ids
  in *gridmodel* (global) numbering -- to `AlgorithmSelector::compute_pf` / `compute_pf_dc`, which
  expect *solver-space* numbering (matching `Ybus_` / `Sbus_` / `bus_pv_` / `bus_pq_`'s own
  indexing), exactly like `LSGrid::ac_pf` correctly passes its `slack_bus_id_ac_solver_`. The
  correct member, `slack_ids_solver_`, was computed by `pre_process_solver` but never actually
  used anywhere. On a grid with any deactivated bus, `id_me_to_solver` compacts the numbering, so
  global bus ids run higher than (and do not correspond to) valid solver-space indices -- silently
  picking the wrong buses as the NR angle reference, up to reading out of bounds on a large grid
  (reproduced as a segfault while isolating this bug). This was invisible on small, fully-connected
  test grids (global id == solver id there, no compaction), which is why the existing test suite
  never caught it. This was the actual cause of the pre-contingency `SolverFactor` divergence
  reported and investigated above: on a real large RTE grid, `ContingencyAnalysisCPP.compute()`'s
  "n" solve diverged while `LSGrid.ac_pf()` called directly on an extracted, verified-bit-identical
  copy (same `Vinit` / `Ybus_` / `Sbus_` / `bus_pv_` / `bus_pq_` / algorithm / `AlgoConfig` / `tol`)
  converged fine; fixing the slack ids resolves it (confirmed convergent, all contingencies, on two
  different real grids).
- [IMPROVED] (cpp) the codebase now compiles warning-free under ``-Wall -Wextra
  -Werror``: fixed 29 signed/unsigned comparison warnings (``int`` / ``Eigen::Index``
  vs ``size_t``) and silenced ~188 intentionally-unused parameters on interface /
  default-virtual-method signatures (kept as ``/*name*/`` comments for documentation,
  matching the convention already used elsewhere in the codebase).
- [IMPROVED] (cpp) the built-in linear solvers now manage their solver memory with RAII
  smart pointers instead of manual allocation. ``CKTSOLinearSolver`` and
  ``NICSLULinearSolver`` hold their CSC index buffers (``ai_`` / ``ap_``) in
  ``std::unique_ptr<T[]>`` allocated with ``std::make_unique`` (was raw ``new[]`` /
  ``delete[]``); ``KLULinearSolver`` holds its ``klu_symbolic`` / ``klu_numeric`` handles in
  ``std::unique_ptr`` with custom deleters that call ``klu_free_symbolic`` /
  ``klu_free_numeric`` (was manual frees in the destructor and ``reset()``). Lifetimes are
  unchanged, but ownership is now RAII-managed, which also removes a latent leak if
  ``analyze()`` / ``factorize()`` were called twice without an intervening ``reset()``.
- [ADDED] `lightsim2grid.network.bake_outer_loops`: rewrites a pypowsybl network's
  input setpoints to the converged PowSyBl OpenLoadFlow (OLF) outer-loop state (tap /
  shunt positions, reactive-limit PV->PQ switches, distributed-slack active power) so
  a subsequent *outer-loop-free* solve in OLF or lightsim2grid reproduces the OLF
  with-loops operating point. Tested in `test_olf_bake` (pure-OLF round trips plus
  lightsim2grid agreement through line / transformer outages).
- [ADDED] `lightsim2grid.network.get_pypowsybl_loopfree_parameters`: a factory for the
  canonical OLF `Parameters` with every outer loop disabled (empty `outerLoopNames`
  allow-list, every trigger forced off). Now the single source of truth used by the
  benchmarks (`compare_lightsim2grid_pypowsybl`), the documentation and the test suite.
- [ADDED] `lightsim2grid.network.compare_baked` (and `ComparisonResult`): a thin helper
  to validate lightsim2grid against OLF on identical inputs (bake, optional outages,
  loop-free solve in both engines, voltage comparison).
- [ADDED] a dedicated `StorageContainer` / `StorageInfo` (exposed through
  `lightsim2grid.elements`) for the storage units, with its convention documented.
- [ADDED] reading the storage units (batteries) from a pypowsybl grid is now tested
  (parity against PowSyBl OpenLoadFlow, see `test_storage_pypowsybl`).
- [ADDED] proper **HVDC** support inside the AC (Newton-Raphson) and DC powerflow:
  HVDC links now participate in the Jacobian, with `ConverterStationInfo` exposing the
  VSC / LCC converter stations and an angle-droop (AC emulation) regime selectable via
  `LSGrid.set_status_droop_hvdc` / `LSGrid.get_status_droop_hvdc`. Tested in DC, in
  batch and for pickling (see `test_hvdc_dc`, `test_hvdc_droop`, `test_hvdc_batch`,
  `test_hvdc_converter_stations`, `test_hvdc_pickle`, `test_hvdc_no_hvdc_bit_identical`).
- [ADDED] reading HVDC lines (and their converter stations) from a pypowsybl grid
  (`init_from_pypowsybl`), tested for parity (see `test_hvdc_pypowsybl`).
- [FIXED] reading a transformer from pypowsybl (`init_from_pypowsybl`) ignored the
  **tap-changer step r / x / g / b corrections** (the per-step percentage deltas):
  only `rho` / `alpha` were folded in, while r/x/g/b were left at their neutral-tap
  value. For phase-shifting transformers whose series impedance varies with the tap
  (common on RTE grids) the through-flow was wrong by tens of MW vs PowSyBl Open Load
  Flow. lightsim2grid has no "tap" concept, so this ``alpha -> r/x correction``
  dependency is stored per transformer and the series impedance is refreshed from the
  **current phase shift** (interpolated) whenever the coefficients are rebuilt -- so an
  *in-place* ``change_shift_trafo`` / ``change_ratio_trafo`` keeps r/x correct too
  (verified to match OLF after a tap change). Enabled only when reading from pypowsybl
  (a flag, off for pandapower). See `LSGrid.set_trafo_shift_dependent_rx`,
  `test_pst_tap_impedance` and `HVDC_OLF_FINDINGS.md`.
- [FIXED] `consider_only_main_component` (and `init_from_pypowsybl(..., only_main_component=True)`,
  the default) used to deactivate an **entire HVDC line** as soon as one of its two
  converters fell outside the main component. For a cross-border / asynchronous HVDC
  link (the two converters are in different *synchronous* components) this silently
  dropped the in-main converter's scheduled injection -- hundreds of MW on real RTE
  grids -- making lightsim2grid disagree with PowSyBl Open Load Flow. The in-main
  converter is now kept injecting (like OLF's boundary injection) and only the
  out-of-main converter is opened; a line with both converters outside is still fully
  dropped. See `test_hvdc_main_component` and `HVDC_OLF_FINDINGS.md`.
- [ADDED] **Static Var Compensators (SVC)**: a dedicated `SvcContainer` / `SvcInfo`
  (exposed through `lightsim2grid.elements`) supporting OFF / VOLTAGE / REACTIVE_POWER
  regulation modes (with an optional voltage / reactive slope), together with the LSGrid
  management methods (`get_svcs`, `deactivate_svc` / `reactivate_svc`, `change_bus_svc`,
  `set_svc_names`).
- [ADDED] **remote voltage control** for generators (and SVCs): a controller can now
  regulate the voltage of a *different* bus through `LSGrid.set_gen_regulated_bus(...)`;
  it is resolved automatically when importing from pypowsybl. Tested against PowSyBl
  OpenLoadFlow (see `test_voltage_control_svc`, `test_voltage_control_remote_gen`,
  `test_voltage_control_pypowsybl`).

  .. warning::
      When importing from pypowsybl, the regulated bus is resolved **once**, at import
      time, and stored by its (fixed) lightsim2grid global bus id. If the regulated
      element is later moved to another bus *inside lightsim2grid* (e.g. through a
      ``change_bus_*`` / topology change), the controller keeps regulating the bus
      resolved at import: the lightsim2grid grid and the original pypowsybl grid then
      desynchronise. Re-import the grid (or call ``set_gen_regulated_bus`` again) to
      follow such a change. This is a known limitation (see the ``TODO`` in the code).
- [ADDED] `LSGrid.get_storages()`, `LSGrid.get_dclines()`, `LSGrid.get_svcs()` and
  `LSGrid.get_voltage_levels()` accessors for the (new) containers.
- [ADDED] multi-threaded contingency analysis: a `nb_thread` attribute on `ContingencyAnalysis`
  and `ContingencyAnalysisCPP` (default ``1``). When set to a value greater than ``1``, the
  contingency list is split into contiguous ranges, each solved by its own OS thread (each
  thread gets its own solver and its own copy of the admittance matrix, and writes to disjoint
  rows of the result matrix). The results do not depend on the number of threads (they match the
  sequential ones up to the solver's convergence tolerance: each thread keeps its own solver
  warm-start state, so converged voltages agree to ~1e-13, far below the powerflow tolerance). It
  is implemented using only the C++ standard library (`std::thread`, no MPI / OpenMP), works for the
  AC (Newton-Raphson), DC and `handle_disconnected_grid` modes, and behaves identically to the
  previous sequential code when `nb_thread == 1`.
- [ADDED] a `handle_disconnected_grid` mode for the contingency analysis (`ContingencyAnalysis`
  and `ContingencyAnalysisCPP`). By default (``False``) a contingency that splits the grid in
  several connected components is skipped (its voltages stay at 0, legacy behaviour). When set
  to ``True``, the largest connected component is solved while the buses of the other
  component(s) are "masked" (their Newton-Raphson equations are forced to identity and their
  voltage reported as 0). This is done **without any extra symbolic factorization** (the
  solver's analyze/factor are not re-triggered): the masked buses keep the Jacobian sparsity
  unchanged. The reference slack is chosen once, up-front, to minimise the number of
  contingencies that have to be skipped (those that would strand the chosen reference), and a
  stranded slack generator has its slack weight zeroed and the remaining weights rescaled.
  Supported both in AC (Newton-Raphson) and in DC: in DC the masked buses' rows of the reduced
  system are forced to the identity (angle 0), the masked injections are dropped and the slack
  imbalance is computed on the live component only, again **without re-triggering the symbolic
  factorization**. An AC non-NR solver (Gauss-Seidel / Fast-Decoupled) raises a clear error.
  Tested in `test_ContingencyAnalysis_split`.
- [IMPROVED] `ContingencyAnalysisCPP.is_grid_connected_after_contingency` now also reports the
  splitting contingencies in DC (it used to always claim "connected"), by labelling the connected
  components of the real `Bbus`.
- [ADDED] (to interpret the new Jacobian layout) the methods `get_theta_to_J_col()`,
  `get_vm_to_J_col()` and `get_q_to_J_col()` on the Newton-Raphson algorithms (and on
  the `AlgorithmSelector`). Each returns a vector, **indexed by the solver bus id**,
  giving the Jacobian column holding that bus' voltage-angle / voltage-magnitude /
  reactive unknown (-1 when the bus owns no such unknown).
- [ADDED] Refactored `ChooseAlgorithm` (used to be ChooseSolver) to a plugin-friendly
 `AlgorithmRegistry` (see doc)
- [ADDED] installing the python package now also comes with the lightsim2grid_core 
  header files.
- [IMPROVED] Removing the "ChooseAlgorithm" API and replacing it with the `AlgorithmRegistry`
- [IMPROVED] clean separation between lightsim2grid_core the main library that could be
  used from cpp and the python bindings.
- [IMPROVED] removed the ".values" and replace them by ".to_numpy()" in pandapower converter.
- [IMPROVED] remove the "typedef" in favor of "using" cpp side (core)
- [IMPROVED] add some "override" and "final" in the algorithm virtual methods.
- [IMPROVED] file names of some example scripts.
- [IMPROVED] documentation: a new "Bus labelling conventions" section (in the LSGrid
  page) explaining the Local / GridModel (global) / Solver bus id conventions, with a
  table of which `change_bus_*` / `get_bus_*` / `update_topo` / `*_solver` method uses
  which convention.
- [IMPROVED] names in the enum for the solver (*eg* AlgorithmType.NR_KLU) now
  matches the names of the solver (*eg* NR_KLU). There is no more differences 
  between the enum and the solver names.
- [IMPROVED] results of get_timers_jacobian() function now measures 13 different timers
  and has is accessible with names (rather than just id). It is still iterable for
  backward compatibility. 
- [IMPROVED] computation speed of DC algorithm in general
- [IMPROVED] dc powerflow is now really separate from AC (different problem, different hypothesis etc.)
- [IMPROVED] the algo controler (used to be solver_control) is now split: there is one
  for DC and one for AC.
- [IMPROVED] a dedicated Jacobian-tester suite (`test_jacobian.py`) now validates the
  refactored Newton-Raphson Jacobian (single and distributed slack, including the HVDC
  and voltage-control extensions).
- [IMPROVED] (cpp) the Newton-Raphson Jacobian layout is now driven by a central
  `NRLedger` registry that assigns each equation (row) and each unknown (column) and
  keeps the bus <-> column mapping introspectable. Columns are now sorted by (solver)
  bus id, which is what changed the ordering w.r.t. the pandapower convention and what
  enables the `get_theta_to_J_col()` / `get_vm_to_J_col()` / `get_q_to_J_col()` mappings
  above. This central registry is also what makes the HVDC / voltage-control Jacobian
  extensions possible.
- [ADDED] a fast, additive binary serialization path: `LSGrid.save_binary(path)` /
  `LSGrid.load_binary(path)`, and the same two methods on every individual element
  container that has its own state (*eg* `gridmodel.get_loads().save_binary(...)`),
  mirroring what is already picklable. This is **not** a replacement for pickle: it
  trades portability for speed, meant for repeatedly re-loading the *same* grid on the
  *same* machine / lightsim2grid build. See the new "Fast binary serialization"
  documentation page and `benchmarks/benchmark_binary_serialization.py` for speed
  comparisons against pickle (the speed up grows with grid size: up to ~17x faster to
  write and ~8x faster to read than pickle on grids with ~9000 buses).
- [IMPROVED] robustness of the powerflow entry points against ill-formed input:
  `LSGrid.ac_pf` / `LSGrid.dc_pf` now validate `max_iter` (>= 0) and `tol`
  (finite and > 0) in addition to the size of `Vinit`, and the python-facing
  `Solver.compute_pf` / `Solver.solve` (Newton-Raphson, single-slack,
  fast-decoupled, Gauss-Seidel -- shared `BaseAlgo::check_pf_inputs`, bound
  via the new `BaseAlgo::compute_pf_with_input_validation`) now validate
  their inputs before touching them: non-square `Ybus`, `V` / `Sbus` /
  `slack_weights` not matching the size of `Ybus`, empty `slack_ids`,
  a negative `max_iter` and a non-positive / non-finite `tol` raise
  `RuntimeError` (`max_iter=0` is accepted: every solver builds its initial
  state / first Jacobian before iterating, so it is a well-defined "return
  the pre-iteration state" call, used by internal tests);
  out-of-range or negative bus ids in `pv` / `pq` / `slack_ids` raise
  `IndexError` (previously they reached raw Eigen indexing: out-of-bounds
  reads/writes in Release builds); a bus listed in more than one of
  slack/pv/pq (or twice in the same one) raises `RuntimeError` instead of
  silently producing a wrong system. This validation is a python-boundary
  concern only: the internal C++ solve path used by `LSGrid` and the batch
  solvers (`ContingencyAnalysis`, `TimeSerie`, security analysis) still calls
  the raw, unchecked `compute_pf` / `compute_pf_dc` in its hot loop, since
  those callers build `Ybus` / `Sbus` / `pv` / `pq` themselves and paying an
  O(n) validation pass on every contingency / timestep would be pure
  overhead. Tested by the new `lightsim2grid/tests/test_pf_input_robustness.py`
  for every solver class available in the build.
- [IMPROVED] robustness of `save_binary` / `load_binary` against ill-formed input:

  - version compatibility is now decided by a dedicated **binary format version**
    (`BINARY_FORMAT_VERSION` in `src/core/BinaryArchive.hpp`, stored in the file
    header) instead of requiring the exact same lightsim2grid version: the format
    number is only bumped when the serialized layout changes, so all lightsim2grid
    versions sharing the same format number read each other's files. An unsupported
    format raises a `RuntimeError` naming both versions and format numbers. A
    committed reference file (`lightsim2grid/tests/binary_format_fixture/`) guards
    against layout changes made without bumping the format number.
  - every count / length field read from a file is now validated against the real
    file size *before* any allocation: a corrupted count raises a clean
    `RuntimeError` instead of a `MemoryError` (or worse, an actual multi-gigabyte
    allocation).
  - the file header now stores the **object type** (*eg* `"LoadContainer"`):
    loading a file into the wrong class raises a `RuntimeError` instead of silently
    succeeding when the layouts happen to match (*eg* `LoadContainer` vs
    `StorageContainer`).
  - trailing bytes after the end of the data are now rejected as corruption.
  - `save_binary` is now **atomic by default**: it writes to a temporary file
    that only replaces the destination once fully written, so an interrupted
    save (crash, full disk, ...) never destroys a previously saved file. Pass
    `atomic=False` for the marginally faster direct write without that
    protection.
  - the integer values of the serialized enums (`AlgorithmType`,
    `SvcContainer.RegulationMode`, `HvdcLineContainer.ConvertersMode`,
    `ConverterStationInfo.ConverterType` -- the last three are now exposed to
    python) are pinned by a test (`TestSerializedEnumValues`) that fails with
    a "bump BINARY_FORMAT_VERSION" message if they are renumbered.
- [ADDED] memory-hardening CI (`.github/workflows/sanitizers.yml`): an
  **ASan + UBSan** job (build with `__SANITIZE=1`, runs the binary
  serialization, solver control, time series and contingency analysis test
  modules under AddressSanitizer + UndefinedBehaviorSanitizer) and a
  **debug-assertions** job (build with `__DEBUG_ASSERTS=1`:
  `-UNDEBUG -D_GLIBCXX_ASSERTIONS`, re-enabling the Eigen bounds assertions
  that Release builds silence, runs the solver-heavy modules). Also a new
  corruption-sweep test (`TestCorruptionSweep`) that corrupts a valid binary
  file at every byte offset and checks `load_binary` never does anything
  worse than raising a clean `RuntimeError`.
- [ADDED] a dedicated C++ unit test suite (Catch2, new git submodule, under
  `src/tests/`) exercising the binary serialization layer (`BinaryArchive`)
  without python or a real grid: synthetic `StateRes` round trips covering
  every serialized field shape, every bounds-check / header-mismatch path,
  the atomic temp-file commit/rollback, and a C++ port of the corruption
  sweep. Built standalone (`cmake -S src/tests`) or via `BUILD_TESTING=ON`,
  and run in CI both through ctest and under `valgrind --error-exitcode=1`
  (`.github/workflows/cpp_unit_tests.yml`) -- practical only because the
  suite is a small plain binary. This is also the first framework for C++
  unit tests of the core (eg future solver-level tests).
- [ADDED] C++ unit tests for the `LSGrid` main API (`src/tests/test_lsgrid.cpp`):
  a 3-bus grid built programmatically through the `init_*` methods and solved
  with the default Eigen SparseLU algorithms -- AC/DC powerflow contract
  (converged => per-bus V, diverged => empty vector), physically-checked
  results (power balance, analytic DC angles), copy / `get_state` /
  `save_binary` round trips, setpoint changes, load deactivation and the
  documented error paths. Also covers every other element type and control
  scenario: shunts, storage units, SVCs (all three regulation modes), HVDC
  (VSC-VSC with and without angle droop, voltage-regulating VSC, LCC power
  factor, LCC+droop rejection), transformers (tap ratio and phase shifter,
  incl. the `change_ratio_trafo` / `change_shift_trafo` setters), distributed
  slack, remote generator voltage control and the rejection of an unfeasible
  local+remote controller pair on one bus. The test target now links
  `lightsim2grid_core`.
- [ADDED] documentation for using `lightsim2grid_core` as a standalone C++
  library (`docs/cpp_library.rst`): building/installing it from source
  (`cmake -S src/core`), consuming the copy shipped inside the python wheel
  (`lightsim2grid.get_cmake_dir()`), linking with CMake via
  `find_package(lightsim2grid_core CONFIG)` / `lightsim2grid::core`, a
  complete build-a-grid-and-solve example, and how to run the C++ unit tests.
- [FIXED] `TrafoContainer` left two bool members (`ignore_tap_side_for_shift_`,
  `shift_dependent_rx_`) uninitialized when `init_trafo` was never called (any
  grid built without trafos): copying or serializing such a grid read
  indeterminate bools (undefined behavior, garbage written into binary files /
  pickles). Found by valgrind over the new C++ LSGrid tests; both members now
  have default initializers.
- [FIXED] `LSGrid.save_binary`/`load_binary` (and pickle, which shares the same
  `LSGrid::get_state()`/`set_state()`/`StateRes` contract) silently dropped the
  per-solver `AlgoConfig` (scaling/refactor policy, line-search tolerances, etc. --
  see `get_ac_algo_config()`/`set_ac_algo_config()` and the DC counterparts): only the
  coarse `AlgorithmType` enum (eg `NR_SparseLU`) round-tripped, so any custom tuning
  applied via `set_ac_algo_config()`/`set_dc_algo_config()` reverted to the solver's
  defaults after a save/load or pickle round trip. Both AC and DC `AlgoConfig` are now
  part of `LSGrid::StateRes` and restored after the algorithm itself is selected on
  `set_state()` (mirroring the existing copy-constructor behavior). `LSGrid::StateRes`
  tuple positions are now named (`LSGrid::SUBSTATION_ID`, `LSGrid::HVDC_ID`, etc.)
  instead of raw integer literals in `get_state()`/`set_state()`.
- [ADDED] `init_from_pypowsybl` now reads and exposes operating limits: per-bus min/max
  voltage (`LSGrid.get_bus_vmin_kv()` / `get_bus_vmax_kv()`, from the source voltage
  levels' `low_voltage_limit` / `high_voltage_limit`) and per-side branch thermal
  (current) limits (`limit_a1_ka` / `limit_a2_ka` on `LineInfo` / `TrafoInfo`, from
  `network.get_operational_limits()`, honoring each side's *selected* limit group and
  taking the minimum finite value across all durations). NaN where a limit is not
  configured; both fields round-trip through pickle and `save_binary`/`load_binary`.
  Grids built from pandapower never populate these (getters return an empty array /
  NaN per element). Not wired into `LightSimBackend.thermal_limit_a`, which keeps its
  existing behaviour.
- [ADDED] `lightsim2grid.network.init_from_powermodels`: a native, feature-complete
  loader that builds an `LSGrid` directly from a PowerModels.jl network data
  dictionary (see https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/):
  bus, branch (line/transformer split), gen, load, shunt, **storage** and dcline
  (HVDC), including distributed-slack handling. Since PowerModels' dict is strictly
  richer than raw MATPOWER (separate load/shunt/storage tables, independent per-side
  line charging, no `tap == 0` sentinel to work around), this is now the shared engine
  other loaders build on: `init_from_matpower` converts its raw `bus`/`gen`/`branch`/
  `dcline` matrices into a PowerModels-style dict and delegates here (its own public
  behavior is unchanged; verified against its existing test suite), rather than the
  other way around.
- [ADDED] `lightsim2grid.network.init_from_pf_delta`: builds an `LSGrid` from one row
  of the PFΔ benchmark dataset (arXiv:2510.22048, a PowerModels.jl-format dict with a
  solved power-flow state attached). Accepts either a parsed row (dict) or a path to
  its `.json` file; a thin wrapper unwrapping the row's `"network"` key and calling
  `init_from_powermodels` directly, since PFΔ rows already are PowerModels dicts.
  Tested in `test_LSGrid_pf_delta.py`, including an end-to-end check that
  lightsim2grid's own AC powerflow reproduces a row's solved `vm`/`va`/`pf`/`qf`/`pt`/
  `qt` to solver tolerance, and synthetic-network checks of the `"storage"` and
  `"dcline"` conversions (PFΔ's own pglib-derived cases never contain either, but the
  schema and `init_from_powermodels` both support them).
- [ADDED] `init_from_pypowsybl(..., fuse_zero_impedance_branches=True)`: a line or
  2-winding transformer whose per-unit impedance is (near-)zero has its two terminal
  buses fused into a single electrical node instead of contributing a `1/Z` admittance
  (`Inf` for an exact zero, which broke the sparse LU factorization outright -- found on
  a real RTE grid with a genuine `r=x=0` line). Mirrors PowSyBl OpenLoadFlow's
  `lowImpedanceBranchMode`/`lowImpedanceThreshold` (new `zero_impedance_threshold_pu`
  parameter, default `1e-8` pu, matching OLF's default). A zero-impedance transformer is
  only fused if it is also at (near-)neutral tap **and** its two sides are at the same
  nominal voltage: PowSyBl's per-unit `rho` is the deviation from the transformer's own
  rated ratio (the tap-changer effect), not its absolute turns ratio, so `rho~=1` alone
  does not mean "no transformation" for a genuine step-up/down transformer. A
  zero-impedance *line* spanning two different nominal voltages raises a `RuntimeError`
  (inconsistent grid data) rather than being silently fused. Off by default (existing
  behaviour unchanged); tested in `test_bus_fusion_pypowsybl.py`.
- [FIXED] `init_from_pypowsybl` no longer crashes (`GeneratorContainer::init: ...min_q
  being above max_q...`) when a generator's reactive capability curve has a malformed
  point (`min_q > max_q` at that active power, a data-entry error) that, after the
  curve is interpolated at the generator's `target_p`, yields an inverted interval.
  OpenLoadFlow tolerates this silently; lightsim2grid now sorts the pair instead of
  hard-rejecting it (found on a real RTE grid snapshot). Tested in
  `test_gen_reactive_curve_swap.py`.
- [FIXED] `init_from_pypowsybl` no longer crashes (`KeyError: "[''] not in index"`) when
  a *connected* generator or SVC remotely regulates the voltage of a terminal that is
  itself disconnected (e.g. a de-energized voltage level). PowSyBl's bus-view id for a
  disconnected element is `''`, not NaN, and could not be resolved to any bus; the
  controller now falls back to local voltage control instead of crashing, matching
  OpenLoadFlow's behaviour (found on a real RTE grid snapshot). Tested in
  `test_remote_control_disconnected_target.py`.
- [FIXED] `TimeSeriesCPP` / `ContingencyAnalysisCPP` (and so the python `TimeSerie`,
  `ContingencyAnalysis` and `SecurityAnalysis` wrappers) could report `NaN` (AC) or
  wildly wrong, physically meaningless flows of several GW / `Inf` (DC) for a line or
  transformer imported as "half-open" (one side open, see `keep_half_open_lines` in
  `init_from_pypowsybl`). `BaseBatchSolverSynch::compute_amps_flows` /
  `compute_active_power_flows` only checked the branch's *global* status before indexing
  `_voltages` / `bus_vn_kv` with each side's bus id, never each side's own status; for a
  half-open branch the open side's bus id is the `_deactivated_bus_id` (`-1`) sentinel,
  so this read out of bounds. AC results happened to often be numerically 0 by
  incidental cancellation (the open side's Kron-reduced `yac_eff_*` coefficient is
  already 0), except for the amps of a branch whose *measured* side ("or"/side 1) was
  itself the open one, which divided by a bogus voltage base and produced `NaN`; the DC
  path has no such Kron reduction (`ydc_11`/`ydc_12` are not reduced) and mixed in an
  arbitrary, disconnected bus's angle unconditionally. Both functions now substitute an
  explicit `0` for an open side's voltage (matching
  `TwoSidesContainer_rxh_A::compute_results_tsc_rxha_no_amps`) instead of indexing with
  `-1`, and force the DC flow to `0` whenever either side is open (matching `fillBdc`'s
  "disco on one side == disco on both sides" DC convention). Tested in
  `test_line_disco_one_side.py`.
- [FIXED] `LSGrid.get_lodf()` indexed `id_me_to_dc_solver_` directly with a branch's raw
  (grid-numbering) bus id, with no status check of any kind (not even the branch's
  *global* one). For a half-open branch (`keep_half_open_lines`) the open side's bus id
  is the `_deactivated_bus_id` (`-1`) sentinel, and `id_me_to_dc_solver_` casts negative
  indices to a huge unsigned offset before indexing its backing `std::vector`, so this
  was a near-certain out-of-bounds crash, reachable from plain Python via
  `init_from_pypowsybl(keep_half_open_lines=True)` -> `dc_pf(...)` -> `get_lodf()`.
  `BaseDCAlgo::get_lodf()` had a second, independent bug in the same spot: its
  `_deactivated_bus_id` guard set the LODF column to `NaN` but never `continue`d, so the
  very next line silently overwrote that with `PTDF.col(-1)` -- again indexing with `-1`.
  Both are fixed: `LSGrid.get_lodf()` now checks each branch's own per-side DC
  connectivity (half-open or fully disconnected both count, matching `fillBdc`'s "disco
  on one side == disco on both sides") before resolving a solver bus id, propagating the
  `_deactivated_bus_id` sentinel instead of indexing with it; `BaseDCAlgo::get_lodf()`
  now gives such a branch the identity treatment (row and column all `0` except a `1` on
  the diagonal) and `continue`s, matching the physical fact that a branch already
  carrying no DC flow has zero impact anywhere, whether it is itself "outaged" or some
  other line is. Tested in `test_line_disco_one_side.py`.
- [FIXED] `ContingencyAnalysis.cpp`'s internal `check_current_violations` (feeds
  `ContingencyAnalysisCPP.compute_limit_violations` / `get_violations` / `get_violations_n`,
  used by the python `ContingencyAnalysis` / `SecurityAnalysis` wrappers) had the same
  class of bug as the `BaseBatchSolverSynch` one above: it only checked a branch's
  *global* status before resolving and indexing with each side's raw bus id (the
  `_deactivated_bus_id` check ran only after that indexing had already happened). Now
  each side's bus id is only resolved/indexed when that side's own status is connected
  (mirroring the `BaseBatchSolverSynch` fix): an open side's voltage is substituted with
  `0`, and for DC (no Kron-reduced coefficients to cancel a bogus open-end angle, unlike
  AC) both sides report `0` current whenever either side is open, matching `fillBdc`'s
  convention. Tested in `test_line_disco_one_side.py`.
- [HARDENING] `HvdcLineContainer::fillSbus` / `compute_results` and
  `LSGrid::fill_hvdc_droop_solver_data` resolve an angle-droop ("AC emulation") HVDC
  line's two converter buses using the *masked* per-side accessor (`-1` if that side is
  individually disconnected) before indexing `id_grid_to_solver`/`id_me_to_solver`, with
  no per-side status check of their own -- only relying on the invariant that
  `LSGrid::deactivate_dcline_side1`/`_side2` and
  `HvdcLineContainer::disconnect_if_not_in_main_component` always call `disable_droop`
  whenever a converter is individually opened (droop across an open converter has no
  physical meaning: the remote angle used by the linear P(theta) relationship is gone).
  This currently always holds (no other code path can set `droop_enabled_` per-element),
  so this was not a live bug, but it depended on every future half-open-creating code
  path remembering to also call `disable_droop`, with no immediate feedback if one
  didn't. All three functions now explicitly check both sides are connected before
  resolving/indexing, and raise a `RuntimeError` instead of silently indexing with `-1`
  if that invariant is ever violated. Existing HVDC test suite (`test_hvdc_*`, 35 tests)
  unaffected.
- [BREAKING] (binary) bumped `BINARY_FORMAT_VERSION` (`src/core/BinaryArchive.hpp`) from
  ``1`` to ``2`` for the new `LSGrid._init_kwargs` field below: a grid ``save_binary``'d
  with a previous lightsim2grid version can no longer be `load_binary`'d. Pickle files are
  unaffected by this specific bump (pickle already requires an exact matching
  MAJOR.medium.patch lightsim2grid version, regardless of any `StateRes` layout change).
- [ADDED] `LSGrid._init_kwargs`: a ``dict[str, str]`` property (get/set, persisted through
  `copy()`, pickle and `save_binary`/`load_binary`) meant for a grid-building function (for
  now only `init_from_pypowsybl`) to record the kwargs it was called with, so that code
  operating on the returned `LSGrid` later does not have to separately remember them. Set
  by `init_from_pypowsybl` to ``{"sort_index": str(sort_index), "buses_for_sub": str(buses_for_sub)}``.
- [ADDED] `lightsim2grid.network.LightsimResultNetwork`: casts a converged `LSGrid` built by
  `init_from_pypowsybl` back into a pypowsybl-``Network``-shaped view, so analysis code
  written against a solved ``pypowsybl.network.Network`` (``get_buses`` / ``get_lines`` /
  ``get_generators`` / ... returning a pandas DataFrame indexed by the pypowsybl element id)
  runs unmodified against a solved lightsim2grid grid. Covers buses, lines, 2-winding
  transformers, generators, loads, shunt compensators, static var compensators, batteries,
  HVDC lines and VSC/LCC converter stations. Relies on `LSGrid._init_kwargs` above (no
  need to separately pass `sort_index` back in) and raises `NotImplementedError` for a grid
  built with the legacy `buses_for_sub=True` mode. See the "Inspecting results in a
  pypowsybl-like way" section of the documentation.
- [FIXED] a build with ``__COMPILE_MARCHNATIVE=1`` (see ``benchmarks/env_compile_all.sh``)
  reliably corrupted the heap and crashed (``free(): invalid size`` / ``double free or
  corruption``) as soon as any C++ function's return value crossed from
  ``lightsim2grid_core`` into ``lightsim2grid_cpp`` (the pybind11 bindings) -- for
  example on the very first ``grid2op.make(...)`` with a pandapower-backed grid, well
  before any solver is even selected. Root cause: ``lightsim2grid_core`` and
  ``lightsim2grid_cpp`` are two separate shared libraries (a recent split), and
  ``-march=native``/``-O3`` were only ever applied (as ``PRIVATE`` compile options) to
  ``lightsim2grid_core``. ``-march=native`` changes which SIMD instruction sets Eigen
  sees enabled (``__AVX__``/``__AVX2__``/...), which changes ``EIGEN_MAX_ALIGN_BYTES``
  and thus how Eigen aligns/allocates/frees its dynamic-size matrices; since both
  libraries instantiate the same ``Eigen::Matrix<...>`` types (every pybind11 Eigen
  return-value cast does), an object allocated under one alignment assumption and freed
  under another silently computed the wrong original-allocation offset. Fixed by
  mirroring ``-march=native``/``-O3`` onto ``lightsim2grid_cpp`` in
  ``src/bindings/python/CMakeLists.txt``, the same way ``__SANITIZE``/``__DEBUG_ASSERTS``
  were already mirrored there. Added a regression test,
  ``.github/workflows/main.yml``'s ``test_march_native`` job, since release wheels do
  not use ``-march=native`` and so never exercised this.
- [FIXED] ``lightsim2grid.compilation_options.compiled_march_native`` /
  ``compiled_o3_optim`` always reported ``False``, regardless of
  ``__COMPILE_MARCHNATIVE``/``__O3_OPTIM``, because ``binding_module.cpp`` (compiled
  into ``lightsim2grid_cpp``) checks the ``__COMPILE_MARCHNATIVE``/``__O3_OPTIM``
  *preprocessor macros*, which were never defined for that target: ``__O3_OPTIM`` was
  only ``target_compile_definitions``'d as ``PRIVATE`` on ``lightsim2grid_core`` (so it
  did not propagate across the library boundary, unlike ``KLU_SOLVER_AVAILABLE`` and
  friends, which are ``PUBLIC``), and ``__COMPILE_MARCHNATIVE`` was never defined as a
  macro anywhere at all -- only used to conditionally add the ``-march=native`` compiler
  *flag* (a separate thing from the macro). This was a pre-existing, purely cosmetic gap
  (the actual ``-march=native``/``-O3`` compiler flags were, and are, applied correctly);
  it just happened to surface right after the fix above, when checking that the build
  driving a benchmark actually had ``-march=native`` active. Fixed alongside the ABI fix
  above, in the same ``src/bindings/python/CMakeLists.txt`` block.
- [ADDED] a runtime guard, on top of the CMake-level ``-march=native``/``-O3`` matching
  above, against the exact class of bug those two ``[FIXED]`` entries describe.
  ``src/core/Ls2gAbiTag.hpp`` defines an "ABI tag" (``EIGEN_MAX_ALIGN_BYTES``, the
  resolved ``EIGEN_VECTORIZE_*`` flags, and the Eigen version) computed independently in
  whichever translation unit calls it; comparing two independently-evaluated tags is the
  only way to observe this kind of drift between separately-compiled binaries (a
  ``static_assert`` cannot, since no single translation unit has visibility into another,
  separately-compiled one's flags). Two checkpoints now compare tags and refuse to
  proceed with a clear error instead of silently corrupting the heap: (1)
  ``AlgorithmRegistry::register_solver`` for third-party solver plugins loaded via
  ``load_algorithm_plugin``/``AlgorithmRegistrar``, and (2) ``lightsim2grid_cpp``'s
  module init, comparing itself against ``lightsim2grid_core``, catching this bug's own
  original failure mode directly at import time.
- [FIXED] ``KLULinearSolver`` never called SuiteSparse's ``klu_defaults()``, so its
  ``klu_common`` control struct was left all-zero (from ``common_ = klu_common();``)
  instead of the library's actual defaults. In particular ``tol=0`` (should be ``0.001``)
  disabled partial-pivoting's diagonal-preference safety, ``scale=0`` (should be ``2``)
  disabled row scaling, ``btf=0`` (should be ``TRUE``) disabled block-triangular
  preordering, and critically ``halt_if_singular=FALSE`` (should be ``TRUE``) made
  ``klu_factor``/``klu_refactor`` silently return ``KLU_OK`` even on a numerically
  singular factorization (found with ``rcond=0``). This caused ``NR_KLU`` to diverge
  (``InifiniteValue``, one Newton step after a degenerate factorization) on a real grid
  (``PtFige-20240807-2300``) where ``NR_SparseLU`` converges fine from the same seed --
  previously misdiagnosed as an inherent SparseLU-vs-KLU numerical-sensitivity artifact.
  ``klu_defaults(&common_)`` is now called in the constructor, ``reset()`` and
  ``analyze()``. Tested in ``test_KLUSolver.py``.
- [BREAKING] (cpp only) every built-in solver's ``LinearSolver`` template parameter
  (``NR_*``/``DC_*``/``FDPF_*``, see ``Solvers.hpp``) is now wrapped in
  ``LinearSolverPolicy<...>`` (``src/core/linear_solvers/LinearSolverPolicy.hpp``): a
  transparent, non-virtual pass-through that counts and times every
  analyze/factorize/refactorize/solve call. ``NRAlgo``/``BaseDCAlgo``/``BaseFDPFAlgo`` no
  longer keep their own ``timer_factor_``/``timer_refactor_``/``timer_initialize_``
  members -- ``get_timers_jacobian()`` now reads those ``TimerJac`` fields live from the
  wrapper instead, with identical externally-observable semantics (still reset every
  ``compute_pf``/``compute_pf_dc`` call). A C++ plugin subclassing ``NRAlgo``/
  ``BaseDCAlgo`` and referencing those protected members directly needs updating; a
  plugin only using the public ``LinearSolver`` API (``analyze``/``factorize``/
  ``refactorize``/``solve``/``reset``, e.g. ``examples/dist_slack_algorithm/``) is
  unaffected.
- [ADDED] ``LinearSolverStats`` (``src/core/linear_solvers/LinearSolverStats.hpp``,
  exported to python): per-call counters (``nb_analyze``/``nb_factorize``/
  ``nb_refactorize``/``nb_refactorize_failed``/``nb_fallback_factorize``/
  ``nb_fallback_factorize_failed``/``nb_solve``/``nb_reset``) and matching durations
  (``timer_initialize_``/``timer_factor_``/``timer_refactor_``/``timer_solve_``) for the
  linear solver backing a solver instance. Counters accumulate over the algorithm's whole
  lifetime (so an occasionally-firing fallback is distinguishable from a systematic one);
  the timer fields reset every call, like ``get_timers_jacobian()``. Available as
  ``solver.get_linear_solver_stats()`` on any solver (``model.get_solver()``, or the
  concrete ``NR_KLU``/``DC_KLU``/... type), and as ``get_linear_solver_stats_bp()`` /
  ``get_linear_solver_stats_bpp()`` on the two-linear-solver ``FDPF_*`` family.
- [ADDED] ``RefactorRetryLinearSolver<LinearSolver>``
  (``src/core/linear_solvers/RefactorRetryLinearSolver.hpp``, ``final``, derives from
  ``LinearSolverPolicy<LinearSolver>``): if a ``refactorize()`` call fails, falls back to
  a full ``factorize()`` (reusing the existing symbolic factorization) before reporting an
  error, tracked separately via ``LinearSolverStats.nb_fallback_factorize`` /
  ``nb_refactorize_failed``. A SuiteSparse-recommended defensive measure for KLU,
  generalized here to any solver with a real factorize/refactorize distinction. New
  built-in algorithms ``NRRefactorRetry_KLU``, ``NRRefactorRetry_CKTSO`` and
  ``NRRefactorRetry_NICSLU`` use it (``SparseLU`` is skipped: its ``factorize()`` and
  ``refactorize()`` are already the same call, so the fallback would be a no-op).
- [FIXED] ``docs/algorithm_names.rst``: the string-only ``NRRefactorRetry_*`` example used
  ``env.backend._grid.change_algorithm(...)`` without warning that this does not survive
  ``env.reset()`` (verified empirically: a reset always re-applies whatever
  ``AlgorithmType`` was last set via ``set_algo_type``), and that ``set_algo_type`` cannot
  be used instead since it requires an actual ``AlgorithmType`` value, which these
  string-only solvers do not have.
- [FIXED] ``docs/solvers.rst``: documented that ``grid.set_ac_algo_config`` /
  ``set_dc_algo_config`` customizations, just like a direct ``change_algorithm`` call,
  do not survive ``env.reset()`` either (verified empirically: ``env.backend._grid`` is a
  new object after reset, with a fresh default ``AlgoConfig``).
- [FIXED] ``docs/physical_law_checker.rst``: removed the misleading suggestion that
  lightsim2grid can only load grids compatible with its own loaders; clarified that
  ``PhysicalLawChecker`` specifically always looks for a pandapower-format ``grid.json``
  in the environment folder regardless of which formats lightsim2grid itself supports
  (see the new ``network-init-formats`` cross-reference in ``docs/network.rst``).
- [FIXED] ``docs/use_solver.rst``: documented the Jacobian row layout (P-mismatch rows
  then Q-mismatch rows, mirroring the column layout) and how to map a solver bus id back
  to the stable GridModel bus id with ``id_ac_solver_to_me``; rewrote the "constraints
  not checked" warning to match what ``compute_pf_with_input_validation`` actually
  validates now (most conditions raise a clean ``RuntimeError``/``IndexError``; only
  "every bus covered by ref/pv/pq" and ``slack_weight`` positivity/sum-to-one are still
  silently unchecked, and the CSC-format requirement is obsolete -- any scipy sparse
  format is accepted); replaced "iteratively update the jacobian matrix J" with "solve
  the linear system J.dx = mismatch" for the 8 Newton-Raphson solver descriptions;
  documented (verified empirically) that ``NRSing_*`` solvers do **not** convert extra
  ``ref`` buses to PV like the Gauss-Seidel and Fast-Decoupled solvers do -- they keep
  every bus in ``ref`` fully fixed (angle and magnitude) with no distributed-slack
  coupling; added a recommendation to remove all but one generator from the slack
  regardless of which solver is used; clarified which 6 (out of 8) solvers are marked as
  possibly-unavailable in each family.
- [FIXED] renamed the ``security_analysis`` example variable to ``contingency_analysis``
  throughout ``lightsim2grid/contingencyAnalysis.py``'s own docstrings, left over from
  before the class was renamed from ``SecurityAnalysis``.
- [FIXED] translated the last remaining French comments in ``pyproject.toml`` to English,
  for consistency with the rest of the (English-only) documentation and codebase.
- [ADDED] ``ContingencyAnalysisCPP.violation_threshold`` (mirrored by the python
  ``ContingencyAnalysis.violation_threshold``): a ``float`` in ``]0., 1.]``, default
  ``1.0``, applied to every check performed when ``compute_limit_violations`` is ``True``,
  so that near-limit situations can be reported before they actually breach 100% of a
  limit. It is the fraction of the usable range still considered acceptable, so lowering it
  makes every check stricter (more violations reported, never fewer). Each of the three
  checks owns one interval, running from a "healthy" anchor to the limit that can be
  violated, and the threshold moves that limit towards its anchor -- one linear
  interpolation, ``effective_limit = threshold * limit + (1 - threshold) * anchor``, for all
  three:

  - ``CURRENT``: anchor ``0``, limit ``limit_a`` -- violates when
    ``value >= threshold * limit_a``. A line's usable range really is ``[0, limit_a]``, so
    this reduces to plainly scaling the limit;
  - ``LOW_VOLTAGE``: anchor ``vn_kv``, limit ``vmin_kv`` -- violates when
    ``v <= threshold * vmin + (1 - threshold) * vn``;
  - ``HIGH_VOLTAGE``: anchor ``vn_kv``, limit ``vmax_kv`` -- violates when
    ``v >= threshold * vmax + (1 - threshold) * vn``.

  A voltage bound has no natural "zero" end, which is why the bus nominal voltage is used as
  its anchor (operating limits are conventionally expressed as +/- x% of it). Anchoring both
  voltage checks there -- rather than on each other -- keeps them independent: the
  ``LOW_VOLTAGE`` verdict never depends on ``vmax_kv`` and vice versa. It also means the two
  effective bounds converge towards the anchor from their own side and can never cross, so
  no bus is ever reported as both too low and too high (the sole exception being a band
  fully collapsed onto nominal, ``vmin == vn == vmax``, where a bus sitting exactly at
  nominal satisfies both and is reported once, as ``LOW_VOLTAGE``). The anchor is the
  nominal voltage *clamped into* ``[vmin_kv, vmax_kv]``: a band is not required to bracket
  it -- a 380 kV level declared with an operating range of ``[390, 450]`` kV is ordinary
  practice on the European 400 kV network -- and clamping keeps each bound moving inwards
  only, without rejecting such a grid. Where the band does bracket the nominal voltage (the
  overwhelmingly common case) the clamp does nothing.

  The reported ``LimitViolation.value`` / ``.limit`` are never rescaled (still the value
  reached and the limit as configured); only the test deciding whether to report is shifted.
  The default ``1.0`` reproduces the previous, threshold-less behaviour (modulo the strict
  ``>`` / ``<`` comparisons becoming ``>=`` / ``<=``, which only differ if a value lands
  exactly on its limit). A value outside ``]0., 1.]`` (or ``NaN``) is rejected. Like
  ``nb_thread`` / ``handle_disconnected_grid`` this is a plain runtime knob affecting only
  the next ``compute()``, with one asymmetry: *lowering* it discards any already-computed
  results (they would under-report), while *raising* it keeps them (a stricter result is a
  superset of a looser one). Unlike ``compute_limit_violations``, neither case clears the
  registered contingencies.
- [ADDED] a bus configured with a minimum voltage above its maximum
  (``vmin_kv > vmax_kv``, see ``LSGrid.set_bus_voltage_limits``) now raises a
  ``RuntimeError`` during a ``ContingencyAnalysis`` limit-violation check, instead of being
  silently reported as an arbitrary one of ``LOW_VOLTAGE`` / ``HIGH_VOLTAGE``. That is a
  genuine input error, and the only way the two effective voltage bounds can end up crossed
  (see ``violation_threshold`` above).
- [ADDED] ``src/tests/test_powerflow_algorithm.cpp``: a C++ (Catch2) test suite for the
  Newton-Raphson system itself (``src/core/powerflow_algorithm/NRSystem*``) and its
  ``MultiSlack`` / ``VoltageControl`` / ``Hvdc`` extensions. It drives an ``NRSystem``
  phase by phase, the way ``NRAlgo::compute_pf`` does, and checks the internal
  consistency the grid-level tests cannot see: the augmented Jacobian is square and
  matches ``total_state_variables()``, every bus -> row / column map stays in range and
  agrees with the ledger's registration pair lists, no Jacobian row or column is
  structurally empty (which is how an unresolved feature entry shows up), a
  ``status_droop`` flip changes values but not the sparsity pattern (so a linear solver
  may reuse its symbolic factorization), bus masking is a pure value-level edit, and
  ``clear_jacobian`` leaves a system that rebuilds bit-identically.
- [ADDED] two jobs in ``.github/workflows/sanitizers.yml`` running the C++ (Catch2) test
  suite under ASan + UBSan and under re-enabled Eigen / libstdc++ assertions
  (``__SANITIZE=1`` / ``__DEBUG_ASSERTS=1``). The suite already ran under valgrind
  (``cpp_unit_tests.yml``), but an index that is wrong while still inside its allocation
  is invisible to both valgrind and ASan -- only the assertion build catches it. The test
  binary itself (not just ``lightsim2grid_core``) is now compiled and linked with the
  sanitizer flags, see ``src/tests/CMakeLists.txt``.
- [ADDED] ``modify_gen_v`` on ``TimeSeriesCPP`` / ``InjectionSweepCPP`` /
  ``ScenarioSweepCPP`` (python: ``TimeSerie.modify_gen_v`` /
  ``ScenarioSweep.modify_gen_v``, the latter inherited "for free" by
  ``InjectionSweep``): a new per-step axis, the target generator voltage magnitude
  (``vm_pu``), shape ``(n_simul, n_gen)``. Unlike ``modify_gen_p`` / ``modify_sgen_p`` /
  ``modify_load_p`` / ``modify_load_q``, this does **not** feed the injection (``Sbus``)
  at all: a PV bus's magnitude is not part of Newton-Raphson's unknown vector, so it
  never moves during a solve once seeded -- ``modify_gen_v`` only re-seeds ``|V|`` at
  each voltage-regulating generator's regulated bus immediately before that row's solve
  (``GeneratorContainer::set_vm``, a new overload taking an explicit per-generator
  vm vector, alongside the existing single-arg overload reading the grid's own
  ``target_vm_pu``). Left unset (the default), every row keeps using the grid's own
  ``target_vm_pu``, exactly as before this setter existed.
- [IMPROVED] speed (DC batch algorithms -- ``TimeSeriesCPP`` / ``InjectionSweepCPP`` /
  ``ContingencyAnalysisCPP`` / ``ScenarioSweepCPP``): a DC solve no longer builds the
  complex voltage every row. ``BaseDCAlgo::compute_pf_dc`` only ever *outputs* the bus
  angles (``theta``, the linear solve's actual result) -- the magnitude ``Vm_`` is a
  pure, unmodified echo of the caller's input, and the complex ``V_`` it used to build
  from both (a ``std::polar`` / hardware ``sincos`` call per bus) existed only so
  downstream code could read it back. That downstream code (``compute_amps_flows`` /
  ``compute_active_power_flows``) was then immediately un-building it again with
  ``.arg()`` (and, for the current-limit checks, ``check_current_violations`` did the
  same) -- a full sin/cos-then-atan2 round trip, per bus, per row, for values that were
  never actually needed in complex form. A new ``BaseAlgo::set_lazy_v`` toggle (DC-only;
  a no-op default for AC / plugin solvers) lets ``BaseBatchSweep`` skip that entirely:
  the per-row sweep now accumulates only ``theta`` (a plain real matrix), and the flow
  computations read it directly, with no ``.arg()``. The (small, per-generator) voltage
  magnitude is reconstructed only if/when something actually asks for it --
  ``get_voltages()`` (lazily, cached on first call) or the amps flows (which do need
  ``|V|`` for the per-unit-to-amps conversion) -- from the grid's own target voltages
  plus ``modify_gen_v``, if set, the exact same (tiny) input the row loop itself already
  used, never from a stored per-row matrix. Bit-for-bit identical results (covered by
  the existing DC / ``modify_gen_v`` / current-limit-violation cases across
  ``src/tests/test_injection_sweep.cpp``, ``test_timeseries_sbus.cpp``,
  ``test_batch_voltage_control.cpp`` and ``test_scenario_sweep_violations.cpp``); the
  "handle disconnected grid" mode (``ContingencyAnalysisCPP`` / ``ScenarioSweepCPP``) is
  unaffected and keeps building the complex voltage eagerly, as it already needs
  per-row current-limit checks whenever that mode is combined with
  ``compute_limit_violations``.
- [FIXED] ``BaseBatchSweep``'s ``_run_range`` aborted the whole remaining range on the first
  row that failed to solve (was non-invertible, or diverged) for **every** instantiation,
  not just ``TimeSeriesCPP``. The abort condition was keyed on
  ``SbusPolicy::supports_vary && !YbusPolicy::supports_contingency``, which happens to be
  true for ``TimeSeriesCPP`` (correctly abort: each row is warm-started from the previous
  one, so nothing past a failure is meaningful) but was *also* true for
  ``InjectionSweepCPP`` (each row is independent of the others and of their order, per its
  own docstring -- a later row must not be abandoned just because an earlier, unrelated one
  failed). ``ContingencyAnalysisCPP`` / ``ScenarioSweepCPP`` were unaffected in practice
  (they go through the masked ``_run_range_masked`` path instead), but the same policy
  combination would have hit them too under ``handle_disconnected_grid``. The abort is now
  keyed directly on the real criterion, ``BatchInitKind::FromPreviousStep`` (chained rows),
  which is true only for ``TimeSeriesCPP``. Every other row still gets its own zero-row /
  ``NOT_SIMULATED``-or-``DIVERGENCE`` sentinel via the existing ``_store_row_status``, same as
  a genuinely diverging row already got; only rows *after* it, which used to be silently
  skipped, are now actually attempted.
- [ADDED] ``nb_converged()`` on every batch algorithm (``TimeSeriesCPP`` /
  ``InjectionSweepCPP`` / ``ContingencyAnalysisCPP`` -- and so ``SecurityAnalysis``, which
  wraps it -- / ``ScenarioSweepCPP``), next to the existing ``nb_solved()``: the count of
  powerflows, among those ``nb_solved()`` attempted, that actually converged. Always
  ``<= nb_solved()`` -- a row skipped outright (a non-invertible / islanding admittance
  matrix, on ``ContingencyAnalysisCPP`` / ``ScenarioSweepCPP``) never reaches the solver at
  all, so it counts towards neither. Implemented once, on the shared
  ``BaseBatchSolverSynch`` base, by threading a new ``int & nb_converged`` out-parameter
  alongside the existing ``nb_solved`` one through ``compute_one_powerflow`` and every
  caller up to ``_compute_threaded`` (including its per-thread accumulate-then-merge
  logic), so all four classes get it identically and cannot drift apart.

[0.13.1]  2026-04-21
--------------------
- [BREAKING] when loading a powergrid from pypowsybl with "use_buses_for_sub" tagged
  and disconnected element on the grid will now raise a RuntimeError. Before there were
  some "automatic" bahaviour to try to find a possible bus which could lead to 
  error afterwards.
- [BREAKING] adding a more precise information about linear solvers. The "refactor" timings
  are now also available in solver.get_timers_jacobian() which now returns a tuple of size 10
- [FIXED] an issue where disconnected powerlines could be tagged as "fakely connected"
- [FIXED] some issues when loading a grid from pypowsybl in case of disconnected elements.
- [FIXED] remove the undefined behaviour while maintaining compile time check to prevent 
  wrong conversion from different bus labelling.
- [IMPROVED] the CI to allow automatic push on pypi on new version tag (introduced in version 0.12.0)
- [IMPROVED] reduce code duplication between ContingencyAnalysis and TimeSerie (cpp side)
- [IMPROVED] handling of branches disconnected at only one side: less code duplication and
  it should be working with TimeSeries and ContingencyAnalysis
- [IMPROVED] speed (DC mode): avoid the systematic call to "refactor" when Ybus is not changed
  when using DC approximation.
- [IMPROVED] simplify the future integration of other linear solvers and the logic when linear_solvers
  are called by decoupling "refactor" steps from "solve" steps (they used to be all under the same
  "solve" method).
- [IMPROVED] added "final" and "override" key-words on some methods in src/element_container
  
[0.13.0] 2026-04-15
--------------------
- [PENDING DEPRECATION] the cpp module (lightsim2grid_cpp) will not be usable directly anymore.
  This means that calls like "from lightsim2grid_cpp import XXX" will not work. To replace them 
  you  need to perform "from lightsim2grid.lightsim2grid_cpp import XXX"
- [FIXED] some compilation issues on some systems (*eg* windows when using c++23 standard)
- [FIXED] some issues with "copy on write" and pandas 3 when init from pandapower grid.
- [IMPROVED] cleaner `cktso_lib` (`from lightsim2grid.compilation_options import cktso_lib`) : the file name and extension are omitted
- [IMPROVED] easier build by relying on cmake and scikit_build_core to build the cpp part
- [IMPROVED] SuiteSparse to version 7.12.2 (2026-02-05)

[0.12.2] 2026-02-05
----------------------
- [FIXED] an issue with shunt buses (was set to 1 even if they were disconnected)
- [FIXED] a warning when applying actions on generator voltage setpoints (due to NaN)
- [FIXED] pandapower grid could be modified when importing a grid from pandapower.
- [IMPROVED] add a test to make sure generator types are available if using
  `dist_slack_non_renew` information.
- [IMPROVED] test coverage on shunts (a test needed to be skipped due to float comparison in grid2op)

[0.12.1]  2026-01-09
---------------------
- [FIXED] phase shift transformers are now properly modeled
  for both pandapower (new in this version) and pypowsybl (already
  the case in previous version)
- [FIXED] a performance issue for all "XXXSingleSlack" (*eg* KLUSingleSlack) algorithm (filling of the initial
  Jacobian matrix was extremly slow due to the massive 'insert' of data in the eigen sparse matrix instead 
  of relying on the "setFromTriplets" method)
- [FIXED] an normal "exception" was not catched in the `close()` method of LightSimBackend in case the
  backend was closed before any grid was loaded.
- [ADDED] possibility to pickle independantly all part of the grid (*eg* gridmodel.get_lines()
  can be pickled independantly from anything else) **NB** pickling and un-pickling 
  lightsim2grid objects can only be used for the same lightsim2grid version.
- [ADDED] the `init_from_n_powerflow` property for ContingencyAnalysis (and 
  ContingencyAnalysisCPP). It allows to chose if the computation of the contingencies
  are initialized with the complex voltages resulting of the powerflow in N 
  (`init_from_n_powerflow=True`) or if they are initialized from the 
  given input vector (`init_from_n_powerflow=False`). Defaults to `False`

[0.12.0] 2026-01-06
--------------------
- [BREAKING] for better consistency, and following pypowsybl convention, trafo and lines "side"
  are now called "1" and "2" instead of "hv" / "lv" (for trafo) or "or" / "ex" for powerlines.
  For example, what used to be accessible with `gridmodel.change_bus_powerline_or(...)` is now called
  `gridmodel.change_bus1_powerline()`. This affects powerlines, transformers and dc powerlines but also
  "LineInfo", "TrafoInfo" and "DCPowerlineInfo". See below for a (should-be exhaustive) list of changes.
- [BREAKING] the `init_pp_backend` public attribute is now private (called now `_init_pp_backend`) and 
  optional, meaning it's `None` when the grid is initialized from pypowsybl (for example).
- [FIXED] some issues with the "load_grid_from_pypowsybl" function (and making sure the graph of the structure
  of the lightsim2grid gridmodel matches the one of the pypowsybl grid).
- [FIXED] an issue with the handling of the slack due to a not correct implementation
  of `update_slack_weights_by_id` cpp side (previous slacks were not removed, slacks get_slack_weights
  were prop to target_p which caused issues when all slacks had targetp==0.)
- [FIXED] an issue with serialization / de serialization caused by an error in serializing the solver types.
- [ADDED] in all "xxxInfo"  (*eg* "LoadInfo") information about subtation and position in the topology
  vector, with the `sub_id` / `pos_topo_vect` (for `LoadInfo`, `SGenInfo`, `GenInfo`, `StorageInfo`, `ShuntInfo`)
  and `sub1_id` / `sub2_id` / `pos1_topo_vect` / `pos2_topo_vect` (for `LineInfo`, `TrafoInfo` and `DCLineInfo`)
- [ADDED] possibility to load the pypowsybl grid with extra key-words arguments (by using `pypowsybl_load_kwargs` in the 
  `loader_kwargs` of LightSimBackend)
- [ADDED] possibility to initialize LightSimBackend with an already loaded grid (by using the `grid` key of the
  `loader_kwargs` of LightSimBackend when loading it with pypowsybl)
- [ADDED] possibility to change the ratio (`rho`) of transformers (`gridmodel.change_ratio_trafo(trafo_id, new_rho)`)
- [ADDED] possibility to change the phase shift (`alpha`) of transformers (`gridmodel.change_shift_trafo(trafo_id, new_alpha)`)
- [ADDED] more consistency checkings to avoid "negative buses" cpp side.
- [ADDED] information about the coefficients assigned on the Ybus matrix for `LineInfo` and `TrafoInfo`: `yac_11`, `yac_12`, 
  `yac_21`, `yac_22`, `ydc_11`, `ydc_12`, `ydc_21`, and `ydc_22`
- [ADDED] possibility to have a powerline / transformer connected on only one side (for DC and Newton-Raphson algorithm, not implemented
  for fast-decoupled yet). This means that powerlines / transformers have 3 statuses: one for each side and one "global".
- [ADDED] possibility to choose the way lightsim2grid will internally treat the powerlines status:

  - `ignore_status_global` (`gridmodel.set_ignore_status_global(True)` or `gridmodel.set_ignore_status_global(False)`). If set
     to `False` (default) the the "global" status is synch with all the others. For example you can deactivate both sides of a line by 
     deactivating `status_global` and conversely if you deactivate both side of a given line, then its "status_global" is set 
     to "activated". If `ignore_status_global` is set to `True` then `global_status` is not updated at all and ignored (*NB* in 
     this case, calling `gridmodel.deactivate_line(...)` will deactivate both sides but not the "global status")
  - `synch_status_both_side` (`self.model.set_synch_status_both_side(XXX)`). If set to `True` (default) then the status of each 
    side of any given line will be synched. Meaning that if an action disconnects one side, it will also disconnect the "status_global"
    and the other side (*NB* in this case if the same actions both connects one side and disconnect another, then the outcome is "undefined").
    If `synch_status_both_side` is `False` then each side of the powerline is independant from the other (which can lead to powerline / 
    transformer being connected at only one side).
  - The complete bahviour is tested in `tests/test_line_disco_one_side.py`. Feel free to have a look if you need more information.

- [IMPROVED] Eigen to version 5.0.1 (2025/11/11)
- [IMPROVED] rename all ".h" file to ".hpp" for cpp headers (cpp side). 
- [IMPROVED] consistency of "bus labelling" cpp side (implement compile time check to prevent accidental conversion from 
  `LocalBusId`, `GlobalBusId` / `GridModelBusId` and / or `SolverBusId`)

Table to upgrade the names:

=============================    =======================  =======================
Class Name                       Old Attribute Name       New Attribute Name
=============================    =======================  =======================
TrafoContainer                   get_bus_from             get_bus_id_side_1
TrafoContainer                   get_bus_to               get_bus_id_side_2
TrafoInfo                        bus_hv_id                bus1_id
TrafoInfo                        bus_lv_id                bus2_id
TrafoInfo                        connected                connected_global
TrafoInfo                        h_pu                     h1_pu or h2_pu
TrafoInfo                        is_tap_hv_side           is_tap_side_1
TrafoInfo                        res_p_hv_mw              res_p1_mw
TrafoInfo                        res_q_hv_mvar            res_q1_mvar      
TrafoInfo                        res_v_hv_kv              res_v1_kv    
TrafoInfo                        res_a_hv_ka              res_a1_ka         
TrafoInfo                        res_p_lv_mw              res_p2_mw      
TrafoInfo                        res_q_lv_mvar            res_q2_mvar    
TrafoInfo                        res_v_lv_kv              res_v2_kv     
TrafoInfo                        res_a_lv_ka              res_a2_ka     
TrafoInfo                        res_theta_hv_deg         res_theta1_deg    
TrafoInfo                        res_theta_lv_deg         res_theta2_deg
LineContainer                    get_bus_from             get_bus_id_side_1
LineContainer                    get_bus_to               get_bus_id_side_2
LineInfo                         connected                connected_global
LineInfo                         bus_or_id                bus1_id
LineInfo                         bus_ex_id                bus2_id
LineInfo                         h_pu                     removed
LineInfo                         h_or_pu                  h1_pu
LineInfo                         h_ex_pu                  h2_pu
LineInfo                         res_p_or_mw              res_p1_mw
LineInfo                         res_q_or_mvar            res_q1_mvar      
LineInfo                         res_v_or_kv              res_v1_kv    
LineInfo                         res_a_or_ka              res_a1_ka         
LineInfo                         res_p_ex_mw              res_p2_mw      
LineInfo                         res_q_ex_mvar            res_q2_mvar    
LineInfo                         res_v_ex_kv              res_v2_kv     
LineInfo                         res_a_ex_ka              res_a2_ka     
LineInfo                         res_theta_or_deg         res_theta1_deg    
LineInfo                         res_theta_ex_deg         res_theta2_deg
DCLineContainer                  get_bus_from             get_bus_id_side_1
DCLineContainer                  get_bus_to               get_bus_id_side_2
DCLineInfo                       connected                connected_global
DCLineInfo                       bus_or_id                bus1_id
DCLineInfo                       bus_ex_id                bus2_id
DCLineInfo                       target_vm_or_pu          target_vm1_pu         
DCLineInfo                       target_vm_ex_pu          target_vm2_pu 
DCLineInfo                       gen_or                   gen1         
DCLineInfo                       gen_ex                   gen2
DCLineInfo                       res_p_or_mw              res_p1_mw
DCLineInfo                       res_q_or_mvar            res_q1_mvar      
DCLineInfo                       res_v_or_kv              res_v1_kv          
DCLineInfo                       res_p_ex_mw              res_p2_mw      
DCLineInfo                       res_q_ex_mvar            res_q2_mvar    
DCLineInfo                       res_v_ex_kv              res_v2_kv      
DCLineInfo                       res_theta_or_deg         res_theta1_deg      
DCLineInfo                       res_theta_ex_deg         res_theta2_deg      
=============================    =======================  =======================


[0.11.0] 2025-12-09
----------------------
- [DEPRECATED] python 3.7 builds will no longer be available
- [FIXED] a bug in the import of the grid from pypowsybl
- [FIXED] a bug with phase shifters in case the tap was on the low voltage side
- [FIXED] a bug with active shunt values (wrong sign in the cpp part) and wrong Ybus diagonal coeff
- [FIXED] a bug in DC computation with shunt active values (wrong sign)
- [FIXED] a bug in DC computation with some phase shifters (when tap was not tagged on correct side)
- [FIXED] a bug in FDPF when phase tap changer was not on high voltage side
- [FIXED] a lots of bug in the conversion of pypowsybl grid when using 
  "old" pypowysbl versions.
- [ADDED] compatibility with python 3.14 and python 3.14 build
- [ADDED] compatibility with pandapower >= 3 version when loading a grid 
  (pandapower changed the way it initilizes the transformers model parameters)
- [ADDED] more kwargs arguments are possible in the LightSimBackend `loader_kwargs`
- [ADDED] name of the substations are now read from the grid when initializing from 
  pypowsybl.
- [ADDED] support for multiple slack when reading a grid from pypowsybl.
- [IMPROVED] the way to initialize the transformers from pypowsybl
- [IMPROVED] possibility to load grid with phase shifters from pypowsybl
- [IMPROVED] function to initialize the grid from pypowsybl has now a 
  basic documentation.

[0.10.3] 2025-04-28
----------------------
- [FIXED] remove deprecated use of numpy<2 function in LightSimBackend
- [FIXED] the "copy.deepcopy" of a lightsim2grid backend does not crash anymore 
  (see issues #36 and #97)
- [FIXED] an issue that could lead to a "segfault" if the "runpf" method of
  `LightSimBackend` was called first with is_dc=True and then with is_dc=False
- [IMPROVED] compat with grid2op 1.11.0
- [IMPROVED] now test proper compilation on clang 20 (was limited to clang 18 before)

[0.10.2] 2025-03-07
----------------------
- [FIXED] an error when changing of bus one of the slack (did not trigger the 
  recompute of pv bus ids)
- [FIXED] an issue when turning off a generator: it was still declared as "slack"
  if it was one.
- [FIXED] could not disconnect a generator when it was a slack bus
- [FIXED] voltage was -1 instead of 0 for disconnected elements (load, generator, storage units etc.)
- [ADDED] an option in `LightSimBackend` automatically disconnect load and generators
  if they are not in the main connected component.
- [IMPROVED] refactoring of the c++ side container element to reduce
  code (for "one end" elements such as loads, generators, static generators and shunts)

[0.10.1] 2025-01-04
----------------------------
- [FIXED] some timings on the benchmarks were not measured at the right time
- [ADDED] more benchmarks especially for DC powerflow
- [ADDED] a `dcpf` function that can replace the pandapower `dcpf` interal function
- [IMPROVED] benchmark on the documentation
  (clarity of what is done)
- [IMPROVED] consistency of the names and measured times accross the different benchmarks

[0.10.0] 2024-12-17
-------------------
- [BREAKING] disconnected storage now raises errors if some power is produced / absorbed, when using legacy grid2op version,
  you can retrieve the previous behaviour by initializing the `LightSimBackend` with
  `backend = LightSimBackend(..., stop_if_storage_disco=False)`
- [BREAKING] with the new `detachment_is_allowed` feature in grid2op, the kwargs `stop_if_load_disco`,
  `stop_if_gen_disco` (and `stop_if_storage_disco`) are now optional. They are set up from the 
  call to `grid2op.make(...)` and are erased by the `allow_detachment` kwargs. In other words,
  you don't need to set `stop_if_load_disco`, `stop_if_gen_disco` or `stop_if_storage_disco`. It is 
  automatically set by `grid2op.make(..., allow_detachment=XXX)` to have the correct bahaviour.
- [FIXED] an issue with the storage units (when asking it to produce / consume 
  but deactivating them with the same action the grid did not diverge)
- [IMPROVED] add the grid2op "detachement" support (loads and generators are allowed
  to be disconnected from the grid)
- [ADDED] a kwargs `stop_if_storage_disco` to control (in legacy grid2op version) the behaviour 
  of the backend when a storage unit is disconnected.

[0.9.2.post2] 2024-11-29
--------------------------
- [FIXED] The attribute `can_output_theta` (of base `Backend` class)
  was not set to `True` if using the pypowsybl loader.
- [FIXED] the github CI to work properly on many linux buit image
- [IMPROVED] build on python 3.13

[0.9.2.post1] 2024-11-28
--------------------------
- [FIXED] There is still a bug with the pypowsybl 1.8.1 version with the 
  tap changer ratio (unconsistency between what needs to be done and the 
  actual documentation). The fix is to set the const variable `PP_BUG_RATIO_TAP_CHANGER`
  to be at least 1.9.0, otherwise results are wrong.
  
[0.9.2] 2024-10-18
--------------------------
- [ADDED] support loading a grid when everything is NOT on the same bus
  (`topo_vect` used to be wrong in this case). This is especially usefull
  for grid loaded with `pypowsybl`
- [ADDED] a file benchmarking the timings for running powerflow on different
  grid sizes.
- [UPDATED] urls to match the new repo location
- [UPDATED] urls to match new grid2op location

[0.9.1] 2024-09-30
--------------------------
- [FIXED] a bug due to wrong type (in a numpy array) for the element name which lead in turn 
  to a fail assertion (equality between two numpy arrays returning a bool and not an array)
- [FIXED] a bug when init a grid from pypowsybl: the wrong value was used for trafos `h` (double)
- [FIXED] a bug when init a grid from pypowsybl: wrong values for `_ls_to_orig` and `_orig_to_ls`
  was set (and later used)
- [FIXED] yet another bug when init a grid from pypowsybl: the voltage in kV (not in pu)
  could be set due to "wrong" labelling of the bus ids
- [FIXED] yet another bug when init a grid from pypowsybl: the ratio of the transformers 
  sent in lightsim2grid did not take into account the "`rated_u1` `rated_u2`" on both side 
  (only used on one side)
- [FIXED] yet another bug when init a grid from pypowsybl: the ratio of the transformers 
  sent in lightsim2grid did not take into account the ratio in the  `pypow_net.get_ratio_tap_changers()`
- [ADDED] a method for the `ContingencyAnalysisCPP` class that returns, for all contingencies
  in the contingency list, which will be simulated and which causes the grid to be disconnected.
- [ADDED] it is now possible to use "one substation" (voltage level) pypowsybl side is
  one substation in lightsim2grid.
- [IMPROVED] removing a weird `1j * h_` when initializing powerlines and transformers. This was 
  part of a pandapower "hack" which is not present anymore (see 
  https://github.com/Grid2Op/lightsim2grid/issues/88#issue-2443299039)

[0.9.0] 2024-07-29
--------------------------
- [BREAKING] installing pandapower lightsim2grid does not require anymore to install
  pandapower (you can initialize `GridModel` with pypowsybl or pandapower if you want). To make it both
  cleaner and clearer the function `lightsim2grid.gridmodel.init` has been removed.
  Please use `lightsim2grid.gridmodel.init_from_pandapower` or 
  `lightsim2grid.gridmodel.init_from_pypowsybl` from now on
- [BREAKING] the previous `gridmodel.get_ptdf()` function was wrongly labelled with the
  "solver" bus id and not the `gridmodel` bus id which could cause issue when it was computed 
  on some grid configuration. It has now been fixed (so the `gridmodel.get_ptdf` returns the
  proper things). If you want the previous behaviour, you need to use `gridmodel.get_ptdf_solver()`
- [BREAKING] similarly, `gridmodel.get_Ybus()`, `gridmodel.get_dcYbus()`, `gridmodel.get_Sbus()`
  and `gridmodel.get_dcSbus()` now return things in the `gridmodel` bus ordering. For the previous
  behaviour you can use `gridmodel.get_Ybus_solver()`, `gridmodel.get_dcYbus_solver()`,
  `gridmodel.get_Sbus_solver()` and `gridmodel.get_dcSbus_solver()`
- [BREAKING] the more rational logic above also extends to all the functions listed in the 
  table below:

===============================    ===================================================
Function with behaviour change      Name of the new function having the same behaviour
===============================    ===================================================
gridmodel.get_ptdf()                gridmodel.get_ptdf_solver()
gridmodel.get_Ybus()                gridmodel.get_Ybus_solver()
gridmodel.get_dcYbus()              gridmodel.get_dcYbus_solver()
gridmodel.get_Sbus()                gridmodel.get_Sbus_solver()
gridmodel.get_dcSbus()              gridmodel.get_dcSbus_solver()
gridmodel.get_pv()                  gridmodel.get_pv_solver()
gridmodel.get_pq()                  gridmodel.get_pq_solver()
gridmodel.get_slack_ids()           gridmodel.get_slack_ids_solver()
gridmodel.get_slack_ids_dc()        gridmodel.get_slack_ids_dc_solver()
gridmodel.get_slack_weights()       gridmodel.get_slack_weights_solver()
gridmodel.get_V()                   gridmodel.get_V_solver()
gridmodel.get_Va()                  gridmodel.get_Va_solver()
gridmodel.get_Vm()                  gridmodel.get_Vm_solver()
gridmodel.get_J()                   gridmodel.get_J_solver()
gridmodel.get_Bf()                  gridmodel.get_Bf_solver()
===============================    ==================================================

- [FIXED] the `change_solver` in the `ContingencyAnalysis` did not work correctly.
  More specifically the solver type used might not be correct if changed which could 
  lead to wrong Ybus being passed to the solver.
- [FIXED] some compatibility mode with python `3.7`
- [FIXED] a bug when "turned off" generator were not PV (slack was 
  "turned off" when its target P was 0. But still the slack so ends up producing something...)
- [FIXED] (consistency with pandapower) when an intial powerflow is run
  to initialized an AC powerflow, the initial voltages are 1 pu (and 
  not `gridmodel.get_init_vm_pu()` as previously).
- [FIXED] `gridmodel.get_ptdf()` now have the 
  normal "gridmodel" bus id representation and not the "solver" bus ordering.
- [FIXED] `gridmodel.get_lodf()` issue wrong results in case of some
  topological modification
- [FIXED] calls to methods such as `gridmodel.get_pv` or `gridmodel.get_V` 
  or `gridmodel.get_Ybus` could lead to severe crashes (segmentation fault)
  on some (rare) cases. Now an exceptions should be thrown in these cases.
- [FIXED] basic backward compatibility is ensured and tested for legacy grid2op >= 0.9.1.post1
  Not all features are tested and only 1.x versions are tested 
  (ie 1.1 or 1.2 but not 1.2.1, 1.2.2, 1.2.3 etc.) and only for python 3.11
- [FIXED] a bug when using `LightSimBackend` with some old (but not too old) grid2op
  versions.
- [FIXED] various compatibility bugs when using old grid2op versions.
- [ADDED] it is now possible to deactivate the support for shunts by 
  subclassing the LightSimBackend class and setting the `shunts_data_available`
  to `False`
- [IMPROVED] in the `ContingencyAnalysis` class, the underlying cpp model will now
  perform an initial powerflow.
- [IMPROVED] distributed wheels are now compiled (whenever possible) with numpy 2. 
  This makes them compatible with both numpy 1.x.y and numpy 2.z.t versions.
- [IMPROVED] tests are now performed when lightsim2grid is compiled with 
  the latest clang (18) and gcc (14) versions on the CI using python 3.11

[0.8.2.post1] 2024-04-xx
--------------------------
- [FIXED] a "forward compatibility" issue with grid2op 1.10.2
  (due to wrong usage of some internal classes when loading a pandapower grid)

[0.8.2] 2024-04-22
--------------------
- [FIXED] CI was broken after migration to artifact v4, set it back to v3 
  (and make the names of the folder clearer)
- [FIXED] CI when using latest pandapower version (2.14) which broke some previous tests
- [ADDED] the computation of the LODF (line outage distribution factor) in 
  lightsim2grid
- [ADDED] some convenience functions to retrieve in a vectorized way the 
  buses to which each elements of a given container is connected 
  (*eg* `gridmodel.get_lines().get_bus_from()`)
- [ADDED] more binaries (windows `arm64` and macos `arm64`)
- [IMPROVED] remove some compilation warnings for clang
- [IMPROVED] possibility to specify generator used as slack by its name when initializing
  from `pypowsybl`.
- [IMPROVED] removing some warnings when grid2op is not installed
  (it should not raise any warning as lightsim2grid does not require grid2op)

[0.8.1] 2024-03-26
--------------------
- [FIXED] a bug with shunts when `nb_busbar_per_sub` >= 2
- [FIXED] some bugs preventing backward compatibility
- [FIXED] an issue in the computation of gen_q when intialized with pypowsybl
  (some overflow cpp side leading to infinite number in gen_q)
- [FIXED] a bug in the "containers" cpp side (wrong bus was assigned)
  when elements was disconnected, which lead to wrong computations for 
  time series or contingency analysis.
- [FIXED] another bug in ContingencyAnalysis (cpp side) leading to wrong computation
  when a powerline was disconnected
- [FIXED] some broken imports when grid2op was not installed
- [FIXED] missing "typing_extension" as required when installation
- [ADDED] some information of compilation directly in the cpp module
- [ADDED] some information of compilation available in the python `compilation_options`
  module python side
- [ADDED] some convenient methods for `ContingencyAnalysis` python side (most 
  notably the possibility to initialize it from a `LightSimBackend` and to
  change the topology of the grid)
- [ADDED] a "reward" module in lightsim2grid with custom reward
  based on lightsim2grid.
- [ADDED] a class `N1ContingencyReward` that can leverage lightsim2grid to 
  assess the number of safe / unsafe N-1.
- [IMPROVED] time measurments in python and c++
- [IMPROVED] now test lightsim2grid with oldest grid2op version
- [IMPROVED] speed, by accelerating the reading back of the data (now read only once and then
  pointers are re used)
- [IMPROVED] c++ side avoid allocating memory (which allow to gain speed python side too)
- [IMPROVED] type hinting in `LightSimBackend` for all 'public' methods (most 
  notably the one used by grid2op)
- [IMPROVED] now the benchmarks are more verbose (detailing some compilation options)

[0.8.0] 2024-03-18
--------------------
- [BREAKING] now able to retrieve `dcSbus` with a dedicated method (and not with the old `get_Sbus`).
  If you previously used `gridmodel.get_Sbus()` to retrieve the Sbus used for DC powerflow, please use
  `gridmodel.get_dcSbus()` instead.
- [DEPRECATED] in the cpp class: the old `SecurityAnalysisCPP` has been renamed `ContingencyAnalysisCPP`
  (you should not import it, but it you do you can `from lightsim2grid.securityAnalysis import ContingencyAnalysisCPP` now)
- [DEPRECATED] in the cpp class: the old `Computers` has been renamed `TimeSerieCPP`
  (you should not import it, but it you do you can `from lightsim2grid.time_serie import TimeSerieCPP` now)
- [FIXED] now voltage is properly set to 0. when shunts are disconnected
- [FIXED] now voltage is properly set to 0. when storage units are disconnected
- [FIXED] a bug where non connected grid were not spotted in DC
- [FIXED] a bug when trying to set the slack for a non existing genererator
- [FIXED] a bug in init from pypowsybl when some object were disconnected. It raises
  an error (because they are not connected to a bus): now this function properly handles
  these cases.
- [FIXED] a bug leading to not propagate correctly the "compute_results" flag when the 
  environment was copied (for example)
- [FIXED] a bug where copying a lightsim2grid `GridModel` did not fully copy it
- [FIXED] a bug in the "topo_vect" comprehension cpp side (sometimes some buses 
  might not be activated / deactivated correctly)
- [FIXED] a bug when reading a grid initialize from pypowsybl (trafo names where put in place 
  of shunt names)
- [FIXED] read the docs was broken
- [FIXED] a bug when reading a grid from pandapower for multiple slacks when slack 
  are given by the "ext_grid" information.
- [FIXED] a bug in "gridmodel.assign_slack_to_most_connected()" that could throw an error if a 
  generator with "target_p" == 0. was connected to the most connected bus on the grid
- [FIXED] backward compat with "future" grid2op version with a 
  better way to copy `LightSimBackend`
- [ADDED] sets of methods to extract the main component of a grid and perform powerflow only on this
  one.
- [ADDED] possibility to set / retrieve the names of each elements of the grid.
- [ADDED] embed in the generator models the "non pv" behaviour. (TODO need to be able to change Q from python side)
- [ADDED] computation of PTPF (Power Transfer Distribution Factor) is now possible
- [ADDED] (not tested) support for more than 2 busbars per substation
- [ADDED] a timer to get the time spent in the gridmodel for the powerflow (env.backend.timer_gridmodel_xx_pf)
  which also include the time 
- [ADDED] support for more than 2 busbars per substation (requires grid2op >= 1.10.0)
- [ADDED] possibility to retrieve the bus id of the original iidm when initializing from pypowsybl 
  (`return_sub_id` kwargs). This is a "beta" feature and will be adressed in a better way
  in a near future.
- [ADDED] possibility to continue the grid2op 'step' when the solver converges but a load or a 
  generator is disconnected from the grid.
- [IMPROVED] now performing the new grid2op `create_test_suite` 
- [IMPROVED] now lightsim2grid properly throw `BackendError`
- [IMPROVED] clean ce cpp side by refactoring: making clearer the difference (linear) solver
  vs powerflow algorithm and move same type of files in the same directory. This change
  does not really affect python side at the moment (but will in future versions)
- [IMPROVED] CI to test on gcc 13 and clang 18 (latest versions to date)
- [IMPROVED] computation speed: grid is not read another time in some cases.
  For example, if load and generators do not change, then Sbus is not
  recomputed. Likewise, if the topology does not change, then the Ybus 
  is not recomputed either see https://github.com/Grid2Op/lightsim2grid/issues/72

[0.7.5.post1] 2024-03-14
-------------------------
- [FIXED] backward compat with "future" grid2op version with a 
  better way to copy `LightSimBackend`
  
[0.7.5] 2023-10-05
--------------------
- [FIXED] a bug in DC powerflow when asking for computation time: it was not reset to 0. when
  multiple powerflows used the same solver
- [FIXED] a bug in AC and DC powerflow when shunts had active values
- [ADDED] possibility to initialize a powergrid based on pypowsybl 
  see https://github.com/Grid2Op/lightsim2grid/issues/53
- [ADDED] some more algorithm to perform powerflow: Fast Decoupled Powerflow (in BX and XB variant)
  see https://github.com/Grid2Op/lightsim2grid/issues/63
- [ADDED] build lightsim2grid for python 3.12
- [ADDED] support for non distributed slack but multiple slack buses
  see https://github.com/Grid2Op/lightsim2grid/issues/50 (ONLY FOR AC powerflow)
- [IMPROVED] now shipping `src` and `eigen` directory in the source of 
  lightsim2grid to allow their installation if wheels are not provided.
- [IMPROVED] in the underlying cpp GridModel powerlines can now have 2
  different values for the `h` parameters (`h_or` and `h_ex`).
- [IMPROVED] now lightsim2grid is able to load a pandapower network with non
  contiguous non starting at 0 bus index

[0.7.3/4] 2023-08-24
--------------------
- [FIXED] a bug where, when you disconnect a load (or gen), the next action cannot be performed
  if it modifies the load (or gen), because you "cannot change the value of a disconnected load (or gen)"
- [FIXED] read-the-docs template is not compatible with latest sphinx version (7.0.0)
  see https://github.com/readthedocs/sphinx_rtd_theme/issues/1463
- [IMPROVED] initialize the underlying "PandaPowerBackend" without numba
- [IMPROVED] grid2op import to be more compliant with renaming of uppercased file names
- [IMPROVED] decoupling of the PandapowerBackend class and the class "internally" used by LightSimBackend
  when loading the grid. This caused some issue, *eg* https://github.com/Grid2Op/grid2op/issues/508

[0.7.2] 2023-06-06
--------------------
- [FIXED] a bug in the `init` function that caused issue when importing a grid with multiple slack
  on some cases
- [FIXED] some bugs in the "SecurityAnalysis" and "TimeSerie" modules especially in DC mode.
- [FIXED] a bug in the DC comptuation: some "divergence" were not catched
- [FIXED] a bug in the "Computer" (cpp) class where the intial voltage could lead to generator not
  participating correctly to the voltage regulation (wrong output voltage level).
- [FIXED] a bug in the "set_bus" of shunt (wrong bus was assigned cpp side)
- [FIXED] an issue when slack bus is added from ext grid (wrong active power value - sign issue)
- [ADDED] support for the CKTSO linear solver (on linux), which is slightly faster than SparseLU, KLU and NICSLU
  (this requires a compilation from source)
- [ADDED] support for distributed slack bus in `LightSimBackend`
- [ADDED] support for "generator with p=0. do not participate in voltage regulation" in `LightSimBackend`
- [ADDED] support for the DC computation for "SecurityAnalysis" and "TimeSerie" modules
- [ADDED] support for DC powerline (in lightsim, they are still not handled in grid2op)
- [IMPROVED] now that multiple slacks is fully supported, the warnings when importing a grid with multiple slacks
  are irrelevant. They have been removed.
- [IMPROVED] the documentation on the "sovlers" part
- [IMPROVED] move the "how to compile" section of the readme in the documentation
- [IMPROVED] `SuiteSparse` is upgraded to version 5.13 (issue with build system based on cmake and BLAS for SuiteSparse >= 6.0)
- [IMPROVED] upgrade to eigen `3.4.0` (stable release)

[0.7.1] 2023-01-11
---------------------
- [BREAKING] drop support for numpy version < 1.20 (to be consistent with grid2op)
- [FIXED] a compatibility issue with grid2op 1.7.2 (missing another backend attribute
  when the environment is copied) see https://github.com/rte-france/Grid2Op/issues/360
- [FIXED] now an error if thrown if the bus indexes in the pandapower grid are not contiguous
  or do not start at 0 (thanks Roman Bolgaryn for spotting this issue)
- [ADDED] automatic build for python 3.11
- [ADDED] support for numpy >= 1.24 (some deprecation *eg** np.str and np.bool are removed)

[0.7.0.post1] 2022-06-20
-------------------------
- [FIXED] a compatibility issue with grid2op 1.7.1 (missing a backend attribute
  when environment is copied)

[0.7.0] 2022-05-30
---------------------
- [ADDED] improved time measurments
- [ADDED] Possibility to set, at creation time, the type of solver used, number
  of iterations and precisions with 
  `LightSimBackend(max_iter=..., tol=..., solver_type=...)`
- [IMPROVED] scripts to load the pandapower grid (json format)
- [IMPROVED] update the automatic tests on more recent compilers.

[0.6.1.post2] 2022-02-08
-------------------------
- [FIXED] add support for python 3.10 now that scipy does (and add proper tests in CI)

[0.6.1.post1] 2022-02-02
-------------------------
- [FIXED] support for python3.7 (and add proper tests in CI)

[0.6.1] 2022-02-01
--------------------
- [BREAKING] the behaviour of the `newton_pf` function is not 
  consistent with pandapower default concerning distributed slack.
- [FIXED] an issue in the distributed slack case spotted by pandapower team 
  thanks to them (see https://github.com/e2nIEE/pandapower/pull/1455)
- [IMPROVED] lightsim2grid will now use the single slack algorithm if the 
  grids counts only one slack bus (performance increase)

[0.6.0] 2021-12-17
-------------------
- [BREAKING] change the interface of the `newton_pf` function to reflect pandapower change in their
  latest version (arguments `ref` has been added). You can still use the old `newton_pf` function, with the
  old signature by importing `newtonpf_old` instead or explicitly importing the new one by importing `newtonpf_new`
- [BREAKING] `SecurityAnalysis` now also returns the active flows when calling `security_analysis.get_flows()`
- [BREAKING] change the file names (python side) to be compliant with pep 8. You can no longer
  do things like `from lightsim2grid.LightSimBackend import LightSimBackend` change it to
  `from lightsim2grid import LightSimBackend` (preferred method)
- [BREAKING] change the file names (python side) to be compliant with pep 8. You can no longer
  do things like `from lightsim2grid.initGridModel import init` change it to
  `from lightsim2grid.gridmodel import init` (preferred method) (same for `GridModel` class)
- [FIXED] a bug that lead to the wrong computation of the dc powerflow in case of `sn_mva != 1.` and phase shifters.
- [FIXED] bug preventing to use the NICSLU linear solver in the `GridModel`
- [FIXED] compilation warnings on clang (missing virtual destructor, unused variables, etc.)
- [FIXED] a bug in the `SecurityAnalysisCPP`: when it diverges for some contingencies, the others were not simulated properly.
- [FIXED] `LightSimBackend` now contains members for `shunts` and `***_theta` as it does for the other quantities. This improves the consistency, but most importantly
  fixes some bugs when used in earlier grid2op versions
- [ADDED] possibility to compute the active flows using the `BaseMultiplePower` 
- [ADDED] possibility to change linear solver used when performing a DC solver
- [ADDED] possibility to make powerflow with distributed slack bus (only for newton raphson at the moment)
- [ADDED] access (read only) to the element of a lightsim2grid grid with the `get_XXX` (*eg* `get_loads()`) methods (see documentation)
- [ADDED] direct access to the solver used in the grid model python side
- [ADDED] unittest in circleci.
- [ADDED] all kind of solvers based on different linear solvers (Eigen sparse LU, KLU or NICSLU) for Newton Raphson and
  DC approximation (9 solvers in total)
- [IMPROVED] use of `steady_clock` to retrieve the ellapse time c++ side
- [IMPROVED] refactoring of the c++ part to use template mecanism instead of inheritance for the
  Newton Raphson and DC solvers.
- [IMPROVED] `GridModel` now contains two different solvers, one for AC powerflow and one for DC powerflow.
- [IMPROVED] error message in the solver are now embedded in an Enum instead of being integers, for better readibility.
- [IMPROVED] error message when the powerflow diverge (error are read from c++ now)

[0.5.5] 2021-11-10
-------------------
- [ADDED] possibility to perform dc powerflow
- [ADDED] a class to compute flows on whole time series when the Ybus does not change (see `TimeSerie`)
- [ADDED] a class to compute flows on multiple contingencies, when Sbus does not change (see `SecurityAnalysis`).
- [IMPROVED] running speed of Newton Raphson solvers with better filling of sparse matrices
- [IMPROVED] upgrade to SuiteSparse `v5.10.1`
- [IMPROVED] upgrade to eigen `3.4.0` (stable release)
- [IMPROVED] clean the compilation warnings on microsoft windows (force the conversion from
  `Eigen::EigenBase<Derived>::Index` to `int` using `static_cast`)
- [IMPROVED] add the proper optimization flag for windows (`/O2` instead of `-03` on linux / macos)
- [IMPROVED] high performance gain when topology is not changed between steps (gain obtained by 
  reusing the previous Ybus)

[0.5.4] 2021-08-20
------------------
- [FIXED] a bug for static generator (wrong signed convention were used in some part of the c++ code). This has
  no impact at all for provided grid2op environments.
- [FIXED] An issue where the backend could get "stuck" in a wrong state because of the way the Vinit was computed (see
  `Issue 30 <https://github.com/Grid2Op/lightsim2grid/issues/30>`_)
- [ADDED] experimental support for the `NICSLU` linear solver (requires a proper license and library, see
  https://github.com/chenxm1986/nicslu for more information. Support does not include multi threaded at the moment).
- [IMPROVED] minor performance improvements for the solvers based on Newton Raphson (faster filling of the Jacobian
  matrix after the first iteration)

[0.5.3] 2021-08-11
-------------------
- [FIXED] minor issues in the benchmark (some time measurments were wrong)
- [ADDED] lightsim2grid package now can be distributed on pypi
- [ADDED] compilation of SuiteSparse using cmake
- [ADDED] compatibility with the KLU linear solver on windows based systems.
- [IMPROVED] the package should now be available on pypi

[0.5.2] 2021-07-26
-------------------
- [FIXED] `GridModel` now properly throw "out_of_range" exception when trying to change the bus of non existing
  elements
- [FIXED] wrong units were displayed for the iterators for lines and transformers.
- [ADDED] now able to retrieve the powerlines parameters python side.
- [IMPROVED] more explicit error messages when the building of the `Ybus` matrix fails.
- [IMPROVED] now the solver is not reset when using the `backend._grid.check_solution`
- [IMPROVED] upgrade SuiteSparse to version `v5.10.1`
- [IMPROVED] upgrade eigen to version `3.4-rc1`

[0.5.1] 2021-04-09
-------------------
- [FIXED] yet another compilation issue with clang (see
  `Issue 22 <https://github.com/Grid2Op/lightsim2grid/issues/22>`_)
- [ADDED] circleci to check compilation for gcc
- [ADDED] circleci to check compilation for clang
- [ADDED] circleci to check compilation for msvc
- [ADDED] function to read the voltage angle from the backend
- [ADDED] compatibility with grid2op 1.5.0 (up to an issue with the storage units)

[0.5.0] 2021-03-01
-------------------
- [FIXED] a compilation issue on macos
- [FIXED] a compilation issue on windows (missing import of vector in `DataConverter.h`)
- [FIXED] an import issue (with `lightsim2grid.SolverType`)
- [FIXED] a bug that lead to the wrong computation of the ratio of the trafo when the tap on hv side.
- [FIXED] wrong timing was measured in the "solver powerflow time" of pandapower in the benchmarks
- [FIXED] a broken handling of shunt modification (wrong bus was assigned)
- [FIXED] an issue in `LightSimBackend.copy` that prevent the copied environment from being reset.
- [FIXED] errors are now raised when pandapower grid cannot be converted in lightsim2grid (*eg.* when
  unsupported elements are present)
- [ADDED] a variant of the Gauss Seidel method which does the update in a "synchronous" fashion
- [ADDED] a function that, given a complex vector is able to check kicchoff's law violation.
- [ADDED] Support for phase shifter (modeled as trafo with an extra parameter `shift`)
- [ADDED] Experimental support for `sn_mva` pandapower parameter.
- [UPDATED] github issue template
- [IMPROVED] warnings are issued when some of the pandapowergrid attributes have been automatically replaced
  when converting to / from pandapower

[0.4.0] - 2020-10-26
---------------------
- [ADDED] the Gauss Seidel method for AC powerflow is now available
- [ADDED] possibility to change easily the solver types from python side

[0.3.0] - 2020-10-06
-------------------------
- [ADDED] Support for pickle for the lightsim Backend.
- [ADDED] LightSim should now be compatible with windows (implementation of a powerflow mode without
  using the SuiteSparse KLU linear solver but rather the Eigen SparseLU one)
- [ADDED] start of the documentation.

[0.2.4] - 2020-08-20
--------------------
- [FIXED] issue for copying environment

[0.2.3] - 2020-08-03
--------------------
- [UPDATED] consistent behaviour between grid2op.PandaPowerBackend and LightSimBackend for action that
  set the bus of only one extremity of a powerline.
- [ADDED] compatibility with grid2op 1.2.0

[0.2.2] - 2020-06-25
---------------------
- [UPDATED] removing the `-march=native` that causes some difficulty for some compilers
- [ADDED] compatibility with grid2op 1.0.0

[0.2.1] - 2020-06-xx
--------------------
- [FIXED] update of the `topo_vect` attribute in class `LightSimBackend` when reset.
- [ADDED] a github issue template

[0.2.0] - 2020-06-15
--------------------
- [ADDED] the changelog
- [FIXED] the import of files when elements where not in service
- [FIXED] a bad catch of a divergence in the solver
- [IMPROVED] the speed to apply the actions
- [FIXED] tests for the backend in grid2op and here are not identical without (too much) duplicates
