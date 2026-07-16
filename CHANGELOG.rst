Change Log
===========

[TODO]
--------
- [refacto] have a structure in cpp for the buses
- [refacto] have the id_grid_to_solver and id_solver_to_grid etc. directly in the solver and NOT in the gridmodel.
- [refacto] put some method in the DataGeneric as well as some attribute (_status for example)
- support 3w trafo (as modeled in pandapower)
- improve speed by not performing internal checks 
  (keep check for boundaries and all for python API instead) [see `TODO DEBUG MODE` in c++ code]
- improve speed
- code parrallelism directly in the `Computer` and `SecurityAnalysisCPP` classes
- a mode to do both `Computer` and `SecurityAnalysisCPP`
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

[0.14.0] 2026-xx-yy
---------------------
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
- [IMPROVED] cleaner `cktso_lib` (`from lightsim2grid.compilation_options import cktso_lib`) : the file name and extension are omitted
- [IMPROVED] easier build by relying on cmake and scikit_build_core to build the cpp part
- [IMPROVED] SuiteSparse to version 7.12.2 (2026-02-05)

[0.12.2] 2026-02-05
----------------------
- [FIXED] an issue with shunt buses (was set to 1 even if they were disconnected)
- [FIXED] a warning when applying actions on generator votlage setpoints (due to NaN)
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
