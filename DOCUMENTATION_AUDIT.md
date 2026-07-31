# Documentation audit — C++ core / Python bindings (prep for 1.0.0)

Scope: `src/core/help_fun_msg.hpp` + `src/core/help_fun_msg.cpp` (the shared
docstring library) and `src/bindings/python/*.cpp` (the pybind11 bindings that
consume it), on branch `n1_full_compute`. Goal, per the request: (1) every
public, actually-used method should have real documentation instead of
`"TODO"` / "internal, do not use", (2) documentation should live in
`help_fun_msg.hpp/.cpp` rather than being duplicated/hand-written per binding,
(3) flag documentation whose content no longer matches the current code.

**Status: all 18 findings in section 3 (accuracy bugs) have been fixed** —
see the commit(s) following this audit. Sections 1, 2, 4 and 5 (placement,
placeholders, missing docs, recommendations) are still just findings, not yet
acted on.

Sphinx pulls these docstrings directly (`autoclass`/`automodule` in
`docs/*.rst`), so anything wrong here is wrong on the public docs site too —
this isn't just an in-REPL `help()` problem.

## 1. Quantitative overview

`help_fun_msg.hpp` declares 217 doc strings across 5 structs (all 217 are
defined in the `.cpp`, none dangling either way):

| Struct               | # entries |
|----------------------|-----------|
| `DocIterator`        | 82        |
| `DocLSGrid`          | 63        |
| `DocSolver`          | 39        |
| `DocSecurityAnalysis`| 18        |
| `DocComputers`       | 15        |

Placeholder / non-informative docstrings directly in the binding files:

| File                       | `"TODO`-ish strings | `_internal_do_not_use` | total bindings (`.def*`) |
|----------------------------|:---:|:---:|:---:|
| `binding_lsgrid.cpp`       | 24  | 155 | ~306 |
| `binding_containers.cpp`   | 47  | 0   | ~302 (incl. `def_readonly`) |
| `binding_misc.cpp`         | 13  | 0   | 21   |
| `binding_solvers.cpp`      | 3   | 6   | ~127 |
| `binding_batch.cpp`        | 2   | 1   | 73   |

`binding_lsgrid.cpp` alone tags roughly half of `LSGrid`'s bound methods with
the single generic `DocLSGrid::_internal_do_not_use` string ("Internal, do
not use unless you know what you're doing"). That string is reused verbatim
~155 times.

**That "internal" label is not accurate for most of them.** Cross-referencing
those 156 distinct method names against the actual `lightsim2grid` Python
package (everything under `lightsim2grid/`, excluding `tests/`) shows:

- **81 of the 156** are called directly by the package itself (mostly
  `lightSimBackend.py` and `lightsim2grid/network/from_*`) — e.g.
  `get_bus_vn_kv`, `get_gen_res_full`, `get_line_res1_full`,
  `get_trafo_res1_full`, `init_bus`, `init_powerlines`, `init_generators`,
  `update_gens_p`, `update_loads_p`, `update_topo`, `change_p_load`,
  `deactivate_gen`, `set_substation_names`, `unset_changes`, `set_n_sub`,
  `set_line_pos1_topo_vect`, etc. These are exactly the "public methods
  actually used" the audit is meant to find — they need real docstrings,
  not the internal-do-not-use boilerplate.
- The remaining 75 are not called from the package outside tests (e.g. the
  new `_side1`/`_side2` half-open-line accessors, `debug_get_Bp_python`,
  `is_grid_connected_after_contingency`'s low-level twin, several
  `get_bus_*`/`change_bus_*` per-element accessors). Some of these look
  genuinely internal/plumbing; others (e.g. `get_bus1_powerline` /
  `get_bus2_powerline`, the counterparts of `change_bus1_powerline` which
  *is* used) are simply the read side of a used read/write pair and are
  probably worth real docs too, since they're reachable public API either
  way.

Full list of the 156 method names checked, split used/not-used, is at the end
of this document (§6) for reference when doing the actual rewrite.

## 2. Structural finding: documentation is not uniformly centralized

The stated goal is "every doc lives in `help_fun_msg`". Today that's true for
a majority of `LSGrid`, `DocIterator`-covered element fields, and the solver
family, but several whole classes/areas are documented **only** inline in the
binding `.cpp` files, which is exactly what would need to be duplicated again
for a second (C, Rust, …) binding:

- `TimerJac` and `LinearSolverStats` (`binding_solvers.cpp`) — full class doc
  + all `def_readonly` fields, no `Doc*` entry at all.
- `bind_nr_algo_policies` / `bind_linear_solver_stats` /
  `bind_fdpf_linear_solver_stats` template helpers in `binding_solvers.cpp` —
  ~25 methods (scaling/refactor policy get/set, `get_config`/`set_config`,
  `get_theta_to_J_col`, etc.), all with good inline prose, none in
  `help_fun_msg`.
- `AlgoConfig`, `AlgoControl` (`binding_misc.cpp`) — `AlgoControl` is 13
  methods, all `"TODO"`.
- `PandaPowerConverter` (`binding_misc.cpp`) — 6 methods, **zero** docstrings
  (not even a placeholder string is passed).
- `SvcContainer`/`SvcInfo`, `ConverterStationInfo` (`binding_containers.cpp`)
  — entirely inline strings, no `DocIterator` entries exist for SVC-specific
  fields at all (regulation mode, slope, susceptance limits, target
  setpoints) even though every other element family has one.
- `LimitViolation`, `ViolationElementType`, `LimitViolationType`
  (`binding_batch.cpp`) — inline only.
- `add_pickle`/`add_binary_serialization` (`pickle_helpers.hpp`,
  `binary_helpers.hpp`) — `save_binary`/`load_binary` have solid inline docs,
  but they live in a template helper header, not `help_fun_msg`; low priority
  since it's centralized in one place already (not copy-pasted per class),
  but still outside the stated single-source-of-truth file.
- ~40 individual methods in `binding_lsgrid.cpp` have good, specific inline
  prose that never made it into `help_fun_msg` (`get_reference_slack_bus`,
  `set_ignore_status_global`, `change_shift_trafo`,
  `set_trafo_shift_dependent_rx`, `get_slack_absorbed_solver`,
  `get_controller_q_solver`, `get_p_buses_solver` / `get_q_buses_solver` /
  `get_theta_buses_solver` / `get_vm_buses_solver`,
  `get_hvdc_droop_data_solver`, `get_status_droop_hvdc*`, and others). These
  are actually **better written** than many `Doc*` entries — moving them
  as-is into `help_fun_msg` is mostly mechanical.

## 3. Confirmed accuracy bugs (docstring contradicts current code)

These are verified against the current implementation, not just suspected —
each is either wired to the wrong `Doc*` constant, self-contradictory, or
references a name that no longer exists.

1. **[FIXED]** **`binding_solvers.cpp:382`** — `AlgorithmSelector.get_error` is bound
   with `DocSolver::get_V.c_str()` (copy-paste from the line above/below);
   it should be `DocSolver::get_error`. Right now calling `help()` on
   `get_error` describes "the complex voltage for each bus."

2. **[FIXED]** **`binding_lsgrid.cpp:102`** — `LSGrid.available_default_algorithms` is
   bound with `DocLSGrid::available_algorithm_names` ("a list of string"),
   while the doc that actually matches what it returns —
   `DocLSGrid::available_default_algorithms` ("a list of
   `lightsim2grid.algorithm.AlgorithmType`") — sits unused right next to it
   in `help_fun_msg.cpp:1900-1913` and is instead attached to nothing. The
   two docstrings exist, they're just cross-wired.

3. **[FIXED]** **`binding_containers.cpp:390`** — `HvdcLineInfo.p2_mw` is documented
   with `DocIterator::target_p_or_mw`, which describes the **side-1/"or"**
   target and its sign convention ("if positive, power flows from ex to
   or"). `p2_mw` is `station_side_2.target_p_mw`
   (`HvdcLineContainer.hpp:32,479`) — the side-2 value, wrong side entirely,
   and the described sign convention doesn't apply to it.

4. **[FIXED]** **`binding_containers.cpp:431-432`** — `SubstationInfo.nb_max_busbars`
   and `SubstationInfo.vn_kv` are both bound with `DocIterator::name` (a
   string-field doc). `nb_max_busbars` is an integer (max busbars per
   substation), `vn_kv` is the substation's nominal voltage in kV — neither
   has anything to do with "name."

5. **[FIXED]** **`binding_batch.cpp:60,106,115,123,133`** — five `def_property`
   docstrings (`TimeSeries.init_from_n_powerflow`,
   `ContingencyAnalysis.compute_limit_violations`, `.init_from_n_powerflow`,
   `.handle_disconnected_grid`, `.nb_thread`) use raw-string literals
   (`R"mydelim(...)"`) that still contain the leftover `"` characters and
   mid-string indentation from what used to be ordinary concatenated string
   literals (`"a" "b" "c"`). Converting to a raw string without stripping
   those artifacts means the *actual* compiled docstring shown to users
   contains stray quote marks and ragged indentation, e.g. (verbatim):
   ```
   Whether to initialize the complex voltages of "
                         "the first time series with the results of a n-powerflow "
                         "(*ie* a powerflow at the start the simulation) or not. "
                         "Default: false
   ```

6. **[FIXED]** **`help_fun_msg.cpp` — all 20 `DocSolver::<Family>` solver descriptions**
   (`NR_SparseLU`, `NRSing_SparseLU`, `DC_SparseLU`, `FDPF_XB_SparseLU`,
   `FDPF_BX_SparseLU`, and the `_KLU`/`_NICSLU`/`_CKTSO` equivalents, lines
   ~158-670) each contain a line like *"In the enum
   `lightsim2grid.solver.AlgorithmType`, it is called `SparseLU`"* (or
   `KLU`, `KLUSingleSlack`, `KLUDC`, `FDPF_KLU`, `NICSLUDC`, `CKTSODC`, …).
   None of those short names exist anymore: `binding_enums.cpp` shows the
   *current* `AlgorithmType` values are `NR_SparseLU`, `NRSing_SparseLU`,
   `DC_SparseLU`, `FDPF_XB_SparseLU`/`FDPF_BX_SparseLU`, `NR_KLU`, … — the
   old short names are literally present in the source only as **commented
   out** deprecated aliases that were never re-enabled. Every one of these
   20 docstrings is telling users to look up an enum member that doesn't
   exist.

7. **[FIXED]** **`help_fun_msg.cpp` — all of `DocComputers` (13 occurrences) and
   `DocSecurityAnalysis` (55 occurrences)** reference the pre-rename Python
   names `lightsim2grid.timeSerie.Computers` and
   `lightsim2grid.securityAnalysis.SecurityAnalysisCPP` /
   `SecurityAnalysis`. Both are now deprecation shims
   (`lightsim2grid/securityAnalysis.py` does
   `warnings.warn("...upgrade to SecurityAnalysisCPP > ContingencyAnalysisCPP...")`
   and just re-exports); the current names are
   `lightsim2grid.timeSerie.TimeSerie` / `TimeSeriesCPP` and
   `lightsim2grid.contingencyAnalysis.ContingencyAnalysis` /
   `ContingencyAnalysisCPP`. This affects essentially the entire
   `DocComputers`/`DocSecurityAnalysis` struct — every cross-reference in it
   points at deprecated names. Given this, renaming the structs themselves
   (`DocComputers` → something TimeSeries-named, `DocSecurityAnalysis` →
   `DocContingencyAnalysis`) as part of the uniformity pass would remove the
   temptation to keep writing "security analysis" in new docs too.

8. **[FIXED]** **`help_fun_msg.cpp:3406-3419` — `DocComputers::get_power_flows`** opens
   with "Get the **active** flows (in **MW**)…" but its own `Returns`
   section says `"The flows (in kA)…"` — copy-pasted from
   `DocComputers::get_flows` a few lines above without updating the unit.
   Self-contradictory within the same docstring.

9. **[FIXED]** **`help_fun_msg.cpp:77` — `DocSolver::get_error`** ships literally as:
   `"Returns the error code (as an integer) encountered by the solver (0 =
   no error). TODO DOC explain better"` — the placeholder text is in the
   compiled, user-facing docstring, not a source comment.

Given the density of issues in items 6-7 (every entry in three of the five
structs affected), it's likely that a full pass would turn up more of the
same pattern once someone starts rewriting — these should be treated as
"the whole struct needs a rename-aware pass," not as three isolated typos.

Deeper cross-checks against the current `NRSystem`/`BaseDCAlgo`/element
implementations (done with fresh, independent re-reads of the code, not
just the bindings) turned up more of the same class of problem,
concentrated in `DocLSGrid::J_description`, `DocIterator`'s hvdc-line block,
and a couple of sign-convention docs that flatly contradict each other:

10. **[FIXED]** **`DocLSGrid::J_description` (`help_fun_msg.cpp:2115-2151`) describes a
    Jacobian shape that predates the current solver architecture.** It
    documents J as a fixed 2×2 block plus one extra row/column for the
    distributed slack. The actual builder is a composable
    `NRSystem<Base, Rest...>` (`powerflow_algorithm/NRSystem.hpp:1238-1239`)
    where `MultiSlack` adds `nb_slack` rows/columns (including a dedicated
    "slack_absorbed" column, not a single row), `VoltageControl` adds a
    per-regulation-group block of reactive-power unknowns/equations, and
    `Hvdc` adds droop-flow derivative entries for angle-droop HVDC lines —
    none of which are mentioned. This propagates into
    `DocLSGrid::get_J_python_solver` (`help_fun_msg.cpp:2221-2240`), which
    appends `J_description` verbatim.

11. **[FIXED]** **The same `J_description` claims "all buses in `ref` will be pv buses
    except the first one… this cannot be changed at the moment."** That's
    no longer true: `NRSystem.hpp`'s `Base::register_in` (~line 210-233)
    explicitly supports slack buses that are *not* locally PV-pinned (a
    `free_vm_slack_buses_` set, e.g. a slack bus that's only remotely/SVC
    regulated), giving them a free Vm unknown like an ordinary PQ bus.

12. **[FIXED]** **`DocLSGrid::get_slack_ids_dc` / `get_slack_ids_dc_solver`
    (`help_fun_msg.cpp:2493`, `2605`) say the distributed slack is
    "currently not taken into account" for DC.** `BaseDCAlgo::compute_pf_dc`
    (`powerflow_algorithm/BaseDCAlgo.tpp:108-132`) explicitly spreads the
    active-power imbalance across all participating slack buses
    proportionally to `slack_weights` — i.e. distributed slack for the DC
    powerflow itself *is* implemented today. (The "no distributed slack"
    caveat is still accurate, separately, for `get_ptdf`/`get_lodf`
    specifically — `BaseDCAlgo.tpp:441` still has a live
    `// TODO PTDF: distributed slack` — so don't just delete the caveat
    everywhere, narrow it to the two functions where it still applies.)

13. **[FIXED]** **`DocSolver::get_timers` / `DocSolver::get_computation_time`
    (`help_fun_msg.cpp:137-156`, `~719-723`) name a return field,
    `timer_total_`, that doesn't exist.** `BaseAlgo::get_timers()`
    (`BaseAlgo.hpp:314-318`) returns the 4-tuple `(timer_Fx_, timer_solve_,
    timer_check_, timer_total_nr_)` — the last element is `timer_total_nr_`.
    (The other three names and the ordering are correct; don't confuse this
    with the separate, correctly-documented 13-field `TimerJac` struct
    returned by `get_timers_jacobian()`.)

14. **[FIXED]** **`DocIterator::bus_ex_id` (`help_fun_msg.cpp:1525-1529`) still says
    "lv" instead of "ex."** Used for `LineInfo.bus2_id`
    (`binding_containers.cpp:299`) and `HvdcLineInfo.bus2_id` (`:388`), the
    text reads *"...at which the 'lv' side of the line is connected"* —
    lines/hvdc-lines have `or`/`ex` sides, not `hv`/`lv` (that's
    transformer-only vocabulary). Looks like a copy from `bus_lv_id`'s
    docstring where "transformer" → "line" was updated but "lv" → "ex" was
    missed.

15. **[FIXED]** **`DocIterator::dc_line_formula` and `DocIterator::target_p_or_mw`
    directly contradict each other on sign convention, on the same
    class.** `dc_line_formula` (`help_fun_msg.cpp:1714-1738`) says the hvdc
    line uses "load convention": positive `or_mw` means power is
    *consumed* at `or` (flows `or`→`ex`). `target_p_or_mw`
    (`:1695-1702`, also used on `HvdcLineInfo`) says positive means power
    is *injected* at `or` (flows `ex`→`or`) — the opposite. The source
    comment in `HvdcLineContainer.hpp:29-30` confirms **generator**
    convention, so `dc_line_formula`'s "load convention" claim is the
    backwards one.

16. **[FIXED]** **`DocIterator::DCLineContainer`/`DCLineInfo` (`help_fun_msg.cpp:1594-1638`,
    `:1641-1693`) describe a retired model, and the example code in them is
    broken.** They still describe the old pandapower-only "dc line = 2
    generators" model; none of the current IIDM converter-station
    (VSC/LCC), angle-droop "AC emulation," or resistive-loss
    (`r_ohm`/`nominal_v_kv`) behavior is mentioned
    (`HvdcLineContainer.hpp:59-90`). Worse, the docstring's own usage
    example calls `dcline.bus_or_id`, which doesn't exist — the real
    attribute is `bus1_id` (`binding_containers.cpp:387`) — so the example
    would raise `AttributeError` if a user copy-pasted it.

17. **[FIXED]** **`DocIterator::min_p_mw`/`max_p_mw` and `min_q_mvar`/`max_q_mvar`
    misdescribe how (and whether) they're used.** `min_p_mw`/`max_p_mw`
    (`help_fun_msg.cpp:1027-1043`, only bound on `SGenInfo`) claim they're
    "used in `check_solution` if `check_q_limits` is set to True" — but
    `LSGrid::check_solution_q_values` (`LSGrid.cpp:1080-1139`) only ever
    checks **reactive** limits and never iterates static generators at all;
    these fields are effectively unused today. `min_q_mvar`/`max_q_mvar`
    (`:1009-1025`) claim the opposite extreme — "not taken into account by
    the solver" — but for `GenInfo`/`ConverterStationInfo` they actively
    drive the runtime Q-mismatch split across multiple PV
    generators/converters on the same bus
    (`GeneratorContainer.cpp:512-593`), changing the reported `res_q_mvar`;
    only for `SGenInfo` is "not used" actually accurate.

18. **[FIXED]** **`DocIterator::is_slack` (`help_fun_msg.cpp:988-999`) references a
    renamed solver class**, `lightsim2grid.solver.SparseLUSingleSlack` —
    same class-of-bug as finding 6 above; current name is `NRSing_SparseLU`
    (`binding_enums.cpp:30`, old name commented out at `:56`).

Corroboration note: two independent research passes over `DocIterator` and
`DocLSGrid` both flagged the exact same `AlgorithmType` rename issue (finding
6) and the same `HvdcLineInfo.p2_mw` mislabeling (finding 3) that were found
manually — increasing confidence in both.

**[FIXED] 19. Follow-up: the "or"/"ex" → "1"/"2" rename wasn't actually complete.**
Findings 3, 15 and 16 only fixed the hvdc-specific corner of a wider problem:
`LineInfo`'s own fields (`bus_or_id`/`bus_ex_id`, `res_p_or_mw`/`res_p_ex_mw`,
`res_q_or_mvar`/`res_q_ex_mvar`, `res_v_or_kv`/`res_v_ex_kv`,
`res_a_or_ka`/`res_a_ex_ka`, `res_theta_or_deg`/`res_theta_ex_deg`) were still
named and worded after the retired "or" (origin) / "ex" (extremity)
convention, even though the bound Python attributes are `bus1_id`/`bus2_id`,
`res_p1_mw`/`res_p2_mw`, etc. — the same mismatch as the hvdc case, just not
flagged the first time because nothing there was *factually wrong* (unlike
`p2_mw`), only stale terminology. The shared `line_model` ASCII schema
(`DocIterator::line_model`, used by `r_pu`/`x_pu`/`h_pu` for both lines and
transformers) had the same problem: it labelled its diagram `or bus`/`ex bus`,
`ior`/`iex`, `vor`/`vex`. Fixed by renaming every affected `DocIterator`
constant to the `_1_`/`_2_` convention (`bus_or_id` → `bus_1_id`,
`res_p_or_mw` → `res_p_1_mw`, the hvdc `_dcline`-suffixed ones too, for
naming consistency with the rest of the file) and rewriting the schema to use
`bus 1`/`bus 2`, `i1`/`i2`, `v1`/`v2`, with a note that transformers keep
`hv`/`lv` instead (a real electrical distinction, not just a naming choice —
left untouched).

Areas that *were* checked closely and found to still be accurate, so they
don't need touching in the accuracy pass (only the reorganization):
`check_grid`, `change_algorithm`'s DC/AC routing, `check_solution`'s KCL-
mismatch description, `ac_pf`/`dc_pf`'s interface description, `get_ptdf`/
`get_lodf`/`get_Bf` (shapes, row ordering, gridmodel-vs-solver labelling,
and their own distributed-slack caveat), the bus-labelling family
(`get_Va`/`Vm`/`V`, `id_*_to_*`, `get_pv`/`pq`/`slack_ids`/`slack_weights`,
`get_Ybus`/`get_Sbus` and their `_solver`/`dc` variants), and the
hv/lv-side reuse for `TrafoInfo` (unlike the or/ex reuse bug above, the
hv/lv one is actually correct).

## 4. Missing documentation entirely

- `PandaPowerConverter` (`binding_misc.cpp:17-24`): `set_f_hz`, `set_sn_mva`,
  `get_line_param_legacy`, `get_line_param`, `get_trafo_param_pp3`,
  `get_trafo_param_pp2` — no docstring argument passed at all (not even a
  placeholder).
- `AlgoControl` (`binding_misc.cpp:36-49`): class + all 13 methods are
  `"TODO"`.
- Many `TrafoInfo`/`LineInfo` fields — the raw/effective admittance
  components (`yac_11`, `yac_12`, `yac_21`, `yac_22`, `yac_eff_*`, `ydc_*`)
  are all `"TODO doc"` (8 fields × 2 classes = 16 placeholders), as are the
  matching `get_yac_eff_*` accessor methods on the containers.
- `get_bus_id`/`get_bus_id_side_1`/`get_bus_id_side_2` on every
  `*Container` class (Generator/SGen/Load/Storage/Shunt/Trafo/Line) —
  `"TODO doc"`.
- `SubstationContainer`/`SubstationInfo` class docstrings are themselves
  `"TODO"` (separately from the field-level bug in §3.4).

## 5. Recommendation for the rewrite pass

1. Treat §3 (accuracy bugs, 18 items) as the first, independent fix — these
   are outright wrong today regardless of the reorganization, and most (the
   swapped-`Doc*`-constant bugs, the enum-rename bugs, the raw-string
   quote artifacts, the self-contradicting units) are cheap, mechanical
   fixes safe to do without touching the uniformity work. Three are not
   mechanical and need someone who actually knows the current solver
   internals to rewrite, not just patch: `J_description` (items 10-11 —
   needs a description of the `NRSystem<Base, MultiSlack, VoltageControl,
   Hvdc>` composable layout, not a touch-up of the old 2×2-block text),
   the DC distributed-slack claim (item 12 — narrow the caveat to
   PTDF/LODF only, don't just delete it), and the hvdc-line family (items
   15-16 — the model description and its example code need a full rewrite
   for the current IIDM/converter-station/angle-droop model, not an edit).
2. Use the §1 used/unused split to prioritize which of the 156
   `_internal_do_not_use` methods get real prose first: the 81 confirmed
   used by `lightSimBackend`/`network/from_*` first, then decide case by
   case on the remaining 75 (some are genuinely internal setup calls only
   ever invoked once from `gridmodel.py`'s `init_*` sequence and can
   reasonably stay terse/internal; others, like the `_side1`/`_side2`
   half-open accessors, are newer public API for this same release and
   deserve full docs, not "internal").
3. When moving the ~40 well-written inline `binding_lsgrid.cpp` docstrings
   and the solver-policy/`TimerJac`/`LinearSolverStats`/`AlgoConfig` docs
   into `help_fun_msg`, this is mostly a cut-and-paste job (content is
   already good) — the main risk is merge conflicts with any concurrent
   work on `n1_full_compute`, not writing new prose.
4. For `DocComputers`/`DocSecurityAnalysis`, decide up front whether to
   rename the C++ structs (`DocContingencyAnalysis`) to match the
   already-renamed Python/C++ classes, since every string in them needs
   touching anyway for the name fix.
5. Add an SVC-specific block to `DocIterator` (currently entirely absent)
   so `SvcContainer`/`SvcInfo` stop being the one element family documented
   ad hoc.

## 6. Appendix — `_internal_do_not_use`-tagged methods, used vs. not

**Used elsewhere in `lightsim2grid/` (outside `tests/`) — 81, prioritize
these:**
`add_gen_slackbus`, `change_bus_gen`, `change_bus_load`, `change_bus_shunt`,
`change_p_gen`, `change_p_load`, `change_p_shunt`, `change_q_load`,
`change_q_shunt`, `change_v_gen`, `deactivate_bus`, `deactivate_dcline`,
`deactivate_gen`, `deactivate_load`, `deactivate_powerline`,
`deactivate_sgen`, `deactivate_shunt`, `deactivate_storage`,
`deactivate_svc`, `deactivate_trafo`, `get_bus_vn_kv`, `get_gen_res_full`,
`get_init_vm_pu`, `get_line_res1_full`, `get_line_res2_full`,
`get_lines_status`, `get_loads_res_full`, `get_n_sub`,
`get_shunts_res_full`, `get_storage_target_p`, `get_storages_res_full`,
`get_substation_names`, `get_trafo_res1_full`, `get_trafo_res2_full`,
`get_trafo_status`, `init_bus`, `init_bus_status`, `init_dclines`,
`init_generators`, `init_generators_full`, `init_hvdc_lines`, `init_loads`,
`init_powerlines`, `init_powerlines_full`, `init_sgens`, `init_shunt`,
`init_storages`, `init_svcs`, `init_trafo`, `init_trafo_pandapower`,
`reactivate_bus`, `reactivate_gen`, `reactivate_load`, `reactivate_shunt`,
`set_gen_pos_topo_vect`, `set_gen_to_subid`, `set_init_vm_pu`,
`set_line_pos1_topo_vect`, `set_line_pos2_topo_vect`,
`set_line_to_sub1_id`, `set_line_to_sub2_id`, `set_load_pos_topo_vect`,
`set_load_to_subid`, `set_n_sub`, `set_shunt_to_subid`, `set_sn_mva`,
`set_storage_pos_topo_vect`, `set_storage_to_subid`,
`set_substation_names`, `set_trafo_pos1_topo_vect`,
`set_trafo_pos2_topo_vect`, `set_trafo_to_sub1_id`, `set_trafo_to_sub2_id`,
`tell_solver_need_reset`, `unset_changes`, `update_gens_p`, `update_gens_v`,
`update_loads_p`, `update_loads_q`, `update_storages_p`, `update_topo`.

**Not found outside `tests/`/bindings — 75, review case by case:**
`change_bus1_dcline`, `change_bus1_powerline`, `change_bus1_trafo`,
`change_bus2_dcline`, `change_bus2_powerline`, `change_bus2_trafo`,
`change_bus_sgen`, `change_bus_storage`, `change_bus_svc`,
`change_p_dcline`, `change_p_sgen`, `change_p_storage`, `change_q_sgen`,
`change_q_storage`, `change_v1_dcline`, `change_v2_dcline`,
`debug_get_Bp_python`, `debug_get_Bpp_python`, `get_all_shunt_buses`,
`get_bus1_dcline`, `get_bus1_powerline`, `get_bus1_trafo`,
`get_bus2_dcline`, `get_bus2_powerline`, `get_bus2_trafo`, `get_bus_gen`,
`get_bus_load`, `get_bus_sgen`, `get_bus_shunt`, `get_bus_status`,
`get_bus_storage`, `get_bus_svc`, `get_dcline_res1_full`,
`get_dcline_res2_full`, `get_gen_res`, `get_gen_status`,
`get_gen_target_p`, `get_gen_theta`, `get_line_res1`, `get_line_res2`,
`get_line_theta1`, `get_line_theta2`, `get_load_target_p`,
`get_load_theta`, `get_loads_res`, `get_loads_status`,
`get_sgen_target_p`, `get_sgens_res`, `get_sgens_res_full`,
`get_sgens_status`, `get_shunt_target_p`, `get_shunt_theta`,
`get_shunts_res`, `get_shunts_status`, `get_sn_mva`, `get_storage_theta`,
`get_storages_res`, `get_storages_status`, `get_trafo_res1`,
`get_trafo_res2`, `get_trafo_theta1`, `get_trafo_theta2`,
`is_grid_connected_after_contingency`, `reactivate_dcline`,
`reactivate_powerline`, `reactivate_sgen`, `reactivate_storage`,
`reactivate_svc`, `reactivate_trafo`, `remove_gen_slackbus`,
`set_max_nb_bus_per_sub`, `tell_recompute_sbus`, `tell_recompute_ybus`,
`tell_ybus_change_sparsity_pattern`, `update_sgens_p`.
