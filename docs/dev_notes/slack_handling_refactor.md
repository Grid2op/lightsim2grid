# Refactoring the slack handling

**Design note, not documentation.** Only what "Already fixed" says is implemented. It
records the root cause behind several slack problems found around the `NROuter_*`
algorithms, the target model, a plan, and the decisions to take before starting. It is not
part of the built documentation (`docs/*.rst`); the matching entry is in `CHANGELOG.rst`'s
`[TODO]` ("The Newton-Raphson conflates two different things under one flag").

## The root cause: two roles under one set

The grid holds one **slack set**: the generators and storage units tagged slack
(`SlackParticipation::slackbus_`, with `weight_`). Being in it means two different things
at once:

1. **the reference role** -- a property of a **bus**: its P is free (it has no P equation),
   its angle is the datum, its magnitude is fixed when a local regulator pins it. A Newton
   needs exactly one per island;
2. **participation in the distributed slack** -- a property of a **unit**: it takes a share
   of the active-power imbalance, by weight and within its limits. Whether that share is
   computed inside the Newton (`slack_absorbed`) or around it (the `DistributedSlack` outer
   loop) is the algorithm's business. The bake (`_olf_bake.py`) is the outer loop's result
   computed *a priori*: the two must not be treated differently.

What the set means then depends on who reads it:

| algorithm | `slack_ids(0)` | the other slack buses | weights |
|---|---|---|---|
| `NR_*` (`MultiSlackNRSystem`) | P equation + `slack_absorbed` share | theta + P equation + share | shared by `slack_absorbed` |
| `NRSing_*` (`SingleSlackNRSystem`) | P-free, angle fixed | **P-free, angle fixed** too | ignored |
| `NROuter_*` | P-free (the single-slack Newton) | PV, with their P equation (see "Already fixed") | the loop shares on the "can participate" units |
| FDPF | P-free | PV (`retrieve_pv_with_slack`) | shared when there is more than one slack |
| Gauss-Seidel | P-free | PV | ignored |
| DC | removed (reference) | PV | the imbalance is spread by weight |

Everyone but `NRSing_*` already reads `slack_ids(0)` as "the reference" and the rest as
"participants". `NRSing_*` gives every slack bus the reference role -- what several
pandapower `ext_grid`s are (`test_multi_slack.py` checks it against `pp.runpp` with two of
them). `docs/use_solver.rst` describes that; `docs/solvers.rst` says the opposite (the
others "are treated as PV buses").

The participants the set could not hold were given a second flag, **"can participate in
the slack"** (`can_participate_weight_`), and that flag has itself ended up with two
meanings:

- the bake's (`_aux_add_slack`): a unit it left out of the slack *only because it sat at an
  active limit* -- it goes back into the slack once moved off that limit;
- the outer loop's (`_aux_olf_slack_participation`, `olf_rules`): every unit the
  `DistributedSlack` loop shares on, kept out of the slack set *by construction* because
  `NROuter_*`'s Newton has a single slack bus -- at any set-point.

## What it causes

1. **`NROuter_*` lost P equations** (fixed, see below): a unit joining the slack set made
   its bus a second P-free bus of the single-slack Newton; a target the loop moved onto it
   changed no equation, the next solve converged in zero iterations, and the loop crept
   towards a wrong dispatch.
2. **The rejoin rule fires too widely.** `change_p_gen` / `change_p_storage`, a
   reconnection and `update_topo` call `rejoin_slackbus_if_able`, which adds to the slack
   set any flagged unit that is out of it, connected and strictly inside its
   `[min_p, max_p]`. Right for the bake's meaning, wrong for the outer loop's: on an
   `olf_rules` grid, every unit whose set-point `change_p_gen` touches (the grid2op backend
   does, every step) joins the slack set, and with `NR_*` the dispatch then depends on
   which set-points were changed before. The other half of the intended rule -- a unit in
   the slack that a set-point change puts on a limit leaves it -- is missing: only the
   redistribution pre-pass (`redistribute_active_power` -> `leave_slackbus`) does that.
3. **The reference moves.** `LSGrid::_slack_bus_id_me` lists the slack buses in generator
   order, then storage. A unit with a lower id joining the slack set becomes `slack_ids(0)`:
   the angle datum moves, and for `NROuter_*` the leftover the loop leaves below its
   threshold lands on that unit.
4. **No reference without a unit.** OpenLoadFlow's slack bus is a bus
   (`slackBusSelectionMode`), which may carry no generator. lightsim2grid's slack is a set
   of units, so the `olf_rules` loader picks a generator as a proxy for OpenLoadFlow's
   choice (its own comments call it a poor one) and the OpenLoadFlow comparisons pin
   `slackBusSelectionMode=NAME` on that generator's bus.
5. **Participation needs a voltage source.** Only generators and storage units can be in
   the set, and a slack unit is PV (`_fillpv` skips slack buses, `is_pseudo_off` never
   holds for one): a load cannot take a share (the `[TODO]` entry).
6. **The reference is not saved.** `_forced_ref_slack_bus_id` (`set_reference_slack_bus`)
   and `dc_slack_on_can_participate_` are copied by the copy constructor but are neither in
   `get_state` / `set_state` nor in the binary format: pickling or `save_binary` loses them.

## Already fixed

`OuterLoopAlgo::compute_pf` hands its Newton `slack_ids(0)` alone and appends the other
slack buses to `pv`: they keep their P equation (and, when no local regulator pins them,
the Vm unknown and Q equation `Base` gives a free-Vm slack bus). Test:
`test_outer_loops.TestDistributedSlack.test_changed_target_shared_as_a_loaded_one`. This is
the target model below, applied locally to one algorithm.

## Target model

- **The reference** is a bus, one per island, held by the grid on its own: set explicitly
  (any bus, with or without a unit), otherwise derived once (today's rule: the first slack
  unit's bus) and kept while the participants change. It is the only P-free bus, for every
  algorithm.
- **The participants** are units with a weight. Each is *active* or *parked at a limit*
  (min or max). Set-point changes follow one rule, whoever makes them (`change_p_*`, a
  reconnection, `update_topo`, the pre-pass, the bake):
  - a parked unit moved strictly inside its limits becomes active again;
  - an active unit put on a limit is parked;
  - both raise `tell_slack_participate_changed` / `tell_slack_weight_changed`.

  The bake marks the units it capped as parked (it knows the side: the overshoot it already
  records). The "can participate" flag disappears: the outer loop shares on the
  participants, active or parked, within their limits; an in-Newton share
  (`slack_absorbed`) goes to the active ones only, since it has no bounds.
- **What the algorithms are handed does not change shape**: `slack_ids` = the reference
  first, then the buses of the active participants, with `slack_weights` (0 at a reference
  that carries none). It is already what every family but `NRSing_*` reads (table above),
  so `BaseAlgo::compute_pf` keeps its signature (plugin ABI); `docs/solver_plugin.rst`
  states the meaning.

## Plan

Test first at every step (a test failing on the current code, then the change).

1. **Core state.** `SlackParticipation`: participants + per-unit parked side, replacing
   `slackbus_` / `can_participate_weight_`; the rejoin / leave rule in the containers'
   `_on_change_p` / `_on_reactivate` / `update_topo`. `LSGrid`: a reference bus field (in
   `get_state`, `set_state` and the binary format, with `dc_slack_on_can_participate_`);
   `_slack_bus_id_me` builds the reference first, then the active participants' buses.
2. **Algorithms.** `NRSing_*` per decision 1 below. `OuterLoopAlgo::compute_pf`'s local fix
   stays correct and becomes redundant.
3. **Results.** A reference bus without a unit: its leftover P is published nowhere today
   (`SlackParticipation::split` only splits onto units). Expose it (decision 2); the
   `SLACK_MISMATCH` detection and `GenPCheck` read it instead of the slack generator's P.
4. **Topology.** `consider_only_main_component` (decision 3), the batch's
   `_select_ref_slack_and_masks` (candidates are slack buses with a positive weight today),
   `set_topo_actions`' refusals for slack participants.
5. **Outer loops.** `TransformerVoltageControlLoop` and `ShuntVoltageControlLoop` count
   every slack bus as generator-controlled: only the buses a unit pins are.
6. **Loaders.** `olf_rules`: the participants go into the participant set, the reference is
   OpenLoadFlow's slack bus (replicate `MostMeshedSlackBusSelector` from OpenLoadFlow
   v2.3.0, or take the bus it reports). The bake: capped units parked. pandapower /
   PowerModels: per decision 1.
7. **Bindings, docs, changelog.** `GenInfo` / `StorageInfo` slack fields;
   `set_gen_can_participate_slack` and friends deprecated or removed (decision 4);
   `docs/solvers.rst` and `docs/use_solver.rst` made consistent (DC and FDPF do distribute
   today, both pages say they do not); `[BREAKING]` entries.

## A distributed-slack Newton inside the outer loops: `NROuterDist_*`

**Implemented**, as a family of its own; `NROuter_*` is untouched and stays the exact
OpenLoadFlow reproduction. It needed none of the refactor above: everything lives in the
outer-loop driver and the Newton's `MultiSlack`.

- **The algorithms.** `NROuterDist_*` = `OuterLoopAlgo<NROuterInner<LS, OuterDistNRSystem>>`,
  `OuterDistNRSystem` = `NRSystem<Base, MultiSlack, VoltageControl, BranchControl,
  ShuntControl, Hvdc>`. The reference (`slack_ids(0)`) stays the only slack bus, and the
  other slack buses are PV with their P equation, as in `NROuter_*`.
- **Participants stay pv / pq buses.** `MultiSlack::set_slack_participant_buses` gives each
  one an entry in the `slack_absorbed` column, at the P row Base already gave it: no new
  row or column, so voltage pinning, PV / PQ switching and the TVC / shunt loops see them as
  before. Their weights are values (`set_slack_weights`), so parking a unit costs no
  `analyze`. With no participant listed, the layout is `NR_*`'s, bit for bit.
- **Who sets what.** `DistributedSlack` reserves the buses (`SlackShareControl`, in
  `declare`): every unit flagged "can participate", on its own bus in the grid and on the
  one the declaring row places it on. It then hands the units' weights over each solve.
  The driver normalises them, and gives the reference weight 1 when there is none (no
  loop, every unit parked): a single slack in the Newton's form, so J is never singular.
- **Publication.** After every solve the driver moves each unit's share into its target
  (`OuterLoopAlgo::_fold_slack_share`), and `get_slack_absorbed()` is 0. The share enters
  the mismatch only through `+ slack_absorbed * weight`, so the mismatch does not change.
  Whatever reads the solve then reads what it reads after `NROuter_*`: the loops, the
  results (`get_outer_target_p`), resume points, `GenPCheck`.
- **The loop.** It shares the same cumulative mismatch as around a single slack, from the
  initial targets, with `slack_redistribution::distribute`. A unit that this saturates
  (above `max_p`, below `min_p`, across 0 MW) is set at its bound with weight 0. One that
  comes back inside its bounds gets its weight back. UNSTABLE when the targets moved, so
  the loop is stable after the first solve when no unit hits a limit. A resumed solve (a
  batch row from its "n" solve, a time series' next step) starts with the units the solve
  it resumes from parked still parked, at their bound: a snapshot's units already at
  `max_p` are then not pushed past it again by every row's first solve. A batch row also
  parks the units the share of its lost power (`_share_lost_p`) put at a bound. That one is
  a bet on the sign of the change in the losses: it saves the re-solve where the
  contingency raises them, and costs one where it lowers them. On real grid snapshots, it
  wins clearly overall. When every unit is
  saturated it is FAILED, or with LEAVE_ON_SLACK_BUS every unit is parked and the
  reference takes the rest.
- **Batch.** The reservation holds for every row. A row's off or stranded unit is not
  among the loop's units, so it gets no weight. Rows are tested against single solves in
  `test_outer_dist_batch.py`.

The dispatch differs from OpenLoadFlow's by at most `slackBusPMaxMismatch`, which OpenLoadFlow
leaves on the reference bus and this variant shares out. The in-Newton alternative to the
loop (projection rows solved by a semismooth Newton) is in `nr_control_limits.md`
("Distributed slack with min_p / max_p").

Left open:

- A batch row that places a participant on a bus no row declared: `SlackShareControl::
  set_units` leaves the unit out of the Newton's share.
- Detection outside a solve (`LSGrid::get_physical_violations`, a context without
  controls) keeps the single-slack rule (SLACK_MISMATCH on the reference).
- A parked unit coming back is handled by construction (each pass shares from the initial
  targets), but no test builds that case.

## Decisions to take first

1. **Several P-free buses.** `NRSing_*` with several slack units is how several pandapower
   `ext_grid`s are solved today. Keep it as an explicit, separate notion ("fixed buses":
   angle and magnitude fixed, P and Q free), or drop it and make `NRSing_*` read the slack
   set like every other algorithm?
2. **The leftover of a reference without a unit**: a getter only (OpenLoadFlow's
   `slackBusActivePowerMismatch`), or also a published pseudo-injection on the bus?
3. **A topology change that strands the reference**: pick a new one in the main component
   (OpenLoadFlow does), or fail / skip the row as a batch does today?
4. **API**: keep `add_gen_slackbus` / `set_gen_can_participate_slack` as deprecated
   wrappers mapping onto the new model, or remove them in the same `[BREAKING]` release?
5. **Loads as participants** (`[TODO]`, problem 5): in this refactor, or later on top of it?
6. ~~The distributed-slack `NROuter_*` variant~~: done before this refactor, as
   `NROuterDist_*` (see above).

Rejected shortcut for problem 4: a dummy 0 MW PQ generator on the bus, tagged slack. It
needs no C++ change -- the leftover lands on it -- but adds an element users and grid2op
see (`n_gen` changes).

## Blast radius

Where the slack set is read or written today, by area (entry points; follow them from
there):

- **`LSGrid`**: `_slack_bus_id_me`, `_build_into_cache` (slack ids, forced reference, raw
  weights, DC on can-participate), `init_slack_bus`, `compute_results` /
  `generator_results` / `_fill_bus_mismatch_*`, `get_gen_slack_shares`,
  `assign_slack_to_most_connected`, `redistribute_active_power`,
  `cap_slack_at_active_limits`, `consider_only_main_component`, `get_physical_violations`,
  `check_grid`, `check_solution_q_values`, `change_p_*` / `reactivate_*` / `update_topo`
  (rejoin), the slack getters.
- **Containers**: `SlackParticipation` (all of it), `GeneratorContainer` /
  `StorageContainer` (slack API, `StateRes` slots, `Info`, `is_pseudo_off`, `_on_*` hooks),
  `VoltageSourceContainer::_fillpv`, `SlackRedistribution` (`Participant.in_slack`),
  `GenPCheck`.
- **`VoltageControlPlan`**: `build_pv_pq`, `build_free_vm_slack`, the controller checks on
  slack buses.
- **Algorithms**: `NRSystem` `Base` / `MultiSlack`, `NRAlgo`, `BaseAlgo` (input checks,
  `retrieve_pv_with_slack`, `distributes_slack`), FDPF, Gauss-Seidel, DC (and the
  single-slack PTDF), `OuterLoopAlgo`, `DistributedSlackLoop`, the TVC and shunt loops.
  `SlackPolicies.hpp/.tpp` is dead code (included nowhere).
- **Batch**: `BaseBatchSweep` (`_select_ref_slack_and_masks`, `_row_slack_weights`,
  `_masked_slack_weights`, the per-row redistribution, `_record_row_*`, the topology-action
  refusals), `BaseBatchSolverSynch`, `ContinuationSweep`.
- **Binary format**: the generator / storage slack slots (`BinaryArchive.hpp`).
- **Bindings**: `binding_lsgrid.cpp` (slack API), `binding_containers.cpp` (`Info` fields),
  `binding_batch.cpp` (`redistribute_slack`, `pick_reference_slack`), `help_fun_msg.cpp`.
- **Python**: `from_pypowsybl/_aux_add_slack.py`, `_aux_handle_slack.py`, `initLSGrid.py`,
  `_olf_bake.py`, `_olf_params.py`, `_olf_compare.py`, `_result_network.py`; the pandapower
  and PowerModels slack converters; `lightSimBackend.py` (`_assign_right_solver`,
  `_handle_dist_slack`); `compare_lsgrid.py`; `pandapower_compat`.
- **Tests**: C++ -- `test_powerflow_algorithm`, `test_slack_redistribution`,
  `test_batch_physical_violations`, `test_lsgrid`, `test_distributed_slack_loop`,
  `test_cache_reuse`, `test_dc_slack_on_can_participate`, the voltage-control and FDPF
  ones. Python -- `test_main_component_slack`, `test_can_participate_slack`,
  `test_batch_redistribute_slack`, `test_multiple_slack`, `test_multi_slack`,
  `test_storage_distributed_slack`, `test_cap_slack_at_active_limits`, `test_dc_dist_slack`,
  `test_dist_slack_jacobian`, `test_ScenarioSweep_topology`, `test_ContingencyAnalysis_split`,
  `test_outer_loops`, `test_solver_control`, `test_physical_violations`,
  `test_voltage_control_pypowsybl`, `test_continuation_powerflow`.
- **Docs**: `use_solver.rst`, `solvers.rst`, `algorithm_names.rst`, `solver_plugin.rst`,
  `security_analysis.rst`, `scenario_sweep.rst`, `comparison_with_pypowsybl.rst`,
  `lightsimbackend.rst`, `network.rst`.
