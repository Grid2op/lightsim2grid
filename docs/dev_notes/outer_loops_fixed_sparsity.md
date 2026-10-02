# OpenLoadFlow outer loops on a fixed Jacobian pattern

**Design note, not documentation.** It describes how lightsim2grid runs OpenLoadFlow's
AC outer loops itself, around a Newton-Raphson that keeps **one symbolic analysis and one
numeric factorization per solve**. It is the reference for the `dev_outerloops` work and is
updated as each phase lands. It is not part of the built documentation (`docs/*.rst`).

Companion notes in this folder:

- `outer_loop_checks_missing_data.md` -- what a post-solve "would this loop have fired?" check
  needs, written when the physical-violation checks were added. Part of it is superseded here
  (the tap / shunt / phase control data is now in scope).
- `nr_control_limits.md` -- the *in-Newton* alternative (limits as complementarity rows). It
  shares the main prerequisite with this note: reserve the union of the layouts up front and
  change values only.

## What is wanted

Two modes, one set of rules.

1. **Outer-loop mode.** A new, non-default algorithm family (`NROuter_KLU`,
   `NROuter_SparseLU`, ...) that behaves like OpenLoadFlow: the same loops, in the same default
   order, with the same triggers and the same actions, around a **single-slack** inner
   Newton-Raphson. The default `NR_*` / `NRSing_*` algorithms are untouched and stay
   bit-identical.
2. **Detection mode.** With any other algorithm, report -- with exactly the rules the loops use
   -- which outer loops would have run. This is what `LSGrid.get_physical_violations` and the
   batch `compute_physical_violations` return. Where today's checks differ from OpenLoadFlow
   they are aligned, which is a breaking change of their output.

Constraints:

- **One `analyze`, one `factorize` per solve**, whatever the loops do. The only accepted extra
  cost is the refactorization fallback: when a `refactorize` fails (a pin or a release can move
  a pivot that KLU's fixed sequence then finds at zero), the linear solver policy runs a fresh
  `factorize`. It is counted in `LinearSolverStats`, never hidden, and tests report it.
- **The loops declare their Jacobian entries up front**, when the solver cache is built and the
  grid data is known: every row, column and switchable bus any of their states may need. A
  generator that may switch between PV and PQ gets both forms reserved; a regulating
  transformer gets its ratio column and every row its control may use. This union pattern is
  what guarantees the single analysis; afterwards a loop only rewrites values.
- **The loops can be registered from Python, in any order** (`algo.clear_outer()`,
  `algo.add_outer(ReactiveLimits())`, ...). Without any call, the list is OpenLoadFlow's
  default list, in its order.
- **Detection is never re-coded.** The trigger of each loop is computed in one place, reused by
  both modes, and built on the existing check headers (see "Reuse" below).

Out of scope for this pass: secondary voltage control, area interchange control, automation
systems, three-winding transformers, loops written in Python (the interface keeps a pybind11
trampoline easy to add later), the batch classes *running* the loops (they keep detection;
`clone()` and per-solve state keep the door open), gpusim2grid.

## Reference

The reference is the OpenLoadFlow bundled in the pypowsybl used for the comparisons: OLF
**2.3.0** (tag `v2.3.0` of `powsybl-open-loadflow`). Read the source at that tag, not at `main`.
Its **default parameters are those of that pypowsybl build**, which differ from OLF's own
source defaults; print them with `pypowsybl.loadflow.Parameters()`. The values below are the
ones this work targets:

| group | parameters |
|---|---|
| loops on | `distributedSlack` (`PROPORTIONAL_TO_GENERATION_P_MAX`, `useActiveLimits`, `slackBusPMaxMismatch` = 1 MW, `slackDistributionFailureBehavior` = FAIL), `hvdcAcEmulation`, `svcVoltageMonitoring`, `useReactiveLimits` (`reactiveLimitsMaxPqPvSwitch` = 3, `voltageRemoteControlRobustMode`) |
| loops off by default, in scope | `transformerVoltageControlOn` (mode `AFTER_GENERATOR_VOLTAGE_CONTROL`, `transformerVoltageControlUseInitialTapPosition`, `generatorVoltageControlMinNominalVoltage` = 120 kV), `shuntCompensatorVoltageControlOn` (mode `WITH_GENERATOR_VOLTAGE_CONTROL`), `phaseShifterRegulationOn` (mode `CONTINUOUS_WITH_DISCRETISATION`) |
| engine | `maxOuterLoopIterations` = 30, `maxNewtonRaphsonIterations` = 30, `newtonRaphsonConvEpsPerEq` = 1e-4 (UNIFORM), `MAX_VOLTAGE_CHANGE` scaling (0.4 pu, 1 rad), `voltageInitMode` = DC_VALUES, realistic voltage [0.8, 1.2] pu checked on buses of nominal voltage >= 180 kV |
| network loading | `extrapolateReactiveLimits`, `forceTargetQInReactiveLimits`, `disableInconsistentVoltageControls`, `generatorsWithZeroMwTargetAreNotStarted`, plausible target voltage [0.8, 1.2] pu above 20 kV, `reactiveRangeCheckMode` = MAX, `plausibleActivePowerLimit` = 10000 MW |
| slack | `slackBusSelectionMode` = MOST_MESHED, `maxSlackBusCount` = 1, `referenceBusSelectionMode` = FIRST_SLACK |

Grid snapshots are loaded with `pypowsybl.network.load(path, {"iidm.die.with-extensions": "all"})`.

Classes to read (paths under `src/main/java/com/powsybl/openloadflow/`):
`lf/outerloop/config/DefaultAcOuterLoopConfig.java` (order), `ac/AcloadFlowEngine.java`
(driver), `ac/outerloop/*` (loops), `network/util/ActivePowerDistribution.java` and
`GenerationActivePowerDistributionStep.java` (slack), `ac/outerloop/tap/*` (transformer
loop helpers), `network/PiModelArray.java` (taps), `network/impl/LfShuntImpl.java` (sections),
`ac/equations/AcEquationSystemCreator.java` (which equations exist), `network/impl/
LfNetworkLoaderImpl.java`, `AbstractLfGenerator.java`, `LfGeneratorImpl.java` (loading rules).

### Accepted differences

- **The reference slack is a generator.** lightsim2grid's single slack stays on a generator;
  OLF picks a bus (MOST_MESHED). Comparisons either pin OLF on the bus of lightsim2grid's slack
  generator (`slackBusSelectionMode = NAME`, the tight comparison) or keep MOST_MESHED and
  rotate the angles so both agree on one common bus. In the second case the residual the
  distributed slack leaves (below `slackBusPMaxMismatch`) sits on a different bus, and the
  tolerance has to allow for it.
- **The inner stopping criterion.** OLF stops on the RMS of the mismatch over all equations,
  lightsim2grid on its infinity norm. Comparisons use tight solves on both sides, so both land
  on the same root; a loop threshold tied to OLF's epsilon (the reactive-limits loop uses
  `newtonRaphsonConvEpsPerEq` as its Q tolerance) is passed explicitly with the same value.

## The loops, in OpenLoadFlow's default order

`isNeeded` drops a loop that has nothing to do (no hvdc line in AC emulation, no stand-by SVC,
no regulating transformer, ...). Loops marked *off* are off with the default parameters but in
scope.

| # | loop | trigger (detection) | action | Jacobian |
|---|---|---|---|---|
| 1 | `DistributedSlack` | slack-bus active mismatch above `slackBusPMaxMismatch` | re-share the **cumulative** mismatch from the units' initial target P, factor `maxP / droop`, clamped to `[minTargetP, maxTargetP]`, never across 0 MW, saturated units leave; FAILED if a residual remains | Sbus values only |
| 2 | `FreezingHvdcACEmulation` | off (`startWithFrozenACEmulation` = false) | -- | -- |
| 3 | `HvdcAcEmulationLimits` | droop flow beyond the per-direction max (to saturation), back inside or reversed (release) | per-line regime: linear / saturated each way | values only; the droop entries are declared in every regime |
| 6 | `VoltageMonitoring` | the voltage a stand-by SVC monitors outside its thresholds | the SVC regulates at the high / low set-point, for good | values only: the SVC is a held controller of its own group from the start |
| 7 | `ReactiveLimits` | PV->PQ: generation Q (load Q included) beyond the bus' summed limits plus epsilon; PQ->PV: voltage back across the target on the right side; robust mode: remote controller outside the realistic band; limit moved since the freeze | freeze at the limit, release, V = 1 in robust mode, keep the strongest PV bus, at most `maxPqPvSwitch` PV->PQ per bus | switchable Vm / Q slots + pinning; remote / group controllers through the held-Q of `VoltageControl` |
| 8 | `PhaseControl` (*off*) | iteration 0: any controller; afterwards a current limiter above its value | continuous alpha in the first solve, rounded to the closest tap; a limiter moves one tap | an alpha column per regulating PST, row `alpha = alpha0` or `P = target` |
| 9 | `TransformerVoltageControl` (*off*) | controlled voltage outside half the deadband | INITIAL / CONTROL / COMPLETE step machine: continuous ratio, generators below 120 kV frozen PQ, bound clamping, initial-tap sharing, rounding | a ratio column per controller, row `rho = rho0`, the controlled bus' `V = target`, or an equal-sharing row |
| 11 | `ShuntVoltageControl` (*off*) | iteration 0: any controller | continuous B in the first solve, then dispatched and rounded to sections | a B column per controller shunt, row `B = B0`, `V = target` or equal sharing |

(4 area interchange, 5 secondary voltage control, 10 transformer reactive power control and 12
automation systems are out of scope.) The unrealistic-voltage check is part of the driver
below; it runs after the last loop able to fix it (`ReactiveLimits`, `TransformerVoltageControl`).

## Architecture

### The algorithm

`NRAlgo::compute_pf` is split into three protected pieces -- setup (rebuild decision,
`update_state`, `init_topology`, `build_J_sparsity`), the Newton loop, and the finalisation --
and `compute_pf` keeps calling them in sequence, so `NR_*` and `NRSing_*` do not move.
`NROuterAlgo<LinearSolver>` is built on the same pieces with the system

```
NRSystem<Base, VoltageControl, Hvdc, BranchControl, ShuntControl>
```

that is `SingleSlackNRSystem` plus two extensions for the continuous tap / phase and shunt
unknowns. The refactorization fallback is always on. The family is registered by name only, as
`NRRefactorRetry_*` is (no new `AlgorithmType` member: that enum is serialized).

The algorithm **owns copies** of the Sbus and Ybus values at stable addresses; the system
reads them through pointers, so a loop edits values in place and nothing goes through
`AlgoControl` (whose flags would force a rebuild, and which `LSGrid::ac_pf` resets after the
solve anyway). The grid is only ever seen `const`.

### The driver

OpenLoadFlow's `AcloadFlowEngine.run`, as of 2.3.0:

1. keep the loops whose `is_needed` holds, then `initialize()` each, in order -- before the
   first solve;
2. first solve; the loops run only if it converged;
3. passes: each loop, in order, repeats `check()`, re-solving after each UNSTABLE result, until
   it is stable, a solve diverges, or the global cap is reached. A pass stops at the loop that
   was the last one to be unstable (a full cycle has then been checked stable). Passes repeat
   while the pass added Newton iterations;
4. FAILED stops everything;
5. the unrealistic-voltage check runs after the last loop able to fix it, and on every later
   solve of that pass;
6. `cleanup()` in reverse order;
7. the outer status is STABLE only if the total iteration count stayed below the cap.

A loop sees its own iteration counter, which counts its UNSTABLE results only ("iteration 0"
means it has not changed anything yet).

Three details of that engine that are easy to get wrong (all checked against the 2.3.0
source and covered by `src/tests/test_outer_loop_driver.cpp`):

- "the last loop that was unstable" is **not** reset between passes: a pass breaks before
  running it, even as the first loop of the pass;
- the unrealistic-voltage check is never off. With robust mode and no loop able to fix an
  unrealistic state -- no loop at all included -- every solve is checked; otherwise the
  check waits for the end of the last such loop. It looks only at buses whose voltage
  magnitude is an unknown of the Newton (OLF's `BUS_V` variables: not a bus held at a
  set-point);
- OLF's Newton tests convergence after a step, so a re-solve is at least one iteration and
  "a pass added Newton iterations" means "a pass re-solved". lightsim2grid's Newton can stop
  at zero iterations, so the driver counts every re-solve as at least one.

On the comparison side: OLF's explicit `outerLoopNames` list cannot name
`AcHvdcAcEmulationLimits` (it is only built from the parameters, by `hvdc_ac_emulation`), so
that loop is isolated through the parameter-driven list, every other loop being off by its
flag (`utils/olf_outer_compare.py`).

### Resuming the Newton

Between two inner solves the state is kept: the voltages, the controllers' reactive
unknowns, the Jacobian pattern and the factorization. The next solve starts from the voltages
of the last one, except where a loop overrode them (the target of a released bus, V = 1 in
robust mode); from the current private Sbus (the slack loop moved targets, a frozen bus
injects its limit); and from the current private Ybus (a tap or a section changed). It
re-evaluates the mismatch there and refactorizes on its first iteration. It never goes through
`update_state`, which would restart from `V_init` and reset the controllers.

### `BaseOuterLoop`

The same non-virtual-interface contract as the element containers
(`element_container/GenericContainer.hpp`): public non-virtual entry points, each forwarding
to one protected `_xxx` hook.

| entry point | role | default |
|---|---|---|
| `name()` | the OLF name | -- |
| `declare(decl)` | claim rows, columns, switchable buses for every reachable state | nothing |
| `is_needed(ctx)` | OLF's `isNeeded` | true |
| `initialize(ctx)` | before the first solve | no-op |
| `detect(ctx, out)` | the trigger, as `LimitViolation`s; pure function of the grid, the solve and the loop state | -- |
| `check(ctx)` | `detect` then act; STABLE / UNSTABLE / FAILED | -- |
| `cleanup(ctx)` | after the last solve, reverse order | no-op |
| `can_fix_unrealistic_state()` | for the driver | false |
| `clone()`, `get_params()` / `set_params()` | copies, threads, persistence | -- |

`declare` is called from `init_topology`, before `build_J_sparsity`: the grid data and the
number of units that may switch, participate or regulate are known there. That is the only
place a loop shapes the pattern.

The context a loop works on holds the private Sbus / Ybus values, the pin and mask sets
(`set_pv_pinned_buses`, `set_masked_buses` -- pushed again after any sparsity rebuild, because
the mask positions are not recomputed by `build_J_sparsity`), the voltage overrides, the held
reactive values of `VoltageControl`, the per-line hvdc regime, the branch and shunt control
values, the solver bus maps, and the grid (const).

**Detection mode** calls the same `detect` on a loop with an empty state, right after a plain
solve. `LSGrid::get_physical_violations` iterates the default loop list; the batch classes do
the same per row.

### The DistributedSlack loop

`outer_loop/DistributedSlackLoop.{hpp,cpp}`. The units it shares on are the ones flagged
"can participate in the slack" (`set_gen_can_participate_slack` /
`set_storage_can_participate_slack`), with that weight as their key, and bounded by their
P limits: `init_from_pypowsybl(olf_rules=True)` fills both with OpenLoadFlow's participation
rule and `activePowerControl` target range, independently of the slack the Newton uses. The
arithmetic is `slack_redistribution::distribute`, called on the units' INITIAL targets with
the cumulative mismatch at every pass, as `ActivePowerDistribution.run` does; the residue is
OpenLoadFlow's `P_RESIDUE_EPS`. The slack bus' mismatch is the quantity
`LSGrid::compute_results` books on the slack generator (its residual, less what an in-Newton
slack carried), so the trigger reads the same after an `NRSing_*` or an `NR_*` solve.

The loop runs (`is_needed`), and `SLACK_MISMATCH` is reported in detection mode, only on a
grid with some unit flagged to share the slack. A grid where no unit is flagged has no
distributed slack set up at all: OpenLoadFlow with `distributedSlack` off, where the loop
does not exist, and its single slack is kept on purpose.

Checked against OpenLoadFlow run with that loop only, on real grid snapshots with a load
increase so that there is something to share: same voltages, same per-unit P, one symbolic
analysis per solve.

**The starting point matters.** OpenLoadFlow's `DC_VALUES` start takes the angles of a DC
load flow, and with the distributed slack on, that DC load flow already shares the
imbalance on the participating units. A single-slack DC instead (the slack generator taking
the whole imbalance) gives initial angles far enough from OpenLoadFlow's for the first
Newton to reach another root of a weak, radial part of a grid: both converge, to visibly
different voltages there, and the loop then keeps sharing from that other root (or a later
re-solve diverges from it). `LSGrid.set_dc_distribute_slack_on_can_participate(True)`
makes `dc_pf` share the imbalance on the units flagged "can participate" with their
weights. Only its angles are kept: `DC_VALUES` starts every magnitude at 1 pu
(`DcValueVoltageInitializer`), and the slack bus keeps its own angle, so a start must be
expressed relative to it -- a start rotated away from the slack is a different start, not
the same one.

With those angles the start is OpenLoadFlow's to a fraction of a degree: lightsim2grid's
DC (`1/x`, the transformer ratio, the same injections and the same sharing) is the DC that
`DcValueVoltageInitializer` solves. The comparison harness can also start from
OpenLoadFlow's own DC angles (`--start olf`), to tell a difference of start from a
difference of solve.

**The linear solver matters too.** A KLU refactorization reuses the pivots chosen at the
last factorization, here the Jacobian at the DC start. After a long first solve and a
large redistribution, one of those pivots can become tiny without being zero:
`klu_refactor` accepts it, the factors are too inaccurate, and the next Newton diverges
(SparseLU, which pivots again every time, does not). powsybl-math-native, hence
OpenLoadFlow, checks the reciprocal pivot growth (`klu_rgrowth`) after every
refactorization and factorizes again below powsybl's `DEFAULT_RGROWTH_THRESHOLD`. The same
check is part of lightsim2grid's refactor fallback (`LinearSolverPolicy::
set_refactor_fallback`, on for `NROuter_*`, `NRRefactorRetry_*` and the batch algorithms
that edit values), and the factorization it costs is counted as a fallback one.

### The VoltageMonitoring loop

With `svcVoltageMonitoring` (`OlfLoadingParameters.svc_voltage_monitoring`, on by default as
in OpenLoadFlow), an SVC whose standby automaton is in standby is a voltage monitor, loaded
idle: off, flagged standby with its thresholds and set-points (`LSGrid.set_svc_standby`).
Every algorithm sees it as off, Q = 0. With the loop in the grid's list and an algorithm
running outer loops, the voltage-control plan also enrols it as a held controller, alone in
a group regulating its bus (`hold_monitors`): the group's row reads "Q = 0", and switching
the SVC on rewrites it into the voltage row at the new set-point, by value. A monitor the
plan cannot enrol (its regulated bus already regulated by a group, or without a Vm unknown)
stays idle.

The automaton's `b0` is a susceptance the SVC carries (`LSGrid.set_svc_b0`), in standby or
not, stamped into Ybus at its bus and part of its reactive output, as OpenLoadFlow's
`LfStandbyAutomatonShunt`. OpenLoadFlow's loading rules around monitors (two on a bus both
regulate, one next to a regulating unit is switched off) are `_olf_rules.voltage_controllers`.

As OpenLoadFlow, the loop only runs when a monitor regulates its own bus, and the
thresholds are in pu of the SVC's own voltage level.

### Publishing the results

`LSGrid::compute_results` splits a generator's P with the grid's own slack weights and its Q
by reactive range. The outer algorithm exports per-element overrides through a new `BaseAlgo`
hook applied there: generator / storage P (the slack loop's targets), generator Q (frozen
limits), and the final PV / PQ state. Taps and sections are published as results
(`res_tap_position`, `res_section_count`), never written back into the inputs -- as pypowsybl's
`solved_tap_position` / `solved_section_count`.

The split of a bus' reactive power between its voltage-controlling units also has to follow
OpenLoadFlow (`AbstractLfBus.updateGeneratorsState` / `dispatchQ`, mode
`Q_EQUAL_PROPORTION`): by reactive keys when every unit has one; otherwise by reactive range,
but only when **every** unit has plausible limits (absolute limits below 1000 MVar, a range
between 1 and 10000 MVar); otherwise equally. With reactive limits on, a unit whose share
crosses a limit is clamped there and the rest is dispatched again over the others. Nothing
is dispatched at all while what is left is at most `Q_DISPATCH_EPSILON` (1e-5 pu, so 1e-3
MVar): a bus whose total is below it publishes 0 on its units.
`LSGrid::_split_q_residual_per_bus` always splits by range today, so a bus mixing a unit with
placeholder limits and a unit with real ones publishes a different per-unit Q (same bus total)
even when no loop runs.

### Registration and persistence

`algo.clear_outer()`, `algo.add_outer(loop)` and `algo.get_outer()` act on the `NROuter_*`
algorithm through a non-const accessor (the existing `get_algo()` binding is const by
convention). Each loop's parameters carry OLF's parameter names and the defaults of the table
above. A grid copy, a pickle or a batch thread only carries the algorithm name and its
`AlgoConfig`, so the loop list is kept as descriptors (type + parameters) next to it and
rebuilt with `clone()`.

## Reuse

Most of the detection, and some of the actions, already exist.

| piece | where | used by |
|---|---|---|
| bus reactive capability (PV->PQ) | `batch_algorithm/BusQCheck.hpp` | `ReactiveLimits.detect` |
| release of a frozen unit (PQ->PV), `can_be_pv` | `GenPvReleaseCheck.hpp`, `GeneratorContainer` | `ReactiveLimits.detect` |
| remote controller outside the realistic band | `RemoteVoltageControlCheck.hpp` | `ReactiveLimits.detect` (robust mode) |
| stand-by SVC thresholds | `SvcStandbyCheck.hpp` | `VoltageMonitoring.detect` |
| hvdc droop saturation / release | `HvdcPCheck.hpp`, `HvdcLineContainer::ac_emulation_frozen` | `HvdcAcEmulationLimits.detect` |
| unit past its p limits (in-Newton slack) | `GenPCheck.hpp` | the `NR_*` algorithms' detection only |
| OLF clamp-and-reshare, no sign change | `element_container/SlackRedistribution.hpp` | `DistributedSlack` action |
| "can take part in the slack" vs "is a slack" | `SlackParticipation::can_participate*` | `DistributedSlack` |
| PV <-> PQ at constant sparsity | `Base::set_switchable_vm_buses`, `NRSystem::set_pv_pinned_buses` | `ReactiveLimits` |
| a controller of a group held at a fixed Q | `LSGrid::set_hold_frozen_regulators`, `VoltageControl` held rows | `ReactiveLimits` (remote / groups), `VoltageMonitoring` (`NRSystem::release_held_svcs`) |
| refactorize -> factorize fallback, counters | `LinearSolverPolicy`, `LinearSolverStats` | the one-factorize accounting |

The check headers depend on no batch state; they take the grid, the solver bus map, the
mismatch and the voltages. Their rules move to OpenLoadFlow's where they differ: epsilon,
load Q in the bus balance, limits at the **current** target P, bus-level limits, no tolerance
on the PQ -> PV test.

New trigger kinds: `SLACK_MISMATCH`, `UNREALISTIC_VOLTAGE`, `REACTIVE_LIMIT_MOVED`,
`TRANSFORMER_VOLTAGE_DEADBAND`, `SHUNT_VOLTAGE_CONTROL`, `PHASE_CONTROL_P`,
`PHASE_LIMITER_CURRENT`.

## What the model lacks

- **Tap changers** -- done in phase 6 (`TapChangers.hpp`; the alpha -> r / x table of
  `set_shift_dependent_rx` still applies between two taps, relative to the current one).
  `TrafoContainer` holds a ratio and a shift, plus an alpha -> r / x
  correction table (`set_shift_dependent_rx`). It needs, per transformer, an optional ratio and
  phase tap changer: step table (rho, alpha, r %, x %, g %, b %), position and range, side,
  regulating flag, target voltage, deadband and regulated bus (ratio), regulation mode and value
  (phase). The pi-model follows OLF's `PiModelArray`: the current tap gives every parameter; a
  continuous ratio or shift overrides rho or alpha only, without interpolating r / x / g / b;
  rounding goes to the closest tap, a strictly better one only.
- **Shunt sections** -- done in phase 6. `ShuntContainer` holds a fixed admittance. It needs linear and
  non-linear sections, the current and maximum count, and the regulation (on, target, deadband,
  regulated bus).
- **Capability curves** -- done in phase 5 (`LSGrid::set_gen_capability_curves`). A generator's reactive limits are read once at its initial target P.
  The slack loop moves that target, and OLF re-evaluates the limits at the current one
  (extrapolating past the curve ends with `extrapolateReactiveLimits`), so the curve points
  have to be in `GeneratorContainer`.
- **`can_be_pv`** currently flags a PQ unit an outer loop froze. It becomes "this unit may be PV
  in some outer iteration": every PV generator, plus the frozen ones.
- **The converter** -- done in phase 6 for the tap changers and the sections. It read transformers at the neutral tap and dropped tap changers, shunt
  sections and regulation; it has to read them (`get_ratio_tap_changers` and steps,
  `get_phase_tap_changers` and steps, `get_shunt_compensators` and the section frames,
  `get_reactive_capability_curve_points`), with r / x / g / b at the current tap.
- **The loading rules.** OpenLoadFlow's network-loading rules are pure predicates in
  `from_pypowsybl/_olf_rules.py`, read by `init_from_pypowsybl(olf_rules=True)` and by
  `bake_outer_loops`. Done for the generators' voltage control (not started, reactive range,
  unresolved regulated bus, implausible target, inconsistent controls on one bus) and their
  target Q (clamped into the limits at target P, curves extrapolated), and for the
  participation in the distributed slack (`participation_weight`, generators and batteries
  alike: the extension's flag and droop, a zero or implausible target, the target range).
  Still to move there: the same voltage-control rules for batteries, VSC stations and SVCs --
  OpenLoadFlow's inconsistency check sees every voltage-controlling unit of a bus, not only
  generators.
  OpenLoadFlow's `POWER_EPSILON_SI` is compared with per-unit values in the not-started rule
  (so 0.01 MW) and with MW in the participation rule (1e-4 MW).

## The continuous controls in the Jacobian

With the default modes above, OLF solves the tap ratio, the phase shift and the shunt
susceptance **inside** the first Newton of their loop, then rounds them. In lightsim2grid terms:

- a new **`BranchControl`** extension adds one column per controlling transformer (rho or
  alpha). The column's entries are the derivatives of the four flow terms of that branch, so
  they sit in the P and Q rows of its two buses (plus the current term of a phase limiter). Its
  row is, by value: `x = x0` while the control is off; the controlled bus' `V = target` for the
  first enabled controller of a group; an equal-sharing row for the others. The union of a
  group's rows and columns is declared up front. Moving the unknown patches that branch's Ybus
  values at known positions, every iteration;
- a new **`ShuntControl`** extension does the same for B (`Q = -B V^2` in the Q row of the
  shunt's bus);
- a generator frozen to PQ by the transformer loop reuses the PV pinning slots.

Rounding to a tap or a section is a value patch on the private Ybus, which only refactorizes.

## Phasing

One loop at a time, in OpenLoadFlow's order. Each phase brings the loop, its detection wired
into `get_physical_violations`, unit tests on small grids, and a comparison against OLF run
with **that loop only** (`outerLoopNames = [name]` plus the flag that creates it) under tight
solves, on real grid snapshots.

| phase | content |
|---|---|
| 0 | this note; the comparison harness (`utils/olf_outer_compare.py`), baseline with no loop on either side |
| 1 | `NRAlgo` split, `NROuterAlgo`, `BaseOuterLoop`, declaration, resume, driver, statistics, result hook, Python registration and persistence. With no loop it must match `NRSing_*`, with one analysis -- **done**, but for the persistence of a custom loop list in a pickle or a binary file (with the format bump of phase 6) |
| 2 | `DistributedSlack`, participation rules unified -- **done** |
| 3 | `HvdcAcEmulationLimits` -- **done**; the snapshots rarely reach their limits, the harness's `--hvdc-limit-factor` lowers them on both engines |
| 4 | `VoltageMonitoring` -- **done**, with OpenLoadFlow's voltage-control loading rules for SVCs, VSC stations and batteries; the harness's `--svc-thresholds-pu` moves the automata's thresholds on both engines |
| 5 | `ReactiveLimits` -- **done**: local buses through switchable Vm / Q slots, group controllers held by value (`VoltageControl` holding, the sharing taken against any active controller), the robust mode; OpenLoadFlow's per-unit reactive split (`set_reactive_dispatch_olf`); capability curves on the grid (`set_gen_capability_curves`), read at the target P `DistributedSlack` gave; a monitor `VoltageMonitoring` switched on is checked like any SVC; a non-regulating unit's target Q clamped into its limits at the target P `DistributedSlack` moved (`set_gen_raw_target_q`, kept apart from the bus residual through `OuterState::Sbus_target`); VSC stations sharing a bus by their widest range, as generators. A plain solve's LOW_Q / HIGH_Q buses are the ones OpenLoadFlow switches PV -> PQ in its first round |
| 6 | tap and section data model, converter, binary format -- **done**: the pi model at the taps (OpenLoadFlow's, step corrections of r, x, g, b included), moving a tap or a section count keeps the Jacobian's pattern, closest-tap rounding; checked against OpenLoadFlow with taps and sections moved on both sides. A check that a control would act is `ViolationCategory::CONTROL` |
| 7 | `PhaseControl` -- **done** (CONTINUOUS_WITH_DISCRETISATION): a `BranchControl` NR extension, a shift column and a row (target alpha or target P, chosen by value) per active power phase shifter, the transformer's Ybus block patched by value in the algorithm's private Ybus; rounding to the closest tap, current limiters moved tap by tap, OpenLoadFlow's connectivity rule; the solved taps published in the results (`res_phase_tap_position`), the inputs untouched. Matches OpenLoadFlow tap for tap, except where its one-tap move undoes itself on the rounding of its Newton (see `PhaseControlLoop.hpp`). **TODO (follow-up)**: decide whether to mirror that rounding or keep it documented |
| 8 | `TransformerVoltageControl` -- **done** (AFTER_GENERATOR_VOLTAGE_CONTROL): `BranchControl` (the former `PhaseShift`) gives each transformer of a visible group a ratio column and the group n rows, each holding by value OpenLoadFlow's BRANCH_TARGET_RHO1, BUS_TARGET_V or DISTR_RHO; groups hidden by a generator's control only round their taps. INITIAL / CONTROL / COMPLETE as OpenLoadFlow, the shared ratio mapping (`useInitialTapPosition`), the generators of buses up to 120 kV frozen at their bus' injection -- OpenLoadFlow passes the injection, load included, as the generation; reproduced, and checked to decide the solved tap -- and released after; `fixTransformerVoltageControls`. Other loops leave a frozen bus alone (`OuterState::suspended_buses`). Matches OpenLoadFlow tap for tap on IEEE 14 and on real grid snapshots as loaded |
| 9 | `ShuntVoltageControl` -- **done** (WITH_GENERATOR_VOLTAGE_CONTROL): a `ShuntControl` NR extension, the regulating shunts of a bus ONE controller with a susceptance column, groups by regulated bus with the same three row forms as the ratios, the Ybus diagonal patched by value; on from the first solve, then `dispatchB` (largest shunt first, each rounded to its closest section) and one more solve. A control hidden by a higher-priority one (generator, SVC, VSC, transformer) never acts -- OpenLoadFlow's `getControllerElements` keeps the visible ones, which also applies to the transformer groups. Matches OpenLoadFlow on IEEE 14; the real grid snapshots have no regulating shunt |
| 10 | all loops together -- **done**: the seven loops in OpenLoadFlow's order match it on every real grid snapshot it converges on, as loaded and with random lines opened (N-1), and with stressed settings that make every loop act; the same taps, one symbolic analysis per solve. Detection cross-check (`utils/olf_outer_compare.py --detection`): a plain solve never misses a loop OpenLoadFlow runs in its first round; a loop that acts only once another changed the state (DistributedSlack after ReactiveLimits moved the losses, TransformerVoltageControl after the voltages moved) is not -- and cannot be -- seen by a plain solve. OpenLoadFlow's report lists only the loops that log something, so the harness reads the silent ones off the positions they leave |

Every outer test asserts a single analysis and reports the fallback factorizations.
