# OpenLoadFlow-style outer loops on a fixed Jacobian pattern

**Assessment note, not documentation. Nothing here is implemented.** It records what it
would take to run OpenLoadFlow's outer loops inside lightsim2grid — for a single solve and
for the batch algorithms — and inside gpusim2grid, under one constraint: **one `analyze` and
one `factorize` per solve (per batch), whatever the outer loops do.** It is not part of the
built documentation (`docs/*.rst`).

Companion to two notes in this folder:

- `outer_loop_checks_missing_data.md` — the *detection* question: what a post-solve "would
  this outer loop have fired?" check needs, and which control data the converters drop;
- `nr_control_limits.md` — the *in-Newton* alternative: enforcing limits as complementarity
  (projection) rows solved by a semismooth Newton.

This note is the third option, the one OpenLoadFlow itself uses: keep the Newton-Raphson as
it is and iterate **around** it. It shares its main prerequisite with `nr_control_limits.md`
(reserve the union of the layouts up front, change values only), so the work is not wasted
if the in-Newton route is taken later.

Reference for the outer loops: PowSyBl OpenLoadFlow,
`lf/outerloop/config/DefaultAcOuterLoopConfig.java` (registration order = nesting order,
innermost first), `ac/AcloadFlowEngine.java` (`runOuterLoop` and the pass loop around it),
and the loop classes in `ac/outerloop/`.

## The short version

It is feasible, and most of the hard part already exists. Three value-level edits that keep
the Jacobian pattern fixed are already in the NR system, and gpusim2grid mirrors all three:

- **PV ↔ PQ at constant sparsity** — `Base::set_switchable_vm_buses` reserves a Vm unknown and
  a Q equation, `NRSystem::set_pv_pinned_buses` pins the Q row to the identity while the bus
  stays PV;
- **bus masking** — `NRSystem::set_masked_buses`, both P and Q rows to identity;
- **row repurposing** — the stranded lone controller of `VoltageControl`, whose voltage row
  becomes `Q_c = 0` by value alone.

What is missing is a driver, a *resume* path for the Newton, a handful of setters, and a
correction to how results are published.

It is also worth knowing that OpenLoadFlow does **not** keep its structure: any equation or
variable activation change marks its `JacobianMatrix` `STRUCTURE_INVALID`, which drops the LU
and runs a full `decomposeLU` again. Every PV → PQ switch pays a symbolic analysis there. A
fixed-pattern design is an improvement on that point, not only parity.

## Loop by loop

| OLF loop | what it changes | fixed-pattern mechanism | state |
|---|---|---|---|
| `ReactiveLimits`, local PV | PV → PQ with Q at its limit, PQ → PV back when the voltage recrosses the set-point | switchable Vm/Q slots + pinning, Q-limit injection | **ready**; needs a Vm setter for the PQ → PV direction |
| `ReactiveLimits`, remote / SVC groups | a controller at its limit leaves its group | lone controller: the existing `(v_row, q_col)` slot, generalising `Q_c = 0` to `Q_c = Q_lim`; group of N: reserve the group's full border block (1 voltage row + N-1 sharing rows × N Q columns + the regulated Vm column) | lone: easy; group: the hard case |
| `DistributedSlack` with `min_p` / `max_p` | a saturated machine leaves the distribution | its weight to zero, its clamped P into Sbus; every participant is already a slack bus of `MultiSlack`, so this is values only | **ready**; `GenPCheck` is the oracle |
| `HvdcAcEmulationLimits` / freezing | a droop line saturates or is frozen | `Hvdc` declares its entries whatever `status_droop` is | **ready**; needs a per-line status setter on the extension |
| `MonitoringVoltage` | a PQ bus becomes PV | pin the Q row the PQ bus already owns | easy |
| tap / phase / shunt voltage control | `Ybus` values (ρ, α, b) | pattern unchanged; value patch as `YbusPolicy::Contingency` does | **blocked**: the control data is not in the model |
| secondary voltage control, area interchange | set-points, P | values only | **blocked**: no pilot points, no `Area` |
| `AutomationSystem` | switches branches | closing a branch absent from the base pattern changes it | hard, out of scope |

The blocked rows are blocked on the converters, not on the solver
(`outer_loop_checks_missing_data.md`); `bake_outer_loops` remains the way to get their effect
until then.

## What lightsim2grid needs

### A new, non-default algorithm

`NRAlgo` is `final`. Split its `compute_pf` into *setup* (phases 1 and 2), *Newton loop* and
*finalise*, and build `NROuterAlgo<LinearSolver, NRSystem>` on the pieces. Register it as
`NROuter_KLU`, `NROuter_SparseLU`, `NROuter_NICSLU`, `NROuter_CKTSO`. The default `NR_*`
algorithms are untouched and stay bit-identical.

The Python API, through `BaseAlgo` virtuals that throw by default (the `supports_cpf`
pattern) forwarded by `AlgorithmSelector`:

```python
algo = grid.get_algo()
algo.clear_outer()
algo.add_outer(DistributedSlackLimits(...))
algo.add_outer(ReactiveLimits(max_pq_pv_switch=3))
```

### An `OuterLoop` interface in the core

Plain C++ in `src/core` (the core stays Python-free), with a pybind11 trampoline in the
bindings so that a loop can be written in Python:

- `name()`;
- `reserve(ReservationSink &)` — called **before** `build_J_sparsity`: switchable buses,
  extra structural entries, custom rows / columns. This is where "one analyze" is decided.
  An `add_outer` after the first solve forces exactly one rebuild, and says so;
- `initialize(ctx)`, `check(ctx) -> STABLE | UNSTABLE | FAILED`, `cleanup(ctx)`;
- `clone()` — one instance per batch thread.

### The driver

Copy OpenLoadFlow's semantics exactly, so that results stay comparable with pypowsybl
through the existing `_olf_compare` tooling: each loop, in registration order, runs until it
is stable, re-solving after each change; the whole list is passed again until a full pass
needed no Newton iteration; a pass stops at the last loop that was unstable; the total number
of outer iterations is capped. A "flat" variant (check all, apply all, re-solve once) is
cheaper but is a different algorithm.

### Resuming the Newton

The core missing piece. `NRSystem::update_state` restarts from `V_init`, re-seeds
`slack_absorbed` and resets the controllers' reactive output to zero. An outer iteration must
instead keep the whole state, re-evaluate the residual at the current point after the
loop's edits, and **refactorize** on its first iteration.

What makes that cheap:

- the system reads Sbus through a raw pointer (`Sbus_data_ptr_`), so an algorithm-owned copy
  can be edited in place between passes;
- nothing goes through `AlgoControl`: its flags feed `need_rebuild` in `NRAlgo.tpp`, which
  would trigger a full rebuild.

What needs a setter: the slack weights (`MultiSlack` keeps a copy), the Vm of a PV bus (the
PQ → PV direction, and secondary voltage control later), the per-line `status_droop` of the
`Hvdc` extension, a per-controller "Q pinned at this value" on `VoltageControl`.

### The edit surface, for Python loops too

The context a loop sees: an Sbus copy, the slack weights, the pinned and masked sets, the Hvdc
statuses, the Vm set-points, the controller pins, and the row / column maps read-only.

Worth adding: generic **J / F value overrides at reserved positions**, the counterpart of
gpusim2grid's `jov` stream. With it, a first version of a loop can pose a custom equation from
Python without any C++ change. A loop never writes `J` directly: `fill_J` zeroes it first.

The GIL: the algorithm's `solve` binding releases it, so the driver re-acquires it around a
Python callback. Python loops are refused in the batch classes.

### One analyze, one factorize

- **One analyze** is guaranteed if `reserve()` asks for the union pattern. For reactive
  limits that is a Vm column and a Q row, with their full dS pattern, for every PV bus with a
  finite reactive range — paid on every solve, pinned or not, because the symbolic analysis
  counts structure, not values. To be measured. The alternative is to reserve a list and
  count any later re-analysis in `LinearSolverStats`.
- **One factorize** has a caveat. A refactorization keeps the pivot sequence, and pinning or
  releasing a row can shrink a pivot; `RefactorRetryLinearSolver` then falls back to a fresh
  numeric factorization. Not an analysis, but a second factorization. It should be rare —
  both states have a strong diagonal, the identity when pinned and `dQ/dVm` when live — but
  it cannot be promised. Either expose the count and assert on it in the tests, or offer a
  strict mode that fails rather than falls back.

### Publishing the results

Easy to miss. `LSGrid::compute_results` splits a generator's active power with the grid's
own slack weights (`ac_cache_.slack_weights`) and the raw per-generator weights, not the
algorithm's. A machine the slack loop saturated would be published with the wrong P, so the
algorithm has to export per-generator overrides through a new hook.

The reactive side works if the injections the loops add are kept **out** of the algorithm's
per-bus mismatch: the mismatch at a bus switched to PQ is then its limit, and
`_split_q_residual_per_bus` publishes it as it does today.

### The batch algorithms

They get it almost for free, since every row calls `compute_pf`. What remains:

- `BaseBatchSweep` already drives `set_switchable_vm_buses` / `set_pv_pinned_buses` for
  generator contingencies (`_push_switchable_to_algo`), and the last caller wins: the
  switchable sets must be **merged**, and a bus is pinned only if both sides keep it PV;
- masked buses are skipped by the checks;
- the outer state is reset per row for the sweeps; whether a `TimeSeries` carries it from one
  row to the next (and inherits the hysteresis) is a decision to take;
- one loop instance per thread (`clone()`);
- `BusQCheck` / `GenPCheck` would then report only what enforcement left over.

## What gpusim2grid needs

The value-level primitives already exist there: the `add_switchable_vm_buses` ledger
extension, the identity-row stream for pinned buses, per-slot slack weights
(`slack_w_stride`), and the stranded-controller J overrides. The gaps:

1. **Masks are fixed per run.** They are host-built sparse streams uploaded once
   (`mask_streams.cuh`). Outer loops need device-resident, per-slot state: pinned flags per
   (slot, reserved bus), a per-slot Sbus (the contingency analysis shares one today), per-slot
   weights — and one kernel writing the identity rows at J positions precomputed once, since
   every slot shares the same pattern.
2. **No convergence test.** `run_nr_loop` runs a fixed number of iterations. The outer driver,
   per chunk: Newton loop, residual (the post-loop kernel exists), per-slot check kernels,
   one device-side "any unstable" reduction copied back per outer iteration. The Q / P checks
   are trivially parallel; per-slot switch counters live on the device.
3. **Resume.** Between passes, skip `init_slack_absorbed_kernel` and the controller-Q reset.
4. **Strategies.** The ones that reuse old factors (`iter0_only`, `direct_base_case_factors`)
   conflict with pin changes: refuse them or force a refactorization. The cuDSS uniform batch
   has the same pivot risk as KLU, with no cheap per-slot fallback; the existing pinned-bus
   paths are the place to measure it.
5. **Coupling.** The `AcPfGPU` and `InjectionSweepGPU` ledgers are never extended today. Any
   row or column `NROuter` adds must be exported explicitly by lightsim2grid: the bridge
   assumes the `VoltageControl` rows come last.

A first Python version is possible there too, at batch granularity: run, check in cupy over
DLPack, update the per-slot state, run again. It needs a warm start from the current V and the
cuDSS analysis kept across `run()` calls — today it is redone at each driver setup.

## A possible phasing

Rough orders of magnitude, for one developer who knows the code, working with Claude Code.
They are estimates, not measurements.

| phase | content | estimate |
|---|---|---|
| A. single solve, Python loops | split `NRAlgo`, `NROuterAlgo`, `OuterLoop` + trampoline, resume, setters, J / F overrides, results hook; `ReactiveLimits` (local), `DistributedSlack` with limits, `HvdcAcEmulationLimits` in Python; tests against pypowsybl and on the analyze / factorize counts | 1–2 weeks |
| B. C++ loops + batch | ports, merged pinning, `clone()`, per-row reset, `TimeSeries` policy, tests on the four batch classes | 1–2 weeks |
| B+ (optional) | reactive limits on multi-controller `VoltageControl` groups | +1 week |
| C. gpusim2grid | device state, per-slot pin kernel, outer driver, check kernels, ledger extension everywhere, export | 2–3 weeks |

What moves the estimates: parity with OpenLoadFlow on `ReactiveLimits`' edge cases (the
PQ → PV direction, the switch cap, robust mode, unrealistic voltages); refactorization
failures, if frequent enough to require another pivoting or ordering strategy; the cost of
the reservation on a plain solve, if it forces the reserve-a-list variant; the GPU test loop,
since CI has no GPU.

The smallest useful milestone is phase A restricted to `DistributedSlack` with limits and
local `ReactiveLimits`: it tells, on real grid snapshots, whether the one-analyze /
one-factorize premise holds and what the reservation costs, before committing to B and C.

## Open decisions

1. OpenLoadFlow's nested ordering, or the flat variant?
2. Reserve every reactive-limited PV bus (analysis guaranteed, J always larger), or reserve a
   list and allow a counted re-analysis?
3. Strict one-factorize, or the counted fallback?
4. Outer loops as the end state, or as the first step towards the projection rows of
   `nr_control_limits.md`, which reuse the same reservation and value-edit surface?
