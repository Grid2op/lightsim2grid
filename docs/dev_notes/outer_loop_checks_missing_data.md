# What lightsim2grid would need to detect an outer loop's trigger

**Status note, not documentation.** Written while adding the bus reactive-capability
check to the batch algorithms (`compute_bus_q_violations`, see
`src/core/batch_algorithm/BusQCheck.hpp`). It records, loop by loop, what a post-solve
"would this outer loop have fired?" check needs, what the `LSGrid` already holds, and what
is missing — so the next piece of this work does not have to re-derive it. It is not part
of the built documentation (`docs/*.rst`).

Reference for the outer loops: PowSyBl OpenLoadFlow,
`src/main/java/com/powsybl/openloadflow/lf/outerloop/config/DefaultAcOuterLoopConfig.java`
(the registration order *is* the nesting order, innermost first) and the loop classes in
`ac/outerloop/`.

## The short version

| OLF outer loop | Result side | Control side | Post-solve check possible? |
|---|---|---|---|
| `ReactiveLimits` | complete | complete | **yes — implemented** for a batch (`compute_bus_q_violations`), for all three families that hold a bus' voltage; a single solve already publishes `res_q_mvar` |
| `TransformerVoltageControl` | complete | **absent** | no |
| `PhaseControl` / `IncrementalPhaseControl` | complete | **absent** | no |
| `ShuntVoltageControl` | complete | **absent** | no |
| `DistributedSlack` | complete | complete | n/a — lightsim2grid distributes the slack inside the Newton (`MultiSlack`), there is no trigger to detect |
| `HvdcAcEmulationLimits` | complete | partial | limits are stored (`status_droop`); the droop itself is in the Newton (`Hvdc`) |
| `VoltageMonitoring` (SVC stand-by automaton) | complete | partial | `b_min` / `b_max` are now **checked** (as part of the reactive capability above) but still never enforced; the stand-by band itself is not modelled |
| `SecondaryVoltageControl` | complete | absent | no control zones / pilot points in the model |
| `AreaInterchangeControl` | complete | absent | no `Area` concept in the model |
| `AutomationSystem` | complete | absent | no overload-management systems in the model |

"Result side complete" everywhere is worth stating plainly: the powerflow results a check
would compare against are all published already — bus voltages, per-branch P/Q/I on both
sides, per-element reactive output, and even the inputs a limit check needs
(`limit_a1_ka` / `limit_a2_ka` per branch, `vmin_kv` / `vmax_kv` per bus). **Every gap
below is on the control-description side, and every one of them is lost in the
converters**, not in the solver.

## 1. `ReactiveLimits` — done

What it needs, and what answers it:

| needed | where it is |
|---|---|
| is this machine regulating voltage? | `is_voltage_controller` on its container (connected + regulating + not treated as off); `HvdcLineContainer::station_is_voltage_controller` for a converter station |
| its reactive capability | a generator: `min_q_` / `max_q_` (`GenInfo.min_q_mvar` / `max_q_mvar`). A converter station: the same, in MVAr (`get_station_min_q_mvar` / `get_station_max_q_mvar`). An SVC: `b_min_` / `b_max_` (`SvcInfo.b_min` / `b_max`), a **susceptance** range in pu, worth `b · \|V\|² · sn_mva` MVAr at the solved voltage |
| the reactive power it (or its bus) produced | `res_q_` (`res_q_mvar` on all three `Info` classes) after a solve; re-derived per row for a batch |
| which bus it holds | `regulated_bus_id_`, local or remote (a station always regulates the bus it stands on) |

**All three families are covered, and they are all of them**: exactly three kinds of element
have a reactive output the solver computes rather than reads — a voltage-regulating
generator, a voltage-regulating hvdc converter station, and a voltage-mode SVC. Everything
else standing on a bus injects reactive power that is *input* data, part of `Sbus`; its own
limits are an input question, not something a solve produced.

The SVC is the only one whose capability is not already in MVAr. `b_min` / `b_max` are
inputs (pu, base `sn_mva`) and the realized susceptance is the result, so the conversion is
one multiplication at the row's own voltage: a shunt `jb` consumes `-j·b·\|V\|²`, hence
injects `+b·\|V\|²` in generator convention, and a positive (capacitive) susceptance
produces reactive power. An asymmetric range (`[0, b]` versus `[-b, 0]`) is what pins that
sign down in the tests — a symmetric one cannot.

**The question is asked per bus, not per machine**, and that is not a shortcut. The
reactive power a bus needs is a fact about the solution; how it is divided between several
machines standing on that bus is not — the solver never decides it,
`LSGrid::_split_q_residual_per_bus` does, by a sharing convention (proportional to each
machine's reactive range). A per-machine check would therefore report the convention: two
20 MVAr machines covering 30 MVAr together is feasible, while the same solution read
machine by machine can show both 5 MVAr "over". So the check compares a bus' reactive
power against the **sum** of its machines' `[min_q, max_q]`.

The single-solve case needed nothing: `LSGrid::compute_results` publishes each generator's
reactive output, so summing it per bus is a comparison in Python. Note that
`LSGrid::check_solution(V, check_q_limits=True)` is **not** that check — it answers "is
this V a solution", and clamps per bus against ONE machine's limits at a time
(`check_solution_q_values_onegen`), which is the per-machine reading this note argues
against.

The batch case is what the new code does, because a batch keeps the voltages and drops
everything else. Per bus it re-derives the raw reactive residual — the algorithm's own
per-bus mismatch (`BaseAlgo::get_bus_mismatch`) with every `VoltageControl` controller's
own injection added back (`get_controller_q()`) — which is exactly what the machines
pinning that bus produced, whether they pin it through the classical PV path or through a
bordered control group.

It is also **a different kind of statement from a voltage or current violation**, and the
report says so: `ViolationCategory::PHYSICAL` vs `OPERATIONAL`. A bus outside its voltage
band or a line above its rating is a state the grid *can* reach and nobody wants to sit
in. A bus needing reactive power its machines do not have is a state the grid *cannot*
reach at all: the converged solution assumes a set-point that could not be held, so it is
a statement about the model's assumptions rather than about how the grid is operated. The
third category, `SOLVER` (`NOT_SIMULATED`, `DIVERGENCE`), is not a limit at all — a
divergence does not distinguish "no solution exists" from "this algorithm did not find
one".

What is still missing around it:

- **The PQ → PV direction does not exist.** OLF's loop also switches a bus *back* to PV
  when its voltage recrosses the set-point on the right side. lightsim2grid never pins a
  machine at a limit, so there is no such state to detect. This only becomes meaningful
  the day the limit is *enforced*.
- **Non-regulating machines are out of scope** by construction: their Q is their own
  setpoint, part of `Sbus`. A violation there is an input error, not something a solve
  produced. Same for a REACTIVE_POWER-mode SVC and a fixed-Q converter station.
- **A machine the algorithm never modelled is left out** rather than counted against a bus
  that never saw its reactive power: a remote-regulating generator absent from the
  controller list, which is what an algorithm that cannot do remote voltage control leaves
  behind.
- **`b_min` / `b_max` are checked, not enforced.** The check reports; nothing clamps an SVC
  to its susceptance range, exactly as nothing clamps a generator to its `[min_q, max_q]`.

## 2. `TransformerVoltageControl` (RTC) — the control data does not exist

`TrafoContainer` keeps `ratio_` and `shift_` and says so outright: *"lightsim2grid has no
'tap' concept"*. The tap is folded into the pi-model at load time and nothing about the
regulation survives.

| needed | present? |
|---|---|
| is there a ratio tap changer, and is it regulating? | **no** |
| its target voltage + deadband | **no** |
| the regulated terminal (which bus/side it watches) | **no** |
| tap position, low/high tap position, step | **no** — only the resulting `ratio` |
| the regulated bus' voltage | yes (`res_v1_kv` / `res_v2_kv`, bus V) |

Where it is lost:

- pandapower (`network/from_pandapower/_aux_add_trafo.py`): reads `tap_step_percent` and
  `tap_pos`, collapses them into `ratio = 1 + 0.01·step·pos` (`TrafoContainer.cpp`).
  `tap_min`, `tap_max`, `tap_side`, `tap_phase_shifter` are never read.
- pypowsybl (`network/from_pypowsybl/_aux_add_buses.py`, `_aux_trafo_rho`): takes the
  current `rho` from `get_ratio_tap_changers()`. `regulating`, `target_v`,
  `target_deadband`, `regulating_bus_id`, `low_tap_position` / `high_tap_position` are
  never read.

## 3. `PhaseControl` / `IncrementalPhaseControl` (PST) — same, plus the regulation mode

| needed | present? |
|---|---|
| is there a phase tap changer, and is it regulating? | **no** |
| its regulation mode (`CURRENT_LIMITER` / `ACTIVE_POWER_CONTROL`) | **no** |
| its regulation value (the P or I target) and the monitored terminal | **no** |
| tap position / range / step | **no** — only the resulting `shift_` (rad) |
| the monitored flow | yes (`res_p1_mw`, `res_a1_ka`) |
| a per-branch current limit | yes (`limit_a1_ka` / `limit_a2_ka`), but that is a *thermal rating*, not the loop's regulation value |

`get_phase_tap_changers()` is not read at all on the pypowsybl path. A
`CURRENT_LIMITER`-mode check is the closest thing to feasible today, and only by abusing
the thermal rating as the regulation value — which is not what OLF regulates against.

## 4. `ShuntVoltageControl` — a shunt is a fixed admittance

`ShuntContainer` stores `target_p_mw` / `target_q_mvar` and stamps
`{p, -q} / sn_mva` on the `Ybus` diagonal. There are no sections and no regulation.

| needed | present? |
|---|---|
| `section_count` / `max_section_count`, per-section b | **no** |
| `voltage_regulation_on`, `target_v`, `target_deadband`, regulated bus | **no** |
| the regulated bus' voltage, the shunt's own reactive output | yes (`ShuntInfo.res_v_kv`, `res_q_mvar`) |

Where it is lost: pandapower (`_aux_add_shunt.py`) reads `p_mw` / `q_mvar` only;
pypowsybl (`_aux_add_shunts.py`) reads `g` / `b` and scales them by the nominal voltage.
`get_shunt_compensators()`'s section and regulation columns are never read.

## The minimum addition that would unblock 2, 3 and 4

One **pure-input control descriptor per element**, which the solver never reads:

- `TapChangerData` on `TrafoContainer`, in two flavours — RTC (regulating flag, target V,
  deadband, regulated bus + side, tap position, low/high position, step) and PST
  (regulation mode, regulation value, monitored terminal, same tap fields);
- `SectionData` on `ShuntContainer` (section, max section, per-section b, regulating flag,
  target V, deadband, regulated bus).

Why it is cheap: it takes no part in `fillYbus` / `fillSbus`, raises no `AlgoControl` flag,
needs no `_on_*` hook, and so cannot break the invalidation contract or the refactorization
path the batch algorithms depend on. It does need the usual leaf paperwork — `StateRes`,
`get_state` / `set_state`, `save_binary` / `load_binary`, an `Info` class — and the two
converters have to read it (everything is available upstream: pypowsybl's tap-changer and
shunt-compensator frames, pandapower's `tap_min` / `tap_max` / `tap_phase_shifter` and
`shunt.step` / `max_step`).

Then the checks themselves are small, and they belong next to the reactive-capability one:
per element, "would have moved, in this direction, by this much, and has / has not a tap
left in that direction". Note that these three are **operational**, not physical: a tap
that should have moved and did not is a control that was not modelled, not an impossible
state — so they report in `ViolationCategory::OPERATIONAL`, unlike the reactive one.

## Why this is worth having beyond curiosity

`network/from_pypowsybl/_olf_params.py` currently **disables** these loops
(`_INLINE_MODE` / `_LEGACY_TRIGGER`) so that OLF solves the same single-shot,
outer-loop-free problem lightsim2grid solves. That is the right way to compare, but it
leaves one thing unknown: whether a given comparison case agrees because both sides were
prevented from acting, or because the loops were genuinely idle. A trigger check answers
exactly that, and is the natural companion to the comparison harness.
