# Enforcing control limits inside the Newton-Raphson

**Design note, not documentation. Nothing here is implemented.** It records a formulation,
why it fits the existing contracts unusually well, and what it would cost — so that the
question can be picked up later without re-deriving it. It is not part of the built
documentation (`docs/*.rst`). The corresponding entry is in `CHANGELOG.rst`'s `[TODO]`.

Companion to `outer_loop_checks_missing_data.md`, which asks the *detection* question —
what a post-solve "would this outer loop have fired?" check needs, and what the model is
missing for it. This note asks the next one. That note ends several sections on the same
sentence, in different words: *"`b_min` / `b_max` are checked, not enforced"*, *"This only
becomes meaningful the day the limit is **enforced**."* This is a proposal for that day.

## The short version

`docs/comparison_with_pypowsybl.rst` ends its "outer loops" section on an open question:

> Whether a given outer loop could instead be folded into the inner NR formulation (like
> distributed slack was) or fundamentally needs iteration around the solve (like discrete
> tap positions, which are not differentiable) is a case-by-case question, and remains
> open for most of OLF's outer loops as far as lightsim2grid is concerned.

It is largely **not** case-by-case. Four of the five are the same mathematical object:

| control | resource | box | regulated quantity |
|---|---|---|---|
| distributed slack | generator P | `[get_min_p, get_max_p]` | system active power balance |
| PV ↔ PQ | generator Q | `[min_q_, max_q_]` | bus voltage magnitude |
| ratio tap changer | ratio ρ | `[rho_min, rho_max]` | bus voltage magnitude (deadband) |
| area interchange | area generation | per-machine limits | area net export |
| secondary voltage control | group Q | per-machine limits | pilot bus voltage |

Every row is a **bounded control resource, complementary to a regulated quantity**:
the resource is strictly inside its box *if and only if* the regulated quantity is exactly
at its target, and a resource pinned at a bound comes with a known *sign* of the
regulation error. That is a complementarity condition, it can be written as an equation,
and equations can be handed to a Newton solver.

So the alternative to "an outer loop" is not "a hand-rolled in-Jacobian special case per
control". It is **one mechanism, instantiated five times**. Only tap *discreteness* (as
opposed to tap *limits*) genuinely needs iteration around the solve.

Caveat up front: none of this is benchmarked, and the cost section is the part to read
sceptically.

## What already exists

Worth stating precisely, because the picture changed recently and the obvious assumption
("lightsim2grid has no active power limits") is now wrong.

| | state |
|---|---|
| generator active limits | **present**, optional: `GeneratorContainer::set_p_limits`, `get_min_p` / `get_max_p` (NaN where never set), `GenInfo.min_p_mw` / `max_p_mw` |
| generator reactive limits | present: `min_q_` / `max_q_` |
| detection, slack past P limits | **done**: `batch_algorithm/GenPCheck.hpp`, via `compute_physical_violations` |
| detection, bus past Q capability | **done**: `batch_algorithm/BusQCheck.hpp`, per bus (see the sibling note on why per bus and not per machine) |
| enforcement, any of them | **absent** — nothing in any algorithm reads a limit |
| transformer voltage control | control side **absent**: no regulated bus, target or deadband stored for a trafo |
| area interchange | **absent**: no `Area` concept in the model |
| secondary voltage control | **absent**: no control zones, no pilot points |

`GeneratorContainer`'s own comment on `set_p_limits` states the gap exactly: *"the
distributed slack is solved INSIDE the Newton system (`MultiSlack`), with fixed
participation factors and no notion of a limit: a participating machine's converged active
power is `target_p + its share of the imbalance`, which can land anywhere."*

The last three rows matter for sequencing: for tap control, area interchange and secondary
voltage control, the **model** is missing before the formulation is. No amount of Jacobian
work helps until a transformer can say which bus it regulates. The first two rows are
where this note has immediate purchase.

## Why the heuristics cycle

Worth spelling out, because it is the whole argument for doing this at all.

The classical PV → PQ heuristic uses **two different tests**:

- switch PV → PQ when `Q > max_q`;
- switch back PQ → PV when `Vm > v_set`.

Each is correct in isolation. The defect is that they are evaluated on *different
iterates*, with a discrete bus-type state carried between them. Nothing forbids the solver
satisfying test 1 at iteration k, test 2 at k+1, test 1 again at k+2, indefinitely. No
single iterate is inconsistent; the *sequence* is.

With two coupled controls it gets qualitatively worse, because there is no correct **order**
to run them in. A tap step changes `Vm`, which un-pins a generator, whose `Q` then exceeds
`max_q`, which re-pins it and frees `Vm`, which moves the tap. Running the tap logic before
or after the reactive-limit logic gives *different* limit cycles, so the answer depends on
an arbitrary implementation choice.

OpenLoadFlow resolves this by fiat, which is the tell: `transformerVoltageControlMode`
offers `AFTER_GENERATOR_VOLTAGE_CONTROL` ("a continuous voltage control is performed after
each generator voltage control outerloop") alongside the default
`CONTINUOUS_WITH_DISCRETISATION`. A parameter that lets the caller choose the interleaving
is an admission that the interleaving matters and that no ordering is right.

In the formulation below the question **cannot be asked**: every control is a row of the
same Jacobian, solved simultaneously by one Newton step. There is no "after".

## The formulation

### NCP functions

An **NCP function** is any `φ: R² → R` with

```
φ(a, b) = 0   <=>   a >= 0,  b >= 0,  a·b = 0
```

Two matter in practice:

```
φ_min(a,b) = min(a, b)                     the natural residual
φ_FB (a,b) = a + b - sqrt(a² + b²)         Fischer-Burmeister
```

`φ_min` is piecewise linear and reproduces exactly the clamp / un-clamp logic one would
write by hand. `φ_FB` is smooth everywhere *except* at `(0,0)`, which matters for
globalisation.

For a **two-sided box** — which is what every control in the table has — the convenient
form is the projection (normal-map) residual. With resource `u`, box `[u_min, u_max]`,
regulation error `e` (target minus actual), `σ = sign(∂(regulated quantity)/∂u)` and any
`c > 0`:

```
Φ = u - mid(u_min, u_max, u + c·σ·e)
```

`mid` is the median of three, i.e. the projection onto the box. Writing `z = u + c·σ·e`:

| branch | Φ reduces to | meaning |
|---|---|---|
| `u_min < z < u_max` | `-c·σ·e` | **regulating** (`e = 0`) |
| `z >= u_max` | `u - u_max` | pinned at the upper bound |
| `z <= u_min` | `u - u_min` | pinned at the lower bound |

**That row is the bus type.** There is no PV / PQ label anywhere in the system.

### Why this kills the cycling

Take the generator case (`u = Q`, `e = v_set - Vm`) and check the two classical tests
against `z = Q + c·(v_set - Vm)`:

- while regulating, `Vm = v_set`, so `z = Q` — pinning when `z > max_q` **is exactly test 1**;
- while pinned, `Q = max_q`, so `z = max_q + c·(v_set - Vm)` — releasing when `z < max_q`
  **is exactly test 2**.

The projection reproduces both classical tests, exactly, but as **one function of one
iterate**. There is no state carried between tests to become inconsistent. Add the line
search below and the iteration is monotone in a merit function, so it cannot return to a
state it has left. Cycling is not suppressed by a latch or an iteration counter — it
becomes structurally impossible.

### Semismooth Newton

`Φ` is not differentiable, so ordinary Newton does not apply. Semismooth Newton does.

A locally Lipschitz, directionally differentiable `F` is *semismooth at x* when, for
**every** `V` in the Clarke generalized Jacobian `∂F(x + d)`,

```
F(x + d) - F(x) - V·d = o(||d||)        as d -> 0
```

and *strongly* semismooth when that is `O(||d||²)`. The quantifier is the content: the
bound must hold for every choice of generalized Jacobian element, which is what lets an
implementation pick an arbitrary one.

The condition is weak in practice. Piecewise smooth functions are semismooth; `min`, `max`,
`|·|`, `mid` and the Euclidean norm are *strongly* semismooth; sums and compositions
preserve it. So the full residual — smooth power-flow equations plus projection rows — is
strongly semismooth.

**Theorem (Qi & Sun 1993).** If `F(x*) = 0`, `F` is semismooth at `x*`, and every
`V ∈ ∂_B F(x*)` is nonsingular (*BD-regularity*), then

```
pick   V_k in the B-subdifferential of F at x_k   (any element)
solve  V_k · d_k = -F(x_k)
       x_{k+1} = x_k + d_k
```

converges Q-superlinearly from any `x_0` near `x*`, Q-quadratically if `F` is strongly
semismooth. One linear solve per iteration — the same cost as ordinary Newton.

### Globalisation

Fischer's result (1992): `ψ_FB = ½·φ_FB²` is **continuously differentiable on all of R²**,
including at the origin where `φ_FB` is not. So

```
Ψ(x) = ½·||Φ(x)||²
```

is C¹ even though `Φ` is not, and a standard Armijo backtracking line search on `Ψ` applies
to the semismooth Newton direction. That is the globalisation, and it is what makes cycling
impossible rather than merely unlikely.

`Ψ` can have stationary points that are not solutions — the complementarity literature
rules these out by assuming the underlying map is P₀, which AC power flow is not. In
practice: line search, with a fallback to an undamped semismooth step when it stalls. Same
posture as the step damping the Newton-Raphson already uses.

## Instantiation per control

### Distributed slack with min_p / max_p

`MultiSlack` (`src/core/powerflow_algorithm/NRSystem.hpp`) is the textbook bordered
formulation: one scalar unknown `slack_absorbed` (call it `k`), one column, one P equation
per slack bus, and feature entries `w_b` at each slack bus' P row.

**Mind the sign.** With the convention in the code (`mis += k * slack_weights` in
`MultiSlack::adjust_mismatch`, and `LSGrid::compute_results` subtracting
`slack_absorbed * slack_weights` back out), the effective injection is `Sbus - k·w`, i.e.
**`k` is negative when generation must rise**. It reads backwards and deserves a comment.

Add one non-negative *shed* variable `s_g` per slack generator:

```
P_g = P_g0 - w_g·k - s_g
0 <= s_g  ⊥  (max_p_g - P_g) >= 0
```

`s_g = 0` with positive headroom is a machine following the slack; `s_g > 0` with zero
headroom is a machine pinned at `max_p` while the rest of the fleet takes its share. A
symmetric block handles `min_p`.

`Φ_g` depends only on `s_g` and `k`, so each row has exactly **two** nonzeros. At bus level
`adjust_mismatch` gains one term, `mis += k·w + E·s`, with `E` the generator-to-bus
incidence (one nonzero per column).

Two points specific to the slack:

- **Granularity.** Limits are per generator; the Jacobian weight is per *bus*
  (`get_slack_weights_solver`), and the generator-level split is post-hoc and proportional
  (`set_p_slack`). With several slack generators on one bus, pinning one means recomputing
  `w_b` from the survivors — so `bus_slack_weight_` becomes an NR-iteration quantity, not a
  per-solve one.
- **Keep-one rule.** If every slack generator pins, the `k` column goes to zero and the
  bordered system is singular. BD-regularity failing *is* the infeasibility detector:
  report it rather than diverging.

This is the one place where detection already exists (`GenPCheck.hpp`), so the enforcement
work has a ready-made oracle: every row `compute_physical_violations` currently flags is a
row this should silence.

### State-dependent slack weights (OLF `balanceType`)

Two different things get conflated here. Using `κ = -k` so that `κ > 0` means "generation
must rise":

**Direction-dependent weights** — different participating sets and factors for up- and
down-regulation — are a kink at `κ = 0` in a function of an unknown. Piecewise smooth,
hence semismooth: evaluate `w` at the current iterate, take either branch at the kink. No
new unknowns, no new sparsity. Essentially free once the rest exists.

**Margin-proportional weights** (OLF's `PROPORTIONAL_TO_GENERATION_REMAINING_MARGIN`) are
implicit: `w_g ∝ max_p_g - P_g`, and `P_g` depends on `w_g`. OLF resolves it by freezing `w`
from the current `P` at the top of each outer iteration. In-NR, make the margin an unknown
with its own local row:

```
P_g = P_g0 + κ·m_g
m_g = max(0, max_p_g - P_g)
```

Two nonzeros per row. Solving the row on the active branch:

```
m_g = (max_p_g - P_g0) / (1 + κ)
P_g = P_g0 + κ·(max_p_g - P_g0) / (1 + κ)
```

Monotone in `κ`, and tends to `max_p_g` **without ever exceeding it**. So
margin-proportional weighting *is* the "reduce the weight as the machine approaches its
limit" idea, and unlike a fixed-steepness sigmoid it enforces the limit exactly, because
the taper is self-consistent rather than imposed. No complementarity needed for this
variant.

Caveats, in decreasing order of nuisance:

- **Conditioning.** Total reachable imbalance is `Σ(max_p_g - P_g0)·κ/(1+κ)`, so exhausting
  the fleet's headroom needs `κ → ∞`. An infeasible case diverges in `κ` rather than failing
  cleanly: cap it explicitly.
- **It is a different physical model.** Margin-proportional means every machine approaches
  its limit asymptotically *together*; AGC participation factors mean machines pin *one at a
  time*. Both are defensible and they give different answers. Do not let "easier to
  formulate" pick the model.
- `max(0, ·)` still matters for a machine starting above `max_p`.

**Do not normalise inside the NR.** `GeneratorContainer::get_slack_weights_solver` does
`res /= sum_res`. Carried into a state-dependent weight, `w_g = m_g / Σ_h m_h` makes every
`∂w_g/∂m_h` nonzero and the slack border **dense** across all participants — which destroys
the batch algorithms' premise. The normalisation is a *gauge*: `κ·m_g/Σm` is just `κ'·m_g`
with `κ' = κ/Σm`, and `κ` is a free unknown, so the absolute scale of `w` is unobservable
inside the NR — only the ratios are. Drop it from the solver and keep it only in
`set_p_slack`'s proportional share, where the sum is a constant computed after the fact.

### Hydro with produce / absorb modes

A machine with a "produce" mode (`0 < P < max_p`) and an "absorb" mode (`min_p < P < 0`),
where switching mode *during* a solve is not wanted.

**If the mode is frozen** — decided by the dispatch, which is the physically defensible
position, since whether a hydro unit turbines or pumps is a commitment decision and not
something a power flow should decide — then the mode boundary is **just another limit**: the
box is `[0, max_p]` in produce mode and `[min_p, 0]` in absorb mode. A produce-mode machine
the slack wants to push negative simply pins at `P = 0`, its weight goes to zero, and its
share is redistributed. The "don't switch" guarantee is not a guard bolted on; it is a
property of the system being solved.

**If the solve should pick the mode**, split the injection: `P = p⁺ - p⁻` with
`0 <= p⁺ ⊥ p⁻ >= 0`. Representable — but this is an MPEC-style complementarity between two
*decision* variables, degenerate exactly at `P = 0` where both vanish, constraint
qualifications failing generically and BD-regularity with them. Machines idling near zero
will be where the solver sticks. Representable, yes; well behaved, no.

**If the mode has a true deadband** (produce mode really starts at some `min_p_turbine > 0`),
the feasible set has a *hole* and is not NCP-representable at all. Freeze the mode, treat
the deadband edge as the machine's limit, and report a commitment infeasibility upward
rather than letting the power flow shut the unit down.

### PV ↔ PQ reactive limits

This is the **original** application of the technique — every reference at the end of this
note is a reactive-limit paper. With resource `Q_j`, box `[min_q, max_q]`, error
`e = v_set - Vm_r`:

```
Φ_j = Q_j - mid(min_q_j, max_q_j, Q_j + c·(v_set - Vm_r))
```

For a group of N machines sharing Q proportionally, the sharing rows take the same shed
treatment as the distributed slack: `Q_j = α_j·Q_grp - t_j` with
`0 <= t_j ⊥ (max_q_j - Q_j) >= 0`, and the group's voltage row regulates while any member
has headroom, releasing `Vm_r` once all are pinned.

Note the interaction with `BusQCheck.hpp`'s central argument, which is that the *detection*
question is a per-**bus** one because the per-machine split is a convention
(`LSGrid::_split_q_residual_per_bus`) and not something the solver decides. Enforcement is
the opposite: it is genuinely per machine, because a pinned machine stops following the
group and the split stops being a convention. So this formulation would *create* the
per-machine truth the check currently cannot assume — which is a point in its favour, and
also means the two must not be expected to agree machine-by-machine before it exists.

### Ratio tap changers

Identical structure, with `u = ρ`:

```
Φ = ρ - mid(rho_min, rho_max, ρ + c·σ·(v_set - Vm_t))
```

**`σ = sign(∂Vm_t/∂ρ)` is not optional and not guessable**: it depends on the ratio
convention and on `is_tap_side1` (`TrafoContainer.hpp`). Getting it backwards selects the
*opposite* branch — the tap runs to the wrong rail and the solve converges, silently, to a
wrong answer. Compute the sensitivity from the model; do not hard-code a sign. This is the
class of bug the `LocalBusId` / `GlobalBusId` / `SolverBusId` tagging exists to prevent
elsewhere in this codebase.

**The deadband makes this better posed, not worse.** A real OLTC regulates to a band, not a
point, and a point target is the ill-posed idealisation — the device moves in steps and can
never sit exactly at `v_set`. The band is also where an otherwise structurally singular `ρ`
column gets an equation, and the device's own answer is a perfectly good one:

| condition | row |
|---|---|
| `Vm_t < v_lo` | ρ rises, up to `rho_max` |
| `Vm_t > v_hi` | ρ falls, down to `rho_min` |
| `v_lo <= Vm_t <= v_hi` | `Φ = ρ - ρ0`, i.e. pin at the incoming tap |

Five branches instead of three, still piecewise linear, still strongly semismooth — and the
deadband branch gives the iteration somewhere to *rest*, which is exactly what the cycling
heuristic lacks.

**Discreteness is the one genuine wall.** Complementarity handles *bounds*, not
*integrality*; no NCP function represents a lattice. The practical answer is the universal
one: solve the continuous relaxation, round to the nearest position, re-solve once with
taps fixed — which is what both OLF modes do.

Here lightsim2grid is unusually well placed: `TrafoContainer.hpp` states that "there is NO
notion of a discrete `tap` in lightsim2grid" — the state is a continuous `ratio_` and the
tap-position arithmetic lives in the converters. So the continuous relaxation is already the
native representation, and the rounding naturally belongs upstream in the pypowsybl path
where tap positions actually exist.

**But the control side does not exist yet.** Per the sibling note, a transformer stores no
regulated bus, no target and no deadband — they are lost in the converters. That is
prerequisite zero here, and it is converter work, not solver work.

### Area interchange and secondary voltage control

Sketched only, because they add no new mathematics and are blocked on the model, not the
formulation:

- **Area interchange**: one shared scalar per area (as `k` is for the system), participating
  machines with limits, regulated quantity = the area's net export. Structurally the
  distributed slack, replicated per area. Blocked: no `Area` concept in the model.
- **Secondary voltage control**: a pilot bus regulated by a group of machines sharing Q by
  participation, each with reactive limits. Structurally the `VoltageControl` group, with the
  shed variables of the reactive-limit section. Blocked: no control zones, no pilot points.

## Fit with the existing architecture

### The batch premise stops being a constraint

`CLAUDE.md` states the rule: *"what is genuinely forbidden is changing the pattern itself
per row — a different pv/pq vector handed to the solver, a different slack set. That raises
`has_pv_changed()` / `has_slack_participate_changed()`, forces a fresh `analyze` on every
row, and the batch loses its reason to exist."* `set_switchable_vm_buses` +
`set_pv_pinned_buses` exist to work around exactly that: reserve the union of layouts, mask
the difference by value.

In this formulation **there is no pv/pq vector**. Every regulated bus permanently owns a Vm
unknown and a Q equation; every controller permanently owns a Q unknown; "PV-ness" is a
*value* of the projection row, recomputed each iteration. So:

- `set_switchable_vm_buses` becomes unnecessary — every regulated bus is switchable by
  construction;
- `set_pv_pinned_buses` becomes unnecessary — pinning is what the projection row does;
- `has_pv_changed()` can never fire from a reactive-limit event, because there is nothing to
  change.

The workaround stops being a workaround and becomes the formulation. Everything that changes
when a control pins or releases is a **value**; positions are fixed at `register_in`. One
`analyze`, `refactorize` thereafter, for every row of a sweep. Turn `set_refactor_fallback`
on for the usual reason — those value edits move pivots.

### `VoltageControl` already has most of the slots

Reading `NRSystem.hpp`, the reactive-limit case needs strikingly little new structure.
`VoltageControl::register_in` already claims, per group:

- `q_cols_[j] = ledger.add_q_unknown(...)` — Q is **already an explicit unknown**;
- `v_rows_[g] = ledger.add_custom_row()` — the voltage constraint already owns a row;
- `h_vm_[g] = sink.add(v_row, vm_cols_[g])` and `h_slope_[j] = sink.add(v_row, q_cols_[j])`
  — **both Jacobian slots the three-way branch needs are already reserved**.

And the `stranded` path is a hand-rolled instance of exactly the row swap this needs.
`fill_feature_values` writes `h_slope_[j] <- stranded ? 1. : slope(j)` and
`h_vm_[g] <- stranded ? 0. : 1.`, while `fill_custom_rows` swaps the residual from
`Vm + Σ s·Q - v_set` to `Q_c`. That is the voltage row being turned into a "`Q_c = 0`" pin
**by value alone, inside a fixed sparsity pattern**. Replace `0` with `max_q` and make the
branch selector the projection test instead of the stranded flag, and the formulation is
there.

For a lone controller this costs **zero new rows, zero new columns, zero new nonzeros** —
the slot exists. The only gating change is the condition that currently reserves
`(v_row, q_col)` for SVCs and for `may_mask_` singletons only; it would become unconditional.

### `Hvdc` is the precedent for a state-dependent Ybus

A continuous ρ makes the transformer's `Ybus` entries state-dependent, and `NRSystem`'s
generic dS pass assumes a fixed `Ybus`. Two ways out:

1. recompute the 2×2 block each iteration and raise `tell_recompute_ybus()` — correct, but it
   fights the whole `AlgoControl` caching design;
2. **stamp the transformer at a nominal ratio ρ0 and let an extension own the difference** as
   an injection correction, `ΔS = V ⊙ conj((Y(ρ) - Y(ρ0))·V)`, applied in `adjust_mismatch`.

Option 2 is not speculative: it is what `Hvdc` already does. `HvdcLineContainer.hpp`
documents it — *"For droop lines, this container does NOT stamp the AC active power in Sbus:
the NR system Hvdc extension owns it (all regimes)."* Same trick one level up. The container
stamps the constant part, the extension owns the state-dependent part, `Ybus` stays constant
and its sparsity untouched, and the correction's Jacobian contributions land in dS positions
that already exist because the branch is already there at ρ0.

Per regulating transformer: **1 column** (ρ), **1 custom row** (2 nonzeros), up to 4 entries
in the ρ column (P and Q rows of both ends). Negligible.

A `TapControl` extension beside `MultiSlack`, `Hvdc` and `VoltageControl` is the obvious
home — the component protocol documented at the top of `NRSystem.hpp` (`update_state` /
`init_topology` / `register_in` / `declare_feature_entries` / `fill_feature_values` /
`adjust_mismatch` / `fill_custom_rows` / `apply_step` / `clear`) already has every hook this
needs.

## Costs and open questions

Stated plainly, because they are the reason this is not obviously the right thing to do.

- **Dimension growth on the common grid.** Today a PV bus has *no* Vm column and *no* Q row.
  Here it has both, plus a Q unknown and a voltage row. On a transmission grid with ~10-15%
  regulated buses, `J` grows meaningfully in dimension and more in nonzeros (the new Q rows
  drag in `dS_dVm` entries). This is paid on **every** solve, including the majority that
  never touch a limit.

  Note where it lands, though: the layout is *identical to the union layout already built*
  when `set_switchable_vm_buses` is populated. So the cost is real for a plain `ac_pf` and
  **zero for the batch algorithms** — which is also where a cycling heuristic hurts most (a
  switching loop on row 4000 of a `ScenarioSweep` costs iterations across the whole sweep and
  is miserable to debug). Keeping the classical pv/pq split for `ac_pf` and using projection
  rows in the batch path is a defensible split.

- **Scaling.** The constant `c` mixes two physical quantities. `max_q - min_q` can be 10 pu
  while the voltage deviation of interest is 0.05 pu; per-machine normalisation
  (`φ(a/σ_g, b/σ_g)`, `σ_g` a rating or a span) is not optional. Textbook implementations
  omit this and fail on real grids for this reason alone. The trap appears once per control
  type.

- **Degeneracy.** BD-regularity can fail exactly at a bound with zero slack — a machine
  dispatched at `max_q` with `Vm` exactly at setpoint, which is not exotic. The cure is
  Kanzow's smoothing, `φ_μ(a,b) = a + b - sqrt(a² + b² + 2μ²)`, solved for a sequence
  `μ_j → 0` with warm starts. This *is* the sigmoid-smoothing idea, but with a homotopy that
  drives the smoothing to zero, so the converged answer respects the limit exactly instead of
  landing at an arbitrary point inside it.

- **Non-uniqueness.** Complementarity fixes the *path*, not the *solution set*. AC power flow
  with voltage controls still has multiple solutions and one can still converge to a
  low-voltage one.

- **Discreteness.** Taps and switched shunts still need round-and-re-solve. This is the only
  part of the outer-loop question that genuinely needs iteration around the solve.

- **Unvalidated here.** None of this has been measured on this codebase. The claims about
  iteration counts and about batch-path cost are structural arguments, not benchmarks.

- **Why nobody does it.** Worth knowing, since "no production load flow works this way" is a
  fair objection. Semismooth Newton is 1993; NR power flow is 1967, and by the time the theory
  existed every utility-grade solver was a validated artefact with regulatory reproducibility
  requirements. The heuristics also work most of the time, and a *naive* implementation —
  without the scaling above — is worse than a heuristic, so early attempts plausibly failed
  quietly. It is not unexplored, though: GRIDOPT / PFNET and the CMU equivalent-circuit line
  both do it.

## A possible incremental path

Each step is independently useful and independently testable. They are the same algorithm at
rising levels of rigour, so this is a climb rather than a commitment.

1. Enforce `get_min_p` / `get_max_p` on the distributed slack by discrete
   clamp-and-renormalise, reserving the shed rows / columns up front. This *is* `φ_min`
   semismooth Newton with a unit step and no merit function.
   `BaseBatchSweep::_masked_slack_weights` already does the zero-and-renormalise on `w`; reuse
   it per NR iteration instead of per row. `GenPCheck` gives the oracle: rows it flags today
   should stop being flagged.
2. Add the Armijo line search on `Ψ = ½||Φ||²`. Cycling disappears. This is now globalised
   semismooth Newton.
3. Apply the projection row to reactive limits, gated to the batch path first (where the
   layout is already reserved). Flip the `(v_row, q_col)` reservation to unconditional. Expect
   `BusQCheck` results to change meaning — see the note in that section.
4. Swap `φ_min → φ_FB` where degenerate machines show up, with a μ-homotopy if needed.
5. `TapControl` extension — but only after the converters carry a transformer's regulated bus,
   target and deadband, which is the actual blocker.

Step 1 alone closes the `min_p` / `max_p` gap on the distributed slack, which is the smallest
useful deliverable and the one with detection already in place.

## References

Nonsmooth analysis and semismooth Newton:

- R. Mifflin, *Semismooth and semiconvex functions in constrained optimization*, SIAM J.
  Control Optim. 15 (1977) — the original definition.
- L. Qi & J. Sun, *A nonsmooth version of Newton's method*, Math. Programming 58 (1993)
  353-367 — the convergence theorem.
- A. Fischer, *A special Newton-type optimization method*, Optimization 24 (1992) 269-284 —
  the FB function and its C¹ merit function.
- T. De Luca, F. Facchinei & C. Kanzow, *A semismooth equation approach to the solution of
  nonlinear complementarity problems*, Math. Programming 75 (1996) 407-439 — the practical
  algorithm, line search included.
- F. Facchinei & J.-S. Pang, *Finite-Dimensional Variational Inequalities and Complementarity
  Problems*, Springer 2003, vol. II ch. 7-9 — the reference text.

Power systems:

- W. Murray, T. Tinoco De Rubira & A. Wigington, *A robust and informative method for solving
  large-scale power flow problems*, Comput. Optim. Appl. 62 (2015) — complementarity for
  reactive limits and PV/PQ switching; the GRIDOPT / PFNET line of work.
- *A modified Newton-Raphson load flow scheme for directly including generator reactive power
  limits using a complementarity framework*, Electric Power Systems Research 105 (2013).
- S. V. Dhople, Y. C. Chen, A. Al-Digs & A. Dominguez-Garcia, *Reexamining the Distributed
  Slack Bus*, IEEE Trans. Power Systems 35(6) (2020) 4870-4879 — what `P0`, `α` and the
  imbalance are *supposed* to mean. Relevant because `GeneratorContainer` derives weights from
  `|target_p|`, which is a heuristic and not an AGC participation factor: enforcing limits on
  top of arbitrary participation makes the *limits* exact and the *distribution* still
  arbitrary.
- A. Pandey, M. Jereminov, G. Hug & L. Pileggi, *Continuously Differentiable Analytical Models
  for Implicit Control within Power Flow*, arXiv:1811.02000 (2018), and A. Agarwal,
  A. Pandey, M. Jereminov & L. Pileggi, *Implicitly Modeling Frequency Control within Power
  Flow*, ISGT-Europe 2019 (arXiv:1908.11778) — the sigmoid-smoothing alternative, motivated
  explicitly by removing the outer loop.
