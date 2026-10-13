# Speed optimisation of the Newton-Raphson with outer loops

**Design note, not documentation.** It records where the time goes in a batch run with
OpenLoadFlow's outer loops (`NROuter_*`), what has been done about it, and what is left. It is
not part of the built documentation (`docs/*.rst`). Unlike the documentation, it keeps the
measured numbers: they are what the next step is decided on.

Companion notes in this folder:

- `static_pivot_refactorization.md` -- the design of `KLUStaticPivotLinearSolver` (the
  `*_KLUStaticPivot` algorithms) and its literature.
- `outer_loops_fixed_sparsity.md` -- why the outer loops keep one Jacobian pattern per solve
  and change values only, which is where the refactorization failures below come from.

The layout facts below (which row lands on which diagonal) were read off `NRLedger`,
`VoltageControl::register_in` and `BranchControl::register_in`, and the zero pivots
classified with temporary instrumentation of those and of the solver, since removed.

## The bench

Everything below is measured the same way, so that the numbers compare:

- one large real grid snapshot, `NROuter_*` with all seven OpenLoadFlow outer loops,
  `tol = 1e-8`, `max_iter = 50`, a DC start, the process pinned to one core (`taskset -c 0`);
- **ContingencyAnalysis**: 100 line N-1 drawn with a fixed seed among the lines closed at both
  ends, `init_from_n_powerflow = True`, `handle_disconnected_grid = True`;
- **InjectionSweep**: 100 rows, every generator's target P scaled by `1 + 0.001 * (k - 50)` and
  every load's by slightly more (a small imbalance for DistributedSlack),
  `init_from_n_powerflow = True`;
- compared: the wall time of `compute()`, the converged rows, the outer and Newton iteration
  counts, the linear-solver counters, and the per-row maximum of dVm / dVa (median, p90, max).

The bench scripts read the snapshot and so live outside the repository.

## Static-pivot refactorization (done)

`NROuter_KLUStaticPivot` against `NROuter_KLU`, same build:

| | `NROuter_KLU` | `NROuter_KLUStaticPivot` |
|---|---|---|
| ContingencyAnalysis, 100 N-1 (3 runs) | 2.98 / 2.94 / 2.93 s | 2.09 / 2.08 / 2.09 s (-29 %) |
| InjectionSweep, 100 rows | 8.00 s | 5.72 s (-28.5 %) |

Converged rows, outer and Newton iterations, and the linear-solver counters (1 analyze, 143
factorizations of which 142 fallbacks, 991 refactorizations for the contingencies) are
identical. Voltages differ at round-off: per-row max dVm 9.7e-13 / 1.6e-12 / 2.1e-12 pu and
dVa 4.1e-11 / 6.1e-11 / 9.6e-11 deg (median / p90 / max) for the contingencies, the same order
for the sweep.

One lesson from getting there: the reciprocal pivot growth guard first computed KLU's
definition literally, on the row-scaled matrix (a pass over the matrix for the row scales,
then a scaled and branchy max on every entry of B and U). That cost about a third of the
kernel's time and held the gain at -17 %. Computed on the matrix the solver actually
factorizes, with a branchless max on values the kernel reads anyway, it costs nothing
measurable (2.06 s with the guard, 2.08 s without).

## Where the time goes now

Callgrind instruction counts (`Ir`) of the ContingencyAnalysis `compute()` alone (collection
toggled on that function), release build:

| | `NROuter_KLU` | `NROuter_KLUStaticPivot` |
|---|---|---|
| total | 29.8 G | 24.6 G |
| refactorization | 14.2 G (`klu_refactor` 9.8, its `klu_scale` ~2.1, `klu_rgrowth` 2.3) | 9.4 G (kernel 8.2, plan freeze + `klu_extract` 1.1) |
| fallback `klu_factor`, 142 calls | 4.4 G | 4.4 G |
| triangular solves | 2.3 G | 2.0 G |
| Newton system (fill J, residual, Ybus . V, sincos, mask positions) | ~4.9 G | ~4.9 G |
| `memset` / `malloc` / `free` | ~1.1 G | ~0.9 G |

Shares of the static-pivot run:

| | share |
|---|---|
| refactorization kernel (991 calls, ~8.2 M `Ir` each) | 33 % |
| fallback factorizations (`klu_kernel`, `klu_scale`, ...) | 18 % |
| plan freeze after each fallback (`freeze_plan` + `klu_extract`) | 4.5 % |
| Newton system | ~20 % |
| triangular solves | 7 % |
| `memset` / allocations | ~3.5 % |

Wall time splits the same way (instrumented run, 2.09 s in all): 849 successful
refactorizations 0.58 s, 143 `klu_factor` 0.48 s, 143 plan freezes 0.12 s, 142 rejected
refactorizations 0.03 s.

### Why there are 142 fallbacks

Every one of them is a refactorization whose frozen pivot became exactly zero. Instrumented:

- 141 of 142: the entry the pivot sits on was written to 0 by a value-level edit of the
  Jacobian (the 142nd pivot sits on a fill-in entry). None comes from cancellation.
- 120 of 142 are on the diagonal (`P[k] == Q[k]`), while KLU pivots off the diagonal in under
  2 % of the columns.
- The bus masking of `handle_disconnected_grid` and the PV pinning (`_recompute_mask_positions`)
  never zero a diagonal entry: they are not the cause.
- What is left are the rows the outer loops switch, by value, between two forms. Each zero
  pivot was classified by the row and the column it sits on:

  | row | column | zero pivots |
  |---|---|---|
  | a voltage-control group's voltage row (remote control by a generator) | the group's first `Q_c` column, its diagonal | 64 |
  | a transformer voltage-control ratio-group row (BranchControl) | a ratio column, on the diagonal | 56 |
  | the same rows | a ratio column, off the diagonal | 21 |
  | the Q row of a voltage-control controller bus | a Vm column, off the diagonal | 1 |

  A voltage row reads "Vm(reg) = target" (one entry, on the regulated bus's Vm column) or,
  held, "Q_c = q" (one entry, on the `Q_c` column); a ratio-group row reads "Vm(bus) = target"
  or "ratio = held", the same way. Each form zeroes the other's entry.
- Where the pivot sits follows from the layout. Rows and columns are numbered by two
  counters that advance together (`NRLedger`), so a bus's Vm column has the bus's Q row on its
  diagonal (dQ/dVm, never switched), while a voltage row lands on the diagonal of its group's
  first `Q_c` column (every group of the snapshot). KLU prefers the diagonal, so it pivots
  the voltage row on `(v_row, Q_c)`: 1 in the held form, 0 in the PV form.
- In every case the Vm columns that would let the pivot fill in come **after** the pivot
  column in KLU's order: for the 64 voltage rows, the regulated bus's Vm column and the
  controller bus's Vm column; for the 77 ratio-group rows, the group bus's Vm
  column. That is what turns the other form into a zero pivot, see "Ordering the regulated Vm
  first" below. It is the necessary condition only: that moving them first makes the pivot
  fill in is the reasoning below, not yet measured.
- Local voltage control by a generator does not go through these groups (a local PV bus is
  the ordinary PV path), so none of the cases is local.
- The pivot orders recur: of the 143 factorizations, 73 chose an order already seen (70
  distinct orders). The sets of exactly-zero entries do not recur (136 distinct out of 143:
  each contingency zeroes its own branch's entries), so they cannot key a cache.

## What is left, by expected gain

1. **Avoid the fallback factorizations** (up to ~24 % of `Ir`: factorization, freeze, rejected
   attempts). The factorization itself is needed only to choose new pivots; when the order it
   would choose was seen before, it can be reused.
   - Keep the last few frozen plans. Before refactorizing, check that every pivot entry of the
     current plan is nonzero in the new values (one load per column, the positions precomputed
     at freeze); if not, take the most recent cached plan that passes, and call `klu_factor`
     only when none does. Upper bound from the recurrence above: about half of the
     factorizations and freezes, ~8 % of `Ir`. Changes the fallback counters; the Newton
     iteration counts must be re-checked, since a cached order was chosen on other values (the
     growth guard stays on).
   - The precheck alone, with `klu_factor` called directly instead of a refactorization that
     fails part-way, saves only the rejected attempts (~1.5 %).

   **Detecting the zero pivot before refactorizing.** The outer loops do know when they switch
   a row (PV -> PQ, a controller held), so the idea is to tell the solver and factorize
   directly. Three points settle how to do it:

   - *Detection alone gains little.* It skips the refactorization that fails part-way (0.03 s
     of 2.09 s, ~1.5 %); the `klu_factor` that follows, which is the 18 %, is still paid.
   - *The outer loop knows the equation, not the pivot.* Whether the entry a switch zeroes is
     a pivot depends on the frozen order, which only the solver holds. The solver can check
     it on its own: store, at freeze, the position in J's value array of every pivot entry
     that is an entry of J, and before each refactorization read those values (one load per
     column, a fraction of a percent of a refactorization). That needs no new interface and
     catches every source of such edits at once (bus masking, PV pinning, every outer loop).
     Passing the switches down instead would cross `NROuterInner` -> `NRAlgo` -> the linear
     solver, and `BaseAlgo` if exposed there, which is the plugin ABI.
   - *No single order can serve both forms of a switched row* (see above): the pivot has to
     move when the row switches. The gain is therefore not in detecting the switch but in not
     paying `klu_factor` to find the new order, i.e. in combining the precheck with the cache
     of plans: precheck the current plan; if it fails, take the most recent cached plan whose
     pivot entries are all nonzero; call `klu_factor` only when none passes. A cache of 4 to 8
     plans, ~8 % of `Ir` expected.
   **Ordering the regulated Vm first (orthogonal to the above, and probably better).** A
   frozen order *can* serve both forms of a switched row, through fill-in. Take a voltage row
   whose regulated bus's Vm column is eliminated before the `Q_c` column, its pivot on the
   bus's Q row (its diagonal, dQ/dVm, never switched):

   - PV form: the voltage row has a 1 on the Vm column; eliminating Vm subtracts
     `1 / (dQ/dVm)` times the Q row from it, and the Q row has an entry on `Q_c` (the
     controller's reactive injection, when the controller sits on the regulated bus). The
     voltage row's entry on `Q_c` is then nonzero: the pivot `(v_row, Q_c)` is valid;
   - held form: that entry is the 1 of "Q_c = q", valid as well.

   So with the order "regulated Vm before the `Q_c` columns", the pivot never moves when the
   row switches, and none of these fallbacks happens -- for `NROuter_KLU` as much as for the
   static-pivot solver. KLU already pivots the Vm column on its diagonal (its diagonal
   preference, dQ/dVm being large); what has to be added is the column order. By the fill
   path theorem, the entry `(v_row, Q_c)` fills in iff a path joins the voltage row to `Q_c`
   through columns eliminated earlier. For remote control -- all the cases seen -- the path
   goes regulated bus -> ... -> controller bus, whose Q row carries `Q_c`: the Vm columns
   along that path (in practice the regulated bus and the controller bus, joined by the
   step-up transformer) must all come before `Q_c`. The ratio-group rows are the same
   picture with the bus's Vm column and the transformer's ratio column (the bus's Q row
   carries dQ/dratio when the transformer touches the bus).

   How it could be done:

   - KLU takes a given ordering (`klu_analyze_given`): run AMD as today, then move each
     control column (`Q_c`, ratio) to just after the last Vm column of its path, keeping the
     permutation symmetric so the diagonal preference still applies.
   - Only the Newton system knows the pairs (`vm_cols_`, `q_cols_`, the ratio groups' columns).
     They reach the linear solver through an optional hook found by template detection, like
     `set_rgrowth_threshold`; `BaseAlgo`, hence the plugin ABI, is not touched.
   - Not "all control columns last": with every `Q_c` at the end, the Schur complement on them
     is dense (the grid connects them all), a dense trailing block that would cost more than it
     saves. The move has to be local.

   To check: the fill it adds (the kernel flop count), the 21 off-diagonal ratio-group
   cases, paths longer than one transformer, the share rows of multi-controller groups and
   SVC slopes, and that the iteration counts do not move. If it works, the remaining
   fallbacks (bus masking, odd cases) can still go through the cache of plans.

   **Measured (throw-away prototype, not committed).** The column pairs went from
   `VoltageControl::register_in` and `BranchControl::register_in` to the static-pivot solver
   through a global registry; its `analyze` ran AMD (`amd_order`, default controls), moved
   indices in a stable sort on the position, kept the permutation symmetric, and called
   `klu_analyze_given`. Five ways of building the order were tried on the contingency bench:

   | order | time | fallbacks | update flops per refactorization |
   |---|---|---|---|
   | AMD, as today | 2.06 s | 142 | 571 k |
   | controls after the Vm of every related bus (ratio groups: group bus and both transformer ends) | 1.75-1.82 s | 0 | 1 126 k |
   | **controls after the Vm they need (voltage groups: regulated and controller buses; ratio groups: the group bus only)** | **1.65 s** | **0** | **720 k** |
   | the same, moving only the ratio groups' rows (not their ratio columns) | 1.66 s | 0 | 720 k |
   | the Vm columns moved up, before the first control column (minimal sets) | 3.39 s | 0 | 4 078 k |
   | the same, full sets | 14.20 s | 0 | 31 819 k |

   Every variant removes every fallback: the fill-in argument holds, for the remote voltage
   controls and for the ratio groups, including the 21 off-diagonal cases. What separates
   them is the fill. Moving the control columns down is local; moving the Vm columns up
   drags dense parts of the grid forward and is ruinous.

   With the third order, against the two other algorithms (same build, one core):

   | | `NROuter_KLU` | `NROuter_KLUStaticPivot`, AMD | `NROuter_KLUStaticPivot`, moved order |
   |---|---|---|---|
   | ContingencyAnalysis, 100 N-1 (3 runs) | 2.95 s | 2.06 / 2.11 / 2.12 s | 1.65 / 1.66 / 1.66 s |
   | InjectionSweep, 100 rows | 8.00 s | 5.83 s | 4.17 s |
   | factorizations (contingencies / sweep) | 143 / 442 | 143 / 442 | 1 / 1 |

   About -20 % on the contingencies and -28 % on the sweep on top of the static pivot (-44 %
   and -48 % against KLU). Converged rows, outer and Newton iteration counts are identical to
   KLU's row by row; voltages differ at round-off (per-row max dVm 1.0e-12 / 1.6e-12 /
   2.7e-12 pu and dVa 3.8e-11 / 6.4e-11 / 9.1e-11 deg, median / p90 / max, for the
   contingencies; the same order for the sweep). The distance to the upper bound (~1.45 s) is
   the extra fill, +26 % of update flops.

   **Implemented (static-pivot solver only).** The third order, as a part of
   `KLUStaticPivotLinearSolver` rather than an option: the static-pivot algorithms are opt-in
   already, and the order only differs from AMD's on a system with rows switched by value.
   - Each Newton component declares its constraints (`declare_ordering_constraints`, in the
     protocol of `NRSystem.hpp`): `VoltageControl` for a group whose voltage row has a slot
     on a `Q_c` column (follow: the `Q_c` columns, the voltage and share rows; lead: the Vm
     of the regulated and controller buses), `BranchControl` and `ShuntControl` for each
     group (follow: its rows; lead: the bus' Vm).
   - `NRAlgo` hands them over before each `analyze`, through
     `detail::set_ordering_constraints_impl` (template detection: a linear solver without
     `set_ordering_constraints` builds none, and `BaseAlgo` is not touched).
   - The solver orders with `constrained_amd_order` and `klu_analyze_given`.

   Measured on the bench, same build, one core:

   | | `NROuter_KLU` | `NROuter_KLUStaticPivot` |
   |---|---|---|
   | ContingencyAnalysis, 100 N-1 | 2.99 s | 1.60 s |
   | InjectionSweep, 100 rows | 8.03 s | 3.96 s |
   | factorizations (contingencies / sweep) | 143 / 442 | 1 / 1 |

   On the contingencies `NROuter_KLU` is bit-identical to the build without the constraints.
   Against it, converged
   rows and outer and Newton iteration counts are identical row by row; voltages differ at
   round-off (per-row max dVm up to 2.7e-12 pu, dVa up to 9.1e-11 deg, both benches).

   Tests: the C++ ones check the order itself, and that on a ring with one switched row, AMD
   alone loses the pivot when the row switches while the constrained order keeps it. A Python
   one checks that a contingency analysis on a public grid with transformer voltage control,
   which has one fallback without the order, has none with it. Remote voltage control on the
   same small grid gets no fallback even without the order (AMD happens to order it well
   there), so that case is only covered by the real grid snapshot.

   Not done: the same order in `KLULinearSolver` (its fallbacks are the same).

2. **A faster kernel** (33 % of `Ir`).
   - The "map" variant of Chen et al. (DATE 2016, see `static_pivot_refactorization.md`):
     precompute the address of every update, no dense work vector, no scatter / gather /
     clearing. Estimated ~4 % of `Ir`.
   - Fewer operations: compare the flops of other fill-reducing orderings (COLAMD; AMD with
     BTF on, which the frozen plan would then have to follow).
3. **Small items** (1-4 % each).
   - `_recompute_mask_positions` (1.5 %) rebuilds its tables, with allocations, for every row;
     the positions per bus could be computed once per batch.
   - `memset` / `malloc` on the hot path (~3.5 %), callers not yet identified (hidden behind PLT
     stubs in a release build).
   - `sincos` (3.8 %) to rebuild V from (Va, Vm) at every iteration.

## Side finding

A grid saved with `save_binary` after its outer-loop setup and reloaded does not converge its
base Newton solve (`TooManyIterations`) from the start voltages the original grid converges
from. Some state the outer loops depend on does not survive the round trip. Not investigated.
