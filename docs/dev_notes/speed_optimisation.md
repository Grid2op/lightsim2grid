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
- What is left is the outer loops' value-level switches of a voltage-control row between its
  two forms, "Vm(reg) = target" (entry on the Vm column) and "Q_c = q" (entry on the Q_c
  column), when ReactiveLimits / VoltageMonitoring switch a bus PV <-> PQ or hold a
  controller. Each form zeroes the other's entry, and such a row has no other nonzero: its
  pivot is on one of the two, so one of the forms always puts it at zero. This is structural:
  **no single pivot order can serve both forms of such a row.** Not yet confirmed loop by loop.
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
