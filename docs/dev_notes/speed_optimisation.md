# Speed optimisation of the Newton-Raphson with outer loops: what is left

**Design note, not documentation.** It lists the performance work still open on a batch run
with OpenLoadFlow's outer loops (`NROuter_*`). It is not part of the built documentation
(`docs/*.rst`).

What is done is documented with the code, not here:

- the static-pivot refactorization (`KLUStaticPivotLinearSolver`, the `*_KLUStaticPivot`
  algorithms): `docs/solvers.rst`, section "static-pivot-refactorization", and the header
  of `KLUStaticPivotSolver.hpp`, with its references;
- the constrained ordering (`constrained_ordering`), which removed the fallback
  factorizations the outer loops caused when they switch a row by value: `docs/solvers.rst`,
  section "constrained-ordering", `OrderingConstraints.hpp`, and the
  `declare_ordering_constraints` step of the component protocol in `NRSystem.hpp`;
- the grid settings `save_binary` used to drop (binary format 22): the `[1.3.0]` changelog.

Companion note: `outer_loops_fixed_sparsity.md` -- why the outer loops keep one Jacobian
pattern per solve and change values only.

## The bench

So that measurements compare from one change to the next:

- one large real grid snapshot, `NROuter_*` with all seven OpenLoadFlow outer loops,
  `tol = 1e-8`, `max_iter = 50`, a DC start, the process pinned to one core (`taskset -c 0`);
- **ContingencyAnalysis**: 100 line N-1 drawn with a fixed seed among the lines closed at both
  ends, `init_from_n_powerflow = True`, `handle_disconnected_grid = True`; and the full N-1
  (`add_all_n1`) with the same settings;
- **InjectionSweep**: 100 rows, every generator's target P scaled by `1 + 0.001 * (k - 50)` and
  every load's by slightly more (a small imbalance for DistributedSlack),
  `init_from_n_powerflow = True`;
- compared: the wall time of `compute()`, the converged rows, the outer and Newton iteration
  counts, the linear-solver counters (`LinearSolverStats`), and the per-row maximum of
  dVm / dVa against the reference algorithm.

The bench scripts read the snapshot and so live outside the repository.

Where things stand on it (wall time of `compute()`, one core):

| | `NROuter_KLU` | `NROuter_KLUStaticPivot` (constrained order) |
|---|---|---|
| ContingencyAnalysis, 100 N-1 | 2.99 s | 1.60 s |
| InjectionSweep, 100 rows | 8.03 s | 3.96 s |
| full N-1 | 301 s | 162 s |

## Open items

1. **Profile again.** The last profile (callgrind, `compute()` only) predates the constrained
   ordering. The fallback factorizations it showed are gone, so the refactorization kernel
   and the Newton system (filling J, the residual, `Ybus . V`, `sincos`, the mask positions)
   now carry most of the time. The items below are ranked on the old profile and should be
   re-ranked first.
2. **The fallbacks left on the full N-1.** A handful of rows still go through a fresh
   factorization with the constrained order. Not classified yet: a contingency stranding a
   bus (`handle_disconnected_grid` masking), or a switched row the constraints do not cover
   (a voltage-control group whose path from the regulated bus to the controller crosses more
   than one transformer, the share rows of multi-controller groups, an SVC slope).
3. **A faster refactorization kernel.**
   - The "map" variant of Chen et al. (DATE 2016, reference in `KLUStaticPivotSolver.hpp`):
     precompute the address of every update, no dense work vector, no scatter / gather /
     clearing.
   - Fewer operations: compare the flops of other fill-reducing orderings (COLAMD; AMD with
     BTF on, which the frozen plan would then have to follow) under the same constraints.
4. **Small items.**
   - `_recompute_mask_positions` (`NRSystem.hpp`) rebuilds its tables, with allocations, for
     every batch row whose mask changed; the positions per bus could be computed once per
     batch.
   - `memset` / `malloc` on the hot path, callers not yet identified (hidden behind PLT
     stubs in a release build).
   - `sincos` to rebuild V from (Va, Vm) at every iteration.
5. **The constrained order for KLU.** On `NROuter_KLU` it removes the same fallbacks but
   gains little, because the fill it adds costs a `klu_refactor` nearly what the fallbacks
   did; hence opt-in. Not profiled further: whether a lighter set of constraints keeps most
   of the benefit at less fill.
