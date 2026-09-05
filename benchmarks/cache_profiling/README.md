# Auditing the cached powerflow path

`lightsim2grid` is fastest on the *second* powerflow: the first one builds the
solver-side data (bus labelling, Ybus, the pv/pq split, the slack weights) and
factorizes the Jacobian, and every subsequent solve on the same grid is supposed
to re-stamp only what actually changed. That "cached path" is what a grid2op
episode, a time series or a contingency sweep spends all of its time in, so it
is the one worth auditing.

This directory measures it, in **instructions retired (Ir)**, with callgrind.
Instruction counts are used rather than wall time because they are exactly
reproducible: the same binary on the same grid gives the same number every run,
so a 2% change is a real 2% and not scheduler noise.

## What is measured, and what is not

The driver (`profile_cached_pf.cpp`) turns callgrind's *collection* on only
around the `ac_pf` / `dc_pf` calls themselves. Loading the grid, the warm-up
solves that fill the cache, and the mutation applied between two steps all run
with collection off. So every number below is the cost of one powerflow, on a
grid that has already solved successfully — no start-up, no measurement
scaffolding.

Phases (one callgrind run each, so nothing is ever mixed):

| phase | what happens between two solves | what it answers |
|---|---|---|
| `cold` | — (the very first solve) | what a full build costs |
| `idem` | nothing at all | the floor: what a solve costs when the answer is already known |
| `inj` | every load's P and Q moved ~2% | the ordinary grid2op step |
| `inj_nores` | idem, results not computed | what `compute_results()` costs (by difference with `inj`) |
| `dcac` | idem, but `dc_pf` then `ac_pf` | what `LightSimBackend`'s default `initdc=True` costs |
| `nocache` | idem, `allow_ac_cache_reuse(false)` | what the cache buys today |
| `topo` | one line opened / closed | a solve whose cache the topology retired |
| `inj:every2`, `inj:every3` | idem `inj` | the same solves refactorizing J every 2 / 3 Newton iterations |
| `idem::NRSing_KLU`, `inj::NRSing_KLU` | idem | the same solves without the distributed-slack extension |

## Running it

```bash
# 1. the grids, dumped from pandapower to lightsim2grid's own binary format
#    (so the profile contains no python and no conversion code)
python make_grids.py grids

# 2. the driver, built straight against src/core -- no python, no pybind11
cmake -S . -B build_profile -DCMAKE_BUILD_TYPE=Release
cmake --build build_profile -j

# 3. one callgrind run per (grid, phase)
./run_profile.sh grids build_profile callgrind_out KLU

# 4. the tables
python3 summarize.py callgrind_out
```

`run_profile.sh` needs `valgrind` and `callgrind_annotate`; the build needs
`valgrind/callgrind.h` (Debian/Ubuntu: `valgrind` + `libsuitesparse-dev` for
KLU). Everything is built `-O3 -g`: `-g` changes no code generation, it is only
what lets `callgrind_annotate` attribute instructions to source lines.

### A/B testing one candidate change

`ab_test.sh` builds `src/core` twice -- once as it is, once with a patch script
applied -- and runs every (grid, phase) under both, reporting the instruction
counts side by side **and** whether the two builds return the same answer. The
driver writes a trace of every solve (its iteration count and its full complex
voltage vector, 17 significant digits) which `compare_traces.py` reads back, so
"same answer" is a measurement and not an assumption. `ab_wallclock.sh` is the
same A/B on wall-clock time (best of N runs), because an instruction count knows
nothing about cache misses or the latency of a libm call.

```bash
./ab_test.sh      grids ab_out patches/bus_mismatch_no_temporary.py
python3 compare_traces.py ab_out
./ab_wallclock.sh grids ab_out patches/bus_mismatch_no_temporary.py 9
```

The tree is restored with `git checkout -- src/core` on exit, including on
failure. `patches/` holds candidates that have NOT been taken; a candidate's
script is deleted once the change lands in `src/core`, and its measurements stay
in the matching `results_ab_*.txt`.

## Results

Measured on this branch (`claude/lightsim2grid-cache-profiling-hweo3g`), gcc 13,
`-O3` (no `-march=native`, which is also what the wheels ship), KLU as the linear
solver, `tol = 1e-8`, `max_iter = 10`. Grids come from `pandapower.networks`.

### Instructions per powerflow

| grid | buses | `cold` | `idem` | `inj` | `inj_nores` | `dcac` | `nocache` | `topo` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| case30 | 30 | 708,845 | 109,300 | **177,116** | 156,678 | 266,865 | 421,090 | 463,484 |
| case118 | 118 | 2,029,113 | 429,240 | **864,166** | 782,137 | 1,176,860 | 1,654,745 | 1,584,581 |
| case1354pegase | 1,354 | 28,444,829 | 5,079,312 | **11,580,775** | 10,725,621 | 19,096,753 | 22,065,963 | 21,264,520 |
| case9241pegase | 9,241 | 292,608,869 | 45,922,151 | **111,291,208** | 104,777,902 | 221,384,553 | 196,271,013 | 190,506,621 |

Read the bold column as "the cost of one ordinary grid2op step". The cache is
already doing its job: `inj` is **2.3x to 4.0x** cheaper than the first solve,
and **1.7x to 2.4x** cheaper than the same solve with `allow_ac_cache_reuse(false)`.
Opening a line (`topo`) costs about what running with no cache at all costs,
which is the expected shape: a sparsity-pattern change retires everything.

### Where a cached solve spends its instructions

| grid | the algorithm (NR + linear solver) | of which KLU refactorization | everything `LSGrid` does around it |
|---|---:|---:|---:|
| case30 | 85.2% | 34.9% | 14.8% |
| case118 | 87.9% | 30.9% | 12.1% |
| case1354pegase | 90.9% | 43.0% | 9.1% |
| case9241pegase | 93.0% | 55.3% | 7.0% |

The bigger the grid, the more completely the answer is "KLU, refactorizing the
Jacobian once per Newton iteration". The `LSGrid` side is a shrinking minority,
and more than half of what is left of it is `compute_results()`.

## What could be saved, measured

Every line below is a real measurement: the same phases, re-run against a
modified build, and compared against the table above. None of these changes is
applied in this branch — this is an audit. Percentages are on the ordinary
cached step (`inj`); row 1 is on `idem`, where it is the whole point. Row 6 is
the one that went through the full A/B (both builds, all four grids, all six
phases, answers compared) — see its section.

| # | change | case30 | case118 | case1354 | case9241 |
|---|---|---:|---:|---:|---:|
| 1 | do not re-derive the distributed-slack state each solve (`idem`, upper bound via `NRSing_KLU`) | **-60%** | **-51%** | **-63%** | **-71%** |
| 2 | `initdc=False` on a warm-started step (`dcac` -> `inj`) | -34% | -27% | -39% | **-50%** |
| 3 | refactorize J every 2 Newton iterations instead of every one | +10% | -17% | -20% | **-23%** |
| 4 | do not recompute the element results (`inj` -> `inj_nores`) | -11.5% | -9.5% | -7.4% | -5.9% |
| 5 | cache the voltage-control controller data instead of rebuilding it per solve | -2.2% | -3.4% | — | -1.1% |
| 6 | take `1/\|V\|` from `Vm_` instead of a `hypot` pass over `V_` (**A/B tested**) | -2.1% | -2.5% | -2.2% | -1.6% |
| 7 | make `compute_results` cheaper: no repeated `vector<bool>` (**applied**), no heap temporary / `std::complex` (superseded) | -1.7% | -1.5% | -1.1% | -0.9% |

### Where each of these landed

An audit measures; it does not decide. What was decided, and why:

* **1 (slack state carry-over) — rejected.** Sbus changes between two solves, so
  a slack state carried over is a coupling between what may be two genuinely
  different scenarios. Starting it from scratch is the safe answer, and the cost
  measured here is the price of that safety, not a bug.
* **2 (`initdc`) — rejected.** Warm-starting the AC solve from the previous
  solution was tried before: when the topology changes a lot between two
  consecutive steps, the last powerflow diverges. The DC init is much more
  robust precisely *because* it introduces no coupling between consecutive
  solves. The 2x on case9241pegase is what that robustness costs.
* **3 (refactorization policy) — out of scope.** `EveryN` is a different
  algorithm, not a cheaper way to run this one.
* **5 (voltage-control data) — to do.** Worth it on its own, and more so because
  the data is currently built twice. It needs the cache to carry more than it
  does today (the `AlgoControl` state and what changed on the elements since),
  so it is its own piece of work.
* **6 (`1/|V|` from `Vm_`) — A/B tested and applied**, see below.
* **7 (a cheaper `compute_results`) — the `std::vector<bool>` half is A/B
  tested and applied**; the mismatch half is superseded by reusing the
  algorithm's own `Ybus . V`. See below. Row 4 is the same work skipped
  entirely; row 7 is the same work done for ~15% less.

### 1. The cached path always runs at least one Newton iteration

An identical re-solve — same grid, same injections, starting from the voltage
the previous solve returned — still runs one **full** Newton iteration: fill
`dS/dVm`, `dS/dVa`, assemble J, refactorize it, solve, step, re-evaluate. On
case9241pegase that iteration is the whole 45.9M-Ir floor, 20.5M of it inside
`klu_refactor`.

It is not the powerflow that needs it. `check_solution` puts the physical
mismatch of the returned voltage at 5.7e-12 MW on case118, five orders of
magnitude below the 1e-8 tolerance. The same solve with the single-slack Newton
system converges in **zero** iterations:

```
NR_KLU         first=4 re-solve=1,1
NRSing_KLU     first=4 re-solve=0,0
```

The difference is `MultiSlack::update_state`, which re-initialises the slack
state on every solve from `real(Sbus.sum())` — the imbalance ignoring losses.
That is a fine cold start and a bad warm start: it puts the slack row's residual
at roughly the grid's total losses, above any sensible tolerance, so the
convergence test before the loop can never pass and one iteration is always
spent rediscovering a number the previous solve already had. Carrying
`slack_absorbed_` (and `VoltageControl::q_`, reset to zero the same way) across
solves would let an unchanged grid converge in zero iterations -- the measured
upper bound of that is the whole `idem` column, -60% to -71%. **That is not
being done**, and for a good reason: Sbus changes between two solves, so a
carried-over slack state couples what may be two genuinely different scenarios.
The number above is therefore the price of starting from scratch, not a saving
on the table.

This also matters for the ordinary step, not just the identical one: on `inj`
the iteration count is the same either way (3 iterations), so the saving there
is the difference in *starting point quality*, which is small — but every
scenario built on repeated solves of a barely-changing grid (time series,
N-1 sweeps, a do-nothing grid2op episode) pays the full floor.

### 2. `initdc=True` doubles a warm-started step

`LightSimBackend` runs a DC powerflow before every AC one (`self.initdc = True`,
commented "does not really hurt computation time"). On a cold solve that is
true and useful. On the cached path it is not: the DC solution is a *worse*
starting point than the AC solution of the previous step, so the AC Newton needs
6 iterations instead of 3, and the DC solve itself is not free. On
case9241pegase a step costs 221.4M Ir with the DC init and 111.3M without —
`initdc` doubles it.

Warm-starting from the previous solution instead was tried before and
**rejected**: when the topology changes a lot between two consecutive steps, the
AC solve diverges. The DC init is robust exactly because it carries nothing over
from the previous solve, and the 2x measured here is what that costs.

### 3. Refactorizing every Newton iteration

`RefactorPolicyType::AlwaysRefactor` is the default, and on case9241pegase
`klu_refactor` alone is 55% of a cached solve. The `EveryN` policy is already
implemented; measured on the ordinary step it is worth **-23%** (every 2) on
case9241pegase, at no cost in iteration count (3 either way), and -17% / -20%
on case118 / case1354pegase. It loses on case30, where the iteration the policy
adds costs more than the factorization it skips — the crossover is somewhere
around a hundred buses. `Chord` (factorize once, never again) does not converge
within 10 iterations on any of these grids.

This is a different algorithm rather than a cheaper way to run the current one,
so it is **out of scope** for this audit -- recorded here only because 55% of a
large cached solve sits in `klu_refactor`, which is worth knowing whichever way
the choice goes.

### 4. Result computation

`compute_results()` re-derives every element's P, Q, V, angle and current from
the solved voltages, and it is the single biggest thing `LSGrid` does around the
algorithm — 6.5M Ir per solve on case9241pegase, 84k on case118. It is 5.9% to
11.5% of a cached solve, and it is *entirely* skippable for a caller that only
wants the voltages (`deactivate_result_computation()`); grid2op is not such a
caller, so this is a knob rather than a saving, but a contingency screen or a
feasibility loop should be using it.

Inside it, the two biggest contributors on case9241pegase are the shunts (542k
Ir for 7,327 of them) and the branches; the per-element costs are in line with
each other, so there is no single outlier to fix -- it is simply a lot of
elements, computed in full on every solve.

### 5. The voltage-control data is rebuilt from scratch every solve

`VoltageControl::update_state` runs on every `compute_pf` (it is in the phase
the code calls "cheap, always") and calls
`LSGrid::fill_voltage_control_solver_data`, which builds five `std::set<int>`
(`get_free_vm_slack_solver_buses` alone builds three, and
`get_group_controlled_buses` is called twice), two `std::vector<bool>` of
`nb_bus_solver`, and walks every generator and SVC. On
case9241pegase that is ~1.1M Ir per solve; on case30 it is 3% of the entire
powerflow.

Nothing it computes depends on the injections or on the voltages. It depends on
the bus labelling, the pv/pq split and the generators' regulation configuration
— all three already tracked by `AlgoControl` (`has_pv_changed`, `has_pq_changed`,
`has_v_changed`, `has_dimension_changed`, `need_reset_solver`). Rebuilding it
only when one of those is raised measures at -2.2% (case30), -3.4% (case118)
and -1.1% (case9241pegase). **To be done**: the data is also built twice today,
which makes it worth more than these numbers, but a correct cache of it needs
the `AlgoControl` state and the element changes since the last solve to be
carried alongside it -- its own piece of work.

### 6. `hypot` where the answer is already known (A/B tested)

`fill_internal_variables` starts with `inv_vm_cache_ = V_.array().abs().inverse()`
— a `std::hypot` per bus per Newton iteration. But `Vm_` *is* `|V_|`. Only two
places write `V_`, and both keep it so: `update_state` sets `Vm_ = |V_init|`
next to `V_ = V_init`, and `apply_step` rebuilds `V_ = Vm_ * exp(i.Va_)` and then
calls `fix_negative_vm`, which flips the sign of a negative `Vm_` and turns
`Va_` by half a turn, leaving `V_` alone and still equal to
`Vm_new * exp(i.Va_new)`. The comment a few lines below in `apply_step` makes
exactly this argument for `Vm_` / `Va_` and declines to recompute them; the
reciprocal pass was not given the same treatment.

`inv_vm_cache_ = Vm_.array().inverse()`, A/B tested with `ab_test.sh` over all
four grids and six phases before it was applied:

| grid | `idem` | `inj` | `dcac` | `topo` | `nocache` | `cold` |
|---|---:|---:|---:|---:|---:|---:|
| case30 | -1.71% | **-2.11%** | -2.10% | -1.21% | -0.89% | -0.93% |
| case118 | -1.73% | **-2.53%** | -2.52% | -1.40% | -1.32% | -1.30% |
| case1354pegase | -1.69% | **-2.22%** | -2.25% | -1.21% | -1.17% | -1.39% |
| case9241pegase | -1.27% | **-1.57%** | -1.59% | -0.92% | -0.89% | -1.12% |

The shape is what it should be: the saving is largest on the phases that are
mostly Newton iterations (`inj`, `dcac`) and smallest on the ones a rebuild
dominates (`cold`, `topo`, `nocache`), because the `hypot` pass is paid once per
iteration.

An instruction count knows nothing about the latency of a libm call, so
`ab_wallclock.sh` ran the same A/B on the clock (best of 9, no valgrind):

| grid | `idem` A -> B | | `inj` A -> B | |
|---|---:|---:|---:|---:|
| case30 | 0.00991 -> 0.00983 ms | -0.84% | 0.01573 -> 0.01553 ms | -1.31% |
| case118 | 0.0357 -> 0.0350 ms | -1.78% | 0.0721 -> 0.0697 ms | **-3.24%** |
| case1354pegase | 0.631 -> 0.614 ms | -2.79% | 1.509 -> 1.469 ms | -2.64% |
| case9241pegase | 5.704 -> 5.581 ms | -2.15% | 14.556 -> 14.072 ms | **-3.33%** |

Time moves slightly *more* than instructions on the large grids, which is what a
removed libm call should do.

**Does it return the same answer?** Every one of the 24 (grid, phase) pairs
converges in **exactly the same number of iterations**, and the returned
voltages agree to:

| grid | max abs. difference | relative | vs. the 1e-8 tolerance |
|---|---:|---:|---|
| case30 | 4.4e-15 pu | 20 ulp | 6 orders of magnitude below |
| case118 | 5.8e-15 pu | 26 ulp | 6 orders below |
| case1354pegase | 2.1e-13 pu | 900 ulp | 5 orders below |
| case9241pegase | 1.4e-12 pu | 6,400 ulp | 4 orders below |

Not bit for bit, and it cannot be: `|cos(a) + i.sin(a)|` is 1 only to within a
rounding, so `|V_|` and `Vm_` differ in the last bits, and a Newton iteration
amplifies that into the iterate. Two things bound what it can do. First, the
value feeds *only* the `dS_dVm` block of the Jacobian — it never enters the
mismatch, which is computed from `V` and `Sbus` and is what the convergence test
reads. A rounding there perturbs the search direction, never the criterion the
answer is accepted on. Second, the differences above are four to six orders of
magnitude below the tolerance both solutions were accepted at: they are not two
different answers, they are the same answer with the Newton having stopped at a
slightly different point inside its own tolerance ball. The same argument
already justifies the sibling shortcut in `apply_step`, whose comment records a
1.9e-16 relative move for the same reason.

Recommendation: take it. One line, no behaviour to reason about beyond the
identity `Vm_ == |V_|`, ~2% off the ordinary cached step and ~1.5% off a large
one, iteration counts untouched.

In the same family, and not measured here: `_reconstruct_V_into` computes
`Va.cos()` and `Va.sin()` as two separate libm passes (`__sincos_fma` shows up
twice in the profile, ~2.5% together on case9241pegase); one `sincos` pass would
halve that.

### 7. `compute_results`, done for ~15% less (A/B tested)

Row 4 measures switching result computation off. Most callers cannot: grid2op
needs every element's P, Q, V, angle and current. So the question is what the
work actually costs, and the profile (the `inj` run minus the `inj_nores` one,
which is exactly `compute_results` and nothing else) splits its 6.51M Ir on
case9241pegase as:

| | Ir/solve | share |
|---|---:|---:|
| `TwoSidesContainer_rxh_A::compute_results_tsc_rxha_no_amps` (branch flows) | 2.87M | 44% |
| — of which `std::vector<bool>` reads and out-of-line bus lookups | 0.98M | 15% |
| `LSGrid::compute_results` itself (the per-bus mismatch) | 1.34M | 21% |
| `OneSideContainer::get_bus_internal` (`std::vector<bool>`) | 0.58M | 9% |
| `GenericContainer::_get_amps` | 0.32M | 5% |
| `GenericContainer::v_kv_theta_from_vpu` | 0.29M | 4% |
| `ShuntContainer::_compute_results` | 0.28M | 4% |

`_get_amps` and `v_kv_theta_from_vpu` are already at ~5 and ~22 instructions per
element and have nothing left in them. Two things do:

**The per-bus mismatch** (`patches/bus_mismatch_no_temporary.py`) -- *superseded,
see below*.
`mismatch = V.array() * (mat * V).conjugate().array() - inj.array()` pays three
things it need not. `mat * V` inside a coefficient-wise expression is a
sparse-times-dense product Eigen can only evaluate into a **heap temporary** --
one allocation and one full-vector zero-fill per solve; `NRSystem::_residual_into`
hit the same wall and answered it with a buffer that outlives the call. The
coefficient-wise multiply then goes through `std::complex::operator*`, whose
NaN-recovery branch after every product is the cost the dS pass and the
branch-flow loop already spell out on real and imaginary parts to avoid. And
`mismatch` is a full complex vector nothing reads: only its two halves are
scaled into the real vectors anyone downstream sees. One pass writing those two
directly removes all three.

**The `std::vector<bool>`** -- *applied*.
`status1[el_id]` and `status2[el_id]` are read five times each per branch, and
each read is a word offset, a shift and a mask rather than a load. On top of
that `get_bus_side_1_internal` re-reads the same status bit to decide whether to
hand back `_deactivated_bus_id`, inside a branch that has just established the
side is connected -- and it does not inline, which is the 0.58M
`get_bus_internal` line above. Read each bit once, and take the bus id straight
from the side's own `bus_id_` where the status is already known.

Measured on its own, it is most of the prize:

| grid | `idem` | `inj` | Ir saved per solve | vs. `compute_results` |
|---|---:|---:|---:|---:|
| case30 | -2.35% | -1.46% | 2,527 | **-12.4%** |
| case118 | -2.42% | -1.21% | 10,195 | **-12.4%** |
| case1354pegase | -2.31% | -1.02% | 115,528 | **-13.5%** |
| case9241pegase | -2.05% | -0.85% | 931,292 | **-14.3%** |

-- 91% of what the two together are worth on case9241pegase, for eleven
bit-vector reads per branch and one call that would not inline.

Measured together, on top of the `1/|V|` change:

| grid | `idem` | `inj` | `dcac` | `topo` | Ir saved per solve | vs. `compute_results` |
|---|---:|---:|---:|---:|---:|---:|
| case30 | -2.77% | -1.72% | -1.16% | -0.65% | 2,981 | **-14.6%** |
| case118 | -3.03% | -1.52% | -1.12% | -0.88% | 12,760 | **-15.6%** |
| case1354pegase | -2.58% | -1.14% | -0.69% | -0.62% | 128,631 | **-15.0%** |
| case9241pegase | -2.25% | -0.93% | -0.47% | -0.54% | 1,020,293 | **-15.7%** |

The absolute saving is the same number in every phase of a given grid --
1,020,293 instructions on case9241pegase, to the instruction, in both `idem` and
`inj` -- which is what removing a fixed per-solve cost looks like, and a check
that the harness measures what it claims. On the clock: -2.5% / -2.2% on `idem`
for case118 / case30, -0.5% to -1.6% elsewhere (best of 9).

**Does it return the same answer?** Yes, and here that is not a judgement call:
the returned voltages are **identical, 0 ulp**, in all 16 (grid, phase) pairs,
with the same iteration counts -- for the pair, and for the `std::vector<bool>`
half on its own. Unlike row 6, which moved a rounding, both of these are exact
rearrangements: the same products in the same grouping and the same order, with
an allocation and a re-read taken out around them. The C++ suite passes
unchanged (229 test cases, 590,862 assertions).

**Where the mismatch half goes next.** It is worth only ~89k Ir per solve on
case9241pegase once the `std::vector<bool>` half has landed, because most of the
1.34M it profiled at is the `Ybus . V` product itself -- real work, not overhead.
But that product is not work `compute_results` has to do at all: the algorithm
computed `Ybus . V` at the accepted voltage on its last residual evaluation, and
`NRSystem::_residual_into` already keeps it in a buffer that outlives the call
(`ybus_v_`, next to the `mis_bus_` that `BaseAlgo::get_bus_mismatch()` exposes).
Reading it back instead of recomputing it is worth the whole product, not just
the temporary around it. The catch to work out is that `mis_bus_` itself is NOT
the raw mismatch `compute_results` wants -- the components fold their own
injections into it (the distributed slack's share, the hvdc droop flows, a
voltage controller's reactive output), so at convergence it is ~0 where
`compute_results` needs the slack's absorbed power and the generators' reactive
output. `ybus_v_` is the reusable half; `mis_bus_` is not, without backing those
adjustments out.

### What is NOT worth attacking

* **The build side of the cache.** `_build_into_cache` costs 1.09M Ir per solve
  on case9241pegase — 1% of the solve — and most of that is the Sbus refill,
  which genuinely has to happen when injections change. The change-flag
  machinery is doing its job.
* **`_get_results_back_to_orig_nodes`.** 176k Ir on case9241pegase (0.2%),
  despite allocating a full `total_bus()` vector per solve.
* **The bus-labelling and pv/pq work.** Zero on the cached path: the flags keep
  `init_converter_bus_id` / `fillpv_pq` from running at all.

## Reproducing the numbers in this file

`results_summary.txt` next to this README is the raw output of `summarize.py`
for the run the tables above were written from (44 callgrind runs: 4 grids x 11
phases). The two candidate fixes in rows 5 and 6 of the savings table were
measured on throwaway builds of `src/core` and are not applied anywhere in the
tree; the patches are three lines each and are described in full in their
sections.
