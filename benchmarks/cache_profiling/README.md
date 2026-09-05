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
applied in this branch — this is an audit.

| # | change | case30 | case118 | case1354 | case9241 |
|---|---|---:|---:|---:|---:|
| 1 | do not re-derive the distributed-slack state each solve (`idem`, upper bound via `NRSing_KLU`) | **-60%** | **-51%** | **-63%** | **-71%** |
| 2 | `initdc=False` on a warm-started step (`dcac` -> `inj`) | -34% | -27% | -39% | **-50%** |
| 3 | refactorize J every 2 Newton iterations instead of every one | +10% | -17% | -20% | **-23%** |
| 4 | do not recompute the element results (`inj` -> `inj_nores`) | -11.5% | -9.5% | -7.4% | -5.9% |
| 5 | cache the voltage-control controller data instead of rebuilding it per solve | -2.2% | -3.4% | — | -1.1% |
| 6 | take `1/\|V\|` from `Vm_` instead of a `hypot` pass over `V_` | — | -2.5% | — | -1.6% |

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
solves when the cache is reusable would let an unchanged grid converge in zero
iterations. The measured upper bound of that is the whole `idem` column: -60% to
-71%.

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

The DC init earns its place on the first step of an episode and after a topology
change. Keeping it there and warm-starting from the previous solution otherwise
is worth half the CPU of a time series on a large grid.

### 3. Refactorizing every Newton iteration

`RefactorPolicyType::AlwaysRefactor` is the default, and on case9241pegase
`klu_refactor` alone is 55% of a cached solve. The `EveryN` policy is already
implemented; measured on the ordinary step it is worth **-23%** (every 2) on
case9241pegase, at no cost in iteration count (3 either way), and -17% / -20%
on case118 / case1354pegase. It loses on case30, where the iteration the policy
adds costs more than the factorization it skips — the crossover is somewhere
around a hundred buses. `Chord` (factorize once, never again) does not converge
within 10 iterations on any of these grids.

This is a numerical-behaviour change, not a free win: it is the one entry in
this table that trades iterations for factorizations. But the trade is currently
not made at all, and on large grids it is clearly favourable.

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
and -1.1% (case9241pegase).

### 6. `hypot` where the answer is already known

`fill_internal_variables` starts with `inv_vm_cache_ = V_.array().abs().inverse()`
— a `std::hypot` per bus per Newton iteration. But `Vm_` *is* `|V_|`:
`apply_step` rebuilds `V_` as `Vm_ * exp(i.Va_)` and then makes `Vm_`
non-negative, and `update_state` sets `Vm_ = |V_init|`. The comment a few lines
below in `apply_step` makes exactly this argument for `Vm_`/`Va_` and declines to
recompute them; the reciprocal pass was not given the same treatment.
`inv_vm_cache_ = Vm_.array().inverse()` measures at -2.5% (case118) and -1.6%
(case9241pegase), with the iteration count unchanged.

In the same family, and not measured here: `_reconstruct_V_into` computes
`Va.cos()` and `Va.sin()` as two separate libm passes (`__sincos_fma` shows up
twice in the profile, ~2.5% together on case9241pegase); one `sincos` pass would
halve that.

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
