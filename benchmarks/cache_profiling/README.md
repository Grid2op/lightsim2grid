# Auditing the cached powerflow path

`lightsim2grid` is fastest on the *second* powerflow: the first one builds the
solver-side data (bus labelling, Ybus, the pv/pq split, the slack weights) and
factorizes the Jacobian, and every subsequent solve on the same grid re-stamps
only what actually changed. That cached path is where a grid2op episode, a time
series or a contingency sweep spends all of its time, so it is the one worth
measuring.

This directory measures it in **instructions retired**, with callgrind.
Instruction counts rather than wall time because they are exactly reproducible:
the same binary on the same grid gives the same number every run, so a 2% change
is a real 2% and not scheduler noise.

## What is measured, and what is not

The driver (`profile_cached_pf.cpp`) turns callgrind's *collection* on only
around the `ac_pf` / `dc_pf` calls themselves. Loading the grid, the warm-up
solves that fill the cache, and the mutation applied between two steps all run
with collection off. Every number below is therefore the cost of one powerflow
on a grid that has already solved successfully.

Phases, one callgrind run each:

| phase | between two solves | answers |
|---|---|---|
| `cold` | — (the very first solve) | what a full build costs |
| `idem` | nothing at all | the floor: a solve whose answer is already known |
| `inj` | every load's P and Q moved ~2% | the ordinary grid2op step |
| `inj_nores` | idem, results not computed | what `compute_results()` costs, by difference |
| `dcac` | idem, but `dc_pf` then `ac_pf` | what `LightSimBackend`'s default `initdc=True` costs |
| `nocache` | idem, `allow_ac_cache_reuse(false)` | what the cache buys |
| `topo` | one line opened / closed | a solve whose cache the topology retired |
| `inj:everyN` | idem `inj` | the same solves refactorizing J every N iterations |
| `inj::NRSing_KLU` | idem | the same solves without the distributed-slack extension |

## Running it

```bash
# 1. the grids, dumped from pandapower to lightsim2grid's binary format, so the
#    profile contains no python and no conversion code
python make_grids.py grids

# 2. the driver, built straight against src/core -- no python, no pybind11
cmake -S . -B build_profile -DCMAKE_BUILD_TYPE=Release
cmake --build build_profile -j

# 3. one callgrind run per (grid, phase), then the tables
./run_profile.sh grids build_profile callgrind_out KLU
python3 summarize.py callgrind_out
```

Needs `valgrind` + `callgrind_annotate`, and `valgrind/callgrind.h` to build.
Everything is `-O3 -g`: `-g` changes no code generation, it is only what lets
`callgrind_annotate` attribute instructions to source lines.

### A/B testing one candidate change

`ab_test.sh` builds `src/core` twice -- once as it is, once with a patch script
applied -- runs every (grid, phase) under both, and reports the instruction
counts side by side **and whether the two builds return the same answer**. The
driver writes a trace of every solve (its iteration count and its full complex
voltage vector, 17 significant digits) which `compare_traces.py` reads back, so
"no behaviour change" is a measurement and not an assumption. `ab_wallclock.sh`
is the same A/B on the clock (best of N), because an instruction count knows
nothing about the latency of a libm call.

```bash
./ab_test.sh      grids ab_out my_patch.py
python3 compare_traces.py ab_out
./ab_wallclock.sh grids ab_out my_patch.py 9
```

A patch script is any python that edits `src/core` in place -- usually a couple
of `str.replace` calls guarded by an assertion that the anchor was found. The
tree is restored with `git checkout -- src/core` on exit, including on failure.

## Results

Measured with gcc 13, `-O3` (no `-march=native`, as the wheels ship), KLU,
`tol = 1e-8`, `max_iter = 10`, on `pandapower.networks` grids.

### Instructions per powerflow

| grid | buses | `cold` | `idem` | `inj` | `inj_nores` | `dcac` | `nocache` | `topo` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| case30 | 30 | 708,845 | 109,300 | **177,116** | 156,678 | 266,865 | 421,090 | 463,484 |
| case118 | 118 | 2,029,113 | 429,240 | **864,166** | 782,137 | 1,176,860 | 1,654,745 | 1,584,581 |
| case1354pegase | 1,354 | 28,444,829 | 5,079,312 | **11,580,775** | 10,725,621 | 19,096,753 | 22,065,963 | 21,264,520 |
| case9241pegase | 9,241 | 292,608,869 | 45,922,151 | **111,291,208** | 104,777,902 | 221,384,553 | 196,271,013 | 190,506,621 |

The bold column is one ordinary grid2op step. The cache is doing its job: it is
**2.3x to 4.0x** cheaper than the first solve and **1.7x to 2.4x** cheaper than
the same solve with `allow_ac_cache_reuse(false)`. Opening a line costs about
what running with no cache costs, which is the expected shape.

### Where a cached solve spends them

| grid | algorithm (NR + linear solver) | of which KLU refactorization | everything `LSGrid` does around it |
|---|---:|---:|---:|
| case30 | 85.2% | 34.9% | 14.8% |
| case118 | 87.9% | 30.9% | 12.1% |
| case1354pegase | 90.9% | 43.0% | 9.1% |
| case9241pegase | 93.0% | 55.3% | 7.0% |

The bigger the grid, the more completely the answer is "KLU, refactorizing the
Jacobian once per Newton iteration". The `LSGrid` side is a shrinking minority.

## What the audit found

Acted on (each A/B'd, answers compared; see the changelog for the details):

| change | case30 | case118 | case1354 | case9241 |
|---|---:|---:|---:|---:|
| `1/\|V\|` from `Vm_` instead of a `hypot` pass over `V_` | -2.1% | -2.5% | -2.2% | -1.6% |
| read each branch status bit once, not five times | -1.5% | -1.2% | -1.0% | -0.9% |
| take the per-bus mismatch off the algorithm | -2.1% | -2.1% | — | -1.4% |

Measured and **declined**, recorded here so they are not re-proposed:

* **Carrying the distributed-slack state across solves.** The cached path always
  runs at least one full Newton iteration, because `MultiSlack::update_state`
  re-derives the slack state from `real(Sbus.sum())` -- the imbalance ignoring
  losses -- rather than reusing the previous solve's. The same solve with
  `NRSing_KLU` converges in *zero* iterations and costs 60-71% less. Declined:
  Sbus changes between two solves, so a carried-over slack state couples what may
  be two genuinely different scenarios. That column is the price of that safety.
* **`initdc=False` on a warm-started step.** `LightSimBackend` runs a DC solve
  before every AC one. On a cached step the DC solution is a *worse* starting
  point than the previous AC one, so the Newton needs 6 iterations instead of 3:
  221.4M instructions against 111.3M on case9241pegase. Declined: warm-starting
  from the previous solution was tried before and diverges when the topology
  changes a lot between two steps. The DC init is robust precisely because it
  carries nothing over.
* **Refactorizing J every 2 Newton iterations** (`RefactorPolicyType::EveryN`,
  already implemented) is worth -23% on case9241pegase at the same iteration
  count, and -17%/-20% on case118/case1354pegase; it loses on case30, where the
  iteration it adds costs more than the factorization it skips. Out of scope: a
  different algorithm, not a cheaper way to run this one.

Measured and **not worth attacking**: the build side of the cache
(`_build_into_cache` is ~1% of a case9241pegase solve, most of it the Sbus refill
that genuinely has to happen), `_get_results_back_to_orig_nodes` (0.2%), and the
bus-labelling and pv/pq work, which the change flags keep from running at all.
