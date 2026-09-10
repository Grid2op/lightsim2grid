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
#    profile contains no python and no conversion code. Two families: the plain
#    pandapower cases, and a `_fancy` variant of the bigger ones carrying remote
#    voltage-control groups and voltage-mode SVCs (see `make_fancy` there).
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
| derive the voltage-control plan once per solve, not four times | -4.9% | -8.3% | -3.2% | -2.2% |

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

### The `_fancy` grids

`make_grids.py` also writes a `_fancy` variant of case118, case1354pegase and
case9241pegase: a few voltage-mode SVCs, plus a few *pairs* of generators
re-pointed at a common neighbouring load bus (a control **group**, solved by the
bordered VoltageControl block). Every setpoint is a magnitude the grid already
holds, so each controller is a no-op at the solution: the SVCs take the plain
case's solved magnitudes, and the groups are pointed at the magnitudes of an AC
solve of the grid *with* the SVCs -- every exotic element in place except the
remote control. Each candidate is kept only if the grid still converges with it
-- the same chunk-and-verify workaround `benchmarks/make_exotic_grid.cpp` uses,
for the same reason (see the remote-voltage-control entry in the changelog's TODO).
(Up to `make_grids.py`'s revision of September 2026 both took the plain case's
magnitudes; the tables below were measured on those grids.)

They exist because the plain pandapower cases have **no** remote voltage control
and no SVC at all: the controller list is empty on every one of them, so nothing
that concerns the bordered block shows up. Two things they measured:

| grid | groups | SVCs | `idem` | `inj` | vs. the plain case (`idem`) |
|---|---:|---:|---:|---:|---:|
| case118_fancy | 8 | 4 | 385,022 | 1,173,068 | +16% |
| case1354pegase_fancy | 20 | 10 | 10,293,316 | 15,583,429 | +138% |
| case9241pegase_fancy | 40 | 30 | 318,321,358 | 596,998,625 | +686% |

That last column is not the plan, and not the bordered rows either: on
case9241pegase_fancy **86% of a cached solve is `klu_refactor`** (against 55% on
the plain case). Forty two-member groups add 80 reactive-injection columns, each
coupling a generator's bus to a regulated bus a branch away, and KLU's ordering
pays for the fill-in that creates. Worth knowing before remote control is enabled
at scale on a large grid; out of scope here.

### Deriving the voltage-control plan once per solve

The "fancy" voltage controllers are described by three derived sets -- which buses
a control GROUP regulates, which slack buses keep a free Vm unknown, and the
controller list itself -- and each of them used to be re-derived where it happened
to be needed: by `fillpv_pq`, by `Base::update_state`, and twice by
`VoltageControl::update_state` (which re-ran the free-Vm slack pass of its own).
Four walks of every generator of the grid per powerflow, each building `std::set`s
as it went. They are now one object (`VoltageControlPlan`), built once into the
solver-side cache and read from there, and kept across a solve that changed none
of its inputs -- which an ordinary grid2op step does not (moving a load's P and Q
raises `need_recompute_sbus` and nothing else).

A/B against the tree as it was, KLU, answers compared bit for bit:

| grid | `inj` (an ordinary step) | `idem` (the floor) | `topo` (a line toggled) |
|---|---:|---:|---:|
| case30 | **-4.9%** | -7.9% | -2.3% |
| case118 | **-8.3%** | -16.9% | -4.9% |
| case118_fancy | **-7.0%** | -18.6% | -4.0% |
| case1354pegase | **-3.2%** | -7.3% | -1.9% |
| case1354pegase_fancy | **-2.6%** | -3.9% | -1.7% |
| case9241pegase | **-2.2%** | -5.4% | -1.4% |
| case9241pegase_fancy | **-0.4%** | -0.8% | -0.5% |

The two phases save the *same* number of instructions, to the last one (8,040 on
case30 up to 2.3M on case9241pegase): what is removed is a fixed per-solve cost,
paid before the Newton loop starts. `inj` reads smaller only because it is a
1.5x-2.7x bigger solve. The proportion falls with grid size because the work
removed is O(generators) while the solve it is measured against is dominated by
KLU -- and it falls furthest on the fancy pegase case for the reason the table
above gives.

About 1.9k of the per-solve figure is a second, unrelated find: `AlgorithmSelector`
took its `error_msg` by `const std::string &` and every call site passes a literal
longer than libstdc++'s small-string buffer, so `get_V` / `get_Va` / `get_Vm` /
`compute_pf` / `tell_solver_control` each did a malloc and a free per solve to build
a string only the error path reads. `const char *` now.

Part of the `topo` column is a third find, and it is the one the benchmark was used
to *settle* rather than to report. `need_recompute_pv_pq()` listed
`ybus_change_sparsity_pattern_` -- "the bus labelling may have moved" -- among the
reasons to rebuild the pv/pq split. It is raised only by branch-side mutations
(reconnecting a line or a trafo, moving one of its ends), and every one of those
goes through `GenericContainer::_apply_and_track_buses`, which raises
`change_dimension_` exactly when the mutation empties or fills a bus -- which is
exactly when the labelling moves. So the term was subsumed. The argument was checked
against the grid before it was believed: `src/tests/test_cache_reuse.cpp` reaches a
state where this flag is the ONLY one of the six raised, solves warm and cold and
compares, and the predicate was then built both ways. Dropping the term leaves those
cases green; dropping `slack_participate_changed_` the same way makes them fail, so
that one stays. What the drop is worth, on its own, on a solve whose cache a
topology change retired:

| grid | `topo`, term kept | `topo`, term dropped |
|---|---:|---:|
| case30 | 440,089 | 436,057 (**-0.9%**) |
| case118 | 1,470,391 | 1,444,534 (**-1.8%**) |
| case118_fancy | 2,443,695 | 2,409,714 (**-1.4%**) |
| case1354pegase | 20,163,544 | 20,013,614 (**-0.7%**) |
| case1354pegase_fancy | 28,649,453 | 28,470,141 (**-0.6%**) |
| case9241pegase | 182,113,949 | 181,393,356 (**-0.4%**) |
| case9241pegase_fancy | 536,874,742 | 536,025,060 (**-0.2%**) |

(the `topo` phase toggles one line, so only every other solve raises the flag at
all; the answers are bit-identical on all seven grids). A note for whoever measures
next: the driver links `liblightsim2grid_core.so` with a RUNPATH into its own build
directory, so two binaries copied out of the SAME build tree load whatever library
that tree holds at run time, not the one they were built with. A/B either from two
separate build directories, or with `LD_LIBRARY_PATH` -- which wins over RUNPATH --
pointing at a saved copy of each library.

### Algorithms that cannot do voltage control

The plan is only ever read by the `Base` and `VoltageControl` components of
NRSystem, i.e. by the Newton-Raphson algorithms. Deriving it once per powerflow --
into the cache, rather than lazily in the extension that wanted it -- therefore made
the fast-decoupled and Gauss-Seidel algorithms pay for two container walks nobody
reads: **+2.20% / +0.57%** of a rebuild on case118 / case9241pegase (`nocache`,
FDPF_XB_KLU), and nothing on an ordinary step, where the plan is reused anyway.

They now build no plan at all -- and, more to the point, they no longer take a
group-regulated bus out of PV, which was a *wrong answer* rather than a slow one
(0.36 pu on case118_fancy). Against the same pre-change baseline:

| grid | phase | before | after |
|---|---|---:|---:|
| case118 | `inj` | 1,017,431 | 1,017,063 (**-0.04%**) |
| case118 | `nocache` | 1,851,440 | 1,849,757 (**-0.09%**) |
| case9241pegase | `inj` | 135,854,646 | 135,886,270 (+0.02%) |
| case9241pegase | `nocache` | 218,035,000 | 217,924,056 (**-0.05%**) |

Profiling those at all needed two fixes of its own: `change_algorithm` by NAME did
not call `init_fdpf_coeffs()` (so an FDPF solver selected the way this driver selects
it threw on the first solve), and the driver's iteration budget was a hard-coded 10
where fast-decoupled needs ~50. Both are fixed; `max_iter` is now derived from
`BaseAlgo::is_fdpf`.

Measured and **not worth attacking**: the build side of the cache
(`_build_into_cache` is ~1% of a case9241pegase solve, most of it the Sbus refill
that genuinely has to happen), `_get_results_back_to_orig_nodes` (0.2%), and the
bus-labelling and pv/pq work, which the change flags keep from running at all.

## The batch algorithms

`profile_batch.cpp` is the same audit for `TimeSeries` and `ContingencyAnalysis`:
one solve per row on a Jacobian whose sparsity is fixed for the whole batch, so a
row should cost what a cached solve costs and not much more. It collects only the
call under audit -- `compute()`, or the two flow computations -- and every number
is **per row**. `run_profile_batch.sh` runs every compute phase twice, over N rows
and over one, so the fixed cost of a `compute()` (the rebuild of the solver input,
the "n" solve with its symbolic analysis) is read apart from the marginal cost of a
row. `ab_test.sh` and `ab_wallclock.sh` drive it with `DRIVER=batch`; the trace
they compare holds every row's voltages (or flows) with 17 significant digits.

| phase | what it measures |
|---|---|
| `ts_ac` / `ts_dc` | `TimeSeries::compute()`, loads moved ~2% per row |
| `ts_flows` | `compute_flows()` + `compute_power_flows()` after it |
| `ca_ac` / `ca_dc` | `ContingencyAnalysis::compute()` over the first N-1 contingencies |
| `ca_ac_mask` / `ca_dc_mask` | idem with `handle_disconnected_grid` |
| `ca_flows` | the flows of a `ca_ac` run |
| `ca_construct` | a fresh `ContingencyAnalysis` from a solved grid, 8 contingencies, `compute()`: what a grid2op loop pays per step |

### Baseline: instructions per row, before any batch-side change

KLU, `-O3`, `tol = 1e-8`; 200 / 50 / 20 rows by grid size. The marginal figures
are the N-row run minus the 1-row run, over N-1.

| grid | `ts_ac` /row | `ts_dc` /row | `ca_ac` /row | `ca_dc` /row | `ca_ac_mask` /row | `ca_dc_mask` /row | `ts_flows` /row | `ca_construct` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| case30 | 176,971 | 10,910 | 50,162 | 3,548 | 54,841 | 5,724 | 13,628 | 2,509,541 |
| case118 | 898,340 | 35,939 | 851,636 | 55,189 | 913,065 | 91,130 | 61,355 | 6,853,738 |
| case1354pegase | 11,902,459 | 324,072 | 7,773,483 | 538,826 | 12,109,772 | 1,140,741 | 679,191 | 41,473,166 |
| case9241pegase | 112,346,159 | 2,439,729 | 124,978,834 | 6,349,065 | 126,462,914 | 9,285,519 | 6,300,964 | 1,108,958,718 |

What a `compute()` pays before its first row (the 1-row run):

| grid | `ts_ac` | `ts_dc` | `ca_ac` | `ca_dc` |
|---|---:|---:|---:|---:|
| case118 | 1,899,303 | 718,719 | 2,134,573 | 692,431 |
| case1354pegase | 27,127,175 | 6,585,505 | 15,221,480 | 6,146,207 |
| case9241pegase | 239,920,456 | 51,417,924 | 239,125,746 | 53,029,954 |

What the baseline says, on case9241pegase:

* a TimeSeries row costs what a cached single solve costs (112M against 111M for
  the `inj` phase above): 98.3% of it is the algorithm, the batch adds 1.7%.
* a ContingencyAnalysis row costs 11% more than a TimeSeries row in AC (more Newton
  iterations from the fixed seed, and a connectivity search per contingency:
  `check_invertible`, 2.25M, 1.7% of the row). In `handle_disconnected_grid` mode a
  second search per contingency is run up front (`_select_ref_slack_and_masks`,
  1.45M per row). In DC, where a row is cheap, that pre-pass is **6.8%** of the row
  (`ca_dc_mask`: 9.29M against 6.35M for `ca_dc`).
* reading the flows back costs 6.3M per row -- 5% of an AC row, but **2.6x a DC
  row**: `compute_amps_flows` walks a column of the row-major voltage matrix per
  branch.
* a DC TimeSeries row is 2.44M, of which the DC solve is ~1.3M: the rest is the
  batch (the Sbus assembly, 0.9M per row on a 20-row run, is one full
  `nb_steps x nb_bus` complex matrix built up front).
* a fresh ContingencyAnalysis per grid2op step costs 1.1G instructions for 8
  contingencies: 8 rows at 125M, plus ~110M of construction (the copy of the grid,
  the rebuild of the solver input, the "n" solve and its analysis).
