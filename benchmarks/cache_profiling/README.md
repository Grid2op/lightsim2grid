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

### One DFS tree per batch instead of a search per contingency

A contingency was checked for connectivity by a breadth-first search over the whole
Ybus pattern -- once per row in AC (`check_invertible`), and once more per row up
front in `handle_disconnected_grid` mode (`_disconnected_buses`, to list what it
strands). One depth-first search of the base graph per `compute()` now answers every
N-1 in constant time (`BusGraph`: Tarjan's bridges, preorder subtree ranges), lists
the stranded side in time proportional to its size, and leaves the search to the
N-k it cannot settle. A plain DC row never asks -- the split is left to the solver
there, as before -- so the tree is not built for it at all.

A/B, KLU, every row's voltages compared bit for bit (`identical` on all 28 rows):

| grid | `ca_ac` | `ca_ac_mask` | `ca_dc_mask` |
|---|---:|---:|---:|
| case30 | -2.2% | -4.0% | -5.3% |
| case118 | -3.0% | -4.7% | -9.0% |
| case1354pegase | -3.8% | -4.1% | -8.4% |
| case9241pegase | -1.6% | -2.6% | -5.8% |
| case118_fancy | -2.5% | -3.9% | -9.0% |
| case1354pegase_fancy | -3.1% | -3.5% | -8.4% |
| case9241pegase_fancy | -1.4% | -2.4% | -5.8% |

The first version of the change built the tree unconditionally and cost a plain DC
row +1.3% to +8.3% (the verdict computed, then never read); that is why the tree is
now built only where AC or the masked mode consumes it: `ca_dc` then reads +0.0% on
the pegase cases and +0.5% on case30 (a bounds check per row).

### The row loop without its copies

Per row, the loop returned the row's injection by value (`nb_bus` complex, from a
row-major matrix whose row is contiguous) and the layout's slack weights by value
(the same vector, every row), re-marked the algorithm's masks and pinning on rows
that strand and flip nothing (each call a pass over the Jacobian's nonzeros at the
next fill), read the clock around Ybus hooks that compile to nothing on a
TimeSeries, and copied every branch's name into a `std::string` while checking its
current. All views and references now; the algorithm's masks are touched only by a
row that needs them; a name is copied into a violation, not per branch.

A/B, KLU, every row bit for bit identical:

| grid | `ts_ac` | `ts_dc` | `ca_ac` | `ca_ac_mask` | `ca_dc_mask` |
|---|---:|---:|---:|---:|---:|
| case30 | -0.3% | -3.9% | -0.1% | -0.3% | -1.1% |
| case118 | -0.2% | -3.0% | -0.1% | -0.1% | -0.5% |
| case1354pegase | -0.1% | -1.9% | -0.2% | -0.3% | -1.0% |
| case9241pegase | -0.05% | -1.1% | -0.1% | -0.2% | -0.7% |

Small, as the audit said they would be: what a row copies is O(nb_bus) and what it
solves is O(nnz(L+U)). Worth having only because they are free, and because the DC
row is small enough for its copies to show.

Found on the way, by the test suite: the "n" solve of a `ScenarioSweep` with
generator contingencies ran with its switchable buses unpinned. The sweep pins them
before `_finish_preprocessing`, which then reset the algorithm -- and a reset drops
the pinning. The per-row loop re-pinned every row, so the rows were right and the
"n" case was not; the loop touching the pinning only where a row flips a bus is what
exposed it. The reset now comes before the preparation hooks.

### The polar form of the seed, once per seed

`NRSystem::update_state` turned the starting voltage into its angle and magnitude on
every solve -- an atan2 and a hypot per bus, 2.1M instructions on case9241pegase,
1.75% of a row -- and every row of a contingency or injection sweep starts from the
same seed. Asked to (`set_start_polar_cache`, which the `FromSeed` sweeps do), the
system keeps the last starting voltage it was handed and its polar form; a call
with the same *bits* (a memcmp, never a comparison of values) copies the polar form
back, which is exact: it was computed from those bits.

A/B, KLU, every row bit for bit identical:

| grid | `ca_ac` | `ca_ac_mask` | `ts_ac` | `idem` / `inj` / `topo` |
|---|---:|---:|---:|---:|
| case30 | -2.0% | -2.1% | +0.00% | +0.01% / +0.01% / -0.05% |
| case118 | -2.6% | -2.7% | +0.00% | +0.00% |
| case1354pegase | -2.2% | -2.4% | +0.00% | +0.00% |
| case9241pegase | -1.5% | -1.5% | +0.00% | +0.00% |
| case118_fancy | -2.2% | -2.2% | +0.00% | +0.00% |
| case1354pegase_fancy | -1.9% | -2.1% | +0.00% | +0.00% |
| case9241pegase_fancy | -1.3% | -1.3% | +0.00% | +0.00% |

Opt-in because the first version cached unconditionally, and a solve that never
hits -- every single solve, every chained row -- paid the three copies into the
cache for nothing: +0.1% to +0.7% on `idem` / `inj` / `topo` and on `ts_ac`. Off,
the code path is the old one, and the table's last two columns say so.

Measured and **declined** on the way: keeping the Newton loop's mismatch vector as
a member instead of allocating it per solve. One malloc less per row, and
+0.1% to +0.35% on every solve: through a member reference the loop cannot keep
the vector's pointer in a register the way it can for a local, and the sparse
product in `_residual_into` compiled worse. The allocation is cheaper than the
aliasing. Found by bisecting the two halves of the change on case1354pegase.

The chained case was measured and **declined**: a TimeSeries row starts from the
system's own last solution, whose polar form the system already holds -- to a
rounding, not to the bit, since the solution is rebuilt from the angle and the
magnitude. Reusing them was worth -0.4% to -2.9% of a row, with the same iteration
counts everywhere, but the rows drift: 5e-15 pu on case30, 1.4e-12 on
case9241pegase, 2.8e-11 on case9241pegase_fancy, growing along the chain. A warm
start that is not the bits the caller handed in is a different contract, and the
A/B compares traces bit for bit for a reason.

### The injection per row, the flows per row

Two things the size of `nb_steps x nb_bus` went away. The injection matrix
(`SbusPolicy::Vary::assemble`) was built up front, complex, one full row per step
-- 1.3 GB for a year of hourly rows on case9241pegase -- and read once; each row is
now built into a buffer the range owns as the loop reaches it (`fill_row`: the same
accumulation, in the same order, with the same operations, so the row is the row of
that matrix bit for bit), and `get_sbuses()` builds the matrix only if asked. The
flows were computed branch by branch, each branch reading a *column* of the
row-major voltage matrix: a strided pass over every row per branch, two heap
temporaries the size of the batch per branch, and on a big grid a matrix re-read
from memory once per branch. They are now computed row by row: the row's voltages
read once, contiguously, the row of flows written once (`_flows_of_row`).

A/B, KLU, 100 rows, best of 5 on the clock, every trace bit for bit identical:

| grid | `ts_flows` | `ca_flows` | `ts_dc` | `ts_ac` | peak rss, `ts_ac` |
|---|---:|---:|---:|---:|---:|
| case30 | -32% | -30% | -26% | +1.4% | 9.4 -> 8.6 MB |
| case118 | -33% | -33% | -44% | -2.6% | 21.8 -> 17.9 MB |
| case1354pegase | -33% | -38% | -46% | -2.6% | 43.2 -> 32.8 MB |
| case9241pegase | -34% | -34% | -28% | +3.2% | 85.2 -> 71.0 MB |
| case118_fancy | -37% | -38% | -47% | -1.8% | 21.8 -> 18.2 MB |
| case1354pegase_fancy | -35% | -43% | -46% | +0.6% | 43.7 -> 33.0 MB |
| case9241pegase_fancy | -32% | -43% | -31% | +3.1% | 89.5 -> 75.0 MB |

The instruction counts tell the other half of the story: the flows cost **more**
instructions on the small grids (+12% to +17% on case30 to case1354pegase, -2.7% on
case9241pegase -- a scalar loop against Eigen's vectorized column operations) and
are a third faster on the clock everywhere. That is the shape of a memory-bound
change, and why the plan asked for the clock on this item.

The `ts_ac` column on the two 9241-bus grids was measured again on a quiet machine
(the table above was taken during a scheduled antivirus scan), the two binaries run
alternately, best of nine, with the driver's split of the row between the solver's
own timer and the rest: on case9241pegase total 8.41 -> 8.44 ms per row (+0.3%),
everything around the solve -25%, preprocessing 33 -> 21 ms; on case9241pegase_fancy
+3.0%, all of it in the solver's own timer (+3.4%) -- code this change does not
touch, whose instruction count fell 0.2%. That is memory placement, not work, and
it was checked: the row's injection now lives in a 148 KB buffer allocated at the
first row instead of inside a 15 MB matrix allocated before the solver's own
buffers, which moves where it lands against the solver's hot vectors. Shifting the
allocator's placement for BOTH binaries (`MALLOC_MMAP_THRESHOLD_=65536`,
`MALLOC_TOP_PAD_=4194304`) takes the solver gap on that grid from +3.2% to +0.6%,
+0.6% and +0.3% under the three settings tried -- the inter-binary noise floor
(the contingency phase, near-identical code in the two builds, reads +0.4% the
same way). Kept: a 1.3 GB matrix on a year of rows and a third off every flow
computation, against a placement effect on one grid that any allocation before
the solver's could flip either way.

A caveat for `-march=native` builds only: the complex division by `sn_mva` is
Eigen's vectorized one, applied to a row of `nb_bus` entries instead of the whole
matrix. With a two-complex packet (AVX) the elements that fall to Eigen's scalar
tail can differ between the two layouts, so a last-bit difference on the last bus
of odd-sized rows is possible there. The shipped build (SSE2, one complex per
packet) has no tail and the traces are identical.

### Registering contingencies, and two steps measured and declined

`ContingencyAnalysis.add_all_n1_contingencies` registered an N-1 sweep one
powerline at a time -- one crossing into C++ and one `clear_results_only()` per
powerline, each resetting the solver and dropping the result buffers even when
nothing had been computed since the last one. `clear_results_only()` is now a
no-op when there is nothing to drop, and the wrapper registers the whole sweep in
one call; a contingency given by name is resolved through a dictionary rather than
a comparison against every name of the grid.

Two items of the plan were sized on the `ca_construct` profile (a fresh
`ContingencyAnalysis` per grid2op step, 8 contingencies, case9241pegase, 1.1G
instructions per step) and **declined**:

* **A move constructor for `LSGrid`.** The copy the batch takes of the grid is
  2.35M instructions of those 1.1G (0.2%); the two solvers it rebuilds and the
  cache it starts cold are what a copy costs, and a move would rebuild them the
  same way. On the Python side the copies are of the *environment*, and a pybind
  constructor cannot take a Python-owned grid by move at all. Nothing to gain that
  the profile can see.
* **Keeping the batch's solver cache and factorization across two `compute()`
  calls.** What it would skip is the input build (14.6M), one symbolic analysis
  (40M) and one factorization (30M): 85M against the 110M of the "n" solve that
  has to run either way, so a third of a one-row repeat and 8% of an eight-row
  one -- and nothing for a fresh object, which is what the Python wrappers build
  per step. Against that, every per-`compute()` configuration that changes the
  Jacobian's sparsity (the switchable buses of a generator-contingency sweep, the
  masked-controller slot) and the DC solver's cached injection would have to be
  invalidated by hand, and a reused factorization refactorizes where a fresh one
  factorizes, so the second `compute()` would no longer match a fresh object bit
  for bit. The natural home for this is an `update_grid()` on the batch, which is
  what would make reuse pay in a grid2op loop; out of scope here.
