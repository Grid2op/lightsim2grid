# Result views that outlive their grid: `reference` vs `reference_internal`

**Design note, not documentation. Switched in 1.3.0: see "Where it stands" at the end.** It
records a lifetime hole in the Python bindings, the policy that closes it, and what was measured
and weighed before switching. It is not part of the built documentation (`docs/*.rst`); the
behaviour users see is described in `docs/result_views.rst`.

## The short version

Most result getters of `LSGrid` (`get_gen_res`, `get_line_res1`, `get_load_res`, ...,
`src/bindings/python/binding_lsgrid.cpp`) are bound with `py::return_value_policy::reference`.
For an Eigen `Ref` / `Map`, pybind11 turns that into a numpy array that **points into the
grid's own buffer and whose `.base` is `None`**: nothing keeps the `LSGrid` alive. When the last
Python reference to the grid goes, the buffer is freed and every array obtained from it reads
freed memory.

The batch classes (`get_voltages`, `get_flows`, ... in `binding_batch.cpp`) use
`py::return_value_policy::reference_internal` instead: still a view, no copy, but `.base` is the
object the method was called on, so it lives as long as any array taken from it.

The proposal: use `reference_internal` for every getter that returns a view into `self`, after
the benchmark and the audit below.

## How it bites

The pattern is a helper that builds a grid, reads a result and returns it:

```python
def _reference():
    grid = ...                       # local
    grid.ac_pf(V0, max_it, tol)
    p, _, _ = grid.get_gen_res()     # view into `grid`, base None
    return p                         # `grid` is collected here: `p` dangles
```

Inside the helper everything checks out — the grid is still alive. Outside, the values are
whatever the allocator put there since: often still the right numbers (the bytes were not
reused yet), sometimes denormal garbage. So a test written that way **passes on one machine
and fails on another**, with no code change in between. This is how
`test_outer_redistribute_slack.py` was red on `test_gcc_latest` only, and green on every
local build tried; it was fixed in the test with a `.copy()`. ASan did not catch it because
the sanitizer job runs a subset of the Python suite that did not include that file.

Copying in the caller is the right fix for a test. It is not a fix for the library: any user
code with the same shape has the same latent bug.

## What pybind11 does with each policy

From `pybind11/eigen/matrix.h` (the `Ref` / `Map` caster, pybind11 3.x):

| policy | numpy array | `.base` | keeps `self` alive |
|---|---|---|---|
| `copy` | owns a copy | — | n/a |
| `reference` | view into the C++ buffer | `None` | **no** |
| `reference_internal` | view into the C++ buffer | `self` | yes |
| `reference` + `py::keep_alive<0, 1>()` | view | depends on the caster | yes (patient list) |

For a getter returning a `std::tuple` of `Ref`s (`tuple3d`, `tuple4d`), the tuple caster passes
the policy and the parent down to each element, so each array of the tuple gets `self` as its
base. Check it rather than trust it: the regression test below asserts it per getter.

## Pros and cons to weigh

**For `reference_internal`:**

- closes the hole for every caller at once, instead of relying on each one to know it must copy;
- keeps the zero-copy view — no new allocation, no memcpy per call — which is why `reference`
  was chosen in the first place;
- the view semantics do not change: an array still shows the next solve's results, exactly as
  today, so code relying on that (if any) keeps working;
- consistent with the batch classes, which already do it.

**Against, or to check:**

- **memory retention.** A single stray result array now keeps the whole `LSGrid` alive — its
  `Ybus`, solver, factorization, every container. On a large real grid, or in a loop building
  many grids (an environment reset, a sweep over snapshots, a test suite), that can turn into
  a real leak-like growth that the current behaviour does not have. Look specifically at
  places that keep result arrays around: `LightSimBackend` caches, grid2op observations,
  `newtonpf`, the batch reference solves in the tests.
- **per-call cost.** Setting the base is an incref on `self`, and a tuple getter does it once
  per element. Expected negligible, but these getters sit on the grid2op step hot path
  (`LightSimBackend` reads every result after every `runpf`), where this repository counts
  microseconds. Measure, do not assume.
- **pickling / copying semantics.** An array with a non-array base pickles as a copy anyway;
  check that nothing downstream (multiprocessing in grid2op, `copy.deepcopy` of a backend)
  starts carrying a grid along, or failing to.
- **`copy` instead.** The safest policy, and the simplest to reason about, but it changes
  behaviour (an array no longer follows the next solve) and adds an allocation per call. Only
  worth considering for getters that are not on a hot path.
- **an unexplained observation.** One local editable build returned `get_gen_res` arrays whose
  base *was* the grid, from the same source and the same pybind11 version as a fresh build that
  returned base `None`. Not understood. Whatever is decided, the regression test must pin the
  behaviour, so that a build can never silently fall back to the other one.

## Benchmark and audit before switching

1. **Micro-benchmark** each result getter of `LSGrid` (and of the containers'
   `get_*` that return views), `reference` vs `reference_internal`, with `timeit`, on a small
   grid and on a large one. Same machine, same build flags, several repeats, report the spread
   and not one number.
2. **Macro-benchmark** the grid2op step loop with `LightSimBackend` (the `benchmarks/`
   scripts), and a `TimeSeries` / `ContingencyAnalysis` run that reads results per row. The
   question is whether the difference is visible at all end to end.
3. **Memory**: a loop that builds and drops grids while keeping one result array per grid,
   watching RSS under both policies; and the same with the arrays dropped, to show the
   retention is only the expected one.
4. **Audit** the Python package for result arrays stored beyond the call (`self._xxx =
   grid.get_..._res()` and the like): each one becomes a strong reference to a grid.
5. **Regression test**: for every view getter, `arr.base is grid` (or `arr.base` reaches it),
   and a helper-returns-a-view test that reads the values after `del grid; gc.collect()` and
   some allocations. It fails today on a CI-style build.

Decide per getter from the numbers: `reference_internal` where the cost is not visible,
`copy` where a getter is cold and the view semantics are not needed, `reference` kept only
with a documented reason. Whatever changes goes with an `[IMPROVED]` (or `[BREAKING]`, if the
view semantics or the retention are judged user-visible) changelog entry.

## Where it stands

Switched in 1.3.0.

- **Point 1** (`benchmarks/view_policy_ab`: the very same getters bound with each policy, timed in
  interleaved rounds and counted under callgrind). `reference_internal` adds a few instructions per
  returned array -- the incref of the grid when the base is set, the decref when the array goes --
  the same at any grid size, and not distinguishable from noise in wall time. `copy` adds an
  allocation and a copy per array, growing with the grid.
- **Point 2: the step loop does not call these getters.** `LightSimBackend._fetch_grid_data`
  fetches each view once and reads it again after every `runpf`; it fetches again only after
  `_reset_res_pointers` (a copy, a reset, `close`, a failed powerflow). The per-call cost is paid
  once per view, not per step, which leaves the end-to-end benchmark nothing per step to see.
- **Point 4, the backend:** its cached views go with `_reset_res_pointers`, which every place
  that replaces or drops `_grid` calls first (`close` included). Under `reference_internal` they
  keep alive only the grid the backend holds anyway.
- **What the audit added: keeping the grid alive is not enough.** A view also needs the grid to
  keep the buffer it shows, and two families did not:
  - the element results. A diverged powerflow reset them by assigning fresh vectors, which freed
    the old buffers (an Eigen move assignment takes over the temporary's storage, and the
    temporary frees the old one). They are now reset in place, to NaN.
  - the solver side (the `*_solver` getters, an algorithm's `get_V` / `get_Va` / `get_Vm`). Later
    powerflows rebuild it in new buffers: a change of topology, a reset, another algorithm, a
    divergence of the DC solver. No policy makes a view of it safe, so these return copies; none
    is on a hot path. `get_Ybus_solver` / `get_dcYbus_solver` had been views since 1.2.0.
- **The decision, per getter** (`binding_lsgrid.cpp`, `binding_containers.cpp`,
  `binding_solvers.cpp`): `reference_internal` for the views of element data (results, setpoints,
  bus ids, admittances, the hvdc droop status) and for the algorithm controls; `copy` for the
  solver side; no policy where the getter returns by value, which pybind11 moves whatever the
  policy says. The containers' `get_bus_id` / `get_yac_eff_*` went from `keep_alive<0, 1>` (safe,
  but invisible in `.base`) to `reference_internal`, and the hvdc lines' `get_bus_id_side_*`,
  which had neither, too.
- **Point 5, the regression tests:** `TestResultViewsKeepTheirGrid` and
  `TestSolverGettersAreSnapshots` (`lightsim2grid/tests/test_LSGrid_out_of_bounds.py`, run under
  ASan by `sanitizers.yml`) sweep every getter of the grid, its containers and its algorithms that
  takes no argument, so that a new getter is checked without being listed; and a C++ test of the
  in-place reset (`src/tests/test_lsgrid.cpp`).
- **Point 3, memory:** the retention is the expected one, and only that. A loop building grids
  and keeping one view of each holds every grid; once the views are dropped no grid is left, and
  the next rounds reuse that memory instead of growing the RSS. A test checks that the grid goes
  with the last array taken from it.
- **Left open:**
  - a view does not survive the re-initialisation of the elements it shows (`init_*`,
    `set_state`), which replaces their buffers. Documented in `docs/result_views.rst`.
  - the batch classes: their views keep the batch alive, but a later `compute` replaces some of
    their buffers -- the flows every time, the voltages when the number of rows changes -- and
    `clear` frees them, so a view held across them reads freed memory. Closing it means copies,
    or result buffers handed over to python at the end of each `compute`; to decide.
