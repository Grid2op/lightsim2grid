# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commits: the DCO signoff is mandatory

Every commit must end with a `Signed-off-by:` trailer naming **a human**, or the DCO check
fails the pull request. This is the single most common reason a PR here is red.

```
Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_...
Signed-off-by: Benjamin Donnot <benjamin.donnot@rte-france.com>
```

Rules:

- **A bot may not sign off.** The DCO is a legal attestation that the contributor has the
  right to submit the work; only a person can make it. `Signed-off-by: Claude
  <noreply@anthropic.com>` is not acceptable even though one or two such commits slipped
  through in the past — do not copy them.
- **The signoff must match the commit author**, so a Claude-written commit is authored by
  the human signing it off, not by Claude. Claude's contribution is recorded in the
  `Co-Authored-By:` trailer, naming the model used.
- Set the identity once per session, then `-s` produces a matching trailer by itself:

  ```
  git config user.name  "Benjamin Donnot"
  git config user.email "benjamin.donnot@rte-france.com"
  git commit -s
  ```

- Adding the trailer to a commit that already exists needs
  `git commit --amend --no-edit -s` (and `--author="..."` if the author is wrong too),
  followed by `git push --force-with-lease`. This is the case that actually bites: the
  mistake is usually noticed only after CI has run.

## Build

The vendored dependencies are git submodules and start empty — a build fails with
`Eigen/Core: No such file or directory` until they are fetched:

```
git submodule update --init --depth 1 eigen SuiteSparse Catch2   # Catch2 only for the C++ tests
```

The build is driven by scikit-build-core. Installing the backend first makes repeated
builds much faster than re-resolving it every time:

```
pip install "scikit-build-core>=0.5.0" pybind11
pip install -e . --no-build-isolation
```

For a quick syntax check of a single C++ change, without a full rebuild:

```
g++ -fsyntax-only -std=c++14 -DLS2G_BUILDING_CORE -I src -I src/core -I eigen \
    -I SuiteSparse/SuiteSparse_config -I SuiteSparse/AMD/Include \
    -I SuiteSparse/BTF/Include -I SuiteSparse/COLAMD/Include -I SuiteSparse/KLU/Include \
    src/core/<file>.cpp
```

Note the language level: the core must compile as **C++14** (CI builds it at both C++14 and
C++26), so no `if constexpr`, structured bindings, `std::optional` or fold expressions.
`NRSystem.hpp` uses the `index_sequence` + dummy-array idiom for variadic folds for exactly
this reason.

## Tests

Python, from `lightsim2grid/tests` — the tests import each other by bare module name, so the
working directory matters:

```
cd lightsim2grid/tests
python -m pytest .                                    # everything
python -m pytest test_ScenarioSweep.py                # one file
python -m pytest test_ScenarioSweep.py::TestScenarioSweepCPP::test_row_count_lock_fails_fast   # one test
python -m pytest . -k "contingency"                   # by name
```

A bare full run reports roughly a hundred failures that have nothing to do with the code:
missing optional dependencies (`pypowsybl`, gymnasium) and grid2op test data
(`data_test/test_PandaPower/test_case14.json`, `data_test/test_multi_chronics`) that a plain
`pip install grid2op` does not ship. Check a failure's message before assuming a change
caused it.

C++ unit tests (Catch2) build standalone — this is the path CI uses, and the one that works.
Configuring the top-level `CMakeLists.txt` with `-DBUILD_TESTING=ON` needs scikit-build-core's
own variables and will not configure on its own:

```
cmake -S src/tests -B build_tests -DCMAKE_BUILD_TYPE=Release
cmake --build build_tests -j$(nproc)
ctest --test-dir build_tests --output-on-failure
ctest --test-dir build_tests -R "some test name"      # one test; they are named by sentence
```

## Architecture

Three layers: a standalone C++ core library (`src/core/`, buildable and usable with no
Python at all — see `docs/cpp_library.rst`), pybind11 bindings (`src/bindings/python/`),
and the Python package (`lightsim2grid/`) whose headline product is `LightSimBackend`, a
grid2op backend.

### `LSGrid` is the hub

`src/core/LSGrid.{hpp,cpp}` owns everything: the element containers, the grid↔solver bus
labelling, and the construction of the solver input (`Ybus`, `Sbus`, the pv/pq split, the
slack). `build_solver_input` / `build_dc_solver_input` produce a `SolverBusLayout` the
algorithms consume; `compute_results` publishes flows and injections back onto the
containers. It is a large file — start from the method you need rather than reading it
through.

### Two bus numbering spaces, and they are different C++ types

A "grid" bus id and a "solver" bus id are not interchangeable, and `TaggedIdVec.hpp` makes
mixing them a **compile error**: `GlobalBusId` / `GlobalBusIdVect` versus `SolverBusId` /
`SolverBusIdVect`, converted through `id_me_to_solver` and `id_solver_to_me`. Deactivated
entries read `BaseConstants::_deactivated_bus_id`. When adding an API, keep the tag — do not
reach for a bare `int` to make a signature compile.

### Element containers

`src/core/element_container/` holds one container per element type (generators, loads,
lines, trafos, shunts, storage, SVCs, HVDC lines and their converter stations), sharing
`GenericContainer` plus the `OneSideContainer` / `TwoSidesContainer` mixins. A container's
job is to stamp itself into `Ybus`/`Sbus` (`fillYbus`, `fillSbus`), declare its contribution
to the pv/pq split (`fillpv`), and receive results.

### `AlgoControl` is the invalidation contract

`src/core/Utils.hpp`. Every element setter that changes something the solver cached must
raise the matching flag — `tell_pv_changed()`, `tell_recompute_ybus()`,
`tell_ybus_change_sparsity_pattern()`, `tell_solver_need_reset()` — and the algorithm reads
those flags to decide what to rebuild (`NRAlgo.tpp`'s `need_rebuild`). This cuts both ways
and both ways are silent: **forget a flag and the solve returns a wrong answer from a stale
cache; raise too much and the performance work in this repo evaporates.** When you add a
setter, find the flag; when you add cached state, find every writer.

### The Newton-Raphson system

`src/core/powerflow_algorithm/` has `BaseAlgo` (the interface), plus DC, fast-decoupled,
Gauss-Seidel and NR implementations. The NR is where the design lives:

- **`NRSystem<Base, Extensions...>`** composes a `Base` block with a variadic tuple of
  extensions (`MultiSlack`, `Hvdc`, `VoltageControl`), each implementing the same component
  protocol — `update_state`, `init_topology`, `register_in`, `declare_feature_entries`,
  `fill_feature_values`, `adjust_mismatch`, `fill_custom_rows`, `apply_step`, `clear`. The
  protocol is documented in a comment block at the top of `NRSystem.hpp`; read it before
  adding a component.
- **`NRLedger`** (`NRLedger.hpp`) is the central registry of equations (Jacobian **rows**)
  and unknowns (Jacobian **columns**). Components claim theirs in registration order, and
  that order *is* the augmented Jacobian's layout. Orientation matters and is easy to get
  backwards: a row is an equation (a P or Q mismatch), a column is an unknown (a Δθ or Δ|V|).
- **`build_J_sparsity`** (`NRSystem.tpp`) is entirely ledger-driven: one pass over the Ybus
  nonzeros emits every dS-derived entry for all components at once. Give a bus a row/column
  in the ledger and the sparsity follows automatically.
- **Value-level row masking** keeps the sparsity fixed across scenarios:
  `set_masked_buses` (P and Q rows → identity, for buses a contingency stranded) and
  `set_pv_pinned_buses` (Q row only, for a bus that is PV in this scenario but PQ in
  another). Both resolve to positions in `J_.valuePtr()` and rewrite values only, so the
  symbolic factorization survives. This is the pattern to reuse for anything that changes a
  bus's role per solve.

### Linear solvers

`src/core/linear_solvers/`: Eigen `SparseLU` (always available), and KLU, NICSLU, CKTSO when
compiled in, behind `LinearSolverPolicy`. The three-way `analyze` (symbolic, expensive) /
`factorize` (numeric) / `refactorize` (numeric, reusing the symbolic analysis and pivot
order) split is the performance story of the whole library — most optimisation work here is
about not triggering `analyze`. `RefactorRetryLinearSolver` falls back to a full `factorize`
when a `refactorize` fails, which is what makes aggressive value-level edits safe.

### Batch algorithms

`src/core/batch_algorithm/BaseBatchSweep.hpp` is one class template with **four**
instantiations, differing only in two policies and an init kind:

| alias | Ybus varies | Sbus varies | each row starts from |
|---|---|---|---|
| `TimeSeries` | no | yes | the previous row's solution |
| `InjectionSweep` | no | yes | the same seed |
| `ContingencyAnalysis` | yes | no | the same seed |
| `ScenarioSweep` | yes | yes | the same seed |

Members that only make sense for some instantiations are SFINAE-gated on
`YbusPolicy::supports_contingency` / `SbusPolicy::supports_vary` rather than split into
subclasses. The performance premise of all four is that **the pv/pq labelling handed to the
solver never changes between rows**, so the batch pays one `analyze` + one `factorize` and
refactorizes thereafter. Anything that would vary the labelling per row has to be expressed
as a value-level mask instead (see the NR section above), or the batch loses its reason to
exist.

### Python package

- `lightSimBackend.py` — the grid2op backend, the main entry point for most users.
- `network/` — grid loaders (pandapower, pypowsybl/iidm, MATPOWER, PowerModels.jl).
- `algorithm/` — algorithm and linear-solver selection, plus the runtime plugin mechanism
  (`AlgorithmRegistry` on the C++ side) that lets an external solver be loaded without
  forking.
- `timeSerie.py`, `injectionSweep.py`, `contingencyAnalysis.py`, `scenarioSweep.py` — thin
  wrappers over the four batch instantiations above.
- `newtonpf/` — a drop-in replacement for pandapower's `newtonpf`.

**Deprecated names still in the tree**, kept working but not for new code: `GridModel` →
`LSGrid`, `lightsim2grid.gridmodel` → `lightsim2grid.network`, `lightsim2grid.solver` →
`lightsim2grid.algorithm`, `SecurityAnalysis` → `ContingencyAnalysis`.

## Changelog

`CHANGELOG.rst`, under the current development version at the top — **not** under `[TODO]`,
which is a list of open problems rather than a release. Entries are prefixed `[ADDED]`,
`[FIXED]`, `[IMPROVED]` and explain *why*, at length where the reasoning is not obvious from
the diff. Match the surrounding style; a one-line entry is out of place here.
