# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commits: the DCO signoff is mandatory

Every commit must end with a `Signed-off-by:` trailer naming **a human**, or the DCO check
fails the pull request. This is the single most common reason a PR here is red.

```
Assisted-by: Claude Code (claude-opus-5)
Claude-Session: https://claude.ai/code/session_...
Signed-off-by: Benjamin Donnot <benjamin.donnot@rte-france.com>
```

Rules:

- **A bot may not sign off.** The DCO is a legal attestation that the contributor has the
  right to submit the work; only a person can make it. `Signed-off-by: Claude
  <noreply@anthropic.com>` is not acceptable even though one or two such commits slipped
  through in the past — do not copy them.
- **The signoff must match the commit author**, so a Claude-written commit is authored by
  the human signing it off, not by Claude.
- The assistant is credited with **`Assisted-by:`**, naming the tool and the model — not
  `Co-Authored-By:`, and never as an author. This follows the kernel's convention for
  coding assistants (https://docs.kernel.org/process/coding-assistants.html).
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

The build is driven by scikit-build-core, which pip fetches on its own:

```
pip install -e .
```

`--no-build-isolation` (with `pip install "scikit-build-core>=0.5.0" pybind11` first) only
buys not re-resolving the build backend on every rebuild. Worth it when iterating on C++;
not needed otherwise.

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

**grid2op has to be installed from source**, not from pypi: a sizeable part of the suite
reads test data (`data_test/test_PandaPower/test_case14.json`,
`data_test/test_multi_chronics`) that the released wheel does not ship, and those tests fail
for that reason alone. If you see ~100 failures complaining about a missing powergrid or
chronics path, that is this, not the code.

The tests are written with **`unittest`**; run them that way. `pytest` does collect and run
them, but it is not what is used here. From `lightsim2grid/tests` — the tests import each
other by bare module name, so the working directory matters:

```
cd lightsim2grid/tests
python -m unittest discover                                     # everything
python -m unittest test_ScenarioSweep                           # one file
python -m unittest test_ScenarioSweep.TestScenarioSweepCPP.test_row_count_lock_fails_fast
```

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
and the Python package (`lightsim2grid/`). `LightSimBackend` (the grid2op backend) is the
best-known entry point but not the only one: the algorithms themselves (`NR_KLU` and the
rest, via `lightsim2grid.algorithm`), and the batch classes `TimeSeries`, `InjectionSweep`,
`ContingencyAnalysis` and `ScenarioSweep`, are used directly and stand on their own.

### `LSGrid` is the hub

A grid is never built by hand: it is loaded, and `lightsim2grid/network/` holds one
converter per source — pandapower (`init_from_pandapower`), pypowsybl / iidm, MATPOWER and
PowerModels.jl — plus `load_binary` for the fast binary format (`docs/network.rst`,
`docs/binary_serialization.rst`). `LightSimBackend` goes through the pandapower or the
pypowsybl path depending on the grid2op environment. When a bug looks like bad input, check
which converter produced the grid before reading the solver.

`src/core/LSGrid.{hpp,cpp}` owns everything: the element containers, the grid↔solver bus
labelling, and the construction of the solver input (`Ybus`, `Sbus`, the pv/pq split, the
slack). `build_solver_input` / `build_dc_solver_input` produce a `SolverBusLayout` the
algorithms consume; `compute_results` publishes flows and injections back onto the
containers. It is a large file — start from the method you need rather than reading it
through.

### Three bus numbering spaces, and they are different C++ types

`TaggedIdVec.hpp` / `Utils.hpp` give each its own type, so mixing them is a **compile
error** rather than a silent wrong answer:

- **`LocalBusId`** — the busbar *within a substation*, **1-based**, from 1 to
  `n_busbar_per_sub`. This is what a grid2op topology action speaks. Converted with
  `SubstationContainer::local_to_gridmodel(sub_id, LocalBusId)`.
- **`GlobalBusId`** (a.k.a. `GridModelBusId` — the same tag) — the grid-wide bus id.
- **`SolverBusId`** — the bus id inside the matrices actually handed to the solver, which
  only covers connected buses. Converted through `id_me_to_solver` / `id_solver_to_me`.

Deactivated entries read `BaseConstants::_deactivated_bus_id`. When adding an API, keep the
tag — do not reach for a bare `int` to make a signature compile. Getting the 1-based local
labelling wrong has caused real bugs (see the note in `SubstationContainer.hpp`).

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
- **Value-level row masking** keeps the sparsity fixed across scenarios. This is a
  *batch-algorithm* concern — a plain solve, or solves called one after another, has no
  need of it — and it is what lets a sweep reuse one symbolic factorization:
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
subclasses.

The performance premise of all four is that **the sparsity pattern of the Jacobian is fixed
for every element of the batch**, so each algorithm pays one `analyze` + one `factorize` and
refactorizes thereafter. Note what that does *not* say: a row is free to change a bus's
role, as long as it does so by rewriting values inside a pattern that was reserved up front.
Three things already work that way, and more are expected along the same line:

- a bus switching **PV ↔ PQ** because a row disconnected the last generator regulating it
  (`set_switchable_vm_buses` + `set_pv_pinned_buses`);
- a trafo disconnection stranding the **single generator** that was regulating through it;
- a contingency **splitting the grid**, where the smaller part's buses are masked out
  entirely (`set_masked_buses`, `handle_disconnected_grid`).

What is genuinely forbidden is changing the pattern itself per row — a different pv/pq
vector handed to the solver, a different slack *set*. That raises `has_pv_changed()` /
`has_slack_participate_changed()`, forces a fresh `analyze` on every row, and the batch
loses its reason to exist. Reserve the slot up front and mask it instead.

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
`[FIXED]`, `[IMPROVED]`, `[BREAKING]`.

**Keep them short: one to four lines, roughly twenty words.** Read the `[0.13.x]` and
earlier blocks for the register to aim at — say what changed and, where it is not obvious,
why, then stop. The recent entries are far too verbose and are not the model to copy; a
changelog is an index, and the reasoning belongs in the commit message and the code
comments, which is where this repository already puts it at length.
