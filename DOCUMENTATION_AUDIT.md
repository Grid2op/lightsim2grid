# Documentation audit — C++ core / Python bindings (prep for 1.0.0)

Scope: `src/core/help_fun_msg.hpp` + `src/core/help_fun_msg.cpp` (the shared
docstring library) and `src/bindings/python/*.cpp` (the pybind11 bindings that
consume it), on branch `n1_full_compute`. Goal, per the request: (1) every
public, actually-used method should have real documentation instead of
`"TODO"` / "internal, do not use", (2) documentation should live in
`help_fun_msg.hpp/.cpp` rather than being duplicated/hand-written per binding,
(3) flag documentation whose content no longer matches the current code.

**Status: sections 2 (centralization), 3 (18 accuracy bugs) and 4 (missing
docs) are all now fixed** — see the commit(s) following this audit.
`DocComputers`/`DocSecurityAnalysis` have been renamed to
`DocTimeSeries`/`DocContingencyAnalysis` (§5 item 4). A handful of new,
not-previously-cataloged findings surfaced along the way (see the §2/§4
sub-sections below) — most fixed in the same pass, two flagged for a
separate follow-up pass: ~15 more `"TODO"` docstrings on other `LSGrid`
methods (out of scope for "the rest of §4" as originally cataloged), and
the SVC family (`SvcContainer`/`SvcInfo`, no `DocIterator` entries at all
today — explicitly deferred, "lots of docs missing there").

Sphinx pulls these docstrings directly (`autoclass`/`automodule` in
`docs/*.rst`), so anything wrong here is wrong on the public docs site too —
this isn't just an in-REPL `help()` problem.

## 0. Full build verification (post section-3 fixes)

Actually built the C++ extension (`pip install -e ".[docs]"`, after
initializing the `eigen` submodule, which wasn't checked out) and ran the
real `sphinx-build -b html` the project uses:

- **`import lightsim2grid` succeeds**, and the compiled docstrings reflect
  the fixes (spot-checked `get_error.__doc__` directly against the running
  module).
- **`sphinx-build -b html . <out>`: build succeeded, 0 warnings, 0
  errors.** Confirmed on a clean run before *and* after the section-3 fixes.
- Additionally ran `sphinx-build -n` (nitpicky: validates every `:class:`/
  `:func:`/`:attr:` cross-reference, which the project's normal build does
  not) purely as an extra check. It produced 1537 warnings — but this mode
  isn't how the project builds its docs, and cross-checking confirmed the
  volume is a **pre-existing, systemic baseline**, not something introduced
  by the section-3 fixes: KLU/NICSLU/CKTSO solver classes aren't compiled
  in this environment (no SuiteSparse) so their `:class:` refs are
  unresolved here only; every `numpy.typing.NDArray`/`numpy.int32`/
  `numpy.complex128` mention in pybind11-generated type stubs is
  unresolved; and, more relevantly, **every single `:class:`/`:attr:`
  reference to `lightsim2grid.solver.*` throughout the *entire*
  `help_fun_msg.cpp` file (hundreds of them, not just the ones touched in
  this pass) is unresolved**, because `lightsim2grid.solver` is a
  deprecated alias module (`lightsim2grid.algorithm` is current) that isn't
  exposed to autodoc under the old path. That's a real, previously
  undiscovered issue, but it's pre-existing, systemic, and out of scope for
  the section-3 accuracy pass — flagged here for a future decision (rename
  every `lightsim2grid.solver.X` reference to `lightsim2grid.algorithm.X`,
  or fix autodoc to resolve the deprecated alias too), not fixed now.
- One nitpicky-only warning did change as a side effect of a section-3 fix:
  `ContingencyAnalysisCPP.compute_limit_violations`'s docstring newly trips
  a numpydoc quirk where a plain-prose property docstring's first line gets
  misparsed as a comma-separated "type" field, generating spurious
  unresolved-class warnings for text fragments like "per contingency".
  Verified by rebuilding with the pre-fix (stray-quote) text reinstated:
  the same quirk still fires on the untouched `nb_thread` property (and on
  untouched prose in `network.rst:374` / `rewards.rst:10`), so it's the
  same pre-existing quirk — my fix just happened to turn a docstring that
  used to be accidentally truncated (by the stray literal newlines) into a
  single clean line long enough to trip it too. Invisible in the actual
  (non-nitpicky) build, both before and after.
- Spot-checked the actual rendered HTML for the `J_description` rewrite,
  the hvdc-line rewrite, and the line/trafo schema: all render correctly
  (bold section headers, `Note`/`Danger`/`Added in version` admonitions,
  the ASCII schema inside its code block, and cross-references), no
  leftover RST artifacts.

### `_READ_THE_DOCS` was silently broken — fixed

Attempting to actually exercise `_READ_THE_DOCS` (the flag meant to let a
docs build bind `NR_KLU`/`NR_NICSLU`/`NR_CKTSO` and friends without the real,
often-proprietary, linear-solver libraries) surfaced a real, pre-existing
bug unrelated to any docstring content: **the resulting module crashed on
import** with `ImportError: generic_type: type "NR_KLU" is already
registered!`.

Root cause: `src/core/Solvers.hpp` handled a missing library, under
`_READ_THE_DOCS`, by making the name a plain C++ type alias of the
corresponding SparseLU type (`using NR_KLU = NR_SparseLU;` etc.) — not a
distinct type. But `src/bindings/python/binding_solvers.cpp` still
registered each one as its own `py::class_<>`, which pybind11 refuses for a
C++ type that's already bound under another Python name (`NR_SparseLU` is
bound first, unconditionally, a few lines earlier in the same function).

A first fix pointed `NR_KLU` at the already-bound `NR_SparseLU` object
(either as a plain Python attribute alias, or by trying a doc-only Python
subclass). Both were rejected: `_READ_THE_DOCS` exists specifically so each
solver shows its **own** documentation, and neither approach gave `NR_KLU`
a distinct docstring, a distinct `__name__`, or (for the subclass route,
which is also blocked by `NRAlgo`/`BaseDCAlgo`/`BaseFDPFAlgo` all being
declared `final`) its own bound methods.

**Fixed** (final design): each of the three optional linear-solver headers
(`src/core/linear_solvers/{KLU,NICSLU,CKTSO}Solver.hpp`) now defines a
second, mutually-exclusive branch — `#elif defined(_READ_THE_DOCS)` — with
a header-only, no-op-bodied stand-in class of the *same name* as the real
one (e.g. `KLULinearSolver`), implementing the same public interface
(`reset()`/`analyze()`/`factorize()`/`refactorize()`/`solve()`/
`CAN_SOLVE_MAT`) with trivial `{ return ErrorType::NoError; }` bodies. This
makes `NR_KLU` (and `NRSing_KLU`/`DC_KLU`/`FDPF_XB_KLU`/`FDPF_BX_KLU`/
`NRRefactorRetry_KLU`) genuinely distinct C++ types from `NR_SparseLU`, with
their own documented, fully-functional-looking public interface — not
aliases and not hollow Python stand-ins — so the original, simple
`binding_solvers.cpp` code (one `#if defined(KLU_SOLVER_AVAILABLE) ||
defined(_READ_THE_DOCS)` guard, one normal `py::class_<>` registration, no
special-casing) needs no changes at all. `src/core/Solvers.hpp` merges what
used to be three separate `#ifdef X_AVAILABLE` / `#elif _READ_THE_DOCS`
alias blocks into one `#if defined(X_AVAILABLE) || defined(_READ_THE_DOCS)`
block per solver family, since the real and stand-in headers now both
define a genuine `XLinearSolver` type either way.

This also surfaced a second, previously-latent linking bug: `SolverInstantiations.cpp`
(part of `lightsim2grid_core`, always compiled) explicitly instantiates
`BaseFDPFAlgo<LinearSolverPolicy<X>, ...>` — needed because `fillBp_Bpp`'s
body lives out-of-line there, not in the header, since it needs the full
`LSGrid` type — but only under the real `X_SOLVER_AVAILABLE` macros, never
`_READ_THE_DOCS`. And `_READ_THE_DOCS`/`BUILD_DOCS_ONLY` was itself only
ever propagated to the `lightsim2grid_cpp` (bindings) target in the
top-level `CMakeLists.txt`, never to `lightsim2grid_core`. Fixed both:
`SolverInstantiations.cpp`'s three `#ifdef` guards are now
`#if defined(X_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)` (matching
`Solvers.hpp`, including its `extern template` declarations for the
bindings side), and the top-level `CMakeLists.txt` now also applies
`target_compile_definitions(lightsim2grid_core PRIVATE _READ_THE_DOCS)`
when `BUILD_DOCS_ONLY` is set and `lightsim2grid_core` was actually built
from source (not a pre-built `IMPORTED` library, which can't be
recompiled).

**Verified both ways**, real build + real import (not just compile):
- With `_READ_THE_DOCS=1` and none of KLU/NICSLU/CKTSO actually available:
  `import lightsim2grid` succeeds; `NR_KLU is NR_SparseLU` is `False`;
  `NR_KLU.__name__` is `'NR_KLU'`; `NR_KLU.__doc__` is genuinely KLU-specific
  text, different from `NR_SparseLU.__doc__`; `NR_KLU()`, `NR_NICSLU()`,
  `NR_CKTSO()` all instantiate and their no-op methods (`converged()`,
  `get_Va()`) run without crashing.
- Without `_READ_THE_DOCS` (the normal build path, still no KLU/NICSLU/CKTSO
  installed): `import lightsim2grid` still works, `NR_SparseLU` binds
  normally, and `NR_KLU` correctly stays absent (clean `ImportError` on
  attribute access) — the fix doesn't change behavior for the normal build
  path.
- Reran the nitpicky Sphinx pass (`sphinx-build -n`) against the
  `_READ_THE_DOCS` build: build succeeds, no warnings/errors related to
  KLU/NICSLU/CKTSO beyond the same pre-existing `lightsim2grid.solver.*`
  cross-reference-resolution warnings that already affect `NR_SparseLU`
  itself (a class that's always compiled, no `_READ_THE_DOCS` involved) —
  confirmed by grepping the log for `lightsim2grid.solver.NR_SparseLU`,
  which shows the identical warning pattern. Nothing new introduced by this
  fix.

### `lightsim2grid.solver.*` cross-references — fixed (stale module rename)

The pre-existing warnings noted just above turned out to be fixable, not
just pre-existing noise: `lightsim2grid.solver` is a deprecated shim module
(`lightsim2grid/solver/__init__.py`, `DeprecationWarning` on import) that
re-exports everything from `lightsim2grid.algorithm`, which is where these
classes are actually documented (`docs/solvers.rst` uses
`.. automodule:: lightsim2grid.algorithm`, registering Sphinx's targets
under that dotted path). But 121 `:class:`/`:attr:`/`:func:` references
across `src/core/help_fun_msg.cpp` (99) and
`src/bindings/python/binding_enums.cpp` (22) still pointed at the old
`lightsim2grid.solver.*` path — the same "rename not fully propagated"
pattern as the earlier or/ex→1/2 and shunt-admittance findings, just for
the `solver`→`algorithm` module rename instead.

**Fixed**: mechanically replaced `lightsim2grid.solver.` with
`lightsim2grid.algorithm.` across both files (every one of the 121
occurrences was a reference to a class/enum/function that's a real,
documented member of `lightsim2grid.algorithm` — confirmed none were
intentional mentions of the deprecated module itself). Also fixed one
adjacent, independent bug surfaced in the same area: `DocSolver::NRSing_SparseLU`
referenced `:class:`lightsim2grid.algorithm.AlgorithmType.NR_KLU`` (an enum
*member*, wrong path entirely) where every parallel sentence elsewhere
references the solver *class* directly (`:class:`lightsim2grid.algorithm.NR_KLU``)
— corrected to match.

**Verified**: nitpicky Sphinx rebuild — `lightsim2grid.solver.*` warnings:
121 occurrences → 0; total warning count 2352 → 2166 (185 resolved: 121 direct
hits plus their cascading duplicates from multiple inclusion points, e.g.
`docs/algorithm_names.rst` autosummary tables). Spot-checked the rendered
HTML for `NR_KLU` in `solvers.html`: the `:class:` link now resolves to the
correct in-page anchor.

### Remaining nitpicky warnings (2166 → 1970) — triaged and mostly fixed

Broke down the 2166 remaining warnings by category:

| Category | Count | Cause |
|---|---|---|
| numpy/scipy/pandas/stdlib type hints (`numpy.typing.NDArray`, `numpy.float64`, `scipy.sparse.csc_matrix`, `pandas.core.frame.DataFrame`, `collections.abc.Sequence`, `os.PathLike`, ...) | ~1700 | `docs/conf.py` loaded `sphinx.ext.intersphinx` but never set `intersphinx_mapping` — no external inventory at all, for anything. |
| `grid2op.*` / `pypowsybl.*` / `pandapower.*` refs | ~340 | Same root cause — no intersphinx mapping for any of them either. |
| `lightsim2grid.LightSimBackend.LightSimBackend` | 163 | Real bug: one shared docstring constant (`DocLSGrid::_internal_do_not_use`, reused across ~150+ bound methods, plus `DocLSGrid::LSGrid` itself) used the wrong path — `lightsim2grid.LightSimBackend.LightSimBackend` instead of the real module `lightsim2grid.lightSimBackend` (lowercase, per the actual file `lightSimBackend.py`). One wrong line, multiplied by every method that shares the boilerplate. |
| `lightsim2grid.network.LSGrid.get_J` | 14 | Real bug: unlike `get_Va`/`get_Vm`, which each have both a gridmodel-labelled (`get_Va`) and solver-labelled (`get_Va_solver`) form on `LSGrid`, the Jacobian only ever got the solver-labelled `get_J_solver` — no gridmodel-labelled `get_J` was ever added. Three docstrings (`DocSolver::get_J_python`, `DocLSGrid::get_J_python_solver` x2) referenced a `get_J` gridmodel counterpart as if it existed. |
| bare `get_timers_jacobian` | 25 | Real bug: `:func:`get_timers_jacobian`` with no module path at all (should be `lightsim2grid.algorithm.AlgorithmSelector.get_timers_jacobian`). |
| `ls2g::AlgoConfig` (in `get_config`/`set_config`) | 24 | Not from our docstrings — pybind11's auto-generated signature line uses the raw C++ namespaced type name. Would need a binding/type-caster change, not a docstring fix; left as-is (out of scope). |

**Fixed**:
- Added `intersphinx_mapping` to `docs/conf.py` for `python`, `numpy`, `scipy`,
  `pandas`, `grid2op`, `pypowsybl`.
- Corrected `lightsim2grid.LightSimBackend.LightSimBackend` → `lightsim2grid.lightSimBackend.LightSimBackend`
  in both source docstrings.
- Reworded the three `get_J`-related passages to state the true (asymmetric)
  behavior accurately instead of claiming a nonexistent gridmodel-labelled
  `get_J` — this was a real "docs describe functionality that isn't there"
  bug, not just a broken link, so text was rewritten, not just re-pointed.
- Fully qualified the bare `get_timers_jacobian` reference.

**Verified**: nitpicky Sphinx rebuild — all three genuine docstring bugs go
to 0 occurrences (`LightSimBackend` 163→0, `get_J` 14→0, bare
`get_timers_jacobian` 25→0); total warnings 2166 → 1970; normal
(non-nitpicky, non-`_READ_THE_DOCS`) build still succeeds with only the
same pre-existing 6 warnings. Spot-checked the rendered HTML anchor for
`lightsim2grid.lightSimBackend.LightSimBackend` resolves correctly.

The `intersphinx_mapping` addition could **not** be verified end-to-end in
this sandbox: its network policy denies outbound access to
`numpy.org`/`scipy.org`/`docs.python.org`/`pandas.pydata.org`/
`grid2op.readthedocs.io`/`pypowsybl.readthedocs.io` (confirmed via
`ProxyError: 403 Forbidden` in the Sphinx build log), so none of those
`objects.inv` inventories could be fetched here — the ~1700-1900 numpy/
scipy/pandas/grid2op/pypowsybl warnings persist in this environment for
that reason, not because the config is wrong. Read the Docs' actual build
servers have normal internet access, so the mapping should resolve there;
this should be re-checked against a real RTD build once available.

## 1. Quantitative overview

`help_fun_msg.hpp` declares 217 doc strings across 5 structs (all 217 are
defined in the `.cpp`, none dangling either way):

| Struct               | # entries |
|----------------------|-----------|
| `DocIterator`        | 82        |
| `DocLSGrid`          | 63        |
| `DocSolver`          | 39        |
| `DocSecurityAnalysis`| 18        |
| `DocComputers`       | 15        |

Placeholder / non-informative docstrings directly in the binding files:

| File                       | `"TODO`-ish strings | `_internal_do_not_use` | total bindings (`.def*`) |
|----------------------------|:---:|:---:|:---:|
| `binding_lsgrid.cpp`       | 24  | 155 | ~306 |
| `binding_containers.cpp`   | 47  | 0   | ~302 (incl. `def_readonly`) |
| `binding_misc.cpp`         | 13  | 0   | 21   |
| `binding_solvers.cpp`      | 3   | 6   | ~127 |
| `binding_batch.cpp`        | 2   | 1   | 73   |

`binding_lsgrid.cpp` alone tags roughly half of `LSGrid`'s bound methods with
the single generic `DocLSGrid::_internal_do_not_use` string ("Internal, do
not use unless you know what you're doing"). That string is reused verbatim
~155 times.

**That "internal" label is not accurate for most of them.** Cross-referencing
those 156 distinct method names against the actual `lightsim2grid` Python
package (everything under `lightsim2grid/`, excluding `tests/`) shows:

- **81 of the 156** are called directly by the package itself (mostly
  `lightSimBackend.py` and `lightsim2grid/network/from_*`) — e.g.
  `get_bus_vn_kv`, `get_gen_res_full`, `get_line_res1_full`,
  `get_trafo_res1_full`, `init_bus`, `init_powerlines`, `init_generators`,
  `update_gens_p`, `update_loads_p`, `update_topo`, `change_p_load`,
  `deactivate_gen`, `set_substation_names`, `unset_changes`, `set_n_sub`,
  `set_line_pos1_topo_vect`, etc. These are exactly the "public methods
  actually used" the audit is meant to find — they need real docstrings,
  not the internal-do-not-use boilerplate.
- The remaining 75 are not called from the package outside tests (e.g. the
  new `_side1`/`_side2` half-open-line accessors, `debug_get_Bp_python`,
  `is_grid_connected_after_contingency`'s low-level twin, several
  `get_bus_*`/`change_bus_*` per-element accessors). Some of these look
  genuinely internal/plumbing; others (e.g. `get_bus1_powerline` /
  `get_bus2_powerline`, the counterparts of `change_bus1_powerline` which
  *is* used) are simply the read side of a used read/write pair and are
  probably worth real docs too, since they're reachable public API either
  way.

Full list of the 156 method names checked, split used/not-used, is at the end
of this document (§6) for reference when doing the actual rewrite.

## 2. Structural finding: documentation is not uniformly centralized

The stated goal is "every doc lives in `help_fun_msg`". Today that's true for
a majority of `LSGrid`, `DocIterator`-covered element fields, and the solver
family, but several whole classes/areas are documented **only** inline in the
binding `.cpp` files, which is exactly what would need to be duplicated again
for a second (C, Rust, …) binding:

- `TimerJac` and `LinearSolverStats` (`binding_solvers.cpp`) — full class doc
  + all `def_readonly` fields, no `Doc*` entry at all.
- `bind_nr_algo_policies` / `bind_linear_solver_stats` /
  `bind_fdpf_linear_solver_stats` template helpers in `binding_solvers.cpp` —
  ~25 methods (scaling/refactor policy get/set, `get_config`/`set_config`,
  `get_theta_to_J_col`, etc.), all with good inline prose, none in
  `help_fun_msg`.
- `AlgoConfig`, `AlgoControl` (`binding_misc.cpp`) — `AlgoControl` is 13
  methods, all `"TODO"`.
- `PandaPowerConverter` (`binding_misc.cpp`) — 6 methods, **zero** docstrings
  (not even a placeholder string is passed).
- `SvcContainer`/`SvcInfo`, `ConverterStationInfo` (`binding_containers.cpp`)
  — entirely inline strings, no `DocIterator` entries exist for SVC-specific
  fields at all (regulation mode, slope, susceptance limits, target
  setpoints) even though every other element family has one.
- `LimitViolation`, `ViolationElementType`, `LimitViolationType`
  (`binding_batch.cpp`) — inline only.
- `add_pickle`/`add_binary_serialization` (`pickle_helpers.hpp`,
  `binary_helpers.hpp`) — `save_binary`/`load_binary` have solid inline docs,
  but they live in a template helper header, not `help_fun_msg`; low priority
  since it's centralized in one place already (not copy-pasted per class),
  but still outside the stated single-source-of-truth file.
- ~40 individual methods in `binding_lsgrid.cpp` have good, specific inline
  prose that never made it into `help_fun_msg` (`get_reference_slack_bus`,
  `set_ignore_status_global`, `change_shift_trafo`,
  `set_trafo_shift_dependent_rx`, `get_slack_absorbed_solver`,
  `get_controller_q_solver`, `get_p_buses_solver` / `get_q_buses_solver` /
  `get_theta_buses_solver` / `get_vm_buses_solver`,
  `get_hvdc_droop_data_solver`, `get_status_droop_hvdc*`, and others). These
  are actually **better written** than many `Doc*` entries — moving them
  as-is into `help_fun_msg` is mostly mechanical.

## 3. Confirmed accuracy bugs (docstring contradicts current code)

These are verified against the current implementation, not just suspected —
each is either wired to the wrong `Doc*` constant, self-contradictory, or
references a name that no longer exists.

1. **[FIXED]** **`binding_solvers.cpp:382`** — `AlgorithmSelector.get_error` is bound
   with `DocSolver::get_V.c_str()` (copy-paste from the line above/below);
   it should be `DocSolver::get_error`. Right now calling `help()` on
   `get_error` describes "the complex voltage for each bus."

2. **[FIXED]** **`binding_lsgrid.cpp:102`** — `LSGrid.available_default_algorithms` is
   bound with `DocLSGrid::available_algorithm_names` ("a list of string"),
   while the doc that actually matches what it returns —
   `DocLSGrid::available_default_algorithms` ("a list of
   `lightsim2grid.algorithm.AlgorithmType`") — sits unused right next to it
   in `help_fun_msg.cpp:1900-1913` and is instead attached to nothing. The
   two docstrings exist, they're just cross-wired.

3. **[FIXED]** **`binding_containers.cpp:390`** — `HvdcLineInfo.p2_mw` is documented
   with `DocIterator::target_p_or_mw`, which describes the **side-1/"or"**
   target and its sign convention ("if positive, power flows from ex to
   or"). `p2_mw` is `station_side_2.target_p_mw`
   (`HvdcLineContainer.hpp:32,479`) — the side-2 value, wrong side entirely,
   and the described sign convention doesn't apply to it.

4. **[FIXED]** **`binding_containers.cpp:431-432`** — `SubstationInfo.nb_max_busbars`
   and `SubstationInfo.vn_kv` are both bound with `DocIterator::name` (a
   string-field doc). `nb_max_busbars` is an integer (max busbars per
   substation), `vn_kv` is the substation's nominal voltage in kV — neither
   has anything to do with "name."

5. **[FIXED]** **`binding_batch.cpp:60,106,115,123,133`** — five `def_property`
   docstrings (`TimeSeries.init_from_n_powerflow`,
   `ContingencyAnalysis.compute_limit_violations`, `.init_from_n_powerflow`,
   `.handle_disconnected_grid`, `.nb_thread`) use raw-string literals
   (`R"mydelim(...)"`) that still contain the leftover `"` characters and
   mid-string indentation from what used to be ordinary concatenated string
   literals (`"a" "b" "c"`). Converting to a raw string without stripping
   those artifacts means the *actual* compiled docstring shown to users
   contains stray quote marks and ragged indentation, e.g. (verbatim):
   ```
   Whether to initialize the complex voltages of "
                         "the first time series with the results of a n-powerflow "
                         "(*ie* a powerflow at the start the simulation) or not. "
                         "Default: false
   ```

6. **[FIXED]** **`help_fun_msg.cpp` — all 20 `DocSolver::<Family>` solver descriptions**
   (`NR_SparseLU`, `NRSing_SparseLU`, `DC_SparseLU`, `FDPF_XB_SparseLU`,
   `FDPF_BX_SparseLU`, and the `_KLU`/`_NICSLU`/`_CKTSO` equivalents, lines
   ~158-670) each contain a line like *"In the enum
   `lightsim2grid.solver.AlgorithmType`, it is called `SparseLU`"* (or
   `KLU`, `KLUSingleSlack`, `KLUDC`, `FDPF_KLU`, `NICSLUDC`, `CKTSODC`, …).
   None of those short names exist anymore: `binding_enums.cpp` shows the
   *current* `AlgorithmType` values are `NR_SparseLU`, `NRSing_SparseLU`,
   `DC_SparseLU`, `FDPF_XB_SparseLU`/`FDPF_BX_SparseLU`, `NR_KLU`, … — the
   old short names are literally present in the source only as **commented
   out** deprecated aliases that were never re-enabled. Every one of these
   20 docstrings is telling users to look up an enum member that doesn't
   exist.

7. **[FIXED]** **`help_fun_msg.cpp` — all of `DocComputers` (13 occurrences) and
   `DocSecurityAnalysis` (55 occurrences)** reference the pre-rename Python
   names `lightsim2grid.timeSerie.Computers` and
   `lightsim2grid.securityAnalysis.SecurityAnalysisCPP` /
   `SecurityAnalysis`. Both are now deprecation shims
   (`lightsim2grid/securityAnalysis.py` does
   `warnings.warn("...upgrade to SecurityAnalysisCPP > ContingencyAnalysisCPP...")`
   and just re-exports); the current names are
   `lightsim2grid.timeSerie.TimeSerie` / `TimeSeriesCPP` and
   `lightsim2grid.contingencyAnalysis.ContingencyAnalysis` /
   `ContingencyAnalysisCPP`. This affects essentially the entire
   `DocComputers`/`DocSecurityAnalysis` struct — every cross-reference in it
   points at deprecated names. Given this, renaming the structs themselves
   (`DocComputers` → something TimeSeries-named, `DocSecurityAnalysis` →
   `DocContingencyAnalysis`) as part of the uniformity pass would remove the
   temptation to keep writing "security analysis" in new docs too.

8. **[FIXED]** **`help_fun_msg.cpp:3406-3419` — `DocComputers::get_power_flows`** opens
   with "Get the **active** flows (in **MW**)…" but its own `Returns`
   section says `"The flows (in kA)…"` — copy-pasted from
   `DocComputers::get_flows` a few lines above without updating the unit.
   Self-contradictory within the same docstring.

9. **[FIXED]** **`help_fun_msg.cpp:77` — `DocSolver::get_error`** ships literally as:
   `"Returns the error code (as an integer) encountered by the solver (0 =
   no error). TODO DOC explain better"` — the placeholder text is in the
   compiled, user-facing docstring, not a source comment.

Given the density of issues in items 6-7 (every entry in three of the five
structs affected), it's likely that a full pass would turn up more of the
same pattern once someone starts rewriting — these should be treated as
"the whole struct needs a rename-aware pass," not as three isolated typos.

Deeper cross-checks against the current `NRSystem`/`BaseDCAlgo`/element
implementations (done with fresh, independent re-reads of the code, not
just the bindings) turned up more of the same class of problem,
concentrated in `DocLSGrid::J_description`, `DocIterator`'s hvdc-line block,
and a couple of sign-convention docs that flatly contradict each other:

10. **[FIXED]** **`DocLSGrid::J_description` (`help_fun_msg.cpp:2115-2151`) describes a
    Jacobian shape that predates the current solver architecture.** It
    documents J as a fixed 2×2 block plus one extra row/column for the
    distributed slack. The actual builder is a composable
    `NRSystem<Base, Rest...>` (`powerflow_algorithm/NRSystem.hpp:1238-1239`)
    where `MultiSlack` adds `nb_slack` rows/columns (including a dedicated
    "slack_absorbed" column, not a single row), `VoltageControl` adds a
    per-regulation-group block of reactive-power unknowns/equations, and
    `Hvdc` adds droop-flow derivative entries for angle-droop HVDC lines —
    none of which are mentioned. This propagates into
    `DocLSGrid::get_J_python_solver` (`help_fun_msg.cpp:2221-2240`), which
    appends `J_description` verbatim.

11. **[FIXED]** **The same `J_description` claims "all buses in `ref` will be pv buses
    except the first one… this cannot be changed at the moment."** That's
    no longer true: `NRSystem.hpp`'s `Base::register_in` (~line 210-233)
    explicitly supports slack buses that are *not* locally PV-pinned (a
    `free_vm_slack_buses_` set, e.g. a slack bus that's only remotely/SVC
    regulated), giving them a free Vm unknown like an ordinary PQ bus.

12. **[FIXED]** **`DocLSGrid::get_slack_ids_dc` / `get_slack_ids_dc_solver`
    (`help_fun_msg.cpp:2493`, `2605`) say the distributed slack is
    "currently not taken into account" for DC.** `BaseDCAlgo::compute_pf_dc`
    (`powerflow_algorithm/BaseDCAlgo.tpp:108-132`) explicitly spreads the
    active-power imbalance across all participating slack buses
    proportionally to `slack_weights` — i.e. distributed slack for the DC
    powerflow itself *is* implemented today. (The "no distributed slack"
    caveat is still accurate, separately, for `get_ptdf`/`get_lodf`
    specifically — `BaseDCAlgo.tpp:441` still has a live
    `// TODO PTDF: distributed slack` — so don't just delete the caveat
    everywhere, narrow it to the two functions where it still applies.)

13. **[FIXED]** **`DocSolver::get_timers` / `DocSolver::get_computation_time`
    (`help_fun_msg.cpp:137-156`, `~719-723`) name a return field,
    `timer_total_`, that doesn't exist.** `BaseAlgo::get_timers()`
    (`BaseAlgo.hpp:314-318`) returns the 4-tuple `(timer_Fx_, timer_solve_,
    timer_check_, timer_total_nr_)` — the last element is `timer_total_nr_`.
    (The other three names and the ordering are correct; don't confuse this
    with the separate, correctly-documented 13-field `TimerJac` struct
    returned by `get_timers_jacobian()`.)

14. **[FIXED]** **`DocIterator::bus_ex_id` (`help_fun_msg.cpp:1525-1529`) still says
    "lv" instead of "ex."** Used for `LineInfo.bus2_id`
    (`binding_containers.cpp:299`) and `HvdcLineInfo.bus2_id` (`:388`), the
    text reads *"...at which the 'lv' side of the line is connected"* —
    lines/hvdc-lines have `or`/`ex` sides, not `hv`/`lv` (that's
    transformer-only vocabulary). Looks like a copy from `bus_lv_id`'s
    docstring where "transformer" → "line" was updated but "lv" → "ex" was
    missed.

15. **[FIXED]** **`DocIterator::dc_line_formula` and `DocIterator::target_p_or_mw`
    directly contradict each other on sign convention, on the same
    class.** `dc_line_formula` (`help_fun_msg.cpp:1714-1738`) says the hvdc
    line uses "load convention": positive `or_mw` means power is
    *consumed* at `or` (flows `or`→`ex`). `target_p_or_mw`
    (`:1695-1702`, also used on `HvdcLineInfo`) says positive means power
    is *injected* at `or` (flows `ex`→`or`) — the opposite. The source
    comment in `HvdcLineContainer.hpp:29-30` confirms **generator**
    convention, so `dc_line_formula`'s "load convention" claim is the
    backwards one.

16. **[FIXED]** **`DocIterator::DCLineContainer`/`DCLineInfo` (`help_fun_msg.cpp:1594-1638`,
    `:1641-1693`) describe a retired model, and the example code in them is
    broken.** They still describe the old pandapower-only "dc line = 2
    generators" model; none of the current IIDM converter-station
    (VSC/LCC), angle-droop "AC emulation," or resistive-loss
    (`r_ohm`/`nominal_v_kv`) behavior is mentioned
    (`HvdcLineContainer.hpp:59-90`). Worse, the docstring's own usage
    example calls `dcline.bus_or_id`, which doesn't exist — the real
    attribute is `bus1_id` (`binding_containers.cpp:387`) — so the example
    would raise `AttributeError` if a user copy-pasted it.

17. **[FIXED]** **`DocIterator::min_p_mw`/`max_p_mw` and `min_q_mvar`/`max_q_mvar`
    misdescribe how (and whether) they're used.** `min_p_mw`/`max_p_mw`
    (`help_fun_msg.cpp:1027-1043`, only bound on `SGenInfo`) claim they're
    "used in `check_solution` if `check_q_limits` is set to True" — but
    `LSGrid::check_solution_q_values` (`LSGrid.cpp:1080-1139`) only ever
    checks **reactive** limits and never iterates static generators at all;
    these fields are effectively unused today. `min_q_mvar`/`max_q_mvar`
    (`:1009-1025`) claim the opposite extreme — "not taken into account by
    the solver" — but for `GenInfo`/`ConverterStationInfo` they actively
    drive the runtime Q-mismatch split across multiple PV
    generators/converters on the same bus
    (`GeneratorContainer.cpp:512-593`), changing the reported `res_q_mvar`;
    only for `SGenInfo` is "not used" actually accurate.

18. **[FIXED]** **`DocIterator::is_slack` (`help_fun_msg.cpp:988-999`) references a
    renamed solver class**, `lightsim2grid.solver.SparseLUSingleSlack` —
    same class-of-bug as finding 6 above; current name is `NRSing_SparseLU`
    (`binding_enums.cpp:30`, old name commented out at `:56`).

Corroboration note: two independent research passes over `DocIterator` and
`DocLSGrid` both flagged the exact same `AlgorithmType` rename issue (finding
6) and the same `HvdcLineInfo.p2_mw` mislabeling (finding 3) that were found
manually — increasing confidence in both.

**[FIXED] 19. Follow-up: the "or"/"ex" → "1"/"2" rename wasn't actually complete.**
Findings 3, 15 and 16 only fixed the hvdc-specific corner of a wider problem:
`LineInfo`'s own fields (`bus_or_id`/`bus_ex_id`, `res_p_or_mw`/`res_p_ex_mw`,
`res_q_or_mvar`/`res_q_ex_mvar`, `res_v_or_kv`/`res_v_ex_kv`,
`res_a_or_ka`/`res_a_ex_ka`, `res_theta_or_deg`/`res_theta_ex_deg`) were still
named and worded after the retired "or" (origin) / "ex" (extremity)
convention, even though the bound Python attributes are `bus1_id`/`bus2_id`,
`res_p1_mw`/`res_p2_mw`, etc. — the same mismatch as the hvdc case, just not
flagged the first time because nothing there was *factually wrong* (unlike
`p2_mw`), only stale terminology. The shared `line_model` ASCII schema
(`DocIterator::line_model`, used by `r_pu`/`x_pu`/`h_pu` for both lines and
transformers) had the same problem: it labelled its diagram `or bus`/`ex bus`,
`ior`/`iex`, `vor`/`vex`. Fixed by renaming every affected `DocIterator`
constant to the `_1_`/`_2_` convention (`bus_or_id` → `bus_1_id`,
`res_p_or_mw` → `res_p_1_mw`, the hvdc `_dcline`-suffixed ones too, for
naming consistency with the rest of the file) and rewriting the schema to use
`bus 1`/`bus 2`, `i1`/`i2`, `v1`/`v2`, with a note that transformers keep
`hv`/`lv` instead (a real electrical distinction, not just a naming choice —
left untouched).

**[FIXED] 20. The schema's shunt-admittance term (`1/2*h`) was itself wrong,
not just stale terminology.** `TwoSidesContainer_rxh_A::compute_yac` (in
`src/core/element_container/TwoSidesContainer_rxh_A.hpp`, around line
1015-1020) builds ``yac_11 = ys + h1``, ``yac_22 = ys + h2`` from two
**independent** per-side values `h_side_1_`/`h_side_2_` (bound as
`h1_pu`/`h2_pu`) -- not one shared `h` split symmetrically in half at each
side, which is what the old `1/2*h` label on both sides of the schema
claimed. They can legitimately differ (eg an asymmetric line/transformer
imported from pypowsybl). Separately, `DocIterator::h_pu`'s own text had `g`
(conductance) and `b` (susceptance) backwards -- it said "capacitance (real
part) and dielectric conductance (imaginary part)", but
`DataConverter.cpp:258-264` builds `h` as ``g + 1j * b`` (conductance from
`branch_g` is the real part; susceptance from `branch_c`, ie the line
charging capacitance, is the imaginary part) -- confirmed independently by
`h_side_1_`'s imaginary part feeding the FDPF `B''` matrix
(`TwoSidesContainer_rxh_A.hpp:929`, a susceptance-only quantity). Fixed the
schema to show `h1`/`h2` (not `1/2*h` twice) with a note that they're
independent, and fixed `h_pu`'s text to state the correct `g`/`b` ↔ real/
imaginary mapping.

Areas that *were* checked closely and found to still be accurate, so they
don't need touching in the accuracy pass (only the reorganization):
`check_grid`, `change_algorithm`'s DC/AC routing, `check_solution`'s KCL-
mismatch description, `ac_pf`/`dc_pf`'s interface description, `get_ptdf`/
`get_lodf`/`get_Bf` (shapes, row ordering, gridmodel-vs-solver labelling,
and their own distributed-slack caveat), the bus-labelling family
(`get_Va`/`Vm`/`V`, `id_*_to_*`, `get_pv`/`pq`/`slack_ids`/`slack_weights`,
`get_Ybus`/`get_Sbus` and their `_solver`/`dc` variants), and the
hv/lv-side reuse for `TrafoInfo` (unlike the or/ex reuse bug above, the
hv/lv one is actually correct).

## 4. Missing documentation entirely

- `PandaPowerConverter` (`binding_misc.cpp:17-24`): `set_f_hz`, `set_sn_mva`,
  `get_line_param_legacy`, `get_line_param`, `get_trafo_param_pp3`,
  `get_trafo_param_pp2` — no docstring argument passed at all (not even a
  placeholder).
- `AlgoControl` (`binding_misc.cpp:36-49`): class + all 13 methods are
  `"TODO"`.
- Many `TrafoInfo`/`LineInfo` fields — the raw/effective admittance
  components (`yac_11`, `yac_12`, `yac_21`, `yac_22`, `yac_eff_*`, `ydc_*`)
  are all `"TODO doc"` (8 fields × 2 classes = 16 placeholders), as are the
  matching `get_yac_eff_*` accessor methods on the containers.
- `get_bus_id`/`get_bus_id_side_1`/`get_bus_id_side_2` on every
  `*Container` class (Generator/SGen/Load/Storage/Shunt/Trafo/Line) —
  `"TODO doc"`.
- `SubstationContainer`/`SubstationInfo` class docstrings are themselves
  `"TODO"` (separately from the field-level bug in §3.4).

### §4 progress: Gen / HVDC / converter-station iterator docs — fixed

First pass at §4, scoped (per request) to `GenInfo`, `HvdcLineInfo` and
`ConverterStationInfo` — the containers most heavily reworked recently — plus
a survey of every other container's bindings to confirm nothing else in that
category was hiding (result: no other container has a field bound to Python
that isn't already covered by an existing `DocIterator` entry or a cataloged
§4/§2 item above; SVC is intentionally deferred, see below).

Every field on these three classes turned out to already be exposed to
Python — this was purely a documentation gap, not a missing-binding one:

- `GenInfo.voltage_regulator_on` / `.target_q_mvar`: were literal `"TODO"`.
  `target_q_mvar` was a bug, not just a placeholder: the real, generic
  `DocIterator::target_q_mvar` (used by SGen/Load/Storage/Shunt) already
  existed and was simply never wired up for `GenInfo`.
- `GenInfo.regulated_bus_id` ("remote voltage control"): had a one-line
  inline string, not centralized, with no explanation of the
  `voltage_regulator_on` interplay or the pypowsybl-import caveat that
  `docs/network.rst` already documented in prose.
- `ConverterStationInfo` (embedded in every HVDC line, one per side): the
  entire class had never been added to `DocIterator` — class doc,
  `converter_type`, `loss_factor`, `voltage_regulator_on` (TODO),
  `target_q_mvar` (TODO), `power_factor` were all inline/terse, even though
  the C++ header already had a good class-level comment (VSC vs LCC
  behaviour) that had just never made it into the Python-facing docs.
- `HvdcLineInfo`'s newer IIDM/droop fields (`converters_mode`,
  `p_setpoint_mw`, `r_ohm`, `nominal_v_kv`, `droop_enabled`, `droop_p0_mw`,
  `droop_k_mw_per_rad`, `pmax_1to2_mw`, `pmax_2to1_mw`, `status_droop`,
  `station1`/`station2`): all inline one-liners, not centralized.

**Fixed**: added 19 new `DocIterator` entries (`help_fun_msg.hpp`/`.cpp`) and
wired all of the above bindings (`binding_containers.cpp`) to them.
`voltage_regulator_on` is one shared entry covering both `GenInfo` and
`ConverterStationInfo` (identical PV/PQ-setpoint semantics); `target_q_mvar`
was extended (not duplicated) with a note on the voltage-regulation
interplay, since it's now reused by 6 different classes. Every
cross-reference inside a shared string is fully qualified
(`lightsim2grid.elements.<Class>.<attr>`) rather than bare, since a bare
`:attr:` role only resolves against the *current* page's class and several
of these strings are rendered on multiple different classes' pages.

While wiring `ConverterStationInfo`'s (pre-existing, generic) `min_q_mvar` /
`max_q_mvar` docs, surfaced a second, unrelated pre-existing bug of the same
species as the `lightsim2grid.solver.*` one fixed earlier: 14 references to
`lightsim2grid.gridmodel.*` (`gridmodel` is *also* a deprecated shim,
`network` is current -- same rename-not-propagated pattern, different
module) plus 8 of those missing the `LSGrid.` class segment entirely
(`lightsim2grid.network.get_Ybus` instead of
`lightsim2grid.network.LSGrid.get_Ybus`) since `get_Ybus` etc. are methods,
not module-level functions. Fixed both in the same pass.

**Verified**: nitpicky Sphinx rebuild, no `_READ_THE_DOCS` (KLU/NICSLU/CKTSO
not in scope here) -- three new dangling-reference warnings surfaced
immediately on the first pass, all the same pre-existing numpydoc quirk
already documented above (a docstring whose first line has an early colon
followed by a comma-separated list gets misparsed as a type field):
`converters_mode`, `status_droop`, `power_factor`. Reworded their first
lines (colon -> em dash) and reran: those three warnings gone, no others
introduced. `gridmodel`/missing-`LSGrid.` fix: 0 occurrences left (was 4 for
`check_solution` alone, more once the other `get_*Ybus*`/`get_*Sbus*`
references across the solver-family docstrings are counted); total warnings
1306 -> 1273 across this whole pass. Spot-checked rendered HTML anchors for
`GenInfo.regulated_bus_id`, `ConverterStationInfo`, and
`HvdcLineInfo.status_droop` -- all present and correctly linked.
Round-tripped both class doc rename usages (`TimeSeriesCPP`,
`ContingencyAnalysisCPP`) through a real import to confirm non-empty
docstrings post-rename.

### §4 completed: the rest of the catalog

The remaining five §4 items are now all fixed:

- **`PandaPowerConverter`** (`set_f_hz`, `set_sn_mva`, `get_line_param_legacy`,
  `get_line_param`, `get_trafo_param_pp3`, `get_trafo_param_pp2` + the class
  itself): all had zero docstring. Corrected two things the class's own
  `.hpp` comment got wrong along the way: it calls itself "provided as
  examples", but it's actually what `lightsim2grid.network.init_from_pandapower`
  uses for every real pandapower grid load (`network/from_pandapower/initLSGrid.py`
  / `_aux_add_line.py` / `_aux_add_trafo.py`); and `get_trafo_param_pp3`'s own
  comment called itself "for legacy (<= 2.14ish) pandapower" when it's
  actually the opposite -- the *newer* pandapower-3 "advanced grid model"
  path (confirmed against `_aux_add_trafo.py`'s actual version-gated call
  site). Also surfaced, and documented honestly rather than glossed over: `get_line_param`
  (the non-legacy, "returns h1/h2 separately" variant) is marked `// TODO` in
  `DataConverter.cpp` and currently just does a plain 50/50 split of the
  legacy unsplit `h` -- it does not yet read a genuinely asymmetric per-side
  split from pandapower, unlike the `h1`/`h2` model itself (which does support
  it, eg for pypowsybl-imported grids that don't go through this converter
  at all). This is a real, pre-existing implementation gap, not just a
  docstring one -- flagged here rather than fixed, since implementing it is
  a separate task from documenting current behavior accurately.
- **`AlgoControl`** (class + its 12 boolean flags): were all literal `"TODO"`.
  Also fixed two adjacent, previously-uncaught `"TODO"` docstrings on
  `LSGrid.get_ac_algo_controler` / `get_dc_algo_controler` (the only way to
  actually obtain an `AlgoControl` from Python) -- without these two, the
  new `AlgoControl` docs would have been unreachable from the class most
  users would start from. `AlgoControl` also was not previously re-exported
  under `lightsim2grid.algorithm` (unlike its neighbor `AlgoConfig`) nor
  autodoc'd anywhere, so `:class:` references to it would have been
  dangling; added it to `lightsim2grid/algorithm/__init__.py`'s imports and
  `__all__` (mirroring `AlgoConfig` exactly) so the existing
  `.. automodule:: lightsim2grid.algorithm` in `docs/solvers.rst` picks it
  up automatically.
- **`TrafoInfo`/`LineInfo` admittance fields** (`yac_11/12/21/22`,
  `yac_eff_11/12/21/22`, `ydc_11/12/21/22`, 8 fields x 2 classes, plus the 4
  `get_yac_eff_*` container accessor methods x 2 containers): all `"TODO
  doc"`. These are the raw / Kron-reduced-effective / DC-approximation
  two-port admittance matrix entries behind the `line_model` schema fixed
  earlier in this audit (§3 items 15-16) -- traced the actual formulas in
  `TwoSidesContainer_rxh_A.hpp` / `TrafoContainer.cpp` (MATPOWER-manual-
  referenced pi-model + Kron reduction for a half-open branch) to document
  them accurately rather than guessing from the field names.
- **`get_bus_id` / `get_bus_id_side_1` / `get_bus_id_side_2`** (container
  accessor methods, 7 one-sided containers + 2 two-sided): all `"TODO doc"`,
  except `HvdcLineContainer`'s two, which had **no docstring argument
  passed at all** (not even a placeholder) -- the most severe gap found in
  this whole pass. All now share one `DocIterator` entry per
  one-sided/two-sided shape.
- **`SubstationContainer`/`SubstationInfo`** class docstrings: were
  `"TODO"`. Also fixed `LSGrid.get_substations` / `get_voltage_levels`
  (`"TODO"`, not in the original catalog but the only way to reach a
  `SubstationContainer` from Python). Surfaced a real, pre-existing gap
  while fixing this: `:class:`lightsim2grid.elements.SubstationContainer` /
  `SubstationInfo`` were referenced from `docs/network.rst` prose but never
  actually `autoclass`'d anywhere, so both references were dangling; added
  the `.. autoclass::` entries (replacing a hand-written field table that's
  now redundant with the real, `DocIterator`-backed field docs).

Also fixed while in the area (not separately cataloged, but the same
"rename not fully propagated" species as the `solver`/`gridmodel` fixes
above, hit by grep while touching the surrounding accessor docstrings): 9
occurrences of a **fully broken** code example, `from lightsim2grid.gridmodel
import init` / `init(pp_net)`, across `get_lines`/`get_trafos`/`get_generators`/
`get_static_generators`/`get_shunts`/`get_storages`/`get_loads`/`get_dclines`
and one `check_grid`-adjacent docstring. `lightsim2grid.gridmodel.init`
does not exist at all (confirmed by import) -- this is not cosmetic, a user
copy-pasting any of these nine examples would hit an `ImportError`. Fixed to
`from lightsim2grid.network import init_from_pandapower` /
`init_from_pandapower(pp_net)` throughout.

**Not fixed, flagged for a future pass** (found while grepping
`binding_lsgrid.cpp` for the `get_substations` fix, but out of scope for
"the rest of §4" as cataloged): roughly 15 more literal `"TODO"` /
`"TODO doc"` docstrings on other `LSGrid` methods (`timer_last_ac_pf`,
`timer_last_dc_pf`, `get_turnedoff_gen_pv`, `update_slack_weights`,
`update_slack_weights_by_id`, `assign_slack_to_most_connected`,
`get_ignore_status_global`, `get_synch_status_both_side`,
`set_line_names`/`set_dcline_names`/`set_trafo_names`/`set_gen_names`/
`set_load_names`/`set_storage_names`/`set_sgen_names`/`set_shunt_names`/
`set_svc_names`, `change_ratio_trafo`). A new, separate finding, not part of
the original 5-item §4 catalog.

**Verified**: nitpicky Sphinx rebuild (no `_READ_THE_DOCS`) -- caught and
fixed two new dangling-reference-causing issues before reporting done: a
bare `:attr:`line_model`` reference (`line_model` is a shared text fragment
appended to `r_pu`/`x_pu`/`h_pu`'s docstrings, not itself a bound
attribute -- reworded to reference `:attr:`r_pu`` instead, a real
attribute whose docstring includes that text), and one more instance of the
already-documented numpydoc colon-quirk in `ydc_11`'s first line (reworded,
colon -> em dash, same fix as before). After both fixes: 0 new dangling
references from any of this section's additions; total warnings
1273 -> 1267 (net decrease despite substantial new content, since fixing a
`"TODO doc"` placeholder removes zero warnings on its own but the several
newly-added cross-references that DO resolve slightly outweigh the small
number of pre-existing numpy-typing warnings now also visible on methods
that previously had so little docstring content they may not have been
fully processed). Confirmed via real import that every newly-documented
class/method has a non-empty `__doc__` (`AlgoControl`, `PandaPowerConverter`,
`SubstationContainer`/`SubstationInfo`, `LSGrid.get_substations` /
`get_ac_algo_controler`, `LineInfo`/`TrafoInfo`'s `yac_11`/`ydc_11`,
`GeneratorContainer.get_bus_id`, `LineContainer.get_yac_eff_11`). Ran the
existing `test_DataConverter.py` suite (untouched functionally, only its
docstrings changed) -- still passes.

## Section 2 (centralization) — fixed

Every item cataloged in §2 has been moved from inline binding-file strings
into `help_fun_msg.hpp`/`.cpp`, plus the two structural gaps §2 didn't
literally call out but which fall in the same bucket:

- `TimerJac` (class + all 13 fields) and `LinearSolverStats` (class + all
  12 fields), `binding_solvers.cpp` — previously had **zero** docstring at
  all, not even a placeholder (`.def_readonly` with no third argument).
  Traced the actual timer semantics through `NRAlgo.tpp` / `BaseDCAlgo.hpp`
  / `BaseFDPFAlgo.hpp` (which fields the base class sets vs. which only
  NR-based solvers fill in, and which linear-solver-stats field each timer
  mirrors) to document them accurately rather than guessing from field
  names.
- `bind_nr_algo_policies` (19 methods: scaling/refactor policy get/set,
  `get_config`/`set_config`, `get_theta_to_J_col`/`get_vm_to_J_col`/
  `get_q_to_J_col`), `bind_linear_solver_stats` (1 method), and
  `bind_fdpf_linear_solver_stats` (2 methods) — moved as-is, content was
  already good.
- `AlgoConfig` (class + `int_params`/`real_params`), `binding_misc.cpp` —
  moved as-is, including the mutate-by-reassignment warning.
- `LimitViolation`/`ViolationElementType`/`LimitViolationType` (class +
  7 fields), `binding_batch.cpp` — moved as-is (per-enum-value docstrings,
  eg `GRID`/`NOT_SIMULATED`/`DIVERGENCE`, were left inline: low duplication
  risk, unlike to be reused elsewhere).
- ~55 `LSGrid` methods/properties in `binding_lsgrid.cpp` (more than the
  audit's original "~40" estimate) — every substantive (non-`"TODO"`,
  non-`_internal_do_not_use`) inline docstring in the file: the
  underscore-prefixed bookkeeping properties (`_ls_to_orig`, `_orig_to_ls`,
  `_init_kwargs`, `_bus_fusion_rep`), algo-config accessors, per-bus voltage
  limits, slack/PV-PQ bookkeeping, names, every half-open (per-side)
  connect/disconnect method, transformer tap/phase-shift, remote voltage
  control, HVDC angle-droop, and the whole family of solver-internal-state
  getters (`get_p_buses_solver` and friends) used by external solvers
  re-deriving the NR system. Left the ~10 short "DEPRECATED: use X
  instead" redirect one-liners inline (not worth centralizing: minimal
  duplication risk, conventionally kept next to the alias they redirect
  from).
- Found and fixed one small, genuinely dead constant along the way:
  `DocLSGrid::available_algorithm_names` already existed with real content
  but was never wired to any binding (`LSGrid.available_algorithm_names`
  still used an inline literal) -- merged into it instead of creating a
  duplicate.
- Not touched: `add_pickle`/`add_binary_serialization` (`pickle_helpers.hpp`/
  `binary_helpers.hpp`) -- per the audit's own note, low priority, already
  centralized in one template header (not copy-pasted per class), just not
  in `help_fun_msg`.

**Verified**: nitpicky Sphinx rebuild after every batch (solver-family,
misc/batch, then the full `LSGrid` sweep) -- caught and fixed several new
dangling references before calling it done, all bare `:func:`/`:attr:`
roles that don't resolve on the class page they're actually rendered on
(eg `get_timers_jacobian` referenced from `TimerJac`'s own page, where it
doesn't exist -- it's only ever bound on `AlgorithmSelector`, not on the
individual solver classes; `get_linear_solver_stats_bp`/`_bpp` referenced
bare from `get_linear_solver_stats`'s docstring, which renders on
`NR_SparseLU`'s page where those FDPF-only methods don't exist), plus two
more instances of the already-documented numpydoc colon-quirk (`"NR-only:
..."` at the start of three `TimerJac` field docstrings). Final
before/after comparison: the exact same *set* of unresolved-reference
targets before and after this pass (confirmed via diff) -- this batch
introduced zero new dangling-reference categories despite touching ~90
docstrings. Round-tripped ~52 of the newly-centralized attributes through
a real import to confirm non-empty, non-placeholder `__doc__`.

## 5. Recommendation for the rewrite pass

1. Treat §3 (accuracy bugs, 18 items) as the first, independent fix — these
   are outright wrong today regardless of the reorganization, and most (the
   swapped-`Doc*`-constant bugs, the enum-rename bugs, the raw-string
   quote artifacts, the self-contradicting units) are cheap, mechanical
   fixes safe to do without touching the uniformity work. Three are not
   mechanical and need someone who actually knows the current solver
   internals to rewrite, not just patch: `J_description` (items 10-11 —
   needs a description of the `NRSystem<Base, MultiSlack, VoltageControl,
   Hvdc>` composable layout, not a touch-up of the old 2×2-block text),
   the DC distributed-slack claim (item 12 — narrow the caveat to
   PTDF/LODF only, don't just delete it), and the hvdc-line family (items
   15-16 — the model description and its example code need a full rewrite
   for the current IIDM/converter-station/angle-droop model, not an edit).
2. Use the §1 used/unused split to prioritize which of the 156
   `_internal_do_not_use` methods get real prose first: the 81 confirmed
   used by `lightSimBackend`/`network/from_*` first, then decide case by
   case on the remaining 75 (some are genuinely internal setup calls only
   ever invoked once from `gridmodel.py`'s `init_*` sequence and can
   reasonably stay terse/internal; others, like the `_side1`/`_side2`
   half-open accessors, are newer public API for this same release and
   deserve full docs, not "internal").
3. When moving the ~40 well-written inline `binding_lsgrid.cpp` docstrings
   and the solver-policy/`TimerJac`/`LinearSolverStats`/`AlgoConfig` docs
   into `help_fun_msg`, this is mostly a cut-and-paste job (content is
   already good) — the main risk is merge conflicts with any concurrent
   work on `n1_full_compute`, not writing new prose.
4. For `DocComputers`/`DocSecurityAnalysis`, decide up front whether to
   rename the C++ structs (`DocContingencyAnalysis`) to match the
   already-renamed Python/C++ classes, since every string in them needs
   touching anyway for the name fix.
5. Add an SVC-specific block to `DocIterator` (currently entirely absent)
   so `SvcContainer`/`SvcInfo` stop being the one element family documented
   ad hoc.

## 6. Appendix — `_internal_do_not_use`-tagged methods, used vs. not

**Used elsewhere in `lightsim2grid/` (outside `tests/`) — 81, prioritize
these:**
`add_gen_slackbus`, `change_bus_gen`, `change_bus_load`, `change_bus_shunt`,
`change_p_gen`, `change_p_load`, `change_p_shunt`, `change_q_load`,
`change_q_shunt`, `change_v_gen`, `deactivate_bus`, `deactivate_dcline`,
`deactivate_gen`, `deactivate_load`, `deactivate_powerline`,
`deactivate_sgen`, `deactivate_shunt`, `deactivate_storage`,
`deactivate_svc`, `deactivate_trafo`, `get_bus_vn_kv`, `get_gen_res_full`,
`get_init_vm_pu`, `get_line_res1_full`, `get_line_res2_full`,
`get_lines_status`, `get_loads_res_full`, `get_n_sub`,
`get_shunts_res_full`, `get_storage_target_p`, `get_storages_res_full`,
`get_substation_names`, `get_trafo_res1_full`, `get_trafo_res2_full`,
`get_trafo_status`, `init_bus`, `init_bus_status`, `init_dclines`,
`init_generators`, `init_generators_full`, `init_hvdc_lines`, `init_loads`,
`init_powerlines`, `init_powerlines_full`, `init_sgens`, `init_shunt`,
`init_storages`, `init_svcs`, `init_trafo`, `init_trafo_pandapower`,
`reactivate_bus`, `reactivate_gen`, `reactivate_load`, `reactivate_shunt`,
`set_gen_pos_topo_vect`, `set_gen_to_subid`, `set_init_vm_pu`,
`set_line_pos1_topo_vect`, `set_line_pos2_topo_vect`,
`set_line_to_sub1_id`, `set_line_to_sub2_id`, `set_load_pos_topo_vect`,
`set_load_to_subid`, `set_n_sub`, `set_shunt_to_subid`, `set_sn_mva`,
`set_storage_pos_topo_vect`, `set_storage_to_subid`,
`set_substation_names`, `set_trafo_pos1_topo_vect`,
`set_trafo_pos2_topo_vect`, `set_trafo_to_sub1_id`, `set_trafo_to_sub2_id`,
`tell_solver_need_reset`, `unset_changes`, `update_gens_p`, `update_gens_v`,
`update_loads_p`, `update_loads_q`, `update_storages_p`, `update_topo`.

**Not found outside `tests/`/bindings — 75, review case by case:**
`change_bus1_dcline`, `change_bus1_powerline`, `change_bus1_trafo`,
`change_bus2_dcline`, `change_bus2_powerline`, `change_bus2_trafo`,
`change_bus_sgen`, `change_bus_storage`, `change_bus_svc`,
`change_p_dcline`, `change_p_sgen`, `change_p_storage`, `change_q_sgen`,
`change_q_storage`, `change_v1_dcline`, `change_v2_dcline`,
`debug_get_Bp_python`, `debug_get_Bpp_python`, `get_all_shunt_buses`,
`get_bus1_dcline`, `get_bus1_powerline`, `get_bus1_trafo`,
`get_bus2_dcline`, `get_bus2_powerline`, `get_bus2_trafo`, `get_bus_gen`,
`get_bus_load`, `get_bus_sgen`, `get_bus_shunt`, `get_bus_status`,
`get_bus_storage`, `get_bus_svc`, `get_dcline_res1_full`,
`get_dcline_res2_full`, `get_gen_res`, `get_gen_status`,
`get_gen_target_p`, `get_gen_theta`, `get_line_res1`, `get_line_res2`,
`get_line_theta1`, `get_line_theta2`, `get_load_target_p`,
`get_load_theta`, `get_loads_res`, `get_loads_status`,
`get_sgen_target_p`, `get_sgens_res`, `get_sgens_res_full`,
`get_sgens_status`, `get_shunt_target_p`, `get_shunt_theta`,
`get_shunts_res`, `get_shunts_status`, `get_sn_mva`, `get_storage_theta`,
`get_storages_res`, `get_storages_status`, `get_trafo_res1`,
`get_trafo_res2`, `get_trafo_theta1`, `get_trafo_theta2`,
`is_grid_connected_after_contingency`, `reactivate_dcline`,
`reactivate_powerline`, `reactivate_sgen`, `reactivate_storage`,
`reactivate_svc`, `reactivate_trafo`, `remove_gen_slackbus`,
`set_max_nb_bus_per_sub`, `tell_recompute_sbus`, `tell_recompute_ybus`,
`tell_ybus_change_sparsity_pattern`, `update_sgens_p`.
