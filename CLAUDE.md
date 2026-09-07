# Working in this repository

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

## Building from a fresh clone

The vendored dependencies are git submodules and start empty — a build fails with
`Eigen/Core: No such file or directory` until they are fetched:

```
git submodule update --init --depth 1 eigen SuiteSparse Catch2   # Catch2 only for the C++ tests
```

The build is driven by scikit-build-core. Installing its backend first makes repeated
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

## Running the tests

Python, from `lightsim2grid/tests` (the tests import each other by bare module name, so the
working directory matters):

```
cd lightsim2grid/tests && python -m pytest .
```

A bare full run reports roughly a hundred failures that have nothing to do with the code:
missing optional dependencies (`pypowsybl`, gymnasium) and grid2op test data
(`data_test/test_PandaPower/test_case14.json`, `data_test/test_multi_chronics`) that a
plain `pip install grid2op` does not ship. Check a failure's message before assuming a
change caused it.

C++ unit tests (Catch2) build standalone — this is the path CI uses, and it is the one
that works; configuring the top-level `CMakeLists.txt` with `-DBUILD_TESTING=ON` needs
scikit-build-core's own variables and will not configure on its own:

```
cmake -S src/tests -B build_tests -DCMAKE_BUILD_TYPE=Release
cmake --build build_tests -j$(nproc)
ctest --test-dir build_tests --output-on-failure
```

## Where the changelog goes

`CHANGELOG.rst`, under the current development version at the top (not under `[TODO]`,
which is a list of open problems rather than a release). Entries are prefixed `[ADDED]`,
`[FIXED]`, `[IMPROVED]` and are written to explain *why*, at length where the reasoning is
not obvious from the diff — match the surrounding style rather than writing a one-liner.
