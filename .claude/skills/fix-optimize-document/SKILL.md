---
name: fix-optimize-document
description: How to carry out and verify a change in lightsim2grid. Covers fixing a bug (failing test first, then the fix, one commit per bug), a speed optimisation (A/B measurement by instruction count or wall clock), and a documentation fix (check the code before editing the doc). Use whenever a potential bug is found or reported, when fixing one, when optimising performance or benchmarking, when assessing whether tests prove a change, or when correcting docs, docstrings or examples.
---

# Fixing, optimising and documenting in lightsim2grid

First decide which kind of change this is: a bug, a speed optimisation or a doc issue.
A doc issue whose code turns out to be wrong is a bug. A speed-up that changes an
answer is a bug too, or at least a behaviour change that must be stated. Follow the
section for that kind of change. `CLAUDE.md` still applies throughout: commit signoff,
changelog, how tests are run, and the rule against running the full suite.

The working tree usually holds the user's own uncommitted work. Stage only the files
your change touched (`git add <paths>`, never `git add -A` / `-a`), and never stash,
reset or check out over their changes.

**Before the first commit, ask the user where each change goes.** List the bugs and
improvements found so far, and ask which should be grouped in a dedicated pull request
(on a branch of their own) and which can be committed directly to the current branch.
Ask again if new ones turn up later that don't fit the answer already given.

## 1. A potential bug

Take each suspected bug through these steps on its own:

1. **Write a failing test.** Write it before touching the code.
   - Python goes in `lightsim2grid/tests/`, using `unittest`, in the file that
     already covers that area if one exists. C++-only behaviour goes in `src/tests/`
     (Catch2). Use the smallest grid that shows the problem.
   - Give the test a docstring or comment that a reviewer can check against the bug
     report. It names the issue if there is one (`# issue #258: ...`, as the
     existing tests do) or the PR it belongs to, and says in a sentence or two what
     was wrong and what the test checks. Describe the behaviour, not numbers
     measured while investigating (see `CLAUDE.md`).
   - Run it against the **unfixed** code and check that it fails, and that it fails
     **for the reason claimed**. That means an assertion on the wrong value or a
     wrong exception, not an `ImportError`, a missing data file or a typo in the test.
   - If no failing test can be written, the bug is **not confirmed**. Do not fix
     it. Record it for the final report (see below) and move on.
2. **Fix the code.** Fix only this bug. Anything else you notice gets its own pass
   through this section.
3. **Run the test again and get green.** Then run the other test files that exercise
   the same code (see "Tests" in `CLAUDE.md`), so you can show the fix broke nothing
   nearby. If a neighbouring test fails, check whether it also fails without your fix
   before blaming the fix: `git stash push -- <only your files>`, run it, then
   `git stash pop`.
4. **Commit it, one commit per bug.** The commit holds the test, the fix and a short
   `[FIXED]` changelog entry. The message says what was wrong, how the test shows
   it, and why the fix is correct. Follow `CLAUDE.md` for the author, `-s` and
   `Assisted-by:`. Commit it to the branch the user chose for it (see above).

**Final report.** End with two lists:
- **Fixed:** for each bug, the commit hash, a one-line description and the test that
  proves it.
- **Suspected but not confirmed:** for each, the suspicion, where in the code, what
  test was tried and why it did not fail (for example: the path is unreachable from
  the public API, the inputs needed cannot be built, or it failed only for an
  unrelated reason). Do not leave this list out to make the result look tidier.

## 2. A speed optimisation

A speed-up is only claimed once it has been **measured A vs B** on the same machine,
with B returning the **same answers** as A.

- **Instruction count** is the preferred measure. It is exactly reproducible, so a
  small difference is real.
  - For a change in `src/core` on the single-powerflow or batch paths, use what
    already exists: `benchmarks/cache_profiling/ab_test.sh`. It builds A and B from
    a patch script and compares instructions retired *and* the voltage vectors to
    17 significant digits. `DRIVER=batch` covers TimeSeries / ContingencyAnalysis.
    The commands are in `benchmarks/cache_profiling/README.md`.
  - Otherwise, use callgrind with collection limited to the code being benchmarked
    (`CALLGRIND_START_INSTRUMENTATION` / `CALLGRIND_STOP_INSTRUMENTATION`, or
    `--collect-atstart=no --toggle-collect=<function>`), so that loading the grid and
    warm-up do not dilute the number. For a coarse whole-program view, use
    `perf stat -e instructions`.
  - Instruction counts do not see cache misses or the latency of a libm call. If a
    change removes arithmetic or changes memory layout, confirm it on the wall clock
    as well.
- **Wall clock** varies more, so treat it with suspicion:
  `benchmarks/cache_profiling/ab_wallclock.sh`, or a script that alternates A and B
  (A B A B…, never all of A and then all of B), repeats each enough times, and
  reports the minimum *and* the spread. Call it a gain only if the difference is
  clearly bigger than the run-to-run spread. Never run two benchmarks, or a
  benchmark and a build, at the same time.
- **Where the benchmark script goes:** the scratchpad by default. If it is reusable
  and covers a path that `benchmarks/` does not, propose adding it there (with the
  licence header the other files carry) and let the user decide.
- **Report:** a table of A, B and the relative change per case, the tool used,
  whether the answers were identical, and the cases that got *slower*. No measured
  numbers go into code comments or docs (see `CLAUDE.md`); they belong in the commit
  message or the PR.
- **Commit:** one commit per optimisation, on the branch the user chose for it, with
  an `[IMPROVED]` changelog entry.

If the A/B shows no gain, or a loss, say so and do not commit the change.

## 3. A documentation issue

Before changing any doc, docstring or example, **check the code it describes.**

- Read the implementation, and grep for the names involved. A name in the doc may
  have been renamed or deprecated (see "Deprecated names" in `CLAUDE.md`).
- If the doc contains a runnable example, run it.
- Then decide:
  - The **doc is wrong** and the code is right: fix the doc.
  - The **code is wrong** and the doc describes the intended behaviour: this is a
    bug, so go to section 1.
  - **Unclear which one is intended:** ask; do not guess.
- If an `.rst` file was touched, build the docs the way the release check does:
  `cd docs && python -m sphinx -b html -W --keep-going . _build/html`. It must
  build with no warnings.
- The `CLAUDE.md` rules on comments and docs apply: no mention of RTE, no measured
  numbers.
- **Commit:** on the branch the user chose for it, with a changelog entry
  (usually `[IMPROVED]`).
