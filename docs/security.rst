.. _security:

Security and trust boundaries
=============================

This page describes what lightsim2grid does, and does **not**, protect you from
when you load a grid or a solver that you did not produce yourself. It is meant
to help you decide which inputs are safe to accept from a third party.

The short version
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Channel
     - Trust required
   * - Unpickling a grid (``pickle.load`` / ``joblib`` / ...)
     - **Trusted input only.** Loading a pickle can execute arbitrary Python
       code *by design* — this is a property of Python's ``pickle`` module, not
       of lightsim2grid. Never unpickle a grid (or anything else) that came from
       an untrusted source.
   * - Loading a solver plugin (:func:`lightsim2grid.load_algorithm_plugin`)
     - **Trusted code only.** A plugin is a native shared library that runs
       in-process; once loaded it can do anything the host process can. Only
       load plugins you built or otherwise trust.
   * - Creating a grid2op environment (``grid2op.make``)
     - **Trusted input only.** An environment is not just data: grid2op reads
       and executes python code shipped with it (its ``config.py``, the classes
       it points to, ...). Only use environments you trust.
   * - Loading a grid from the fast binary format (``save_binary`` /
       ``load_binary``)
     - **Least dangerous, still not a safe channel.** See below.
   * - Building a grid from a source file (pandapower, pypowsybl, matpower,
       powermodels — including through a grid2op environment)
     - **Validated on load.** See below.

Pickle, plugins and grid2op environments are trusted-input-only by design
-------------------------------------------------------------------------

Unpickling arbitrary data, loading a native plugin or using a grid2op
environment are all, fundamentally, ways of running code chosen by whoever
produced the input. No amount of validation inside lightsim2grid changes that,
so the only safe rule is: **do not load a pickle or a plugin, nor call**
``grid2op.make``\ **, on something coming from a source you do not trust.** This
is the same posture as ``numpy``, ``pandas`` (``read_pickle``), ``torch.load``
and any other library that supports pickling.

The version string embedded in a pickled grid is a *compatibility* check (it
refuses to load a grid saved by a different lightsim2grid version); it is **not**
a security control.

The binary format is the least dangerous channel
-------------------------------------------------

.. warning::
    "Least dangerous" is not "safe": loading data you do not trust is a bad
    idea whatever the format, and none of the validation described here makes
    it a good one. The point below is only that this format gives an attacker
    strictly less to work with than pickle does — not that it is meant for
    untrusted input.

The fast binary format (:ref:`binary-serialization`) was designed as an
additive alternative to pickle. Loading a binary file never executes code from
the file, and it is validated at two levels:

* **Byte level.** The reader checks the file's magic number, format version and
  object type, validates every length field against the actual file size
  *before* allocating (so a corrupted size cannot trigger a huge allocation),
  and requires the file to be fully consumed. Malformed, truncated or corrupted
  bytes raise ``RuntimeError``.
* **Semantic level.** After the bytes are read back into a grid, ``check_grid()``
  runs automatically (see below): a file that is byte-wise well-formed but
  *inconsistent* — for example one with an out-of-range bus id — raises a clean
  exception instead of leaving the grid in a state that would cause an
  out-of-bounds memory access on the next powerflow.

``check_grid()``: whole-grid consistency validation
---------------------------------------------------

:py:meth:`lightsim2grid.network.LSGrid.check_grid` verifies that a grid
is internally consistent and safe to run a powerflow on. It checks that every
index the grid carries is in range — the bus id of every element, the substation
id and the position in the topology vector (both optional), and the generator and
SVC slack / remote-regulated bus references. The topology-vector positions are
additionally required to form a permutation of ``[0, dim_topo)`` over exactly the
elements :py:meth:`lightsim2grid.network.LSGrid.update_topo` drives (loads, generators, storage units,
powerlines and transformers): that is what makes them safe to index the
``dim_topo``-sized arrays that method is given, so a shunt, static generator, SVC
or hvdc line claiming a position in that vector is rejected. It checks the two
grid-wide per-unit scalars, ``sn_mva`` (the base power the whole grid is
normalised against) and ``init_vm_pu``, are finite and strictly positive — a
degenerate value there does not make a powerflow fail, it makes it silently
*wrong*. It also checks the *shape* of the grid
itself: that the substation container's own arrays agree with each other (the bus
count, the per-bus status vector and ``n_sub × nmax_busbar_per_sub`` all describe
the same set of buses), and that the bus-id mapping vectors carried alongside the
grid (``_ls_to_orig`` / ``_orig_to_ls`` / ``_bus_fusion_rep``) are in range. That
part matters as much as the per-element checks: the bus count is the *bound* the
element checks are expressed against, so validating elements against a grid whose
own arrays disagree would prove nothing. It raises ``IndexError`` (an out-of-range
index) or ``RuntimeError`` (a structural inconsistency), and returns ``None`` for
a consistent grid.

You normally do not need to call it yourself:

* it runs automatically when a grid is loaded from a pickle or the binary format
  (via ``set_state``);
* every grid loader (``init_from_pandapower``, ``init_from_pypowsybl``,
  ``init_from_matpower``, ``init_from_powermodels``) calls it before returning.

It is exposed so you can validate a grid you built or modified by hand. It runs
in time proportional to the number of elements in the grid, so it is relatively
cheap compared to a powerflow.

What it deliberately does **not** check is the *value* of the per-element
electrical parameters: a line's ``r`` / ``x``, a transformer's tap ratio, an
injection's ``p`` / ``q``. Values that look degenerate are usually legitimate —
``r = 0`` is a lossless branch, a **negative reactance is a series capacitor**,
and negative ``r`` / ``x`` also come out of three-winding-transformer star
equivalents (pandapower's own ``case300`` and ``case9241pegase`` carry them).
Non-finite values are legitimate too: pypowsybl encodes "no thermal limit" as a
non-finite ``limit_a_ka``, and unbounded generator reactive limits are
conventionally ``±Inf``. None of these can corrupt anything either — they are
never used as indices, so they cannot cause an out-of-bounds access; they flow
into ``Ybus`` / ``Bbus`` and the solvers' own guards report the solve as
**non-converged**.

What *is* checked are the two grid-wide scalars ``sn_mva`` and ``init_vm_pu``.
They escape the safety net above: ``sn_mva`` scales results back to MW / MVar
*outside* the solver, so a bad value is invisible to the solver's guards — a NaN
there produced a DC powerflow that reported convergence and returned NaN flows.

.. warning::
    The solver's guards are not a complete safety net for the *derived* results.
    They check the bus voltages, not the branch flows computed from them, so a
    grid whose model is degenerate in a way that leaves the voltages finite can
    still produce a **converged** solve with non-finite flows. A connected branch
    with ``r == 0`` *and* ``x == 0`` is the known case: its admittance
    ``1 / (r + j.x)`` is infinite, and in DC the flows come out non-finite while
    the solve reports success. Such a branch is not rejected on load, because a
    lossless branch can be perfectly legitimate (an ideal step-up transformer),
    and ``init_from_pypowsybl(fuse_zero_impedance_branches=True)`` exists exactly
    to fuse the buses of the ones that are not. **Check that the flows you get
    back are finite** if your grid may contain zero-impedance branches.

.. note::

   Release wheels are compiled with ``-O3 -DNDEBUG``, which removes Eigen's and
   the standard library's own bounds checks. ``check_grid()`` (and the size /
   finiteness check applied to the voltages an *external* solver returns) are
   deliberately *not* gated behind that flag: they always run, because they are
   what turns a malformed input into a clean exception rather than an
   out-of-bounds access.

The memory-safety invariant: raise, never corrupt
--------------------------------------------------

Underneath all of the above is a single design rule the C++ core holds itself
to: **a** :py:class:`~lightsim2grid.network.LSGrid` **must never corrupt process
memory.** Whatever it is handed — a malformed pickle or binary file, a grid
loaded from a source file, or a perfectly ordinary Python API call made in the
wrong order or with an out-of-range id — the worst it is allowed to do is raise
a clean Python exception (``IndexError`` / ``RuntimeError`` / ``ValueError`` /
...). It must not read or write out of bounds, use freed memory, or leave the
grid in a state that makes the *next* call do so. C++ is not a memory-safe
language and release wheels are built ``-O3 -DNDEBUG`` (so neither Eigen's nor
the standard library's own bounds checks are compiled in), so this is enforced
by explicit validation on every entry point rather than by the language: an id
or size that reaches the core from outside is range-checked *before* it is used
to index anything, and a grid restored from a serialized state is put through
``check_grid()`` before it is usable.

This invariant is about *safety*, not *correctness of the answer*: a grid you
built with nonsensical (but in-range and consistent) parameters can still return
a wrong or non-converged powerflow — see the ``check_grid()`` notes above on the
values it deliberately does not police. It is also independent of the trust
boundaries above: those tell you which *inputs* are safe to accept; this tells
you that whatever you do accept, a bug in it surfaces as an exception you can
catch, not as silent heap corruption. If you ever observe a genuine crash
(segfault, ``malloc``/``free`` abort, ASan/UBSan report) instead of an
exception, that is a bug in lightsim2grid — please report it.

The solver hot path holds the same invariant
---------------------------------------------

Everything above is about what enters the grid. The Newton-Raphson solver
(``src/core/powerflow_algorithm/``) holds the same rule on its own inputs.

The arrays a powerflow is called with — ``Ybus``, ``Sbus``, the initial
voltages, ``slack_ids`` / ``slack_weights`` / ``pv`` / ``pq`` — are validated
before the solve: the matrix must be square, the vectors must all have one entry
per bus, at least one slack is required, and every bus id must be in ``[0,
n_bus)`` with ``slack_ids``, ``pv`` and ``pq`` pairwise disjoint. A violation
raises ``RuntimeError`` or ``IndexError``.

The **extensions** — hvdc angle droop, remote voltage control, SVCs,
distributed slack — need more than that, because their data does not arrive
through those arguments at all: the solver pulls it from the grid itself at the
start of each solve. Two properties are therefore checked explicitly, once per
solve:

* every bus id that data carries is in range *before* it is used as an index
  (the solver indexes its bus-keyed tables with these ids directly, and a
  release wheel has no bounds check to fall back on);
* the per-element row / column bookkeeping the solver cached when it last
  rebuilt its Jacobian still matches the data it just pulled. These can only
  disagree if a grid change altered the *number* of hvdc droop lines or voltage
  controllers without the solver being told to rebuild — in which case the
  solve would otherwise index stale tables and write outside the Jacobian.

Both raise ``RuntimeError`` naming the block that is inconsistent. The second
one reports an internal error: reaching it means lightsim2grid failed to signal
a topology change to its own solver, so please report it with the grid that
triggered it.

Like ``check_grid()``, these checks are **not** compiled out by ``-DNDEBUG``;
they run in release wheels. They are placed so that nothing on the per-solve
path scales with the size of the grid:

* the range check on the bus classification (``pv`` / ``pq``) runs only when the
  topology is rebuilt, because that is the only time those ids are used as
  indices — and a rebuild already pays for a symbolic refactorization, which
  dominates it by orders of magnitude;
* what runs on *every* solve is proportional to the number of hvdc droop lines
  and voltage controllers only. On a grid with none of them — which includes
  every grid in the benchmark suite — that is a few dozen integer comparisons,
  independent of the number of buses, and every loop is empty.

Nothing runs per Newton iteration.

Note the split, which mirrors ``check_grid()``'s: what is checked are the
*indices and sizes*, never the electrical *values*. A droop coefficient or a
voltage setpoint that is degenerate is not rejected here — it flows into the
Jacobian and the solve is reported as non-converged, exactly as described in the
``check_grid()`` section above.

A note for solver-plugin authors
---------------------------------

A solver plugin is trusted code, so lightsim2grid does not try to sandbox it.
It does, however, sanity-check the voltages your ``compute_pf`` returns before
using them — for *external* solvers only, so the built-in ones pay nothing for
it: if you report convergence but return a voltage vector of the wrong
size, the powerflow raises ``RuntimeError`` (this would otherwise be an
out-of-bounds read); if you return non-finite values, the solve is reported as
non-converged rather than propagating ``NaN`` / ``Inf`` into the results. Writing
a correct plugin therefore means returning ``V`` / ``Va`` / ``Vm`` sized to the
solver's bus count, and reporting non-convergence yourself when appropriate.
