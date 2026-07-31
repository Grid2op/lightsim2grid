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
