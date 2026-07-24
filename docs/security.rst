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
   * - Loading a grid from the fast binary format (``save_binary`` /
       ``load_binary``)
     - **Hardened against untrusted input.** See below.
   * - Building a grid from a source file (pandapower, pypowsybl, matpower,
       powermodels)
     - **Validated on load.** See below.

Pickle and plugins are trusted-input-only by design
---------------------------------------------------

Unpickling arbitrary data and loading a native plugin are both, fundamentally,
ways of running code chosen by whoever produced the input. No amount of
validation inside lightsim2grid changes that, so the only safe rule is: **do not
load a pickle or a plugin from a source you do not trust.** This is the same
posture as ``numpy``, ``pandas`` (``read_pickle``), ``torch.load`` and any other
library that supports pickling.

The version string embedded in a pickled grid is a *compatibility* check (it
refuses to load a grid saved by a different lightsim2grid version); it is **not**
a security control.

The binary format is the channel to use for untrusted data
----------------------------------------------------------

The fast binary format (:ref:`binary-serialization`) was designed as a safe,
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

:py:meth:`LSGrid.check_grid` (also available as ``GridModel.check_grid``) verifies that a grid
is internally consistent and safe to run a powerflow on. It checks that every
index the grid carries is in range — the bus id of every element, the substation
id and the position in the topology vector (both optional), and the generator
slack and remote-regulated bus references — and that the physical input arrays
(line / transformer ``r``, ``x`` and shunt values, thermal limits, and the
``p`` / ``q`` / voltage set-points) contain no ``NaN`` or infinite value. It
raises ``IndexError`` (an out-of-range index) or ``RuntimeError`` (a structural
inconsistency or a non-finite value), and returns ``None`` for a consistent grid.

You normally do not need to call it yourself:

* it runs automatically when a grid is loaded from a pickle or the binary format
  (via ``set_state``);
* every grid loader (``init_from_pandapower``, ``init_from_pypowsybl``,
  ``init_from_matpower``, ``init_from_powermodels``) calls it before returning.

It is exposed so you can validate a grid you built or modified by hand. It runs
in time proportional to the number of elements in the grid, so it is cheap
compared to a powerflow.

.. note::

   Release wheels are compiled with ``-O3 -DNDEBUG``, which removes Eigen's and
   the standard library's own bounds checks. ``check_grid()`` (and the size /
   finiteness check applied to the voltages a solver returns) are deliberately
   *not* gated behind that flag: they always run, because they are what turns a
   malformed input into a clean exception rather than an out-of-bounds access.

A note for solver-plugin authors
---------------------------------

A solver plugin is trusted code, so lightsim2grid does not try to sandbox it.
It does, however, sanity-check the voltages your ``compute_pf`` returns before
using them: if you report convergence but return a voltage vector of the wrong
size, the powerflow raises ``RuntimeError`` (this would otherwise be an
out-of-bounds read); if you return non-finite values, the solve is reported as
non-converged rather than propagating ``NaN`` / ``Inf`` into the results. Writing
a correct plugin therefore means returning ``V`` / ``Va`` / ``Vm`` sized to the
solver's bus count, and reporting non-convergence yourself when appropriate.
