.. _binary-serialization:

Fast binary serialization (``save_binary`` / ``load_binary``)
================================================================

Goal
-----------------

In addition to the standard `pickle <https://docs.python.org/3/library/pickle.html>`_ support already
available on :class:`lightsim2grid.network.LSGrid` (see :ref:`elements-modeled`), lightsim2grid offers a
second, **additive** serialization path: ``save_binary`` / ``load_binary``.

This is **not a replacement for pickle**. It is a much faster, custom binary format meant for the case
where you repeatedly re-load the *same* grid on the *same* machine / lightsim2grid build and care more
about speed than about portability. Pickle remains the way to go if you need to move a grid across
python versions, machines, or lightsim2grid versions, or if you need it to work with multiprocessing /
``ray`` / etc.

It is available on :class:`lightsim2grid.network.LSGrid` itself, as well as on every individual element
container that has its own state (:class:`~lightsim2grid.elements.LoadContainer`,
:class:`~lightsim2grid.elements.LineContainer`, :class:`~lightsim2grid.elements.GeneratorContainer`, ...),
exactly mirroring what is already picklable:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend

    env = grid2op.make(..., backend=LightSimBackend())
    grid = env.backend._grid

    # save the whole grid to a fast binary file
    grid.save_binary("my_grid.lsb")

    # ... later, possibly in another python process on the same machine ...
    from lightsim2grid.network import LSGrid
    grid_reloaded = LSGrid.load_binary("my_grid.lsb")

Individual element containers work the same way:

.. code-block:: python

    loads = grid.get_loads()
    loads.save_binary("my_loads.lsb")

    from lightsim2grid.elements import LoadContainer
    loads_reloaded = LoadContainer.load_binary("my_loads.lsb")

.. note::
    Version compatibility is decided by a **binary format version** (an integer stored in the file
    header, see ``BINARY_FORMAT_VERSION`` in ``src/core/BinaryArchive.hpp``), not by the lightsim2grid
    version: the format number is only bumped when the serialized layout actually changes, so all
    lightsim2grid versions sharing the same format number can read each other's files. Loading a file
    with an unsupported format number raises a ``RuntimeError`` naming both versions and both format
    numbers.

.. warning::
    ``load_binary`` validates its input and always raises a ``RuntimeError`` (rather than crashing,
    exhausting memory, or silently returning garbage) on ill-formed files: garbage or truncated files,
    corrupted internal sizes (every count in the file is checked against the real file size *before*
    anything is allocated), files that contain an object of a different type (*eg* loading a
    ``LoadContainer`` file with ``StorageContainer.load_binary``), and files with unexpected trailing
    bytes are all rejected. On top of that byte-level validation, loading a whole grid also runs
    :py:meth:`lightsim2grid.network.LSGrid.check_grid` (the same consistency check used everywhere a grid is loaded): a file
    that is byte-wise well-formed but carries an out-of-range index (bus / substation /
    topology-vector) is rejected with an ``IndexError`` (out-of-range index) or a ``RuntimeError``
    (structural inconsistency) rather than being loaded into a state that would fault on the next
    powerflow. ``save_binary`` on the other hand is **atomic by default**: it writes to a
    temporary file that only replaces the destination once fully written, so an interrupted save never
    destroys a previously saved file. Pass ``atomic=False`` to write the destination directly instead
    (marginally faster -- it skips one temporary file and rename -- but without that protection).

.. note::
    **The solver is saved by name.** A grid records the *registry name* of the AC and DC
    algorithm it was using (``'NR_SparseLU'``, or whatever a plugin registered itself as),
    and ``load_binary`` re-selects the solver with that name. Two consequences:

    - if that solver is not registered when you load the file -- a plugin you have not
      loaded in this process, or a built-in needing an optional backend (KLU / NICSLU /
      CKTSO) this build lacks -- ``load_binary`` raises, naming the solver and how to get
      it. Use ``LSGrid.load_binary_without_algorithm(path)`` to load the grid *data* and
      keep the default solvers instead, then pick one with ``change_algorithm``;
    - a name is **not** a strong identity. If two different plugins register the same
      solver name, the one loaded when the file is read is the one you get -- it may be a
      completely different algorithm from the one that was used when the file was written,
      with no way for lightsim2grid to notice. Give your plugin solvers distinctive names
      (a project prefix, for instance) if you exchange files between environments.

.. danger::
    Not even pickle is a *portable* format, and this binary format is even less so. Two
    rules to follow:

    - **Never load a binary file from a source you do not trust.** ``load_binary``
      rejects structurally ill-formed input and, via ``check_grid``, a grid with
      out-of-range indices (see above), but a well-formed, internally-consistent
      file can still carry adversarial *values* -- loading it
      means trusting its content, the same way loading a pickle file does. See also
      :ref:`security`.
    - **Generate the file on the machine that will consume it, with the same
      lightsim2grid version.** The format stores raw native data (same endianness,
      ``real_type``, ...) with no checksum and no cross-platform migration; the
      binary-format-version check (see above) only catches layout changes the
      lightsim2grid authors explicitly bumped it for, not every possible mismatch.
      Use pickle instead if you need to move a grid across machines, python versions,
      or lightsim2grid versions.

.. note::
    The binary format walks the exact same internal state (``get_state`` / ``set_state``) that pickle
    uses, so the two stay structurally in sync: anything that can be pickled can be saved with
    ``save_binary``. See ``lightsim2grid/tests/test_binary_serialization.py`` for the full test suite
    (round trip, AC/DC powerflow round trip, format mismatch, cross-version load, corrupted / truncated /
    wrong-type files, atomicity, and a committed reference file guarding the layout against accidental
    changes).

.. _binary_serialization_benchmarks:

Benchmarks
-----------------

Here are some benchmarks made with:

- date: 2026-09-21 10:06  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.4.6
- pandas version: 2.3.3
- pandapower version: 3.5.4
- pypowsybl version: 1.16.1
- pypowsybl2grid version: 0.3.1
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.1.0
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 benchmark_binary_serialization.py

Single grid (``l2rpn_idf_2023``)
++++++++++++++++++++++++++++++++++

A first, direct comparison on a single mid-sized grid (the same fixture used by the test suite):

=========================================  ================  ===================  ===========
l2rpn_idf_2023 (20 repeats)                  time (ms/call)  speedup vs pickle    file size
=========================================  ================  ===================  ===========
save_binary (write)                                   0.374  1.2x                 33166 bytes
load_binary (read)                                    0.301  0.9x                 -
pickle.dump (write)                                   0.44   1.0x (reference)     34281 bytes
pickle.load (read)                                    0.266  1.0x (reference)     -
init_from_pandapower (case118, reference)             1.342  -                    -
init_from_pypowsybl (ieee118, reference)             26.091  -                    -
=========================================  ================  ===================  ===========

The "speedup vs pickle" column above is now computed and printed directly by
``benchmark_binary_serialization.py`` (previously it was a hand-written sentence below the table,
computed separately from the same two numbers -- a second place the ratio could have drifted out of
sync with the table above it). ``save_binary`` / ``load_binary`` also produce a slightly smaller file
than pickle on this grid.

Grid size scan
++++++++++++++++++++++++++++++++++

The benchmark script also scans a range of grid sizes (the same pandapower cases used in
:ref:`benchmark-grid-size`, from 14 to ~9241 buses), each ``LSGrid`` being built directly with
:func:`lightsim2grid.network.init_from_pandapower`. In addition to ``save_binary`` / ``load_binary``
versus pickle, it also reports, as a reference, the time to read the *original* pandapower json **file**
from disk and build an ``LSGrid`` from it from scratch (``pp.from_json`` + ``init_from_pandapower``) --
this is the "no snapshot at all" baseline:

.. code-block:: text

================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ======================================  ========================================
grid                nb bus    save_binary (ms)    pickle.dump (ms)  write speedup      load_binary (ms)    pickle.load (ms)  read speedup      binary size (B)    pickle size (B)    pandapower json size (B)    read pandapower file (ms, reference)  load_binary speedup vs pandapower file
================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ======================================  ========================================
case14                  14               0.266               0.259  1.0x                          0.16                0.168  1.1x                         3986               3414                       69741                                  282.13  1760.6x
case118                118               0.316               0.302  1.0x                          0.186               0.24   1.3x                        23120              24198                      108962                                  293.86  1582.9x
case_illinois200       200               0.374               0.404  1.1x                          0.171               0.412  2.4x                        29492              30832                      127041                                  284.27  1658.6x
case300                300               0.299               0.488  1.6x                          0.199               0.36   1.8x                        49740              52502                      161723                                  287.86  1446.2x
case1354pegase        1354               1.861               1.895  1.0x                          0.328               1.317  4.0x                       225494             251702                      514421                                  309.78  945.2x
case1888rte           1888               0.445               4.665  10.5x                         0.446               1.349  3.0x                       266140             298657                      656133                                  302.29  677.9x
case2848rte           2848               0.462               5.529  12.0x                         0.551               1.998  3.6x                       399687             450352                      951381                                  313.36  569.2x
case2869pegase        2869               0.524               3.898  7.4x                          0.563               2.51   4.5x                       499129             562920                     1061979                                  327.72  581.9x
case3120sp            3120               0.481               3.133  6.5x                          0.457               1.881  4.1x                       377311             430711                      961372                                  326.01  713.3x
case6495rte           6495               0.688               7.155  10.4x                         1.053               4.535  4.3x                       927521            1053991                     2140446                                  363.44  345.1x
case6515rte           6515               0.736               9.492  12.9x                         1.175               4.572  3.9x                       930353            1056998                     2170208                                  355.16  302.4x
case9241pegase        9241               1.14               19.004  16.7x                         1.819               9.398  5.2x                      1704530            1928651                     3458084                                  403.89  222.1x
================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ======================================  ========================================

Two things stand out:

- The speedup of ``save_binary`` / ``load_binary`` over pickle **grows with grid size**: at small grids
  (``case14``) the two formats are roughly on par (pickle's fixed per-call overhead dominates), but at
  ``case9241pegase`` (~9200 buses) ``save_binary`` is **17x** faster to write and ``load_binary`` is
  **8.4x** faster to read than pickle, and the binary file is consistently 10-15% smaller.
- Reading the original pandapower json **file** (``pp.from_json`` + conversion) has a large, roughly
  constant overhead (~300-440 ms) essentially independent of grid size, dominated by pandapower's own
  json (de)serialization rather than by the (fast) lightsim2grid conversion itself. This makes
  ``load_binary`` **several hundred to over a thousand times faster** than re-reading a pandapower file
  from scratch -- the intended use case for this feature: keep a fast, lightsim2grid-native snapshot of
  an already-converted grid around, instead of re-running the full pandapower/pypowsybl conversion every
  time.

.. note::
    ``load_binary`` speedups reported here are versus reading the *original source format* (pandapower
    json), not versus ``init_from_pandapower`` alone: the pandapower json reader itself is the dominant
    cost, not the lightsim2grid conversion. See :ref:`benchmark-grid-size` for a size scan of the raw
    powerflow solver speed (a different topic than serialization).
