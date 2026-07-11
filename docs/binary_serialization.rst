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
    bytes are all rejected. ``save_binary`` on the other hand is **atomic**: it writes to a temporary
    file that only replaces the destination once fully written, so an interrupted save never destroys
    a previously saved file. Note that the format stores raw native data: files are meant to be written
    and read by builds sharing the same data layout (same endianness, ``real_type``, ...) -- there is
    no checksum nor cross-platform migration, pickle remains the portable format.

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

- date: 2026-07-02 10:21 CEST
- system: Linux 6.8.0-60-generic
- OS: Ubuntu 22.04.5 LTS
- processor: x86_64
- python version: 3.12.8
- numpy version: 2.0.2
- pandas version: 2.3.3
- pandapower version: 3.4.0
- pypowsybl version: 1.15.0
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.0.0.rc1
- lightsim2grid extra information:

	- klu_solver_available: True
	- nicslu_solver_available: False
	- cktso_solver_available: False
	- compiled_march_native: False
	- compiled_o3_optim: False

This benchmark is available by running, from the root of the lightsim2grid repository:

.. code-block:: bash

    cd benchmarks
    python3 benchmark_binary_serialization.py

Single grid (``l2rpn_idf_2023``)
++++++++++++++++++++++++++++++++++

A first, direct comparison on a single mid-sized grid (the same fixture used by the test suite):

====================================  ================  ===========
l2rpn_idf_2023 (15 repeats)             time (ms/call)    file size
====================================  ================  ===========
save_binary (write)                             0.284      31606 B
load_binary (read)                              0.168      -
pickle.dump (write)                             0.418      32681 B
pickle.load (read)                              0.293      -
====================================  ================  ===========

``save_binary`` is **1.5x** faster than ``pickle.dump`` and ``load_binary`` is **1.7x** faster than
``pickle.load`` on this grid, while also producing a slightly smaller file.

Grid size scan
++++++++++++++++++++++++++++++++++

The benchmark script also scans a range of grid sizes (the same pandapower cases used in
:ref:`benchmark-grid-size`, from 14 to ~9241 buses), each ``LSGrid`` being built directly with
:func:`lightsim2grid.network.init_from_pandapower`. In addition to ``save_binary`` / ``load_binary``
versus pickle, it also reports, as a reference, the time to read the *original* pandapower json **file**
from disk and build an ``LSGrid`` from it from scratch (``pp.from_json`` + ``init_from_pandapower``) --
this is the "no snapshot at all" baseline:

.. code-block:: text

    ================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ================================  ================================
    grid                nb bus    save_binary (ms)    pickle.dump (ms)  write speedup      load_binary (ms)    pickle.load (ms)  read speedup      binary size (B)    pickle size (B)    pandapower json size (B)    read pandapower file (ms, ref)  load_binary speedup vs pp file
    ================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ================================  ================================
    case14                  14               0.315               0.3    1.0x                          0.206               0.273  1.3x                         3368               3054                       69741                            312.12  1515.3x
    case118                118               0.381               0.35   0.9x                          0.271               0.231  0.9x                        19502              21654                      108962                            288.89  1065.5x
    case_illinois200       200               0.356               0.625  1.8x                          0.243               0.279  1.1x                        25936              28600                      127041                            311.49  1280.0x
    case300                300               0.41                0.649  1.6x                          0.212               0.415  2.0x                        43648              48446                      161723                            312.43  1475.3x
    case1354pegase        1354               0.566               1.812  3.2x                          0.324               1.025  3.2x                       196544             229176                      514421                            320.61  989.2x
    case1888rte           1888               0.559               2.286  4.1x                          0.436               1.393  3.2x                       235928             274971                      656133                            326.26  748.4x
    case2848rte           2848               0.698               3.265  4.7x                          0.355               2.119  6.0x                       355939             415935                      951381                            339.28  955.4x
    case2869pegase        2869               0.812               4.005  4.9x                          0.413               2.45   5.9x                       436746             513498                     1061979                            346.86  839.9x
    case3120sp            3120               0.598               3.311  5.5x                          0.427               1.863  4.4x                       334367             398073                      961372                            346.68  812.3x
    case6495rte           6495               1.271              11.11   8.7x                          0.643               5.317  8.3x                       833860             981887                     2140446                            365.18  568.1x
    case6515rte           6515               1.09               11.233  10.3x                         0.66                4.668  7.1x                       836436             984720                     2170208                            385.7   584.7x
    case9241pegase        9241               1.06               18.023  17.0x                         1.061               8.923  8.4x                      1497483            1763983                     3458084                            438.61  413.6x
    ================  ========  ==================  ==================  ===============  ==================  ==================  ==============  =================  =================  ==========================  ================================  ================================

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
