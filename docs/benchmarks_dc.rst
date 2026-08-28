.. _benchmark-dc-solvers:

Benchmarks (dc solvers)
========================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings when
performing DC powerflow.

If you are interested in other type of benchmark, let us know !

.. note::
  This page is really similar to the page ":ref:`benchmark-solvers`" and if some explanation are missing, they
  can probably be found there or in the more detailed page ":ref:`benchmark-deep-dive`"

  Only summary will be posted here.

TL;DR
---------

.. danger::
    If you want to perform only DC powerflow (a linear model as long as the topology is not modified)
    then you should probably avoid doing some powerflow directly, but rather use linear algebra and the PTDF and LODF
    matrices (see :ref:`ptdf-lodf-section` for what they are and how to get them with lightsim2grid).

When using grid2op, for these small environment, the difference in computation time for an AC or a DC powerflow
is neglectible. Because the Newton-Raphson algorithm has been much more optimized, it is even faster to run
AC powerflow than DC powerflow for the case 14 (~2940 steps per second for AC, see :ref:`benchmark-solvers`,
and ~2750 steps per second for DC, see table below).
For the bigger case 118 the DC environment is slightly faster (~2380 steps per second for DC vs ~2180 steps
per second for AC).

.. note::
  If you want to be faster in grid2op, switching to DC powerflow instead of AC will probably not be
  a good solution if you use lightsim2grid.

Lightsim2grid is still much faster than pandapower (*eg* for case 118, ~2380 steps / s for lightsim2grid and
~129 for pandapower, see table below). A comparison against the pypowsybl DC backend is not part of the current
benchmark tables (unlike the AC comparison in :ref:`benchmark-solvers`).

Last, but not least, if you want to perform DC computations and knows in advance the generations and loads
and the topology of the grid, then you probably should use the PTDF and LODF matrices (:ref:`ptdf-lodf-section`).
With them, using a matrix multiplication (and numpy) you can run (on one CPU core) multiple millions of
DC powerflows each second.

.. note::
  As for :ref:`benchmark-solvers`, ``benchmark_dc_solvers.py`` (run twice by ``benchmarks_dc.sh``, once per
  environment) now prints, right after the tables, the descriptive text ("Description (computation time)" /
  "Description (differences)") computed from the numbers actually measured during that run. Update this page
  (including the TL;DR above) by copy / pasting the tables and that generated text after a new run, instead of
  re-deriving the numbers by hand -- this is what let the TL;DR above drift out of sync with its own tables in
  the first place.

Machine used on the benchmarks
-------------------------------

In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
compared with pandapower and pypowsybl when using grid2op.

All of them has been run on a computer with a the following characteristics:

- date: 2026-08-28 16:51  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.3.5
- pandas version: 2.3.3
- pandapower version: 3.4.0
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.0.0
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

Command to replicate the benchmark on your machine
----------------------------------------------------

To run the benchmark, ``cd`` into the ``benchmarks`` folder and install the dependencies
(we suppose here that you have already installed lightsim2grid):

.. code-block:: bash

  pip install -r req_benchmarks.txt

This will install the required packages to run the benchmark smoothly (most notably `grid2op` and `numba`)
and then you can start the benchmark with the following commands:

.. code-block:: bash

    python3 benchmark_dc_solvers.py --env_name l2rpn_case14_sandbox --no_test --number 8000
    python3 benchmark_dc_solvers.py --env_name l2rpn_neurips_2020_track2_small --no_test --number 8000


Results
---------

For an environment based on the IEEE case 14:

===========================  ======================  ========================================  ==========================
case14_sandbox                 grid2op speed (it/s)    grid2op 'backend.runpf' time (ms / pf)    time in 'algo' (ms / pf)
===========================  ======================  ========================================  ==========================
PP DC                                           135                               6.06                        0.888
pypowsybl                                       652                               1.2                         0.865
DC (SparseLU)                                  2760                               0.0531                      0.00181
DC (KLU)                                       2760                               0.0528                      0.00144
DC (NICSLU\*)                                  2750                               0.0531                      0.00145
DC (CKTSO\*)                                   2770                               0.0527                      0.00146
time serie \*\*                                                                   0.000815316                 0.000301582
PTDF \*\*                                                                         0.000141541                 0.000140672
contingency analysis \*\*\*                                                       0.00423725                  0.00102835
LODF \*\*\*                                                                       0.00060045                  0.0005572
===========================  ======================  ========================================  ==========================

And for an environment based on the IEEE case 118:

===========================  ======================  ========================================  ==========================
neurips_2020_track2            grid2op speed (it/s)    grid2op 'backend.runpf' time (ms / pf)    time in 'algo' (ms / pf)
===========================  ======================  ========================================  ==========================
PP DC                                           125                               6.56                        1.06
DC (SparseLU)                                  2390                               0.0645                      0.00528
DC (KLU)                                       2430                               0.0621                      0.00349
DC (NICSLU\*)                                  2430                               0.0621                      0.00347
DC (CKTSO\*)                                   2420                               0.0625                      0.00365
time serie \*\*                                                                   0.00557626                  0.000930515
PTDF \*\*                                                                         0.000809913                 0.000790683
contingency analysis \*\*\*                                                       0.00447826                  0.0026454
LODF \*\*\*                                                                       0.000408027                 0.000256081
===========================  ======================  ========================================  ==========================

(see the section "Comments" below for details and especially the meaning of \*, \*\* and \*\*\*)

Descriptions
--------------

The tables in the previous sections are a condensed report of different figures more or less comparable (sorry for that...):

The rows:

- **PP DC** reports the computation time when using the pandapower backend of grid2op
- **pypowsybl** reports the timings when using the pypowsybl backend 
- **DC** uses the lightsim2grid DC algorithm with the default Eigen "SparseLU" linear solver
- **DC (KLU)** uses the lightsim2grid DC algorithm with the KLU linear solver
- **DC (NICSLU)** uses the lightsim2grid DC algorithm with the NICSLU linear solver
- **DC (CKTSO)** uses the lightsim2grid DC algorithm with the CKTSO linear solver
- **time serie** uses the lightsim2grid `TimeSerie` module to perform the same computation 
  as the one done with grid2op but in c++ only (this is why there is nothing in the column "grid2op speed (it/s)")
- **PTDF**: uses lightsim2grid to get the PTDF matrix and then numpy to perform the 
  same computation as all of the above from the PTDF matrix. grid2op is not involved either hence the absence of value
  for the "grid2op speed (it/s)" column
- **contingency analysis** reports a different kind of computation, when all the powerlines are disconnected one
  after the other (for given value of loads and generators). There are as many computation here as the number of 
  powerlines (and transformers) on the grid. It does not use grid2op either.
- **LODF** also performs a contingency analysis but it uses the "Line Outage Distribution Factor" matrix to 
  compute it. Just like PTDF it uses lightsim2grid to retrieve the LODF and then uses numpy to perform the
  flows computation from this LODF.

The columns:

- **grid2op speed (it/s)** reports the number of iteration per second that can be performed for each given methods
  (when applicable). It is measured counting only the time of the grid2op environment
- **grid2op 'backend.runpf' time (ms / pf)**:
  
  - for **PP DC**, **DC**, **pypowsybl**, **DC (KLU)**, **DC (NICSLU)** and **DC (CKTSO)** it reports the time
    spent in the grid2op backend
  - for **time serie** and **contingency analysis** : it reports the time to do all the powerflows (including pre processing, 
    post processing, etc.) and the time to compute, from these, the current flows
  - for **PTDF**, and **LODF** it reports the time to perform the 
    closest thing to the above, which in this case would be the time to compute the PTDF / LODF matrix
    using lightsim2grid and the time to compute the flows from these matrix (this last part is only matrix multiplication
    done in numpy)
- **time in 'algo' (ms / pf)**: 

  - for **PP DC**, **DC**, **pypowsybl**, **DC (KLU)**, **DC (NICSLU)** and **DC (CKTSO)** it reports the time
    spent in the algorithm that compute the flows (discarding everything not related in the backend)
  - for **time serie** and **contingency analysis** : it reports the time spent in the algorithm that performs the 
    powerflows
  - for **PTDF**, and **LODF** it reports the time to compute the flows from the PTDF / LODF matrix 
    (this last part is only matrix multiplication done in numpy)


The rows **DC (NICSLU \*)** and **DC (CKTSO \*)** requires lightsim2grid to be built from source.

The rows **time serie \*\*** and **PTDF \*\*** perform the same computation as the above but withtout
the use of grid2op. It is less flexible (in grid2op you could change the topology, apply redispatching etc.)
here you would not be able to. But it is also much faster (especially when using the PTDF)

The rows **contingency analysis \*\*\*** and **LODF \*\*\*** perform a different computation which is 
often denoted by "contingency analysis" or "security analysis" or "N-1" in the power system
community. It consists in disconnecting line one after the other and compute the flows.

Comments
--------

This is the text printed by ``benchmark_dc_solvers.py`` (see the note above) for the two tables above.

For the IEEE case 14:
+++++++++++++++++++++++++

From a grid2op perspective, lightsim2grid allows to compute up to ~2765 DC steps each second (column `grid2op speed`, row `DC (CKTSO\*)`) on the case14_sandbox and "only" ~135 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~21** (2765 / 135) in this case.

When compared to powsybl (with the pypowsybl backend), lightsim2grid is around **~4.2** times faster (652 vs 2765).

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2759 vs 2756 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.362ms (`1. / 2765. * 1000.`) on average and on this 362 µs (or 0.362ms), only 1 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~360 µs are spent in the grid2op extra layer or in the backend implementation in this case (`100%` of the computation time - `=360 / 362`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.000815 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0527 ms for the fastest grid2op DC backend (`DC (CKTSO\*)`), a **~65x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00424 ms on average (row `contingency analysis`), a **~12x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 0.000142 ms per powerflow, a **~372x** speed up compared to the fastest grid2op DC backend.

Likewise, using the LODF matrix (row `LODF`) to perform the contingency analysis takes 0.0006 ms per contingency, a **~88x** speed up compared to the fastest grid2op DC backend.


For the IEEE case 118:
+++++++++++++++++++++++++

From a grid2op perspective, lightsim2grid allows to compute up to ~2432 DC steps each second (column `grid2op speed`, row `DC (KLU)`) on the neurips_2020_track2 and "only" ~125 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~19** (2432 / 125) in this case.

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2392 vs 2432 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.411ms (`1. / 2432. * 1000.`) on average and on this 411 µs (or 0.411ms), only 3 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~408 µs are spent in the grid2op extra layer or in the backend implementation in this case (`99%` of the computation time - `=408 / 411`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.00558 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0621 ms for the fastest grid2op DC backend (`DC (KLU)`), a **~11x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00448 ms on average (row `contingency analysis`), a **~14x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 0.00081 ms per powerflow, a **~77x** speed up compared to the fastest grid2op DC backend.

Likewise, using the LODF matrix (row `LODF`) to perform the contingency analysis takes 0.000408 ms per contingency, a **~152x** speed up compared to the fastest grid2op DC backend.

See TL;DR section at the top of the file.
