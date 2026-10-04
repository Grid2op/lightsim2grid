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
PP DC                                           105                               8.19                        0.876
pypowsybl                                       444                               1.93                        1.58
DC (SparseLU)                                  2820                               0.0525                      0.00199
DC (KLU)                                       2780                               0.0519                      0.00154
DC (NICSLU\*)                                  2840                               0.0517                      0.00157
DC (CKTSO\*)                                   2840                               0.0519                      0.00162
time serie \*\*                                                                   0.00057823                  0.000300405
PTDF \*\*                                                                         1.52969e-05                 1.44983e-05
contingency analysis \*\*\*                                                       0.0043565                   0.0010985
LODF \*\*\*                                                                       0.00054865                  0.00046915
===========================  ======================  ========================================  ==========================

And for an environment based on the IEEE case 118:

===========================  ======================  ========================================  ==========================
neurips_2020_track2            grid2op speed (it/s)    grid2op 'backend.runpf' time (ms / pf)    time in 'algo' (ms / pf)
===========================  ======================  ========================================  ==========================
PP DC                                            98                               8.77                        1.05
pypowsybl                                       333                               2.64                        2.2
DC (SparseLU)                                  2490                               0.0614                      0.00524
DC (KLU)                                       2500                               0.0595                      0.00371
DC (NICSLU\*)                                  2510                               0.0592                      0.00365
DC (CKTSO\*)                                   2510                               0.0593                      0.00374
time serie \*\*                                                                   0.00272882                  0.00105585
PTDF \*\*                                                                         0.000605272                 0.000586309
contingency analysis \*\*\*                                                       0.00412654                  0.00251818
LODF \*\*\*                                                                       0.000488484                 0.000331645
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

From a grid2op perspective, lightsim2grid allows to compute up to ~2843 DC steps each second (column `grid2op speed`, row `DC (NICSLU\*)`) on the case14_sandbox and "only" ~105 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~27** (2843 / 105) in this case.

When compared to powsybl (with the pypowsybl backend), lightsim2grid is around **~6.4** times faster (444 vs 2843).

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2824 vs 2785 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.352ms (`1. / 2843. * 1000.`) on average and on this 352 µs (or 0.352ms), only 2 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~350 µs are spent in the grid2op extra layer or in the backend implementation in this case (`100%` of the computation time - `=350 / 352`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.000578 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0517 ms for the fastest grid2op DC backend (`DC (NICSLU\*)`), a **~89x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00436 ms on average (row `contingency analysis`), a **~12x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 1.53e-05 ms per powerflow, a **~3383x** speed up compared to the fastest grid2op DC backend.

Using the LODF matrix (row `LODF`) to perform the contingency analysis is faster: 0.000549 ms per contingency, a **~94x** speed up compared to the fastest grid2op DC backend.


For the IEEE case 118:
+++++++++++++++++++++++++

From a grid2op perspective, lightsim2grid allows to compute up to ~2508 DC steps each second (column `grid2op speed`, row `DC (NICSLU\*)`) on the neurips_2020_track2 and "only" ~98 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~26** (2508 / 98) in this case.

When compared to powsybl (with the pypowsybl backend), lightsim2grid is around **~7.5** times faster (333 vs 2508).

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2487 vs 2495 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.399ms (`1. / 2508. * 1000.`) on average and on this 399 µs (or 0.399ms), only 4 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~395 µs are spent in the grid2op extra layer or in the backend implementation in this case (`99%` of the computation time - `=395 / 399`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.00273 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0592 ms for the fastest grid2op DC backend (`DC (NICSLU\*)`), a **~22x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00413 ms on average (row `contingency analysis`), a **~14x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 0.000605 ms per powerflow, a **~98x** speed up compared to the fastest grid2op DC backend.

Using the LODF matrix (row `LODF`) to perform the contingency analysis is faster: 0.000488 ms per contingency, a **~121x** speed up compared to the fastest grid2op DC backend.

See TL;DR section at the top of the file.
