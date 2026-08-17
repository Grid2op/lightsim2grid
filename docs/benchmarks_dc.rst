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

- date: 2026-08-06 18:17  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.3.5
- pandas version: 2.3.3
- pandapower version: 3.4.0
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.0.0.rc3
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
PP DC                                           138                               5.91                        0.866
pypowsybl                                       713                               1.08                        0.741
DC (SparseLU)                                  2850                               0.0516                      0.00178
DC (KLU)                                       2780                               0.0525                      0.00155
DC (NICSLU\*)                                  2850                               0.0513                      0.00146
DC (CKTSO\*)                                   2840                               0.0514                      0.00149
time serie \*\*                                                                   0.00120099                  0.000382153
PTDF \*\*                                                                         7.15695e-05                 7.07186e-05
contingency analysis \*\*\*                                                       0.00633635                  0.0011499
LODF \*\*\*                                                                       0.00042515                  0.0003688
===========================  ======================  ========================================  ==========================

And for an environment based on the IEEE case 118:

===========================  ======================  ========================================  ==========================
neurips_2020_track2            grid2op speed (it/s)    grid2op 'backend.runpf' time (ms / pf)    time in 'algo' (ms / pf)
===========================  ======================  ========================================  ==========================
PP DC                                           127                               6.46                        1.05
DC (SparseLU)                                  2460                               0.0623                      0.00498
DC (KLU)                                       2430                               0.0616                      0.00328
DC (NICSLU\*)                                  2410                               0.0619                      0.00329
DC (CKTSO\*)                                   2460                               0.0607                      0.00334
time serie \*\*                                                                   0.00954453                  0.00151285
PTDF \*\*                                                                         0.000715833                 0.000696339
contingency analysis \*\*\*                                                       0.00806741                  0.0030212
LODF \*\*\*                                                                       0.000335581                 0.000229532
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

From a grid2op perspective, lightsim2grid allows to compute up to ~2851 DC steps each second (column `grid2op speed`, row `DC (NICSLU\*)`) on the case14_sandbox and "only" ~138 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~21** (2851 / 138) in this case.

When compared to powsybl (with the pypowsybl backend), lightsim2grid is around **~4.0** times faster (713 vs 2851).

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2848 vs 2781 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.351ms (`1. / 2851. * 1000.`) on average and on this 351 µs (or 0.351ms), only 1 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~349 µs are spent in the grid2op extra layer or in the backend implementation in this case (`100%` of the computation time - `=349 / 351`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.0012 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0513 ms for the fastest grid2op DC backend (`DC (NICSLU\*)`), a **~43x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00634 ms on average (row `contingency analysis`), a **~8x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 7.16e-05 ms per powerflow, a **~717x** speed up compared to the fastest grid2op DC backend.

Likewise, using the LODF matrix (row `LODF`) to perform the contingency analysis takes 0.000425 ms per contingency, a **~121x** speed up compared to the fastest grid2op DC backend.


For the IEEE case 118:
+++++++++++++++++++++++++

From a grid2op perspective, lightsim2grid allows to compute up to ~2460 DC steps each second (column `grid2op speed`, row `DC (CKTSO\*)`) on the neurips_2020_track2 and "only" ~127 for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), leading to a speed up of **~19** (2460 / 127) in this case.

For this environment there is no sensible difference in using `KLU` linear solver (row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) (2456 vs 2433 iterations on the reported runs, might slightly vary across runs).

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, at least we think the observed differences are within error margins.

For this environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) takes 0.407ms (`1. / 2460. * 1000.`) on average and on this 407 µs (or 0.407ms), only 3 µs are spent in the backend (column `time in 'algo' (ms / pf)`). Meaning that ~403 µs are spent in the grid2op extra layer or in the backend implementation in this case (`99%` of the computation time - `=403 / 407`- is external to the powerflow algorithm)

The `TimeSerie` module performs one DC powerflow in 0.00954 ms on average (row `time serie`, column `grid2op 'backend.runpf' time`), compared to 0.0607 ms for the fastest grid2op DC backend (`DC (CKTSO\*)`), a **~6x** speed up.

Similarly, the `ContingencyAnalysis` module performs one DC contingency in 0.00807 ms on average (row `contingency analysis`), a **~8x** speed up compared to the fastest grid2op DC backend.

Using the PTDF matrix directly (row `PTDF`) is even faster: 0.000716 ms per powerflow, a **~85x** speed up compared to the fastest grid2op DC backend.

Likewise, using the LODF matrix (row `LODF`) to perform the contingency analysis takes 0.000336 ms per contingency, a **~181x** speed up compared to the fastest grid2op DC backend.

See TL;DR section at the top of the file.
