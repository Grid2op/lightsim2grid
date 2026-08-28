.. _benchmark-grid-size:

Benchmarks (grid size)
======================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package in the ``benchmarks`` folder.

If you are interested in other type of benchmarks, let us know !

TL;DR
-------

In summary, lightsim2grid (when using KLU linear solver) perfomances are:

================  ===============  ==================  =====================  ====================  ==============================
grid                size (nb bus)    time (recycling)    time (no recycling)    time (`TimeSerie`)    time (`ContingencyAnalysis`)
================  ===============  ==================  =====================  ====================  ==============================
case14                         14           0.0180543              0.0393453            0.00691772                       0.0141547
case118                       118           0.0925261              0.245769             0.034242                         0.0501795
case_illinois200              200           0.163772               0.418551             0.06606                          0.112861
case300                       300           0.314592               0.688923             0.183331                         0.235984
case1354pegase               1354           1.63302                3.10312              0.883601                         1.09353
case1888rte                  1888           2.44727                4.32304              1.11172                          1.40246
case2848rte                  2848           3.8338                 6.75225              1.7368                           2.20093
case2869pegase               2869           3.82374                7.33523              2.05757                          2.41441
case3120sp                   3120           4.30735                7.57317              1.61617                          2.42757
case6495rte                  6495          11.8292                18.8783               4.82027                          5.84038
case6515rte                  6515          13.4284                20.365                4.96732                          5.87362
case9241pegase               9241          17.6527                29.9001               8.38969                          9.64497
================  ===============  ==================  =====================  ====================  ==============================
   

All timings reported above are in milliseconds (ms) for one powerflow (in all cases lots of powerflow are carried out, up to a thousands
and the timings here are averaged accross all the powerflows performed)

For detailed explanation about each column as well as the hardware used, please refer to the section below, but in summary:

- benchmark were run on python 3.12 with a laptop (see section :ref:`bench_grid_size_hardware`
  and page :ref:`benchmark-deep-dive` for more information about the exact definition of the timers ):
- `time (recycling)` indicates the average time it took to run 1 powerflow (with consecutive run of 288 powerflows)
  while allowing lighsim2grid to re use some basic previous computation from one powerflow to another. This is the most common
  usecase in grid2op for example (default behaviour). See :ref:`bench_grid_size_glop` for more information
- `time (no recycling)` indicates the same average time as aboved but lightsim2grid is forced to restart the 
  computation from scratch each time, as if it was a completely different grid on a completely different computers. 
  See :ref:`bench_grid_size_glop` for more information.
- `time (TimeSerie)` reports the time it takes to run one powerflow using the lightsim2grid `TimeSerie` module, were 
  everything is in c++ and some care has been taken to improve the performance (reuse of as many things as possible, 
  carefull memory allocation, etc.). See :ref:`bench_grid_size_ts` for more information.
- `time (ContingencyAnalysis)` reports the time it takes to run one powerflow using the lightsim2grid `ContingencyAnalysis` module, were
  everything is in c++ and some care has been taken to improve the performance (reuse of as many things as possible,
  carefull memory allocation, etc.). See :ref:`bench_grid_size_ca` for more information. **NB** on this settings,
  as opposed to the others, the grid production / generations stay the same, but the grid topology changes by the
  connection and disconnection of powerlines.

.. note::
  Unlike the other benchmark pages, the TL;DR table above is not hand-copied from the detailed tables further
  down this page: ``benchmark_grid_size.py`` prints it directly, computed from the exact same run that produces
  the 4 detailed tables below, so the two cannot drift apart from each other. The script also now prints a
  "Description" paragraph (see the "Comments" section at the end of this page) commenting on these numbers,
  generated the same way as on the other benchmark pages.

.. _bench_grid_size_hardware:

Using a grid2op environment
----------------------------
In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
on different grid sizes varying from the ieee case 14 grid (14 buses) up to the pegase 9241 grid (case9241 from pandapower
counting 9241 buses).

All of them has been run on a computer with a the following characteristics:

- date: 2026-08-06 18:23  CEST
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


Solver used for linear algebra: NR single (KLU)


To run the benchmark, ``cd`` into the ``benchmarks`` folder and type:

.. code-block:: bash

    python benchmark_grid_size.py

(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

.. _bench_grid_size_glop:

Computation time using grid2op
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This benchmark in doing repeat calls to `env.step(do_nothing)` (usually 288 or 1000) for a given environment build 
on a grid coming from data available in pandapower.

Then we compare different measurments:

- `avg step duration (ms)` is the average time it takes to perform the `grid2op.step`. It is given in milliseconds (ms).
  It takes into account the time to read the data, to feed the data to the underlying c++ model, to run the powerflow
  and to read back the data from the c++ model.
- `time [DC + AC] (ms / pf)` is the time it takes to perform the entire powerflow, which consists in first 
  providing an initial guess (DC approximation) and then to compute the powerflow. As compared to the 
  above timings, it only take into account the time to run the powerflow. This "time to run the powerflow" 
  can be at this stage decomposed in:

  - converting the provided data into valid matrix / vector to run a DC powerflow
  - computing a DC powerflow (used to initialize the AC powerflow)
  - converting again the provided data into valid matrix / vector to run an AC powerflow
  - computint the AC Powerflow
  - post processing the internal data (which includes *eg* the flows on the lines in amps, the reactive value
    produced / absorbed by each generator etc.)

- `time in 'solver' (ms / pf)` gives the time it takes to only perform the AC powerflow:

  - converting the provided data into valid matrix / vector to run an AC powerflow
  - computing the AC Powerflow
  - post processing the internal data (which includes *eg* the flows on the lines in amps, the reactive value
    produced / absorbed by each generator etc.)
    
- `time in 'algo' (ms / pf)` gives the time spent in the algorithm that computes the AC powerflow only

.. warning::
  For more information about what is actually done and the wordings used in this section, 
  you can consult the page :ref:`benchmark-deep-dive`
  
The results are given in two tables:

- the first one corresponds to the default settings were lightsim2grid is allowed to "recycle" previous
  results, which is the default in grid2op and lightsim2grid. This corresponds to a generic grid2op usecase.
- the second one is the same run for the same environment, but this time lightsim2grid recreate everything from
  scratch each time, the "recycling" is deactivated.

The main impact on "recycling" is that, when activated (default), lightsim2grid can skip some of its internal 
computation, especially in the steps:

- "converting the provided data into valid matrix / vector to run a DC powerflow"
- "converting again the provided data into valid matrix / vector to run an AC powerflow"
- also the computation of the DC and AC powerflows can be a little bit faster (depending on the linear solver used)

The "no recycling" strategy is closer to a situation were you would simulate different powerflows on 
different cores or even  on different computers and cannot share the internal state of the solvers (for example). 
It can also represent a situation were you would run powerflows for vastly different grids one after 
the other.


Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (recyling allowed, default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid                size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.416912                   0.0294001        34013.5                        0.0180543                   0.0141478
case118                       118                  0.456141                   0.107554          9297.69                       0.0925261                   0.0824107
case_illinois200              200                  0.529129                   0.180569          5538.04                       0.163772                    0.151446
case300                       300                  0.722134                   0.335469          2980.9                        0.314592                    0.296417
case1354pegase               1354                  2.49597                    1.68992            591.744                      1.63302                     1.55694
case1888rte                  1888                  3.16992                    2.5168             397.329                      2.44727                     2.35697
case2848rte                  2848                  4.74195                    3.93506            254.126                      3.8338                      3.69879
case2869pegase               2869                  5.24078                    3.93682            254.012                      3.82374                     3.65829
case3120sp                   3120                  5.283                      4.42263            226.11                       4.30735                     4.16654
case6495rte                  6495                 13.6542                    12.0822              82.7664                    11.8292                     11.4815
case6515rte                  6515                 15.2783                    13.6811              73.0938                    13.4284                     13.0776
case9241pegase               9241                 22.0619                    18.0824              55.3025                    17.6527                     16.9855
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.382921                   0.0647661        15440.2                        0.0393453                   0.0315124
case118                       118                  0.716991                   0.335443          2981.13                       0.245769                    0.218422
case_illinois200              200                  1.11168                    0.538909          1855.6                        0.418551                    0.385138
case300                       300                  1.27529                    0.8701            1149.29                       0.688923                    0.635676
case1354pegase               1354                  4.63179                    3.81146            262.367                      3.10312                     2.8725
case1888rte                  1888                  5.87786                    5.20158            192.249                      4.32304                     4.06026
case2848rte                  2848                  8.95185                    8.10464            123.386                      6.75225                     6.35489
case2869pegase               2869                 10.4034                     9.0233             110.824                      7.33523                     6.81353
case3120sp                   3120                  9.96229                    9.06944            110.26                       7.57317                     7.18675
case6495rte                  6495                 23.8855                    22.1442              45.1586                    18.8783                     17.904
case6515rte                  6515                 25.3544                    23.6233              42.3312                    20.365                      19.3901
case9241pegase               9241                 39.9119                    35.7023              28.0094                    29.9001                     27.9824
================  ===============  ========================  ==========================  ================  ============================  ==========================

.. _bench_grid_size_ts:

Computation time using the lightsim2grid `TimeSerie` module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As opposed to the experiment above, the `TimeSerie` lightsim2grid module allows to perform sequential computation
of varying productions and loads with the exact same grid topology.

This does not rely on grid2op and is coded in "pure c++" still using one single CPU core. It should be faster than 
the timings reported on the above sequence because:

- the loop is made in c++ instead of python
- the code has been optimize to run faster and "recycle" as many things as possible: the 
  matrices representing the grid is computed only once, it is factorized only once, 
  conversion from the internal solver representation to MW, MVAr and A is done in 
  a vectorized way etc.

This rapidity has a cost, it is much less flexible. With the grid2op framework an "agent"
can do a lot of different actions (even though "do nothing" was used for the benchmark). Here
on the other hand, only a "*do nothing*" action can be performed (and without emulation of
any kind of protections).

The column `time (ms / pf)` can be compared with the column `time [DC + AC] (ms / pf)` of the 
table in the previous benchmark.

================  ===============  ================  ================
grid                size (nb bus)    time (ms / pf)    speed (pf / s)
================  ===============  ================  ================
case14                         14        0.00691772        144556
case118                       118        0.034242           29203.9
case_illinois200              200        0.06606            15137.7
case300                       300        0.183331            5454.61
case1354pegase               1354        0.883601            1131.73
case1888rte                  1888        1.11172              899.509
case2848rte                  2848        1.7368               575.772
case2869pegase               2869        2.05757              486.01
case3120sp                   3120        1.61617              618.746
case6495rte                  6495        4.82027              207.457
case6515rte                  6515        4.96732              201.316
case9241pegase               9241        8.38969              119.194
================  ===============  ================  ================

.. _bench_grid_size_ca:

Computation time using the lightsim2grid `ContingencyAnalysis` module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As opposed to the benchmarks reported in the previous two sections, this benchmark 
is focused on the `ContingencyAnalysis` lightsim2grid module.

A "contingency analysis" is often carried out in power system. The objective is
to assess whether or not the current grid state is safe if one (or more)
powerline would be disconnected. It uses the same 
productions / consumptions for each computation. Each time it disconnects
one or more powerlines, run the powerflow and then stores the results.

For this benchmark we focus on disconnecting only one powerline (though 
lightsim2grid offers the possibility to disconnect as many as you want) with 
a limit on 1000 contingency simulated (even for grid were there would be 
more than 1000 powerlines / trafos to disconnect we limit the computation to 
only 1000).

================  ===============  ===================  ===================
grid                size (nb bus)    time (ms / cont.)    speed (cont. / s)
================  ===============  ===================  ===================
case14                         14            0.0141547            70648
case118                       118            0.0501795            19928.4
case_illinois200              200            0.112861              8860.47
case300                       300            0.235984              4237.57
case1354pegase               1354            1.09353                914.474
case1888rte                  1888            1.40246                713.035
case2848rte                  2848            2.20093                454.354
case2869pegase               2869            2.41441                414.18
case3120sp                   3120            2.42757                411.934
case6495rte                  6495            5.84038                171.222
case6515rte                  6515            5.87362                170.253
case9241pegase               9241            9.64497                103.681
================  ===============  ===================  ===================

Comments
--------

This is the text printed by ``benchmark_grid_size.py`` (see the note above the TL;DR table) for the tables
above, computed from the numbers actually measured during that run.

Allowing lightsim2grid to "recycle" previous computation (column `avg step duration (ms)`, default behaviour) instead of restarting from scratch at every step makes grid2op between **~0.9x** (on `case14`) and **~2.1x** (on `case_illinois200`) faster, depending on the grid size.

Compared to a regular grid2op step (with recycling), the `TimeSerie` module is between **~2.5x** (on `case2869pegase`) and **~60.3x** (on `case14`) faster.

Similarly, the `ContingencyAnalysis` module is between **~2.2x** (on `case2848rte`) and **~29.5x** (on `case14`) faster than a regular grid2op step (with recycling) to evaluate one contingency.
