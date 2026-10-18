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
case14                         14           0.0123822              0.0323205            0.00451251                      0.00960011
case118                       118           0.0699949              0.19534              0.0289472                       0.0403928
case_illinois200              200           0.146258               0.34833              0.0620114                       0.0974685
case300                       300           0.269318               0.575865             0.151003                        0.188785
case1354pegase               1354           1.42593                2.63143              0.818186                        0.938555
case1888rte                  1888           2.1658                 3.62574              1.02095                         1.19667
case2848rte                  2848           3.42232                5.70926              1.58737                         1.88313
case2869pegase               2869           3.40024                6.06937              1.88812                         2.08635
case3120sp                   3120           3.96382                6.61183              1.4962                          2.14265
case6495rte                  6495          10.7373                16.3327               4.49756                         5.1188
case6515rte                  6515          12.2724                17.8264               4.52332                         5.09823
case9241pegase               9241          16.0023                25.8985               7.76783                         8.49665
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

- date: 2026-10-12 19:54  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.4.6
- pandas version: 2.3.3
- pandapower version: 3.5.4
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.2.0
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
case14                         14                  0.310173                   0.0214094        46708.5                        0.0123822                  0.00969786
case118                       118                  0.422622                   0.0840015        11904.5                        0.0699949                  0.0630786
case_illinois200              200                  0.516392                   0.163776          6105.9                        0.146258                   0.138146
case300                       300                  0.681292                   0.291091          3435.35                       0.269318                   0.258211
case1354pegase               1354                  2.27568                    1.48391            673.895                      1.42593                    1.38161
case1888rte                  1888                  2.89869                    2.24131            446.168                      2.1658                     2.10899
case2848rte                  2848                  4.31387                    3.52781            283.462                      3.42232                    3.34242
case2869pegase               2869                  4.79877                    3.51743            284.298                      3.40024                    3.30119
case3120sp                   3120                  4.94227                    4.08751            244.648                      3.96382                    3.87784
case6495rte                  6495                 12.4769                    11.0023              90.8897                    10.7373                    10.5486
case6515rte                  6515                 14.1533                    12.5477              79.6956                    12.2724                    12.0579
case9241pegase               9241                 20.4096                    16.4566              60.766                     16.0023                    15.5775
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.363198                   0.0554112        18046.9                        0.0323205                   0.0245302
case118                       118                  0.639341                   0.274908          3637.58                       0.19534                     0.170539
case_illinois200              200                  0.826866                   0.459052          2178.4                        0.34833                     0.319771
case300                       300                  1.13999                    0.737171          1356.54                       0.575865                    0.533276
case1354pegase               1354                  4.07029                    3.25836            306.903                      2.63143                     2.44203
case1888rte                  1888                  5.08718                    4.40123            227.209                      3.62574                     3.41198
case2848rte                  2848                  7.7166                     6.89612            145.009                      5.70926                     5.39332
case2869pegase               2869                  8.8433                     7.52683            132.858                      6.06937                     5.66095
case3120sp                   3120                  8.86395                    7.98812            125.186                      6.61183                     6.29116
case6495rte                  6495                 20.827                     19.2179              52.0347                    16.3327                     15.6061
case6515rte                  6515                 22.4607                    20.7357              48.226                     17.8264                     17.0543
case9241pegase               9241                 35.1783                    31.0138              32.2437                    25.8985                     24.3234
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
case14                         14        0.00451251        221606
case118                       118        0.0289472          34545.6
case_illinois200              200        0.0620114          16126.1
case300                       300        0.151003            6622.4
case1354pegase               1354        0.818186            1222.22
case1888rte                  1888        1.02095              979.481
case2848rte                  2848        1.58737              629.972
case2869pegase               2869        1.88812              529.627
case3120sp                   3120        1.4962               668.362
case6495rte                  6495        4.49756              222.343
case6515rte                  6515        4.52332              221.076
case9241pegase               9241        7.76783              128.736
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
case14                         14           0.00960011           104166
case118                       118           0.0403928             24756.9
case_illinois200              200           0.0974685             10259.7
case300                       300           0.188785               5297.04
case1354pegase               1354           0.938555               1065.47
case1888rte                  1888           1.19667                 835.654
case2848rte                  2848           1.88313                 531.031
case2869pegase               2869           2.08635                 479.306
case3120sp                   3120           2.14265                 466.713
case6495rte                  6495           5.1188                  195.358
case6515rte                  6515           5.09823                 196.147
case9241pegase               9241           8.49665                 117.693
================  ===============  ===================  ===================

Comments
--------

Allowing lightsim2grid to "recycle" previous computation (column `avg step duration (ms)`, default behaviour) instead of restarting from scratch at every step makes grid2op between **~1.2x** (on `case14`) and **~1.9x** (on `case2869pegase`) faster, depending on the grid size.

Compared to a regular grid2op step (with recycling), the `TimeSerie` module is between **~2.6x** (on `case2869pegase`) and **~69.5x** (on `case14`) faster.

Similarly, the `ContingencyAnalysis` module is between **~2.3x** (on `case3120sp`) and **~31.8x** (on `case14`) faster than a regular grid2op step (with recycling) to evaluate one contingency.
