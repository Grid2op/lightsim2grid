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
case14                         14           0.0151524              0.0393676            0.00546965                       0.0120966
case118                       118           0.0907367              0.239201             0.0381798                        0.055575
case_illinois200              200           0.163444               0.423087             0.0683476                        0.112028
case300                       300           0.311167               0.704807             0.176352                         0.229164
case1354pegase               1354           1.59221                3.11083              0.900187                         1.11285
case1888rte                  1888           2.41214                4.32001              1.12627                          1.42144
case2848rte                  2848           3.7678                 6.71873              1.74596                          2.21609
case2869pegase               2869           3.78651                7.26164              2.08606                          2.42732
case3120sp                   3120           4.26542                7.49865              1.6428                           2.48339
case6495rte                  6495          11.5761                18.7157               4.78278                          5.74913
case6515rte                  6515          12.9611                20.1715               4.91258                          5.81299
case9241pegase               9241          17.4095                29.7788               8.40066                          9.6462
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

- date: 2026-08-28 16:56  CEST
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
case14                         14                  0.320099                   0.0244134        40961.2                        0.0151524                   0.0117544
case118                       118                  0.451427                   0.105576          9471.88                       0.0907367                   0.080653
case_illinois200              200                  0.528523                   0.180514          5539.72                       0.163444                    0.150956
case300                       300                  0.725959                   0.333054          3002.51                       0.311167                    0.292701
case1354pegase               1354                  2.45353                    1.65307            604.936                      1.59221                     1.51582
case1888rte                  1888                  3.15837                    2.48821            401.895                      2.41214                     2.32146
case2848rte                  2848                  4.70391                    3.87832            257.843                      3.7678                      3.63089
case2869pegase               2869                  5.2151                     3.90676            255.967                      3.78651                     3.61854
case3120sp                   3120                  5.30337                    4.39481            227.541                      4.26542                     4.11986
case6495rte                  6495                 13.403                     11.845               84.4239                    11.5761                     11.2339
case6515rte                  6515                 14.7969                    13.225               75.6144                    12.9611                     12.6123
case9241pegase               9241                 21.7958                    17.866               55.9721                    17.4095                     16.7472
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.375508                   0.0636576        15709.1                        0.0393676                   0.0312537
case118                       118                  0.701734                   0.327901          3049.7                        0.239201                    0.211855
case_illinois200              200                  0.920673                   0.546093          1831.19                       0.423087                    0.387856
case300                       300                  1.29638                    0.888698          1125.24                       0.704807                    0.650414
case1354pegase               1354                  4.64983                    3.82789            261.24                       3.11083                     2.87081
case1888rte                  1888                  5.89917                    5.20745            192.033                      4.32001                     4.04502
case2848rte                  2848                  8.93703                    8.088              123.64                       6.71873                     6.3059
case2869pegase               2869                 10.2971                     8.9443             111.803                      7.26164                     6.71196
case3120sp                   3120                  9.92825                    9.01849            110.883                      7.49865                     7.09966
case6495rte                  6495                 23.9327                    22.0092              45.4356                    18.7157                     17.7245
case6515rte                  6515                 25.1584                    23.4604              42.625                     20.1715                     19.1766
case9241pegase               9241                 39.761                     35.6363              28.0612                    29.7788                     27.8505
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
case14                         14        0.00546965        182827
case118                       118        0.0381798          26191.9
case_illinois200              200        0.0683476          14631.1
case300                       300        0.176352            5670.46
case1354pegase               1354        0.900187            1110.88
case1888rte                  1888        1.12627              887.888
case2848rte                  2848        1.74596              572.75
case2869pegase               2869        2.08606              479.372
case3120sp                   3120        1.6428               608.716
case6495rte                  6495        4.78278              209.083
case6515rte                  6515        4.91258              203.559
case9241pegase               9241        8.40066              119.038
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
case14                         14            0.0120966            82668
case118                       118            0.055575             17993.7
case_illinois200              200            0.112028              8926.3
case300                       300            0.229164              4363.69
case1354pegase               1354            1.11285                898.597
case1888rte                  1888            1.42144                703.511
case2848rte                  2848            2.21609                451.244
case2869pegase               2869            2.42732                411.976
case3120sp                   3120            2.48339                402.676
case6495rte                  6495            5.74913                173.939
case6515rte                  6515            5.81299                172.029
case9241pegase               9241            9.6462                 103.668
================  ===============  ===================  ===================

Comments
--------

This is the text printed by ``benchmark_grid_size.py`` (see the note above the TL;DR table) for the tables
above, computed from the numbers actually measured during that run.
Allowing lightsim2grid to "recycle" previous computation (column `avg step duration (ms)`, default behaviour) instead of restarting from scratch at every step makes grid2op between **~1.2x** (on `case14`) and **~2.0x** (on `case2869pegase`) faster, depending on the grid size.

Compared to a regular grid2op step (with recycling), the `TimeSerie` module is between **~2.5x** (on `case2869pegase`) and **~58.5x** (on `case14`) faster.

Similarly, the `ContingencyAnalysis` module is between **~2.1x** (on `case2848rte`) and **~26.5x** (on `case14`) faster than a regular grid2op step (with recycling) to evaluate one contingency.
