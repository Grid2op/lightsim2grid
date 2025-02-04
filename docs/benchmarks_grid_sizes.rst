.. _benchmark-grid-size:

Benchmarks (grid size)
======================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

TODO DOC in progress

If you are interested in other type of benchmarks, let us know !

TL;DR
-------

In summary, lightsim2grid (when using KLU linear solver) perfomances are:

================  ===============  ==================  =====================  ====================  ==============================
grid                size (nb bus)    time (recycling)    time (no recycling)    time (`TimeSerie`)    time (`ContingencyAnalysis`)
================  ===============  ==================  =====================  ====================  ==============================
case14                         14           0.0130073              0.0343035            0.00457648                      0.00986142
case118                       118           0.0727236              0.209858             0.031906                        0.0466558
case_illinois200              200           0.151148               0.363769             0.0633617                       0.102013
case300                       300           0.264309               0.59775              0.129651                        0.174466
case1354pegase               1354           1.50257                2.7327               0.826231                        1.05762
case1888rte                  1888           2.37322                3.92464              1.06475                         1.37497
case2848rte                  2848           3.7093                 6.13028              1.62492                         2.19232
case2869pegase               2869           3.6351                 6.42046              1.92647                         2.3468
case3120sp                   3120           4.13678                6.85112              1.52145                         2.32498
case6495rte                  6495          11.4654                17.3329               4.96104                         5.75883
case6515rte                  6515          13.1227                18.832                4.77071                         5.86091
case9241pegase               9241          16.9394                27.053                8.11946                         9.34644
================  ===============  ==================  =====================  ====================  ==============================
   

All timings reported above are in milliseconds (ms) for one powerflow (in all cases lots of powerflow are carried out, up to a thousands
and the timings here are averaged accross all the powerflows performed)

For detailed explanation about each column as well as the hardware used, please refer to the section below, but in summary:

- benchmark were run on python 3.12 with an old laptop (i7-6820HQ CPU @ 2.70GHz from 2015) (see section :ref:`bench_grid_size_hardware`
  and page :ref:`benchmark-deep-dive` for more information about the exact definition of the timers ):
- `time (recycling)` indicates the average time it took to run 1 powerflow (with consecutive run of 288 powerflows)
  while allowing lighsim2grid to re use some basic previous computation from one powerflow to another. This is the most consommations
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

.. _bench_grid_size_hardware:

Using a grid2op environment
----------------------------
In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
on different grid sizes varying from the ieee case 14 grid (14 buses) up to the pegase 9241 grid (case9241 from pandapower
counting 9241 buses).

All of them has been run on a computer with a the following characteristics:

- date: 2025-01-09 11:59  CET
- system: Linux 6.5.0-1024-oem
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.9.21.final.0 (64 bit)
- numpy version: 1.26.4
- pandas version: 2.2.3
- pandapower version: 2.14.10
- grid2op version: 1.10.5
- lightsim2grid version: 0.10.0
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 


Solver used for linear algebra: NR single (KLU)


To run the benchmark `cd` in the [benchmark](./benchmarks) folder and type:

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
case14                         14                  0.350895                   0.0245781        40686.6                        0.0130073                   0.0105687
case118                       118                  0.438527                   0.0921714        10849.4                        0.0727236                   0.0640808
case_illinois200              200                  0.531842                   0.177022          5649                          0.151148                    0.139818
case300                       300                  0.692534                   0.298294          3352.4                        0.264309                    0.247054
case1354pegase               1354                  2.54428                    1.61281            620.037                      1.50257                     1.42742
case1888rte                  1888                  3.1374                     2.50807            398.713                      2.37322                     2.27984
case2848rte                  2848                  4.66414                    3.90836            255.862                      3.7093                      3.56542
case2869pegase               2869                  5.45635                    3.87341            258.171                      3.6351                      3.4594
case3120sp                   3120                  5.16431                    4.37043            228.81                       4.13678                     3.99066
case6495rte                  6495                 13.3672                    11.9835              83.4479                    11.4654                     11.1138
case6515rte                  6515                 15.0186                    13.6416              73.305                     13.1227                     12.7565
case9241pegase               9241                 22.7308                    17.9356              55.755                     16.9394                     16.294
================  ===============  ========================  ==========================  ================  ============================  ==========================


Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.394679                   0.0589122        16974.4                        0.0343035                    0.028186
case118                       118                  0.66635                    0.292747          3415.92                       0.209858                     0.187849
case_illinois200              200                  0.851082                   0.476794          2097.34                       0.363769                     0.336049
case300                       300                  1.17444                    0.764839          1307.46                       0.59775                      0.554902
case1354pegase               1354                  4.3213                     3.37901            295.945                      2.7327                       2.52633
case1888rte                  1888                  5.38228                    4.7376             211.077                      3.92464                      3.68506
case2848rte                  2848                  8.18769                    7.40336            135.074                      6.13028                      5.75152
case2869pegase               2869                  9.52221                    7.92512            126.181                      6.42046                      5.94842
case3120sp                   3120                  9.07648                    8.25089            121.199                      6.85112                      6.4863
case6495rte                  6495                 21.8053                    20.3641              49.1061                    17.3329                      16.4422
case6515rte                  6515                 23.2821                    21.8367              45.7945                    18.832                       17.9478
case9241pegase               9241                 37.3876                    32.5509              30.7211                    27.053                       25.3161
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
case14                         14        0.00457648        218508
case118                       118        0.031906           31342.1
case_illinois200              200        0.0633617          15782.4
case300                       300        0.129651            7713.03
case1354pegase               1354        0.826231            1210.32
case1888rte                  1888        1.06475              939.19
case2848rte                  2848        1.62492              615.415
case2869pegase               2869        1.92647              519.085
case3120sp                   3120        1.52145              657.267
case6495rte                  6495        4.96104              201.571
case6515rte                  6515        4.77071              209.613
case9241pegase               9241        8.11946              123.161
================  ===============  ================  ================


.. _bench_grid_size_ca:

Computation time using the lightsim2grid `ContingencyAnalysis` module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As opposed to the benchmarks reported in the previous two sections, this benchmark 
is focused on the `ContingencyAnalysis` lightsim2grid module.

A "contingency analysis" is often carried out in power system. The objective is
to assess whether or not the current grid state is safe if one (or more)
powerline would be disconnected. It uses the same 
productions / consommations for each computation. Each time it disconnects
one or more powerlines, run the powerflow and then stores the results.

For this benchmark we focus on disconnecting only one powerline (though 
lightsim2grid offers the possibility to disconnect as many as you want) with 
a limit on 1000 contingency simulated (even for grid were there would be 
more than 1000 powerlines / trafos to disconnect we limit the computation to 
only 1000).

================  ===============  ===================  ===================
grid                size (nb bus)    time (ms / cont.)    speed (cont. / s)
================  ===============  ===================  ===================
case14                         14           0.00986142           101405
case118                       118           0.0466558             21433.6
case_illinois200              200           0.102013               9802.67
case300                       300           0.174466               5731.77
case1354pegase               1354           1.05762                 945.523
case1888rte                  1888           1.37497                 727.287
case2848rte                  2848           2.19232                 456.137
case2869pegase               2869           2.3468                  426.112
case3120sp                   3120           2.32498                 430.111
case6495rte                  6495           5.75883                 173.646
case6515rte                  6515           5.86091                 170.622
case9241pegase               9241           9.34644                 106.993
================  ===============  ===================  ===================


