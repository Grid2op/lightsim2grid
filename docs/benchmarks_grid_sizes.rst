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
case14                         14           0.014566               0.0360995            0.00546858                       0.0109907
case118                       118           0.0835847              0.221692             0.038204                         0.0535954
case_illinois200              200           0.167528               0.393047             0.0718412                        0.112911
case300                       300           0.277824               0.637875             0.1498                           0.201345
case1354pegase               1354           1.58167                2.97318              0.909525                         1.14373
case1888rte                  1888           2.46013                4.16834              1.19669                          1.4888
case2848rte                  2848           3.82885                6.47039              1.81877                          2.34153
case2869pegase               2869           3.73823                6.84596              2.1224                           2.49936
case3120sp                   3120           4.26664                7.29163              1.65787                          2.47875
case6495rte                  6495          11.8336                18.459                5.33566                          6.15428
case6515rte                  6515          13.6057                20.0196               5.22459                          6.27699
case9241pegase               9241          17.3575                28.5716               8.76378                          9.95424
================  ===============  ==================  =====================  ====================  ==============================
   

All timings reported above are in milliseconds (ms) for one powerflow (in all cases lots of powerflow are carried out, up to a thousands
and the timings here are averaged accross all the powerflows performed)

For detailed explanation about each column as well as the hardware used, please refer to the section below, but in summary:

- benchmark were run on python 3.12 with a laptop (see section :ref:`bench_grid_size_hardware`
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

- date: 2026-04-21 09:51  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.13.5.final.0 (64 bit)
- numpy version: 2.3.5
- pandas version: 2.3.3
- pandapower version: 3.4.0
- grid2op version: 1.12.4.dev0
- lightsim2grid version: 0.13.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: False 
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
case14                         14                  0.317006                   0.0236869        42217.4                        0.014566                    0.0110552
case118                       118                  0.439654                   0.0977845        10226.6                        0.0835847                   0.0714621
case_illinois200              200                  0.52879                    0.184574          5417.89                       0.167528                    0.152389
case300                       300                  0.687898                   0.300009          3333.23                       0.277824                    0.255286
case1354pegase               1354                  2.44809                    1.64378            608.353                      1.58167                     1.48944
case1888rte                  1888                  3.22291                    2.54091            393.56                       2.46013                     2.34192
case2848rte                  2848                  4.75012                    3.94042            253.78                       3.82885                     3.65101
case2869pegase               2869                  5.20068                    3.86489            258.74                       3.73823                     3.52338
case3120sp                   3120                  5.25079                    4.3943             227.568                      4.26664                     4.08523
case6495rte                  6495                 13.6749                    12.1174              82.5257                    11.8336                     11.3858
case6515rte                  6515                 15.4622                    13.8893              71.9976                    13.6057                     13.1395
case9241pegase               9241                 21.9036                    17.8509              56.0197                    17.3575                     16.5355
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.368089                   0.0590322        16939.9                        0.0360995                   0.0290161
case118                       118                  0.6709                     0.305026          3278.4                        0.221692                    0.194384
case_illinois200              200                  0.87917                    0.509671          1962.05                       0.393047                    0.357242
case300                       300                  1.21275                    0.808751          1236.47                       0.637875                    0.584015
case1354pegase               1354                  4.48071                    3.65709            273.441                      2.97318                     2.71566
case1888rte                  1888                  5.67926                    4.99325            200.27                       4.16834                     3.88297
case2848rte                  2848                  8.57918                    7.73595            129.267                      6.47039                     6.03774
case2869pegase               2869                  9.76003                    8.38388            119.277                      6.84596                     6.29586
case3120sp                   3120                  9.6056                     8.71158            114.79                       7.29163                     6.85212
case6495rte                  6495                 23.2027                    21.5196              46.4693                    18.459                      17.391
case6515rte                  6515                 24.7994                    23.0958              43.2979                    20.0196                     18.9434
case9241pegase               9241                 38.38                      34.1234              29.3054                    28.5716                     26.5461
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
case14                         14        0.00546858        182863
case118                       118        0.038204           26175.2
case_illinois200              200        0.0718412          13919.6
case300                       300        0.1498              6675.56
case1354pegase               1354        0.909525            1099.48
case1888rte                  1888        1.19669              835.64
case2848rte                  2848        1.81877              549.823
case2869pegase               2869        2.1224               471.164
case3120sp                   3120        1.65787              603.184
case6495rte                  6495        5.33566              187.418
case6515rte                  6515        5.22459              191.402
case9241pegase               9241        8.76378              114.106
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
case14                         14            0.0109907            90986.1
case118                       118            0.0535954            18658.3
case_illinois200              200            0.112911              8856.57
case300                       300            0.201345              4966.6
case1354pegase               1354            1.14373                874.333
case1888rte                  1888            1.4888                 671.68
case2848rte                  2848            2.34153                427.072
case2869pegase               2869            2.49936                400.102
case3120sp                   3120            2.47875                403.429
case6495rte                  6495            6.15428                162.488
case6515rte                  6515            6.27699                159.312
case9241pegase               9241            9.95424                100.46
================  ===============  ===================  ===================
