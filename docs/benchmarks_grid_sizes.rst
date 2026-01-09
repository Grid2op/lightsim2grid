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
case14                         14           0.0147596              0.0362591            0.00550065                       0.0107742
case118                       118           0.0846219              0.227989             0.0378691                        0.0545186
case_illinois200              200           0.162386               0.401702             0.0709988                        0.115073
case300                       300           0.29452                0.659511             0.168504                         0.218641
case1354pegase               1354           1.59691                3.0413               0.914828                         1.14361
case1888rte                  1888           2.47497                4.31092              1.18708                          1.48549
case2848rte                  2848           3.86355                6.59429              1.81397                          2.33079
case2869pegase               2869           3.78899                6.98899              2.16268                          2.52204
case3120sp                   3120           4.2851                 7.31591              1.66029                          2.47509
case6495rte                  6495          11.9117                18.7309               5.30932                          6.15445
case6515rte                  6515          13.4341                20.1974               5.21714                          6.25699
case9241pegase               9241          17.3358                28.8583               8.75372                          9.89653
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
case14                         14                  0.309347                   0.0242351        41262.4                        0.0147596                   0.0113769
case118                       118                  0.443979                   0.103476          9664.06                       0.0846219                   0.0725103
case_illinois200              200                  0.533608                   0.188029          5318.33                       0.162386                    0.14703
case300                       300                  0.710634                   0.328086          3047.98                       0.29452                     0.271714
case1354pegase               1354                  2.50252                    1.70885            585.19                       1.59691                     1.50292
case1888rte                  1888                  3.27753                    2.61607            382.253                      2.47497                     2.35618
case2848rte                  2848                  4.87324                    4.071              245.64                       3.86355                     3.68091
case2869pegase               2869                  5.31515                    4.03376            247.907                      3.78899                     3.57292
case3120sp                   3120                  5.37515                    4.52247            221.118                      4.2851                      4.10012
case6495rte                  6495                 14.0517                    12.4536              80.2983                    11.9117                     11.4436
case6515rte                  6515                 15.5669                    13.9739              71.5619                    13.4341                     12.9532
case9241pegase               9241                 22.2166                    18.2666              54.7448                    17.3358                     16.4956
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.359185                   0.0586501        17050.3                        0.0362591                   0.0294048
case118                       118                  0.676157                   0.311162          3213.76                       0.227989                    0.20073
case_illinois200              200                  0.884864                   0.518262          1929.53                       0.401702                    0.365925
case300                       300                  1.23709                    0.834281          1198.64                       0.659511                    0.604085
case1354pegase               1354                  4.53825                    3.73396            267.812                      3.0413                      2.78292
case1888rte                  1888                  5.84507                    5.16601            193.573                      4.31092                     4.01593
case2848rte                  2848                  8.70574                    7.87036            127.059                      6.59429                     6.16275
case2869pegase               2869                  9.88912                    8.55548            116.884                      6.98899                     6.43026
case3120sp                   3120                  9.61628                    8.73533            114.478                      7.31591                     6.89268
case6495rte                  6495                 23.5506                    21.8496              45.7675                    18.7309                     17.6645
case6515rte                  6515                 24.9907                    23.2837              42.9486                    20.1974                     19.1217
case9241pegase               9241                 38.539                     34.4091              29.062                     28.8583                     26.8374
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
case14                         14        0.00550065        181797
case118                       118        0.0378691          26406.8
case_illinois200              200        0.0709988          14084.7
case300                       300        0.168504            5934.57
case1354pegase               1354        0.914828            1093.1
case1888rte                  1888        1.18708              842.403
case2848rte                  2848        1.81397              551.278
case2869pegase               2869        2.16268              462.39
case3120sp                   3120        1.66029              602.305
case6495rte                  6495        5.30932              188.348
case6515rte                  6515        5.21714              191.676
case9241pegase               9241        8.75372              114.237
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
case14                         14            0.0107742            92814.2
case118                       118            0.0545186            18342.4
case_illinois200              200            0.115073              8690.1
case300                       300            0.218641              4573.7
case1354pegase               1354            1.14361                874.428
case1888rte                  1888            1.48549                673.177
case2848rte                  2848            2.33079                429.039
case2869pegase               2869            2.52204                396.504
case3120sp                   3120            2.47509                404.026
case6495rte                  6495            6.15445                162.484
case6515rte                  6515            6.25699                159.821
case9241pegase               9241            9.89653                101.046
================  ===============  ===================  ===================
