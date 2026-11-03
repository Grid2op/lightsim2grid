.. _benchmark-grid-size:

Benchmarks (grid size)
======================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package in the ``benchmarks`` folder.

If you are interested in other type of benchmarks, let us know !

TL;DR
-------

In summary, lightsim2grid (when using **NR single (KLU static pivot)** linear solver) perfomances are:

================  ===============  ==================  =====================  ====================  ==============================
grid                size (nb bus)    time (recycling)    time (no recycling)    time (`TimeSerie`)    time (`ContingencyAnalysis`)
================  ===============  ==================  =====================  ====================  ==============================
case14                         14           0.0104296              0.0325648            0.00345285                      0.00836953
case118                       118           0.0464082              0.183256             0.0199597                       0.0263036
case_illinois200              200           0.0945836              0.324233             0.0411208                       0.0659203
case300                       300           0.163388               0.531766             0.0816072                       0.0999404
case1354pegase               1354           0.971325               2.37345              0.534323                        0.604285
case1888rte                  1888           1.37489                3.16793              0.651767                        0.745152
case2848rte                  2848           2.23207                4.86451              1.05791                         1.22029
case2869pegase               2869           2.32699                5.28641              1.31411                         1.41216
case3120sp                   3120           2.78034                5.94121              1.03495                         1.47456
case6495rte                  6495           7.46306               13.9053               3.02447                         3.49925
case6515rte                  6515           8.13535               14.4648               3.09987                         3.49521
case9241pegase               9241          11.1062                22.0007               5.46248                         5.86093
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

- date: 2026-11-02 15:46  CET
- system: Linux 6.8.0-111-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.4.6
- pandas version: 2.3.3
- pandapower version: 3.5.4
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.3.0
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
case14                         14                  0.313829                   0.0201867         49537.6                       0.0104296                  0.00721816
case118                       118                  0.409336                   0.0616762         16213.7                       0.0464082                  0.0382189
case_illinois200              200                  0.464154                   0.112736           8870.26                      0.0945836                  0.0854076
case300                       300                  0.574369                   0.185849           5380.72                      0.163388                   0.150938
case1354pegase               1354                  1.90407                    1.04146             960.19                      0.971325                   0.919179
case1888rte                  1888                  2.14263                    1.45727             686.216                     1.37489                    1.31099
case2848rte                  2848                  3.18645                    2.34842             425.818                     2.23207                    2.14288
case2869pegase               2869                  3.80467                    2.45582             407.197                     2.32699                    2.21592
case3120sp                   3120                  4.0576                     2.94449             339.618                     2.78034                    2.6715
case6495rte                  6495                  9.5597                     7.77764             128.574                     7.46306                    7.23404
case6515rte                  6515                 10.1586                     8.43499             118.554                     8.13535                    7.9111
case9241pegase               9241                 15.4995                    11.577                86.378                    11.1062                    10.6844
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.369143                   0.0569078        17572.3                        0.0325648                   0.0236052
case118                       118                  0.639279                   0.265991          3759.52                       0.183256                    0.155578
case_illinois200              200                  0.809087                   0.437373          2286.38                       0.324233                    0.29225
case300                       300                  1.11389                    0.696982          1434.76                       0.531766                    0.487273
case1354pegase               1354                  4.01283                    3.05624            327.199                      2.37345                     2.15298
case1888rte                  1888                  4.7751                     3.99719            250.175                      3.16793                     2.91133
case2848rte                  2848                  7.0342                     6.11279            163.592                      4.86451                     4.5043
case2869pegase               2869                  8.19752                    6.77525            147.596                      5.28641                     4.84499
case3120sp                   3120                  8.66598                    7.42988            134.592                      5.94121                     5.57504
case6495rte                  6495                 18.8975                    16.9526              58.9881                    13.9053                     13.0838
case6515rte                  6515                 19.3087                    17.4531              57.2964                    14.4648                     13.6614
case9241pegase               9241                 31.2481                    27.1468              36.8368                    22.0007                     20.3864
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
case14                         14        0.00345285        289616
case118                       118        0.0199597          50101.1
case_illinois200              200        0.0411208          24318.6
case300                       300        0.0816072          12253.8
case1354pegase               1354        0.534323            1871.53
case1888rte                  1888        0.651767            1534.29
case2848rte                  2848        1.05791              945.264
case2869pegase               2869        1.31411              760.97
case3120sp                   3120        1.03495              966.233
case6495rte                  6495        3.02447              330.636
case6515rte                  6515        3.09987              322.594
case9241pegase               9241        5.46248              183.067
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
case14                         14           0.00836953           119481
case118                       118           0.0263036             38017.7
case_illinois200              200           0.0659203             15169.8
case300                       300           0.0999404             10006
case1354pegase               1354           0.604285               1654.85
case1888rte                  1888           0.745152               1342.01
case2848rte                  2848           1.22029                 819.479
case2869pegase               2869           1.41216                 708.133
case3120sp                   3120           1.47456                 678.167
case6495rte                  6495           3.49925                 285.775
case6515rte                  6515           3.49521                 286.106
case9241pegase               9241           5.86093                 170.621
================  ===============  ===================  ===================

Comments
--------

Allowing lightsim2grid to "recycle" previous computation (column `avg step duration (ms)`, default behaviour) instead of restarting from scratch at every step makes grid2op between **~1.2x** (on `case14`) and **~2.2x** (on `case1888rte`) faster, depending on the grid size.

Compared to a regular grid2op step (with recycling), the `TimeSerie` module is between **~2.8x** (on `case9241pegase`) and **~90.9x** (on `case14`) faster.

Similarly, the `ContingencyAnalysis` module is between **~2.6x** (on `case2848rte`) and **~37.5x** (on `case14`) faster than a regular grid2op step (with recycling) to evaluate one contingency.
