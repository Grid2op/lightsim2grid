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
case14                         14           0.01213                0.0321324            0.00441931                      0.00964653
case118                       118           0.0660303              0.190708             0.0289332                       0.0397558
case_illinois200              200           0.129221               0.345833             0.0571606                       0.0883376
case300                       300           0.2383                 0.561618             0.135418                        0.164824
case1354pegase               1354           1.35702                2.59791              0.770959                        0.890424
case1888rte                  1888           2.01847                3.57255              0.965818                        1.11261
case2848rte                  2848           3.21817                5.57038              1.49533                         1.75219
case2869pegase               2869           3.19276                5.98939              1.79082                         1.95654
case3120sp                   3120           3.73936                6.33616              1.4179                          2.03015
case6495rte                  6495          10.1448                15.8663               4.26454                         4.83098
case6515rte                  6515          11.4255                17.1227               4.28735                         4.83759
case9241pegase               9241          15.0171                24.9988               7.34957                         8.05451
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

- date: 2026-09-21 10:16  CEST
- system: Linux 6.8.0-60-generic
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.12.8.final.0 (64 bit)
- numpy version: 2.4.6
- pandas version: 2.3.3
- pandapower version: 3.5.4
- grid2op version: 1.12.5.dev0
- lightsim2grid version: 1.1.0
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
case14                         14                  0.307161                   0.0212539        47050.2                        0.01213                    0.00977312
case118                       118                  0.416393                   0.0802482        12461.3                        0.0660303                  0.0593461
case_illinois200              200                  0.486859                   0.146183          6840.74                       0.129221                   0.121794
case300                       300                  0.644224                   0.259927          3847.23                       0.2383                     0.227654
case1354pegase               1354                  2.19068                    1.41384            707.293                      1.35702                    1.31353
case1888rte                  1888                  2.72162                    2.08776            478.983                      2.01847                    1.96435
case2848rte                  2848                  4.10201                    3.32159            301.061                      3.21817                    3.13888
case2869pegase               2869                  4.57518                    3.30607            302.474                      3.19276                    3.09604
case3120sp                   3120                  4.69402                    3.85867            259.157                      3.73936                    3.65639
case6495rte                  6495                 11.8413                    10.4046              96.1115                    10.1448                     9.96245
case6515rte                  6515                 13.1832                    11.6824              85.5987                    11.4255                    11.2268
case9241pegase               9241                 19.3348                    15.4619              64.6753                    15.0171                    14.604
================  ===============  ========================  ==========================  ================  ============================  ==========================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ============================  ==========================
grid name           size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'solver' (ms / pf)    time in 'algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ============================  ==========================
case14                         14                  0.354886                   0.0548783        18222.2                        0.0321324                   0.0247829
case118                       118                  0.633091                   0.269854          3705.7                        0.190708                    0.166857
case_illinois200              200                  0.820722                   0.457051          2187.94                       0.345833                    0.316472
case300                       300                  1.11828                    0.71989           1389.1                        0.561618                    0.521425
case1354pegase               1354                  4.02602                    3.22855            309.736                      2.59791                     2.40328
case1888rte                  1888                  5.00038                    4.3515             229.806                      3.57255                     3.35274
case2848rte                  2848                  7.56474                    6.75564            148.025                      5.57038                     5.2524
case2869pegase               2869                  8.76664                    7.45625            134.116                      5.98939                     5.55479
case3120sp                   3120                  8.51066                    7.65408            130.649                      6.33616                     6.04839
case6495rte                  6495                 20.3295                    18.7453              53.3467                    15.8663                     15.1643
case6515rte                  6515                 21.6117                    20.005               49.9874                    17.1227                     16.4003
case9241pegase               9241                 34.1935                    30.0892              33.2345                    24.9988                     23.4816
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
case14                         14        0.00441931        226280
case118                       118        0.0289332          34562.3
case_illinois200              200        0.0571606          17494.6
case300                       300        0.135418            7384.55
case1354pegase               1354        0.770959            1297.09
case1888rte                  1888        0.965818            1035.39
case2848rte                  2848        1.49533              668.749
case2869pegase               2869        1.79082              558.402
case3120sp                   3120        1.4179               705.271
case6495rte                  6495        4.26454              234.492
case6515rte                  6515        4.28735              233.244
case9241pegase               9241        7.34957              136.062
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
case14                         14           0.00964653           103664
case118                       118           0.0397558             25153.5
case_illinois200              200           0.0883376             11320.2
case300                       300           0.164824               6067.09
case1354pegase               1354           0.890424               1123.06
case1888rte                  1888           1.11261                 898.786
case2848rte                  2848           1.75219                 570.714
case2869pegase               2869           1.95654                 511.105
case3120sp                   3120           2.03015                 492.575
case6495rte                  6495           4.83098                 206.997
case6515rte                  6515           4.83759                 206.715
case9241pegase               9241           8.05451                 124.154
================  ===============  ===================  ===================

Comments
--------

Allowing lightsim2grid to "recycle" previous computation (column `avg step duration (ms)`, default behaviour) instead of restarting from scratch at every step makes grid2op between **~1.2x** (on `case14`) and **~1.9x** (on `case2869pegase`) faster, depending on the grid size.

Compared to a regular grid2op step (with recycling), the `TimeSerie` module is between **~2.6x** (on `case2869pegase`) and **~69.5x** (on `case14`) faster.

Similarly, the `ContingencyAnalysis` module is between **~2.3x** (on `case3120sp`) and **~31.8x** (on `case14`) faster than a regular grid2op step (with recycling) to evaluate one contingency.
