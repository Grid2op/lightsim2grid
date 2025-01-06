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

================  ===============  ===================  =====================  =====================  ================================    
grid  name         size (nb bus)    time (recycling)    time (no recycling)     time (`TimeSerie`)       time (`ContingencyAnalysis`)    
================  ===============  ===================  =====================  =====================  ================================    
case14                         14         0.0303807          0.0688201               0.00972245                   0.0344761       
case118                       118         0.167014           0.383771                0.0651537                    0.0940448  
case_illinois200              200         0.366178           0.747475                0.152984                     0.251852    
case300                       300         0.592916           1.21181                 0.379875                     0.467905   
case1354pegase               1354         3.13735            5.17859                 1.58152                      2.01299         
case1888rte                  1888         4.78187            7.58089                 2.08743                      2.72081      
case2848rte                  2848         7.49326           12.0294                  3.22694                      4.19178    
case2869pegase               2869         7.06508           12.0486                  3.64617                      4.62894         
case3120sp                   3120         8.54887           13.4784                  3.04654                      4.64494         
case6495rte                  6495        26.4778            37.8204                 10.9002                      12.5037              
case6515rte                  6515        29.8737            42.66                   11.38                        12.9684         
case9241pegase               9241        36.0544            55.4857                 16.6537                      20.0572      
================  ===============  ===================  =====================  =====================  ================================        

All timings reported above are in milliseconds (ms) for one powerflow (in all cases lots of powerflow are carried out, up to a thousands
and the timings here are averaged accross all the powerflows performed)

For detailed explanation about each column as well as the hardware used, please refer to the section below, but in summary:

- benchmark were run on python 3.12 with an old laptop (i7-6820HQ CPU @ 2.70GHz from 2015) see section :ref:`bench_grid_size_hardware`
  for more information
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

- date: 2024-10-18 09:35  CEST
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.12.4.final.0 (64 bit)
- numpy version: 1.26.4
- pandas version: 2.2.3
- pandapower version: 2.14.10
- grid2op version: 1.10.5
- lightsim2grid version: 0.9.2
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

- `time in 'gridmodel' (ms / pf)` gives the time it takes to only perform the AC powerflow:

  - converting the provided data into valid matrix / vector to run an AC powerflow
  - computing the AC Powerflow
  - post processing the internal data (which includes *eg* the flows on the lines in amps, the reactive value
    produced / absorbed by each generator etc.)
    
- `time in 'pf algo' (ms / pf)` gives the time spent in the algorithm that computes the AC powerflow only

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

================  ===============  ========================  ==========================  ================  ===============================  =============================
grid                size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'gridmodel' (ms / pf)    time in 'pf algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ===============================  =============================
case14                         14                  0.758799                   0.0597669        16731.7                           0.0303807                      0.0250171
case118                       118                  0.913219                   0.211025          4738.78                          0.167014                       0.149728
case_illinois200              200                  1.18555                    0.424583          2355.25                          0.366178                       0.340139
case300                       300                  1.44624                    0.661998          1510.58                          0.592916                       0.557392
case1354pegase               1354                  5.26387                    3.37046            296.695                         3.13735                        2.9635
case1888rte                  1888                  6.32057                    5.04453            198.234                         4.78187                        4.58628
case2848rte                  2848                  9.52315                    7.88586            126.809                         7.49326                        7.19927
case2869pegase               2869                 10.428                      7.51632            133.044                         7.06508                        6.70432
case3120sp                   3120                 10.6149                     9.01426            110.935                         8.54887                        8.24586
case6495rte                  6495                 30.5814                    27.5533              36.2933                       26.4778                        25.6759
case6515rte                  6515                 34.0398                    30.9591              32.3007                       29.8737                        29.0781
case9241pegase               9241                 46.1182                    37.7921              26.4606                       36.0544                        34.7085
================  ===============  ========================  ==========================  ================  ===============================  =============================

Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (**no recycling allowed**, non default)

================  ===============  ========================  ==========================  ================  ===============================  =============================
grid                size (nb bus)    avg step duration (ms)    time [DC + AC] (ms / pf)    speed (pf / s)    time in 'gridmodel' (ms / pf)    time in 'pf algo' (ms / pf)
================  ===============  ========================  ==========================  ================  ===============================  =============================
case14                         14                  0.777772                    0.119986         8334.27                          0.0688201                      0.0567457
case118                       118                  1.26015                     0.531649         1880.94                          0.383771                       0.343062
case_illinois200              200                  1.77514                     0.961583         1039.95                          0.747475                       0.688786
case300                       300                  2.39949                     1.52385           656.232                         1.21181                        1.12254
case1354pegase               1354                  8.08618                     6.32786           158.031                         5.17859                        4.75853
case1888rte                  1888                 10.3294                      9.00365           111.066                         7.58089                        7.0991
case2848rte                  2848                 16.0491                     14.2892             69.9832                       12.0294                        11.2664
case2869pegase               2869                 17.6752                     14.6977             68.0376                       12.0486                        11.0712
case3120sp                   3120                 17.6044                     15.9006             62.8906                       13.4784                        12.7485
case6495rte                  6495                 46.697                      43.6531             22.9079                       37.8204                        35.8113
case6515rte                  6515                 51.8558                     48.7368             20.5184                       42.66                          40.588
case9241pegase               9241                 74.1648                     65.6422             15.2341                       55.4857                        51.7239
================  ===============  ========================  ==========================  ================  ===============================  =============================


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
case14                         14        0.00972245       102855
case118                       118        0.0651537         15348.3
case_illinois200              200        0.152984           6536.64
case300                       300        0.379875           2632.45
case1354pegase               1354        1.58152             632.305
case1888rte                  1888        2.08743             479.059
case2848rte                  2848        3.22694             309.891
case2869pegase               2869        3.64617             274.26
case3120sp                   3120        3.04654             328.241
case6495rte                  6495       10.9002               91.7417
case6515rte                  6515       11.38                 87.8737
case9241pegase               9241       16.6537               60.0467
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
case14                         14            0.0344761           29005.6
case118                       118            0.0940448           10633.2
case_illinois200              200            0.251852             3970.58
case300                       300            0.467905             2137.18
case1354pegase               1354            2.01299               496.774
case1888rte                  1888            2.72081               367.537
case2848rte                  2848            4.19178               238.562
case2869pegase               2869            4.62894               216.032
case3120sp                   3120            4.64494               215.288
case6495rte                  6495           12.5037                 79.9763
case6515rte                  6515           12.9684                 77.1104
case9241pegase               9241           20.0572                 49.8575
================  ===============  ===================  ===================

