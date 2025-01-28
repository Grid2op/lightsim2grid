.. _benchmark-dc-solvers:

Benchmarks (dc solvers)
========================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings when
performing DC powerflow.

TODO DOC in progress

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
    matrices. They can be obtained with lightsim2grid and allow to perform much more powerflow.

When using grid2op, for these small environment, the difference in computation time for an AC or a DC powerflow 
is neglectible. Because the Newton-Raphson algorithm has been much more optimized, it is even faster to run
AC powerflow than DC powerflow for the case 14 (2600 steps per second for AC and 2400 steps per second for DC).
For the bigger case 118 the DC environment is slightly faster (2100 steps per second for DC vs 2000 steps per second for AC).

.. note::
  If you want to be faster in grid2op, switching to DC powerflow instead of AC will probably not be 
  a good solution if you use lightsim2grid.

Lightsim2grid is still much faster than pandapower (*eg* for case 118, 2000 steps / s for lightsim2grid and 180 for
pandapower) and pypowsybl (*eg* for case 118: 650 steps per second for pypowsybl and 2000 for lightsim2grid).

Last, but not least, if you want to perform DC computations and knows in advance the generations and loads
and the topology of the grid, then you probably should use the PTDF and LODF matrices. With them, 
using a matrix multiplication (and numpy) you can run (on one CPU core) multiple millions of 
DC powerflows each second.

Machine used on the benchmarks
-------------------------------

In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
compared with pandapower and pypowsybl when using grid2op.

All of them has been run on a computer with a the following characteristics:

- system: Linux 6.5.0-1024-oem
- OS: ubuntu 22.04
- processor: 13th Gen Intel(R) Core(TM) i7-13700H
- python version: 3.9.21.final.0 (64 bit)
- numpy version: 1.26.4
- pandas version: 2.2.3
- pandapower version: 2.14.10
- pywposybl version: 1.9.0.dev1
- pypowsybl2grid version: 0.1.0
- grid2op version: 1.10.5
- lightsim2grid version: 0.10.0
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 

Command to replicate the benchmark on your machine
----------------------------------------------------

To run the benchmark `cd` in the [benchmark](./benchmarks) folder and install the dependencies
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
PP DC                                        204                               3.58                        0.624
pypowsybl                                   1020                               0.609                       0.558
DC                                          2330                               0.0481                      0.0057
DC (KLU)                                    2380                               0.0437                      0.00174
DC (NICSLU \*)                              2370                               0.0437                      0.00178
DC (CKTSO \*)                               2370                               0.0436                      0.00165
time serie \*\*                       NA                                       0.00122758                  0.000456796
PTDF \*\*                             NA                                       7.46892e-05                 7.3936e-05
contingency analysis \*\*\*           NA                                       0.00455775                  0.0010859
LODF \*\*\*                           NA                                       0.0004457                   0.0003768
===========================  ======================  ========================================  ==========================

And for an environment based on the IEEE case 118:

===========================  ======================  ========================================  ==========================
neurips_2020_track2            grid2op speed (it/s)    grid2op 'backend.runpf' time (ms / pf)    time in 'algo' (ms / pf)
===========================  ======================  ========================================  ==========================
PP DC                                        184                               4.01                        0.8
pypowsybl                                    655                               1.1                         1.02
DC                                          1850                               0.0919                      0.0423
DC (KLU)                                    2050                               0.054                       0.00677
DC (NICSLU \*)                              2050                               0.0538                      0.00658
DC (CKTSO \*)                               2060                               0.0529                      0.00576
time serie \*\*                        NA                                      0.0113052                   0.00314729
PTDF \*\*                              NA                                      0.000755997                 0.000739758
contingency analysis \*\*\*            NA                                      0.00826261                  0.00303673
LODF \*\*\*                            NA                                      0.000479167                 0.000297538
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

TODO


See TL;DR section at the top of the file.
