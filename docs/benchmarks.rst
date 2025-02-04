.. _benchmark-solvers:

Benchmarks (solvers)
======================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

TODO DOC in progress

If you are interested in other type of benchmark, let us know !

.. note::
  Benchmarks performed here does not reflect performance of your usage of grid2op. In particular, there is no "agent",
  the time to "reset" the environment is not taken into account etc.

  The objective is to assess the performance of lightsim2grid (compared to pandapower and pypowsybl) when used as
  a grid2op backend. Hence all timings exclude part of the code not related with this goal.

TL;DR
---------

The default pandapower backend is the slowest one.

Pypowsybl is, for this benchmark, always 3-4 times faster than pandapower (*ie* it allows to perform
3 times more grid2op iteration in the same amount of time).

Lightsim2grid is ~20 times faster than pandapower (it allows to perform 20 times more interation
for the same amount of time).

When diving into lightsim2grid, the faster resolution method is "Newton-Raphson" (with or without distributed slack),
closely followed by the "Fast Decoupled Powerflow" method. 

Concerning the linear solver used for the lightsim2grid backends, it has little to no impact if you use KLU, NICSLU or 
CKTSO. To avoid license issues, we then recommend to use the KLU linear solver (if you compiled lightsim2grid from source,
we strongly recommend you to perform its installation too).

The default linear solver in Eigen (called here Sparse LU) is much slower than KLU and you should avoid to use it if possible
(it's probably better to take some more time to build lightsim2grid with KLU).
But if you don't have the choice, SparseLU still give better performance than pypowsybl or pandapower.

We do not recommend to use Gauss Seidel method (much slower than everything else).

.. note::
  More time have been spent to optimize Newton-Raphson algorithm in lightsim2grid than to optimize Gauss Seidel or
  Fast Decoupled method. The results here should not be extrapolate outside of lightsim2grid usage.

  At time of writing only a prelimary version of the backend based on pypowsybl was avaialble, speed improvement
  might be achieved for the official release.

Using a grid2op environment
----------------------------

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

To run the benchmark `cd` in the [benchmark](./benchmarks) folder and install the dependencies
(we suppose here that you have already installed lightsim2grid):

.. code-block:: bash
  pip install -r req_benchmarks.txt

This will install the required packages to run the benchmark smoothly (most notably `grid2op` and `numba`)
and then you can start the benchmark with the following commands:

.. code-block:: bash

    python3 benchmark_solvers.py --env_name l2rpn_case14_sandbox --no_test --number 1000
    python3 benchmark_solvers.py --env_name l2rpn_neurips_2020_track2_small --no_test --number 1000

.. note::
  The first time you execute them, some data might be downloaded. These data comes from the grid2op
  package and contains the time series (how each load and generation behaves) for different steps.

(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

We compare up to 20 different "solvers" (combination of "linear solver used" (*eg* Eigen, KLU, CKTSO, NICSLU)
and powerflow algorithm (*eg* "Newton Raphson", or "Fast Decoupled")):

- **PP**: PandaPowerBackend (default grid2op backend) which is the reference in our benchmarks (uses the numba
  acceleration). It is our reference solver.
- **PP (no numba)** : is still the `PandaPowerBackend` but it does not use the numba accelaration (see 
  pandapower documentation for more information)
- **PP (with lightsim)**: is again the `PandaPowerBackend` which uses lightsim2grid to compute the powerflow (drr
  pandapower documentation for more information)
- **pypowsybl**: benchmarks the "pypowsybl2grid" backend based on the python implementation of the powsybl framework.
- **GS** (Gauss Seidel): the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows.
- **GS synch** (Gauss Seidel synch version): the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows.
- **NR single (SLU)** (Newton Raphson -single slack- with SparseLU): the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform). This solver supports distributed slack bus.
- **NR (SLU)** (Newton Raphson -distributed slack- with SparseLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **NR (KLU)** (Newton Raphson -distributed slack- with KLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "KLU" from the `SuiteSparse` C package. This solver supports distributed slack bus.
- **NR single (KLU)** (Newton Raphson -single slack- with KLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **NR (NICSLU *)** (Newton Raphson -distributed slack- with NICSLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "NICSLU". [**NB** NICSLU is a free software but not open source, in order to use
  it with lightsim2grid, you need to install lightsim2grid from source for such solver]
- **NR single (NICSLU *)** (Newton Raphson -single slack- with NICSLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **NR (CKTSO *)** (Newton Raphson -distributed slack- with CKTSO): the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "CKTSO". [**NB** CKTSO is a free software but not open source, in order to use
  it with lightsim2grid, you need to install lightsim2grid from source for such solver]
- **NR single (CKTSO *)** (Newton Raphson -single slack- with CKTSO): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **FDPF XB (SLU)** (Fast Decoupled Powerflow, XB variant - with SparseLU linear solver): It is the lightsim2grid 
  implementation of the Fast Decoupled powerflow (in its "XB" variant) that uses the native linear solver in 
  Eigen (called SparseLU in this documentation)
- **FDPF BX (SLU)** (Fast Decoupled Powerflow, BX variant - with SparseLU linear solver): It is the lightsim2grid 
  implementation of the Fast Decoupled powerflow (in its "BX" variant) that uses the native linear solver in 
  Eigen (called SparseLU in this documentation)
- **FDPF XB (KLU)** (Fast Decoupled Powerflow, XB variant - with KLU linear solver) same as `FDPF XB (SLU)` but using KLU instead 
  of SparseLU
- **FDPF BX (KLU)** (Fast Decoupled Powerflow, BX variant - with KLU linear solver) same as `FDPF BX (SLU)` but using KLU instead 
  of SparseLU
- **FDPF XB (NICSLU *)** (Fast Decoupled Powerflow, XB variant - with NICSLU linear solver) same as `FDPF XB (SLU)` but using NICSLU instead 
  of SparseLU
- **FDPF BX (NICSLU *)** (Fast Decoupled Powerflow, BX variant - with NICSLU linear solver) same as `FDPF BX (SLU)` but using NICSLU instead 
  of SparseLU
- **FDPF XB (CKTSO *)** (Fast Decoupled Powerflow, XB variant - with CKTSO linear solver) same as `FDPF XB (SLU)` but using CKTSO instead 
  of SparseLU
- **FDPF BX (CKTSO *)** (Fast Decoupled Powerflow, BX variant - with CKTSO linear solver) same as `FDPF BX (SLU)` but using CKTSO instead 
  of SparseLU

**NB** all backends above (except pandapower) are implemented in lightsim2grid.

**NB** solver with \* are available provided that lightsim2grid is installed from source and following the instructions 
in the documentation. Some license might be needed.

All benchmarks where done with all the customization (for speed, *eg* `-O3` and `-march=native` for linux). 
See the readme for more information.

.. warning::
  At time of writing only a development version of the powsybl backend was available. We will update these figures when 
  the first version will be available.


Computation time
~~~~~~~~~~~~~~~~~~~

In this first subsection we compare the computation times (detailed explanation on page :ref:`benchmark-deep-dive`):

- **grid2op speed** from a grid2op point of view
  (this include the time to compute the powerflow, plus the time to modify 
  the powergrid plus the
  time to read back the data once the powerflow has run plus the time to update 
  the environment and
  the observations etc.). It is reported in "iteration per second" (`it/s`) and 
  represents the number of grid2op "step"
  that can be computed per second.
- **grid2op 'backend.runpf' time** corresponds to the time the grid2op backend take 
  to perform a powerflow
  as seen from grid2op (counting the resolution time of the powerflow solver and some time to check 
  the validity of the results but
  not the time to update the grid nor the grid2op environment), for lightsim2grid 
  it includes the time to read back the data
  from c++ to python. It is reported in milli seconds (ms).
- **solver powerflow time** corresponds only to the time spent in the solver 
  itself. It does not take into
  account any of the checking, nor the transfer of the data python side, nor the 
  computation of the initial state (done through DC approximation) etc. 
  It is reported in milli seconds (ms) as well.

There are two major differences between **grid2op 'backend.runpf' time** and **solver powerflow time**. In **grid2op 'backend.runpf' time**
the time to initialize the solver (usually with the DC approximation) is counted (it is not in **solver powerflow time**). Secondly,
in **grid2op 'backend.runpf' time** the time to read back the data is also included. This explain why **grid2op 'backend.runpf' time** is
stricly greater, for all benchmarks, than **solver powerflow time** (the closer it is, the better the implementation of the LightSimBackend)


First on an environment based on the IEEE case 14 grid:

====================  ======================  ===================================  ==========================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    time in 'algo' (ms / pf)
====================  ======================  ===================================  ==========================
PP                                     130                                 6.37                       2.83
PP (no numba)                           86.3                              10.2                        6.69
PP (with lightsim)                     121                                 6.97                       1.47
pypowsybl                              722                                 0.996                      0.94
GS                                    1590                                 0.311                      0.261
GS synch                              1500                                 0.348                      0.299
NR single (SLU)                       2420                                 0.0847                     0.0338
NR (SLU)                              2400                                 0.0867                     0.0353
NR single (KLU)                       2650                                 0.059                      0.00993
NR (KLU)                              2640                                 0.0591                     0.00959
NR single (NICSLU *)                  2630                                 0.0596                     0.0101
NR (NICSLU *)                         2640                                 0.0592                     0.00973
NR single (CKTSO *)                   2640                                 0.0589                     0.00957
NR (CKTSO *)                          2640                                 0.0586                     0.00915
FDPF XB (SLU)                         2610                                 0.0642                     0.0152
FDPF BX (SLU)                         2560                                 0.0714                     0.0222
FDPF XB (KLU)                         2640                                 0.0613                     0.0125
FDPF BX (KLU)                         2590                                 0.0675                     0.0186
FDPF XB (NICSLU *)                    2630                                 0.0613                     0.0124
FDPF BX (NICSLU *)                    2570                                 0.0679                     0.0186
FDPF XB (CKTSO *)                     2630                                 0.0616                     0.0127
FDPF BX (CKTSO *)                     2580                                 0.0681                     0.0189
====================  ======================  ===================================  ==========================


From a grid2op perspective, lightsim2grid allows to compute up to ~2600 steps each second (column `grid2op speed`, rows `NR XXX`) 
on the case 14 and
"only" ~130 for the default PandaPower Backend (column `grid2op speed`, row `PP`), leading to a speed up of 
**~20** (2600 / 130) in this case
(lightsim2grid Backend is ~20 times faster than pandapower Backend when comparing grid2op speed). 

When compared to powsybl (with the pypowsybl backend), lightsim2grid (with newton raphson) is around 4 times faster (720 vs 2600).

For such a small environment, there is no sensible
difference in using `KLU` linear solver  (rows `NR single (KLU)` or  `NR (KLU)`) compared to using the SparseLU solver of Eigen 
(rows `NR single (SLU)` or  `NR (SLU)`)
(2420 vs 2650 iterations on the reported runs, might slightly vary across runs). 

Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, 
at least we think the observed differences are within error margins.

There are also very little differences between non distributed slack (`NR Single` rows) 
and distributed slack (`NR` rows) for all of the tested linear solvers.

Finally, the "fast decoupled" methods also leads to equivalent performances for almost all linear solvers and are slightly 
less performant than the Newton Raphson one.

For this small environment, for lightsim2grid backend (and if we don't take into account the "agent time"), the computation time 
is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) 
takes 0.38ms (`1. / 2640. * 1000.`) on average and on this 380 ns (or 0.38ms), only 
9.5 ns are spent in the backend. Meaning that 370 ns are spent in the grid2op extra layer or
in the backend implementation in this case (97% of the computation time 
- `=370 / 380`- is external to the powerflow solver)

Then on an environment based on the IEEE case 118:

=====================  ======================  ===================================  ==========================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    time in 'algo' (ms / pf)
=====================  ======================  ===================================  ==========================
PP                                     108                                   7.76                       4.07
PP (no numba)                           73.7                                12.2                        8.47
PP (with lightsim)                     104                                   8.27                       1.93
pypowsybl                              275                                   3.18                       3.09
GS                                       7.23                              138                        138
GS synch                                43.8                                22.5                       22.4
NR single (SLU)                       1120                                   0.496                      0.42
NR (SLU)                              1090                                   0.516                      0.439
NR single (KLU)                       1960                                   0.137                      0.0679
NR (KLU)                              1980                                   0.132                      0.0628
NR single (NICSLU *)                  1970                                   0.133                      0.0647
NR (NICSLU *)                         1990                                   0.129                      0.0608
NR single (CKTSO *)                   1990                                   0.129                      0.0608
NR (CKTSO *)                          2000                                   0.126                      0.0571
FDPF XB (SLU)                         1810                                   0.183                      0.117
FDPF BX (SLU)                         1750                                   0.202                      0.135
FDPF XB (KLU)                         1890                                   0.163                      0.0967
FDPF BX (KLU)                         1830                                   0.178                      0.112
FDPF XB (NICSLU *)                    1890                                   0.162                      0.0965
FDPF BX (NICSLU *)                    1830                                   0.178                      0.112
FDPF XB (CKTSO *)                     1890                                   0.162                      0.0965
FDPF BX (CKTSO *)                     1830                                   0.178                      0.112
=====================  ======================  ===================================  ==========================

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU) is **~18** time faster than
using the default `PandaPower` backend (~2000 it/s vs ~110 for pandapower with numba). 

When compared to powsybl (with the pypowsybl backend), lightsim2grid (with newton raphson) is around **7** times faster (280 vs 2000).

The speed up of lightsim + SparseLU (`1100` it / s) is a bit lower (than using KLU, CKTSO or NICSLU), but it is still **~10**
times faster than using the default backend (1100 / 110). 

For this environment the `NR (KLU)` solver (solver powerflow time) is ~6-7 times faster than the `NR (SLU)` solver 
(`0.439` ms per powerflow for `NR (SLU)` compared to `0.0628` ms for `NR (KLU)`), but it only translates to `NR (KLU)` 
providing ~70% more iterations per second in the total program (`1100` vs `1980`) 
mainly because grid2op itself takes some times to modify the
grid and performs some consistency cheks. 

For this test case once again there is no noticeable difference between `NICSLU`, `CKTSO` and `KLU`.

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc. --  *ie* looking at the column "`solver powerflow time (ms)`") 
we can notice that it takes on average (over 1000 different states) approximately **0.0628**
to compute a powerflow with the LightSimBackend (if using the Newton-Raphson algorithm with the `KLU` linear solver) 
compared to the **4.1 ms** when using
the `PandaPowerBackend` (with numba, but without lightsim2grid) (speed up of **~65** times)

For this small environment, once again, for lightsim2grid backend (and if we don't take into account the "agent time"), 
the computation time 
is vastly dominated by factor external to the powerflow solver. Indeed, doing a 'env.step' (column `grid2op speed (it/s)`) 
takes 0.50ms (`1. / 2000. * 1000.`) on average and on this 500 ns (or 0.50ms), only 
63 ns are spent in the powerflow algorithm. Meaning that 440 ns are spent 
in the grid2op extra layer or in the grid2op backend in this case (~90% of the computation time 
- `=440 / 500`- is external to the powerflow solver)

.. note:: The "solver powerflow time" reported for pandapower is obtained by summing, over the 1000 powerflow performed
    the `pandapower_backend._grid["_ppc"]["et"]` (the "estimated time" of the pandapower newton raphson computation)

    For the lightsim backend, the "solver powerflow time" corresponds to the average of the results of
    `gridmodel.get_computation_time()` function that, for each powerflow, returns the time spent in the solver
    uniquely (time inside the `basesolver.compute_pf()` function. In particular it do not count the time
    to initialize the vector `V` with the DC approximation nor the time taken to convert the lightsim2grid external modeling 
    to something that can be processed by the powerflow algorithm). This is a behaviour similar to the one of pandapower.

    For pypowsybl backend, the column `solver powerflow time (ms)` also counts the time to convert the data from the 
    iidm modeling to the data model used by the "Open Load Flow" powerflow solver.

More information
~~~~~~~~~~~~~~~~~

For all the benchmarks, the stopping criteria (of the solver, for each steps) is based on the maximum discrepency (at each bus) of the Kirchhoff Current Law's. 
The solver stops if this `tolerance` is under `1e-8` "per unit".

For all "steps", the model is intialized with the "DC" (direct current) approximation before the AC powerflow is run. In these benchmark, only the timings
of the AC powerflow is reported in the column "`solver powerflow time (ms)`". This time is however included in the other columns.

These benchmarks mimic a typical behaviour in grid2op where in most cases the "agent" does nothing: between two consecutive steps, 
there is no modification of the topology. Only the injection (active and reactive power consumed and active generation as well as target voltage 
setpoint of generators are modified between two steps). The topology does not change, the tap of the transformers stay the same etc.

Finally, as opposed to pandapower Backend, lightsim2grid Backend is able to "recycle" partially some of its comptuation. Concretely this means, in this 
case, that the Ybus matrix is not recomputed at each steps (but computed only at the first one) for example. This can lead to some time savings in these
cases.

.. warning::
  For more information about what is actually done and the wordings used in this section, 
  you can consult the page :ref:`benchmark-deep-dive`


Differences
~~~~~~~~~~~~~~~~~~~
Using the same command, we report the maximum value of the differences when compared with the reference 
implementation which in this case is pandapower (this explain why pandapower always have a difference of 0.) for different 
output values :

- `aor` : the current flow (in Amps) at the origin side of each powerline
- `gen_p` : the generators active production values
- `gen_q`: the generators reactive production values

Note that only the maximum values (of the absolute differences) across all the steps (1000 for the IEEE case 14 and
1000 for the IEEE case 118)
and across all the lines (or generators) is displayed.

We report only the difference compared with the baseline which is pandapower (PP).

Here are the results for the IEEE case 14 (max over 1000 powerflows):

============================  ==============  ==============  ================
case14_sandbox (1000 iter)      Δ aor (amps)    Δ gen_p (MW)    Δ gen_q (MVAr)
============================  ==============  ==============  ================
PP (ref)                            0               0                 0
GS                                  0.000122        7.63e-06          7.63e-06
GS synch                            0.000122        7.63e-06          7.63e-06
NR single (SLU)                     0.000122        7.63e-06          7.63e-06
NR (SLU)                            0.000122        7.63e-06          7.63e-06
NR single (KLU)                     0.000122        7.63e-06          7.63e-06
NR (KLU)                            0.000122        7.63e-06          7.63e-06
NR single (NICSLU *)                0.000122        7.63e-06          7.63e-06
NR (NICSLU *)                       0.000122        7.63e-06          7.63e-06
NR single (CKTSO *)                 0.000122        7.63e-06          7.63e-06
NR (CKTSO *)                        0.000122        7.63e-06          7.63e-06
FDPF XB (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (SLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF BX (KLU)                       0.000122        7.63e-06          7.63e-06
FDPF XB (NICSLU *)                  0.000122        7.63e-06          7.63e-06
FDPF BX (NICSLU *)                  0.000122        7.63e-06          7.63e-06
FDPF XB (CKTSO *)                   0.000122        7.63e-06          7.63e-06
FDPF BX (CKTSO *)                   0.000122        7.63e-06          7.63e-06
============================  ==============  ==============  ================

.. note::

    Flows are here measured in amps (and not kA). The maximum difference of flows is approximately 0.1mA
    or 1e-4 A. This difference is totally neglectible on power transportation side where the current is usually
    around 1kA (1e3 A).

Here are the results for the IEEE case 118 (max over 1000 powerflows):

=================================  ==============  ==============  ================
neurips_2020_track2 (1000 iter)      Δ aor (amps)    Δ gen_p (MW)    Δ gen_q (MVAr)
=================================  ==============  ==============  ================
PP (ref)                                  0              0                 0
GS                                        6.1e-05        3.81e-06          1.53e-05
GS synch                                  6.1e-05        3.81e-06          1.53e-05
NR single (SLU)                           6.1e-05        0                 9.54e-07
NR (SLU)                                  6.1e-05        0                 9.54e-07
NR single (KLU)                           6.1e-05        0                 9.54e-07
NR (KLU)                                  6.1e-05        0                 9.54e-07
NR single (NICSLU *)                      6.1e-05        0                 9.54e-07
NR (NICSLU *)                             6.1e-05        0                 9.54e-07
NR single (CKTSO *)                       6.1e-05        0                 9.54e-07
NR (CKTSO *)                              6.1e-05        0                 9.54e-07
FDPF XB (SLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (SLU)                             6.1e-05        0                 9.54e-07
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        0                 9.54e-07
FDPF XB (NICSLU *)                        6.1e-05        1.91e-06          1.53e-05
FDPF BX (NICSLU *)                        6.1e-05        0                 9.54e-07
FDPF XB (CKTSO *)                         6.1e-05        1.91e-06          1.53e-05
FDPF BX (CKTSO *)                         6.1e-05        0                 9.54e-07
=================================  ==============  ==============  ================

As we can see on all the tables above, the difference when using lightsim and pandapower is rather
small, even when using a different algorithm to solve the powerflow (LS + GS corresponds to
using Gauss Seidel as opposed to using Newton Raphson solver)

When using Newton Raphson solvers, the difference in absolute values when using lightsim2grid compared
with using PandaPowerBackend is neglectible: less than 1e-06 in all cases (and 0.00 when comparing the
flows on the powerline for both environments).

.. note::
  The differences reported here are in comparison with pandapower. This is why there is 0. to 
  all the columns corresponding to the `PP (ref)` row.

Other benchmark
----------------

We have at our disposal different computers with different software / hardware.

From time to time, we benchmark grid2op and lightsim2grid. 
The results can be found in:

.. toctree::
  :maxdepth: 1
  :caption: For a laptop with a i7 of 2015 wth a frequency of 2.70 GHz

  benchmark_solver/ubuntu_2004_dell/ls0.8.1_glop1.10.1
  benchmark_solver/ubuntu_2004_dell/ls0.8.0_glop1.10.0
  benchmark_solver/ubuntu_2004_dell/ls0.8.0_glop1.9.8

.. toctree::
  :maxdepth: 1
  :caption: For a laptop with a i7 of 2014 wth a frequency of 3.0 GHz

  benchmark_solver/ubuntu_2004_server/ls0.8.1_glop1.10.1
  benchmark_solver/ubuntu_2004_server/ls0.8.1_glop1.10.1_py311
  benchmark_solver/ubuntu_2004_server/ls0.8.1_glop1.10.1_py312

.. toctree::
  :maxdepth: 1
  :caption: For a laptop with a ryzen 7 of 2020 wth a frequency of 4.2 GHz

  benchmark_solver/windows_10_portable/ls0.8.1_glop1.9.7_py38
  benchmark_solver/windows_10_portable/ls0.8.1_glop1.9.6_py38
  benchmark_solver/windows_10_portable/ls0.8.1_glop1.9.5_py38
  
.. toctree::
  :maxdepth: 1
  :caption: For a desktop with a i7 of 2014 with a frequency of 4.00GHz

  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.10.0_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.10.0_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.10.1_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.10.1_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.0_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.0_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.1_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.1_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.2_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.2_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.3_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.3_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.4_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.4_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.5_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.5_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.6_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.6_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.7_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.7_py3.9
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.8_py3.8
  benchmark_solver/ubuntu_2004_desktop/ls0.8.2_glop1.9.8_py3.9


.. note::
  Any contribution here is more than welcomed. You can write a github discussion here 
  https://github.com/Grid2Op/lightsim2grid/discussions/new?category=show-and-tell 
  and describe rapidly your setup and we'll make sure to include your benchmark in future release.

  Thanks !


Benchmarks of other lightsim2grid functions
--------------------------------------------

With lightsim2grid 0.5.5 some new feature has been introduced, which are the "security analysis" and the "comptuation 
of time series". 

The respective benchmarks are put in their respective section :ref:`sa_benchmarks` and :ref:`ts_benchmarks`. These 
classes allow to achieve a *15x* and even *100x* speed ups over grid2op (using lightsim2grid), for example 
allowing to perform 186 powerflow on the IEEE 118 in less than 3 ms.
