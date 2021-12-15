
Benchmarks
============

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

TODO DOC in progress

If you are interested in other type of benchmark, let us know !

Using a grid2op environment
----------------------------
In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
compared with pandapower when using grid2op.

All of them has been run on a computer with a the following characteristics:

- system: Linux 5.11.0-38-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.18.5
- pandas version: 1.1.4
- pandapower version: 2.6.0
- lightsim2grid version: 0.5.5
- grid2op version: 1.6.4


To run the benchmark `cd` in the [benchmark](./benchmarks) folder and type:

.. code-block:: bash

    python3 benchmark_solvers.py --name l2rpn_case14_sandbox --no_test --number 1000
    python3 benchmark_solvers.py --name l2rpn_neurips_2020_track2_small --no_test --number 1000

(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

We compare up to 9 different solvers:

- **PP**: PandaPowerBackend (default grid2op backend) which is the reference in our benchmarks (uses the numba
  acceleration). It is our reference solver.
- **LS+GS** (LightSimBackend+Gauss Seidel): the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows.
- **LS+GS S** (LightSimBackend+Gauss Seidel Synchronous): the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows).
- **LS+SLU** (Newton Raphson+SparseLU): the grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform). This solver supports distributed slack bus.
- **LS+SLU (single)** (Newton Raphson+SparseLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **LS+KLU** (Newton Raphson+KLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "KLU" from the `SuiteSparse` C package. This solver supports distributed slack bus.
- **LS+KLU (single)** (Newton Raphson+KLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.
- **LS+NICSLU** (Newton Raphson+NICSLU): he grid2op backend based on lightsim2grid that uses the 
  "Newton Raphson" algorithm coupled with the linear solver 
  "NICSLU". [**NB** NICSLU is a free software but not open source, in order to use
  it with lightsim2grid, you need to install lightsim2grid from source for such solver]
- **LS+NICSLU (single)** (Newton Raphson+NICSLU): same as above but this solver does not support distributed slack bus and
  can thus be slightly faster.

All benchmarks where done with all the customization (for speed, *eg* `-O3` for linux). See the readme for more information.

Computation time
~~~~~~~~~~~~~~~~~~~

In this first subsection we compare the computation times:

- **grid2op speed** from a grid2op point of view
  (this include the time to compute the powerflow, plus the time to modify the powergrid plus the
  time to read back the data once the powerflow has run plus the time to update the environment and
  the observations etc.). It is reported in "iteration per second" (`it/s`) and represents the number of grid2op "step"
  that can be computed per second.
- **grid2op 'backend.runpf' time** corresponds to the time the solver take to perform a powerflow
  as seen from grid2op (counting the resolution time and some time to check the validity of the results but
  not the time to update the grid nor the grid2op environment), for lightsim2grid it includes the time to read back the data
  from c++ to python. It is reported in milli seconds (ms).
- **solver powerflow time** corresponds only to the time spent in the solver itself. It does not take into
  account any of the checking, nor the transfer of the data python side etc. It is reported in milli seconds (ms) as well.

There are two major differences between **grid2op 'backend.runpf' time** and **solver powerflow time**. In **grid2op 'backend.runpf' time**
the time to initialize the solver (usually with the DC approximation) is counted (it is not in **solver powerflow time**). Secondly,
in **grid2op 'backend.runpf' time** the time to read back the data is also included. This explain why **grid2op 'backend.runpf' time** is
stricly greater, for all benchmarks, than **solver powerflow time** (the closer it is, the better the implementation of the LightSimBackend)


First on an environment based on the IEEE case 14 grid:

================  ======================  ===================================  ============================
case14_sandbox      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
================  ======================  ===================================  ============================
PP                                  68.8                               11.3                          4.36
LS+GS                              866                                  0.41                         0.308
LS+GS S                            859                                  0.416                        0.314
LS+SLU                            1090                                  0.168                        0.0632
LS+KLU                            1160                                  0.12                         0.0188
LS+NICSLU                         1150                                  0.121                        0.0186
================  ======================  ===================================  ============================

From a grid2op perspective, lightsim2grid allows to compute up to ~1200 steps each second on the case 14 and
"only" 69 for the default PandaPower Backend, leading to a speed up of **~17** in this case
(lightsim2grid is ~16 times faster than Pandapower). For such a small environment, there is no sensible
difference in using
KLU linear solver compared to using the SparseLU solver of Eigen (1160 vs 1090 iterations on the reported
runs, might slightly vary across runs). KLU and NICSLU achieve almost identical performances.

Then on an environment based on the IEEE case 118:

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      39.4                                13.5                           5.65
LS+GS                                    5.12                              194                           194
LS+GS S                                 35.2                                27.5                          27.3
LS+SLU                                 621                                   0.775                         0.591
LS+KLU                                 883                                   0.301                         0.12
LS+NICSLU                              881                                   0.302                         0.121
=====================  ======================  ===================================  ============================

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU) is **~22** time faster than
using the default PandaPower backend. The speed up of lightsim + SparseLU is a bit lower, but it is still **~16**
times faster than using the default backend [the `LS+KLU` solver is ~4-5 times faster than the `LS+SLU` solver 
(`0.12` ms per powerflow for `L2+KLU`  compared to `0.59` ms for `LS+SLU`), but it only translates to `LS+KLU` 
providing ~30-40% more
iterations per second in the total program (`880` vs `650`) mainly because grid2op itself takes some times to modify the
grid and performs all the check it does.] For this testcase once again there is no noticeable difference between
`NICSLU` and `KLU`.

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc.) we can notice that it takes on average (over 1000 different states) approximately **0.12ms**
to compute a powerflow with the LightSimBackend (if using the KLU linear solver) compared to the **5.6 ms** when using
the PandaPowerBackend (speed up of **~46** times)

**NB** pandapower performances heavily depends on the pandas version used, we used here a version of pandas which
we found gave the best performances on our machine.

.. note:: The "solver powerflow time" reported for pandapower is obtained by summing, over the 1000 powerflow performed
    the `pandapower_backend._grid["_ppc"]["et"]` (the "estimated time" of the pandapower newton raphson computation
    with the numba accelaration enabled)

    For the lightsim backend, the "solver powerflow time" corresponds to the sum of the results of
    `gridmodel.get_computation_time()` function that, for each powerflow, returns the time spent in the solver
    uniquely (time inside the `basesolver.compute_pf()` function. In particular it do not count the time
    to initialize the vector V with the DC approximation)

Differences
~~~~~~~~~~~~~~~~~~~
Using the same command, we report the maximum value of the differences between different quantities:

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
LS+GS                               0.000122        7.63e-06          7.63e-06
LS+GS S                             0.000122        7.63e-06          7.63e-06
LS+SLU                              0.000122        7.63e-06          7.63e-06
LS+KLU                              0.000122        7.63e-06          7.63e-06
LS+NICSLU                           0.000122        7.63e-06          7.63e-06
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
LS+GS                                     6.1e-05        3.81e-06          1.53e-05
LS+GS S                                   6.1e-05        3.81e-06          1.53e-05
LS+SLU                                    6.1e-05        0                 9.54e-07
LS+KLU                                    6.1e-05        0                 9.54e-07
LS+NICSLU                                 6.1e-05        0                 9.54e-07
=================================  ==============  ==============  ================

As we can see on all the tables above, the difference when using lightsim and pandapower is rather
small, even when using a different algorithm to solve the powerflow (LS + GS corresponds to
using Gauss Seidel as opposed to using Newton Raphson solver)

When using Newton Raphson solvers, the difference in absolute values when using lightsim2grid compared
with using PandaPowerBackend is neglectible: less than 1e-06 in all cases (and 0.00 when comparing the
flows on the powerline for both environments).

Other benchmarks
----------------------------

With lightsim2grid 0.5.5 some new feature has been introduced, which are the "security analysis" and the "comptuation 
of time series". 

The respective benchmarks are put in their respective section :ref:`sa_benchmarks` and :ref:`ts_benchmarks`. These 
classes allow to achieve a *15x* and even *100x* speed ups over grid2op (using lightsim2grid), for example 
allowing to perform 186 powerflow on the IEEE 118 in less than 3 ms.