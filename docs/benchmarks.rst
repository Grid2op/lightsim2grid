
Benchmarks (solvers)
======================

In this paragraph we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

TODO DOC in progress

If you are interested in other type of benchmark, let us know !

Using a grid2op environment
----------------------------
In this section we perform some benchmark of a `do nothing` agent to test the raw performance of lightsim2grid
compared with pandapower when using grid2op.

All of them has been run on a computer with a the following characteristics:

- system: Linux 5.11.0-40-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.18.5
- pandas version: 1.1.4
- pandapower version: 2.7.0
- lightsim2grid version: 0.6.0
- grid2op version: 1.6.4


To run the benchmark `cd` in the [benchmark](./benchmarks) folder and type:

.. code-block:: bash

    python3 benchmark_solvers.py --env_name l2rpn_case14_sandbox --no_test --number 1000
    python3 benchmark_solvers.py --env_name l2rpn_neurips_2020_track2_small --no_test --number 1000

(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

We compare up to 19 different "solvers" (combination of "linear solver used" (*eg* Eigen, KLU, CKTSO, NICSLU)
and powerflow algorithm (*eg* "Newton Raphson", or "Fast Decoupled")):

- **PP**: PandaPowerBackend (default grid2op backend) which is the reference in our benchmarks (uses the numba
  acceleration). It is our reference solver.
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

**NB** all backend above are implemented in lightsim2grid.

**NB** solver with \* are available provided that lightsim2grid is installed from source and following the instructions 
in the documentation.

All benchmarks where done with all the customization (for speed, *eg* `-O3` and `-march=native` for linux). 
See the readme for more information.

Computation time
~~~~~~~~~~~~~~~~~~~

In this first subsection we compare the computation times:

- **grid2op speed** from a grid2op point of view
  (this include the time to compute the powerflow, plus the time to modify 
  the powergrid plus the
  time to read back the data once the powerflow has run plus the time to update 
  the environment and
  the observations etc.). It is reported in "iteration per second" (`it/s`) and 
  represents the number of grid2op "step"
  that can be computed per second.
- **grid2op 'backend.runpf' time** corresponds to the time the solver take 
  to perform a powerflow
  as seen from grid2op (counting the resolution time and some time to check 
  the validity of the results but
  not the time to update the grid nor the grid2op environment), for lightsim2grid 
  it includes the time to read back the data
  from c++ to python. It is reported in milli seconds (ms).
- **solver powerflow time** corresponds only to the time spent in the solver 
  itself. It does not take into
  account any of the checking, nor the transfer of the data python side etc. 
  It is reported in milli seconds (ms) as well.

There are two major differences between **grid2op 'backend.runpf' time** and **solver powerflow time**. In **grid2op 'backend.runpf' time**
the time to initialize the solver (usually with the DC approximation) is counted (it is not in **solver powerflow time**). Secondly,
in **grid2op 'backend.runpf' time** the time to read back the data is also included. This explain why **grid2op 'backend.runpf' time** is
stricly greater, for all benchmarks, than **solver powerflow time** (the closer it is, the better the implementation of the LightSimBackend)


First on an environment based on the IEEE case 14 grid:

==================  ======================  ===================================  ============================
case14_sandbox        grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
==================  ======================  ===================================  ============================
PP                                    70.5                               11                            4.27
LS+GS                                881                                  0.447                        0.327
LS+GS S                              877                                  0.446                        0.327
LS+SLU (single)                     1110                                  0.191                        0.0655
LS+SLU                              1120                                  0.195                        0.0683
LS+KLU (single)                     1200                                  0.138                        0.0176
LS+KLU                              1180                                  0.141                        0.0188
LS+NICSLU (single)                  1200                                  0.139                        0.0179
LS+NICSLU                           1200                                  0.139                        0.0184
==================  ======================  ===================================  ============================

From a grid2op perspective, lightsim2grid allows to compute up to ~1200 steps each second on the case 14 and
"only" 70 for the default PandaPower Backend, leading to a speed up of **~17** in this case
(lightsim2grid is ~17 times faster than `Pandapower`). For such a small environment, there is no sensible
difference in using `KLU` linear solver compared to using the SparseLU solver of Eigen (1120 vs 1200 iterations on the reported
runs, might slightly vary across runs). `KLU` and `NICSLU` achieve almost identical performances.

Then on an environment based on the IEEE case 118:

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                       39.6                               13.3                           5.58
LS+GS                                     5.3                              188                           188
LS+GS S                                  36.5                               26.6                          26.4
LS+SLU (single)                         642                                  0.775                         0.607
LS+SLU                                  588                                  0.932                         0.769
LS+KLU (single)                         945                                  0.277                         0.116
LS+KLU                                  918                                  0.306                         0.144
LS+NICSLU (single)                      947                                  0.274                         0.11
LS+NICSLU                               929                                  0.298                         0.134
=====================  ======================  ===================================  ============================

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU) is **~24** time faster than
using the default `PandaPower` backend (~950 it/s vs ~40). 

The speed up of lightsim + SparseLU (`0.11`) is a bit lower, but it is still **~16**
times faster than using the default backend [the `LS+KLU` solver is ~5-6 times faster than the `LS+SLU` solver 
(`0.11` ms per powerflow for `L2+KLU`  compared to `0.6` ms for `LS+SLU`), but it only translates to `LS+KLU` 
providing ~40-50% more
iterations per second in the total program (`950` vs `640`) mainly because grid2op itself takes some times to modify the
grid and performs all the check it does.] For this testcase once again there is no noticeable difference between
`NICSLU` and `KLU`.

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc. -- column "solver powerflow time (ms)") we can notice that it takes on average (over 1000 different states) approximately **0.12ms**
to compute a powerflow with the LightSimBackend (if using the `KLU` linear solver) compared to the **5.6 ms** when using
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
LS+SLU (single)                     0.000122        7.63e-06          7.63e-06
LS+SLU                              0.000122        7.63e-06          7.63e-06
LS+KLU (single)                     0.000122        7.63e-06          7.63e-06
LS+KLU                              0.000122        7.63e-06          7.63e-06
LS+NICSLU (single)                  0.000122        7.63e-06          7.63e-06
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
LS+SLU (single)                           6.1e-05        0                 9.54e-07
LS+SLU                                    6.1e-05        0                 9.54e-07
LS+KLU (single)                           6.1e-05        0                 9.54e-07
LS+KLU                                    6.1e-05        0                 9.54e-07
LS+NICSLU (single)                        6.1e-05        0                 9.54e-07
LS+NICSLU                                 6.1e-05        0                 9.54e-07
=================================  ==============  ==============  ================

As we can see on all the tables above, the difference when using lightsim and pandapower is rather
small, even when using a different algorithm to solve the powerflow (LS + GS corresponds to
using Gauss Seidel as opposed to using Newton Raphson solver)

When using Newton Raphson solvers, the difference in absolute values when using lightsim2grid compared
with using PandaPowerBackend is neglectible: less than 1e-06 in all cases (and 0.00 when comparing the
flows on the powerline for both environments).

Other benchmark
----------------

We have at our disposal different computers with different software / hardware.

From time to time, we benchmark grid2op and lightsim2grid. 
The results can be found in:

.. toctree::
  :maxdepth: 1
  :caption: For a laptop with a i7 of 2015 wth a frequency of 2.70 GHz

  benchmark_solver/ubuntu_2004_dell/ls0.8.0_glop1.10.0
  benchmark_solver/ubuntu_2004_dell/ls0.8.0_glop1.9.8

Benchmarks of other lightsim2grid functions
--------------------------------------------

With lightsim2grid 0.5.5 some new feature has been introduced, which are the "security analysis" and the "comptuation 
of time series". 

The respective benchmarks are put in their respective section :ref:`sa_benchmarks` and :ref:`ts_benchmarks`. These 
classes allow to achieve a *15x* and even *100x* speed ups over grid2op (using lightsim2grid), for example 
allowing to perform 186 powerflow on the IEEE 118 in less than 3 ms.