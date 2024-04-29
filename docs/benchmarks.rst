
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

- date: 2024-03-25 17:53  CET
- system: Linux 5.15.0-56-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
- python version: 3.10.13.final.0 (64 bit)
- numpy version: 1.23.5
- pandas version: 2.2.1
- pandapower version: 2.13.1
- grid2op version: 1.10.1
- lightsim2grid version: 0.8.1
- lightsim2grid extra information: 

	- klu_solver_available: True 
	- nicslu_solver_available: True 
	- cktso_solver_available: True 
	- compiled_march_native: True 
	- compiled_o3_optim: True 


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

====================  ======================  ===================================  ============================
case14_sandbox          grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
====================  ======================  ===================================  ============================
PP                                      46.3                               18.4                          6.57
GS                                     757                                  0.474                        0.378
GS synch                               769                                  0.445                        0.348
NR single (SLU)                        960                                  0.184                        0.0831
NR (SLU)                               952                                  0.189                        0.0819
NR single (KLU)                       1030                                  0.12                         0.0221
NR (KLU)                              1030                                  0.118                        0.0202
NR single (NICSLU *)                  1020                                  0.121                        0.022
NR (NICSLU *)                         1020                                  0.119                        0.02
NR single (CKTSO *)                   1020                                  0.119                        0.0211
NR (CKTSO *)                           989                                  0.121                        0.0192
FDPF XB (SLU)                         1010                                  0.13                         0.032
FDPF BX (SLU)                         1010                                  0.143                        0.0451
FDPF XB (KLU)                         1020                                  0.124                        0.0263
FDPF BX (KLU)                         1010                                  0.134                        0.0377
FDPF XB (NICSLU *)                    1010                                  0.126                        0.0267
FDPF BX (NICSLU *)                    1020                                  0.134                        0.0383
FDPF XB (CKTSO *)                     1010                                  0.125                        0.0268
FDPF BX (CKTSO *)                     1000                                  0.136                        0.0381
====================  ======================  ===================================  ============================

From a grid2op perspective, lightsim2grid allows to compute up to ~1200 steps each second on the case 14 and
"only" 70 for the default PandaPower Backend, leading to a speed up of **~17** in this case
(lightsim2grid is ~17 times faster than `Pandapower`). For such a small environment, there is no sensible
difference in using `KLU` linear solver compared to using the SparseLU solver of Eigen (1120 vs 1200 iterations on the reported
runs, might slightly vary across runs). `KLU` and `NICSLU` achieve almost identical performances.

Then on an environment based on the IEEE case 118:

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      41.5                                20.7                           8.6
GS                                       3.74                              266                           266
GS synch                                35.8                                26.9                          26.8
NR single (SLU)                        536                                   0.897                         0.767
NR (SLU)                               505                                   0.959                         0.818
NR single (KLU)                        811                                   0.268                         0.144
NR (KLU)                               820                                   0.256                         0.131
NR single (NICSLU *)                   813                                   0.259                         0.134
NR (NICSLU *)                          827                                   0.243                         0.118
NR single (CKTSO *)                    814                                   0.257                         0.131
NR (CKTSO *)                           829                                   0.24                          0.116
FDPF XB (SLU)                          762                                   0.352                         0.232
FDPF BX (SLU)                          749                                   0.373                         0.252
FDPF XB (KLU)                          786                                   0.307                         0.188
FDPF BX (KLU)                          776                                   0.327                         0.206
FDPF XB (NICSLU *)                     786                                   0.308                         0.188
FDPF BX (NICSLU *)                     771                                   0.324                         0.204
FDPF XB (CKTSO *)                      784                                   0.309                         0.19
FDPF BX (CKTSO *)                      773                                   0.329                         0.209
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
FDPF BX (SLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (KLU)                             6.1e-05        1.91e-06          1.53e-05
FDPF BX (KLU)                             6.1e-05        1.91e-06          7.63e-06
FDPF XB (NICSLU *)                        6.1e-05        1.91e-06          1.53e-05
FDPF BX (NICSLU *)                        6.1e-05        1.91e-06          7.63e-06
FDPF XB (CKTSO *)                         6.1e-05        1.91e-06          1.53e-05
FDPF BX (CKTSO *)                         6.1e-05        1.91e-06          7.63e-06
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
  https://github.com/BDonnot/lightsim2grid/discussions/new?category=show-and-tell 
  and describe rapidly your setup and we'll make sure to include your benchmark in future release.

  Thanks !


Benchmarks of other lightsim2grid functions
--------------------------------------------

With lightsim2grid 0.5.5 some new feature has been introduced, which are the "security analysis" and the "comptuation 
of time series". 

The respective benchmarks are put in their respective section :ref:`sa_benchmarks` and :ref:`ts_benchmarks`. These 
classes allow to achieve a *15x* and even *100x* speed ups over grid2op (using lightsim2grid), for example 
allowing to perform 186 powerflow on the IEEE 118 in less than 3 ms.
