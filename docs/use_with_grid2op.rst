.. _use_with_g2op:

Use with grid2op
==================

The preferred way to use light2im simulator is with Grid2op. And in this case, you can simply use it
this way.

The integration with grid2op is rather easy. You simply need to provide the key-word argument
`backend=LightSimBackend()` when building your environment using the `grid2op.make` function and you
can use it transparently.

**NB** By default, the fastest resolution method (among KLU and SparseLU linear solver is used). For now it
is not easy to change it. If this is of any interest for you, please let us know with a
`feature request <https://github.com/BDonnot/lightsim2grid/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=>`_

*************************
Regular environments
*************************
For standard grid2op environment, you can use it like:

.. code-block:: python

    import grid2op
    from lightsim2grid.LightSimBackend import LightSimBackend
    from grid2op.Agent import RandomAgent

    # create an environment
    env_name = "rte_case14_realistic"  # for example, other environments might be usable
    env = grid2op.make(env_name,
                       backend=LightSimBackend()  # this is the only change you have to make!
                       )

    # create an agent
    my_agent = RandomAgent(env.action_space)

    # proceed as you would any open ai gym loop
    nb_episode = 10
    for _ in range(nb_episde):
        # you perform in this case 10 different episodes
        obs = env.reset()
        reward = env.reward_range[0]
        done = False
        while not done:
            # here you loop on the time steps: at each step your agent receive an observation
            # takes an action
            # and the environment computes the next observation that will be used at the next step.
            act = agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            # the `LightSimBackend` will be used to carry out the powerflow computation instead
            # of the default grid2op `PandaPowerBackend`

This also works transparently if you use the grid2op `Runner` or the grid2op `MultiEnvironment` paradigm.
Some working examples are given bellow (do not hesitate to consult the documentation of grid2Op
for more information about these features `there <https://grid2op.readthedocs.io/en/latest/>`_

*************************
Runner
*************************
For the grid2op `Runner <https://grid2op.readthedocs.io/en/latest/runner.html>`_ you can use it as
follow:

.. code-block:: python

    import grid2op
    from grid2op.Runner import Runner
    from grid2op.Agent import RandomAgent
    from lightsim2grid.LightSimBackend import LightSimBackend

    env_name = "rte_case14_realistic"  # for example, other environments might be usable
    env = grid2op.make(env_name,
                       backend=LightSimBackend()  # this is the only change you have to make!
                       )

    NB_EPISODE = 10  # assess the performance for 10 episodes, for example
    NB_CORE = 2  # do it on 2 cores, for example
    PATH_SAVE = "agents_log"  # and store the results in the "agents_log" folder
    runner = Runner(**env.get_params_for_runner(), agentClass=RandomAgent)
    runner.run(nb_episode=NB_EPISODE, nb_process=NB_CORE, path_save=PATH_SAVE)

*************************
MultiProcessing
*************************
Finally if you want to use multi processing to get even faster, you can use the grid2op multi processing
class transparently too (see their documentation
`here <https://grid2op.readthedocs.io/en/latest/environment.html#grid2op.Environment.MultiEnvMultiProcess>`_)


.. code-block:: python

    import grid2op
    from grid2op.Environment import MultiEnvMultiProcess
    from lightsim2grid.LightSimBackend import LightSimBackend
    env_name = "rte_case14_realistic"  # for example, other environments might be usable
    # create an environment
    env0 = grid2op.make(env_name,
                       backend=LightSimBackend()
                       )
    # create a second environment, that can be similar, or not
    env1 = grid2op.make(env_name,
                       backend=LightSimBackend()
                       )
    # it is recommended to filter or create the environment with different parameters, otherwise this class
    # is of little interest
    envs = [env0, env1]  # list of all environments created
    nb_envs = [1, 7]  # number of "copies" of each environment that will be made.
    # in this case the first one will be copied only once, and the second one 7 times.
    # the total number of environments used in the multi env will be the sum(nb_envs), here 8.

    multi_env = MultiEnvMultiProcess(envs=envs, nb_envs=nb_envs)
    # and now you can use it like any other grid2op environment (almost)
    observations = multi_env.reset()

Benchmarks
************

In this section we will expose some brief benchmarks about the use of lightsim2grid in the grid2op settings.
The code to run these benchmarks are given with this package int the [benchmark](./benchmarks) folder.

All of them has been run on a computer with a the following characteristics:

- system: Linux 5.8.0-63-generic
- OS: ubuntu 20.04
- processor: Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz
- python version: 3.8.10.final.0 (64 bit)
- numpy version: 1.18.5
- pandas version: 1.1.4
- pandapower version: 2.6.0
- grid2op version: 1.6.3
- lightsim2grid version: 0.5.4


To run the benchmark `cd` in the [benchmark](./benchmarks) folder and type:

.. code-block:: bash

    python3 benchmark_solvers.py --name l2rpn_case14_sandbox --no_test --number 1000
    python3 benchmark_solvers.py --name l2rpn_neurips_2020_track2_small --no_test --number 1000

(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

We compare 5 different backends:

- **PP**: PandaPowerBackend (default grid2op backend) which is the reference in our benchmarks (uses the numba
  acceleration). It is our reference solver.
- **LS+GS** (LightSimBackend+Gauss Seidel): the grid2op backend based on lightsim2grid that uses the "Gauss Seidel"
  solver to compute the powerflows It is implemented in `GaussSeidelSolver.h`.
- **LS+GS S** (LightSimBackend+Gauss Seidel Synchronous): the grid2op backend based on lightsim2grid that uses a
  variant of the "Gauss Seidel" method to compute the powerflows. It is implemented in `GaussSeidelSynchSolver.h`.
- **LS+SLU** (Newton Raphson+SparseLU): the grid2op backend based on lightsim2grid that uses the
  "Newton Raphson" algorithm coupled with the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform) and is implemented in `SparseLUSolver.h`.
- **LS+KLU** (Newton Raphson+KLU): he grid2op backend based on lightsim2grid that uses the
  "Newton Raphson" algorithm coupled with the linear solver
  "KLU" from the `SuiteSparse` c package implemented in `KLUSolver.h`.
- **LS+NICSLU** (Newton Raphson+NICSLU): he grid2op backend based on lightsim2grid that uses the
  "Newton Raphson" algorithm coupled with the linear solver
  "NICSLU" implemented in
  `[NICSLUSolver.h`. [**NB** NICSLU is a free software but not open source, in order to use
  it with lightsim2grid, you need to check section
  `(optional) Include NICSLU linear solver (experimental)` of the Readme file.
  It is required to install lightsim2grid from source for such solver]

Computation time
~~~~~~~~~~~~~~~~~~~

In this first subsection we compare the computation times:

- **grid2op speed** from a grid2op point of view
  (this include the time to compute the powerflow, plus the time to modify the powergrid plus the
  time to read back the data once the powerflow has run plus the time to update the environment and
  the observations etc.). It is reported in "iteration per second" (`it/s`) and represents the number of grid2op "step"
  that can be computed per second.
- **grid2op powerflow time** corresponds to the time the solver take to perform a powerflow
  as seen from grid2op (counting the resolution time and some time to check the validity of the results but
  not the time to update the grid nor the grid2op environment). It is reported in milli seconds (ms).
- **solver powerflow time** corresponds to only the time spend in the solver itself. It does not take into
  account any of the checking, nor the reading back of the data etc. It is reported in milli seconds (ms).


First on an environment based on the IEEE case 14 grid:

================  ======================  ===================================  ============================
case14_sandbox      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
================  ======================  ===================================  ============================
PP                                  60.8                               12.8                          4.91
LS+GS                              800                                  0.458                        0.334
LS+GS S                            790                                  0.47                         0.347
LS+SLU                             992                                  0.214                        0.0886
LS+KLU                            1030                                  0.177                        0.0507
LS+NICSLU                         1020                                  0.178                        0.0513
================  ======================  ===================================  ============================

From a grid2op perspective, lightsim2grid allows to compute up to ~1000 steps each second on the case 14 and
"only" 61 for the default PandaPower Backend, leading to a speed up of **~16** in this case
(lightsim2grid is ~16 times faster than Pandapower). For such a small environment, there is no sensible
difference in using
KLU linear solver compared to using the SparseLU solver of Eigen (1030 vs 992 iterations on the reported
runs, might slightly vary across runs). KLU and NICSLU achieve almost identical performances.

Then on an environment based on the IEEE case 118:

=====================  ======================  ===================================  ============================
neurips_2020_track2      grid2op speed (it/s)    grid2op 'backend.runpf' time (ms)    solver powerflow time (ms)
=====================  ======================  ===================================  ============================
PP                                      38.6                                13.6                           5.67
LS+GS                                    5.39                              184                           184
LS+GS S                                 31.8                                30.5                          30.2
LS+SLU                                 533                                   1.06                          0.791
LS+KLU                                 706                                   0.602                         0.329
LS+NICSLU                              704                                   0.603                         0.33
=====================  ======================  ===================================  ============================

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU) is **~18** time faster than
using the default PandaPower backend. The speed up of lightsim + SparseLU is a bit lower, but it is still **~10**
times faster than using the default backend [the `LS+KLU` solver is ~2-3 times faster than the `LS+SLU` solver
(`0.33` ms per powerflow for `L2+KLU`  compared to `0.79` ms for `LS+SLU`), but it only translates to `LS+KLU`
providing ~30-40% more
iterations per second in the total program (`706` vs `533`) mainly because grid2op itself takes some times to modify the
grid and performs all the check it does.] For this testcase once again there is no noticeable difference between
`NICSLU` and `KLU`.

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc.) we can notice that it takes on average (over 1000 different states) approximately **0.33ms**
to compute a powerflow with the LightSimBackend (if using the KLU linear solver) compared to the **6.6 ms** when using
the PandaPowerBackend (speed up of **~18** times)

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

.. info::

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

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
