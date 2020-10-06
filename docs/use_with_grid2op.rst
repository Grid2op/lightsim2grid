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

All of them has been run on a computer with a `Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz` processor. The command
to run the benchmark is (once `cd` in the [benchmark](./benchmarks) folder folder):

.. code-block:: bash

    python3 do_nothing.py --name l2rpn_case14_sandbox  # for environment based on IEEE 14 grid
    python3 do_nothing.py --name l2rpn_neurips_2020_track2  --no_test  # for environment based on IEEE 118 grid

We compare 3 different backends:

- **PP**: PandaPowerBackend (default grid2op backend)
- **LS+SLU** (LightSimBackend+SparseLU): the grid2op backend based on lightsim2grid that uses the linear solver "SparseLU" from the
  Eigen c++ library (available on all platform)
- **LS+KLU** (LightSimBackend+KLU): the grid2op backend based on lightsim2grid that uses the linear solver
  "KLU" from the SuiteSparse c package, available only (for now) on Linux and Mac Os.

Computation time
~~~~~~~~~~~~~~~~~~~

In this first subsection we compare the computation time from a grid2op point of view (number of step
per second) and the time the solver take to perform a powerflow (counting the resolution time and
not the time to update the grid or run all the checks perform by grid2op).

+------------+--------------------+--------------------+-------------------------------+-------------------------------+
|            | IEEE 14 (it / s)   |  IEEE 118 (it / s) |  IEEE 14 (powerflow time, ms) | IEEE 118 (powerflow time, ms) |
+============+====================+====================+===============================+===============================+
| PP         |   68.5             |  42.2              |   12.1                        |  14.3                         |
+------------+--------------------+--------------------+-------------------------------+-------------------------------+
| LS+SLU     |   880              |  458               |   0.20                        |  1.08                         |
+------------+--------------------+--------------------+-------------------------------+-------------------------------+
| LS+KLU     |   880              |  596               |   0.17                        |  0.62                         |
+------------+--------------------+--------------------+-------------------------------+-------------------------------+


(results may vary depending on the hard drive, the ram etc. and are presented here for illustration only)

(we remind that these simulations correspond to simulation on one core of the CPU. Of course it is possible to
make use of all the available cores, which would increase the number of steps that can be performed)

From a grid2op perspective, lightsim2grid allows to compute 880 steps each second on the case 14 and "only" 68.5
for the default PandaPower Backend, leading to a speed up of **~13** in this case (lightsim2grid is 13 times faster
than Pandapower). For such a small environment, there is no difference in using KLU linear solver (not available on
Windows based machine) compared to using the SparseLU solver of Eigen.

For an environment based on the IEEE 118, the speed up in using lightsim + KLU (LS+KLU)
[for now only available on linux and MacOS] is **~14** time faster than
using the default PandaPower backend. The speed up of lightsim + SparseLU is a bit lower, but it is still **~11**
times faster than using the default backend [using sparseLU linear solver is approximately 30% slower than using KLU.]

If we look now only at the time to compute one powerflow (and don't take into account the time to load the data, to
initialize the solver, to modify the grid, read back the results, to perform the other update in the
grid2op environment etc.) we can notice that it takes on average (over 575 different states) approximately **0.62 ms**
to compute a powerflow with the LightSimBackend (if using the KLU linear solver) compared to the **14.3 ms** when using
the PandaPowerBackend.

Differences
~~~~~~~~~~~~~~~~~~~
Using the same command, we report the maximum value of the differences between different quantities:

- `aor` : the current flow (in Amps) at the origin side of each powerline
- `gen_p` : the generators active production values
- `gen_q`: the generators reactive production values

Note that only the maximum values across all the steps (1000 for the IEEE case 14 and 525 for the IEEE case 118)
and across all the lines (or generators) is displayed. This is then an upper bound.

We report only the difference compared with the baseline which is pandapower (PP).

Here are the results for the IEEE case 14 (max over 1000 powerflows):
+------------+--------------------+--------------------+-------------------------------+
|            | IEEE 14 aor        |  IEEE 14 gen_p     |  IEEE 14 gen_q                |
+============+====================+====================+===============================+
| PP         |   0.00             |  0.00              |  0.00                         |
+------------+--------------------+--------------------+-------------------------------+
| LS+SLU     |   0.00             |  0.00              |   0.00                        |
+------------+--------------------+--------------------+-------------------------------+
| LS+KLU     |   4.48e-12         |  0.00              |   9.54e-7                     |
+------------+--------------------+--------------------+-------------------------------+

Here are the results for the IEEE case 118 (max over 575 powerflows):
+------------+--------------------+--------------------+-------------------------------+
|            | IEEE 118 aor       |  IEEE 118 gen_p    |  IEEE 118 gen_q               |
+============+====================+====================+===============================+
| PP         |   0.00             |  0.00              |  0.00                         |
+------------+--------------------+--------------------+-------------------------------+
| LS+SLU     |   0.00             |  0.00              |   0.00                        |
+------------+--------------------+--------------------+-------------------------------+
| LS+KLU     |   4.57e-12         |  0.00              |   9.54e-7                     |
+------------+--------------------+--------------------+-------------------------------+


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
