.. _use_with_g2op:

Use with grid2op
==================================================

The preferred way to use light2im simulator is with Grid2op. And in this case, you can simply use it
this way.

The integration with grid2op is rather easy. You simply need to provide the key-word argument
`backend=LightSimBackend()` when building your environment using the `grid2op.make` function and you
can use it transparently.

**NB** By default, the fastest resolution method (among KLU and SparseLU linear solver is used). For now it
is not easy to change it. If this is of any interest for you, please let us know with a
`feature request <https://github.com/Grid2Op/lightsim2grid/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=>`_


Regular environments
----------------------------------

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
    # and use this as any grid2op environment, for example

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

Runner
-----------------

For the grid2op `Runner <https://grid2op.readthedocs.io/en/latest/runner.html>`_ using lightsim2grid is
also perfectly transparent once the environment is loaded.


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

MultiProcessing
-----------------

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

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
