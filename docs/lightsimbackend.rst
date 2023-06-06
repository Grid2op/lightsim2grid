LightSimBackend
===================================
This is an implementation of a grid2op `Backend <https://grid2op.readthedocs.io/en/latest/backend.html>`_ that uses lightsim2grid simulator coded in c++.

The integration with grid2op is rather easy. You simply need to provide the key-word argument
`backend=LightSimBackend()` when building your environment using the `grid2op.make` function and you
can use it transparently.

Example
--------
See the section :ref:`use_with_g2op` for more information and more examples.

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

Customization of the backend
-------------------------------
.. warning::
    Use grid2op > 1.7.1 for this feature to work properly. Otherwise some bugs (hard to detect) will occur.

You can customize the way the backend behaves in different ways:

- `max_iter`: maximum number of iterations you allow the solver to perform. If a valid solution to the Kirchhoff Current Laws (KCL)
  is not found after this number of iterations, then the backend will "diverge". Default is 10 which is a good value for medium size
  powergrid if you use a Newton Raphson based method (default)
- `tol`: During its internal iterations, the underlying solver will say the Kirchhoff Current Laws (KCL) are matched if the 
  maximum value of the difference is lower than this. Default is `1e-8`.
- `solver_type`: which type of "solver" you want to use. See :ref:`solvers_doc` for more information. By default it uses 
  what it considers the fastest solver available which is likely to be :class:`lightsim2grid.solver.SolverType.KLUSolverSingleSlack`
- `turned_off_pv` : by default (set to `turned_off_pv=True`) all generators partipate in the voltage regulation, which is not completely realistic.
  When you initialize a backend with `turned_off_pv=False` then the generators that do not produce power (*eg* "p=0.") or that are
  turned off are excluded from the voltage regulation.
- `dist_slack_non_renew`: by default in most grid2op environment, the slack bus is "centralize" / "single slack". This parameters
  allows to bypass this restriction and use all non renewable generators (and turned on and with  > 0.) in a distributed
  slack bus setting. It might change the default `solver_type` used.
- \* `detailed_infos_for_cascading_failures`: for exhaustivity, do not modify.
- \* `can_be_copied`: for exhaustivity, do not modify.

The easiest way to customize your backend is when you create the grid2op environment, like this:

.. code-block:: python

    import grid2op
    import lightsim2grid
    from lightsim2grid import LightSimBackend

    env_name = ...
    env = grid2op.make(env_name,
                       backend=LightSimBackend(
                        max_iter=15,
                        tol=1e-9,
                        solver_type=lightsim2grid.solver.SolverType.KLUSolverSingleSlack,
                        # etc.
                        )
                      )

Detailed documentation
--------------------------

.. automodule:: lightsim2grid.lightSimBackend
    :members:
    :autosummary:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
