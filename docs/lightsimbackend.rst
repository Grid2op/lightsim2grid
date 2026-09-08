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
    from lightsim2grid import LightSimBackend
    from grid2op.Agent import RandomAgent

    # create an environment
    env_name = "l2rpn_case14_sandbox"  # for example, other environments might be usable
    env = grid2op.make(env_name,
                       backend=LightSimBackend()  # this is the only change you have to make!
                       )

    # create an agent
    my_agent = RandomAgent(env.action_space)

    # proceed as you would any open ai gym loop
    nb_episode = 10
    for _ in range(nb_episode):
        # you perform in this case 10 different episodes
        obs = env.reset()
        reward = env.reward_range[0]
        done = False
        while not done:
            # here you loop on the time steps: at each step your agent receive an observation
            # takes an action
            # and the environment computes the next observation that will be used at the next step.
            act = my_agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            # the `LightSimBackend` will be used to carry out the powerflow computation instead
            # of the default grid2op `PandaPowerBackend`

Customization of the solver
-------------------------------
.. warning::
    Use grid2op > 1.7.1 for this feature to work properly. Otherwise some bugs (hard to detect) will occur.

You can customize the way the backend behaves in different ways:

- `max_iter`: maximum number of iterations you allow the solver to perform. If a valid solution to the Kirchhoff Current Laws (KCL)
  is not found after this number of iterations, then the backend will "diverge". Default is 10 which is a good value for medium size
  powergrid if you use a Newton Raphson based method (default)
- `tol`: During its internal iterations, the underlying solver will say the Kirchhoff Current Laws (KCL) are matched if the 
  maximum value of the difference is lower than this. Default is `1e-8`.
- `algo_type`: which type of powerflow algorithm (combined with which linear solver) you want to use. See
  :ref:`solvers_doc` for more information. By default it uses what it considers the fastest one available, which
  is likely to be :class:`lightsim2grid.algorithm.NRSing_KLU`.

  .. deprecated:: 1.0.0
      This kwarg used to be called `solver_type`. That name is kept for backward
      compatibility (it still works and is mapped to `algo_type`), but it is deprecated: "solver" now
      refers specifically to the *linear* solver (KLU, SparseLU, NICSLU, CKTSO), not the powerflow
      algorithm nor the combination of both that `algo_type` selects. Passing both `solver_type` and
      `algo_type` with different values raises. See :ref:`algorithm_names` for the full naming
      rationale.
- `turned_off_pv` : by default (set to `turned_off_pv=True`) all generators partipate in the voltage regulation, which is not completely realistic.
  When you initialize a backend with `turned_off_pv=False` then the generators that do not produce power (*eg* "p=0.") or that are
  turned off are excluded from the voltage regulation.
- `dist_slack_non_renew`: by default in most grid2op environment, the slack bus is "centralize" / "single slack". This parameters
  allows to bypass this restriction and use all non renewable generators (and turned on and with  > 0.) in a distributed
  slack bus setting. It might change the default `algo_type` used.
- \* `use_static_gen`: bool=False, DO NOT USE AT THE MOMENT. When it will be available, you will be able to load
  both "static" generators (pq generators) and "regular" (pv generators) as generators in lightsim2grid. It does
  not work at the moment and has no effect.
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
                        algo_type=lightsim2grid.algorithm.AlgorithmType.NRSing_KLU,
                        # etc.
                        )
                      )

Customization of the input format
----------------------------------

For a few versions now, we try to extend the capability of lightsim2grid and make it work
with other data "reader". We started this process by allowing to initialize a 
lightsim2grid `LSGrid` from a pypowsybl network.

For example, if you environment contains a grid in the iidm format (native format of pypowsybl networks), 
you can load it with:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend
    from grid2op.Agent import RandomAgent

    # create an environment
    env_with_iidm_as_the_grid_description = ...  # eg a path to a directory containing an iidm file
    env = grid2op.make(env_with_iidm_as_the_grid_description,
                       backend=LightSimBackend(loader_method="pypowsybl")
                       )

.. _lightsimbackend_matpower:

Using a MATPOWER case as the grid of an environment
++++++++++++++++++++++++++++++++++++++++++++++++++++

.. versionadded:: 1.0.1

The third reader is MATPOWER: if your environment ships a `grid.m` (the `.m` script a
MATPOWER case is distributed as) or a `grid.mat` (the binary MATPOWER saves) instead of a
`grid.json` or a `grid.xiidm`, load it with:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend

    # a path to a directory containing a "grid.m" (or a "grid.mat") file
    env_with_matpower_as_the_grid_description = ...
    env = grid2op.make(env_with_matpower_as_the_grid_description,
                       backend=LightSimBackend(loader_method="matpower")
                       )

This goes through :func:`lightsim2grid.network.init_from_matpower` and never builds a
pandapower or a pypowsybl network on the way. Reading a `.m` file needs the optional
`matpowercaseframes` package, a `.mat` file the optional `scipy` one (both come with
``pip install lightsim2grid[matpower]``); an already parsed case needs neither, and can
be handed over directly:

.. code-block:: python

    from pypower.api import case118

    env = grid2op.make(env_path,
                       backend=LightSimBackend(loader_method="matpower",
                                               loader_kwargs={"grid": case118()})
                       )

What MATPOWER does not say, and what lightsim2grid does about it:

- **substations**: MATPOWER has no notion of a busbar section within a bus, so there is
  exactly one grid2op substation per MATPOWER bus. How many busbar sections each of them
  gets is up to you (`n_busbar` of `grid2op.make`, or the `n_busbar_per_sub`
  `loader_kwargs`); the extra ones start empty and deactivated, waiting for a topology
  action.
- **names**: MATPOWER numbers its buses and names nothing, so grid2op's own default names
  are used -- `sub_0`, `load_1_0`, `gen_5_3`, `0_4_1` (a powerline is named after the two
  substations it joins and its own id), ... They are made by grid2op itself
  (`Backend._fill_names_obj`), so they are exactly what any other nameless grid gets. A
  chronics folder addressing the elements by name has to use those.
- **thermal limits**: MATPOWER's `RATE_A` is a branch MVA rating, not the ampere limit
  grid2op works with, and is 0 ("unlimited") in a fair share of the published cases. The
  limits are therefore left open by the loader and are read from the environment's
  ``config.py``, which is where a grid2op environment declares them anyway.
- **nominal voltages**: a case that never leaves per unit leaves its `BASE_KV` column at 0,
  and several of the published ones do. Since lightsim2grid reports voltages in kV, such a bus
  is given a nominal voltage of 1 kV -- so every voltage it reports is numerically its
  per-unit value -- and a warning says so. Set the `BASE_KV` column of the case if you
  want actual kV.
- **generators**: several `mpc.gen` rows on the same bus stay independent generators (they
  are *not* aggregated), and each one is a grid2op generator of its own.

A full example of such an environment lives in `lightsim2grid/tests/case_14_matpower`.

You can also customize the way lightsim2grid works with some extra options:


- `loader_method`: Literal["pandapower", "pypowsybl", "matpower"]: from which grid "file description" 
  the grid will be loaded. If you use `pandapower` then pandapower needs to be installed.
  If you specified `pypowsybl` then pypowsybl needs to be installed on your machine.
  If you specified `matpower` then, unless you hand over an already parsed case, reading a
  ".m" file needs `matpowercaseframes` and reading a ".mat" file needs `scipy`.
- `loader_kwargs` : ``dict``: some customization to use when loading the grid. It is not
  not used when loading the grid from `pandapower`. Please refer to the documentation of
  :attr:`LightSimBackend._loader_kwargs` for more information. 

Other Customization
--------------------

- `stop_if_load_disco`, `stop_if_gen_disco`, `stop_if_storage_disco` : ``Optional[bool] = None``:
  whether to raise a `BackendError` if a load / generator / (producing or absorbing) storage unit
  ends up disconnected. The default, ``None``, defers to grid2op's own ``allow_detachment`` setting
  (grid2op >= 1.11.0 and lightsim2grid >= 0.10.0): ``False`` (do not raise) if detachment is allowed,
  ``True`` (raise) otherwise -- matching the legacy, pre-``allow_detachment`` behaviour. Passing an
  explicit ``True``/``False`` that contradicts what ``allow_detachment`` would otherwise select is
  overridden (with a warning) to stay consistent with it.
- `automatically_disconnect` : ``bool = False``: if ``True``, automatically disconnects any load /
  generator that ends up outside the grid's main connected component instead of raising a "grid not
  connected" error. This should only be used together with grid2op's ``allow_detachment``.
- `gen_slack_id` : ``Optional[int] = None``: id (or name, or a collection of either) of the
  generator(s) that should participate to the slack. Only supported when
  `loader_method="pypowsybl"`, and mutually exclusive with `dist_slack_non_renew` (pick one or
  the other).

Detailed documentation
--------------------------

.. automodule:: lightsim2grid.lightSimBackend
    :members:
    :autosummary:
    :private-members:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
