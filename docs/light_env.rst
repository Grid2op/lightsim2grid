Light environment
=======================================

Goal
--------------------------

:class:`lightsim2grid.lightEnv.LightEnv` is a (very) limited grid2op environment written
entirely in C++: the loop ``reset`` / ``step`` / protections / cooldowns / observation runs
without going back to python, and without grid2op's bookkeeping around it. It is meant for
the workloads where the environment itself is the bottleneck: training an agent that only
acts on the topology, or a search (tree search, beam search, ...) that evaluates many
candidate actions one step ahead.

It is not a replacement of grid2op. It models one thing -- a grid replaying a time series of
injections, on which an agent changes the topology -- and it follows grid2op's semantics for
that one thing, so an agent can be trained or searched on the light environment and run on
grid2op afterwards. See :ref:`light-env-differences` for what it leaves out.

On ``l2rpn_case14_sandbox`` a step is **25 to 50 times faster** than a grid2op step with
:class:`lightsim2grid.LightSimBackend`, and evaluating a candidate action is **about 20 times
faster** than ``obs.simulate`` (see :ref:`light-env-benchmark`).

.. note::

    The light environment is new and its interface may still change. It is also usable from
    C++ with no python at all: it is header-only, in ``src/core/light_env/`` (installed with
    the core library, see :doc:`cpp_library`).

Quick start
--------------------------

The easiest way to build one is from a grid2op environment using
:class:`lightsim2grid.LightSimBackend`: the grid of its backend is the initial state of every
episode, and the thermal limits, the protections and the cooldowns are taken from its
parameters. This is the setup used by the benchmark (``benchmarks/light_env.py``):

.. code-block:: python

    import numpy as np
    import grid2op
    from lightsim2grid import LightSimBackend
    from lightsim2grid.lightEnv import LightEnv, Protections

    env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend())
    env.reset()

    # the grid the episode starts from (the topology is restored at every reset)
    light_env = LightEnv(env.backend._grid)

    # the injections: one row per step, one column per element, NaN meaning "unchanged"
    data = env.chronics_handler.real_data.data
    nb_row = data.load_p.shape[0]
    light_env.assign_time_series(
        data.load_p.astype(float),
        data.load_q.astype(float),
        data.prod_p.astype(float),
        (data.prod_v / env.backend.prod_pu_to_kv).astype(float),  # in pu, not kV
        np.full((nb_row, env.n_storage), np.nan),
        np.full((nb_row, env.n_shunt), np.nan),
        np.full((nb_row, env.n_shunt), np.nan),
        np.full((nb_row, 0), np.nan),
        np.full((nb_row, 0), np.nan),
    )

    # the protections: thermal limits in kA (grid2op uses A), on both sides of each line
    th_lim_ka = env.get_thermal_limit() * 1e-3
    protections = Protections()
    protections.set_thermal_limit_or(th_lim_ka)
    protections.set_thermal_limit_ex(th_lim_ka * env.backend.lines_or_pu_to_kv / env.backend.lines_ex_pu_to_kv)
    protections.set_max_line_time_step_overflow(
        np.full(env.n_line, env.parameters.NB_TIMESTEP_OVERFLOW_ALLOWED, dtype=np.int32))
    light_env.protections = protections

    # the cooldowns (grid2op parameters of the same name)
    light_env.nb_timestep_cooldown_sub = env.parameters.NB_TIMESTEP_COOLDOWN_SUB
    light_env.nb_timestep_cooldown_line = env.parameters.NB_TIMESTEP_COOLDOWN_LINE
    light_env.nb_timestep_reconnection = env.parameters.NB_TIMESTEP_RECONNECTION

    # the actions the agent can take: step(i) plays actions[i]
    actions = [env.action_space()] + list(env.action_space.get_all_unitary_topologies_set(env.action_space))
    light_env.init_actions(actions)

    obs, info = light_env.reset()
    done = False
    while not done:
        act_id = 0  # your agent here
        obs, reward, done, truncated, info = light_env.step(act_id)

Setting it up
--------------------------

**The grid.** ``LightEnv(grid)`` takes an :class:`lightsim2grid.network.LSGrid` and copies it:
this copy is the state every episode starts from (``reset`` restores its topology). A grid
loaded by :class:`lightsim2grid.LightSimBackend` also carries the position of every element in
grid2op's topology vector, which ``obs.topo_vect`` needs; a grid loaded directly (eg with
``init_from_pandapower``) works, but has no ``topo_vect``.

**The time series** (``assign_time_series``). Nine matrices, one row per step, one column per
element: ``load_p``, ``load_q``, ``gen_p`` (MW / MVAr), ``gen_v`` (**pu**), ``storage_p``,
``shunt_p``, ``shunt_q``, ``sgen_p``, ``sgen_q``. A NaN leaves the value of the grid unchanged.
They must all have the same number of rows, which becomes ``max_step``; the sizes are checked
at the next ``reset``. For now only ``load_p``, ``load_q``, ``gen_p`` and ``gen_v`` are
applied, the other five are only checked.

**The protections** (:class:`lightsim2grid.lightEnv.Protections`). The thermal limits, in kA,
on both sides of every line (grid2op numbering: powerlines then transformers, the "or" side
of a transformer being its hv side), and for each line the number of steps it can stay in
overflow. A line whose overflow counter exceeds that number is disconnected, the powerflow is
run again, and so on until no more line trips (a cascade, within one step). This is grid2op's
"soft" overflow protection (``NB_TIMESTEP_OVERFLOW_ALLOWED``).

**The cooldowns.** ``nb_timestep_cooldown_sub``, ``nb_timestep_cooldown_line`` and
``nb_timestep_reconnection`` behave as the grid2op parameters of the same name: after an
action touches a substation (resp. changes the status of a line) that substation (resp. line)
cannot be acted on for that many steps, and a line disconnected by the protections cannot be
reconnected before ``nb_timestep_reconnection`` steps.

**The actions** (``init_actions``). A list of grid2op actions or of
:class:`lightsim2grid.lightEnv.TopoAction`, mixing both is fine. Only ``set_bus`` and
``set_line_status`` are modelled; a grid2op action using anything else that affects the grid
(``change_bus``, ``redispatch``, ``curtail``, storage...) is refused. Every action is checked
against the initial grid (the element exists, the busbar exists, no contradiction) and an
invalid one raises a ``ValueError`` naming it; nothing is registered in that case. Without
any action registered, ``step(0)`` (do nothing) is the only valid call.

A :class:`lightsim2grid.lightEnv.TopoAction` can also be built by hand, with grid2op semantics:

.. code-block:: python

    from lightsim2grid.lightEnv import TopoAction, ElementType

    act = TopoAction()
    act.add_element(ElementType.load, 0, 2)      # set_bus: load 0 on busbar 2 of its substation
    act.add_element(ElementType.line_or, 2, 2)   # origin side of line 2 on busbar 2
    act.set_line_status(3, -1)                   # set_line_status: disconnect line 3

Playing an episode
--------------------------

``reset()`` restores the initial topology, the overflow counters and the cooldowns, applies
the first row of the time series, runs a powerflow and returns ``(obs, info)``. It must be
called before the first step, and again after ``assign_time_series`` or a change of the
protections.

``step(act_id)`` returns ``(obs, reward, done, truncated, info)``:

- the action ``act_id`` is played, unless it touches a substation or a line still in
  cooldown: it is then replaced by "do nothing" and ``info["is_illegal"]`` is ``"true"``;
- the next row of the time series is applied, and the powerflow is run with the protections
  (the cascade described above);
- ``reward`` is the fraction of the episode survived so far, ``current_step / max_step``;
- ``done`` is ``True`` at the end of the time series (``info["success"] == "true"``, reward 1)
  or when the powerflow diverges (``info["failure"] == "true"``, reward 0, ``truncated`` is
  ``True``). ``info["survival_time"]`` is then the fraction of the episode survived.

The values in ``info`` are strings.

The observation
--------------------------

``obs`` is a :class:`lightsim2grid.lightEnv.LightEnvObservation`, a **view** on the state of
the environment: every attribute is a read-only numpy array on the environment's memory,
nothing is copied. ``reset`` and every ``step`` return the same object, whose values follow
the environment. Use ``np.array(obs.p_or)`` to keep a snapshot.

======================================== ==================================================================
attribute                                meaning
======================================== ==================================================================
``rho``                                  max over both sides of current / thermal limit, per line
``p_or``, ``q_or``, ``a_or``             flows at the origin side (MW, MVAr, **kA**), per line
``p_ex``, ``q_ex``, ``a_ex``             flows at the extremity side (MW, MVAr, **kA**), per line
``load_p``, ``gen_p``                    active power of the loads / generators (MW)
``topo_vect``                            grid2op topology vector (local busbar, -1 if disconnected)
``time_before_cooldown_line``            steps before the status of each line can change again
``time_before_cooldown_sub``             steps before each substation can be acted on again
``current_step``                         step of the environment (0 at reset)
======================================== ==================================================================

Lines follow grid2op numbering. Currents are in kA, where grid2op uses A.

A view is valid until the next ``reset``, or a step ending the episode by a divergence (the
grid then drops its results); once ``done`` is ``True`` the values are not meaningful.

Copying: lookahead and search
------------------------------

A :class:`lightsim2grid.lightEnv.LightEnv` can be copied (``env.copy()``, ``copy.copy``,
``copy.deepcopy``, or ``LightEnv(env)``, all the same). The copy is an independent environment
at the same point of the same episode, with its own observation. What does not change during
an episode -- the initial grid, the time series and the registered actions -- is shared
read-only between the copies, so a copy is cheap (about 0.015 ms on ``l2rpn_case14_sandbox``).
The first step of a copy is slower than a step of the original, about 0.05 ms instead of
0.02 ms: a copied grid starts with a cold solver (no cached ``Ybus`` nor factorisation) and
rebuilds it. This is the light environment's counterpart of grid2op's ``obs.simulate``:

.. code-block:: python

    def best_action(light_env, candidates):
        """the candidate that leaves the lowest max(rho) one step ahead"""
        best, best_rho = 0, np.inf
        for act_id in candidates:
            sim = light_env.copy()
            obs, reward, done, truncated, info = sim.step(act_id)
            if done or info["is_illegal"] == "true":
                continue
            if obs.rho.max() < best_rho:
                best, best_rho = act_id, obs.rho.max()
        return best

Unlike ``obs.simulate``, the copy steps on the *actual* next row of the time series, not on a
forecast: it looks at the future the environment will really see.

.. _light-env-differences:

What it does not model
--------------------------

Compared with a grid2op environment:

- the only actions are ``set_bus`` and ``set_line_status``: no ``change_bus`` /
  ``change_line_status``, no redispatching, curtailment, storage, or setpoint change;
- only ``load_p``, ``load_q``, ``gen_p`` and ``gen_v`` are replayed: storage units, shunts and
  static generators keep the values of the initial grid;
- no maintenance, no opponent, no forecasts, no alarms / alerts;
- no limit on the number of substations or lines an action touches (``MAX_SUB_CHANGED``,
  ``MAX_LINE_STATUS_CHANGED``), only the cooldowns;
- only the "soft" overflow protection: no instantaneous disconnection above
  ``HARD_OVERFLOW_THRESHOLD``, and a line is in overflow when ``rho >= 1`` (grid2op:
  above ``SOFT_OVERFLOW_THRESHOLD``, 1 by default);
- a single reward, the fraction of the episode survived;
- every powerflow starts from a flat DC initialisation, rather than from the previous state.

Within these limits it follows grid2op step by step: ``lightsim2grid/tests/test_LightEnv.py``
plays the same actions on both and compares the flows, the topology and the cooldowns at every
step.

.. _light-env-benchmark:

Benchmark
--------------------------

``benchmarks/light_env.py`` plays the same chronics of ``l2rpn_case14_sandbox`` on a grid2op
environment with :class:`lightsim2grid.LightSimBackend` and on a light environment built as
in the quick start, with the same thermal limits and protections. grid2op's
``HARD_OVERFLOW_THRESHOLD`` is raised out of reach, so both sides model the same thing (see
above). Three workloads:

- **do nothing**: ``env.step(do_nothing)`` vs ``light_env.step(0)``;
- **topology**: every 10 steps a random unitary ``set_bus`` action (the same sequence of
  actions on both sides), do nothing otherwise. Such an agent ends its episode quickly, so 50
  episodes are played;
- **lookahead**: at every step, 10 candidate actions are evaluated one step ahead
  (``obs.simulate(act)`` vs ``light_env.copy().step(act_id)``), then the episode goes on with
  do nothing. Only the evaluation is timed.

To run it (``cd`` into the ``benchmarks`` folder, grid2op installed):

.. code-block:: bash

    python light_env.py            # the test chronics shipped with grid2op (3 x 575 steps)
    python light_env.py --no_test  # the full dataset, downloaded by grid2op

The results below were obtained on the test chronics, on a shared cloud machine (Intel Xeon
@ 2.10GHz, python 3.11, grid2op 1.12.5, lightsim2grid 1.0.1rc0 with KLU, compiled without
``-O3`` / ``-march=native``). Timings vary by up to a factor 2 from one run to the next on such
a machine; the ratio between the two columns is what matters.

"ms / step (in env)" is the time the environment accounts for itself (grid2op
``env._time_step``, the light env ``step_time``) and "powerflow" the time spent in the
powerflows, protections included (grid2op ``env._time_powerflow``, the light env
``protections.powerflow_time``, which also counts the powerflow of ``reset``). Both sides play
the same number of steps and end the same episodes: they agree on the physics.

Do nothing:

=========================== ======= =========== =========== =========== ==================== =======================
\                           steps   game over   steps / s   ms / step   ms / step (in env)   ms / step (powerflow)
=========================== ======= =========== =========== =========== ==================== =======================
grid2op + LightSimBackend   1725    0           1153        0.868       0.790                0.187
LightEnv                    1725    0           54503       0.018       0.017                0.015
=========================== ======= =========== =========== =========== ==================== =======================

Topology (a random unitary ``set_bus`` action every 10 steps, 50 episodes):

=========================== ======= =========== =========== =========== ==================== =======================
\                           steps   game over   steps / s   ms / step   ms / step (in env)   ms / step (powerflow)
=========================== ======= =========== =========== =========== ==================== =======================
grid2op + LightSimBackend   1246    50          1294        0.773       0.706                0.176
LightEnv                    1246    50          31620       0.032       0.029                0.030
=========================== ======= =========== =========== =========== ==================== =======================

Lookahead (10 candidate actions per step, one step ahead):

=========================== ======= =========== =========== ===========
\                           steps   game over   steps / s   ms / step
=========================== ======= =========== =========== ===========
grid2op ``obs.simulate``    17220   0           634         1.577
LightEnv copy + step        17220   0           12160       0.082
=========================== ======= =========== =========== ===========

In short, over three runs: a step is 25 to 50 times faster, and a one-step lookahead about 20
times faster. Most of a grid2op step is spent outside the powerflow (0.19 ms of 0.87 ms here):
applying the action, the rules, building the observation. The light environment spends
almost all of its time in the powerflow, and the python binding adds about 1 µs per step (the
gap between "ms / step" and "ms / step (in env)"). A step with a topological change is
slower than a do nothing one (0.03 ms vs 0.02 ms), as the admittance matrix has to be rebuilt,
and the lookahead is dominated by the first step of each copy, which starts from a cold
solver (see `Copying: lookahead and search`_).

Detailed usage
--------------------------

.. automodule:: lightsim2grid.lightEnv
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
