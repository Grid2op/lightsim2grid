Continuation Power Flow
=======================================

Goal
--------------------------

A continuation power flow (CPF) traces the solution curve of a grid as its injections
move from where they are now towards a target state, and reports where that curve turns
back on itself -- the "nose", i.e. the maximum loading the grid can serve at all. It
answers "how much more load can this grid take, and which bus gives way first", which a
sequence of ordinary powerflows cannot: near the collapse point the Newton-Raphson simply
stops converging, and a divergence does not tell you whether the grid is at its limit or
your starting point was bad.

The parameter along the curve is called :math:`\lambda`: :math:`\lambda = 0` is the grid's
own state, :math:`\lambda = 1` the target state. This is MATPOWER's ``runcpf``
formulation, and the option names and defaults below follow its ``cpf.*`` options.

.. note::
    This is a *natural* parameterisation: the corrector solves at a fixed
    :math:`\lambda`, so the curve stops **at** the nose and does not trace the lower
    branch. Rounding the nose needs an arc-length parameterisation (MATPOWER's
    ``cpf.parameterization`` 2 / 3), which is not implemented yet.

Basic usage
--------------------------

.. code-block:: python

    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    from lightsim2grid import run_cpf

    grid = init_from_pandapower(pn.case118())

    res = run_cpf(grid, loading_factor=3.0)
    print(res.msg)
    print(f"collapse at {res.lam_max:.3f} of the requested increase")
    print(f"ie a load {1 + res.lam_max * (3.0 - 1.0):.3f} times the base one")

``res`` is a :class:`lightsim2grid.continuationPowerflow.CPFResult`: ``lam``, ``V``,
``Vm``, ``Va_deg``, ``tangent_lam``, ``lam_max``, ``success`` and ``msg``. The helper
``plot_pv_curve(res)`` draws the PV curve of the worst buses if matplotlib is installed.

Steering which elements move
--------------------------------

By default every load grows together and generation follows, so the direction is
"the whole grid, uniformly". Two vectors change that, one coefficient per element in
``[0, 1]``: **0 holds that element at its base value, 1 moves it fully**.

.. code-block:: python

    import numpy as np

    # only the loads of one area grow; everything else stays put
    steering = np.zeros(len(grid.get_loads()))
    steering[my_area_load_ids] = 1.0
    res = run_cpf(grid, loading_factor=3.0, load_steering=steering)

    # ... and this one grows at only 40% of the rate of the others
    steering[another_load_id] = 0.4

``gen_steering`` is the same, per generator. Its default is 1 for every non-slack
generator and 0 for the slack ones, so generation follows the load and the slack picks up
only the incremental losses -- MATPOWER requires exactly that of a target case ("same as
base case w.r.t. Qg and slack Pg"). Passing ``gen_steering=0`` instead holds generation
fixed and makes the slack supply the whole increase; that is a legitimate study, but it
traces a **different curve with a different nose**, so do not compare the two margins.

``scale_q=True`` (the default) scales each load's reactive power by the same coefficient
as its active power, i.e. at constant power factor. For a direction that no combination
of these expresses, ``direction=`` takes explicit per-element deltas in MW / MVAr:

.. code-block:: python

    res = run_cpf(grid, direction={"load_p": delta_load_p, "gen_p": delta_gen_p})

Performance
--------------------------

The continuation is a batch algorithm
(:class:`lightsim2grid.lightsim2grid_cpp.ContinuationSweepCPP`), not a Python loop around
``ac_pf``: the topology is fixed for the whole curve, so the Jacobian's sparsity pattern
is too, and the **whole curve costs one symbolic factorization** whatever the number of
points. Every corrector after the first only refactorizes, and each predictor is a single
triangular solve reusing the factorization the corrector left standing -- which is what
makes tracing a curve of several hundred points affordable.

You can check this on any run::

    cpf = ContinuationPowerFlow(grid)
    res = cpf.run(loading_factor=3.0)
    assert cpf.cpp.get_linear_solver_stats().nb_analyze == 1

Only the Newton-Raphson algorithms are supported (the predictor needs a Jacobian);
``AlgorithmType.GaussSeidel``, the fast-decoupled family and the DC algorithms are
refused with an explicit message.

Options
--------------------------

Named and defaulted after MATPOWER's ``cpf.*``:

===========================  =========  ==============================================
option                       default    meaning
===========================  =========  ==============================================
``step``                     0.05       nominal step, an arc length along the tangent
``step_min``                 1e-4       failing at this step ends the curve
``step_max``                 0.2        largest step the adaptation may grow to
``adapt_step``               ``False``  adapt the step to the predictor's error
``adapt_step_damping``       0.7        damping of that adaptation
``adapt_step_tol``           1e-3       predictor error the adaptation aims at
``nose_tol``                 1e-5       ``tangent_lam`` below this means "at the nose"
``stop_at_lam``              ``None``   stop at this lambda instead of at the nose
``max_steps``                1000       hard cap on the number of traced points
``exact_tangent``            ``False``  refactorize J at each point before its tangent
===========================  =========  ==============================================

``tangent_lam`` is the tangent's :math:`\lambda` component. With this parameterisation it
is strictly positive and *tends to zero* at the nose -- it is a threshold, never a sign
change.

Detailed usage
--------------------------

.. automodule:: lightsim2grid.continuationPowerflow
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
