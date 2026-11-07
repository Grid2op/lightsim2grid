Differentiable powerflow (pytorch)
===================================

.. _differentiable:

A batch of powerflows, differentiable with respect to its injections, as a
pytorch operation: a loss computed on the results of thousands of AC powerflows
can be back-propagated to the loads and generators that produced them. A single
grid's powerflow is differentiable the same way (:class:`GridPowerFlow`), and both
work with the outer loops (``NROuter_*``, ``NROuterDist_*``).

.. code-block:: python

    import torch
    from lightsim2grid.network import init_from_pandapower
    from lightsim2grid.differentiable import BatchCPUPowerFlow

    grid = init_from_pandapower(some_pandapower_net)
    pf = BatchCPUPowerFlow.init_from_grid(grid, nb_thread=8)

    load_p = torch.tensor(..., requires_grad=True)    # (n_scenarios, n_load), MW
    V = pf(load_p=load_p)                             # complex (n_scenarios, n_bus)

    loss = pf.compute_flows(V).i_or_a.relu().sum()    # eg an overload penalty
    loss.backward()
    load_p.grad                                       # (n_scenarios, n_load)

``pytorch`` is **not** a dependency of lightsim2grid: this subpackage is the only
thing that needs it, and it is never imported by ``import lightsim2grid``.
Install it with ``pip install lightsim2grid[torch]``, or just ``pip install torch``.

What is differentiable
----------------------

Inputs are ``(n_scenarios, n_elements)`` tensors; ``None`` means "the grid's own
value, in every row". Row *i* is one independent scenario, with its own
injections and its own set of disconnected elements.

================== ===================================== ==================
input              meaning                               gradient
================== ===================================== ==================
``load_p``         active load, MW                       yes
``load_q``         reactive load, MVAr                   yes
``gen_p``          generator active setpoint, MW         yes
``sgen_p``         static generator active power, MW     yes
``gen_v``          generator voltage setpoint, pu        yes [#genv]_
``line_status``    True = connected                      no (discrete)
``trafo_status``   True = connected                      no (discrete)
``gen_status``     True = connected                      no (discrete)
``shunt_status``   True = connected                      no (discrete)
``load_status``    True = connected                      no (discrete)
================== ===================================== ==================

.. [#genv] Generators regulating one bus must be given the same set-point: they share
   the one derivative that exists (see the module documentation of
   :class:`BatchCPUPowerFlow`).

The statuses follow grid2op's convention (**True means connected**), and a row
that trips a generator gets the whole treatment the C++ already implements: the
lost MW go to the slack, the distributed slack is re-weighted without the
machine, and a bus whose last voltage-regulating generator is gone becomes PQ
for that row. A tripped shunt leaves the row's admittance matrix, a tripped load
its injection; an element a row disconnects has no gradient in that row (its
value there changes nothing).

Why it is not expensive
-----------------------

The gradient comes from the adjoint method (the implicit function theorem), not
from differentiating the Newton-Raphson iterations. Each row costs **one
transposed solve** on top of the forward -- and that cost does not depend on how
many inputs are being differentiated, nor on how wide ``V`` is.

Two things make it cheap here, both of them properties the batch algorithms
already had:

* every row of a batch shares one Jacobian **sparsity pattern**, so the whole
  batch keeps one symbolic factorization;
* KLU answers the transposed system out of the factorization of :math:`J`
  itself (``klu_tsolve``), so no transposed matrix is ever built.

The price is memory: the converged Jacobian of every row is kept, which is
``n_scenarios * nnz(J)`` floats. On a large grid with a large batch that is
gigabytes -- :func:`BatchCPUPowerFlow.adjoint_memory_bytes` reports it, and it
is worth reading before scaling a batch up.

One grid
--------

:class:`GridPowerFlow` is the same operation for one :class:`LSGrid`: its inputs are
``(n_elements,)`` tensors, its forward is the grid's own ``ac_pf`` (on a copy of the
grid), and ``V`` is ``(n_bus,)``.

.. code-block:: python

    from lightsim2grid.differentiable import GridPowerFlow

    pf = GridPowerFlow(grid)                  # the grid's algorithm, and its outer loops
    V = pf(load_p=load_p, gen_v=gen_v)
    loss = pf.compute_flows(V).i_or_a.relu().sum()
    loss.backward()

Nothing is stored by the forward: the backward reads the system the grid's algorithm
still holds. A graph whose forward is not the last one any more has its powerflow run
again first. The same is available without pytorch, on the grid itself:
``LSGrid.solve_JT``, ``get_theta_col_of_bus`` and the other maps,
``gen_v_indirect_grad`` and ``gen_p_indirect_grad``, read after an ``ac_pf`` that
converged, the grid unchanged since.

With outer loops
----------------

With an ``NROuter_*`` or ``NROuterDist_*`` algorithm, the forward runs the outer loops
(per row for a batch), and the gradient is that of the **last Newton solve, every
decision of the loops held as it ended**: a bus they switched PQ stays PQ, a unit held at
a reactive limit stays held, a tap stays where they left it. The same adjoint formula
as without loops, on the Jacobian of the system the loops left, refreshed at the
solution. It is the derivative of the whole run wherever a small enough change of the
inputs changes none of these decisions, which is what the tests check against finite
differences.

What the loops leave is taken into account:

* a bus still holding its voltage is pinned PV in the Jacobian: its ``gen_v`` gradient
  goes through that pinned row, and a bus switched PQ gives its generators no ``gen_v``
  gradient (``get_gen_v_share_rows``, per row);
* a voltage-control group whose controllers are all held at a limit has lost its
  set-point: no ``gen_v`` gradient for its generators;
* the dS/d|V| column a fixed-magnitude bus (the slack) is differentiated with is that of
  the solve's own admittance matrix, at the taps and sections the loops settled on; and
  ``compute_flows`` uses those taps too;
* a unit frozen at a limit of its capability curve injects a reactive power that follows
  its active set-point: ``gen_p_indirect_grad`` adds that slope to its ``gen_p``
  gradient. The curve is read at the target P a batch row or an outer loop gave the unit,
  so on a single grid without a distributed slack its fixed limits apply, and the slope
  is zero;
* ``NROuterDist_*`` shares the slack inside the Newton: its Jacobian carries the
  distributed slack, and the gradient with it. A unit the distributed slack holds at a
  bound gets no ``gen_p`` gradient (its set-point moves nothing); a unit already beyond
  its bound stays at its own set-point, and keeps its gradient.

What it does not cover:

* ``NROuter_*`` with ``DistributedSlack`` moves the units' targets **between** Newton
  solves: its last Jacobian is single-slack, and its gradient puts every change of the
  balance on the reference bus, the shares the loop settled on held fixed. Use
  ``NROuterDist_*`` for the gradient of the distributed slack.
* A value a loop froze from an intermediate solve is held, not differentiated: the
  output of an SVC held at a limit (``B . V^2`` at the voltage it was frozen at), and the
  share of the distributed slack moving the point a capability curve is read at.
* A decision that would flip: the gradient is one-sided information at a point where a
  loop is about to switch, as it is for any piecewise-smooth function.

Rows that fail
--------------

A row that did not converge, or that a contingency islanded, comes back as
``NaN`` rather than as zeros -- zeros are a plausible-looking voltage, ``NaN`` is
not -- and its gradient is zero. :func:`BatchCPUPowerFlow.converged` says which
rows those are, so a loss can mask them out.

Detailed documentation
----------------------

.. automodule:: lightsim2grid.differentiable
    :members:
    :autosummary:

