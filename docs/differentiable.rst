Differentiable powerflow (pytorch)
===================================

.. _differentiable:

A batch of powerflows, differentiable with respect to its injections, as a
pytorch operation: a loss computed on the results of thousands of AC powerflows
can be back-propagated to the loads and generators that produced them.

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
``gen_v``          generator voltage setpoint, pu        not yet [#genv]_
``line_status``    True = connected                      no (discrete)
``trafo_status``   True = connected                      no (discrete)
``gen_status``     True = connected                      no (discrete)
================== ===================================== ==================

.. [#genv] ``gen_v`` can be used as a fixed setpoint; a tensor that requires a
   gradient is rejected rather than silently given none. Its gradient needs a
   term the Jacobian does not store for a voltage-fixed bus.

The statuses follow grid2op's convention (**True means connected**), and a row
that trips a generator gets the whole treatment the C++ already implements: the
lost MW go to the slack, the distributed slack is re-weighted without the
machine, and a bus whose last voltage-regulating generator is gone becomes PQ
for that row.

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

