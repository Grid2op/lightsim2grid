.. _nr_control_limits:

Control limits inside the Newton-Raphson (design note)
=======================================================

.. warning::
    **Nothing described here is implemented.** This is a design note: it records a
    formulation, why it fits lightsim2grid's architecture unusually well, and what it
    would cost. It is written so that the question can be picked up later without
    re-deriving it. The corresponding entry is in ``CHANGELOG.rst``'s ``[TODO]``.

Scope
--------------------------

:ref:`comparison_with_pypowsybl` ends its "outer loops" section on an open question:

    *Whether a given outer loop could instead be folded into the inner NR formulation
    (like distributed slack was) or fundamentally needs iteration around the solve (like
    discrete tap positions, which are not differentiable) is a case-by-case question,
    and remains open for most of OLF's outer loops as far as lightsim2grid is
    concerned.*

This note argues that it is **not** case-by-case: four of the five outer loops named
there are the same mathematical object, one formulation covers all four, and the
remaining one (tap *discreteness*, as opposed to tap *limits*) is the only part that
genuinely needs iteration around the solve.

It also covers two questions that come up alongside: enforcing ``pmin`` / ``pmax`` on
distributed-slack participants, and making the slack weights depend on the state (OLF's
``balanceType``).

.. contents:: Contents
    :local:
    :depth: 2


The observation
--------------------------

Every one of these controls is a **bounded resource, complementary to a regulated
quantity**:

==============================  =====================  ======================  ==============================
control                         resource               box                     regulated quantity
==============================  =====================  ======================  ==============================
distributed slack               generator :math:`P`    ``[pmin, pmax]``        system active power balance
PV :math:`\leftrightarrow` PQ   generator :math:`Q`    ``[min_q, max_q]``      bus voltage magnitude
ratio tap changer               ratio :math:`\rho`     ``[rho_min, rho_max]``  bus voltage magnitude (deadband)
area interchange                area generation        per-machine limits      area net export
secondary voltage control       group :math:`Q`        per-machine limits      pilot bus voltage
==============================  =====================  ======================  ==============================

"Complementary" is meant literally, in the mathematical-programming sense: the resource
is strictly inside its box **if and only if** the regulated quantity is exactly at its
target; a resource pinned at a bound comes with a known *sign* of the regulation error.
That is a complementarity condition, and complementarity conditions can be written as
equations and handed to a Newton solver.

So the alternative to "an outer loop" is not "a hand-rolled in-Jacobian special case per
control". It is **one mechanism, instantiated five times**.


Why the heuristics cycle
--------------------------

The classical PV :math:`\rightarrow` PQ heuristic uses **two different tests**:

- switch PV :math:`\rightarrow` PQ when ``Q > max_q``;
- switch back PQ :math:`\rightarrow` PV when ``Vm > v_set``.

Each is correct in isolation. The defect is that they are evaluated on *different
iterates*, with a discrete bus-type state carried between them. Nothing forbids the
solver satisfying test 1 at iteration :math:`k`, test 2 at :math:`k+1`, test 1 again at
:math:`k+2`, indefinitely. No single iterate is inconsistent; the *sequence* is.

With two coupled controls this gets qualitatively worse, because there is no correct
**order** to run them in. A tap step changes ``Vm``, which un-pins a generator, whose
``Q`` then exceeds ``max_q``, which re-pins it and frees ``Vm``, which moves the tap.
Running the tap logic before or after the reactive-limit logic gives *different* limit
cycles, so the answer depends on an arbitrary implementation choice.

OpenLoadFlow resolves this by fiat, which is the tell: ``transformerVoltageControlMode``
offers ``AFTER_GENERATOR_VOLTAGE_CONTROL`` ("a continuous voltage control is performed
after each generator voltage control outerloop") alongside the default
``CONTINUOUS_WITH_DISCRETISATION``. A parameter that lets the caller choose the
interleaving is an admission that the interleaving matters and that no ordering is
right.

In the formulation below the question **cannot be asked**: every control is a row of the
same Jacobian, solved simultaneously by one Newton step. There is no "after".


The formulation
--------------------------

NCP functions
****************

An **NCP function** is any :math:`\varphi : \mathbb{R}^2 \to \mathbb{R}` with

.. math::
    \varphi(a, b) = 0 \quad \Longleftrightarrow \quad a \geq 0,\; b \geq 0,\; ab = 0

Two matter in practice:

.. math::
    \varphi_{\min}(a,b) &= \min(a, b) \\
    \varphi_{FB}(a,b)   &= a + b - \sqrt{a^2 + b^2} \qquad \text{(Fischer-Burmeister)}

:math:`\varphi_{\min}` is piecewise linear and reproduces exactly the clamp / un-clamp
logic one would write by hand. :math:`\varphi_{FB}` is smooth everywhere *except* at the
single point :math:`(0,0)`, which matters for globalisation (below).

For a **two-sided box** -- which is what every control in the table above has -- the
convenient form is the projection (normal-map) residual. With resource :math:`u`, box
:math:`[u_{\min}, u_{\max}]`, regulation error :math:`e` (target minus actual), and
:math:`\sigma = \operatorname{sign}(\partial(\text{regulated quantity}) / \partial u)`:

.. math::
    \Phi = u - \operatorname{mid}\big(u_{\min},\; u_{\max},\; u + c\,\sigma\, e\big),
    \qquad c > 0

``mid`` is the median of three, i.e. the projection onto the box. Writing
:math:`z = u + c\,\sigma\,e` for the quantity being projected, the three branches are:

.. list-table::
    :header-rows: 1
    :widths: 30 30 40

    * - branch
      - :math:`\Phi` reduces to
      - meaning
    * - :math:`u_{\min} < z < u_{\max}`
      - :math:`-c\,\sigma\,e`
      - **regulating** (:math:`e = 0`)
    * - :math:`z \geq u_{\max}`
      - :math:`u - u_{\max}`
      - pinned at the upper bound
    * - :math:`z \leq u_{\min}`
      - :math:`u - u_{\min}`
      - pinned at the lower bound

**That row is the bus type.** There is no PV / PQ label anywhere in the system.

Why this kills the cycling
*****************************

Take the generator case (:math:`u = Q`, :math:`e = v_{set} - V_m`) and check the two
classical tests against :math:`z = Q + c\,(v_{set} - V_m)`:

- while regulating, :math:`V_m = v_{set}`, so :math:`z = Q` -- pinning when
  :math:`z > Q_{\max}` **is exactly test 1**;
- while pinned, :math:`Q = Q_{\max}`, so :math:`z = Q_{\max} + c(v_{set} - V_m)` --
  releasing when :math:`z < Q_{\max}` **is exactly test 2**.

The projection reproduces both classical tests, exactly, but as **one function of one
iterate**. There is no state carried between tests to become inconsistent. Add the line
search below and the iteration is monotone in a merit function, so it cannot return to a
state it has left. Cycling is not suppressed by a latch or an iteration counter -- it
becomes structurally impossible.

Semismooth Newton
********************

:math:`\Phi` is not differentiable, so ordinary Newton does not apply. **Semismooth
Newton** does.

A locally Lipschitz, directionally differentiable :math:`F` is *semismooth at* :math:`x`
when, for **every** :math:`V` in the Clarke generalized Jacobian
:math:`\partial F(x + d)`,

.. math::
    F(x + d) - F(x) - V d = o(\lVert d \rVert) \qquad (d \to 0)

and *strongly* semismooth when the right-hand side is :math:`O(\lVert d \rVert^2)`. The
quantifier is the content: the bound must hold for every choice of generalized Jacobian
element, which is what lets an implementation pick an arbitrary one.

The condition is weak in practice. Piecewise smooth functions are semismooth;
:math:`\min`, :math:`\max`, :math:`\lvert\cdot\rvert`, ``mid`` and the Euclidean norm are
*strongly* semismooth; sums and compositions preserve it. So the full residual -- smooth
power-flow equations plus projection rows -- is strongly semismooth.

**Theorem (Qi & Sun 1993).** If :math:`F(x^\star) = 0`, :math:`F` is semismooth at
:math:`x^\star`, and every :math:`V \in \partial_B F(x^\star)` is nonsingular
(*BD-regularity*), then

.. math::
    V_k &\in \partial_B F(x_k) \quad \text{(any element)} \\
    V_k\, d_k &= -F(x_k) \\
    x_{k+1} &= x_k + d_k

converges Q-superlinearly from any :math:`x_0` near :math:`x^\star`, and Q-quadratically
if :math:`F` is strongly semismooth. One linear solve per iteration -- the same cost as
ordinary Newton.

Globalisation
****************

Fischer's result (1992): :math:`\psi_{FB} = \tfrac{1}{2}\varphi_{FB}^2` is
**continuously differentiable on all of** :math:`\mathbb{R}^2`, including at the origin
where :math:`\varphi_{FB}` is not. So

.. math::
    \Psi(x) = \tfrac{1}{2}\lVert \Phi(x) \rVert^2

is :math:`C^1` even though :math:`\Phi` is not, and a standard Armijo backtracking line
search on :math:`\Psi` applies to the semismooth Newton direction. That is the
globalisation, and it is what makes cycling impossible rather than merely unlikely.

.. note::
    :math:`\Psi` can have stationary points that are not solutions. The complementarity
    literature rules these out by assuming the underlying map is :math:`P_0`, which AC
    power flow is not. In practice: line search, with a fallback to an undamped
    semismooth step when it stalls. This is the same posture as the step damping the
    Newton-Raphson already uses.


Instantiation per control
--------------------------

Distributed slack with ``pmin`` / ``pmax``
********************************************

Today ``MultiSlack`` (``src/core/powerflow_algorithm/NRSystem.hpp``) is the textbook
bordered formulation: one scalar unknown ``slack_absorbed`` (call it :math:`k`), one
column, one P equation per slack bus, and feature entries :math:`w_b` at each slack bus'
P row. With the sign convention in the code (``mis += k * slack_weights``, see
``_evaluate_Fx`` in ``NRSystem.tpp`` and ``LSGrid::compute_results``), the effective
injection is :math:`S_{bus} - k\,w`, i.e. **:math:`k` is negative when generation must
rise**. Worth a comment in the code -- it reads backwards.

Add one non-negative *shed* variable :math:`s_g` per slack generator:

.. math::
    P_g &= P_g^0 - w_g k - s_g \\
    0 &\leq s_g \;\perp\; (P^{\max}_g - P_g) \geq 0

Semantics: :math:`s_g = 0` with positive headroom is a machine following the slack;
:math:`s_g > 0` with zero headroom is a machine pinned at ``pmax`` while the rest of the
fleet takes its share. A symmetric block handles ``pmin``.

:math:`\Phi_g` depends only on :math:`s_g` and :math:`k`, so each row has exactly **two**
nonzeros. At bus level ``adjust_mismatch`` gains one term, ``mis += k*w + E*s``, with
:math:`E` the generator-to-bus incidence (one nonzero per column).

.. warning::
    ``GeneratorContainer`` has ``min_q_`` / ``max_q_`` but **no active power limits at
    all** -- see the ``// TODO add pmin and pmax here !`` in ``GeneratorContainer.hpp``.
    That is prerequisite zero for any of this.

Two further points specific to the slack:

- **Granularity.** Limits are per generator; the Jacobian weight is per *bus*
  (``get_slack_weights_solver``), and the generator-level split is post-hoc and
  proportional (``set_p_slack``). With several slack generators on one bus, pinning one
  means recomputing :math:`w_b` from the survivors -- so ``bus_slack_weight_`` becomes an
  NR-iteration quantity, not a per-solve one.
- **Keep-one rule.** If every slack generator pins, the :math:`k` column goes to zero and
  the bordered system is singular. BD-regularity failing *is* the infeasibility detector:
  report it rather than diverging.

State-dependent slack weights (OLF ``balanceType``)
*****************************************************

Two different things get conflated here. Using :math:`\kappa = -k` so that
:math:`\kappa > 0` means "generation must rise":

**Direction-dependent weights** -- different participating sets and factors for up- and
down-regulation -- are a kink at :math:`\kappa = 0` in a function of an unknown.
Piecewise smooth, hence semismooth: evaluate :math:`w` at the current iterate, take
either branch at the kink. No new unknowns, no new sparsity. Essentially free once the
rest exists.

**Margin-proportional weights** (OLF's ``PROPORTIONAL_TO_GENERATION_REMAINING_MARGIN``)
are implicit: :math:`w_g \propto P^{\max}_g - P_g`, and :math:`P_g` depends on
:math:`w_g`. OLF resolves it by freezing :math:`w` from the current :math:`P` at the top
of each outer iteration. In-NR, make the margin an unknown with its own local row:

.. math::
    P_g &= P_g^0 + \kappa\, m_g \\
    m_g &= \max\big(0,\; P^{\max}_g - P_g\big)

Two nonzeros per row. Solving the row on the active branch:

.. math::
    m_g = \frac{P^{\max}_g - P_g^0}{1 + \kappa}
    \qquad\Longrightarrow\qquad
    P_g = P_g^0 + \kappa\,\frac{P^{\max}_g - P_g^0}{1 + \kappa}

which is monotone in :math:`\kappa` and tends to :math:`P^{\max}_g` **without ever
exceeding it**. So margin-proportional weighting *is* the "reduce the weight as the
machine approaches pmax" idea, and unlike a fixed-steepness sigmoid it enforces the limit
exactly, because the taper is self-consistent rather than imposed. No complementarity
needed for this variant.

Caveats, in decreasing order of nuisance:

- **Conditioning.** Total reachable imbalance is
  :math:`\sum_g (P^{\max}_g - P^0_g)\cdot\kappa/(1+\kappa)`, so exhausting the fleet's
  headroom needs :math:`\kappa \to \infty`. An infeasible case diverges in :math:`\kappa`
  rather than failing cleanly: cap it explicitly.
- **It is a different physical model.** Margin-proportional means every machine
  approaches its limit asymptotically *together*; AGC participation factors mean machines
  pin *one at a time*. Both are defensible and they give different answers. Do not let
  "easier to formulate" pick the model.
- :math:`\max(0,\cdot)` still matters for a machine starting above ``pmax``.

.. warning::
    **Do not normalise inside the NR.** ``GeneratorContainer::get_slack_weights_solver``
    does ``res /= sum_res``. Carried into a state-dependent weight,
    :math:`w_g = m_g / \sum_h m_h` makes every :math:`\partial w_g/\partial m_h` nonzero
    and the slack border **dense** across all participants -- which destroys the batch
    algorithms' premise.

    The normalisation is a *gauge*: :math:`\kappa \cdot m_g/\sum m` is just
    :math:`\kappa' m_g` with :math:`\kappa' = \kappa/\sum m`, and :math:`\kappa` is a free
    unknown, so the absolute scale of :math:`w` is unobservable inside the NR -- only the
    ratios are. Drop it from the solver and keep it only in ``set_p_slack``'s proportional
    share, where the sum is a constant computed after the fact.

Hydro with produce / absorb modes
************************************

A machine with a "produce" mode (:math:`0 < P < p_{\max}`) and an "absorb" mode
(:math:`p_{\min} < P < 0`), where switching mode *during* a solve is not wanted.

**If the mode is frozen** (decided by the dispatch, which is the physically defensible
position -- whether a hydro unit turbines or pumps is a commitment decision, not
something a power flow should decide), then the mode boundary is **just another limit**:
the box is ``[0, pmax]`` in produce mode and ``[pmin, 0]`` in absorb mode. A produce-mode
machine that the slack wants to push negative simply pins at :math:`P = 0`, its weight
goes to zero, and its share is redistributed. The "don't switch" guarantee is not a guard
bolted on; it is a property of the system being solved.

**If the solve should pick the mode**, split the injection:
:math:`P = p^+ - p^-` with :math:`0 \leq p^+ \perp p^- \geq 0`. Representable -- but this
is an MPEC-style complementarity between two *decision* variables, degenerate exactly at
:math:`P = 0` where both vanish, constraint qualifications fail generically and
BD-regularity with them. Machines idling near zero will be where the solver sticks.
Representable, yes; well behaved, no.

**If the mode has a true deadband** (produce mode really starts at some
:math:`p_{\min}^{turbine} > 0`), the feasible set has a *hole* and is not
NCP-representable at all. Freeze the mode, treat the deadband edge as the machine's
limit, and report a commitment infeasibility upward rather than letting the power flow
shut the unit down.

PV :math:`\leftrightarrow` PQ reactive limits
***********************************************

This is the **original** application of the technique -- every reference at the end of
this note is a reactive-limit paper. With resource :math:`Q_j`, box
``[min_q, max_q]``, error :math:`e = v_{set} - V_{m,r}`:

.. math::
    \Phi_j = Q_j - \operatorname{mid}\big(Q^{\min}_j,\; Q^{\max}_j,\;
             Q_j + c\,(v_{set} - V_{m,r})\big)

For a group of :math:`N` machines sharing :math:`Q` proportionally, the sharing rows take
the same shed treatment as the distributed slack:
:math:`Q_j = \alpha_j Q_{grp} - t_j` with :math:`0 \leq t_j \perp (Q^{\max}_j - Q_j) \geq 0`,
and the group's voltage row regulates while any member has headroom, releasing
:math:`V_{m,r}` once all are pinned.

Ratio tap changers
********************

Identical structure, with :math:`u = \rho`:

.. math::
    \Phi = \rho - \operatorname{mid}\big(\rho_{\min},\; \rho_{\max},\;
           \rho + c\,\sigma\,(v_{set} - V_{m,t})\big)

.. warning::
    :math:`\sigma = \operatorname{sign}(\partial V_{m,t} / \partial \rho)` is **not
    optional and not guessable**: it depends on the ratio convention and on
    ``is_tap_side1`` (``TrafoContainer.hpp``). Getting it backwards selects the *opposite*
    branch -- the tap runs to the wrong rail and the solve converges, silently, to a wrong
    answer. Compute the sensitivity from the model; do not hard-code a sign. This is the
    class of bug the ``LocalBusId`` / ``GlobalBusId`` / ``SolverBusId`` tagging exists to
    prevent elsewhere in this codebase.

**The deadband makes this better posed, not worse.** A real OLTC regulates to a band, not
a point, and a point target is the ill-posed idealisation -- the device moves in steps and
can never sit exactly at :math:`v_{set}`. The band is also where an otherwise
structurally singular :math:`\rho` column gets an equation, and the device's own answer
("don't move unless the voltage leaves the band") is a perfectly good one:

.. list-table::
    :widths: 40 60

    * - :math:`V_{m,t} < v_{lo}`
      - :math:`\rho` rises, up to :math:`\rho_{\max}`
    * - :math:`V_{m,t} > v_{hi}`
      - :math:`\rho` falls, down to :math:`\rho_{\min}`
    * - :math:`v_{lo} \leq V_{m,t} \leq v_{hi}`
      - :math:`\Phi = \rho - \rho_0`, i.e. pin at the incoming tap

Five branches instead of three, still piecewise linear, still strongly semismooth -- and
the deadband branch gives the iteration somewhere to *rest*, which is exactly what the
cycling heuristic lacks.

**Discreteness is the one genuine wall.** Complementarity handles *bounds*, not
*integrality*; no NCP function represents a lattice. The practical answer is the
universal one: solve the continuous relaxation, round to the nearest position, re-solve
once with taps fixed -- which is what both OLF modes do.

Here lightsim2grid is unusually well placed: ``TrafoContainer.hpp`` states that "there is
NO notion of a discrete ``tap`` in lightsim2grid" -- the state is a continuous ``ratio_``
and the tap-position arithmetic lives in the converters. So the continuous relaxation is
already the native representation, and the rounding naturally belongs upstream in the
pypowsybl path where tap positions actually exist.

Area interchange and secondary voltage control
*************************************************

Sketched only, because they add no new mathematics:

- **Area interchange**: one shared scalar per area (as :math:`k` is for the system),
  participating machines with limits, regulated quantity = the area's net export.
  Structurally the distributed slack, replicated per area.
- **Secondary voltage control**: a pilot bus regulated by a group of machines sharing
  :math:`Q` by participation, each with reactive limits. Structurally the
  ``VoltageControl`` group, with the shed variables of the reactive-limit section.


Fit with lightsim2grid's architecture
----------------------------------------

The batch premise stops being a constraint
*********************************************

``CLAUDE.md`` states the rule: *"what is genuinely forbidden is changing the pattern
itself per row -- a different pv/pq vector handed to the solver, a different slack set.
That raises* ``has_pv_changed()`` */* ``has_slack_participate_changed()``\ *, forces a
fresh* ``analyze`` *on every row, and the batch loses its reason to exist."*
``set_switchable_vm_buses`` + ``set_pv_pinned_buses`` exist to work around exactly that:
reserve the union of layouts, mask the difference by value.

In this formulation **there is no pv/pq vector**. Every regulated bus permanently owns a
:math:`V_m` unknown and a :math:`Q` equation; every controller permanently owns a
:math:`Q` unknown; "PV-ness" is a *value* of the projection row, recomputed each
iteration. So:

- ``set_switchable_vm_buses`` becomes unnecessary -- every regulated bus is switchable by
  construction;
- ``set_pv_pinned_buses`` becomes unnecessary -- pinning is what the projection row does;
- ``has_pv_changed()`` can never fire from a reactive-limit event, because there is
  nothing to change.

The workaround stops being a workaround and becomes the formulation. Everything that
changes when a control pins or releases is a **value** (:math:`\varphi_a`,
:math:`\varphi_b`); positions are fixed at ``register_in``. One ``analyze``,
``refactorize`` thereafter, for every row of a sweep. Turn ``set_refactor_fallback`` on
for the usual reason -- those value edits move pivots.

``VoltageControl`` already has most of the slots
***************************************************

Reading ``NRSystem.hpp``, the reactive-limit case needs strikingly little new structure.
``VoltageControl::register_in`` already claims, per group:

- ``q_cols_[j] = ledger.add_q_unknown(...)`` -- :math:`Q` is **already an explicit
  unknown**;
- ``v_rows_[g] = ledger.add_custom_row()`` -- the voltage constraint already owns a row;
- ``h_vm_[g] = sink.add(v_row, vm_cols_[g])`` and
  ``h_slope_[j] = sink.add(v_row, q_cols_[j])`` -- **both Jacobian slots the three-way
  branch needs are already reserved**.

And the ``stranded`` path is a hand-rolled instance of exactly the row swap this needs.
``fill_feature_values`` writes ``h_slope_[j] <- stranded ? 1. : slope(j)`` and
``h_vm_[g] <- stranded ? 0. : 1.``, while ``fill_custom_rows`` swaps the residual from
``Vm + sum s.Q - v_set`` to ``Q_c``. That is the voltage row being turned into a
"``Q_c = 0``" pin **by value alone, inside a fixed sparsity pattern**. Replace ``0`` with
``max_q`` and make the branch selector the projection test instead of the stranded flag,
and the formulation is there.

For a lone controller this costs **zero new rows, zero new columns, zero new nonzeros** --
the slot exists. The only gating change is the condition that currently reserves
``(v_row, q_col)`` for SVCs and for ``may_mask_`` singletons only; it would become
unconditional.

``Hvdc`` is the precedent for a state-dependent ``Ybus``
***********************************************************

A continuous :math:`\rho` makes the transformer's ``Ybus`` entries state-dependent, and
``NRSystem``'s generic dS pass assumes a fixed ``Ybus``. Two ways out:

1. recompute the 2x2 block each iteration and raise ``tell_recompute_ybus()`` -- correct,
   but it fights the whole ``AlgoControl`` caching design;
2. **stamp the transformer at a nominal ratio** :math:`\rho_0` **and let an extension own
   the difference** as an injection correction,
   :math:`\Delta S = V \odot \overline{(Y(\rho) - Y(\rho_0)) V}`, applied in
   ``adjust_mismatch``.

Option 2 is not speculative: it is what ``Hvdc`` already does. ``HvdcLineContainer.hpp``
documents it -- *"For droop lines, this container does NOT stamp the AC active power in
Sbus: the NR system Hvdc extension owns it (all regimes)."* Same trick one level up. The
container stamps the constant part, the extension owns the state-dependent part, ``Ybus``
stays constant and its sparsity untouched, and the correction's Jacobian contributions
land in dS positions that already exist because the branch is already there at
:math:`\rho_0`.

Per regulating transformer: **1 column** (:math:`\rho`), **1 custom row** (2 nonzeros),
up to 4 entries in the :math:`\rho` column (P and Q rows of both ends). Negligible.

A ``TapControl`` extension beside ``MultiSlack``, ``Hvdc`` and ``VoltageControl`` is the
obvious home -- the component protocol documented at the top of ``NRSystem.hpp``
(``update_state`` / ``init_topology`` / ``register_in`` / ``declare_feature_entries`` /
``fill_feature_values`` / ``adjust_mismatch`` / ``fill_custom_rows`` / ``apply_step`` /
``clear``) already has every hook this needs.


Costs and open questions
--------------------------

Stated plainly, because they are the reason this is not obviously the right thing to do.

- **Dimension growth on the common grid.** Today a PV bus has *no* :math:`V_m` column and
  *no* :math:`Q` row. Here it has both, plus a :math:`Q` unknown and a voltage row. On a
  transmission grid with ~10-15% regulated buses, ``J`` grows meaningfully in dimension
  and more in nonzeros (the new :math:`Q` rows drag in ``dS_dVm`` entries). This is paid
  on **every** solve, including the majority that never touch a limit.

  Note where it lands, though: the layout is *identical to the union layout already built*
  when ``set_switchable_vm_buses`` is populated. So the cost is real for a plain
  ``ac_pf`` and **zero for the batch algorithms** -- which is also where a cycling
  heuristic hurts most (a switching loop on row 4000 of a ``ScenarioSweep`` costs
  iterations across the whole sweep and is miserable to debug). Keeping the classical
  pv/pq split for ``ac_pf`` and using projection rows in the batch path is a defensible
  split.

- **Scaling.** The constant :math:`c` mixes two physical quantities. ``max_q - min_q``
  can be 10 pu while the voltage deviation of interest is 0.05 pu; per-machine
  normalisation (:math:`\varphi(a/\sigma_g, b/\sigma_g)`, :math:`\sigma_g` a rating or a
  span) is not optional. Textbook implementations omit this and fail on real grids for
  this reason alone. The same trap appears once per control type.

- **Degeneracy.** BD-regularity can fail exactly at a bound with zero slack -- a machine
  dispatched at ``max_q`` with :math:`V_m` exactly at setpoint, which is not exotic. The
  cure is Kanzow's smoothing,
  :math:`\varphi_\mu(a,b) = a + b - \sqrt{a^2+b^2+2\mu^2}`, solved for a sequence
  :math:`\mu_j \to 0` with warm starts. Note this *is* the sigmoid-smoothing idea, but
  with a homotopy that drives the smoothing to zero, so the converged answer respects the
  limit exactly instead of landing at an arbitrary point inside it.

- **Non-uniqueness.** Complementarity fixes the *path*, not the *solution set*. AC power
  flow with voltage controls still has multiple solutions and one can still converge to a
  low-voltage one.

- **Discreteness.** Taps and switched shunts still need round-and-re-solve. This is the
  only part of the outer-loop question that genuinely needs iteration around the solve.

- **Unvalidated at scale, here.** None of this has been measured on this codebase. The
  claims about iteration counts and about batch-path cost are structural arguments, not
  benchmarks.


A possible incremental path
------------------------------

Each step is independently useful and independently testable. They are the same algorithm
at rising levels of rigour, so this is a climb rather than a commitment.

1. Add ``pmin`` / ``pmax`` to ``GeneratorContainer`` (prerequisite for anything on the
   slack side; see the existing ``TODO`` there).
2. Implement discrete clamp-and-renormalise for slack ``pmax``, reserving the shed
   rows / columns up front. This *is* :math:`\varphi_{\min}` semismooth Newton with a unit
   step and no merit function. ``BaseBatchSweep::_masked_slack_weights`` already does the
   zero-and-renormalise on :math:`w`; reuse it per NR iteration instead of per row.
3. Add the Armijo line search on :math:`\Psi = \tfrac{1}{2}\lVert\Phi\rVert^2`. Cycling
   disappears. This is now globalised semismooth Newton.
4. Apply the projection row to reactive limits, gated to the batch path first (where the
   layout is already reserved). Flip the ``(v_row, q_col)`` reservation to unconditional.
5. Swap :math:`\varphi_{\min} \to \varphi_{FB}` where degenerate machines show up, with a
   :math:`\mu`-homotopy if needed.
6. ``TapControl`` extension, continuous :math:`\rho` with deadband, discretisation left to
   the converter layer.

Steps 1-3 alone close the ``pmin`` / ``pmax`` gap on the distributed slack, which is the
smallest useful deliverable.


References
--------------------------

Nonsmooth analysis and semismooth Newton:

- R. Mifflin, *Semismooth and semiconvex functions in constrained optimization*, SIAM J.
  Control Optim. 15 (1977) -- the original definition.
- L. Qi & J. Sun, *A nonsmooth version of Newton's method*, Math. Programming 58 (1993)
  353-367 -- the convergence theorem.
- A. Fischer, *A special Newton-type optimization method*, Optimization 24 (1992) 269-284
  -- the FB function and its :math:`C^1` merit function.
- T. De Luca, F. Facchinei & C. Kanzow, *A semismooth equation approach to the solution of
  nonlinear complementarity problems*, Math. Programming 75 (1996) 407-439 -- the
  practical algorithm, line search included.
- F. Facchinei & J.-S. Pang, *Finite-Dimensional Variational Inequalities and
  Complementarity Problems*, Springer 2003, vol. II ch. 7-9 -- the reference text.

Power systems:

- W. Murray, T. Tinoco De Rubira & A. Wigington, *A robust and informative method for
  solving large-scale power flow problems*, Comput. Optim. Appl. 62 (2015)
  -- complementarity for reactive limits and PV/PQ switching; the GRIDOPT / PFNET line of
  work.
- *A modified Newton-Raphson load flow scheme for directly including generator reactive
  power limits using a complementarity framework*, Electric Power Systems Research 105
  (2013).
- S. V. Dhople, Y. C. Chen, A. Al-Digs & A. Dominguez-Garcia, *Reexamining the Distributed
  Slack Bus*, IEEE Trans. Power Systems 35(6) (2020) 4870-4879 -- what :math:`P^0`,
  :math:`\alpha` and the imbalance are *supposed* to mean. Relevant because
  ``GeneratorContainer`` currently derives weights from ``|target_p|``, which is a
  heuristic and not an AGC participation factor: enforcing limits on top of arbitrary
  participation makes the *limits* exact and the *distribution* still arbitrary.
- A. Pandey, M. Jereminov, G. Hug & L. Pileggi, *Continuously Differentiable Analytical
  Models for Implicit Control within Power Flow*, arXiv:1811.02000 (2018), and
  A. Agarwal, A. Pandey, M. Jereminov & L. Pileggi, *Implicitly Modeling Frequency Control
  within Power Flow*, ISGT-Europe 2019 (arXiv:1908.11778) -- the sigmoid-smoothing
  alternative, motivated explicitly by removing the outer loop.

.. seealso::
    :ref:`comparison_with_pypowsybl` for the outer-loop architecture this note responds
    to, and ``examples/dist_slack_algorithm/`` for the opposite approach (an explicit
    OLF-style outer loop around a single-slack inner solve) expressed as a solver plugin.
