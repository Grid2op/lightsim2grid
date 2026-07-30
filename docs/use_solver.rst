.. _use-solver:

Use as Pandapower Solver
=========================

Advanced usage
-----------------

.. versionchanged:: 0.6.0

    As of version 0.6.0 lightsim2grid implements the new API of "newtonpf" required by pandapower. This means that
    it is asked specifically to provide `ref` a vector of integer identifying the slack buses (implements a distributed
    slack)

LightSim2grid can be used as a specific implementation of the pandapower "newtonpf" function.

Suppose you somehow get:

- `Ybus` the admittance matrix of your powersystem, for example given by pandapower
  (will be converted to a scipy `sparse.csc_matrix` )
- `V0` the (complex) voltage vector at each bus, for example given by pandapower
- `Sbus` the (complex) power absorb at each bus, for example as given by pandapower
- `ref` Ids of the slack buses (added in version 0.5.6 to match recent pandapower changes)
- `pv` list of PV buses
- `pq` list of PQ buses
- `ppci` a ppc internal pandapower test case (or dictionary, is used to retrieve the coefficients associated to each slack bus)
- `options` list of pandapower "options" (or dictionary with keys `max_iteration` and `tolerance_mva`)

You can define replace the `newtonpf` function of `pandapower.pandapower.newtonpf` function with the following
piece of code:

.. code-block:: python

    from lightsim2grid.newtonpf import newtonpf

    # when pandapower version <= 2.7.0
    # V, converged, iterations, J, Vm_it, Va_it = newtonpf(Ybus, Sbus, V0, pv, pq, ppci, options)

    # when pandapower version > 2.7.0
    V, converged, iterations, J, Vm_it, Va_it = newtonpf(Ybus, Sbus, V0, ref, pv, pq, ppci, options)

This function uses the KLU algorithm (or the solver provided in Eigen if KLU has not been installed)
and a c++ implementation of a Newton solver for speed.

.. note::

  The oldest `newtonpf` function compatible with older version of pandapower (*eg* <=2.6.0) can still be accessed with
  `from lightsim2grid.newtonpf import newtonpf_old`

.. note::

  `lightsim2grid.newtonpf` is a thin re-export of `lightsim2grid.pandapower_compat.newtonpf`, which is
  its actual home. `lightsim2grid.pandapower_compat` also provides the DC counterpart, `dcpf`, a
  drop-in replacement for pandapower's own `dcpf` function (same use case as `newtonpf` above, but
  for DC powerflows): `from lightsim2grid.pandapower_compat import dcpf`.

Even more advanced usage
--------------------------
You can customize even more the solvers that you want to use.

Lightsim2grid comes with 22 available solvers (see :class:`lightsim2grid.algorithm.AlgorithmType`, which also has a
``Custom`` sentinel value used for solvers loaded via a plugin, see :ref:`solver_plugin`) that can solve either AC or
DC powerflows. We can cluster them into 4 categories.

If you want to stay "relatively high level", once you have a grid model, you can change the solver using
the "enum" of the solvers you want to use as showed bellow:

.. code-block:: python

    from lightsim2grid.algorithm import AlgorithmType
    # init the grid model
    from lightsim2grid.network import init_from_pandapower
    pp_net = ...  # any pandapower grid
    lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

    # change the solver used for the powerflow
    lightsim_grid_model.change_algorithm(AlgorithmType.NR_KLU)  # change the NR solver that uses KLU
    # you can replace `AlgorithmType.NR_KLU` by any of the available algorithms described bellow,
    # for example `AlgorithmType.NRSing_KLU`, `AlgorithmType.NR_SparseLU`
    # or even `AlgorithmType.NR_NICSLU`

All algorithms can be accessed with the same API (if you want to use the raw python class, not recommended):

.. code-block:: python

  from lightsim2grid.algorithm import NR_KLU  # or any of the names above; here NR_KLU is used as example
  Ybus = ...  # a csc sparse matrix (it's really important that it is a csc and not a csr !)
  V0 = ...  # a complex vector (initial guess)
  Sbus = ...  # a complex vector (power injected at each bus)
  ref = ...  # a vector of integer giving the bus id of the slack buses
  slack_weight = ...  # a vector of real number giving the weights associated to each slack bus
  pv = ...  # a vector of integers giving the bus id of the PV bus
  pq = ...  # a vector of integers giving the bus id of the PQ bus
  max_it = ...  # a > 0 integer maximum number of iterations the solver is allowed to perform
  tol = ...  # a > 0. real number giving the maximum KCL violation allowed for a all nodes

  solver = NR_KLU()
  converged = solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)

  # to retrieve the voltages related information (in case converged is True)
  solver.get_Va()  # voltage angle
  solver.get_Vm()  # voltage magnitude
  solver.get_V()  # complex voltage
  # for compatible solvers
  solver.get_J()  # see documentation of the `newtonpf` function for more information about the shape of J.

  # since lightsim2grid 1.0.0, the Jacobian rows / columns no longer follow the
  # pandapower convention. To know which column of `J` holds which unknown, use
  # the mappings below. Each returns a vector **indexed by the (solver) bus id**
  # -- ie the same indexing as `Ybus` / `V0` above -- giving the Jacobian column
  # of that bus' unknown, or -1 when the bus has no such unknown (eg. the
  # voltage-angle column of a slack bus, or the voltage-magnitude of a PV bus):
  theta_col = solver.get_theta_to_J_col()  # bus id -> column of its voltage-angle (theta) unknown
  vm_col = solver.get_vm_to_J_col()         # bus id -> column of its voltage-magnitude (vm) unknown
  q_col = solver.get_q_to_J_col()           # bus id -> column of its reactive (q) unknown
  # eg. the angle unknown of (solver) bus 4 is column `theta_col[4]` of `J`

  # the rows of `J` follow the same layout as its columns: the first rows are the
  # active power (P) mismatch equations, indexed exactly like the theta / vm columns
  # above (one row per pvpq bus), followed by the reactive power (Q) mismatch
  # equations, one row per pq bus. So `J[theta_col[4], :]` is the P-mismatch equation
  # of (solver) bus 4, and `J[q_col[4], :]` (when not -1) is its Q-mismatch equation.
  # See the ``J_description`` note of :func:`lightsim2grid.network.LSGrid.get_J_solver`
  # for the full block structure (including how the distributed-slack coupling row /
  # column is inserted).

  # all the ids above are in the *solver* bus numbering (see :ref:`bus-labelling`),
  # which is compact (0 ... nb_connected_bus() - 1) and can change with the topology.
  # To map a solver bus id back to the stable GridModel (global) bus id -- eg to know
  # which substation / bus a given row or column of `J` actually corresponds to --
  # use :func:`~lightsim2grid.network.LSGrid.id_ac_solver_to_me` (and its converse
  # :func:`~lightsim2grid.network.LSGrid.id_me_to_ac_solver`) on the `LSGrid` object
  # used to build the solver (not on the raw solver object itself):
  #
  #     lightsim_grid_model.id_ac_solver_to_me()
  #
  # eg. the angle unknown of (solver) bus 4 above is for GridModel bus
  # `lightsim_grid_model.id_ac_solver_to_me()[4]`

  # some other usefull information
  solver.get_nb_iter()  # return the number of iteration performed
  solver.get_timers()  # timer_Fx_ / timer_solve_ / timer_check_ (seconds spent in each stage)
  solver.get_error()  # an `ErrorType` value, eg ErrorType.NoError, .SingularMatrix, .TooManyIterations, ...
  solver.converged()  # equal to the boolean `converged` above

``solver.solve(...)`` validates most of its inputs and raises a regular python exception
(``RuntimeError`` or ``IndexError``) instead of crashing when they are not met. lightsim2grid
expects:

- `tol` to be a finite, strictly positive number (raises otherwise).
- `max_it` to be a non-negative integer (raises otherwise; `max_it = 0` is accepted -- it just
  means the solver stops immediately without converging).
- `Ybus` to be a square matrix, with the same size as `Sbus` and `V0` (raises otherwise). `Ybus`
  no longer needs to specifically be a `scipy.sparse.csc_matrix`: any 2D sparse format accepted
  by `scipy.sparse` is converted internally, so this is only a (minor) performance consideration,
  not a correctness one, if you already have it in CSC format.
- every id appearing in `ref`, `pv` or `pq` to be a valid bus id (in ``[0, len(V0))``), and
  `ref`, `pv` and `pq` to be pairwise disjoint (raises ``IndexError`` / ``RuntimeError`` otherwise).
- `ref` (the list of slack buses) to be non-empty (raises otherwise).
- `slack_weight` to have the same size as `Ybus` (raises otherwise).

The following are **not** validated, and can silently produce an incorrect (but not crashing)
result if violated:

- all buses should appear in the union `ref U pv U pq` (a bus missing from all three is silently
  dropped from the system of equations instead of raising an error).
- `slack_weight[node id] > 0.` for all node id in `ref`, and `sum(slack_weight) = 1.` (an
  inconsistent `slack_weight` will not raise, it will just distribute the slack power
  differently than you might expect).

.. warning::

    Because the two conditions above are not checked, we still do not recommend using these
    raw solver classes directly unless you know exactly what you are doing: prefer building
    the grid model from an :func:`~lightsim2grid.network.init_from_pandapower` (or any other
    ``init_from_*``) call and letting :class:`~lightsim2grid.network.LSGrid` /
    :class:`~lightsim2grid.lightSimBackend.LightSimBackend` derive `ref` / `pv` / `pq` /
    `slack_weight` for you.


AC solvers using Newton Raphson
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are 8 solvers in this categorie. They can in turn, be split into two main sub categories. The first one (`NR_*`) allows
for a distributed slack bus (but can be a bit slower); the other one (`NRSing_*`) does not distribute the slack according to
`slack_weight` at all.

.. warning::

   For `NRSing_*`, if `ref` contains more than one bus id, they are **not** converted to PV buses (unlike what a
   quick reading of "single slack" might suggest): every bus in `ref` keeps both its voltage angle **and**
   magnitude fixed (exactly like a lone slack bus would), `slack_weight` is entirely ignored, and there is no
   equation coupling the buses in `ref` together. In practice this means each bus in `ref` independently absorbs
   whatever active / reactive mismatch shows up locally at that bus -- it behaves like an (uncoordinated)
   multi-reference powerflow, not like a single-slack powerflow with the extra buses silently downgraded to PV.
   Verified empirically: with a second bus added to `ref` for a `NRSing_KLU` solve, that bus' resulting voltage
   angle is exactly ``0.`` (the same as the "real" slack bus), whereas the Gauss-Seidel and Fast-Decoupled solvers
   below, running on the very same grid, converge with that second bus' angle at a genuine, non-zero, solved
   value (*ie* they do actually convert it to a PV bus). If you rely on "single slack" behaviour, make sure `ref`
   truly contains a single bus id -- see the general recommendation at the end of the Gauss-Seidel section below
   about removing all but one generator from the slack in the first place.

The list is:

- `NR_KLU` \*: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the
  fast `KLU` implementation is used to solve the linear system `J.dx = mismatch` at each iteration.
- `NR_NICSLU` \*: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the
  fast `NICSLU` implementation is used to solve the linear system `J.dx = mismatch` at each iteration.
- `NR_CKTSO` \*: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the
  fast `CKTSO` implementation is used to solve the linear system `J.dx = mismatch` at each iteration.
- `NR_SparseLU`: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the
  Eigen default implementation is used to solve the linear system `J.dx = mismatch` at each iteration (instead of the faster `KLU`, `NICSLU` or `CKTSO`)
- `NRSing_KLU` \*: implementation of the Newton Raphson algorithm ignoring `slack_weight` (see warning above), where the
  fast `KLU` implementation is used to solve the linear system `J.dx = mismatch` at each iteration
- `NRSing_NICSLU` \*: implementation of the Newton Raphson algorithm ignoring `slack_weight` (see warning above), where the
  fast `NICSLU` implementation is used to solve the linear system `J.dx = mismatch` at each iteration.
- `NRSing_CKTSO` \*: implementation of the Newton Raphson algorithm ignoring `slack_weight` (see warning above), where the
  fast `CKTSO` implementation is used to solve the linear system `J.dx = mismatch` at each iteration.
- `NRSing_SparseLU`: implementation of the Newton Raphson algorithm ignoring `slack_weight` (see warning above), where the
  Eigen default implementation is used to solve the linear system `J.dx = mismatch` at each iteration (instead of the faster `KLU`, `NICSLU` or `CKTSO`)

You can use them as:

.. code-block:: python

  from lightsim2grid.algorithm import NR_KLU  # or any of the names above

  # retrieve some Ybus, V0, etc. as above
  solver = NR_KLU()
  converged = solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)
  # process the results as above

.. note::
  \* the 6 solvers marked above (the `KLU`, `NICSLU` and `CKTSO` variants of both `NR_*` and `NRSing_*`, the 2
  `SparseLU` ones are not marked) might not be available on all platforms (`KLU` is available if installed from
  pypi, but not necessarily when installed from source). The solvers based on `NICSLU` or `CKTSO` also require an
  installation from source.

AC solvers using Gauss Seidel method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are 2 solvers in this categorie. Neither of them supports distributed slack bus [they both ignore `slack_weight` and
assign all elements of `ref` into `pv` except the first one]. If a grid with more
than 1 slack bus is provided, only the first one will be used as a slack bus, the others will be considered as "PV" nodes
(verified empirically: unlike the `NRSing_*` solvers above, the "extra" buses in `ref` do get a genuine, non-zero solved
voltage angle here, confirming they are indeed treated as PV and not kept fixed).

.. warning::

   In any case (whichever solver you use, `NRSing_*`, Gauss-Seidel or Fast-Decoupled), if your grid has more than
   one generator tagged as slack, we advise removing all but one of them from `ref` beforehand (typically by
   picking, among the generators tagged as slack, the one connected to a high voltage level, well connected, and
   as "central" as possible in the powerflow graph) rather than relying on whichever of the two behaviours above
   the algorithm you picked happens to implement.

These solvers use the Gauss Seidel method to compute powerflows. This method will iteratively update the component
of a bus based on the mismatch of the KCL. The "Gauss Seidel Synch" method is a custom implementation of this method
that updates every components at once intead of updating them one by one for each iterations.

The two solvers there are `GaussSeidelAlgo` and `GaussSeidelSynchAlgo`. Unless for some particular use case, we
do not recommend to use them as they often are slower than the Newton Raphson based solvers above.

DC solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is another type of solvers available in lightsim2grid, they use a DC modeling of the powergrid equation and
are often really fast compared to full AC powerflow.

The DC equations comes from the linearization of the AC equation, and solving a DC powerflow is basically equivalent to 
inverting a sparse matrix (or said differently solving for an equation of the sort `Ybus * Theta = Sbus` - strictly speaking 
it's not exactly this equation as we need a slack bus, for various reasons out of the scope of this documentation). 
In the current implementation it does not uses `slack_weight` and does not model distributed slack.

There are 4 solvers of this type that are different in the way they solve `Ybus * Theta = Sbus`:

- `DC_SparseLU` uses the default Eigen sparse LU implementation
- `DC_KLU` uses the fast `KLU` solver
- `DC_NICSLU` uses the fast `NICSLU` solver
- `DC_CKTSO` uses the fast `CKTSO` solver

.. code-block:: python

  from lightsim2grid.algorithm import DC_SparseLU  # or any of the names above

  # retrieve some Ybus, V0, etc. as above
  dc_solver = DC_SparseLU()
  converged = dc_solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)
  # process the results as above

Fast Decoupled solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are 8 solvers of this type. The "Fast Decoupled Power Flow" (FDPF) method exploits the (approximate)
decoupling between active power / voltage angle on one side, and reactive power / voltage magnitude on the
other, to avoid rebuilding the full Newton-Raphson jacobian at each iteration. Two variants are available,
`XB` and `BX`, that differ in which simplifying assumption is made on the B' / B'' matrices; neither supports
distributed slack bus [ignores `slack_weight`, assign all elements of `ref` into `pv` except the first one].

The list is:

- `FDPF_XB_SparseLU` / `FDPF_BX_SparseLU`: XB / BX variant, using the default Eigen sparse LU implementation.
- `FDPF_XB_KLU` \* / `FDPF_BX_KLU` \*: XB / BX variant, using the fast `KLU` implementation.
- `FDPF_XB_NICSLU` \* / `FDPF_BX_NICSLU` \*: XB / BX variant, using the fast `NICSLU` implementation.
- `FDPF_XB_CKTSO` \* / `FDPF_BX_CKTSO` \*: XB / BX variant, using the fast `CKTSO` implementation.

.. code-block:: python

  from lightsim2grid.algorithm import FDPF_XB_SparseLU  # or any of the names above

  # retrieve some Ybus, V0, etc. as above
  solver = FDPF_XB_SparseLU()
  converged = solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)
  # process the results as above

.. note::
  \* the 6 solvers marked above (the `KLU`, `NICSLU` and `CKTSO` variants of both `XB` and `BX`, out of the 8 total
  -- the 2 `SparseLU` ones are not marked) might not be available on all platforms (`KLU` is available if
  installed from pypi, but not necessarily when installed from source). The solvers based on `NICSLU` or `CKTSO`
  also require an installation from source.


Detailed documentation
--------------------------

.. automodule:: lightsim2grid.newtonpf
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`