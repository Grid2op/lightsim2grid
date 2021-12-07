Use as Pandapower Solver
=========================

.. versionchanged:: 0.5.6

    As of version 0.5.6 lightsim2grid implements the new API of "newtonpf" required by pandapower. This means that
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

    from lighsim2grid.newtonpf import newtonpf
    V, converged, iterations, J = newtonpf(Ybus, V, Sbus, ref, pv, pq, ppci, options)

This function uses the KLU algorithm (or the solver provided in Eigen if KLU has not been instealld) 
and a c++ implementation of a Newton solver for speed.

.. note::

  The oldest `newtonpf` function compatible with older version of pandapower (*eg* <=2.6.0) can still be accessed with
  `from lightsim2grid.newtonpf import newtonpf_old`

.. _available-powerflow-solvers: 

Even more advanced usage
########################
You can customize even more the solvers that you want to use.

Lightsim2grid comes with 11 available solvers that can solver either AC or DC powerflows. We can cluster them into 3 categories.

If you want to stay "relatively high level", once you have a grid model, you can change the solver using
the "enum" of the solvers you want to use as showed bellow:

.. code-block:: python

    from lightsim2grid.solver import SolverType
    # init the grid model
    from lightsim2grid.initGridModel import init
    pp_net = ...  # any pandapower grid
    lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

    # change the solver used for the powerflow
    lightsim_grid_model.change_solver(SolverType.KLUSolver)  # change the NR solver that uses KLU
    # you can replace `SolverType.KLUSolver` by any of the 11 available solvers described bellow, 
    # for example (and we will not write the 11...) `SolverType.KLUSolverSingleSlack`, `SolverType.SparseLUSolver` 
    # or even `SolverType.NICSLUSolver`
        
All solvers can be accessed with the same API (if you want to use the raw python class, not recommended):

.. code-block:: python

  from lightsim2grid.solver import ASolverAvailable
  Ybus = ...  # a csc sparse matrix (it's really important that it is a csc and not a csr !)
  V0 = ...  # a complex vector (initial guess)
  Sbus = ...  # a complex vector (power injected at each bus)
  ref = ...  # a vector of integer giving the bus id of the slack buses
  slack_weight = ...  # a vector of real number giving the weights associated to each slack bus
  pv = ...  # a vector of integers giving the bus id of the PV bus
  pq = ...  # a vector of integers giving the bus id of the PQ bus
  max_it = ...  # a > 0 integer maximum number of iterations the solver is allowed to perform
  tol = ...  # a > 0. real number giving the maximum KCL violation allowed for a all nodes

  solver = ASolverAvailable()
  converged = solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)

  # to retrieve the voltages related information (in case converged is True)
  solver.get_Va()  # voltage magnitude
  solver.get_Vm()  # voltage angle
  solver.get_V()  # complex voltage
  # for compatible solvers
  solver.get_J()  # see documentation of the `newton_pf` function for more information about the shape of J.

  # some other usefull information
  solver.get_nb_iter()  # return the number of iteration performed
  solver.get_timers()  # some execution times for some function (TODO DOC)
  sovler.get_error()  # the id of the error encountered  (TODO DOC)
  sovler.converged()  # equal to the boolean `converged` above

Be carefull, there are some constraints on the data that are not necessarily checked, and might lead to hugly crash of the
python virtual machine at execution time. So we encourage you to check that:

- tol > 0.
- maxt_it > 0
- Ybus is a squared sparse matrix, in CSC format (see documentation of scipy sparse for more information) **It is 
  really important that this matrix is in CSC format**
- `Sbus` and `V0` have the same size which corresponds to the size (number of rows or columns) of `Ybus`
- for all node id in `ref`, `slack_weight[node id] > 0.`
- `sum(slack_weight) = 1.` and all elements of `slack_weight` are > 0.
- all the buses are on `ref` (for slack buses) or on `pv` (for PV buses) or on `pq` (for PQ buses)
  [informatically, this means that the ensemble `[0, len(V0) - 1]` is included in the union `ref U pv U pq` ]
- all buses are only in one of `ref`, `pv` and `pq` [informatically an element of `[0, len(V0) - 1]` cannot be at the 
  same time in `pv` and `pq` or in `ref` and `pv` or in `ref` and `pq`
- there should be at least one element in `ref` (`len(ref) > 0`)
  
.. warning::

    Just to emphasize that if any of the condition above is not met, this can result in crash of the python
    virtual machine without any exception thrown (segfault). 

    This is why we do not recommend to use these solvers directly !


AC solvers using Newton Raphson
+++++++++++++++++++++++++++++++

There are 6 solvers in this categorie. They can in turn, be split into two main sub categories. The first one allows for a
distributed slack bus (but can be a bit slower) as the other one does not allow for such (in case of multiple slack bus, only 
the first one is used as a real slack bus, the other ones are converted silently to PV buses)

The list is:

- `KLUSolver` \*: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the 
  fast `KLU` implementation is used to iteratively update the jacobian matrix `J`.
- `NICSLUSolver` \*: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the 
  fast `NICSLU` implementation is used to iteratively update the jacobian matrix `J`.
- `SparseLUSolver`: implementation of the Newton Raphson algorithm supporting the distributed slack bus, where the 
  Eigen default implementation is used to iteratively update the jacobian matrix `J` (instead of the faster `KLU` or `NICSLU`)
- `KLUSolverSingleSlack` \*: implementation of the Newton Raphson algorithm only supporting single slack bus [ignores `slack_weight`, assign 
  all elements of `ref` into `pv` except the first one], where the 
  fast `KLU` implementation is used to iteratively update the jacobian matrix `J`
- `NICSLUSolverSingleSlack` \*: implementation of the Newton Raphson algorithm only supporting single slack bus [ignores `slack_weight`, assign 
  all elements of `ref` into `pv` except the first one], where the 
  fast `NICSLU` implementation is used to iteratively update the jacobian matrix `J`.
- `SparseLUSolverSingleSlack`: implementation of the Newton Raphson algorithm only supporting single slack bus [ignores `slack_weight`, assign 
  all elements of `ref` into `pv` except the first one], where the 
  Eigen default implementation is used to iteratively update the jacobian matrix `J` (instead of the faster `KLU` or `NICSLU`)

You can use them as:

.. code-block:: python

  from lightsim2grid.solver import KLUSolver  # or any of the names above

  # retrieve some Ybus, V0, etc. as above
  solver = KLUSolver()
  converged = solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)
  # process the results as above

.. note::
  \* these 4 solvers might not be available on all platforms (KLU is available if installed from pypi, but not
  necessarily when installed from source). The solvers based on `NICSLU` also requires an installation from
  source.

AC solvers using Gauss Seidel method
+++++++++++++++++++++++++++++++++++++

There are 2 solvers in this categorie. Neither of them supports distributed slack bus [they both ignore `slack_weight` and
assign all elements of `ref` into `pv` except the first one]. If a grid with more
more than 1 slack bus is provided, only the first one will be used as a slack bus, the others will be considered as "PV" nodes.

These solvers use the Gauss Seidel method to compute powerflows. This method will iteratively update the component
of a bus based on the mismatch of the KCL. The "Gauss Seidel Synch" method is a custom implementation of this method
that updates every components at once intead of updating them one by one for each iterations.

The two solvers there are `GaussSeidelSolver` and `GaussSeidelSynchSolver`. Unless for some particular use case, we
do not recommend to use them as they often are slower than the Newton Raphson based solvers above.

DC solvers
+++++++++++
This is another type of solvers available in lightsim2grid, they use a DC modeling of the powergrid equation and
are often really fast compared to full AC powerflow.

The DC equations comes from the linearization of the AC equation, and solving a DC powerflow is basically equivalent to 
inverting a sparse matrix (or said differently solving for an equation of the sort `Ybus * Theta = Sbus` - strictly speaking 
it's not exactly this equation as we need a slack bus, for various reasons out of the scope of this documentation). 
In the current implementation it does not uses `slack_weight` and does not model distributed slack.

There are 3 solvers of this type that are different in the way they solve `Ybus * Theta = Sbus`:

- `DCSolver` uses the default Eigen sparse LU implementation
- `KLUDCSolver` uses the fast `KLU` solver
- `NICSLUDCSolver` uses the fast `NICSLU` solver    

.. code-block:: python

  from lightsim2grid.solver import DCSolver  # or any of the names above

  # retrieve some Ybus, V0, etc. as above
  dc_solver = DCSolver()
  converged = dc_solver.solve(Ybus, V0, Sbus, ref, slack_weights, pv, pq, max_it, tol)
  # process the results as above

Detailed documentation
######################

.. automodule:: lightsim2grid.newtonpf
    :members:
    :autosummary:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`