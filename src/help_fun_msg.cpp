// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// containst he help message of some common functions (not to write them dozens of time)

#include "help_fun_msg.h"

const std::string DocSolver::get_J_python = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.
)mydelimiter";

const std::string DocSolver::get_Va = R"mydelimiter(
    Returns the voltage angles for each buses as a numpy vector of real number.
)mydelimiter";

const std::string DocSolver::get_Vm = R"mydelimiter(
    Returns the voltage magnitude for each buses as a numpy vector of real number.
)mydelimiter";

const std::string DocSolver::get_V = R"mydelimiter(
    Returns the complex voltage for each buses as a numpy vector of complex number.
)mydelimiter";

const std::string DocSolver::get_error = R"mydelimiter(
    Returns the error code (as an integer) encountered by the solver (0 = no error). TODO DOC explain better
)mydelimiter";

const std::string DocSolver::get_nb_iter = R"mydelimiter(
    Returns the number of iterations effectively performed by the solver (> 0 integer).
)mydelimiter";
const std::string DocSolver::reset = R"mydelimiter(
    Reset the solver. In this context this will clear all data used by the solver. It is mandatory to do it each time the `Ybus` matrix 
    (or any of the `pv`, or `pq` or `ref` indices vector are changed).
)mydelimiter";
const std::string DocSolver::converged = R"mydelimiter(
    Returns whether or not the solver has converged or not.
)mydelimiter";

const std::string DocSolver::compute_pf = R"mydelimiter(
    Function used to perform a powerflow.

    see section :ref:`available-powerflow-solvers` for more information about these. 

    .. warning::
        There are strong assumption made on the validity of the parameters provided. Please have a look at the section :ref:`available-powerflow-solvers`
        for more details about this.

        If a non compliant state is provided, it might result in crash of the python virtual machine (session terminates with error like 
        `segfault`) and there is absolutely no way to do anything and any data not saved on the hard drive will be lost.

    Parameters
    ------------
    Ybus: ``scipy.sparse`` matrix, CSC format
        The admittance matrix of the system
    V: ``numpy.ndarray``, vector of complex numbers
        The initial guess (and final result) for the complex angle at each bus (it is modified during the computation :)
    Sbus: ``numpy.ndarray``, vector of complex numbers
        Complex power injected at each bus
    slack_ids: ``numpy.ndarray``, vector of integers
        Gives all the ids of the buses participating to the distributed slack bus. [might be ignore by some solvers]
    slack_weights: ``numpy.ndarray``, vector of real numbers
        For each bus taking part in the distributed slack, it gives its coefficient
    pv: ``numpy.ndarray``, vector of integers
        Index of the pv buses
    pq: ``numpy.ndarray``, vector of integers
        Index of the pq buses
    max_iter: ``int``
        Maximum number of iterations performed by the solver. [might be ignore by some solvers]
    tol: ``float``
        Solver tolerance (*eg* 1e-8) [might be ignore by some solvers]

    Examples
    ---------
    Some detailed examples are provided in section :ref:`available-powerflow-solvers` of the documentation.

)mydelimiter";

const std::string DocSolver::get_timers = R"mydelimiter(
    Returns information about the time taken by some part of the solvers (in seconds)

    Times are measured in seconds using the c++ `steady_clock <https://www.cplusplus.com/reference/chrono/steady_clock/>`_ clock.

    Returns
    ---------
    timer_Fx_: ``float``
        Time spent to compute the mismatch at the KCL for each bus (both for active and reactive power)

    timer_solve_: ``float``
        Total time spent in the underlying linear solver

    timer_check_: ``float``
        Time spent in checking whether or not the mismatch of the KCL met the specified tolerance

    timer_total_: ``float``
        Total time spent in the solver

)mydelimiter";
    
const std::string DocSolver::SparseLUSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the default Eigen sparse solver available in Eigen
    for the linear algebra. 

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `SparseLU` member (*eg* `env_lightsim.backend.set_solver_typelightsim2grid.solver.SolverType.SparseLU)`).

)mydelimiter";

const std::string DocSolver::SparseLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, using the default Eigen sparse solver available in Eigen
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.SparseLUSolver` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `SparseLUSolverSingleSlack` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.SparseLUSolverSingleSlack)`).

)mydelimiter";

const std::string DocSolver::DCSolver =  R"mydelimiter(
    Default implementation of the DC solver, it uses the default Eigen sparse lu decomposition to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `DC` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.DC)`).

)mydelimiter";

const std::string DocSolver::KLUSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster KLU solver available in the SuiteSparse library
    for the linear algebra (can be unavailable if you build lightsim2grid from source). It is usually faster than the :class:`lightsim2grid.solver.SparseLUSolver`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `KLU` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.KLU)`).

)mydelimiter";

const std::string DocSolver::KLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm,the faster KLU solver available in the SuiteSparse library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.KLUSolver`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `KLUSingleSlack` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.KLUSingleSlack)`).

)mydelimiter";

const std::string DocSolver::KLUDCSolver = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster KLU solver available in the SuiteSparse library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (can be unavailable if you build lightsim2grid from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `KLUDC` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.KLUDC)`).

)mydelimiter";

const std::string DocSolver::NICSLUSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It is usually faster than the :class:`lightsim2grid.solver.SparseLUSolver`. (requires a build from source)
    
    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `NICSLU` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.NICSLU)`).

)mydelimiter";

const std::string DocSolver::NICSLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.NICSLUSolver` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `NICSLUSingleSlack` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.NICSLUSingleSlack)`).

)mydelimiter";

const std::string DocSolver::NICSLUDCSolver = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster NICSLU solver available in the NICSLU library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (requires a build from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `NICSLUSingleSlack` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.NICSLUSingleSlack)`).

)mydelimiter";

const std::string DocSolver::GaussSeidelSolver = R"mydelimiter(
    Default implementation of the "Gauss Seidel" powerflow solver. We do not recommend to use it as the Newton Raphson based solvers
    are usually much faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `GaussSeidel` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.GaussSeidel)`).

)mydelimiter";

const std::string DocSolver::GaussSeidelSynchSolver = R"mydelimiter(
    Variant implementation of the "Gauss Seidel" powerflow solver, where every buses are updated at once (can be significantly faster than the 
    :class:`lightsim2grid.solver.GaussSeidelSolver` for larger grid). We still do not recommend to use it as the Newton Raphson based solvers
    are usually much faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is referred to by the `GaussSeidelSynch` member (*eg* `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SolverType.GaussSeidelSynch)`).

)mydelimiter";

const std::string DocSolver::AnySolver = R"mydelimiter(
    This is a "wrapper" class that allows the user to perform some powerflow using the same API using different solvers. It is not recommended
    to use this wrapper directly. It is rather a class exported to be compatible with the `env_lightsim2grid.backend._grid.get_solver()` method.

    Examples
    ---------
    This class is built to be used like this:

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        anysolver = env.backend._grid.get_solver()

        anysolver.get_type()  # type of solver currently used
        anysolver.get_J()  # current Jacobian matrix, if available by the method

)mydelimiter";

const std::string DocSolver::get_type = R"mydelimiter(
    Retrieve the current solver used. This will return an instance of :class:`lightsim2grid.solver.SolverType` indicating which
    is the underlying solver in use.

    This should be equivalent to :func:`lightsim2grid.initGridModel.GridModel.get_solver_type()`

)mydelimiter";
const std::string DocSolver::chooseSolver_get_J_python = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.

    .. note::
        Depending on the underlying solver used (*eg* :class:`lightsim2grid.solver.DCSolver` or :class:`lightsim2grid.solver.GaussSeidelSolver`)
        the jacobian matrix might be irrelevant and an attempt to use this function will throw a RuntimeError. 

)mydelimiter";
const std::string DocSolver::get_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver when performing a powerflow.

    This is equivalent to the `timer_total_` returned value of the`***.get_timers()` function.
)mydelimiter";

const std::string DocIterator::id = R"mydelimiter(
    Get the ideas of the element. Ids are integer from 0 to n-1 (if `n` denotes the number of such elements on the grid.)

    Examples
    --------
    We give the example only for generators, but it works similarly for every other types of objects
    in a :class:`lightsim2grid.initGridModel.GridModel`.
    
    This gives something like:

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        first_gen = grid_model.get_generators()[0]  # or get_loads for loads, etc.
        first_gen.id  # should be 0

)mydelimiter";

const std::string DocIterator::connected = R"mydelimiter(
    Get the status (True = connected, False = disconnected) of each element of a :class:`lightsim2grid.initGridModel.GridModel`

)mydelimiter";

const std::string DocIterator::bus_id = R"mydelimiter(
    Get the bus id (as an integer) at which each element of a :class:`lightsim2grid.initGridModel.GridModel` is connected. If `-1` is returned it means
    that the object is disconnected.

)mydelimiter";

const std::string DocIterator::target_p_mw = R"mydelimiter(
    Get the active production (or consumption) setpoint in MW for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

)mydelimiter";

const std::string DocIterator::target_q_mvar = R"mydelimiter(
    Get the reactive production (or consumption) setpoint in MVAr for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

)mydelimiter";

const std::string DocIterator::line_model = R"mydelimiter(
    The "line model" (also valid for transformers) is: 
    
    .. code-block :: none

                    ior                      ________            iex
        `or bus` o------>   -----------------|r + j.x|---------<-------o `ex bus`
                 |       ) (            |                  |           |
                 |       ) (         |     |            |     |        |
                 | vor   ) ( n:1     |1/2*h|            |1/2*h|        | vex
                 |       ) (         |     |            |     |        |
                 \/      ) (            |                  |           \/
        ground---o-------   -------------------------------------------o---- ground

    (fyi: `ior`, `iex`, `n` and `y` are all complex numbers. `r` and `x` are real numbers. `j` is a complex number such that `j^2 = -1`)

)mydelimiter";

const std::string DocIterator::r_pu = R"mydelimiter(
    Retrieve the resistance (given in pair unit system, and not in Ohm) of the powerlines or the transformers. This is a real number
    and is represented by the number `r` in the line model.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::x_pu = R"mydelimiter(
    Retrieve the reactance (given in pair unit system, and not in Ohm) of the powerlines or the transformers. This is a real number
    and is represented by the number `x` in the line model.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::h_pu = R"mydelimiter(
    Retrieve the capacitance (real part) and dielectric conductance (imaginary part)
    (given in pair unit system) of the powerlines or the transformers. 
    
    This is a complex number and is represented by the number `h` in the line model.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::res_p_mw = R"mydelimiter(
    Get the active production (or consumption) in MW for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_q_mvar = R"mydelimiter(
    Get the reactive production (or consumption) in MVAr for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_theta_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this object is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter";

const std::string DocIterator::res_v_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this object is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter";

const std::string DocIterator::target_vm_pu = R"mydelimiter(
    Get the voltage magnitude setpoint (in pair unit and NOT in kV) for each element of the grid supporting this feature.

    .. warning::
        This is given in "pair unit" (pu) system and not in kilo Volt (kV) !

)mydelimiter";

const std::string DocIterator::has_res = R"mydelimiter(
    This property specify whether or not a given element contains some "result" information. If set to ``True`` then the fields
    starting with `res_` (*eg* `res_p_mw`) are filled otherwise they are initialized with an arbitrary (and meaningless) value.

)mydelimiter";

const std::string DocIterator::DataGen = R"mydelimiter(
    This class allows to iterate through the generators of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    In lightsim2grid they are modeled as "pv" meanings you give the active production setpoint and voltage magnitude setpoint
    (see :attr:`lightsim2grid.elements.DataSGen` for more exotic PQ generators).

    The active production value setpoint are modified only for the generators participating to the slack buses
    (see :attr:`lightsim2grid.elements.GenInfo.is_slack` and :attr:`lightsim2grid.elements.GenInfo.slack_weight`).

    Generators are modeled as in pandapower and can be represented a the 
    `pandapower generators <https://pandapower.readthedocs.io/en/latest/elements/gen.html#electric-model>`_ .

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        for gen in grid_model.get_generators():
            # do something with gen !
            gen.bus_id

        print(f"There are {len(grid_model.get_generators())} generators on the grid.")

        first_generator = grid_model.get_generators()[0]

    You can have a look at :class:`lightsim2grid.elements.GenInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::GenInfo = R"mydelimiter(
    This class represents what you get from retrieving some elements from :class:`lightsim2grid.elements.DataGen`

    It allows to read information from each generator of the powergrid.

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        first_generator = grid_model.get_generators()[0]  # first generator is a `GenInfo`

        for gen in grid_model.get_generators():
            # gen is a `GenInfo`
            gen.bus_id

)mydelimiter";

const std::string DocIterator::is_slack = R"mydelimiter(
    Tells whether or not this generator paticipated to the distributed slack bus.

    .. note:: 
        Depending on the solver used, it is possible that a generator we asked to participate to the distributed slack bus
        do not participate to it (for example if there is a more than one generator where `is_slack` is ``True`` but the model used
        to computed the powerflow do not support distributed slack buses - **eg** :class:`lightsim2grid.solver.SparseLUSolverSingleSlack`)

        This is why we recommend to use the (slower) but more accurate :class:`lightsim2grid.solver.SparseLUSolver` or 
        :class:`lightsim2grid.solver.KLUSolver` for example.

)mydelimiter";

const std::string DocIterator::slack_weight = R"mydelimiter(
    For each generators, gives the participation (for the distributed slack) of this particular generator.
    
    .. note::
        Weights do not scale to one for this variable thus this number has no meaning by itself and should be compared with the others.

)mydelimiter";

const std::string DocIterator::min_q_mvar = R"mydelimiter(
    Minimum reactive value that can be produced / absorbed by this generator given MVAr.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.initGridModel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::max_q_mvar = R"mydelimiter(
    Maximum reactive value that can be produced / absorbed by this generator given MVAr.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.initGridModel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::min_p_mw = R"mydelimiter(
    Minimum active value that can be produced / absorbed by this generator given in MW.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.initGridModel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::max_p_mw = R"mydelimiter(
    Maximum active value that can be produced / absorbed by this generator given in MW.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.initGridModel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::DataSGen = R"mydelimiter(
    This class allows to iterate through the static generators of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    In lightsim2grid they are two types of generators the more standard PV generators (see :attr:`lightsim2grid.elements.DataGen`). These
    are more exotic generators known as PQ, where you give the active production value and reactive production value. It's basically like loads,
    but using the generator convention (if the value is positive, it means power is taken from the grid to the element)

    They cannot participate to the distributed slack bus.

    Static generators are modeled as in pandapower and can be represented a the 
    `pandapower static generators <https://pandapower.readthedocs.io/en/latest/elements/sgen.html#electric-model>`_ .

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the static generators
        for sgen in grid_model.get_static_generators():
            # do something with sgen !
            sgen.bus_id

        print(f"There are {len(grid_model.get_static_generators())} static generators on the grid.")

        first_static_generator = grid_model.get_static_generators()[0]

    You can have a look at :class:`lightsim2grid.elements.SGenInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::SGenInfo = R"mydelimiter(
    This class represents what you get from retrieving some elements from :class:`lightsim2grid.elements.DataSGen`

    It allows to read information from each static generator of the powergrid.

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # do something with the static generators
        first_static_generator = grid_model.get_static_generators()[0]  # first static generator is a `SGenInfo`

        for sgen in grid_model.get_static_generators():
            # sgen is a `SGenInfo`
            sgen.bus_id

)mydelimiter";

const std::string DocIterator::DataLoad = R"mydelimiter(
    This class allows to iterate through the loads **and storage units** of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    They cannot participate to the distributed slack bus yet. If you want this feature, fill free to send us a github issue.

    Loads are modeled as in pandapower and can be represented a the 
    `pandapower loads <https://pandapower.readthedocs.io/en/latest/elements/load.html#electric-model>`_ .

    .. note::
        lightsim2grid Storages are modeled as load.
    
    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the load
        for load in grid_model.get_loads():
            # do something with load !
            load.bus_id

        print(f"There are {len(grid_model.get_loads())} loads on the grid.")

        first_load = grid_model.get_loads()[0]

        # or the storage units
        for storage in grid_model.get_storages():
            # do something with storage !
            storage.bus_id

        print(f"There are {len(grid_model.get_storages())} storage units on the grid.")

        first_storage_unit = grid_model.get_storages()[0]

    You can have a look at :class:`lightsim2grid.elements.LoadInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::LoadInfo = R"mydelimiter(
    This class represents what you get from retrieving some elements from :class:`lightsim2grid.elements.DataLoad`.
    We remind the reader that storage units are also modeled as load in lightsim2grid.

    It allows to read information from each load / storage unit of the powergrid.

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    .. note::
        lightsim2grid Storages are modeled as load.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # for loads
        first_load = grid_model.get_loads()[0]  # first static generator is a `LoadInfo`
        for load in grid_model.get_loads():
            # load is a `LoadInfo`
            load.bus_id


        # for loads
        first_storage_unit = grid_model.get_storages()[0]  # first static generator is a `LoadInfo`
        for storage in grid_model.get_storages():
            # storage is a `LoadInfo`
            storage.bus_id

)mydelimiter";

const std::string DocIterator::DataShunt = R"mydelimiter(
    This class allows to iterate through the load of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    Shunts are modeled as in pandapower and can be represented a the 
    `pandapower shunts <https://pandapower.readthedocs.io/en/latest/elements/shunt.html#electric-model>`_ .
    
    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the load
        for shunt in grid_model.get_shunts():
            # do something with shunt !
            shunt.bus_id

        print(f"There are {len(grid_model.get_shunts())} shunts on the grid.")

        first_shunt = grid_model.get_shunts()[0]

    You can have a look at :class:`lightsim2grid.elements.ShuntInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::ShuntInfo = R"mydelimiter(
    This class represents what you get from retrieving the shunts from :class:`lightsim2grid.elements.DataShunt`.

    It allows to read information from each shunt of the powergrid.

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # for shunts
        first_shunt = grid_model.get_shunts()[0]  # first shunt, this is a `ShuntInfo`
        for shunt in grid_model.get_shunts():
            # shunt is a `ShuntInfo`
            shunt.bus_id

)mydelimiter";

const std::string DocIterator::DataTrafo = R"mydelimiter(
    This class allows to iterate through the transformers of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    Transformers are modeled as in pandapower and can be represented a the 
    `pandapower transformers <https://pandapower.readthedocs.io/en/latest/elements/trafo.html#electric-model>`_ .
    
    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the tranformers
        for trafo in grid_model.get_trafos():
            # do something with trafo !
            trafo.bus_hv_id

        print(f"There are {len(grid_model.get_trafos())} transformers on the grid.")

        first_transformer = grid_model.get_trafos()[0]

    You can have a look at :class:`lightsim2grid.elements.TrafoInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::TrafoInfo = R"mydelimiter(
    This class represents what you get from retrieving the transformers from :class:`lightsim2grid.elements.DataTrafo`.

    It allows to read information from each transformer of the powergrid.

    Transformers have two sides, one is "hv" for "high voltage" and one is "lv" for "low voltage" that are connected and linked to each other
    by some equations.

    For accessing the results, it's basically the same as having two "elements" (so you get two "voltage_magnitude" `res_v_kv`,
    two "injected power" `res_p_mw` etc.)

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # for transformers
        first_transformer = grid_model.get_trafos()[0]  # first transformer, this is a `TrafoInfo`
        for trafo in grid_model.get_trafos():
            # trafo is a `TrafoInfo`
            trafo.bus_hv_id

    Notes
    -----
    Transformer are modeled using the "line model".

    Usually, the "or" side is the "hv" side and the "ex" side is the "lv" side.

    The tap ratio `n` bellow is a complex number with its magnitude corresponding to the tap ratio and its
    angle to the phase shifter.

    For more information about the model and the equations linking all the quantities, please visit 
    `matpower manual <https://matpower.org/docs/MATPOWER-manual.pdf>`_ , especially the "3. Modeling" and the
    "3.2 Branches" subsection, as well as the equation 3.1, 3.2 and 3.3 therein.
    
)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::ratio = R"mydelimiter(
    Retrieve the ratio (absolute value of the complex coefficient `n` in the powerline model). It has no units

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::shift_rad = R"mydelimiter(
    Retrieve the shift angle (angle of the complex coefficient `n` in the powerline model). It is given in radian (and not in degree)

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::is_tap_hv_side = R"mydelimiter(
    Gives whether the tap (both for the ratio and the phase shifter) is located "hv" side (default, when ``True``) or 
    "lv" side (when ``False``).

)mydelimiter";

const std::string DocIterator::bus_hv_id = R"mydelimiter(
    Get the bus id (as an integer) at which the "hv" side of the transformer is connected. If `-1` is returned it means
    that the transformer is disconnected.

)mydelimiter";

const std::string DocIterator::bus_lv_id = R"mydelimiter(
    Get the bus id (as an integer) at which the "lv" side of the transformer is connected. If `-1` is returned it means
    that the transformer is disconnected.

)mydelimiter";

const std::string DocIterator::res_p_hv_mw = R"mydelimiter(
    Get the active power in MW for at the "hv" side of the transformer. If it is positive it means power is absorbed by the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_p_lv_mw = R"mydelimiter(
    Get the active power in MW for at the "lv" side of the transformer. If it is positive it means power is absorbed by the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_q_hv_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "hv" side of the transformer. If it is positive it means power is absorbed by the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_q_lv_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "lv" side of the transformer. If it is positive it means power is absorbed by the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_theta_hv_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "hv" side of the transformer is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter";

const std::string DocIterator::res_theta_lv_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "lv" side of the transformer is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter";

const std::string DocIterator::res_v_hv_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "hv" side of the transformer is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter";

const std::string DocIterator::res_v_lv_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "lv" side of the transformer is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter";

const std::string DocIterator::res_a_lv_ka = R"mydelimiter(
    Get the current flows (in kA) at the "lv" side of the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_a_hv_ka = R"mydelimiter(
    Get the current flows (in kA) at the "hv" side of the transformer.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::DataLine = R"mydelimiter(
    This class allows to iterate through the powerlines of the :class:`lightsim2grid.initGridModel.GridModel` easily, as if they were
    in a python list.

    Powerlines are modeled as in pandapower and can be represented a the 
    `pandapower lines <https://pandapower.readthedocs.io/en/latest/elements/line.html#electric-model>`_ .
    
    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the powerlines
        for line in grid_model.get_lines():
            # do something with line !
            line.bus_or_id

        print(f"There are {len(grid_model.get_lines())} lines on the grid.")

        first_line = grid_model.get_lines()[0]

    You can have a look at :class:`lightsim2grid.elements.LineInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::LineInfo = R"mydelimiter(
    This class represents what you get from retrieving the powerlines from :class:`lightsim2grid.elements.DataLine`.

    It allows to read information from each powerline of the powergrid.

    Powerlines have two sides, one is "or" for "origin" and one is "ex" for "extremity" that are connected and linked to each other
    by some equations.

    For accessing the results, it's basically the same as having two "elements" (so you get two "voltage_magnitude" `res_v_kv`,
    two "injected power" `res_p_mw` etc.)

    .. warning::
        Data ca only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # for powerlines
        first_line = grid_model.get_lines()[0]  # first line, this is a `LineInfo`
        for line in grid_model.get_lines():
            # line is a `LineInfo`
            line.bus_or_id

    Notes
    -----
    Line are modeled using the "line model" as shown in the schema at the end of the paragraph.

    The tap ratio `n` on this schema will be 1.0 for all powerline. If you want to model phase shifters, please
    model them as Trafo (see :class:`lightsim2grid.elements.TrafoInfo`)

    For more information about the model and the equations linking all the quantities, please visit 
    `matpower manual <https://matpower.org/docs/MATPOWER-manual.pdf>`_ , especially the "3. Modeling" and the
    "3.2 Branches" subsection, as well as the equation 3.1, 3.2 and 3.3 therein.
    
)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::bus_or_id = R"mydelimiter(
    Get the bus id (as an integer) at which the "or" side of the line is connected. If `-1` is returned it means
    that the line is disconnected.

)mydelimiter";

const std::string DocIterator::bus_ex_id = R"mydelimiter(
    Get the bus id (as an integer) at which the "lv" side of the line is connected. If `-1` is returned it means
    that the line is disconnected.

)mydelimiter";

const std::string DocIterator::res_p_or_mw = R"mydelimiter(
    Get the active power in MW for at the "or" side of the line. If it is positive it means power is absorbed by the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_p_ex_mw = R"mydelimiter(
    Get the active power in MW for at the "ex" side of the line. If it is positive it means power is absorbed by the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_q_or_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "or" side of the line. If it is positive it means power is absorbed by the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_q_ex_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "ex" side of the line. If it is positive it means power is absorbed by the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_theta_or_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "or" side of the line is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter";

const std::string DocIterator::res_theta_ex_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "ex" side of the line is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter";

const std::string DocIterator::res_v_or_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "or" side of the line is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter";

const std::string DocIterator::res_v_ex_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "ex" side of the line is connected.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter";

const std::string DocIterator::res_a_or_ka = R"mydelimiter(
    Get the current flows (in kA) at the "or" side of the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_a_ex_ka = R"mydelimiter(
    Get the current flows (in kA) at the "ex" side of the line.

    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";
