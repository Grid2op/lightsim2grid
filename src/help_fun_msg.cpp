// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
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

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `SparseLU`.
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.SparseLU)` at creation time    
    
    .. note::
        Available on all plateform, this is the default solver used when :class:`lightsim2grid.solver.KLUSolverSingleSlack`
        is not found (when a "single slack" is detected).

)mydelimiter";

const std::string DocSolver::SparseLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, using the default Eigen sparse solver available in Eigen
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.SparseLUSolver` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `SparseLUSingleSlack` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.SparseLUSingleSlack)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.SparseLUSingleSlack)` at creation time

    .. note::
        Available on all plateform, this is the default solver used when a distributed slack bus is detected and :class:`lightsim2grid.solver.SolverType.KLUSolver`
        is not found.

)mydelimiter";

const std::string DocSolver::DCSolver =  R"mydelimiter(
    Default implementation of the DC solver, it uses the default Eigen sparse lu decomposition to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `DC` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.DC)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.DC)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

)mydelimiter";

const std::string DocSolver::FDPF_XB_SparseLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the default Eigen sparse lu decomposition for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_SparseLU)` at creation time

)mydelimiter";

const std::string DocSolver::FDPF_BX_SparseLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the default Eigen sparse lu decomposition for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_SparseLU)` at creation time

)mydelimiter";

const std::string DocSolver::KLUSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster KLU solver available in the SuiteSparse library
    for the linear algebra (can be unavailable if you build lightsim2grid from source). It is usually faster than the :class:`lightsim2grid.solver.SparseLUSolver`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.KLU)` at creation time

    .. note::
        This is the default solver used when a distributed slack bus is detected (when it's available, otherwise see :class:`lightsim2grid.solver.SparseLUSolver`).

)mydelimiter";

const std::string DocSolver::KLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm,the faster KLU solver available in the SuiteSparse library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.KLUSolver`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `KLUSingleSlack` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.KLUSingleSlack)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.KLUSingleSlack)` at creation time

    .. note::
        This is the default solver used when available.

)mydelimiter";

const std::string DocSolver::KLUDCSolver = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster KLU solver available in the SuiteSparse library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (can be unavailable if you build lightsim2grid from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `KLUDC` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.KLUDC)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.KLUDC)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

)mydelimiter";

const std::string DocSolver::FDPF_XB_KLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdbx"  in pypower / pandapower), it uses the fast KLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_KLU)` at creation time

)mydelimiter";

const std::string DocSolver::FDPF_BX_KLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 3" / "fdxb"  in pypower / pandapower), it uses the fast KLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_KLU)` at creation time

)mydelimiter";

const std::string DocSolver::NICSLUSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It is usually faster than the :class:`lightsim2grid.solver.SparseLUSolver`. (requires a build from source)
    
    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.NICSLU)` at creation time

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu
        
)mydelimiter";

const std::string DocSolver::NICSLUSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.NICSLUSolver` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `NICSLUSingleSlack` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.NICSLUSingleSlack)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.NICSLUSingleSlack)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::NICSLUDCSolver = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster NICSLU solver available in the NICSLU library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (requires a build from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `NICSLUDC` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.NICSLUDC)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.NICSLUDC)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu
 
)mydelimiter";

const std::string DocSolver::FDPF_XB_NICSLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the fast NICSLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_NICSLU)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::FDPF_BX_NICSLUSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the fast NICSLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_NICSLU)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::CKTSOSolver = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster CKTSO solver available in the CKTSO library
    for the linear algebra (requires a build from source)
    
    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.CKTSO)` at creation time

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso
 
)mydelimiter";

const std::string DocSolver::CKTSOSolverSingleSlack = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, the faster CKTSO solver available in the CKTSO library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.solver.CKTSOSolver` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `CKTSOSingleSlack` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.CKTSOSingleSlack)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.CKTSOSingleSlack)` at creation time

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso
 
)mydelimiter";

const std::string DocSolver::CKTSODCSolver = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster CKTSO solver available in the CKTSO library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (requires a build from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `CKTSODC` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.CKTSODC)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.CKTSODC)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::FDPF_XB_CKTSOSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the fast CKTSO library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_CKTSO)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for cktso.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::FDPF_BX_CKTSOSolver =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the fast CKTSO library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `FDPF_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.FDPF_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.FDPF_CKTSO)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for cktso.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::GaussSeidelSolver = R"mydelimiter(
    Default implementation of the "Gauss Seidel" powerflow solver. We do not recommend to use it as the Newton Raphson based solvers
    are usually much (much) faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it is called `GaussSeidel` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.GaussSeidel)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.GaussSeidel)` at creation time

    .. warning::
        It currently does not support distributed slack.

)mydelimiter";

const std::string DocSolver::GaussSeidelSynchSolver = R"mydelimiter(
    Variant implementation of the "Gauss Seidel" powerflow solver, where every buses are updated at once (can be significantly faster than the 
    :class:`lightsim2grid.solver.GaussSeidelSolver` for larger grid). We still do not recommend to use it as the Newton Raphson based solvers
    are usually much (much) faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.solver.SolverType`, it called `GaussSeidelSynch` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_solver_type(lightsim2grid.solver.GaussSeidelSynch)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.solver.GaussSeidelSynch)` at creation time

    .. warning::
        It currently does not support distributed slack.
        
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

    This should be equivalent to :func:`lightsim2grid.gridmodel.GridModel.get_solver_type()`

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
    Get the id of the element. Ids are integer from 0 to n-1 (if `n` denotes the number of such elements on the grid.)

    Examples
    --------
    We give the example only for generators, but it works similarly for every other types of objects
    in a :class:`lightsim2grid.gridmodel.GridModel`.
    
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

const std::string DocIterator::name = R"mydelimiter(
    Get the name of the element. Names are string that should be unique. But if you really want things unique, use the `id`

    .. warning::
        Names are optional and might not be set when reading the grid. 

    Examples
    --------
    We give the example only for generators, but it works similarly for every other types of objects
    in a :class:`lightsim2grid.gridmodel.GridModel`.
    
    This gives something like:

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        first_gen = grid_model.get_generators()[0]  # or get_loads for loads, etc.
        first_gen.name 

)mydelimiter";

const std::string DocIterator::connected = R"mydelimiter(
    Get the status (True = connected, False = disconnected) of each element of a :class:`lightsim2grid.gridmodel.GridModel`

)mydelimiter";

const std::string DocIterator::bus_id = R"mydelimiter(
    Get the bus id (as an integer) at which each element of a :class:`lightsim2grid.gridmodel.GridModel` is connected. If `-1` is returned it means
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


const std::string DocIterator::only_avail_res = R"mydelimiter(
    
    .. warning::
        This feature is only relevant if the results have been computed (for example if a powerflow has successfully run)

)mydelimiter";

const std::string DocIterator::res_p_mw = R"mydelimiter(
    Get the active production (or consumption) in MW for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_mvar = R"mydelimiter(
    Get the reactive production (or consumption) in MVAr for element of the grid supporting this feature.

    For generators (and static generators) it is given following the "generator convention" (positive = power is injected to the grid)
    
    For loads (and storage units) it is given following the "load convention" (positive = power is absorbed from the grid)

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this object is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this object is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::target_vm_pu = R"mydelimiter(
    Get the voltage magnitude setpoint (in pair unit and NOT in kV) for each element of the grid supporting this feature.

    .. warning::
        This is given in "pair unit" (pu) system and not in kilo Volt (kV) !

)mydelimiter";

const std::string DocIterator::has_res = R"mydelimiter(
    This property specify whether or not a given element contains some "result" information. If set to ``True`` then the fields
    starting with `res_` (*eg* `res_p_mw`) are filled otherwise they are initialized with an arbitrary (and meaningless) value.

)mydelimiter";

const std::string DocIterator::GeneratorContainer = R"mydelimiter(
    This class allows to iterate through the generators of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
    in a python list.

    In lightsim2grid they are modeled as "pv" meanings you give the active production setpoint and voltage magnitude setpoint
    (see :attr:`lightsim2grid.elements.SGenContainer` for more exotic PQ generators).

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
    This class represents what you get from retrieving some elements from 
    :class:`lightsim2grid.elements.GeneratorContainer`

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
        to computed the powerflow do not support distributed slack buses - **eg** :class:`lightsim2grid.solver.SparseLUSingleSlack`)

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
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.gridmodel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::max_q_mvar = R"mydelimiter(
    Maximum reactive value that can be produced / absorbed by this generator given MVAr.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.gridmodel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::min_p_mw = R"mydelimiter(
    Minimum active value that can be produced / absorbed by this generator given in MW.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.gridmodel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::max_p_mw = R"mydelimiter(
    Maximum active value that can be produced / absorbed by this generator given in MW.

    .. note:: 
        This is for now not taken into account by the solver. It is only used in :func:`lightsim2grid.gridmodel.check_solution` if `check_q_limits` is
        set to ``True``

)mydelimiter";

const std::string DocIterator::SGenContainer = R"mydelimiter(
    This class allows to iterate through the static generators of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
    in a python list.

    In lightsim2grid they are two types of generators the more standard PV generators (see 
    :attr:`lightsim2grid.elements.GeneratorContainer`). These
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
    This class represents what you get from retrieving some elements from 
    :class:`lightsim2grid.elements.SGenContainer`

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

const std::string DocIterator::LoadContainer = R"mydelimiter(
    This class allows to iterate through the loads **and storage units** of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
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
    This class represents what you get from retrieving some elements from 
    :class:`lightsim2grid.elements.LoadContainer`.
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

const std::string DocIterator::ShuntContainer = R"mydelimiter(
    This class allows to iterate through the load of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
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
    This class represents what you get from retrieving the shunts from 
    :class:`lightsim2grid.elements.ShuntContainer`.

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

const std::string DocIterator::TrafoContainer = R"mydelimiter(
    This class allows to iterate through the transformers of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
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
    This class represents what you get from retrieving the transformers from 
    :class:`lightsim2grid.elements.TrafoContainer`.

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

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_p_lv_mw = R"mydelimiter(
    Get the active power in MW for at the "lv" side of the transformer. If it is positive it means power is absorbed by the transformer.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_hv_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "hv" side of the transformer. If it is positive it means power is absorbed by the transformer.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_lv_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "lv" side of the transformer. If it is positive it means power is absorbed by the transformer.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_hv_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "hv" side of the transformer is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_lv_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "lv" side of the transformer is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_hv_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "hv" side of the transformer is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_lv_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "lv" side of the transformer is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_lv_ka = R"mydelimiter(
    Get the current flows (in kA) at the "lv" side of the transformer.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_hv_ka = R"mydelimiter(
    Get the current flows (in kA) at the "hv" side of the transformer.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::LineContainer = R"mydelimiter(
    This class allows to iterate through the powerlines of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
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
    This class represents what you get from retrieving the powerlines from 
    :class:`lightsim2grid.elements.LineContainer`.

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

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_p_ex_mw = R"mydelimiter(
    Get the active power in MW for at the "ex" side of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_or_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "or" side of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_ex_mvar = R"mydelimiter(
    Get the reactive power in MVAr for at the "ex" side of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_or_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "or" side of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_ex_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "ex" side of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_or_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "or" side of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_ex_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "ex" side of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_or_ka = R"mydelimiter(
    Get the current flows (in kA) at the "or" side of the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_ex_ka = R"mydelimiter(
    Get the current flows (in kA) at the "ex" side of the line.

)mydelimiter" + DocIterator::only_avail_res;


const std::string DocIterator::DCLineContainer = R"mydelimiter(
    This class allows to iterate through the dc lines of the :class:`lightsim2grid.gridmodel.GridModel` easily, as if they were
    in a python list.

    DC lines are modeled as in pandapower and can be represented a the 
    `pandapower dclines <https://pandapower.readthedocs.io/en/latest/elements/dcline.html#electric-model>`_ . Basically
    a dc line is made of 2 generators (one at each side `or` or `ex`). These
    two generators are linked together: if one produces xx MW the other one consume yy MW (and there 
    exists a relation between xx and yy).

    To dive a bit into the modelling, the two underlying generators can be controlled independantly for the voltage
    setpoint. But they are linked together for their active value. A dc powerline has some losses (both in MW and in 
    percent) and the formula for computing the power injected / produced by each generator is:

    - if `xx` is positive, then `yy = -1.0 * (xx - loss_mw) * (1.0 - 0.01 * loss_percent)`
    - if `xx` is negative, then `yy = -1.0 * xx / (1.0 - 0.01 * loss_percent) + loss_mw

    The first formula directly comes from pandapower. The second one ensures that if the direction of the flow is 
    inverted, then flows should also be inverted (`xx` becomes `yy` and reciprocally).

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the dc powerlines (usually there are none...)
        for dcline in grid_model.get_dclines():
            # do something with line !
            dcline.bus_or_id

        print(f"There are {len(grid_model.get_dclines())} lines on the grid.")

        // first_dcline = grid_model.get_dclines()[0]

    You can have a look at :class:`lightsim2grid.elements.DCLineInfo` for properties of these elements.

)mydelimiter";


const std::string DocIterator::DCLineInfo = R"mydelimiter(
    This class represents what you get from retrieving the dc powerlines from 
    :class:`lightsim2grid.elements.DCLineContainer`.

    It allows to read information from each dc powerline of the powergrid.

    DC Powerlines have two sides, one is "or" for "origin" and one is "ex" for "extremity" that are connected and linked to each other
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
        first_line = grid_model.get_dclines()[0]  # first dcline, this is a `DCLineInfo`
        for dcline in grid_model.get_dclines():
            # dcline is a `LineInfo`
            dcline.bus_or_id

    Notes
    -----
    DC lines are modeled by two underlying generators.
    
    Each of them can be controlled independantly for the voltage setpoint. 
    
    But they are linked together for their active value. A dc powerline has some losses (both in MW `loss_mw` and in 
    percent `loss_percent`) and the formula for computing the power injected / produced by each generator is:

    - if `xx` is positive, then `yy = -1.0 * (xx - loss_mw) * (1.0 - 0.01 * loss_percent)`
    - if `xx` is negative, then `yy = -1.0 * xx / (1.0 - 0.01 * loss_percent) + loss_mw

    The first formula directly comes from pandapower. The second one ensures that if the direction of the flow is 
    inverted, then flows should also be inverted (`xx` becomes `yy` and reciprocally).

    For the sake of simplicity, you can only control the active value at the `or` side of the dc powerline. The
    active value at the `ex` side 
)mydelimiter";

const std::string DocIterator::target_p_or_mw = R"mydelimiter(
    The target active production setpoint at the `or` side of the dc powerline (in MW).

    .. note:: 
        If it's positive it means that power is actually injected at the `or` side, so the
        power flows from `ex` to `or`.
        
)mydelimiter";

const std::string DocIterator::target_vm_or_pu = R"mydelimiter(
    The target active voltage setpoint at the `or` side of the powerline (in pu NOT in kV).
)mydelimiter";

const std::string DocIterator::target_vm_ex_pu = R"mydelimiter(
    The target active voltage setpoint at the `ex` side of the powerline (in pu NOT in kV).
)mydelimiter";

const std::string DocIterator::dc_line_formula = R"mydelimiter(
    .. note::
        A DC line is modeled by two connected generators and some losses to convert the power from one to the other.

        Two types of losses are considered:
        
        - `flat` losses: `loss_mw` (in MW) 
        - `relative` losses: `loss_percent` (no unit) 
        
        The formula for computing the power injected / produced by each generator is:

        - if `or_mw` is positive, then `ex_mw = -1.0 * (or_mw - loss_mw) * (1.0 - 0.01 * loss_percent)`
        - if `or_mw` is negative, then `ex_mw = -1.0 * or_mw / (1.0 - 0.01 * loss_percent) + loss_mw

        Where `or_mw` denotes the power injected at the origin side and `ex_mw` the power injected at the `extremity`
        side.

    .. note::
        By convention, a dc powerline adopts the `load convention`.

        This means that if `or_mw` is positive then power is consumed at the `or` side, so the
        power flows from `or` to `ex` and vice versa.

)mydelimiter";

const std::string DocIterator::loss_pct = R"mydelimiter(
    The value of the `loss percent` parameter for the dc line.

)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocIterator::loss_mw = R"mydelimiter(
    The value of the `loss_mw` parameter for the dc line.

)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocIterator::res_p_or_mw_dcline = R"mydelimiter(
    The amount of active power injected at the `or` side of the dc powerline (in MW).

)mydelimiter" + DocIterator::only_avail_res + DocIterator::dc_line_formula;

const std::string DocIterator::res_p_ex_mw_dcline = R"mydelimiter(
    The amount of active power injected at the `ex` side of the dc powerline (in MW).

)mydelimiter" + DocIterator::only_avail_res + DocIterator::dc_line_formula;

const std::string DocIterator::res_q_or_mvar_dcline = R"mydelimiter(
    The amount of reactive power injected at the `or` side of the dc powerline (in MVAr).

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_ex_mvar_dcline = R"mydelimiter(
    The amount of reactive power injected at the `ex` side of the dc powerline (in MVAr).

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_or_kv_dcline = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "or" side of the dc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_ex_kv_dcline = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which this "ex" side of the dc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_or_deg_dcline = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "or" side of the dc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_ex_deg_dcline = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which this "ex" side of the dc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::gen_or = R"mydelimiter(
    Direct access to the "or" generators, directly returns a :class:`lightsim2grid.elements.GenInfo`
)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocIterator::gen_ex = R"mydelimiter(
    Direct access to the "ex" generators, directly returns a :class:`lightsim2grid.elements.GenInfo`
)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocGridModel::GridModel = R"mydelimiter(
    This class represent a lightsim2grid power network. All the elements that can be manipulated by
    lightsim2grid are represented here.

    We do not recommend to use this class directly, but rather to use a :class:`lightsim2grid.LightSimBackend.LightSimBackend`.

    Examples
    ---------

    We **DO NOT** recommend to do:

    .. code-block:: python

        import lightsim2grid
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower network for example pp_net = pn.case118() 

        grid_model = init(pp_net)

    It's better to do:

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # any grid2op environment
        grid2op_env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = grid2op_env.backend._grid

    The best way to use this class is through the `LightSimBackend` and not to use it directly !

)mydelimiter";

const std::string DocGridModel::change_solver =  R"mydelimiter(
    This function allows to control which solver is used during the powerflow. See the section :ref:`available-powerflow-solvers` for 
    more information about them.

    .. seealso:: :attr:`lightsim2grid.solver.SolverType` for a list of the available solver (NB: some solvers might not be available on all platform)

    .. note::
        If the solver type entered is a `DC` solver (**eg** from :attr:`lightsim2grid.solver.SolverType`, 
        `DC`, `KLUDC` or `NICSLUDC`), it will change the `_dc_solver` otherwise the regular `_solver` 
        is modified.

    Examples
    ---------

    .. code-block:: python
        
        from lightsim2grid.solver import SolverType
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # change the solver used for the powerflow
        # to use internally a solver based on Newton Raphson algorithme using Eigen sparse LU
        lightsim_grid_model.change_solver(SolverType.SparseLUSolver)  

)mydelimiter";

const std::string DocGridModel::available_solvers =  R"mydelimiter(
    Return the list of solver available on the current lightsim2grid installation.

    This is a list of :attr:`lightsim2grid.solver.SolverType`.

)mydelimiter";
const std::string DocGridModel::get_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver when performing a powerflow.

    This is equivalent to the `get_computation_time` of the :func:`lightsim2grid.solver.AnySolver.get_computation_time` of
    the solver used (:func:`lightsim2grid.gridmodel.GridModel.get_solver`)
    
)mydelimiter";
const std::string DocGridModel::get_dc_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver (used to perform DC approximation) when performing a DC powerflow.

    This is equivalent to the `get_computation_time` of the :func:`lightsim2grid.solver.AnySolver.get_computation_time` of
    the DC solver used (:func:`lightsim2grid.gridmodel.GridModel.get_dc_solver`)
    
)mydelimiter";
const std::string DocGridModel::get_solver_type = R"mydelimiter(
    Return the type of the solver currently used.

    This is equivalent to the `get_type` of the :func:`lightsim2grid.solver.AnySolver.get_type` of
    the solver used.

)mydelimiter";
const std::string DocGridModel::get_dc_solver_type = R"mydelimiter(
    Return the type of the solver currently used to compute DC powerflow.

)mydelimiter";
const std::string DocGridModel::get_solver = R"mydelimiter(
    Return the solver currently in use as a :func:`lightsim2grid.solver.AnySolver` instance.

)mydelimiter";
const std::string DocGridModel::get_dc_solver = R"mydelimiter(
    Return the solver currently in use as a :func:`lightsim2grid.solver.AnySolver` instance for the dc powerflow.

)mydelimiter";

const std::string DocGridModel::get_lines = R"mydelimiter(
    This function allows to retrieve the powerlines (as a 
    :class:`lightsim2grid.elements.LineContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the powerlines
        print([el.x_pu for el in lightsim_grid_model.get_lines()]) # to print the "x" for each powerlines

)mydelimiter";
const std::string DocGridModel::get_trafos = R"mydelimiter(
    This function allows to retrieve the transformers (as a 
    :class:`lightsim2grid.elements.LineContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the trafos
        print([el.x_pu for el in lightsim_grid_model.get_trafos()]) # to print the "x" for each transformer

)mydelimiter";
const std::string DocGridModel::get_generators = R"mydelimiter(
    This function allows to retrieve the (standard) generators (as a 
    :class:`lightsim2grid.elements.GeneratorContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the generators
        print([el.target_p_mw for el in lightsim_grid_model.get_generators()]) # to print the active production setpoint for each generators

)mydelimiter";
const std::string DocGridModel::get_static_generators = R"mydelimiter(
    This function allows to retrieve the (more exotic) static generators (as a 
    :class:`lightsim2grid.elements.SGenContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the static generators
        print([el.target_p_mw for el in lightsim_grid_model.get_static_generators()]) # to print the active production setpoint for each static generator

)mydelimiter";
const std::string DocGridModel::get_shunts = R"mydelimiter(
    This function allows to retrieve the shunts (as a 
    :class:`lightsim2grid.elements.ShuntContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the shunts
        print([el.target_q_mvar for el in lightsim_grid_model.get_shunts()]) # to print the reactive consumption for each shunts

)mydelimiter";
const std::string DocGridModel::get_storages = R"mydelimiter(
    This function allows to retrieve the storage units (as a 
    :class:`lightsim2grid.elements.LoadContainer` object,
    see :ref:`elements-modeled` for more information)

    .. note::
        We want to emphize that, as far as lightsim2grid is concerned, the storage units are modeled as loads. This is why
        this function will return a :class:`lightsim2grid.elements.LoadContainer`.

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # print the target consumption of each storage units
        print([el.target_p_mw for el in lightsim_grid_model.get_storages()]) # to print the active consumption for each storage unit

)mydelimiter";
const std::string DocGridModel::get_loads = R"mydelimiter(
    This function allows to retrieve the loads (as a :class:`lightsim2grid.elements.LoadContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # print the target consumption of each loads
        print([el.target_p_mw for el in lightsim_grid_model.get_loads()]) # to print the active consumption for each load

)mydelimiter";

const std::string DocGridModel::get_dclines = R"mydelimiter(
    This function allows to retrieve the dc powerlines (as a 
    :class:`lightsim2grid.elements.DCLineContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.gridmodel import init
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the powerlines
        print([el.x_pu for el in lightsim_grid_model.get_dclines()]) # to print the "x" for each powerlines

)mydelimiter";

const std::string DocGridModel::_internal_do_not_use = R"mydelimiter(
        INTERNAL

        .. warning:: /!\\ Internal, do not use unless you know what you are doing /!\\

        This is used as part of a dedicated code for :class:`lightsim2grid.LightSimBackend.LightSimBackend`

)mydelimiter";

const std::string DocGridModel::J_description = R"mydelimiter(
    J has the shape::
    
        | s | slack_bus |               | (pvpq+1,1) |   (1, pvpq)  |  (1, pq)   |
        | l |  -------  |               |            | ------------------------- |
        | a | J11 | J12 | = dimensions: |            | (pvpq, pvpq) | (pvpq, pq) |
        | c | --------- |               |   ------   | ------------------------- |
        | k | J21 | J22 |               |  (pq, 1)   |  (pq, pvpq)  | (pq, pq)   |
        

    With:
    
    - `J11` = dS_dVa[array([pvpq]).T, pvpq].real (= real part of dS / dVa for all pv and pq buses)
    - `J12` = dS_dVm[array([pvpq]).T, pq].real
    - `J21` = dS_dVa[array([pq]).T, pvpq].imag
    - `J22` = dS_dVm[array([pq]).T, pq].imag (= imaginary part of dS / dVm for all pq buses)
    - `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real 
      and dS_dVm[slack_bus_id, pq].real
    - `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
      the remaining pq components are all 0.

    .. note::
        By default (and this cannot be changed at the moment), all buses in `ref` will be pv buses except the first one.

    .. note::
        the notation `pvpq` above means "the concatenation of the pv vector and the pq vector" (after the distributed slack is taken into account - see note just above)
)mydelimiter";

const std::string DocGridModel::get_J_python = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.

    .. note::
        Some powerflows (*eg* DC or Gauss Seidel) do not rely on jacobian matrix, in this case, calling this function will return an exception. 
)mydelimiter" + DocGridModel::J_description;

const std::string DocGridModel::get_Va = R"mydelimiter(
    Returns the voltage angles for each buses as a numpy vector of real number. This vector have the size of the total number of active buses on the system.

    You can use the :attr:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` (or :attr:`lightsim2grid.gridmodel.GridModel.id_dc_solver_to_me`) to know at which bus
    (on the grid) they corresponds.

)mydelimiter";

const std::string DocGridModel::get_Vm = R"mydelimiter(
    Returns the voltage magnitude for each buses as a numpy vector of real number. This vector have the size of the total number of active buses on the system.

    You can use the :attr:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` (or :attr:`lightsim2grid.gridmodel.GridModel.id_dc_solver_to_me`) to know at which bus
    (on the grid) they corresponds.
)mydelimiter";

const std::string DocGridModel::get_V = R"mydelimiter(
    Returns the complex voltage for each buses as a numpy vector of complex number. This vector have the size of the total number of active buses on the system.

    You can use the :attr:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` (or :attr:`lightsim2grid.gridmodel.GridModel.id_dc_solver_to_me`) to know at which bus
    (on the grid) they corresponds.
)mydelimiter";


const std::string DocGridModel::id_me_to_ac_solver = R"mydelimiter(
    In lightsim2grid, buses are labelled from `0` to `n-1` (if `n` denotes the total number of buses on the grid) [this is called "**grid model bus id**"]

    At any given point in time, some buses might be deactivated (for example because nothing is connected to them).

    On the other end, the solvers need a contiguous list of only active buses (otherwise they might run into divergence issue) [this will be called 
    "**solver bus id**" later on]

    This function allows, for all buses of the :class:`lightsim2grid.gridmodel.GridModel` to know on which "solver bus" they are affected. It
    has the same size as the total number of buses on the grid. And for each of them it tells to which "solver bus" it is connected (unless there is a `-1`,
    meaning the associated bus is deactivated).

    Examples
    ---------

    .. code-block:: python

        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimbackend())
        grid_model = env.backend._grid

        id_me_to_ac_solver = grid.id_me_to_ac_solver()
        # is [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]

        # put everything to bus 2 on substation O
        _ = env.step(env.action_space({"set_bus": {"substations_id": [(0, (2, 2, 2))]}}))

        id_me_to_ac_solver2 = grid.id_me_to_ac_solver()
        # is [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    
    .. seealso:: :class:`lightsim2grid.gridmodel.GridModel.id_me_to_dc_solver` for its counterpart when a dc powerflow is used
    
    .. seealso:: :class:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for the "reverse" operation (given a "solver bus" id, returns 
        the "gridmodel bus id")

    Notes
    -----

    For all steps, you have the propertie that, if `id_ac_solver_to_me = gridmodel.id_ac_solver_to_me()` and `id_me_to_ac_solver = gridmodel.id_me_to_ac_solver()`
    and by denoting `gridmodel_bus_id = np.arange(gridmodel.total_bus())` and `solver_bus_id = np.arange(gridmodel.nb_bus())`:

    - `solver_bus_id` and `id_ac_solver_to_me` have the same shape
    - `gridmodel_bus_id` and `id_me_to_ac_solver` have the same shape
    - `solver_bus_id` is shorter (or of the same length) than `gridmodel_bus_id`
    - the connected bus (in the grid model) are given by `gridmodel_bus_id[id_ac_solver_to_me]`, and it gives their order

)mydelimiter";

const std::string DocGridModel::id_ac_solver_to_me = R"mydelimiter(
    In lightsim2grid, buses are labelled from `0` to `n-1` (if `n` denotes the total number of buses on the grid) [this is called "**grid model bus id**"]

    At any given point in time, some buses might be deactivated (for example because nothing is connected to them).

    On the other end, the solvers need a contiguous list of only active buses (otherwise they might run into divergence issue) [this will be called 
    "**solver bus id**" later on]

    This function allows, for all buses exported in the solver, to retrieve which was the initial bus in the :class:`lightsim2grid.gridmodel.GridModel`. It
    has the same size as the number of active buses on the grid.

    Examples
    ---------

    .. code-block:: python

        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimbackend())
        grid_model = env.backend._grid

        id_ac_solver_to_me = grid.id_ac_solver_to_me()
        # is [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

        # put everything to bus 2 on substation O
        _ = env.step(env.action_space({"set_bus": {"substations_id": [(0, (2, 2, 2))]}}))

        id_ac_solver_to_me2 = grid.id_ac_solver_to_me()
        # is [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    
    .. seealso:: :class:`lightsim2grid.gridmodel.GridModel.id_dc_solver_to_me` for its counterpart when a dc powerflow is used
    
    .. seealso:: :class:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` for the "reverse" operation (given a "solver bus" id, returns 
        the "gridmodel bus id")

    Notes
    -----

    For all steps, you have the propertie that, if `id_ac_solver_to_me = gridmodel.id_ac_solver_to_me()` and `id_me_to_ac_solver = gridmodel.id_me_to_ac_solver()`
    and by denoting `gridmodel_bus_id = np.arange(gridmodel.total_bus())` and `solver_bus_id = np.arange(gridmodel.nb_bus())`:

    - `solver_bus_id` and `id_ac_solver_to_me` have the same shape
    - `gridmodel_bus_id` and `id_me_to_ac_solver` have the same shape
    - `solver_bus_id` is shorter (or of the same length) than `gridmodel_bus_id`
    - the connected bus (in the grid model) are given by `gridmodel_bus_id[id_ac_solver_to_me]`, and it gives their order

)mydelimiter";

const std::string DocGridModel::id_me_to_dc_solver = R"mydelimiter(
    Same as :class:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` but only used for the DC approximation.
)mydelimiter";

const std::string DocGridModel::id_dc_solver_to_me = R"mydelimiter(
    Same as :class:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` but only used for the DC approximation.
)mydelimiter";

const std::string DocGridModel::total_bus = R"mydelimiter(
    Returns (>0 integer) the total number of buses in the powergrid (both connected and disconnected)
)mydelimiter";

const std::string DocGridModel::nb_bus = R"mydelimiter(
    Returns (>0 integer) the number of connected buses on the powergrid (ignores the disconnected bus).
)mydelimiter";

const std::string DocGridModel::get_pv = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PV" (ie the buses on which at least a generator is connected.).

    It returns a vector of integer.

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
        
)mydelimiter";

const std::string DocGridModel::get_pq = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PQ".

    It returns a vector of integer.

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it will might be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocGridModel::get_slack_ids = R"mydelimiter(
    Returns the ids of the buses that are part of the distributed slack.

    It returns a vector of integer.

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocGridModel::get_slack_weights = R"mydelimiter(
    For each bus used by the solver, it outputs its participation to the distributed slack.

    It's 0 if the current bus does not participate to it, otherwise it is made of > 0. real numbers.

    This vector sums to 1 and has the same size as the number of active buses on the grid.

    .. warning:: 
        This vector represents "solver buses" and not "original grid model buses".

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocGridModel::get_Ybus = R"mydelimiter(
    This function returns the (complex) `Ybus` matrix used to compute the powerflow.

    The resulting matrix is a CSC scipy sparse matrix of complex number.

    It is a square matrix, as many rows (columns) as there are connected buses on the grid.

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system !

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Ybus[id_me_to_ac_solver[k],:]` (rows of this bus), `Ybus[:, id_me_to_ac_solver[k]]` (column for this bus) 

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter";

const std::string DocGridModel::get_dcYbus = R"mydelimiter(
    It is the equivalent of :func:`lightsim2grid.gridmodel.GridModel.get_Ybus` but for the dc solver.

    .. warning::
        As opposed to some other librairies (for example Matpower of pandapower), the Ybus for the dc approximation in lightsim2grid has no
        imaginary components. 
        
        It could have returned a real matrix, but we choose (out of consistency with other solvers) to keep the representation
        as a complex numbers.
    
)mydelimiter"; 

const std::string DocGridModel::get_Sbus = R"mydelimiter(
    This function returns the (complex) `Sbus` vector, which is the vector of active / reactive power injected at each active bus

    The resulting vector is a vector of complex number having the size of the number of connected buses on the grid.

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system and in load convention (so generation will be negative)

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
    
    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Sbus[id_me_to_ac_solver[k]]` is the total power injected at the grid model bus solver `k`.

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter"; 

const std::string DocGridModel::check_solution = R"mydelimiter(
    This function allows to check that a given complex voltage vector satisfies the KCL or not, given the state of the sytem.

    .. note::
        It is expected that you provide a complex number even for the buses that are disconnected in the grid model. They will not be ignored
        so you can put anything you want. We keep the public interface this way to avoid headaches with the bus order between
        the grid model and the solver (you can refer to :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and 
        :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` if you still want to have a look)

    .. seealso:: :class:`lightsim2grid.physical_law_checker.PhysicalLawChecker` for an easier to use, more pythonic function !

    Parameters
    ------------
    V:
      It expects a complex voltage vector (having as many components as the total number of buses in the grid.) representing the
      vector you want to test.

    check_q_limits: ``bool``
      whether you want to take into account the reactive limit of generators when performing the check 

    Returns
    -------
    mismatch: 
        A complex vector having the size of the number of total buses on the grid, given, for each of them, the active / reactive power mismatch
        at each bus (ie the power you would need to take from the grid and have the input vector `V` checking the KCL given the current state of
        the grid)

)mydelimiter"; 

const std::string DocGridModel::deactivate_result_computation = R"mydelimiter(
    Allows to deactivate the computation of the flows, reactive power absorbed by generators etc. to gain a bit of time when it is not needed.

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.reactivate_result_computation`
)mydelimiter";     

const std::string DocGridModel::reactivate_result_computation = R"mydelimiter(
    Allows to reactivate the computation of the flows, reactive power absorbed by generators etc. when they are needed again after having been
    deactivated.

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.deactivate_result_computation`
)mydelimiter";     

const std::string DocGridModel::ac_pf = R"mydelimiter(
    Allows to perform an AC (alternating current) powerflow.

    .. note::
        It is expected that you provide a complex number even for the buses that are disconnected in the grid model. They will not be affected (if the powerflow converges)
        and you can put anything you want there. We keep the public interface this way to avoid headaches with the bus order between
        the grid model and the solver (you can refer to :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and 
        :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` if you still want to have a look)

    .. seealso:: :func:`lightsim2grid.gridmodel.GridModel.dc_pf` if you want to perform DC powerflow (same interface, same results, same behaviour)

    .. warning::
        The input vector `V` is modified (and is equal to the resulting vector `V`)

    Parameters
    ------------
    V:
      It expects a complex voltage vector (having as many components as the total number of buses in the grid.) representing the
      initial guess of the resulting flows. This vector will be modified !

    max_iter: ``int``
        Maximum number of iterations allowed (this might be ignored) and should be a >= 0 integer
    
    tol: ``float``
        Tolerance criteria to stop the computation. This should be > 0 real number.

    Returns
    -------
    V:
        A complex vector given the complex voltage at each buses of the grid model. Will be empty when the powerflow diverged.

    Examples
    --------

    .. code-block:: python

        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimbackend())
        grid_model = env.backend._grid

        V = grid_model.ac_pf(V, 10, 1e-8)
        # if the powerflow has converged, V.shape > 0 otherwise V is empty (size 0)
        # the original V is modified in the process !

)mydelimiter";     

const std::string DocGridModel::dc_pf = R"mydelimiter(
    This function has the same interface, inputs, outputs, behaviour, etc. as the :func:`lightsim2grid.gridmodel.GridModel.ac_pf`.
)mydelimiter";       

const std::string DocComputers::Computers = R"mydelimiter(
    Allows the computation of time series, that is, the same grid topology is used while the active / reactive power injected
    at each buse vary. The grid topology is fixed, the injections vary.

    This is a "raw" c++ class, for an easier to use interface, please refer to the python documentation of the 
    :class:`lightsim2grid.timeSerie.TimeSerie` class.

)mydelimiter";

const std::string DocComputers::total_time = R"mydelimiter(
    Total time spent in solving the powerflows, pre processing the data, post processing them, initializing everything etc.
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocComputers::solver_time = R"mydelimiter(
    Total time spent only in solving the powerflows (excluding pre processing the data, post processing them, initializing everything etc.)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocComputers::amps_computation_time = R"mydelimiter(
    Time spent in computing the flows (in amps) after the voltages have been computed at each nodes
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocComputers::preprocessing_time = R"mydelimiter(
    Time spent in pre processing the data (this involves, but is not limited to the computation of the Sbus)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocComputers::nb_solved = R"mydelimiter(
    Total number of powerflows solved.

)mydelimiter";

const std::string DocComputers::get_status = R"mydelimiter(
    Status of the solvers (1: success, 0: failure).

    .. note::
        Even if the solver failed at some point, some results might still be available (but not totally).

)mydelimiter";

const std::string DocComputers::compute_Vs = R"mydelimiter(
    Compute the voltages (at each bus of the grid model) for some time series of injections (productions, loads, storage units, etc.)

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.Computers.compute_flows` and 
        :func:`lightsim2grid.timeSerie.Computers.get_flows`, :func:`lightsim2grid.timeSerie.Computers.get_voltages` or
        :func:`lightsim2grid.timeSerie.Computers.get_sbuses`.

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

    Parameters
    -----------
    gen_p:  ``numy.ndarray``, float
        Active generation for each generators. Its counts as many column as the number of generators on the grid and as many rows as
        the number of steps to compute.

    sgen_p:  ``numy.ndarray``, float
        Active generation for each static generator. Its counts as many column as the number of static generators on the grid and as many rows as
        the number of steps to compute.

    load_p:``numy.ndarray``, float
        Active consumption for each loads. Its counts as many column as the number of loads on the grid and as many rows as
        the number of steps to compute.

    load_q: ``numy.ndarray``, float
        Reactive consumption for each loads. Its counts as many column as the number of loads on the grid and as many rows as
        the number of steps to compute.

    Vinit:  ``numy.ndarray``, complex
        First voltage at each bus of the grid model (including the disconnected buses)

    max_iter:  ``int``
        Total number of iteration (>0 integer)

    tol: ``float``
        Solver tolerance (> 0. float)

    Returns
    ----------
    status: ``int``
        The status of the computation. 1 means "success": all powerflows were computed sucessfully, 0 means there were some errors and that 
        the computation stopped after a certain number of steps.

)mydelimiter";

const std::string DocComputers::compute_flows = R"mydelimiter(
    Retrieve the flows (in amps, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.Computers.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocComputers::compute_power_flows = R"mydelimiter(
    Retrieve the active flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.Computers.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocComputers::get_flows = R"mydelimiter(
    Get the current flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a time step, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.timeSerie.Computers.compute_Vs` has been caleed)

    Returns
    -------
    As: ``numpy.ndarry`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocComputers::get_power_flows = R"mydelimiter(
    Get the active flows (in MW) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a time step, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_power_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.timeSerie.Computers.compute_Vs` has been caleed)

    Returns
    -------
    As: ``numpy.ndarry`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocComputers::get_voltages = R"mydelimiter(
    Get the complex voltage angles at each bus of the powergrid.

    Each rows correspond to a time step, each column to a bus.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_Vs`.

    Returns
    -------
    Vs: ``numpy.ndarry`` (matrix)
        The complex voltage angles at each bus of the powergrid.

)mydelimiter";

const std::string DocComputers::get_sbuses = R"mydelimiter(
    Get the complex power injected at each (solver id) bus of the powergrid. Results are given in pair unit.
    We do not recommend to use it as it uses the solver id and NOT the powergrid bus id (you can refer to 
    :func:`lightsim2grid.gridmodel.GridModel.id_me_to_ac_solver` and 
    :func:`lightsim2grid.gridmodel.GridModel.id_ac_solver_to_me` for more information)

    Each rows correspond to a time step, each column to a bus (bus are identified by their solver id !)

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.Computers.compute_Vs`.

    Returns
    -------
    Sbuses: ``numpy.ndarry`` (matrix)
        The complex power injected at each bus (pair unit, load sign convention)

)mydelimiter";

const std::string DocComputers::clear = R"mydelimiter(
    Clear the solver and to as if the class never performed any powerflow.

)mydelimiter";

const std::string DocSecurityAnalysis::SecurityAnalysis = R"mydelimiter(
    Allows the computation of "security analysis", that consists in computing the flows that would result from the disconnection of one or multiple
    disconnections of some powerlines.

    This is a "raw" c++ class, for an easier to use interface, please refer to the python documentation of the 
    :class:`lightsim2grid.securityAnalysis.SecurityAnalysis` class.

    .. warning::
        This function might give wrong result for lightsim2grid version 0.5.5 were they were a bug : when some contingencies made the grid
        non connex, it made all the other contingencies diverge. This bug has been fixed in version 0.6.0 and this is why we do not recommend
        to use this feature with lightsim2grid version < 0.6.0 !

    .. note::
        Even if you instruct it to simulate the same contingency multiple times, it will only do it once.

    .. note::
        You can only simulate disconnection of powerlines / transformers

    At a glance, this class should be used in three steps:

    1) Modify the list of contingencies to simulate, with the functions:

    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_n1`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_all_n1`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_nk`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_multiple_n1`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.clear`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.remove_n1`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.remove_nk`
    - :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.remove_multiple_n1`

    2) Then you can start the computation of the security analysis with 
    :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute` then optionally
    :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute_flows` .

    3) And finally inspect the results with :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_flows` and 
    :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_voltages` .

)mydelimiter";

const std::string DocSecurityAnalysis::preprocessing_time = R"mydelimiter(
    Time spent in pre processing the data (this involves, the checking whether the grid would be still connex after the contingency for example)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocSecurityAnalysis::modif_Ybus_time = R"mydelimiter(
    Time spent to modify the Ybus matrix before simulating each contingency.
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocSecurityAnalysis::add_all_n1 = R"mydelimiter(
    This allows to add all the "n-1" in the contingency list to simulate.

    .. seealso:: :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_n1` to add only a single line

    .. seealso::
        :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call 
        (but not necessarily all)

)mydelimiter";

const std::string DocSecurityAnalysis::add_n1 = R"mydelimiter(
    This allows to add a single  "n-1" in the contingency list to simulate.

    .. seealso:: :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_all_n1` to add all contingencies at the same time

    .. seealso::
        :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call.

    Parameters
    ----------
    line_id: ``int``
        The line id you would like to see disconnected
        
)mydelimiter";

const std::string DocSecurityAnalysis::add_nk = R"mydelimiter(
    This allows to add a single  "n-k" in the contingency list to simulate (it will only add at most one contingency)

    .. warning::
        A "n-k" will disconnect multiple powerlines at the same time. It's not the same as adding muliple "n-1" contingencies, where
        powerlines will be disconnected one after the other.

    Parameters
    ----------
    vect_nk: ``list`` (of ``int``)
        The lines id you want to add in the single contingency added.

)mydelimiter";

const std::string DocSecurityAnalysis::add_multiple_n1 = R"mydelimiter(
    This allows to add a multiple "n-1" in the contingency list to simulate (it will add as many contingency as the size of the list)
    and is equivalent to call multiple times :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_n1`

    .. seealso::
        :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.add_all_n1` to add all the "n-1" contingencies.
    
    .. warning::
        A "n-k" will disconnect multiple powerlines at the same time. It's not the same as adding muliple "n-1" contingencies, where
        powerlines will be disconnected one after the other.

    Parameters
    ----------
    vect_n1: ``list`` (of ``int``)
        The lines id you want to add to the contingency list

)mydelimiter";

const std::string DocSecurityAnalysis::clear = R"mydelimiter(
    Clear the list of all contingencies. After a call to this method, you will need to re add some contingencies with

)mydelimiter";

const std::string DocSecurityAnalysis::remove_n1 = R"mydelimiter(
    Remove a single "n-1" contingency from the contingency list to simulate.

    Parameters
    ----------
    line_id: ``int``
        The line id you would like to remove from contingency list (will remove a single "n-k" contingencies)
    
    Returns
    -------
    success: ``bool``
        Whether or not the contingency has been properly removed

)mydelimiter";

const std::string DocSecurityAnalysis::remove_nk = R"mydelimiter(
    Remove a single "n-k" contingency from the contingency list to simulate. This removes at much one single contingency

    Parameters
    ----------
    vect_nk: ``list`` (of ``int``)
        The lines id you want to remove from contingency list.

    Returns
    -------
    nb_removed: ``int``
        The total number of contingencies removed from the contingency list

)mydelimiter";

const std::string DocSecurityAnalysis::remove_multiple_n1 = R"mydelimiter(
    Remove multiple "n-1" contingency from the contingency list to simulate. This can remove up to `len(vect_n1)` single contingencies
    from the contingency list.

    Parameters
    ----------
    vect_n1: ``list`` (of ``int``)
        The lines id you want to remove from contingency list (will remove multiple "n-1" single contingency)

    Returns
    -------
    success: ``bool``
        Whether or not the contingency has been properly removed

)mydelimiter";

const std::string DocSecurityAnalysis::my_defaults_vect = R"mydelimiter(
    Allows to inspect the contingency list that will be simulated.

    Returns
    -------
    my_defaults_vect: ``list``
        The list (of list) of all the current contingencies. Its length corresponds to the number of contingencies simulated.
        For each contingency, it gives which powerline will be disconnected.

)mydelimiter";

const std::string DocSecurityAnalysis::compute = R"mydelimiter(
    Compute the voltages (at each bus of the grid model) for some time series of injections (productions, loads, storage units, etc.)

    .. note::
        This function must be called before :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute_flows` and 
        :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_flows` or
        :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_voltages` .

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

    Parameters
    -----------
    Vinit:  ``numy.ndarray``, complex
        First voltage at each bus of the grid model (including the disconnected buses)

    max_iter:  ``int``
        Total number of iteration (>0 integer)

    tol: ``float``
        Solver tolerance (> 0. float)

)mydelimiter";

const std::string DocSecurityAnalysis::compute_flows = R"mydelimiter(
    Compute the current flows (in amps, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocSecurityAnalysis::compute_power_flows = R"mydelimiter(
    Compute the current flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocSecurityAnalysis::get_flows = R"mydelimiter(
    Get the flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a contingency, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute` has been caleed)

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.securityAnalysis.SecurityAnalysis.get_flows` method.

    Returns
    -------
    As: ``numpy.ndarray`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocSecurityAnalysis::get_voltages = R"mydelimiter(
    Get the complex voltage angles at each bus of the powergrid.

    Each rows correspond to a contingency, each column to a bus.

    .. warning::
        This function must be called after :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute`.

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.securityAnalysis.SecurityAnalysis.get_flows` method.

    Returns
    -------
    Vs: ``numpy.ndarray`` (matrix)
        The complex voltage angles at each bus of the powergrid.

)mydelimiter";

const std::string DocSecurityAnalysis::get_power_flows = R"mydelimiter(
    Get the active flows (in MW) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a contingency, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute_power_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.compute` has been caleed)

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.securityAnalysis.SecurityAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.securityAnalysis.SecurityAnalysis.get_flows` method.

    Returns
    -------
    As: ``numpy.ndarray`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";
