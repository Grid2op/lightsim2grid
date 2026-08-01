// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// containst he help message of some common functions (not to write them dozens of time)

#include "help_fun_msg.hpp"

namespace ls2g {

const std::string DocSolver::get_J_python = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.

    The "jacobian" matrix is only available for some powerflow (the one based on the Newton Raphson algorithm)
    and we provide it only for the last computed iteration.

    .. note::
        It is using the "solver" labelling, as this is accessed from the solvers. Unlike
        :func:`get_Va`/:func:`get_Vm`, the Jacobian has no "gridmodel" labelled equivalent on
        :class:`lightsim2grid.network.LSGrid` -- only :func:`lightsim2grid.network.LSGrid.get_J_solver`,
        which keeps the solver labelling.

    .. seealso::
        This function should be equal to :func:`lightsim2grid.network.LSGrid.get_J_solver`

)mydelimiter";

const std::string DocSolver::get_Va = R"mydelimiter(
    Returns the voltage angles for each buses as a numpy vector of real number.

    .. note::
        It is using the "solver" labelling, as this is accessed from the solvers.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.get_Va` for the same things, but rather using 
        the "gridmodel" labelling.

    .. seealso::
        This function should be equal to :func:`lightsim2grid.network.LSGrid.get_Va_solver`

)mydelimiter";

const std::string DocSolver::get_Vm = R"mydelimiter(
    Returns the voltage magnitude for each buses as a numpy vector of real number.

    .. note::
        It is using the "solver" labelling, as this is accessed from the solvers.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.get_Vm` for the same things, but rather using 
        the "gridmodel" labelling.

    .. seealso::
        This function should be equal to :func:`lightsim2grid.network.LSGrid.get_Vm_solver`
)mydelimiter";

const std::string DocSolver::get_V = R"mydelimiter(
    Returns the complex voltage for each buses as a numpy vector of complex number.    
    
    .. note::
        It is using the "solver" labelling, as this is accessed from the solvers.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.get_V` for the same things, but rather using 
        the "gridmodel" labelling.

    .. seealso::
        This function should be equal to :func:`lightsim2grid.network.LSGrid.get_V_solver`
)mydelimiter";

const std::string DocSolver::get_error = R"mydelimiter(
    Returns the error encountered by the solver during the last ``compute_pf`` / ``solve`` call,
    as a :class:`lightsim2grid.algorithm.ErrorType` value (``ErrorType.NoError``, ie 0, when nothing
    went wrong).

    .. note::
        Reaching ``max_iter`` without meeting the requested tolerance is itself reported as an
        error here (``ErrorType.TooManyIterations``), so :func:`converged` (which is exactly
        ``get_error() == ErrorType.NoError``) is ``False`` in that case too.

    See :class:`lightsim2grid.algorithm.ErrorType` for the full list of possible values and what each
    one means.
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

    .. note::
        This python-facing method (also available as ``solve``) validates its inputs before doing
        anything else: a non-square ``Ybus``, a size mismatch between ``Ybus`` / ``V`` / ``Sbus`` /
        ``slack_weights``, an out-of-range id in ``slack_ids`` / ``pv`` / ``pq``, a bus listed in more
        than one of them, an empty ``slack_ids``, a negative ``max_iter`` (0 is accepted: it returns
        the pre-iteration state, before any Newton-Raphson / Gauss-Seidel step), or a non-finite or
        non-positive ``tol`` all raise a clean ``RuntimeError`` (or ``IndexError`` for out-of-range
        ids) instead of touching the underlying solver. This validation is skipped on the internal
        C++ code path used by :class:`lightsim2grid.gridmodel.LSGrid` and the batch solvers
        (:class:`~lightsim2grid.contingencyAnalysis.ContingencyAnalysis`,
        :class:`~lightsim2grid.timeSerie.TimeSerie`, security analysis), which build these arrays
        themselves and call the solver many times in a loop: paying this check on every call there
        would be pure overhead, so it is only performed at this python entry point.

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

    .. note::
        This is returned as a plain ``(float, float, float, float)`` tuple, in the order below
        (there are no named attributes on it) -- for named access to a wider set of timers, see
        :func:`lightsim2grid.algorithm.AlgorithmSelector.get_timers_jacobian` instead, which returns
        a :class:`lightsim2grid.algorithm.TimerJac`.

    Returns
    ---------
    timer_Fx_: ``float``
        Time spent to compute the mismatch at the KCL for each bus (both for active and reactive power)

    timer_solve_: ``float``
        Total time spent in the underlying linear solver

    timer_check_: ``float``
        Time spent in checking whether or not the mismatch of the KCL met the specified tolerance

    timer_total_nr_: ``float``
        Total time spent in the solver

)mydelimiter";
    
const std::string DocSolver::NR_SparseLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the default Eigen sparse solver available in Eigen
    for the linear algebra. 

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NR_SparseLU`.
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NR_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NR_SparseLU)` at creation time    
    
    .. note::
        Available on all plateform, this is the default solver used when :class:`lightsim2grid.algorithm.NRSing_KLU`
        is not found (when a "single slack" is detected).

)mydelimiter";

const std::string DocSolver::NRSing_SparseLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, using the default Eigen sparse solver available in Eigen
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.algorithm.NR_SparseLU` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NRSing_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NRSing_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NRSing_SparseLU)` at creation time

    .. note::
        Available on all plateform, this is the default solver used when a distributed slack bus is detected and :class:`lightsim2grid.algorithm.NR_KLU`
        is not found.

)mydelimiter";

const std::string DocSolver::DC_SparseLU =  R"mydelimiter(
    Default implementation of the DC solver, it uses the default Eigen sparse lu decomposition to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `DC_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.DC_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.DC_SparseLU)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

)mydelimiter";

const std::string DocSolver::FDPF_XB_SparseLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the default Eigen sparse lu decomposition for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_XB_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_XB_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_XB_SparseLU)` at creation time

)mydelimiter";

const std::string DocSolver::FDPF_BX_SparseLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (BX version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the default Eigen sparse lu decomposition for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_BX_SparseLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_BX_SparseLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_BX_SparseLU)` at creation time

)mydelimiter";

const std::string DocSolver::NR_KLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster KLU solver available in the SuiteSparse library
    for the linear algebra (can be unavailable if you build lightsim2grid from source). It is usually faster than the :class:`lightsim2grid.algorithm.NR_SparseLU`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NR_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NR_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NR_KLU)` at creation time

    .. note::
        This is the default solver used when a distributed slack bus is detected (when it's available, otherwise see :class:`lightsim2grid.algorithm.NR_SparseLU`).

)mydelimiter";

const std::string DocSolver::NRSing_KLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm,the faster KLU solver available in the SuiteSparse library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.algorithm.NR_KLU`.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NRSing_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NRSing_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NRSing_KLU)` at creation time

    .. note::
        This is the default solver used when available.

)mydelimiter";

const std::string DocSolver::NRRefactorRetry_KLU = R"mydelimiter(
    Same as :class:`lightsim2grid.algorithm.NR_KLU` (Newton Raphson, distributed slack, KLU
    linear solver), except that if a Jacobian refactorize() fails it falls back to a full
    factorize() (reusing the existing symbolic factorization) before giving up, rather
    than reporting an error immediately.

    Use `get_linear_solver_stats()` on the solver to inspect how often factor/refactor
    calls happen and how often the fallback fires (see :class:`lightsim2grid.algorithm.LinearSolverStats`).

)mydelimiter";

const std::string DocSolver::DC_KLU = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster KLU solver available in the SuiteSparse library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (can be unavailable if you build lightsim2grid from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `DC_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.DC_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.DC_KLU)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

)mydelimiter";

const std::string DocSolver::FDPF_XB_KLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the fast KLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_XB_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_XB_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_XB_KLU)` at creation time

)mydelimiter";

const std::string DocSolver::FDPF_BX_KLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (BX version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the fast KLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_BX_KLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_BX_KLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_BX_KLU)` at creation time

)mydelimiter";

const std::string DocSolver::NR_NICSLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It is usually faster than the :class:`lightsim2grid.algorithm.NR_SparseLU`. (requires a build from source)
    
    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NR_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NR_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NR_NICSLU)` at creation time

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::NRRefactorRetry_NICSLU = R"mydelimiter(
    Same as :class:`lightsim2grid.algorithm.NR_NICSLU` (Newton Raphson, distributed slack,
    NICSLU linear solver), except that if a Jacobian refactorize() fails it falls back to
    a full factorize() before giving up, rather than reporting an error immediately. For
    NICSLU, factorize() and refactorize() call the same underlying routine, so this
    fallback is effectively a no-op retry -- included mainly for API symmetry with
    :class:`lightsim2grid.algorithm.NRRefactorRetry_KLU` and
    :class:`lightsim2grid.algorithm.NRRefactorRetry_CKTSO`.

    Use `get_linear_solver_stats()` on the solver to inspect how often factor/refactor
    calls happen and how often the fallback fires (see :class:`lightsim2grid.algorithm.LinearSolverStats`).

)mydelimiter";

const std::string DocSolver::NRSing_NICSLU = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, the faster NICSLU solver available in the NICSLU library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.algorithm.NR_NICSLU` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NRSing_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NRSing_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NRSing_NICSLU)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::DC_NICSLU = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster NICSLU solver available in the NICSLU library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (requires a build from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `DC_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.DC_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.DC_NICSLU)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu
 
)mydelimiter";

const std::string DocSolver::FDPF_XB_NICSLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the fast NICSLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_XB_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_XB_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_XB_NICSLU)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::FDPF_BX_NICSLU =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (BX version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the fast NICSLU library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_BX_NICSLU` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_BX_NICSLU)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_BX_NICSLU)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for nicslu.

    .. note::

        NICSLU is available at https://github.com/chenxm1986/nicslu

)mydelimiter";

const std::string DocSolver::NR_CKTSO = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, allowing for distributed slack and using the faster CKTSO solver available in the CKTSO library
    for the linear algebra (requires a build from source)
    
    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NR_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NR_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NR_CKTSO)` at creation time

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::NRRefactorRetry_CKTSO = R"mydelimiter(
    Same as :class:`lightsim2grid.algorithm.NR_CKTSO` (Newton Raphson, distributed slack,
    CKTSO linear solver), except that if a Jacobian refactorize() fails it falls back to
    a full factorize() (reusing the existing symbolic factorization) before giving up,
    rather than reporting an error immediately.

    Use `get_linear_solver_stats()` on the solver to inspect how often factor/refactor
    calls happen and how often the fallback fires (see :class:`lightsim2grid.algorithm.LinearSolverStats`).

)mydelimiter";

const std::string DocSolver::NRSing_CKTSO = R"mydelimiter(
    This classes implements the Newton Raphson algorithm, the faster CKTSO solver available in the CKTSO library
    for the linear algebra. It does not support the distributed slack, but can be slightly faster than the :class:`lightsim2grid.algorithm.NR_CKTSO` .

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `NRSing_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.NRSing_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.NRSing_CKTSO)` at creation time

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso
 
)mydelimiter";

const std::string DocSolver::DC_CKTSO = R"mydelimiter(
    Alternative implementation of the DC solver, it uses the faster CKTSO solver available in the CKTSO library to solve for the DC voltage given the DC admitance matrix and
    the power injected at each nodes (requires a build from source).

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `DC_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.DC_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.DC_CKTSO)` at creation time

    .. warning::
        This is a DC solver that uses the DC approximation. If you want to use this approximation, you need to specified
        it when you create the grid2op environment, for example with "param.ENV_DC=True".

        Otherwise, it is used internally to find good starting point to intialize the real AC solver.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::FDPF_XB_CKTSO =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (XB version: "alg 2" / "fdxb"  in pypower / pandapower), it uses the fast CKTSO library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_XB_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_XB_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_XB_CKTSO)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for cktso.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::FDPF_BX_CKTSO =  R"mydelimiter(
    Default implementation of the Fast Decoupled Powerflow solver (BX version: "alg 3" / "fdbx"  in pypower / pandapower), it uses the fast CKTSO library for 
    its underlying sparse matrix manipulation.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `FDPF_BX_CKTSO` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.FDPF_BX_CKTSO)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.FDPF_BX_CKTSO)` at creation time    

    .. warning::
        
        Use this solver requires a compilation of lightsim2grid from source (see readme) AND an appropriate license for cktso.

    .. note::

        CKTSO is available at https://github.com/chenxm1986/cktso

)mydelimiter";

const std::string DocSolver::GaussSeidelAlgo = R"mydelimiter(
    Default implementation of the "Gauss Seidel" powerflow solver. We do not recommend to use it as the Newton Raphson based solvers
    are usually much (much) faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it is called `GaussSeidel` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.GaussSeidel)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.GaussSeidel)` at creation time

    .. warning::
        It currently does not support distributed slack.

)mydelimiter";

const std::string DocSolver::GaussSeidelSynchAlgo = R"mydelimiter(
    Variant implementation of the "Gauss Seidel" powerflow solver, where every buses are updated at once (can be significantly faster than the 
    :class:`lightsim2grid.algorithm.GaussSeidelAlgo` for larger grid). We still do not recommend to use it as the Newton Raphson based solvers
    are usually much (much) faster.

    See :ref:`available-powerflow-solvers` for more information on how to use it.

    .. note::

        In the enum :attr:`lightsim2grid.algorithm.AlgorithmType`, it called `GaussSeidelSynch` 
        
        You can use it with:
        
        - `env_lightsim.backend.set_algo_type(lightsim2grid.algorithm.GaussSeidelSynch)` after creation
        - `LightSimBackend(solver_type=lightsim2grid.algorithm.GaussSeidelSynch)` at creation time

    .. warning::
        It currently does not support distributed slack.
        
)mydelimiter";

const std::string DocSolver::AlgorithmSelector = R"mydelimiter(
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
    Retrieve the current solver used. This will return an instance of :class:`lightsim2grid.algorithm.AlgorithmType` indicating which
    is the underlying solver in use.

    This should be equivalent to :func:`lightsim2grid.network.LSGrid.get_algo_type()`

)mydelimiter";
const std::string DocSolver::chooseSolver_get_J_python = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.

    .. note::
        Depending on the underlying solver used (*eg* :class:`lightsim2grid.algorithm.DC_SparseLU` or :class:`lightsim2grid.algorithm.GaussSeidelAlgo`)
        the jacobian matrix might be irrelevant and an attempt to use this function will throw a RuntimeError. 

)mydelimiter";
const std::string DocSolver::get_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver when performing a powerflow.

    This is equivalent to the last (4th) element, ``timer_total_nr_``, of the tuple returned by
    ``***.get_timers()``.
)mydelimiter";

const std::string DocIterator::id = R"mydelimiter(
    Get the id of the element. Ids are integer from 0 to n-1 (if `n` denotes the number of such elements on the grid.)

    Examples
    --------
    We give the example only for generators, but it works similarly for every other types of objects
    in a :class:`lightsim2grid.network.LSGrid`.
    
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
    in a :class:`lightsim2grid.network.LSGrid`.
    
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

const std::string DocIterator::sub_id = R"mydelimiter(
    Get the substation id of the element.

    .. note::
        In pypowsybl, this is called "voltage levels".

    .. warning::
        Substation ids are optional and might not be set when reading the grid. In that case -1 is set for this attribute.

)mydelimiter";

const std::string DocIterator::pos_topo_vect = R"mydelimiter(
    Get the position of the element in the grid2op "topo_vect" vector.

    .. warning::
        Position in the "topo vector" are optional and might not be set when reading the grid. In that case -1 is set for this attribute.

)mydelimiter";

const std::string DocIterator::connected = R"mydelimiter(
    Get the status (True = connected, False = disconnected) of each element of a :class:`lightsim2grid.network.LSGrid`

)mydelimiter";

const std::string DocIterator::bus_id = R"mydelimiter(
    Get the bus id (as an integer) at which each element of a :class:`lightsim2grid.network.LSGrid` is connected. If `-1` is returned it means
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

                     i1                       ________             i2
         `bus 1` o------>   -----------------|r + j.x|---------<-------o `bus 2`
                 |       ) (            |                  |           |
                 |       ) (         |     |            |     |        |
                 | v1    ) ( n:1      | h1  |            | h2  |        | v2
                 |       ) (         |     |            |     |        |
                 \/      ) (            |                  |           \/
        ground---o-------   -------------------------------------------o---- ground

    (fyi: `i1`, `i2`, `n`, `h1` and `h2` are all complex numbers. `r` and `x` are real numbers. `j` is a complex number such that `j^2 = -1`)

    .. note::
        `h1` and `h2` are independent per-side shunt admittances, NOT necessarily one half of a single
        total value each (they can differ, eg for an asymmetric line/transformer coming from pypowsybl):
        the admittance matrix contribution of one branch is ``[[ys + h1, -ys], [-ys, ys + h2]]`` with
        ``ys = 1 / (r + j.x)`` (see :func:`lightsim2grid.elements.LineContainer.get_yac_eff_11` and
        friends for the coefficients actually used, including any tap-side / phase-shift correction for
        transformers).

    .. note::
        For a powerline, side 1 / side 2 used to be called `or` (origin) / `ex` (extremity) in older
        lightsim2grid versions; for a transformer they are `hv` (high voltage) / `lv` (low voltage)
        instead, since which physical side is tap-side matters there (see `is_tap_side_1`).

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
    Retrieve the shunt admittance (in pair unit system) of one side of the powerline / transformer:
    conductance `g` as the real part, susceptance `b` (related to the line charging capacitance) as the
    imaginary part, ie ``h = g + 1j * b``.

    This is a complex number, represented by `h1` (side 1) or `h2` (side 2) in the line model -- they are
    independent values, not half of a single shared `h` each (see the note in the line model below).

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::nb_max_busbars = R"mydelimiter(
    Maximum number of busbars (independent buses) allowed at this substation (``int``, > 0).

    This is the per-substation value of what :func:`lightsim2grid.network.LSGrid.set_max_nb_bus_per_sub`
    sets grid-wide: the substation has exactly this many candidate buses, some of which may be
    unused (disconnected) at any given time.

)mydelimiter";

const std::string DocIterator::vn_kv = R"mydelimiter(
    Nominal voltage of this substation, in kV (``float``).

)mydelimiter";


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
    This class allows to iterate through the generators of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
        to computed the powerflow do not support distributed slack buses - **eg** :class:`lightsim2grid.algorithm.NRSing_SparseLU`)

        This is why we recommend to use the (slower) but more accurate :class:`lightsim2grid.algorithm.NR_SparseLU` or 
        :class:`lightsim2grid.algorithm.NR_KLU` for example.

)mydelimiter";

const std::string DocIterator::slack_weight = R"mydelimiter(
    For each generators, gives the participation (for the distributed slack) of this particular generator.
    
    .. note::
        Weights do not scale to one for this variable thus this number has no meaning by itself and should be compared with the others.

)mydelimiter";

const std::string DocIterator::min_q_mvar = R"mydelimiter(
    Minimum reactive value that can be produced / absorbed by this generator, in MVAr.

    .. note::
        On a :class:`lightsim2grid.elements.GenInfo` or :class:`lightsim2grid.elements.ConverterStationInfo`
        that is locally voltage-regulating (``voltage_regulator_on`` is ``True`` and it does not regulate a
        remote bus), this is genuinely used at every solve: when several such units share the same bus, their
        reactive-power mismatch is split between them proportionally to ``max_q_mvar - min_q_mvar``. It is
        also used, in the same case, by :func:`lightsim2grid.gridmodel.LSGrid.check_solution` when
        ``check_q_limits`` is ``True``, to report any part of the mismatch that falls outside
        ``[min_q_mvar, max_q_mvar]`` instead of masking it.

        On a "PQ" generator (``voltage_regulator_on`` is ``False``), a remotely-regulating one, or a
        :class:`lightsim2grid.elements.SGenInfo` (static generators never regulate voltage), this value is
        NOT used anywhere by lightsim2grid: it is pure metadata carried over from the source model.

)mydelimiter";

const std::string DocIterator::max_q_mvar = R"mydelimiter(
    Maximum reactive value that can be produced / absorbed by this generator, in MVAr. See `min_q_mvar` for
    when (and how) this is actually used.

)mydelimiter";

const std::string DocIterator::min_p_mw = R"mydelimiter(
    Minimum active value that can be produced / absorbed by this static generator, in MW.

    .. note::
        This is NOT used anywhere by lightsim2grid today: it is not enforced by the solver, and
        :func:`lightsim2grid.gridmodel.LSGrid.check_solution` does not examine static generators at all
        (only :class:`lightsim2grid.elements.GenInfo` / :class:`lightsim2grid.elements.ConverterStationInfo`,
        see `min_q_mvar`). It is pure metadata carried over from the source model.

)mydelimiter";

const std::string DocIterator::max_p_mw = R"mydelimiter(
    Maximum active value that can be produced / absorbed by this static generator, in MW. See `min_p_mw`.

)mydelimiter";

const std::string DocIterator::SGenContainer = R"mydelimiter(
    This class allows to iterate through the static generators of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
    This class allows to iterate through the loads **and storage units** of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
    This class allows to iterate through the load of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
    This class allows to iterate through the transformers of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
    This class allows to iterate through the powerlines of the :class:`lightsim2grid.network.LSGrid` easily, as if they were
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
            line.bus1_id

        print(f"There are {len(grid_model.get_lines())} lines on the grid.")

        first_line = grid_model.get_lines()[0]

    You can have a look at :class:`lightsim2grid.elements.LineInfo` for properties of these elements.

)mydelimiter";

const std::string DocIterator::LineInfo = R"mydelimiter(
    This class represents what you get from retrieving the powerlines from 
    :class:`lightsim2grid.elements.LineContainer`.

    It allows to read information from each powerline of the powergrid.

    Powerlines have two sides, "1" and "2" (called "or" for "origin" and "ex" for "extremity" in older
    lightsim2grid versions), that are connected and linked to each other by some equations.

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
            line.bus1_id

    Notes
    -----
    Line are modeled using the "line model" as shown in the schema at the end of the paragraph.

    The tap ratio `n` on this schema will be 1.0 for all powerline. If you want to model phase shifters, please
    model them as Trafo (see :class:`lightsim2grid.elements.TrafoInfo`)

    For more information about the model and the equations linking all the quantities, please visit 
    `matpower manual <https://matpower.org/docs/MATPOWER-manual.pdf>`_ , especially the "3. Modeling" and the
    "3.2 Branches" subsection, as well as the equation 3.1, 3.2 and 3.3 therein.
    
)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::bus_1_id = R"mydelimiter(
    Get the bus id (as an integer) at which side 1 of the line is connected. If `-1` is returned it means
    that the line is disconnected.

)mydelimiter";

const std::string DocIterator::bus_2_id = R"mydelimiter(
    Get the bus id (as an integer) at which side 2 of the line is connected. If `-1` is returned it means
    that the line is disconnected.

)mydelimiter";

const std::string DocIterator::res_p_1_mw = R"mydelimiter(
    Get the active power in MW at side 1 of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_p_2_mw = R"mydelimiter(
    Get the active power in MW at side 2 of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_1_mvar = R"mydelimiter(
    Get the reactive power in MVAr at side 1 of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_2_mvar = R"mydelimiter(
    Get the reactive power in MVAr at side 2 of the line. If it is positive it means power is absorbed by the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_1_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which side 1 of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_2_deg = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which side 2 of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_1_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which side 1 of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_2_kv = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which side 2 of the line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_1_ka = R"mydelimiter(
    Get the current flows (in kA) at side 1 of the line.

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_a_2_ka = R"mydelimiter(
    Get the current flows (in kA) at side 2 of the line.

)mydelimiter" + DocIterator::only_avail_res;


const std::string DocIterator::DCLineContainer = R"mydelimiter(
    This class allows to iterate through the hvdc lines of the :class:`lightsim2grid.network.LSGrid` easily,
    as if they were in a python list. (Kept under the historical name `DCLineContainer` / `get_dclines()` for
    backward compatibility; the legacy pandapower dc line is now just a special case of the model below.)

    The model follows powsybl IIDM / open-loadflow: the hvdc line itself owns the active power
    (`p_setpoint_mw`, drawn at the rectifier, and `converters_mode`, saying which side rectifies), while its
    two embedded converter stations (`station1` / `station2`, one on each side of the line, see
    :class:`lightsim2grid.elements.ConverterStationInfo`) own the reactive power / voltage behaviour and can
    be controlled independently for their voltage setpoint. See :class:`lightsim2grid.elements.HvdcLineInfo`
    for the loss model turning the setpoint at the rectifier side into the active power actually injected at
    the other side, and for the angle-droop ("AC emulation") alternative to a fixed setpoint.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # manipulate the hvdc lines (usually there are none...)
        for hvdc_line in grid_model.get_dclines():
            # do something with the line !
            hvdc_line.bus1_id

        print(f"There are {len(grid_model.get_dclines())} hvdc lines on the grid.")

    You can have a look at :class:`lightsim2grid.elements.HvdcLineInfo` for properties of these elements.

)mydelimiter";


const std::string DocIterator::DCLineInfo = R"mydelimiter(
    This class represents what you get from retrieving the hvdc lines from
    :class:`lightsim2grid.elements.HvdcLineContainer`.

    It allows to read information from each hvdc line of the powergrid.

    Hvdc lines have two sides, "1" and "2", each with its own converter station (`station1` / `station2`,
    see :class:`lightsim2grid.elements.ConverterStationInfo`) -- the equivalent of the "origin" / "extremity"
    naming used in older lightsim2grid versions for AC powerlines and transformers.

    For accessing the results, it's basically the same as having two "elements" (so you get two "voltage
    magnitude" `res_v1_kv` / `res_v2_kv`, two "injected power" `res_p1_mw` / `res_p2_mw`, etc.)

    .. warning::
        Data can only be read from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        # create a lightsim2grid "gridmodel"
        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # for hvdc lines
        first_hvdc_line = grid_model.get_dclines()[0]  # first hvdc line, this is an `HvdcLineInfo`
        for hvdc_line in grid_model.get_dclines():
            # hvdc_line is an `HvdcLineInfo`
            hvdc_line.bus1_id

    Notes
    -----
    See :func:`lightsim2grid.elements.HvdcLineInfo.target_p1_mw` for the active-power loss model turning the
    setpoint given at the rectifier side into the active power actually injected at the other side, and
    `droop_enabled` / `droop_p0_mw` / `droop_k_mw_per_rad` for the angle-droop ("AC emulation") alternative,
    where the active power follows the angle difference between the two sides instead of a fixed setpoint.

)mydelimiter";

const std::string DocIterator::target_p_1_mw_dcline = R"mydelimiter(
    The active power target (in MW) of the converter station on side 1 of the hvdc line, generator sign
    convention (positive = power injected into the AC grid at side 1).

    For a line NOT in angle-droop mode, this is derived from `p_setpoint_mw` / `converters_mode` through the
    loss model described in :class:`lightsim2grid.elements.HvdcLineInfo` -- it is the target for the AC
    powerflow, not necessarily equal to `p_setpoint_mw` itself (which is always >= 0 and lives at the
    rectifier side, whichever side that is). See `target_p2_mw` for the side-2 counterpart.

    .. note::
        In angle-droop mode (`droop_enabled` is True), the active power actually used by the solver instead
        follows `p0 + k * (theta1 - theta2)` and is NOT read from this field; see `droop_p0_mw` /
        `droop_k_mw_per_rad`.

)mydelimiter";

const std::string DocIterator::target_p_2_mw_dcline = R"mydelimiter(
    The active power target (in MW) of the converter station on side 2 of the hvdc line, generator sign
    convention (positive = power injected into the AC grid at side 2). See `target_p1_mw` for the side-1
    counterpart and the loss model turning one into the other.

)mydelimiter";

const std::string DocIterator::target_vm_1_pu_dcline = R"mydelimiter(
    The target voltage setpoint (in pu, NOT in kV) of the converter station on side 1 of the hvdc line.

)mydelimiter";

const std::string DocIterator::target_vm_2_pu_dcline = R"mydelimiter(
    The target voltage setpoint (in pu, NOT in kV) of the converter station on side 2 of the hvdc line.

)mydelimiter";

const std::string DocIterator::dc_line_formula = R"mydelimiter(

    .. note::
        The active power actually injected at one side of an hvdc line is derived from the active power
        setpoint at the rectifier side through a loss model (mirrors open-loadflow's
        `HvdcUtils.getConverterStationTargetP`, extended with the legacy pandapower fixed-loss term)::

            line_in   = (1 - lf_rect) * (1 - loss_pct / 100) * p_setpoint_mw
            line_loss = r_ohm * line_in^2 / nominal_v_kv^2        (0 when nominal_v_kv == 0)
            received  = (1 - lf_inv) * (line_in - line_loss) - loss_mw

        where `p_setpoint_mw` (>= 0) is drawn at the rectifier side (`converters_mode` says which side that
        is), `lf_rect` / `lf_inv` are the rectifier / inverter converter stations' own loss factors, and
        `received` is the target active power (generator convention) at the non-rectifier side. The legacy
        pandapower dc line maps onto this exactly with station loss factors = 0 and `r_ohm` = 0.

    .. note::
        Both `target_p1_mw` and `target_p2_mw` use the **generator** sign convention: a positive value means
        power is injected into the AC grid at that side (so a positive `target_p1_mw` means power flows from
        side 2 to side 1 through the line).

    .. note::
        In angle-droop mode (`droop_enabled` is True), none of the above applies: the active power instead
        follows `p0 + k * (theta1 - theta2)`; see `droop_p0_mw` / `droop_k_mw_per_rad` / `status_droop`.

)mydelimiter";

const std::string DocIterator::loss_pct = R"mydelimiter(
    The `loss_pct` (relative loss, in percent) parameter of the hvdc line, used in the active-power loss
    model below.

)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocIterator::loss_mw = R"mydelimiter(
    The `loss_mw` (flat loss, in MW) parameter of the hvdc line, used in the active-power loss model below.

)mydelimiter" + DocIterator::dc_line_formula;

const std::string DocIterator::res_p_1_mw_dcline = R"mydelimiter(
    The active power actually injected at side 1 of the hvdc line (in MW, generator convention).

)mydelimiter" + DocIterator::only_avail_res + DocIterator::dc_line_formula;

const std::string DocIterator::res_p_2_mw_dcline = R"mydelimiter(
    The active power actually injected at side 2 of the hvdc line (in MW, generator convention).

)mydelimiter" + DocIterator::only_avail_res + DocIterator::dc_line_formula;

const std::string DocIterator::res_q_1_mvar_dcline = R"mydelimiter(
    The reactive power actually injected at side 1 of the hvdc line (in MVAr, generator convention).

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_q_2_mvar_dcline = R"mydelimiter(
    The reactive power actually injected at side 2 of the hvdc line (in MVAr, generator convention).

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_1_kv_dcline = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which side 1 of the hvdc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_v_2_kv_dcline = R"mydelimiter(
    Get the magnitude of the complex voltage (in kV) of the bus at which side 2 of the hvdc line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_v_kv"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_1_deg_dcline = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which side 1 of the hvdc
    line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;

const std::string DocIterator::res_theta_2_deg_dcline = R"mydelimiter(
    Get the angle of the complex voltage (in degree, not in radian) of the bus at which side 2 of the hvdc
    line is connected.

    .. note::
        All elements (load, generators, side of powerline etc.) connected at the same bus have the same "res_theta_deg"

)mydelimiter" + DocIterator::only_avail_res;


const std::string DocLSGrid::LSGrid = R"mydelimiter(
    This class represent a lightsim2grid power network. All the elements that can be manipulated by
    lightsim2grid are represented here.

    We do not recommend to use this class directly, but rather to use a :class:`lightsim2grid.lightSimBackend.LightSimBackend`.

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

const std::string DocLSGrid::check_grid = R"mydelimiter(
    Check that the grid is internally consistent and safe to run a powerflow on.

    It verifies that every index the grid carries is in range: the bus id of each
    element (load, generator, static generator, storage, shunt, line, transformer,
    hvdc line, static var compensator), the substation id and the position in the
    topology vector (both optional), and the generator slack / remote-regulated bus
    references.

    This is called automatically when a grid is loaded (from a pickle or from the
    fast binary format), and by the grid loaders (from pandapower, pypowsybl,
    matpower or powermodels). You normally do not need to call it yourself; it is
    exposed so you can validate a grid you built or modified by hand.

    Raises
    ------
    An ``IndexError`` (C++ ``std::out_of_range``) if an index is out of range, or a
    ``RuntimeError`` (C++ ``std::runtime_error``) on a structural inconsistency.
    Returns ``None`` if the grid is consistent.

    Notes
    -----
    Runs in time proportional to the number of elements in the grid, so it is cheap
    compared to a powerflow.

)mydelimiter";

const std::string DocLSGrid::change_algorithm =  R"mydelimiter(
    This function allows to control which solver is used during the powerflow. See the section :ref:`available-powerflow-solvers` for
    more information about them.

    .. seealso:: :attr:`lightsim2grid.algorithm.AlgorithmType` for a list of the available algorithms (NB: some algorithms might not be available on all platform)

    .. note::
        If the algorithm type entered is a `DC` algorithm (**eg** from :attr:`lightsim2grid.algorithm.AlgorithmType`,
        `DC_SparseLU`, `DC_KLU` or `DC_NICSLU`), it will change the `_dc_solver` otherwise the regular `_solver`
        is modified.

    .. rubric:: Examples

    .. code-block:: python

        from lightsim2grid.algorithm import AlgorithmType
        # init the grid model
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

        # change the algorithm used for the powerflow
        # to use internally a Newton Raphson algorithm with the Eigen sparse LU linear solver
        lightsim_grid_model.change_algorithm(AlgorithmType.NR_SparseLU)

)mydelimiter";

const std::string DocLSGrid::available_default_algorithms =  R"mydelimiter(
    Return the list of the names of the algorithm available on the current lightsim2grid installation.

    This is a list of :attr:`lightsim2grid.algorithm.AlgorithmType`.

)mydelimiter";


const std::string DocLSGrid::available_algorithm_names =  R"mydelimiter(
    Return the list of the names of the algorithm available on the current lightsim2grid installation.

    This is a list of string.

)mydelimiter";

const std::string DocLSGrid::get_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver when performing a powerflow.

    This is equivalent to the `get_computation_time` of the :func:`lightsim2grid.algorithm.AlgorithmSelector.get_computation_time` of
    the solver used (:func:`lightsim2grid.network.LSGrid.get_solver`)
    
)mydelimiter";
const std::string DocLSGrid::get_dc_computation_time = R"mydelimiter(
    Return the total computation time (in second) spend in the solver (used to perform DC approximation) when performing a DC powerflow.

    This is equivalent to the `get_computation_time` of the :func:`lightsim2grid.algorithm.AlgorithmSelector.get_computation_time` of
    the DC solver used (:func:`lightsim2grid.network.LSGrid.get_dc_solver`)
    
)mydelimiter";
const std::string DocLSGrid::get_algo_type = R"mydelimiter(
    Return the type of the solver currently used.

    This is equivalent to the `get_type` of the :func:`lightsim2grid.algorithm.AlgorithmSelector.get_type` of
    the solver used.

)mydelimiter";
const std::string DocLSGrid::get_dc_algo_type = R"mydelimiter(
    Return the type of the solver currently used to compute DC powerflow.

)mydelimiter";
const std::string DocLSGrid::get_algo = R"mydelimiter(
    Return the solver currently in use as a :func:`lightsim2grid.algorithm.AlgorithmSelector` instance.

)mydelimiter";
const std::string DocLSGrid::get_dc_algo = R"mydelimiter(
    Return the solver currently in use as a :func:`lightsim2grid.algorithm.AlgorithmSelector` instance for the dc powerflow.

)mydelimiter";

const std::string DocLSGrid::get_lines = R"mydelimiter(
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
const std::string DocLSGrid::get_trafos = R"mydelimiter(
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
const std::string DocLSGrid::get_generators = R"mydelimiter(
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
const std::string DocLSGrid::get_static_generators = R"mydelimiter(
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
const std::string DocLSGrid::get_shunts = R"mydelimiter(
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
const std::string DocLSGrid::get_storages = R"mydelimiter(
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
const std::string DocLSGrid::get_loads = R"mydelimiter(
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

const std::string DocLSGrid::get_dclines = R"mydelimiter(
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

const std::string DocLSGrid::_internal_do_not_use = R"mydelimiter(
        INTERNAL

        .. warning:: /!\\ Internal, do not use unless you know what you are doing /!\\

        This is used as part of a dedicated code for :class:`lightsim2grid.lightSimBackend.LightSimBackend`

)mydelimiter";

const std::string DocLSGrid::J_description = R"mydelimiter(
    J is NOT a fixed-shape 2x2 block anymore: it is assembled by composing a common "Base" block
    with independent extensions, and which extensions are active depends on the solver
    (:class:`lightsim2grid.algorithm.NRSing_SparseLU` and friends use only ``Base`` + ``VoltageControl``
    + ``Hvdc``; :class:`lightsim2grid.algorithm.NR_SparseLU` and friends additionally use
    ``MultiSlack``). Each part below claims its own rows (equations) / columns (unknowns); nothing
    below is claimed twice.

    **Base** (always present) is the usual decoupled-looking core::

        | J11 | J12 | = dimensions: | (pvpq, pvpq) | (pvpq, pq) |
        | --------- |               | ------------------------ |
        | J21 | J22 |               |  (pq, pvpq)  |  (pq, pq)  |

    with:

    - `J11` = dS_dVa[array([pvpq]).T, pvpq].real (= real part of dS / dVa for all pv and pq buses)
    - `J12` = dS_dVm[array([pvpq]).T, pq].real
    - `J21` = dS_dVa[array([pq]).T, pvpq].imag
    - `J22` = dS_dVm[array([pq]).T, pq].imag (= imaginary part of dS / dVm for all pq buses)

    .. note::
        A slack bus that is NOT locally pinned by its own directly-connected voltage-regulating
        generator (a PQ distributed-slack participant, or a slack bus regulated only remotely / by
        an SVC -- see the ``VoltageControl`` extension below) also gets a free vm unknown and a Q
        equation added by ``Base``, exactly like an ordinary pq bus. This is NOT restricted to "all
        but the first ref bus": it depends on how each individual slack bus is actually voltage-pinned.

    **MultiSlack** (distributed-slack solvers only, ie ``NR_*``, not ``NRSing_*``): for every slack
    bus it adds one P-equation row (including the reference bus'); for every slack bus OTHER than
    the reference, it additionally adds a theta unknown column. On top of that, it adds exactly ONE
    extra column, shared by the whole system, for the "slack_absorbed" unknown (the distributed
    slack's total absorbed mismatch) -- not one extra row/column pair per slack bus.

    **VoltageControl** (always present; covers both a generator remotely regulating another bus'
    voltage and an SVC in voltage-control mode): controllers sharing the same regulated bus are
    grouped; each group of N controllers adds N reactive-power (Q) unknown columns (one per
    controller -- a plain, non-regulating "PQ" generator gets none), 1 voltage-setpoint row, and
    N-1 reactive-power-sharing rows.

    **Hvdc** (always present, angle-droop / "AC emulation" hvdc lines only): claims no row or column
    of its own; it only adds extra dP/dtheta terms into the P-mismatch rows / theta columns that
    ``Base``/``MultiSlack`` already registered for the two buses of each droop-controlled hvdc line.

    .. note::
        the notation `pvpq` above means "the concatenation of the pv vector and the pq vector".
        Slack buses (participating in the distributed slack or not) are NOT part of `pvpq`: they are
        registered by ``Base``/``MultiSlack`` as described above, independently of it.

    .. note::
        All notation here are notation for the solver. You should use `gridmodel.get_pq_solver()` and
        `gridmodel.get_pv_solver()` to retrieve their  value.

)mydelimiter";

const std::string DocLSGrid::get_J_python = R"mydelimiter(
    The jacobian matrix is an internal object to the solver and should only be used
    when on knows how exactly it is filled.
)mydelimiter";

const std::string DocLSGrid::get_Va = R"mydelimiter(
    Returns the voltage angles for each buses as a numpy vector of real number.
    This vector have the size of the total number buses on the system, including the 
    disconnected bus. It adopts the "gridmodel" labelling.

    .. versionchanged:: 0.9.0
        They are labelled with the `grimodel` labelling. To retrieve the
        previous behaviour (solver labelling) you can use the current
        :func:`lightsim2grid.network.LSGrid.get_Va_solver` (before version 0.9.0)

    .. danger:: 
        Some breaking change have been introduced in lighsim2grid 0.9.0.
        You can :func:`lightsim2grid.network.LSGrid.get_Va_solver` to get the
        previous (before 0.9.0) behaviour.

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` 
        (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.

)mydelimiter";

const std::string DocLSGrid::get_Vm = R"mydelimiter(
    Returns the voltage magnitude for each buses as a numpy vector of real number. 
    This vector have the size of the total number buses on the system, including the 
    disconnected bus. It adopts the "gridmodel" labelling.

    .. versionchanged:: 0.9.0
        They are labelled with the `grimodel` labelling. To retrieve the
        previous behaviour (solver labelling) you can use the current
        :func:`lightsim2grid.network.LSGrid.get_Vm_solver` (before version 0.9.0)

    .. danger:: 
        Some breaking change have been introduced in lighsim2grid 0.9.0.
        You can :func:`lightsim2grid.network.LSGrid.get_Vm_solver` to get the previous 
        (before 0.9.0) behaviour.

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` 
        (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.
)mydelimiter";

const std::string DocLSGrid::get_V = R"mydelimiter(
    Returns the complex voltage for each buses as a numpy vector of complex number. 
    This vector have the size of the total number buses on the system, including the 
    disconnected bus. It adopts the "gridmodel" labelling.

    .. versionchanged:: 0.9.0
        They are labelled with the `grimodel` labelling. To retrieve the
        previous behaviour (solver labelling) you can use the current
        :func:`lightsim2grid.network.LSGrid.get_V_solver` (before version 0.9.0)

    .. danger:: 
        Some breaking change have been introduced in lighsim2grid 0.9.0.
        You can :func:`lightsim2grid.network.LSGrid.get_V_solver` to get the previous 
        (before 0.9.0) behaviour.

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.
)mydelimiter";

const std::string DocLSGrid::get_J_python_solver = R"mydelimiter(
    Returns the Jacobian matrix used for solving the powerflow as a scipy sparse CSC matrix matrix of real number.

    The "jacobian" matrix is only available for some powerflow algorithms
    (the one based on the Newton Raphson algorithm) and we provide it only for the last computed iteration.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous (before
        version 0.9.0) behaviour of this function, which used to be called ``get_J``.

    .. versionadded:: 0.9.0
        This function is the renamed ``get_J`` of earlier lightsim2grid versions. Unlike
        :func:`lightsim2grid.network.LSGrid.get_Va`/:func:`lightsim2grid.network.LSGrid.get_Vm`,
        no gridmodel-labelled ``get_J`` was added: the Jacobian is only ever exposed with the
        solver labelling, through this function.

    .. note::
        Some powerflows (*eg* DC or Gauss Seidel) do not rely on jacobian matrix, in this case, 
        calling this function will return an exception. 

)mydelimiter" + DocLSGrid::J_description;

const std::string DocLSGrid::get_Va_solver = R"mydelimiter(
    Returns the voltage angles for each buses as a numpy vector of real number. 
    This vector have the size of the number of active buses on the system and
    adopts the "solver" labelling.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_Va` (before version 0.9.0)

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_Va` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_Va`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` 
        (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.

)mydelimiter";

const std::string DocLSGrid::get_Vm_solver = R"mydelimiter(
    Returns the voltage magnitude for each buses as a numpy vector of real number. 
    This vector have the size of the number of active buses on the system and
    adopts the "solver" labelling.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_Vm` (before version 0.9.0)

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_Vm` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_Vm`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.

)mydelimiter";

const std::string DocLSGrid::get_V_solver = R"mydelimiter(
    Returns the complex voltage for each buses as a numpy vector of complex number. 
    This vector have the size of the number of active buses on the system and
    adopts the "solver" labelling.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_V` (before version 0.9.0)

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_V` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_V`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. note::
        You can use the :attr:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` (or :attr:`lightsim2grid.network.LSGrid.id_dc_solver_to_me`) to know at which bus
        (on the grid) they corresponds.

)mydelimiter";


const std::string DocLSGrid::id_me_to_ac_solver = R"mydelimiter(
    In lightsim2grid, buses are labelled from `0` to `n-1` (if `n` denotes the total number of buses on the grid) [this is called "**grid model bus id**"]

    At any given point in time, some buses might be deactivated (for example because nothing is connected to them).

    On the other end, the solvers need a contiguous list of only active buses (otherwise they might run into divergence issue) [this will be called 
    "**solver bus id**" later on]

    This function allows, for all buses of the :class:`lightsim2grid.network.LSGrid` to know on which "solver bus" they are affected. It
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
    
    .. seealso:: :class:`lightsim2grid.network.LSGrid.id_me_to_dc_solver` for its counterpart when a dc powerflow is used
    
    .. seealso:: :class:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for the "reverse" operation (given a "solver bus" id, returns 
        the "gridmodel bus id")

    Notes
    -----

    For all steps, you have the propertie that, if `id_ac_solver_to_me = gridmodel.id_ac_solver_to_me()` and `id_me_to_ac_solver = gridmodel.id_me_to_ac_solver()`
    and by denoting `gridmodel_bus_id = np.arange(gridmodel.total_bus())` and `solver_bus_id = np.arange(gridmodel.nb_connected_bus())`:

    - `solver_bus_id` and `id_ac_solver_to_me` have the same shape
    - `gridmodel_bus_id` and `id_me_to_ac_solver` have the same shape
    - `solver_bus_id` is shorter (or of the same length) than `gridmodel_bus_id`
    - the connected bus (in the grid model) are given by `gridmodel_bus_id[id_ac_solver_to_me]`, and it gives their order

)mydelimiter";

const std::string DocLSGrid::id_ac_solver_to_me = R"mydelimiter(
    In lightsim2grid, buses are labelled from `0` to `n-1` (if `n` denotes the total number of buses on the grid) [this is called "**grid model bus id**"]

    At any given point in time, some buses might be deactivated (for example because nothing is connected to them).

    On the other end, the solvers need a contiguous list of only active buses (otherwise they might run into divergence issue) [this will be called 
    "**solver bus id**" later on]

    This function allows, for all buses exported in the solver, to retrieve which was the initial bus in the :class:`lightsim2grid.network.LSGrid`. It
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
    
    .. seealso:: :class:`lightsim2grid.network.LSGrid.id_dc_solver_to_me` for its counterpart when a dc powerflow is used
    
    .. seealso:: :class:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` for the "reverse" operation (given a "solver bus" id, returns 
        the "gridmodel bus id")

    Notes
    -----

    For all steps, you have the propertie that, if `id_ac_solver_to_me = gridmodel.id_ac_solver_to_me()` and `id_me_to_ac_solver = gridmodel.id_me_to_ac_solver()`
    and by denoting `gridmodel_bus_id = np.arange(gridmodel.total_bus())` and `solver_bus_id = np.arange(gridmodel.nb_connected_bus())`:

    - `solver_bus_id` and `id_ac_solver_to_me` have the same shape
    - `gridmodel_bus_id` and `id_me_to_ac_solver` have the same shape
    - `solver_bus_id` is shorter (or of the same length) than `gridmodel_bus_id`
    - the connected bus (in the grid model) are given by `gridmodel_bus_id[id_ac_solver_to_me]`, and it gives their order

)mydelimiter";

const std::string DocLSGrid::id_me_to_dc_solver = R"mydelimiter(
    Same as :class:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` but only used for the DC approximation.
)mydelimiter";

const std::string DocLSGrid::id_dc_solver_to_me = R"mydelimiter(
    Same as :class:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` but only used for the DC approximation.
)mydelimiter";

const std::string DocLSGrid::total_bus = R"mydelimiter(
    Returns (>0 integer) the total number of buses in the powergrid (both connected and disconnected)
)mydelimiter";

const std::string DocLSGrid::nb_connected_bus = R"mydelimiter(
    Returns (>0 integer) the number of connected buses on the powergrid (ignores the disconnected bus).
)mydelimiter";

const std::string DocLSGrid::get_pv = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PV" (ie the buses on which at least a generator is connected.).

    It returns a vector of integer. 

    .. danger::
        From lightsim2grid 0.9.0 they are labelled with the `gridmodel` labelling.

        This behaviour is now accessible with the
        :func:`lightsim2grid.network.LSGrid.get_pv` before version 0.9.0

    .. versionchanged:: 0.9.0
        The new version of this function returns the id labelled with the gridmodel convention (for consistency).

        Earlier version returned the labelling in the "solver" convention. To access the earlier function, please 
        use the :func:`lightsim2grid.network.LSGrid.get_pv` function.

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
        
)mydelimiter";

const std::string DocLSGrid::get_pq = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PQ".

    It returns a vector of integer.

    .. danger::
        From lightsim2grid 0.9.0 they are labelled with the `gridmodel` labelling.

        This behaviour is now accessible with the
        :func:`lightsim2grid.network.LSGrid.get_pq` before version 0.9.0

    .. versionchanged:: 0.9.0
        The new version of this function returns the id labelled with the gridmodel convention (for consistency).

        Earlier version returned the labelling in the "solver" convention. To access the earlier function, please 
        use the :func:`lightsim2grid.network.LSGrid.get_pq` function.

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it will might be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_ids = R"mydelimiter(
    Returns the ids of the buses that are part of the distributed slack.

    It returns a vector of integer.

    .. danger::
        From lightsim2grid 0.9.0 they are labelled with the `gridmodel` labelling.

        This behaviour is now accessible with the
        :func:`lightsim2grid.network.LSGrid.get_slack_ids_solver` before version 0.9.0

    .. versionchanged:: 0.9.0
        The new version of this function returns the id labelled with the gridmodel convention (for consistency).

        Earlier version returned the labelling in the "solver" convention. To access the earlier function, please 
        use the :func:`lightsim2grid.network.LSGrid.get_slack_ids_solver` function.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_ids_dc = R"mydelimiter(
    Returns the ids of the buses that are part of the distributed slack. For DC, the active-power mismatch is spread across these buses proportionally to their `slack_weights` (see :func:`lightsim2grid.network.LSGrid.dc_pf`) -- distributed slack IS taken into account for the DC powerflow itself; it is only `get_ptdf` / `get_lodf` that still assume a single slack bus.

    It returns a vector of integer.

    .. danger::
        From lightsim2grid 0.9.0 they are labelled with the `gridmodel` labelling.

        This behaviour is now accessible with the
        :func:`lightsim2grid.network.LSGrid.get_slack_ids_dc_solver` before version 0.9.0

    .. versionchanged:: 0.9.0
        The new version of this function returns the id labelled with the gridmodel convention (for consistency).

        Earlier version returned the labelling in the "solver" convention. To access the earlier function, please 
        use the :func:`lightsim2grid.network.LSGrid.get_slack_ids_dc_solver` function.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_weights = R"mydelimiter(
    For each bus in the gridmodel solver, it outputs its participation to the distributed slack.

    It's 0 if the current bus does not participate to it, otherwise it is made of > 0. real numbers.

    This vector sums to 1 and has the same size as the number of active buses on the grid.

    .. danger::
        From lightsim2grid 0.9.0 they are labelled with the `gridmodel` labelling.

        This behaviour is now accessible with the
        :func:`lightsim2grid.network.LSGrid.get_slack_weights_solver` before version 0.9.0

    .. versionchanged:: 0.9.0
        The new version of this function returns the id labelled with the gridmodel convention (for consistency).

        Earlier version returned the labelling in the "solver" convention. To access the earlier function, please 
        use the :func:`lightsim2grid.network.LSGrid.get_slack_weights_solver` function.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_pv_solver = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PV" (ie the buses on which at least a generator is connected.).

    It returns a vector of integer. 

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_pv` before version 0.9.0

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_pv` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_pv`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
        
)mydelimiter";

const std::string DocLSGrid::get_pq_solver = R"mydelimiter(
    Returns the ids of the buses that are labelled as "PQ".

    It returns a vector of integer.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_pq` before version 0.9.0

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_pq` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_pq`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it will might be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_ids_solver = R"mydelimiter(
    Returns the ids of the buses that are part of the distributed slack.

    It returns a vector of integer.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_slack_ids` before version 0.9.0

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_slack_ids` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_slack_ids`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_ids_dc_solver = R"mydelimiter(
    Returns the ids of the buses that are part of the distributed slack. For DC, the active-power mismatch is spread across these buses proportionally to their `slack_weights` (see :func:`lightsim2grid.network.LSGrid.dc_pf`) -- distributed slack IS taken into account for the DC powerflow itself; it is only `get_ptdf` / `get_lodf` that still assume a single slack bus.

    It returns a vector of integer.

    .. versionadded:: 0.9.0
        Only what is now :func:`lightsim2grid.network.LSGrid.get_slack_ids_solver` (that used 
        to be called :func:`lightsim2grid.network.LSGrid.get_slack_ids`) was available. 

        There were no possibility to retrieve that for DC powerflow.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_slack_ids` before version 0.9.0

    .. warning:: 
        The index are given in the "solver bus" convention. This means that it might not be the bus of the original grid model.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_slack_weights_solver = R"mydelimiter(
    For each bus used by the solver, it outputs its participation to the distributed slack.

    It's 0 if the current bus does not participate to it, otherwise it is made of > 0. real numbers.

    This vector sums to 1 and has the same size as the number of active buses on the grid.

    .. danger::
        They are labelled with the `solver` labelling, which corresponds to the previous behaviour in 
        :func:`lightsim2grid.network.LSGrid.get_slack_weights` before version 0.9.0

    .. versionadded:: 0.9.0
        This function replace the :func:`lightsim2grid.network.LSGrid.get_slack_weights_solver` of earlier
        lightsim2grid version. The new version of :func:`lightsim2grid.network.LSGrid.get_slack_weights_solver`
        now returns the id labelled with the gridmodel convention (for consistency).

    .. warning:: 
        This vector represents "solver buses" and not "original grid model buses".

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

)mydelimiter";

const std::string DocLSGrid::get_Ybus_solver = R"mydelimiter(
    This function returns the (complex) `Ybus` matrix used to compute the AC powerflow.

    The resulting matrix is a CSC scipy sparse matrix of complex number.

    It is a square matrix, as many rows (columns) as there are **connected** buses on the grid.

    .. seealso::
        If you want to retrieve the Ybus adopting the "gridmodel" bus labelling, you can use 
        :func:`lightsim2grid.gridmodel.get_Ybus`

    .. versionadded:: 0.9.0
        It was named `get_Ybus` before this version, but the name has been changed to avoid confusing AND a new
        function (this one) has been made with the proper `gridmodel` labelling.

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system !

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Ybus[id_me_to_ac_solver[k],:]` (rows of this bus), `Ybus[:, id_me_to_ac_solver[k]]` (column for this bus) 

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter";

const std::string DocLSGrid::get_dcYbus_solver = R"mydelimiter(
    This function returns the (complex) `Ybus` matrix used to compute the DC powerflow
    (its imaginary part should be 0.).

    The resulting matrix is a CSC scipy sparse matrix of complex number.

    It is a square matrix, as many rows (columns) as there are **connected** buses on the grid.

    .. seealso::
        If you want to retrieve the Ybus adopting the "gridmodel" bus labelling, you can use 
        :func:`lightsim2grid.gridmodel.get_dcYbus`

    .. versionadded:: 0.9.0
        It was named `get_dcYbus` before this version, but the name has been changed to avoid confusing AND a new
        function (this one) has been made with the proper `gridmodel` labelling.

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system !

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Ybus[id_me_to_ac_solver[k],:]` (rows of this bus), `Ybus[:, id_me_to_ac_solver[k]]` (column for this bus) 

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter";


const std::string DocLSGrid::get_Sbus_solver = R"mydelimiter(
    This function returns the (complex) `Sbus` vector used by the AC solver. 
    It is the vector of active / reactive power injected at each active bus

    The resulting vector is a vector of complex number having the size of the number of connected buses on the grid.

    .. seealso::
        If you want to retrieve the Sbus with the "gridmodel" convention, you can use :func:`lightsim2grid.gridmodel.get_Sbus`

    .. versionadded:: 0.9.0
        It was named `get_Sbus` before this version, but the name has been changed to avoid confusing AND a new
        function (this one) has been made with the proper `gridmodel` labelling.

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system and in load convention (so generation will be negative)

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
    
    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Sbus[id_me_to_ac_solver[k]]` is the total power injected at the grid model bus solver `k`.

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter"; 


const std::string DocLSGrid::get_dcSbus_solver = R"mydelimiter(
    This function returns the (complex) `Sbus` vector used by the DC sovler. 
    It is the vector of active / reactive power injected at each active bus

    The resulting vector is a vector of complex number having the size of the number of **connected** buses on the grid.

    .. seealso::
        If you want to retrieve the Sbus with the "gridmodel" convention, you can use :func:`lightsim2grid.gridmodel.get_dcSbus`

    .. versionadded:: 0.9.0
        
    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system and in load convention (so generation will be negative)

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
    
    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Sbus[id_me_to_ac_solver[k]]` is the total power injected at the grid model bus solver `k`.

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter"; 


const std::string DocLSGrid::get_Ybus = R"mydelimiter(
    This function returns the (complex) `Ybus` matrix (for the AC powerflow)
    with the gridmodel convention.

    The resulting matrix is a CSC scipy sparse matrix of complex number.

    It is a square matrix, as many rows (columns) as there are **total** buses on the grid.

    .. seealso::
        If you want to retrieve the Ybus adopting the "solver" bus labelling (old behaviour), you can use 
        :func:`lightsim2grid.gridmodel.get_Ybus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_Ybus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.gridmodel.get_Ybus_solver`

    .. warning:: 
        Each row / columns of this matrix represents a "solver bus" (and not a "grid model bus"). In other word, the first row / column of this
        matrix is not necessarily the first bus of the grid model.

    .. warning::
        This is given in the pair unit system !

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Ybus[id_me_to_ac_solver[k],:]` (rows of this bus), `Ybus[:, id_me_to_ac_solver[k]]` (column for this bus) 

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter";


const std::string DocLSGrid::get_dcYbus = R"mydelimiter(
    This function returns the (complex) `Ybus` matrix (for the DC powerflow)
    (its imaginary part should be 0.) with the gridmodel convention.

    The resulting matrix is a CSC scipy sparse matrix of complex number.

    It is a square matrix, as many rows (columns) as there are **total** buses on the grid.

    .. seealso::
        If you want to retrieve the Ybus adopting the "solver" bus labelling (old behaviour), you can use 
        :func:`lightsim2grid.gridmodel.get_dcYbus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_dcYbus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.gridmodel.get_dcYbus_solver`

    .. warning::
        This is given in the pair unit system !

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.

    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Ybus[id_me_to_ac_solver[k],:]` (rows of this bus), `Ybus[:, id_me_to_ac_solver[k]]` (column for this bus) 

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter";


const std::string DocLSGrid::get_Sbus = R"mydelimiter(
    This function returns the (complex) `Sbus` vector of the gridmodel. It is build
    using the "Sbus" passed to the AC solver for which the buses have been properly relabelled
    in the gridmodel convention.

    The resulting vector is a vector of complex number having the size of the number of **total** buses on the grid.

    .. seealso::
        If you want to retrieve the Sbus with the "solver" convention, you can use 
        :func:`lightsim2grid.gridmodel.get_Sbus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_Sbus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.gridmodel.get_Sbus_solver`

    .. warning::
        This is given in the pair unit system and in load convention (so generation will be negative)

    .. seealso:: :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
    
    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Sbus[id_me_to_ac_solver[k]]` is the total power injected at the grid model bus solver `k`.

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter"; 


const std::string DocLSGrid::get_dcSbus = R"mydelimiter(
    This function returns the (complex) `Sbus` vector of the gridmodel for the DC solver (imaginary part should be 0.). 
    It is build using the "dcSbus" passed to the DC solver for which the buses have been properly relabelled
    in the gridmodel convention.

    The resulting vector is a vector of complex number having the size of the number of **total** buses on the grid.

    .. seealso::
        If you want to retrieve the Sbus with the "sovler" convention, you can use 
        :func:`lightsim2grid.gridmodel.get_dcSbus_solver`

    .. versionadded:: 0.9.0

    .. warning::
        This is given in the pair unit system and in load convention (so generation will be negative)

    .. seealso:: 
        :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for
        ways to link the "grid model" bus id to the "solver" bus id.
    
    Notes
    -----

    Suppose that the grid model bus of id k is connected. Then the row / column `id_me_to_ac_solver[k]` (will be >= 0) and will represent this bus:
    `Sbus[id_me_to_ac_solver[k]]` is the total power injected at the grid model bus solver `k`.

    .. warning:: 
        The above only holds when the bus of id `k` is connected which is when `id_me_to_ac_solver[k] >= 0` !

)mydelimiter"; 


const std::string DocLSGrid::check_solution = R"mydelimiter(
    This function allows to check that a given complex voltage vector satisfies the KCL or not, given the state of the sytem.

    .. note::
        It is expected that you provide a complex number even for the buses that are disconnected in the grid model. They will not be ignored
        so you can put anything you want. We keep the public interface this way to avoid headaches with the bus order between
        the grid model and the solver (you can refer to :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and 
        :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` if you still want to have a look)

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

const std::string DocLSGrid::deactivate_result_computation = R"mydelimiter(
    Allows to deactivate the computation of the flows, reactive power absorbed by generators etc. to gain a bit of time when it is not needed.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.reactivate_result_computation`
)mydelimiter";     

const std::string DocLSGrid::reactivate_result_computation = R"mydelimiter(
    Allows to reactivate the computation of the flows, reactive power absorbed by generators etc. when they are needed again after having been
    deactivated.

    .. seealso:: :func:`lightsim2grid.network.LSGrid.deactivate_result_computation`
)mydelimiter";     


const std::string DocLSGrid::ac_pf = R"mydelimiter(
    Allows to perform an AC (alternating current) powerflow.

    .. note::
        It is expected that you provide a complex number even for the buses that are disconnected in the grid model. They will not be affected (if the powerflow converges)
        and you can put anything you want there. We keep the public interface this way to avoid headaches with the bus order between
        the grid model and the solver (you can refer to :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and 
        :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` if you still want to have a look)

    .. seealso:: :func:`lightsim2grid.network.LSGrid.dc_pf` if you want to perform DC powerflow (same interface, same results, same behaviour)

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
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # have an initial guess for the complex voltage at each bus
        Vinit = np.ones(grid_model.total_bus(), dtype=complex)
        
        # maximum number of iteration
        nb_iter = 10  # a good default

        # tolerance
        tol = 1e-8

        V = grid_model.ac_pf(Vinit, nb_iter, tol)
        # if the powerflow has converged, V.shape > 0 otherwise V is empty (size 0)
        # the original V is modified in the process !

)mydelimiter";     

const std::string DocLSGrid::dc_pf = R"mydelimiter(
    This function has the same interface, inputs, outputs, behaviour, etc.
    as the :func:`lightsim2grid.network.LSGrid.ac_pf`.
)mydelimiter";       


const std::string DocLSGrid::get_ptdf = R"mydelimiter(
    This function returns the PTDF (Power Transfer Distribution Factor) which tells you
    how much the flows on each powerline / tranformer will vary if some given power
    is injected on each bus of the grid.

    It adopts the `gridmodel` bus labelling.

    It is a dense matrix, with (nb lines + nb tranformers) rows and (nb **total** bus) columns.

    .. note::
        You need to run a DC powerflow before calling this method (otherwise 
        an exception is raised.)

    It is an alternative to compute DC powerflows (provided that the topology of the 
    grid is not modified). You can do it with:

    .. code-block:: python

        import numpy as np
        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # have an initial guess for the complex voltage at each bus
        Vinit = np.ones(grid_model.total_bus(), dtype=complex)

        Vdc = grid_model.dc_pf(Vinit, 1, 1e-8)

        PTDF = grid_model.get_ptdf()

        new_Sbus = 1.7 * grid_model.get_dcSbus()

        new_flows = np.dot(PTDF, new_Sbus * grid_model.get_sn_mva())
        # the flows on the grid if every injection is multiplied by 1.7

    .. note::
        If a bus is disconnected, then the associated columns is full of 0.

    .. note::
        If the vector Sbus does not sum to 0. the "slack" used is the first slack of
        the slack vector. No distributed slack is used for DC at the moment.

        If you want distributed slack in this case, please open a feature request on 
        github.

    .. note::
        The 'power' "injected" at disconnected buses (buses with colums of PTDF full of 0.)
        is completely discarded (multiplied by 0.)

)mydelimiter";       

const std::string DocLSGrid::get_ptdf_solver = R"mydelimiter(
    This function returns the PTDF (Power Transfer Distribution Factor) which tells you
    how much the flows on each powerline / tranformer will vary if some given power
    is injected on each bus of the grid.

    It adopts the `solver` bus labelling.

    It is a dense matrix, with (nb lines + nb tranformers) rows and (nb **activated** bus) columns.

    Each rows represents a powerline (or a transformer) and each columns represent a bus.

    So the coefficient at row `i` and column `j` of this matrix represents the increase of
    flow (in MW) of powerline `i` if the power on bus `j` is increased of 1MW.

    .. note::
        First `len(gridmodel.get_lines())` rows represent the powerlines, the remaining
        `len(gridmodel.get_trafos())` represent transformers.

    .. note::
        You need to run a DC powerflow before calling this method (otherwise 
        an exception is raised.)

    It is an alternative to compute DC powerflows (provided that the topology of the 
    grid is not modified). You can do it with:

    .. code-block:: python

        import numpy as np
        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # have an initial guess for the complex voltage at each bus
        Vinit = np.ones(grid_model.total_bus(), dtype=complex)

        Vdc = grid_model.dc_pf(Vinit, 1, 1e-8)

        PTDF = grid_model.get_ptdf_solver()

        new_Sbus = 1.7 * grid_model.get_dcSbus_solver()

        new_flows = np.dot(PTDF, new_Sbus * grid_model.get_sn_mva())
        # the flows on the grid if every injection is multiplied by 1.7
        # spoiler: it will be multiplied by 1.7, but you get the idea, 
        # you can change Sbus in a different ways...

    .. note::
        If a bus is disconnected, then the associated columns is full of 0.

    .. note::
        If the vector Sbus does not sum to 0. the "slack" used is the first slack of
        the slack vector. No distributed slack is used for DC at the moment.

        If you want distributed slack in this case, please open a feature request on 
        github.

    .. note::
        With this convention, the disconnected bus are not modeled.

)mydelimiter";     

const std::string DocLSGrid::get_lodf = R"mydelimiter(
    This function returns the LODF (Line Outage Distribution Factor) which tells you
    how much the flows on each powerline / tranformer will vary if some given 
    powerline / transformer is disconnected.

    It is a dense matrix, with (nb lines + nb tranformers) rows and (nb lines + nb tranformers)
    columns.

    Each rows / columns represent a powerline / transformers. More concretely, the coefficient
    at row `i` and column `j` represents how much the flows on line / transformer `i` will vary
    if line / transformer `j` is disconnected. 

    .. note::
        First `len(gridmodel.get_lines())` rows / columns represent the powerlines, the remaining
        `len(gridmodel.get_trafos())` represent transformers.

    .. note::
        You need to run a DC powerflow before calling this method (otherwise 
        an exception is raised.)

        Internally, this method will compute the PTDF

    It is an alternative to compute DC powerflows when powerlines are disconnected.

    .. code-block:: python

        import numpy as np
        # create a grid model
        import grid2op
        from lightsim2grid import LightSimBackend
        env_name = ...  # eg "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        # have an initial guess for the complex voltage at each bus
        Vinit = np.ones(grid_model.total_bus(), dtype=complex)

        Vdc = grid_model.dc_pf(Vinit, 1, 1e-8)

        LODF_mat = 1. * grid_model.get_lodf()

        lor_p, *_ = grid_model.get_lineor_res()
        tor_p, *_ = grid_model.get_trafohv_res()
        init_powerflow = np.concatenate((lor_p, tor_p))

        # if you want to see the impact of a single line disconnected
        l_id = 0 # (or anything between 0 and n_line + n_trafo)
        por_lodf = init_powerflow + LODF_mat[:, l_id] * init_powerflow[l_id]

        # the effect when disconnecting all powerlines (one powerline disconnected each steps)
        mat_flow = np.tile(init_powerflow, LODF_mat.shape[0]).reshape(LODF_mat.shape)
        por_lodf = mat_flow + LODF_mat.T * mat_flow.T

)mydelimiter";     

const std::string DocLSGrid::get_Bf = R"mydelimiter(
    Returns the "Bus from" matrix, with the bus having the 
    `gridmodel` id (sparse matrix).

    More specifically, it is a matrix with `(nb line + nb trafo)` rows and
    (nb **total** bus) columns.

    For each powerline / transformer (row `i`), there is a `+1` for the 
    "origin side" bus and a `-1` for the "extremity side" bus if the 
    line / trafo is connected. If it is disconnected then the associated 
    row will be full of 0.

    .. note::
        First `len(gridmodel.get_lines())` rows represent the powerlines, the remaining
        `len(gridmodel.get_trafos())` represent transformers.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.get_Bf_solver` which will give
        the same matrix but with buses with the "solver" labelling (thus having
        no columns of 0)

)mydelimiter";

const std::string DocLSGrid::get_Bf_solver = R"mydelimiter(
    Returns the "Bus from" matrix, with the bus having the 
    `solver` id (sparse matrix).

    More specifically, it is a matrix with `(nb line + nb trafo)` rows and
    (nb **connected** bus) columns.

    For each powerline / transformer (row `i`), there is a `+1` for the 
    "origin side" bus and a `-1` for the "extremity side" bus if the 
    line / trafo is connected. If it is disconnected then the associated 
    row will be full of 0.

    .. note::
        First `len(gridmodel.get_lines())` rows represent the powerlines, the remaining
        `len(gridmodel.get_trafos())` represent transformers.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.get_Bf` which will give
        the same matrix but with the buses having the "gridmodel" labelling

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
    Total time spent only in solving the powerflows 
    (excluding pre processing the data, post processing them, initializing everything etc.)
    
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
        This function must be called before :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_flows` and 
        :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_flows`, :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_voltages` or
        :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_sbuses`.

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
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocComputers::compute_power_flows = R"mydelimiter(
    Retrieve the active flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocComputers::get_flows = R"mydelimiter(
    Get the current flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a time step, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been caleed)

    Returns
    -------
    As: ``numpy.ndarry`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocComputers::get_power_flows = R"mydelimiter(
    Get the active flows (in MW) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a time step, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_power_flows` has been called.
        (`compute_power_flows` also requires that :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been caleed)

    Returns
    -------
    Ps: ``numpy.ndarry`` (matrix)
        The active flows (in MW) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocComputers::get_voltages = R"mydelimiter(
    Get the complex voltage angles at each bus of the powergrid.

    Each rows correspond to a time step, each column to a bus.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs`.

    Returns
    -------
    Vs: ``numpy.ndarry`` (matrix)
        The complex voltage angles at each bus of the powergrid.

)mydelimiter";

const std::string DocComputers::get_sbuses = R"mydelimiter(
    Get the complex power injected at each (solver id) bus of the powergrid. Results are given in pair unit.
    We do not recommend to use it as it uses the solver id and NOT the powergrid bus id (you can refer to 
    :func:`lightsim2grid.network.LSGrid.id_me_to_ac_solver` and 
    :func:`lightsim2grid.network.LSGrid.id_ac_solver_to_me` for more information)

    Each rows correspond to a time step, each column to a bus (bus are identified by their solver id !)

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs`.

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
    :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis` class.

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

    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_n1`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_all_n1`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_nk`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_multiple_n1`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.clear`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.remove_n1`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.remove_nk`
    - :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.remove_multiple_n1`

    2) Then you can start the computation of the security analysis with 
    :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` then optionally
    :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute_flows` .

    3) And finally inspect the results with :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows` and 
    :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_voltages` .

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

    .. seealso:: :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_n1` to add only a single line

    .. seealso::
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call 
        (but not necessarily all)

)mydelimiter";

const std::string DocSecurityAnalysis::add_n1 = R"mydelimiter(
    This allows to add a single  "n-1" in the contingency list to simulate.

    .. seealso:: :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_all_n1` to add all contingencies at the same time

    .. seealso::
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call.

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
    and is equivalent to call multiple times :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_n1`

    .. seealso::
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_all_n1` to add all the "n-1" contingencies.
    
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
        This function must be called before :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute_flows` and 
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows` or
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_voltages` .

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
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocSecurityAnalysis::compute_power_flows = R"mydelimiter(
    Compute the current flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocSecurityAnalysis::get_flows = R"mydelimiter(
    Get the flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a contingency, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been caleed)

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.get_flows` method.

    Returns
    -------
    As: ``numpy.ndarray`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

const std::string DocSecurityAnalysis::get_voltages = R"mydelimiter(
    Get the complex voltage angles at each bus of the powergrid.

    Each rows correspond to a contingency, each column to a bus.

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute`.

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.get_flows` method.

    Returns
    -------
    Vs: ``numpy.ndarray`` (matrix)
        The complex voltage angles at each bus of the powergrid.

)mydelimiter";

const std::string DocSecurityAnalysis::get_power_flows = R"mydelimiter(
    Get the active flows (in MW) at the origin side / high voltage side of each transformers / powerlines.

    Each rows correspond to a contingency, each column to a powerline / transformer

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute_power_flows` has been called.
        (`compute_flows` also requires that :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been caleed)

    .. warning::
        The order in which the contingencies are computed is **NOT** (in this c++ class) the order in which you enter them. They are computed
        in the order given by :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.my_defaults`. For an easier, more "human readable" please
        use the :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis.get_flows` method.

    Returns
    -------
    As: ``numpy.ndarray`` (matrix)
        The flows (in kA) at the origin side / high voltage side of each transformers / powerlines.

)mydelimiter";

} // namespace ls2g
