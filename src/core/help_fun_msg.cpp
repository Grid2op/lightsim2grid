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
        C++ code path used by :class:`lightsim2grid.network.LSGrid` and the batch solvers
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

const std::string DocSolver::TimerJac = R"mydelimiter(
    Named timer record returned by
    :func:`~lightsim2grid.algorithm.AlgorithmSelector.get_timers_jacobian`, breaking down the plain
    ``(timer_Fx, timer_solve, timer_check, timer_total_nr)`` tuple of
    :func:`~lightsim2grid.algorithm.NR_SparseLU.get_timers` into a finer-grained, named set of
    phases.

    All fields default to ``-1.`` when not measured by the active solver: only the
    Newton-Raphson family (``NR_*``/``NRSing_*``/``NRRefactorRetry_*``) fills in every field;
    Gauss-Seidel and DC solvers only ever set the handful of phases they actually go through
    (see each field's own doc for which).

    Supports tuple-style iteration, indexing and unpacking, in the field declaration order above.

)mydelimiter";

const std::string DocSolver::timer_Fx = R"mydelimiter(
    Time spent computing the KCL mismatch (both active and reactive power) for each bus.

)mydelimiter";

const std::string DocSolver::timer_solve = R"mydelimiter(
    Total time spent in the underlying linear solver's ``solve()`` step (same value as
    :attr:`lightsim2grid.algorithm.LinearSolverStats.timer_solve` for NR-based solvers).

)mydelimiter";

const std::string DocSolver::timer_factor = R"mydelimiter(
    Total time spent in the underlying linear solver's ``factorize()`` step (same value as
    :attr:`lightsim2grid.algorithm.LinearSolverStats.timer_factor`). NR-only: ``-1.`` for
    Gauss-Seidel and DC solvers.

)mydelimiter";

const std::string DocSolver::timer_refactor = R"mydelimiter(
    Total time spent in the underlying linear solver's ``refactorize()`` step (same value as
    :attr:`lightsim2grid.algorithm.LinearSolverStats.timer_refactor`). NR-only: ``-1.`` for
    Gauss-Seidel and DC solvers.

)mydelimiter";

const std::string DocSolver::timer_initialize = R"mydelimiter(
    Total time spent in the underlying linear solver's ``analyze()`` (symbolic factorization)
    step (same value as :attr:`lightsim2grid.algorithm.LinearSolverStats.timer_initialize`).
    NR-only: ``-1.`` for Gauss-Seidel and DC solvers.

)mydelimiter";

const std::string DocSolver::timer_check = R"mydelimiter(
    Time spent checking whether the KCL mismatch met the convergence tolerance.

)mydelimiter";

const std::string DocSolver::timer_dSbus = R"mydelimiter(
    NR-only -- time spent computing the bus-injection sensitivities the Jacobian is built from
    (``-1.`` for Gauss-Seidel and DC solvers).

)mydelimiter";

const std::string DocSolver::timer_fillJ = R"mydelimiter(
    NR-only -- time spent assembling the Jacobian matrix itself, from the sensitivities measured
    by :attr:`timer_dSbus` (``-1.`` for Gauss-Seidel and DC solvers).

)mydelimiter";

const std::string DocSolver::timer_Va_Vm = R"mydelimiter(
    NR-only -- time spent updating the voltage angles / magnitudes from the linear solver's solved
    increments (``-1.`` for Gauss-Seidel and DC solvers).

)mydelimiter";

const std::string DocSolver::timer_pre_proc = R"mydelimiter(
    Time spent in pre-processing (setup before the main iteration loop starts).

)mydelimiter";

const std::string DocSolver::timer_scale = R"mydelimiter(
    NR-only -- time spent applying the active step-scaling policy (see
    :func:`~lightsim2grid.algorithm.NR_SparseLU.get_scaling_policy_type`), eg the line-search
    backtracking of the ``LineSearch`` policy (``-1.`` for Gauss-Seidel and DC solvers).

)mydelimiter";

const std::string DocSolver::timer_mismatch = R"mydelimiter(
    Time spent (re)computing the mismatch, including any post-processing done once the main loop
    has finished.

)mydelimiter";

const std::string DocSolver::timer_total_nr = R"mydelimiter(
    Total time spent in the solver (the same value as
    :func:`~lightsim2grid.algorithm.NR_SparseLU.get_timers`'s ``timer_total_nr``).

)mydelimiter";

const std::string DocSolver::LinearSolverStats = R"mydelimiter(
    Per-call counters and timings for a linear solver, as tracked by every built-in solver
    (``LinearSolverPolicy``) and by the ``NRRefactorRetry_*`` solvers' extra fallback bookkeeping
    (``RefactorRetryLinearSolver``).

    Returned by :func:`~lightsim2grid.algorithm.NR_SparseLU.get_linear_solver_stats` (or, for the
    fast-decoupled ``FDPF_*`` family, which holds two independent linear solvers, by
    :func:`~lightsim2grid.algorithm.FDPF_XB_SparseLU.get_linear_solver_stats_bp` /
    :func:`~lightsim2grid.algorithm.FDPF_XB_SparseLU.get_linear_solver_stats_bpp`).

    The ``nb_*`` counters accumulate over the solver's whole lifetime (across every
    :func:`~lightsim2grid.algorithm.NR_SparseLU.compute_pf` call, not reset in between); the
    ``timer_*`` fields reset every ``compute_pf`` call, like
    :func:`~lightsim2grid.algorithm.AlgorithmSelector.get_timers_jacobian`'s
    :class:`~lightsim2grid.algorithm.TimerJac`.

)mydelimiter";

const std::string DocSolver::nb_reset = R"mydelimiter(
    Number of times ``reset()`` was called on the underlying linear solver (discarding any cached
    factorization).

)mydelimiter";

const std::string DocSolver::nb_analyze = R"mydelimiter(
    Number of times the underlying linear solver's ``analyze()`` (symbolic factorization) was
    called.

)mydelimiter";

const std::string DocSolver::nb_factorize = R"mydelimiter(
    Number of times the underlying linear solver's ``factorize()`` (full numeric factorization)
    was called.

)mydelimiter";

const std::string DocSolver::nb_refactorize = R"mydelimiter(
    Number of times the underlying linear solver's ``refactorize()`` (cheaper, reusing the
    existing symbolic factorization / pivot order) was called.

)mydelimiter";

const std::string DocSolver::nb_refactorize_failed = R"mydelimiter(
    Number of times a ``refactorize()`` call failed (eg the matrix became too ill-conditioned for
    the reused pivot order).

    .. seealso::
        :attr:`nb_fallback_factorize`, on the ``NRRefactorRetry_*`` solvers: a failed refactorize
        there is retried with a full ``factorize()`` before giving up.

)mydelimiter";

const std::string DocSolver::nb_fallback_factorize = R"mydelimiter(
    ``NRRefactorRetry_*`` solvers only: number of times a failed ``refactorize()`` was retried
    with a full ``factorize()`` (``0`` for every other solver, which does not retry).

)mydelimiter";

const std::string DocSolver::nb_fallback_factorize_failed = R"mydelimiter(
    ``NRRefactorRetry_*`` solvers only: number of times that fallback ``factorize()`` retry (see
    :attr:`nb_fallback_factorize`) itself failed (``0`` for every other solver).

)mydelimiter";

const std::string DocSolver::nb_solve = R"mydelimiter(
    Number of times the underlying linear solver's ``solve()`` was called.

)mydelimiter";

const std::string DocSolver::get_linear_solver_stats = R"mydelimiter(
    Per-call counters and timings for the underlying linear solver, as a
    :class:`lightsim2grid.algorithm.LinearSolverStats`.

    .. seealso::
        :func:`~lightsim2grid.algorithm.FDPF_XB_SparseLU.get_linear_solver_stats_bp` /
        :func:`~lightsim2grid.algorithm.FDPF_XB_SparseLU.get_linear_solver_stats_bpp`, the
        equivalent for the fast-decoupled ``FDPF_*`` family, which holds two independent linear
        solvers (this method does not exist there).

)mydelimiter";

const std::string DocSolver::get_linear_solver_stats_bp = R"mydelimiter(
    ``FDPF_*`` solvers only: per-call counters and timings for the B' linear solver, as a
    :class:`lightsim2grid.algorithm.LinearSolverStats`.

    .. seealso::
        :func:`get_linear_solver_stats_bpp` for the B'' linear solver;
        :func:`~lightsim2grid.algorithm.NR_SparseLU.get_linear_solver_stats` for the
        single-linear-solver equivalent used by every other solver family.

)mydelimiter";

const std::string DocSolver::get_linear_solver_stats_bpp = R"mydelimiter(
    ``FDPF_*`` solvers only: per-call counters and timings for the B'' linear solver, as a
    :class:`lightsim2grid.algorithm.LinearSolverStats`.

    .. seealso::
        :func:`get_linear_solver_stats_bp` for the B' linear solver;
        :func:`~lightsim2grid.algorithm.NR_SparseLU.get_linear_solver_stats` for the
        single-linear-solver equivalent used by every other solver family.

)mydelimiter";

const std::string DocSolver::get_scaling_policy_type = R"mydelimiter(
    Return the current step-scaling policy (:class:`lightsim2grid.algorithm.ScalingPolicyType`):
    how the Newton-Raphson step is scaled down before being applied, if at all.

)mydelimiter";

const std::string DocSolver::set_scaling_policy = R"mydelimiter(
    Set the step-scaling policy (:class:`lightsim2grid.algorithm.ScalingPolicyType`). The
    per-policy parameters below (:func:`set_max_dVa`/:func:`set_max_dVm`, :func:`set_ls_c`/
    :func:`set_ls_rho`/:func:`set_ls_max_iter`, :func:`set_iw_mu_min`/:func:`set_iw_mu_max`) are
    only read by their corresponding policy; changing them has no effect while a different
    policy is active.

)mydelimiter";

const std::string DocSolver::get_refactor_policy = R"mydelimiter(
    Return the current Jacobian refactorization policy
    (:class:`lightsim2grid.algorithm.RefactorPolicyType`): when the linear solver does a cheaper
    ``refactorize()`` instead of a full ``factorize()``.

)mydelimiter";

const std::string DocSolver::set_refactor_policy = R"mydelimiter(
    Set the Jacobian refactorization policy (:class:`lightsim2grid.algorithm.RefactorPolicyType`).
    :func:`set_refactor_every_n` is only read by the ``EveryN`` policy.

)mydelimiter";

const std::string DocSolver::get_max_dVa = R"mydelimiter(
    Maximum voltage angle step (radian) allowed per iteration, for the ``MaxVoltageChange``
    scaling policy. Only read while that policy is active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_max_dVa = R"mydelimiter(
    Set :func:`get_max_dVa`.

)mydelimiter";

const std::string DocSolver::get_max_dVm = R"mydelimiter(
    Maximum voltage magnitude step (pu) allowed per iteration, for the ``MaxVoltageChange``
    scaling policy. Only read while that policy is active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_max_dVm = R"mydelimiter(
    Set :func:`get_max_dVm`.

)mydelimiter";

const std::string DocSolver::get_ls_c = R"mydelimiter(
    Armijo sufficient-decrease constant ``c`` for the ``LineSearch`` scaling policy. Only read
    while that policy is active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_ls_c = R"mydelimiter(
    Set :func:`get_ls_c`.

)mydelimiter";

const std::string DocSolver::get_ls_rho = R"mydelimiter(
    Backtracking factor ``rho`` (in ``(0, 1)``) for the ``LineSearch`` scaling policy. Only read
    while that policy is active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_ls_rho = R"mydelimiter(
    Set :func:`get_ls_rho`.

)mydelimiter";

const std::string DocSolver::get_ls_max_iter = R"mydelimiter(
    Maximum number of backtracking iterations for the ``LineSearch`` scaling policy. Only read
    while that policy is active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_ls_max_iter = R"mydelimiter(
    Set :func:`get_ls_max_iter`.

)mydelimiter";

const std::string DocSolver::get_iw_mu_min = R"mydelimiter(
    Minimum optimal multiplier for the ``Iwamoto`` scaling policy. Only read while that policy is
    active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_iw_mu_min = R"mydelimiter(
    Set :func:`get_iw_mu_min`.

)mydelimiter";

const std::string DocSolver::get_iw_mu_max = R"mydelimiter(
    Maximum optimal multiplier for the ``Iwamoto`` scaling policy. Only read while that policy is
    active (see :func:`set_scaling_policy`).

)mydelimiter";

const std::string DocSolver::set_iw_mu_max = R"mydelimiter(
    Set :func:`get_iw_mu_max`.

)mydelimiter";

const std::string DocSolver::get_refactor_every_n = R"mydelimiter(
    Refactorize (full ``factorize()``, not the cheaper ``refactorize()``) every N-th iteration,
    for the ``EveryN`` refactor policy. Only read while that policy is active (see
    :func:`set_refactor_policy`).

)mydelimiter";

const std::string DocSolver::set_refactor_every_n = R"mydelimiter(
    Set :func:`get_refactor_every_n`.

)mydelimiter";

const std::string DocSolver::get_config = R"mydelimiter(
    Return a :class:`lightsim2grid.algorithm.AlgoConfig` capturing every scaling/refactor policy
    type and parameter above, as a single serializable object.

    .. seealso::
        :func:`set_config` to restore it; going through a
        :class:`~lightsim2grid.lightSimBackend.LightSimBackend` instead of a raw solver object,
        see :func:`lightsim2grid.lightSimBackend.LightSimBackend.get_ac_algo_config`.

)mydelimiter";

const std::string DocSolver::set_config = R"mydelimiter(
    Restore every scaling/refactor policy type and parameter from a
    :class:`lightsim2grid.algorithm.AlgoConfig` previously obtained from :func:`get_config`.

)mydelimiter";

const std::string DocSolver::get_theta_to_J_col = R"mydelimiter(
    ``bus_id -> Jacobian column`` for that bus's voltage-angle (theta) unknown, or ``-1`` if that
    bus has none (eg the slack bus, or a PQ-only DC solve). Only valid after a powerflow has been
    run.

)mydelimiter";

const std::string DocSolver::get_vm_to_J_col = R"mydelimiter(
    ``bus_id -> Jacobian column`` for that bus's voltage-magnitude (Vm) unknown, or ``-1`` if that
    bus has none (eg a PV bus). Only valid after a powerflow has been run.

)mydelimiter";

const std::string DocSolver::get_q_to_J_col = R"mydelimiter(
    ``bus_id -> Jacobian column`` for that bus's reactive-power (Q) unknown -- currently always
    ``-1``: no solver in this version stamps a reactive-power unknown as its own Jacobian column.

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

const std::string DocIterator::get_bus_id = R"mydelimiter(
    ``bus_id`` (see the field of the same name on this container's element type, eg
    :class:`lightsim2grid.elements.GenInfo`) for every element of this container, as a single
    array: element ``i`` of the result is that element's bus id, ``-1`` if disconnected.

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

    .. note::
        For elements that can regulate a voltage instead of applying a fixed reactive setpoint
        (see :attr:`lightsim2grid.elements.GenInfo.voltage_regulator_on` /
        :attr:`lightsim2grid.elements.ConverterStationInfo.voltage_regulator_on`), this value is
        only actually used when voltage regulation is OFF. When it is ON, the reactive power is
        computed by the powerflow instead and this setpoint is ignored.

)mydelimiter";

const std::string DocIterator::voltage_regulator_on = R"mydelimiter(
    Whether this element tries to regulate a bus voltage (PV-like behaviour, following
    :attr:`target_vm_pu`) or applies a fixed reactive setpoint instead (PQ-like behaviour,
    following :attr:`target_q_mvar`).

    When ``True``, the reactive power is not an independent input: it is computed by the
    powerflow so that the regulated bus's voltage magnitude matches ``target_vm_pu`` (within
    ``min_q_mvar`` / ``max_q_mvar``). When ``False``, ``target_q_mvar`` is used directly and
    ``target_vm_pu`` / ``min_q_mvar`` / ``max_q_mvar`` are ignored.

    .. note::
        On a :class:`lightsim2grid.elements.GenInfo`, the regulated bus is not necessarily this
        generator's own bus -- see :attr:`lightsim2grid.elements.GenInfo.regulated_bus_id`
        ("remote voltage control").

        On a :class:`lightsim2grid.elements.ConverterStationInfo`, this is only meaningful for
        VSC stations (:attr:`lightsim2grid.elements.ConverterStationInfo.converter_type` ``== 0``):
        LCC stations (``converter_type == 1``) always have it ``False`` and instead consume
        reactive power following
        :attr:`lightsim2grid.elements.ConverterStationInfo.power_factor`.

)mydelimiter";

const std::string DocIterator::regulated_bus_id = R"mydelimiter(
    The grid bus id whose voltage this element regulates, when :attr:`voltage_regulator_on` is
    ``True``.

    Defaults to this element's own :attr:`bus_id` ("local" voltage control). When it differs
    from ``bus_id``, the element performs "remote voltage control": instead of behaving as an
    ordinary PV bus itself, it acts as a controller contributing (jointly with any other element
    regulating the same bus) to holding that bus's voltage magnitude at ``target_vm_pu``.

    .. seealso::
        :func:`lightsim2grid.network.LSGrid.set_gen_regulated_bus` to change it for a generator.

    .. warning::
        When the grid is read from pypowsybl, the regulated bus is resolved **once**, at import
        time, and stored by its (fixed) lightsim2grid global bus id. If the regulated element is
        later moved to another bus *inside lightsim2grid* (e.g. through a ``change_bus_*`` /
        topology change), the controller keeps regulating the bus resolved at import: the
        lightsim2grid grid and the original pypowsybl grid then desynchronise. Re-import the grid
        (or call ``set_gen_regulated_bus`` again) if you need to follow such a topology change.

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

const std::string DocIterator::yac_11 = R"mydelimiter(
    One entry of this branch's raw two-port AC admittance matrix, computed as if both sides were
    connected (see :attr:`yac_eff_11` for the version that accounts for the actual connection
    status).

    With ``ys = 1 / (r_pu + 1j * x_pu)``, for a plain powerline (``ratio == 1``,
    ``shift_rad == 0``): ``yac_11 = ys + h1``, ``yac_22 = ys + h2``, ``yac_12 = yac_21 = -ys`` --
    see the note in :attr:`r_pu`'s line model. For a transformer, the tap
    :attr:`~lightsim2grid.elements.TrafoInfo.ratio` and
    :attr:`~lightsim2grid.elements.TrafoInfo.shift_rad` additionally fold into all four entries.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::yac_12 = DocIterator::yac_11;
const std::string DocIterator::yac_21 = DocIterator::yac_11;
const std::string DocIterator::yac_22 = DocIterator::yac_11;

const std::string DocIterator::yac_eff_11 = R"mydelimiter(
    One entry of this branch's *effective* two-port AC admittance matrix -- :attr:`yac_11` and
    friends, corrected for the actual connection status. This is exactly what is stamped into the
    grid's Ybus.

    - Both sides connected: equal to :attr:`yac_11` (etc) unchanged.
    - Exactly one side connected (a "half-open" branch): Kron-reduced to a single self-admittance
      at the connected end (the open end is eliminated); the three other entries are ``0``.
    - Neither side connected (or the branch itself disconnected): all four entries are ``0``.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::yac_eff_12 = DocIterator::yac_eff_11;
const std::string DocIterator::yac_eff_21 = DocIterator::yac_eff_11;
const std::string DocIterator::yac_eff_22 = DocIterator::yac_eff_11;

const std::string DocIterator::get_yac_eff_11 = R"mydelimiter(
    ``yac_eff_11`` (etc, see :class:`lightsim2grid.elements.LineInfo` /
    :class:`lightsim2grid.elements.TrafoInfo`) for every element of this container, as a single
    array.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::get_yac_eff_12 = DocIterator::get_yac_eff_11;
const std::string DocIterator::get_yac_eff_21 = DocIterator::get_yac_eff_11;
const std::string DocIterator::get_yac_eff_22 = DocIterator::get_yac_eff_11;

const std::string DocIterator::ydc_11 = R"mydelimiter(
    One entry of this branch's two-port DC admittance matrix -- the DC powerflow linearization
    only keeps the series susceptance (``1 / x_pu``), so ``ydc_11 = ydc_22 = 1 / x_pu`` and
    ``ydc_12 = ydc_21 = -1 / x_pu`` for a plain powerline (a transformer's tap
    :attr:`~lightsim2grid.elements.TrafoInfo.ratio` additionally divides it in). Real numbers,
    unlike the AC :attr:`yac_11` family.

    .. note::
        Unlike :attr:`yac_eff_11`, there is no status-aware "effective" counterpart exposed for
        the DC admittance: a disconnected side is instead handled directly by the DC solver /
        Ybus construction.

)mydelimiter" + DocIterator::line_model;

const std::string DocIterator::ydc_12 = DocIterator::ydc_11;
const std::string DocIterator::ydc_21 = DocIterator::ydc_11;
const std::string DocIterator::ydc_22 = DocIterator::ydc_11;

const std::string DocIterator::SubstationContainer = R"mydelimiter(
    This class allows to iterate through the substations of the :class:`lightsim2grid.network.LSGrid`
    easily, as if they were in a python list.

    A substation is not itself an electrical element: it is the group of candidate buses
    (busbars) that the elements connected "at" a given site can be assigned to (see
    :ref:`bus-labelling` and :attr:`lightsim2grid.elements.SubstationInfo.nb_max_busbars`).

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        for sub in grid_model.get_substations():
            # sub is a `SubstationInfo`
            sub.vn_kv

)mydelimiter";

const std::string DocIterator::SubstationInfo = R"mydelimiter(
    This class represents what you get from retrieving some elements from
    :class:`lightsim2grid.elements.SubstationContainer`.

    It allows to read information from each substation of the powergrid.

    .. warning::
        Data can only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())

        grid_model = env.backend._grid

        first_substation = grid_model.get_substations()[0]  # first substation is a `SubstationInfo`

)mydelimiter";

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
        also used, in the same case, by :func:`lightsim2grid.network.LSGrid.check_solution` when
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
        :func:`lightsim2grid.network.LSGrid.check_solution` does not examine static generators at all
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

const std::string DocIterator::SvcContainer = R"mydelimiter(
    This class allows to iterate through the Static Var Compensators (SVC) of the
    :class:`lightsim2grid.network.LSGrid` easily, as if they were in a python list.

    An SVC injects reactive power only (its active power is always ``0``). It follows the IIDM
    model of powsybl, with three regulation modes (see :attr:`~lightsim2grid.elements.SvcInfo.regulation_mode`):

    - ``VOLTAGE``: regulates the voltage of a bus (local or remote), optionally with a
      voltage/reactive slope ("droop"). Stamps nothing directly in Sbus: it is never a PV bus,
      and is always a controller of a ``VoltageControl`` group (the bordered formulation), even
      for the local, non-sloped case.
    - ``REACTIVE_POWER``: a fixed reactive injection (behaves like a non-regulating generator, or
      a load): stamps Q (and P = 0) into Sbus directly.
    - ``OFF``: behaves as if disconnected.

    :attr:`~lightsim2grid.elements.SvcInfo.b_min` / :attr:`~lightsim2grid.elements.SvcInfo.b_max`
    are stored for introspection only: they are never enforced by the powerflow (no outer loop,
    no limit check), mirroring how a generator's ``min_q_mvar`` / ``max_q_mvar`` is handled.

    Voltage regulation equations (``VOLTAGE`` mode)
    ------------------------------------------------

    A ``VOLTAGE``-mode SVC never becomes a PV bus. Instead it is solved as a "controller" of a
    bordered Newton-Raphson block, exactly like a remote-regulating generator (see
    :attr:`~lightsim2grid.elements.GenInfo.regulated_bus_id`): its reactive injection
    :math:`Q_c` (generator sign convention, per unit) becomes an extra unknown of the powerflow,
    solved for jointly with the bus voltages and angles.

    All controllers (generators and/or SVCs) that regulate the same bus form one "group". For a
    group regulating bus ``reg`` at setpoint :math:`v_{set}`, with controllers :math:`c = 1..N`:

    - voltage constraint (one equation for the whole group):

      .. math::

          V_m(reg) + \sum_{c=1}^{N} s_c \, Q_c = v_{set}

      where :math:`s_c` is the slope (:attr:`~lightsim2grid.elements.SvcInfo.slope_pu`) of
      controller :math:`c`, ``0`` for a generator or a non-sloped (``slope_pu = 0``) SVC. With a
      single non-sloped controller in the group this reduces to the usual PV-like
      :math:`V_m(reg) = v_{set}`, only enforced through this bordered :math:`Q_c` unknown rather
      than by reclassifying the bus.
    - reactive sharing (:math:`N-1` equations, only when the group has more than one controller):

      .. math::

          \frac{Q_1}{w_1} = \frac{Q_2}{w_2} = \dots = \frac{Q_N}{w_N}

      i.e. controllers share the group's total reactive effort in proportion to their weight
      :math:`w_c`, with :math:`w_c` = :attr:`~lightsim2grid.elements.SvcInfo.b_max` :math:`-`
      :attr:`~lightsim2grid.elements.SvcInfo.b_min` for an SVC (``qmax - qmin`` for a generator).

    :attr:`~lightsim2grid.elements.SvcInfo.slope_pu` is expressed directly in per-unit (a pu
    voltage deviation per pu of :math:`Q_c`). When importing from a pypowsybl grid, the slope is
    read from the ``voltagePerReactivePowerControl`` extension in kV/MVar and converted as

    .. math::

          s_{pu} = slope_{kV/MVar} \cdot \frac{s_{n,mva}}{v_{n,kv}(reg)}

    with :math:`v_{n,kv}(reg)` the nominal voltage of the regulated bus.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        for svc in grid_model.get_svcs():
            # svc is a `SvcInfo`
            svc.bus_id

)mydelimiter";

const std::string DocIterator::SvcInfo = R"mydelimiter(
    This class represents what you get from retrieving some elements from
    :class:`lightsim2grid.elements.SvcContainer`.

    It allows to read information from each Static Var Compensator (SVC) of the powergrid.

    .. warning::
        Data can only be accessed from this element. You cannot modify (yet) the grid using this class.

    Examples
    --------

    .. code-block:: python

        import grid2op
        from lightsim2grid import LightSimBackend

        env_name = ... # eg. "l2rpn_case14_test"
        env = grid2op.make(env_name, backend=LightSimBackend())
        grid_model = env.backend._grid

        first_svc = grid_model.get_svcs()[0]  # first SVC is a `SvcInfo`

        for svc in grid_model.get_svcs():
            # svc is a `SvcInfo`
            svc.bus_id

)mydelimiter";

const std::string DocIterator::RegulationMode = R"mydelimiter(
    The regulation mode of a Static Var Compensator (values follow the IIDM model of powsybl):
    ``OFF`` (``0``), ``VOLTAGE`` (``1``), or ``REACTIVE_POWER`` (``2``) -- see
    :class:`lightsim2grid.elements.SvcContainer` for what each means.

)mydelimiter";

const std::string DocIterator::regulation_mode = R"mydelimiter(
    This SVC's regulation mode, as a :class:`~lightsim2grid.elements.SvcContainer.RegulationMode`
    (``0`` = ``OFF``, ``1`` = ``VOLTAGE``, ``2`` = ``REACTIVE_POWER``) -- see
    :class:`lightsim2grid.elements.SvcContainer` for what each means.

)mydelimiter";

const std::string DocIterator::svc_target_vm_pu = R"mydelimiter(
    Voltage setpoint (pu of the regulated bus). Only meaningful in ``VOLTAGE`` mode (see
    :attr:`regulation_mode`).

)mydelimiter";

const std::string DocIterator::svc_target_q_mvar = R"mydelimiter(
    Reactive power setpoint (MVAr, generator sign convention -- positive injects into the grid).
    Only meaningful in ``REACTIVE_POWER`` mode (see :attr:`regulation_mode`).

)mydelimiter";

const std::string DocIterator::slope_pu = R"mydelimiter(
    Voltage/reactive slope ("droop", pu) -- in ``VOLTAGE`` mode, ``0.`` means the SVC holds
    :attr:`target_vm_pu` exactly; a non-zero slope lets the regulated voltage deviate from
    the setpoint in proportion to the reactive power delivered. Unused (but still stored) outside
    ``VOLTAGE`` mode. See :class:`~lightsim2grid.elements.SvcContainer` for the exact voltage
    regulation equations.

)mydelimiter";

const std::string DocIterator::b_min = R"mydelimiter(
    Minimum susceptance (pu) -- stored for introspection only, it is never enforced by the
    powerflow (no outer loop, no limit check).

)mydelimiter";

const std::string DocIterator::b_max = R"mydelimiter(
    Maximum susceptance (pu) -- stored for introspection only, it is never enforced by the
    powerflow (no outer loop, no limit check).

)mydelimiter";

const std::string DocIterator::svc_regulated_bus_id = R"mydelimiter(
    The grid bus id whose voltage this SVC regulates, when :attr:`regulation_mode` is
    ``VOLTAGE``.

    Defaults to this element's own :attr:`bus_id` ("local" voltage control). When it differs
    from ``bus_id``, the SVC performs "remote voltage control": instead of behaving as an
    ordinary PV bus itself, it acts as a controller contributing (jointly with any other element
    regulating the same bus, eg a remote-regulating generator) to holding that bus's voltage
    magnitude at :attr:`target_vm_pu`. Same mechanism as
    :attr:`lightsim2grid.elements.GenInfo.regulated_bus_id`.

    .. warning::
        When the grid is read from pypowsybl, the regulated bus is resolved **once**, at import
        time, and stored by its (fixed) lightsim2grid global bus id. If the regulated element is
        later moved to another bus *inside lightsim2grid* (e.g. through a ``change_bus_*`` /
        topology change), the controller keeps regulating the bus resolved at import: the
        lightsim2grid grid and the original pypowsybl grid then desynchronise. Re-import the grid
        if you need to follow such a topology change.

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

const std::string DocIterator::get_bus_id_side_1 = R"mydelimiter(
    ``bus_1_id`` for every element of this container, as a single array: element ``i`` of the
    result is that element's side-1 bus id, ``-1`` if disconnected on that side.

)mydelimiter";

const std::string DocIterator::get_bus_id_side_2 = R"mydelimiter(
    ``bus_2_id`` for every element of this container, as a single array: element ``i`` of the
    result is that element's side-2 bus id, ``-1`` if disconnected on that side.

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

const std::string DocIterator::converters_mode = R"mydelimiter(
    Which side of the HVDC line rectifies -- ``0`` means side 1 is the rectifier (side 2 the
    inverter), ``1`` means side 2 is the rectifier (side 1 the inverter).

    Active power flows from the rectifier side to the inverter side, minus losses -- see
    :attr:`p_setpoint_mw` and the loss model described in
    :class:`lightsim2grid.elements.HvdcLineInfo`.

)mydelimiter";

const std::string DocIterator::p_setpoint_mw = R"mydelimiter(
    Active power drawn at the rectifier side of the HVDC line (MW, always ``>= 0``): which
    physical side that is depends on :attr:`converters_mode`.

    The power actually delivered at the other (inverter) side is this value minus the resistive
    (:attr:`r_ohm`) and converter
    (:attr:`~lightsim2grid.elements.ConverterStationInfo.loss_factor`) losses -- see the loss
    model described in :class:`lightsim2grid.elements.HvdcLineInfo`.

    .. note::
        When :attr:`droop_enabled` is ``True``, this setpoint is not used: the active power
        instead follows the angle-droop equation, see :attr:`droop_enabled`.

)mydelimiter";

const std::string DocIterator::r_ohm = R"mydelimiter(
    DC line resistance (Ohm), used in the resistive loss term of the loss model described in
    :class:`lightsim2grid.elements.HvdcLineInfo`.

    ``0.`` for lines that do not model a resistive loss (e.g. the legacy pandapower-shaped hvdc
    lines).

)mydelimiter";

const std::string DocIterator::nominal_v_kv = R"mydelimiter(
    DC nominal voltage (kV) of the line, used together with :attr:`r_ohm` in the resistive loss
    term described in :class:`lightsim2grid.elements.HvdcLineInfo`.

    The resistive loss term is ``0.`` when this is ``0.`` (e.g. the legacy pandapower-shaped
    hvdc lines, which do not model DC resistive losses).

)mydelimiter";

const std::string DocIterator::droop_enabled = R"mydelimiter(
    Whether angle-droop control ("AC emulation", IIDM ``HvdcAngleDroopActivePowerControl``) is
    enabled for this HVDC line.

    When ``True``, the active power is not fixed at :attr:`p_setpoint_mw` but instead follows
    the angle difference between the two sides::

        raw_mw = droop_p0_mw + droop_k_mw_per_rad * (theta_1 - theta_2)

    saturated at :attr:`pmax_1to2_mw` / :attr:`pmax_2to1_mw` -- see :attr:`status_droop` for the
    regime currently in effect.

    .. note::
        Angle-droop cannot run once either converter station is individually disconnected while
        the line stays otherwise connected (the remote angle is no longer available): it then
        falls back to the fixed :attr:`p_setpoint_mw` for that line, regardless of this flag.

)mydelimiter";

const std::string DocIterator::droop_p0_mw = R"mydelimiter(
    Angle-droop set point ``p0`` (MW): the active power flow (side 1 to side 2) when the two
    sides' voltage angles are equal.

    Only meaningful when :attr:`droop_enabled` is ``True`` -- see :attr:`droop_enabled` for the
    full equation.

)mydelimiter";

const std::string DocIterator::droop_k_mw_per_rad = R"mydelimiter(
    Angle-droop slope ``k`` (MW per radian of angle difference between the two sides).

    Only meaningful when :attr:`droop_enabled` is ``True`` -- see :attr:`droop_enabled` for the
    full equation.

)mydelimiter";

const std::string DocIterator::pmax_1to2_mw = R"mydelimiter(
    Maximum active power (MW) the angle-droop equation is allowed to deliver from side 1 to
    side 2 before saturating -- see :attr:`droop_enabled` and :attr:`status_droop`.

    Only meaningful when :attr:`droop_enabled` is ``True``.

)mydelimiter";

const std::string DocIterator::pmax_2to1_mw = R"mydelimiter(
    Maximum active power (MW) the angle-droop equation is allowed to deliver from side 2 to
    side 1 before saturating -- see :attr:`droop_enabled` and :attr:`status_droop`.

    Only meaningful when :attr:`droop_enabled` is ``True``.

)mydelimiter";

const std::string DocIterator::status_droop = R"mydelimiter(
    The angle-droop regime currently in effect -- ``0`` means the raw droop equation applies
    unsaturated (linear), ``+1`` means it is saturated at :attr:`pmax_1to2_mw` (flow forced from
    side 1 to side 2), ``-1`` means it is saturated at :attr:`pmax_2to1_mw` (flow forced from
    side 2 to side 1).

    .. note::
        This is an INPUT to the powerflow, not something it decides on its own: switching regime
        changes which equation is stamped in the jacobian, so which regime applies is decided by
        an outer loop (typically in Python, between two solves), not by this solve itself. Use
        :func:`lightsim2grid.network.LSGrid.set_status_droop_hvdc` /
        :func:`lightsim2grid.network.LSGrid.get_status_droop_hvdc` to set / read it at the grid
        level.

    Only meaningful when :attr:`droop_enabled` is ``True``.

)mydelimiter";

const std::string DocIterator::station_side_1 = R"mydelimiter(
    The converter station on side 1 of this HVDC line, as a
    :class:`lightsim2grid.elements.ConverterStationInfo`.

)mydelimiter";

const std::string DocIterator::station_side_2 = R"mydelimiter(
    The converter station on side 2 of this HVDC line, as a
    :class:`lightsim2grid.elements.ConverterStationInfo`.

)mydelimiter";

const std::string DocIterator::ConverterStationInfo = R"mydelimiter(
    This class represents what you get from retrieving one side's converter station of an
    :class:`lightsim2grid.elements.HvdcLineInfo` (:attr:`~lightsim2grid.elements.HvdcLineInfo.station1`
    / :attr:`~lightsim2grid.elements.HvdcLineInfo.station2`).

    It follows the IIDM model of powsybl: a station is either a VSC (voltage source converter --
    behaves like a generator, either regulating voltage or with a fixed reactive setpoint, see
    :attr:`voltage_regulator_on`) or a LCC (line commutated converter -- behaves like a load,
    always consuming ``Q = abs(P) * tan(acos(power_factor))``), see :attr:`converter_type`.

    The active power of a station (:attr:`target_p_mw`, generator sign convention) is not an
    independent input: it is derived from the owning HVDC line's active power setpoint (or its
    angle-droop behaviour) and the loss model -- see :class:`lightsim2grid.elements.HvdcLineInfo`.

    .. warning::
        Data can only be read from this element. You cannot modify (yet) the grid using this
        class directly (see :class:`lightsim2grid.elements.HvdcLineInfo` for how to act on the
        owning HVDC line).

)mydelimiter";

const std::string DocIterator::converter_type = R"mydelimiter(
    Whether this converter station is a VSC (``0``, voltage source converter) or a LCC (``1``,
    line commutated converter).

    See :class:`lightsim2grid.elements.ConverterStationInfo` for the behaviour of each.

)mydelimiter";

const std::string DocIterator::loss_factor = R"mydelimiter(
    Converter loss factor (fraction, between ``0.`` and ``1.``) applied when deriving this
    station's active power from the owning HVDC line's power flow -- see the loss model
    described in :class:`lightsim2grid.elements.HvdcLineInfo`.

)mydelimiter";

const std::string DocIterator::power_factor = R"mydelimiter(
    LCC power factor -- the reactive power consumed by the station is
    ``Q = abs(P) * tan(acos(power_factor))``.

    Only meaningful when :attr:`converter_type` is ``1`` (LCC); always ``1.`` (unused) for VSC
    stations.

)mydelimiter";


const std::string DocLSGrid::LSGrid = R"mydelimiter(
    This class represent a lightsim2grid power network. All the elements that can be manipulated by
    lightsim2grid are represented here.

    We do not recommend to use this class directly, but rather to use a :class:`lightsim2grid.lightSimBackend.LightSimBackend`.

    Examples
    ---------

    We **DO NOT** recommend to do:

    .. code-block:: python

        import lightsim2grid
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower network for example pp_net = pn.case118() 

        grid_model = init_from_pandapower(pp_net)

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
    IndexError
        (C++ ``std::out_of_range``) if an index is out of range.
    RuntimeError
        (C++ ``std::runtime_error``) on a structural inconsistency.

    Returns
    -------
    None
        If the grid is consistent.

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
    Returns the names of all registered algorithms, including any loaded plugins, as a list of
    string.

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

const std::string DocLSGrid::get_ac_algo_controler = R"mydelimiter(
    Return the AC solver family's change-tracking flags, as a
    :class:`lightsim2grid.algorithm.AlgoControl` instance.

    A grid modification (eg. disconnecting a line, changing a setpoint) sets one or more of these
    flags; the AC solver reads and resets them the next time it runs an AC powerflow, so it only
    recomputes what actually changed since the last one. Mostly useful for debugging /
    introspecting exactly what a given modification invalidated.

    .. seealso::
        :func:`get_dc_algo_controler` for the independent set of flags tracked for the DC solver
        family.

)mydelimiter";

const std::string DocLSGrid::get_dc_algo_controler = R"mydelimiter(
    Return the DC solver family's change-tracking flags, as a
    :class:`lightsim2grid.algorithm.AlgoControl` instance.

    Same as :func:`get_ac_algo_controler`, but for the DC solver family: the two are tracked
    independently since a DC powerflow does not consume (and reset) the AC flags, and vice versa.

)mydelimiter";

const std::string DocLSGrid::_ls_to_orig = R"mydelimiter(
    Has the size of the number of possible buses in lightsim2grid (*ie*
    ``n_sub_ * max_nb_bus_per_sub_``) and gives the id of the corresponding bus in the original
    grid (pandapower or pypowsybl).

    If ``-1`` is present, then this bus does not exist in the original grid, it is only present
    in the lightsim2grid gridmodel.

)mydelimiter";

const std::string DocLSGrid::_orig_to_ls = R"mydelimiter(
    Opposite to :attr:`_ls_to_orig`. Has the size of the number of buses in the original grid
    (pandapower or pypowsybl) and tells to which bus of lightsim2grid it corresponds. It should
    be an integer ``>= 0`` and ``< n_sub_ * max_nb_bus_per_sub_``.

)mydelimiter";

const std::string DocLSGrid::_max_nb_bus_per_sub = R"mydelimiter(
    Do not modify it after loading!

)mydelimiter";

const std::string DocLSGrid::_init_kwargs = R"mydelimiter(
    ``dict`` (``str`` -> ``str``) of the relevant kwargs this grid was built with (eg by
    :func:`lightsim2grid.network.init_from_pypowsybl`), for example
    ``{"sort_index": "True", "buses_for_sub": "False"}``. Empty for a grid not built that way, or
    a default-constructed one.

)mydelimiter";

const std::string DocLSGrid::_bus_fusion_rep = R"mydelimiter(
    Fused-bus lookup, size :func:`total_bus` (empty if unset / not built with bus fusion). For
    each ls bus id, gives the ls bus id of the "representative" bus it was merged into by
    ``fuse_zero_impedance_branches`` (identity for a bus not involved in any fusion). Set by
    :func:`lightsim2grid.network.init_from_pypowsybl`; only ever read by downstream Python result
    views (eg ``LightsimResultNetwork``), never by any C++ powerflow logic.

)mydelimiter";

const std::string DocLSGrid::get_ac_algo_config = R"mydelimiter(
    Return the AC solver's :class:`lightsim2grid.algorithm.AlgoConfig` (scaling/refactor policy
    type and parameters).

)mydelimiter";

const std::string DocLSGrid::set_ac_algo_config = R"mydelimiter(
    Apply a :class:`lightsim2grid.algorithm.AlgoConfig` to the AC solver (restores scaling/refactor
    policy and parameters).

)mydelimiter";

const std::string DocLSGrid::get_dc_algo_config = R"mydelimiter(
    Return the DC solver's :class:`lightsim2grid.algorithm.AlgoConfig` (no-op for non-NR solvers,
    returns an empty config).

)mydelimiter";

const std::string DocLSGrid::set_dc_algo_config = R"mydelimiter(
    Apply a :class:`lightsim2grid.algorithm.AlgoConfig` to the DC solver.

)mydelimiter";

const std::string DocLSGrid::change_algorithm_by_name = R"mydelimiter(
    Change the AC (or DC) algorithm by registry name. Accepts built-in names and plugin names
    registered via ``load_solver_plugin()``.

    .. seealso::
        :func:`change_algorithm` to change it by :class:`lightsim2grid.algorithm.AlgorithmType`
        instead.

)mydelimiter";

const std::string DocLSGrid::set_bus_voltage_limits = R"mydelimiter(
    Set the per-bus min/max operating voltage (in kV), one value per bus (see :func:`get_bus_vn_kv`).

)mydelimiter";

const std::string DocLSGrid::get_bus_vmin_kv = R"mydelimiter(
    Per-bus min operating voltage, in kV (``NaN`` if not provided for a given bus, empty array if
    never set).

)mydelimiter";

const std::string DocLSGrid::get_bus_vmax_kv = R"mydelimiter(
    Per-bus max operating voltage, in kV (``NaN`` if not provided for a given bus, empty array if
    never set).

)mydelimiter";

const std::string DocLSGrid::get_svcs = R"mydelimiter(
    Get the container of all the Static Var Compensators (SVC), as a
    :class:`lightsim2grid.elements.SvcContainer`.

)mydelimiter";

const std::string DocLSGrid::turnedoff_no_pv = R"mydelimiter(
    Turned-off generators (or generators with ``target_p_mw == 0``) will not be PV buses: they
    will not maintain voltage.

)mydelimiter";

const std::string DocLSGrid::turnedoff_pv = R"mydelimiter(
    Turned-off generators (or generators with ``target_p_mw == 0``) will be PV buses: they will
    maintain voltage. This is the default.

)mydelimiter";

const std::string DocLSGrid::set_reference_slack_bus = R"mydelimiter(
    Force a (gridmodel) bus to be the angle reference among the slack buses (reordered to
    ``slack_ids[0]``) without changing the slack set / weights; ``-1`` clears it.

)mydelimiter";

const std::string DocLSGrid::get_reference_slack_bus = R"mydelimiter(
    Forced angle-reference slack bus (gridmodel id), or ``-1`` if none.

)mydelimiter";

const std::string DocLSGrid::set_ignore_status_global = R"mydelimiter(
    Ignore the ``global_status`` flags for powerlines and transformers (set to ``True`` if you
    want to control each side of a powerline / transformer independently). Default: ``False``.

)mydelimiter";

const std::string DocLSGrid::set_synch_status_both_side = R"mydelimiter(
    Synchronize the status of each side of a powerline / transformer: if you disconnect one side,
    the other side is also disconnected. Default: ``True``.

)mydelimiter";

const std::string DocLSGrid::get_line_names = R"mydelimiter(
    Names of the powerlines, as set by ``set_line_names``; empty if never set.

)mydelimiter";

const std::string DocLSGrid::get_trafo_names = R"mydelimiter(
    Names of the transformers, as set by ``set_trafo_names``; empty if never set.

)mydelimiter";

const std::string DocLSGrid::set_line_current_limit_side1 = R"mydelimiter(
    Set the side-1 current limit of each powerline, in kA (see
    :attr:`lightsim2grid.elements.LineInfo.limit_a1_ka`).

)mydelimiter";

const std::string DocLSGrid::set_line_current_limit_side2 = R"mydelimiter(
    Set the side-2 current limit of each powerline, in kA (see
    :attr:`lightsim2grid.elements.LineInfo.limit_a2_ka`).

)mydelimiter";

const std::string DocLSGrid::set_trafo_current_limit_side1 = R"mydelimiter(
    Set the side-1 current limit of each transformer, in kA (see
    :attr:`lightsim2grid.elements.TrafoInfo.limit_a1_ka`).

)mydelimiter";

const std::string DocLSGrid::set_trafo_current_limit_side2 = R"mydelimiter(
    Set the side-2 current limit of each transformer, in kA (see
    :attr:`lightsim2grid.elements.TrafoInfo.limit_a2_ka`).

)mydelimiter";

const std::string DocLSGrid::deactivate_powerline_side1 = R"mydelimiter(
    Disconnect only side 1 of a powerline (half-open). Needs
    ``set_synch_status_both_side(False)`` to keep side 2 connected.

)mydelimiter";

const std::string DocLSGrid::deactivate_powerline_side2 = R"mydelimiter(
    Disconnect only side 2 of a powerline (half-open). Needs
    ``set_synch_status_both_side(False)`` to keep side 1 connected.

)mydelimiter";

const std::string DocLSGrid::reactivate_powerline_side1 = R"mydelimiter(
    Reconnect only side 1 of a powerline.

)mydelimiter";

const std::string DocLSGrid::reactivate_powerline_side2 = R"mydelimiter(
    Reconnect only side 2 of a powerline.

)mydelimiter";

const std::string DocLSGrid::deactivate_trafo_side1 = R"mydelimiter(
    Disconnect only side 1 of a transformer (half-open). Needs
    ``set_synch_status_both_side(False)`` to keep side 2 connected.

)mydelimiter";

const std::string DocLSGrid::deactivate_trafo_side2 = R"mydelimiter(
    Disconnect only side 2 of a transformer (half-open). Needs
    ``set_synch_status_both_side(False)`` to keep side 1 connected.

)mydelimiter";

const std::string DocLSGrid::reactivate_trafo_side1 = R"mydelimiter(
    Reconnect only side 1 of a transformer.

)mydelimiter";

const std::string DocLSGrid::reactivate_trafo_side2 = R"mydelimiter(
    Reconnect only side 2 of a transformer.

)mydelimiter";

const std::string DocLSGrid::get_lines_status_side1 = R"mydelimiter(
    Per-side status of each powerline's side 1 (relevant for half-open lines: ``get_lines_status()``
    is ``True`` as soon as either side is connected).

)mydelimiter";

const std::string DocLSGrid::get_lines_status_side2 = R"mydelimiter(
    Per-side status of each powerline's side 2, see :func:`get_lines_status_side1`.

)mydelimiter";

const std::string DocLSGrid::get_trafo_status_side1 = R"mydelimiter(
    Per-side status of each transformer's side 1, see :func:`get_lines_status_side1`.

)mydelimiter";

const std::string DocLSGrid::get_trafo_status_side2 = R"mydelimiter(
    Per-side status of each transformer's side 2, see :func:`get_lines_status_side1`.

)mydelimiter";

const std::string DocLSGrid::deactivate_dcline_side1 = R"mydelimiter(
    Disconnect only converter station 1 of an HVDC line; station 2 stays active (injecting /
    regulating).

)mydelimiter";

const std::string DocLSGrid::deactivate_dcline_side2 = R"mydelimiter(
    Disconnect only converter station 2 of an HVDC line; station 1 stays active (injecting /
    regulating).

)mydelimiter";

const std::string DocLSGrid::change_shift_trafo = R"mydelimiter(
    Change the phase-shift angle for a given transformer.

    .. warning::
        It should be expressed in radian (not in degree) -- see :func:`change_shift_trafo_deg`
        for the degree variant.

    If the flag ``ignore_tap_side_for_shift`` (*eg*
    ``lightsim_grid_model.get_trafos().ignore_tap_side_for_shift``) is ``False`` (the default),
    the angle is given at the tap side (side 1 or side 2). If this flag is ``True`` (*eg* the
    grid comes from pandapower) the phase-shift angle should instead be given at side 1 (the hv
    side in pandapower).

)mydelimiter";

const std::string DocLSGrid::change_shift_trafo_deg = R"mydelimiter(
    Same as :func:`change_shift_trafo` but the phase-shift angle is expressed in degree, not in
    radian.

)mydelimiter";

const std::string DocLSGrid::set_trafo_shift_dependent_rx = R"mydelimiter(
    Declare that (some) transformers have a series impedance (``r``, ``x``) that depends on their
    phase-shift angle ``alpha``, supplied as a per-transformer table of sample points
    ``alpha (rad) -> r/x correction (%)`` (the per-step r/x deltas of a pypowsybl phase-tap-changer;
    ``r%`` == ``x%``).

    The effective ``r`` / ``x`` is then ``base * (1 + corr(shift) / 100)``, interpolated on the
    current shift and refreshed whenever :func:`change_shift_trafo` / ``change_ratio_trafo`` is
    called. There is NO "tap" concept here: the dependency is purely on the (continuous) shift.

    Pass an empty list for a transformer without such a dependency; ``enable`` should be kept
    ``False`` for pandapower grids, which have no such data.

)mydelimiter";

const std::string DocLSGrid::set_gen_regulated_bus = R"mydelimiter(
    Set the grid bus whose voltage a generator regulates ("remote voltage control", see
    :attr:`lightsim2grid.elements.GenInfo.regulated_bus_id`; ``bus == own bus`` for local control).

)mydelimiter";

const std::string DocLSGrid::set_status_droop_hvdc = R"mydelimiter(
    Set the angle-droop regime of an HVDC line (see
    :attr:`lightsim2grid.elements.HvdcLineInfo.status_droop`): ``0`` = linear, ``+1`` = saturated
    side 1 to side 2, ``-1`` = saturated side 2 to side 1.

    This is an INPUT of the solver, constant across one solve: the saturation logic is meant to be
    run between two solves (a Python outer loop).

)mydelimiter";

const std::string DocLSGrid::get_status_droop_hvdc = R"mydelimiter(
    Angle-droop regime of one HVDC line, see :func:`set_status_droop_hvdc`.

)mydelimiter";

const std::string DocLSGrid::get_status_droop_hvdc_vect = R"mydelimiter(
    Angle-droop regimes of every HVDC line, see :func:`set_status_droop_hvdc`.

)mydelimiter";

const std::string DocLSGrid::get_slack_col_solver = R"mydelimiter(
    Jacobian column of the ``MultiSlack`` ``slack_absorbed`` unknown (``-1`` when distributed
    slack is inactive).

)mydelimiter";

const std::string DocLSGrid::get_slack_absorbed_solver = R"mydelimiter(
    Converged value (pu) of the ``MultiSlack`` ``slack_absorbed`` unknown (``0`` when distributed
    slack is inactive). This is the ground truth after convergence -- not the ``0`` initial guess
    an external solver's own linearized derivation starts from.

)mydelimiter";

const std::string DocLSGrid::get_controller_q_solver = R"mydelimiter(
    Converged reactive injection (pu) per ``VoltageControl`` controller (a remote-regulating
    generator or a voltage-mode SVC), in controller registration order. Empty when the extension
    is inactive.

)mydelimiter";

const std::string DocLSGrid::get_controller_kind_solver = R"mydelimiter(
    Kind of each ``VoltageControl`` controller (``0`` = generator, ``1`` = SVC), same order as
    :func:`get_controller_q_solver`.

)mydelimiter";

const std::string DocLSGrid::get_controller_elem_id_solver = R"mydelimiter(
    Element id of each ``VoltageControl`` controller (generator id if a generator, SVC id if an
    SVC), same order as :func:`get_controller_q_solver`.

)mydelimiter";

const std::string DocLSGrid::get_controller_q_col_solver = R"mydelimiter(
    Jacobian column of each ``VoltageControl`` controller's own Q unknown, same order as
    :func:`get_controller_q_solver`.

    NOT the same as the bus-keyed ``get_q_to_J_col_solver``: that map only keeps the LAST
    controller registered at a given bus, so it silently collides whenever two controllers
    regulate reactive power from the same bus. External solvers rebuilding this bordered block
    must use this instead.

)mydelimiter";

const std::string DocLSGrid::get_p_buses_solver = R"mydelimiter(
    Compact ``(bus, row)`` pair list for the P equations -- the row/col counterpart of
    ``get_p_to_J_row_solver()``, preserving EVERY registration (a bus may appear more than once;
    see ``NRLedger``'s "Multiplicity rules"). Same length as :func:`get_p_rows_solver`.

)mydelimiter";

const std::string DocLSGrid::get_p_rows_solver = R"mydelimiter(
    Jacobian row of each entry in :func:`get_p_buses_solver`, same order.

)mydelimiter";

const std::string DocLSGrid::get_q_buses_solver = R"mydelimiter(
    Compact ``(bus, row)`` pair list for the Q equations, see :func:`get_p_buses_solver`.

)mydelimiter";

const std::string DocLSGrid::get_q_rows_solver = R"mydelimiter(
    Jacobian row of each entry in :func:`get_q_buses_solver`, same order.

)mydelimiter";

const std::string DocLSGrid::get_theta_buses_solver = R"mydelimiter(
    Compact ``(bus, col)`` pair list for the theta unknowns, see :func:`get_p_buses_solver`.

)mydelimiter";

const std::string DocLSGrid::get_theta_cols_solver = R"mydelimiter(
    Jacobian column of each entry in :func:`get_theta_buses_solver`, same order.

)mydelimiter";

const std::string DocLSGrid::get_vm_buses_solver = R"mydelimiter(
    Compact ``(bus, col)`` pair list for the Vm unknowns, see :func:`get_p_buses_solver`.

)mydelimiter";

const std::string DocLSGrid::get_vm_cols_solver = R"mydelimiter(
    Jacobian column of each entry in :func:`get_vm_buses_solver`, same order.

)mydelimiter";

const std::string DocLSGrid::get_hvdc_droop_data_solver = R"mydelimiter(
    ``(bus1, bus2, status, p0, k, lf1, lf2, r, pmax12, pmax21)``, one entry per CONNECTED
    droop-enabled HVDC line (solver bus numbering, pu). Ground truth for external solvers
    re-deriving the theta-dependent droop-flow contribution to F independently -- see
    ``HvdcDroopSolverData`` for the flow formula.

)mydelimiter";

const std::string DocLSGrid::copy = R"mydelimiter(
    Return a full, independent deep copy of this grid.

)mydelimiter";

const std::string DocLSGrid::timer_last_ac_pf = R"mydelimiter(
    Wall-clock time (seconds) of the last :func:`ac_pf` call, from pre-processing through
    result storage -- the whole call, not just the solver's own internal timers (see
    :func:`lightsim2grid.algorithm.NR_SparseLU.get_timers` /
    :func:`~lightsim2grid.algorithm.AlgorithmSelector.get_timers_jacobian` for those). ``0.`` if
    :func:`ac_pf` was never called.

)mydelimiter";

const std::string DocLSGrid::timer_last_dc_pf = R"mydelimiter(
    Same as :attr:`timer_last_ac_pf`, but for the last :func:`dc_pf` call.

)mydelimiter";

const std::string DocLSGrid::get_turnedoff_gen_pv = R"mydelimiter(
    Whether a turned-off generator (or one with ``target_p_mw == 0``) counts as a PV bus, as set
    by :func:`turnedoff_pv` / :func:`turnedoff_no_pv` (default: ``True``, ie :func:`turnedoff_pv`).

)mydelimiter";

const std::string DocLSGrid::update_slack_weights = R"mydelimiter(
    Recompute the distributed-slack weight of every generator, restricted to the ones for which
    ``could_be_slack`` is ``True`` (a boolean array, one entry per generator): each such
    generator's weight becomes proportional to its ``abs(target_p_mw)`` (or, if every candidate's
    ``target_p_mw`` is ``0.``, an equal split among them). Every other generator stops
    participating in the slack.

    .. seealso::
        :func:`update_slack_weights_by_id`, the same but taking a list of generator ids instead
        of a boolean mask.

)mydelimiter";

const std::string DocLSGrid::update_slack_weights_by_id = R"mydelimiter(
    Same as :func:`update_slack_weights`, but ``slack_ids`` is a list of candidate generator ids
    instead of a per-generator boolean mask.

)mydelimiter";

const std::string DocLSGrid::assign_slack_to_most_connected = R"mydelimiter(
    Pick a single new slack generator automatically: among the buses with at least one generator
    producing (``target_p_mw > 0``), the one with the most powerline / transformer ends connected
    to it: then, at that bus, the generator with the highest ``abs(target_p_mw)``.

    Clears every existing slack assignment first, so the result is always a single slack
    generator, not a distributed one.

    Returns ``(bus_id, gen_id)`` (gridmodel ids) of the bus and generator picked.

)mydelimiter";

const std::string DocLSGrid::consider_only_main_component = R"mydelimiter(
    Restrict the grid to its main synchronous component: starting a breadth-first search from the
    slack bus(es) over the branch graph (powerlines, transformers, and any other connecting
    element), find every bus reachable from them, then disconnect every element with no bus in
    that component.

    An HVDC line with only one converter station in the main component is not fully disconnected:
    the in-main-component converter stays active (still injecting / regulating its scheduled
    power), and only the out-of-component one is opened -- see
    :class:`lightsim2grid.elements.HvdcLineContainer`.

    Requires at least one slack bus to already be defined (see
    :func:`assign_slack_to_most_connected`); raises otherwise.

)mydelimiter";

const std::string DocLSGrid::get_ignore_status_global = R"mydelimiter(
    Current value of the ``ignore_status_global`` flag, see :func:`set_ignore_status_global`.

)mydelimiter";

const std::string DocLSGrid::get_synch_status_both_side = R"mydelimiter(
    Current value of the ``synch_status_both_side`` flag, see :func:`set_synch_status_both_side`.

)mydelimiter";

const std::string DocLSGrid::set_line_names = R"mydelimiter(
    Set the powerlines' names, one per powerline (raises if the length does not match the number
    of powerlines).

    .. seealso::
        :func:`get_line_names` to read them back.

)mydelimiter";

const std::string DocLSGrid::set_dcline_names = R"mydelimiter(
    Set the HVDC lines' names, one per HVDC line (raises if the length does not match the number
    of HVDC lines).

)mydelimiter";

const std::string DocLSGrid::set_trafo_names = R"mydelimiter(
    Set the transformers' names, one per transformer (raises if the length does not match the
    number of transformers).

    .. seealso::
        :func:`get_trafo_names` to read them back.

)mydelimiter";

const std::string DocLSGrid::set_gen_names = R"mydelimiter(
    Set the generators' names, one per generator (raises if the length does not match the number
    of generators).

)mydelimiter";

const std::string DocLSGrid::set_load_names = R"mydelimiter(
    Set the loads' names, one per load (raises if the length does not match the number of loads).

)mydelimiter";

const std::string DocLSGrid::set_storage_names = R"mydelimiter(
    Set the storage units' names, one per storage unit (raises if the length does not match the
    number of storage units).

)mydelimiter";

const std::string DocLSGrid::set_sgen_names = R"mydelimiter(
    Set the static generators' names, one per static generator (raises if the length does not
    match the number of static generators).

)mydelimiter";

const std::string DocLSGrid::set_shunt_names = R"mydelimiter(
    Set the shunts' names, one per shunt (raises if the length does not match the number of
    shunts).

)mydelimiter";

const std::string DocLSGrid::set_svc_names = R"mydelimiter(
    Set the Static Var Compensators' names, one per SVC (raises if the length does not match the
    number of SVCs).

)mydelimiter";

const std::string DocLSGrid::change_ratio_trafo = R"mydelimiter(
    Change the tap ratio of a given transformer (see
    :attr:`lightsim2grid.elements.TrafoInfo.ratio`).

    .. seealso::
        :func:`change_shift_trafo` / :func:`change_shift_trafo_deg` to change its phase-shift
        angle instead.

)mydelimiter";

const std::string DocLSGrid::get_lines = R"mydelimiter(
    This function allows to retrieve the powerlines (as a 
    :class:`lightsim2grid.elements.LineContainer` object,
    see :ref:`elements-modeled` for more information)

    Examples
    ---------

    .. code-block:: python
        
        # init the grid model
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

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
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the powerlines
        print([el.x_pu for el in lightsim_grid_model.get_dclines()]) # to print the "x" for each powerlines

)mydelimiter";

const std::string DocLSGrid::get_substations = R"mydelimiter(
    This function allows to retrieve the substations (as a
    :class:`lightsim2grid.elements.SubstationContainer` object, see :ref:`elements-modeled` for
    more information). Also available as ``get_voltage_levels`` (its powsybl / IIDM name).

    Examples
    ---------

    .. code-block:: python

        # init the grid model
        from lightsim2grid.network import init_from_pandapower
        pp_net = ...  # any pandapower grid
        lightsim_grid_model = init_from_pandapower(pp_net)  # some warnings might be issued as well as some warnings

        # usage example: print some information about the substations
        print([el.vn_kv for el in lightsim_grid_model.get_substations()]) # to print the nominal voltage of each substation

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
        :func:`lightsim2grid.network.LSGrid.get_Ybus`

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
        :func:`lightsim2grid.network.LSGrid.get_dcYbus`

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
        If you want to retrieve the Sbus with the "gridmodel" convention, you can use :func:`lightsim2grid.network.LSGrid.get_Sbus`

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
        If you want to retrieve the Sbus with the "gridmodel" convention, you can use :func:`lightsim2grid.network.LSGrid.get_dcSbus`

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
        :func:`lightsim2grid.network.LSGrid.get_Ybus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_Ybus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.network.LSGrid.get_Ybus_solver`

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
        :func:`lightsim2grid.network.LSGrid.get_dcYbus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_dcYbus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.network.LSGrid.get_dcYbus_solver`

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
        :func:`lightsim2grid.network.LSGrid.get_Sbus_solver`

    .. danger:: 
        Major change in version 0.9.0 of lightsim2grid (see versionchanged below)

    .. versionchanged:: 0.9.0
        It has not the same definition as the "old" behaviour. In the old behaviour, the `get_Sbus` used the
        solver convention. To get the "old" behaviour, you need to use :func:`lightsim2grid.network.LSGrid.get_Sbus_solver`

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
        :func:`lightsim2grid.network.LSGrid.get_dcSbus_solver`

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

const std::string DocTimeSeries::TimeSeries = R"mydelimiter(
    Allows the computation of time series, that is, the same grid topology is used while the active / reactive power injected
    at each buse vary. The grid topology is fixed, the injections vary.

    This is a "raw" c++ class, for an easier to use interface, please refer to the python documentation of the 
    :class:`lightsim2grid.timeSerie.TimeSerie` class.

)mydelimiter";

const std::string DocTimeSeries::total_time = R"mydelimiter(
    Total time spent in solving the powerflows, pre processing the data, post processing them, initializing everything etc.
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocTimeSeries::solver_time = R"mydelimiter(
    Total time spent only in solving the powerflows 
    (excluding pre processing the data, post processing them, initializing everything etc.)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocTimeSeries::amps_computation_time = R"mydelimiter(
    Time spent in computing the flows (in amps) after the voltages have been computed at each nodes
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocTimeSeries::preprocessing_time = R"mydelimiter(
    Time spent in pre processing the data (this involves, but is not limited to the computation of the Sbus)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocTimeSeries::nb_solved = R"mydelimiter(
    Total number of powerflows solved.

)mydelimiter";

const std::string DocTimeSeries::get_status = R"mydelimiter(
    Status of the solvers (1: success, 0: failure).

    .. note::
        Even if the solver failed at some point, some results might still be available (but not totally).

)mydelimiter";

const std::string DocTimeSeries::compute_Vs = R"mydelimiter(
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

const std::string DocTimeSeries::compute_flows = R"mydelimiter(
    Retrieve the flows (in amps, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocTimeSeries::compute_power_flows = R"mydelimiter(
    Retrieve the active flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.timeSerie.TimeSeriesCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocTimeSeries::get_flows = R"mydelimiter(
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

const std::string DocTimeSeries::get_power_flows = R"mydelimiter(
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

const std::string DocTimeSeries::get_voltages = R"mydelimiter(
    Get the complex voltage angles at each bus of the powergrid.

    Each rows correspond to a time step, each column to a bus.

    .. warning::
        This function must be called after :func:`lightsim2grid.timeSerie.TimeSeriesCPP.compute_Vs`.

    Returns
    -------
    Vs: ``numpy.ndarry`` (matrix)
        The complex voltage angles at each bus of the powergrid.

)mydelimiter";

const std::string DocTimeSeries::get_sbuses = R"mydelimiter(
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

const std::string DocTimeSeries::clear = R"mydelimiter(
    Clear the solver and to as if the class never performed any powerflow.

)mydelimiter";

const std::string DocContingencyAnalysis::ContingencyAnalysis = R"mydelimiter(
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

const std::string DocContingencyAnalysis::preprocessing_time = R"mydelimiter(
    Time spent in pre processing the data (this involves, the checking whether the grid would be still connex after the contingency for example)
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocContingencyAnalysis::modif_Ybus_time = R"mydelimiter(
    Time spent to modify the Ybus matrix before simulating each contingency.
    
    It is given in seconds (``float``).

)mydelimiter";

const std::string DocContingencyAnalysis::add_all_n1 = R"mydelimiter(
    This allows to add all the "n-1" in the contingency list to simulate.

    .. seealso:: :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_n1` to add only a single line

    .. seealso::
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call 
        (but not necessarily all)

)mydelimiter";

const std::string DocContingencyAnalysis::add_n1 = R"mydelimiter(
    This allows to add a single  "n-1" in the contingency list to simulate.

    .. seealso:: :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_all_n1` to add all contingencies at the same time

    .. seealso::
        :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.add_multiple_n1` to add multiple single contingencies in the same call.

    Parameters
    ----------
    line_id: ``int``
        The line id you would like to see disconnected
        
)mydelimiter";

const std::string DocContingencyAnalysis::add_nk = R"mydelimiter(
    This allows to add a single  "n-k" in the contingency list to simulate (it will only add at most one contingency)

    .. warning::
        A "n-k" will disconnect multiple powerlines at the same time. It's not the same as adding muliple "n-1" contingencies, where
        powerlines will be disconnected one after the other.

    Parameters
    ----------
    vect_nk: ``list`` (of ``int``)
        The lines id you want to add in the single contingency added.

)mydelimiter";

const std::string DocContingencyAnalysis::add_multiple_n1 = R"mydelimiter(
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

const std::string DocContingencyAnalysis::clear = R"mydelimiter(
    Clear the list of all contingencies. After a call to this method, you will need to re add some contingencies with

)mydelimiter";

const std::string DocContingencyAnalysis::remove_n1 = R"mydelimiter(
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

const std::string DocContingencyAnalysis::remove_nk = R"mydelimiter(
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

const std::string DocContingencyAnalysis::remove_multiple_n1 = R"mydelimiter(
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

const std::string DocContingencyAnalysis::my_defaults_vect = R"mydelimiter(
    Allows to inspect the contingency list that will be simulated.

    Returns
    -------
    my_defaults_vect: ``list``
        The list (of list) of all the current contingencies. Its length corresponds to the number of contingencies simulated.
        For each contingency, it gives which powerline will be disconnected.

)mydelimiter";

const std::string DocContingencyAnalysis::compute = R"mydelimiter(
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

const std::string DocContingencyAnalysis::compute_flows = R"mydelimiter(
    Compute the current flows (in amps, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocContingencyAnalysis::compute_power_flows = R"mydelimiter(
    Compute the current flows (in MW, at the origin of each powerlines / high voltage size of each transformers.

    .. warning::
        This function must be called after :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.compute` has been called.

    .. note::
        This function must be called before :func:`lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_flows`

    .. note::
        During this computation, the GIL is released, allowing easier parrallel computation

)mydelimiter";

const std::string DocContingencyAnalysis::get_flows = R"mydelimiter(
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

const std::string DocContingencyAnalysis::get_voltages = R"mydelimiter(
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

const std::string DocContingencyAnalysis::get_power_flows = R"mydelimiter(
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

const std::string DocContingencyAnalysis::ViolationElementType = R"mydelimiter(
    The kind of element on which a limit was violated: ``BUS`` (a voltage limit), ``LINE`` /
    ``TRAFO`` (a current limit), or ``GRID`` -- the whole grid / contingency rather than a
    specific element, used when the contingency itself could not be simulated at all (see
    :class:`LimitViolationType`'s ``NOT_SIMULATED`` / ``DIVERGENCE``).

)mydelimiter";

const std::string DocContingencyAnalysis::LimitViolationType = R"mydelimiter(
    The kind of limit that was violated: ``LOW_VOLTAGE`` / ``HIGH_VOLTAGE`` (a bus voltage
    magnitude limit) or ``CURRENT`` (a line / transformer thermal limit) for an ordinary,
    element-level violation; ``NOT_SIMULATED`` or ``DIVERGENCE`` for a contingency-level one (see
    :class:`ViolationElementType`'s ``GRID``):

    - ``NOT_SIMULATED``: a pre-check (eg graph connectivity) skipped this contingency -- the
      solver was never invoked for it.
    - ``DIVERGENCE``: the solver was invoked for this contingency but did not converge.

)mydelimiter";

const std::string DocContingencyAnalysis::LimitViolation = R"mydelimiter(
    A single limit violation, as detected by
    :class:`~lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP`.

    .. seealso::
        :func:`~lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_violations` /
        :func:`~lightsim2grid.contingencyAnalysis.ContingencyAnalysisCPP.get_violations_n`, which
        return a list of these per contingency.

)mydelimiter";

const std::string DocContingencyAnalysis::element_type = R"mydelimiter(
    The kind of element this violation is about, as a :class:`ViolationElementType`.

)mydelimiter";

const std::string DocContingencyAnalysis::element_id = R"mydelimiter(
    Which element this violation is about -- the grid-model bus id for ``BUS``; the local (0-based,
    own-type) line / transformer id for ``LINE`` / ``TRAFO``; unused (``-1``) for ``GRID``.

)mydelimiter";

const std::string DocContingencyAnalysis::side = R"mydelimiter(
    ``1`` or ``2`` for a ``LINE`` / ``TRAFO`` violation (which side's current limit was
    violated); unused (``0``) for ``BUS`` / ``GRID``.

)mydelimiter";

const std::string DocContingencyAnalysis::violation_type = R"mydelimiter(
    The kind of limit that was violated, as a :class:`LimitViolationType`.

)mydelimiter";

const std::string DocContingencyAnalysis::value = R"mydelimiter(
    The value actually reached (the voltage magnitude or the current, matching
    :attr:`violation_type`); unused (``NaN``) for ``NOT_SIMULATED`` / ``DIVERGENCE``.

)mydelimiter";

const std::string DocContingencyAnalysis::limit = R"mydelimiter(
    The limit that was violated; unused (``NaN``) for ``NOT_SIMULATED`` / ``DIVERGENCE``.

)mydelimiter";

const std::string DocContingencyAnalysis::violation_name = R"mydelimiter(
    The violating element's name -- for ``LINE`` / ``TRAFO``, as set by
    :func:`lightsim2grid.network.LSGrid.set_line_names` /
    :func:`lightsim2grid.network.LSGrid.set_trafo_names`; for ``BUS``, the name of the substation
    the violating bus belongs to (see
    :func:`lightsim2grid.network.LSGrid.set_substation_names` -- there is no per-bus name in
    ``LSGrid``, only per-substation ones). Empty string if names were never set on the grid for
    the relevant kind, or for ``GRID``.

)mydelimiter";

const std::string DocMisc::PandaPowerConverter = R"mydelimiter(
    Converts electrical parameters given in the "pandapower" format (percentages, kA, kW, ...)
    into the per-unit ``r`` / ``x`` / ``h`` (and ``r`` / ``x`` / ``h1`` / ``h2``) representation
    lightsim2grid uses internally -- see :attr:`lightsim2grid.elements.LineInfo.r_pu` and friends
    for that representation.

    This is what :func:`lightsim2grid.network.init_from_pandapower` uses under the hood to build
    the powerlines and transformers of the grid it returns; it is not just illustrative, and is
    not needed if you build (or read) a grid any other way.

    Before calling any of the conversion methods below, :func:`set_sn_mva` and :func:`set_f_hz`
    must both have been called with the grid's base power (MVA) and frequency (Hz): every
    conversion method raises if either is unset, non-finite or not strictly positive.

    .. note::
        Every conversion method reads all of its array arguments together, element by element:
        they must all have the same length (one entry per line / transformer), or the call
        raises.

)mydelimiter";

const std::string DocMisc::set_f_hz = R"mydelimiter(
    Set the grid's frequency (Hz, eg ``50.`` or ``60.``), used by the shunt-admittance (``h``)
    computation of :func:`get_line_param` / :func:`get_line_param_legacy`.

    Must be called (with a finite, strictly positive value) before any conversion method.

)mydelimiter";

const std::string DocMisc::set_sn_mva = R"mydelimiter(
    Set the grid's base power (MVA), used to convert every electrical parameter to the per-unit
    system.

    Must be called (with a finite, strictly positive value) before any conversion method.

)mydelimiter";

const std::string DocMisc::get_line_param_legacy = R"mydelimiter(
    Convert a powerline's electrical parameters from the pandapower format to lightsim2grid's
    per-unit ``r`` / ``x`` / ``h`` -- the legacy variant, for pandapower versions where a line's
    shunt admittance is not split between its two sides (see
    :attr:`lightsim2grid.elements.LineInfo.h_pu` for the split ``h1`` / ``h2`` representation
    used internally): the single ``h`` returned here is used, unsplit, for both sides.

    Every argument is a per-line array (one entry per powerline), read together element by
    element: ``branch_r`` / ``branch_x`` (resistance / reactance, Ohm), ``branch_g`` / ``branch_c``
    (shunt conductance in micro-Siemens per km-equivalent / capacitance in nF-equivalent, already
    aggregated over the line's length), ``branch_from_kv`` / ``branch_to_kv`` (nominal voltage of
    the two end buses, kV).

    Returns a ``(r_pu, x_pu, h_pu)`` tuple of per-line arrays.

    .. seealso::
        :func:`get_line_param` for the non-legacy variant, which returns the shunt admittance
        already split per side.

)mydelimiter";

const std::string DocMisc::get_line_param = R"mydelimiter(
    Convert a powerline's electrical parameters from the pandapower format to lightsim2grid's
    per-unit ``r`` / ``x`` / ``h1`` / ``h2`` (see
    :attr:`lightsim2grid.elements.LineInfo.h_pu`): the non-legacy variant, splitting the shunt
    admittance between the two sides.

    Takes the same arguments as :func:`get_line_param_legacy`. Returns a
    ``(r_pu, x_pu, h1_pu, h2_pu)`` tuple of per-line arrays.

    .. warning::
        As of this version, the split is always a plain 50 / 50 share of the total shunt
        admittance (``h1_pu == h2_pu``, each half the legacy, unsplit ``h``): this converter does
        not yet read an asymmetric per-side split from pandapower. Grids with a genuinely
        asymmetric line (eg imported from pypowsybl / IIDM instead) are not affected, since they
        do not go through this converter.

)mydelimiter";

const std::string DocMisc::get_trafo_param_pp3 = R"mydelimiter(
    Convert a transformer's electrical parameters from the pandapower format to lightsim2grid's
    per-unit ``r`` / ``x`` / ``h``, for the "advanced grid model" transformer parametrization
    introduced in pandapower 3 (``pp_net.trafo`` when the grid was built with pandapower >= 3 and
    originates from a pandapower-3-native file).

    ``tap_step_pct`` / ``tap_pos`` / ``tap_angles`` (tap angle in **radian**) / ``is_tap_hv_side``
    describe the tap changer; ``vn_hv`` / ``vn_lv`` are the nominal voltages of the hv / lv buses
    (kV); ``trafo_vk_percent`` / ``trafo_vkr_percent`` / ``trafo_sn_trafo_mva`` /
    ``trafo_pfe_kw`` / ``trafo_i0_pct`` are the transformer's short-circuit / no-load test values,
    as in pandapower. ``trafo_model_is_t`` selects the equivalent-circuit model (``True`` for the
    "T" model, ``False`` for "pi"). Every argument other than ``trafo_model_is_t`` is a
    per-transformer array, read together element by element.

    Returns a ``(r_pu, x_pu, h_pu)`` tuple of per-transformer arrays (a single, unsplit ``h`` --
    see :func:`get_line_param_legacy` for the equivalent caveat on powerlines).

    .. seealso::
        :func:`get_trafo_param_pp2` for the pre-pandapower-3 transformer parametrization.

)mydelimiter";

const std::string DocMisc::get_trafo_param_pp2 = R"mydelimiter(
    Convert a transformer's electrical parameters from the pandapower format to lightsim2grid's
    per-unit ``r`` / ``x`` / ``h``, for the transformer parametrization used before pandapower 3
    (always the "T" equivalent-circuit model; there is no ``trafo_model_is_t`` argument).

    Takes the same arguments as :func:`get_trafo_param_pp3`, minus ``trafo_model_is_t``. Returns
    the same ``(r_pu, x_pu, h_pu)`` tuple of per-transformer arrays.

)mydelimiter";

const std::string DocMisc::AlgoControl = R"mydelimiter(
    Change-tracking flags for one solver family (AC or DC), read via
    :func:`lightsim2grid.network.LSGrid.get_ac_algo_controler` /
    :func:`lightsim2grid.network.LSGrid.get_dc_algo_controler`.

    A grid modification (eg. disconnecting a line, changing a setpoint) sets one or more of these
    flags; the corresponding solver family reads and resets them the next time it runs a
    powerflow, so it only recomputes / re-stamps what actually changed since its last run (a
    plain change in dimension does not, by itself, imply the sparsity pattern changed, for
    instance). Each flag answers one narrow question about what kind of change happened -- none
    of them say *what* changed, only that *something in that category* did.

    .. warning::
        This is read-only introspection: there is no way to set these flags from Python, and
        nothing in lightsim2grid resets them for you except the solver itself, on its next run of
        the family this instance tracks.

)mydelimiter";

const std::string DocMisc::has_dimension_changed = R"mydelimiter(
    Whether the number of buses (so the size of Ybus / Sbus) changed since the last powerflow of
    this solver family -- eg after a topology change that changes how many buses are in use.

)mydelimiter";

const std::string DocMisc::has_pv_changed = R"mydelimiter(
    Whether the set of PV buses (voltage-regulated, see
    :attr:`lightsim2grid.elements.GenInfo.voltage_regulator_on`) changed since the last powerflow
    of this solver family.

)mydelimiter";

const std::string DocMisc::has_pq_changed = R"mydelimiter(
    Whether the set of PQ buses (fixed reactive setpoint) changed since the last powerflow of
    this solver family -- the complement of :func:`has_pv_changed`.

)mydelimiter";

const std::string DocMisc::has_slack_participate_changed = R"mydelimiter(
    Whether the set of generators participating in the distributed slack (see
    :attr:`lightsim2grid.elements.GenInfo.is_slack`) changed since the last powerflow of this
    solver family.

)mydelimiter";

const std::string DocMisc::need_reset_solver = R"mydelimiter(
    Whether the solver needs a full reset (discarding any cached factorization / matrix) before
    its next powerflow -- set for changes too disruptive to recompute incrementally.

)mydelimiter";

const std::string DocMisc::need_recompute_sbus = R"mydelimiter(
    Whether the bus injection vector (Sbus) needs recomputing before the next powerflow of this
    solver family -- eg a setpoint changed, but not necessarily the grid's topology.

)mydelimiter";

const std::string DocMisc::need_recompute_ybus = R"mydelimiter(
    Whether the admittance matrix (Ybus) needs recomputing before the next powerflow of this
    solver family -- some of its coefficients changed, though not necessarily its sparsity
    pattern (see :func:`ybus_change_sparsity_pattern`).

)mydelimiter";

const std::string DocMisc::ybus_change_sparsity_pattern = R"mydelimiter(
    Whether Ybus's sparsity pattern changed (which non-zero entries it has, not just their
    values) since the last powerflow of this solver family -- eg after a topology change. This is
    a stronger condition than :func:`need_recompute_ybus`: a changed sparsity pattern requires the
    linear solver to redo its symbolic factorization, not just its numeric one.

)mydelimiter";

const std::string DocMisc::has_slack_weight_changed = R"mydelimiter(
    Whether the distributed-slack weight (see
    :attr:`lightsim2grid.elements.GenInfo.slack_weight`) of at least one participating generator
    changed since the last powerflow of this solver family.

)mydelimiter";

const std::string DocMisc::has_v_changed = R"mydelimiter(
    Whether at least one generator's voltage setpoint (see
    :attr:`lightsim2grid.elements.GenInfo.target_vm_pu`) changed since the last powerflow of this
    solver family.

)mydelimiter";

const std::string DocMisc::has_ybus_some_coeffs_zero = R"mydelimiter(
    Whether at least one coefficient of Ybus may have been set to exactly ``0.`` (and Ybus
    re-compressed to drop it) since the last powerflow of this solver family. Some solvers
    (notably the Newton-Raphson family) may need to recompute some cached state when this
    happens, even though it does not by itself change Ybus's sparsity pattern.

)mydelimiter";

const std::string DocMisc::has_one_el_changed_bus = R"mydelimiter(
    Whether at least one element (a generator, a load, one side of a line, ...) changed which bus
    it is connected to -- including being reconnected or disconnected -- since the last powerflow
    of this solver family.

)mydelimiter";

const std::string DocMisc::AlgoConfig = R"mydelimiter(
    Serializable configuration blob for Newton-Raphson algorithm parameters: stores the
    scaling / refactor policy type and every associated parameter (see
    :func:`lightsim2grid.algorithm.NR_SparseLU.get_scaling_policy_type` and friends) as a single
    object.

    Obtain via ``solver.get_config()`` or, going through a
    :class:`~lightsim2grid.lightSimBackend.LightSimBackend`,
    :func:`lightsim2grid.lightSimBackend.LightSimBackend.get_ac_algo_config`.

    .. warning::
        :attr:`int_params` / :attr:`real_params` are plain lists returned **by value**:
        ``config.int_params[0] = ...`` silently does nothing, because it mutates a temporary
        copy, not the object's actual state. You must build the new list and reassign the whole
        attribute (``int_params = new_list; config.int_params = int_params``) for the change to
        take effect.

)mydelimiter";

const std::string DocMisc::int_params = R"mydelimiter(
    Integer parameters, as a plain list -- ``[ScalingPolicyType, RefactorPolicyType, ls_max_iter,
    refactor_every_n]`` -- see :class:`~lightsim2grid.algorithm.AlgoConfig`'s warning about
    reassigning the whole list to mutate it.

)mydelimiter";

const std::string DocMisc::real_params = R"mydelimiter(
    Real-valued parameters, as a plain list -- ``[max_dVa, max_dVm, ls_c, ls_rho, iw_mu_min,
    iw_mu_max]`` -- see :class:`~lightsim2grid.algorithm.AlgoConfig`'s warning about reassigning
    the whole list to mutate it.

)mydelimiter";

} // namespace ls2g
