// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "Solvers.hpp"
#include "AlgorithmSelector.hpp"
#include "help_fun_msg.hpp"
#include "powerflow_algorithm/ScalingPolicies.hpp"
#include "powerflow_algorithm/RefactorPolicies.hpp"

using namespace ls2g;

// Bind the common method set shared by ALL solver types.
// Covers: get_Va/Vm/V, get_error, get_nb_iter, reset, converged, compute_pf, get_timers, solve.
// Does NOT include get_J (NR-only) or debug methods (FDPF-only).
template<typename Solver>
void bind_algo_methods(py::class_<Solver>& cls) {
    cls
        .def("get_Va",      &Solver::get_Va,      DocSolver::get_Va.c_str())
        .def("get_Vm",      &Solver::get_Vm,      DocSolver::get_Vm.c_str())
        .def("get_V",       &Solver::get_V,       DocSolver::get_V.c_str())
        .def("get_error",   &Solver::get_error,   DocSolver::get_error.c_str())
        .def("get_nb_iter", &Solver::get_nb_iter, DocSolver::get_nb_iter.c_str())
        .def("reset",       &Solver::reset,       DocSolver::reset.c_str())
        .def("converged",   &Solver::converged,   DocSolver::converged.c_str())
        // python-facing compute_pf / solve are bound to the *checked* wrapper
        // (BaseAlgo::compute_pf_with_input_validation), not the raw virtual
        // compute_pf: the raw one is the hot, unchecked path used internally
        // by LSGrid / BaseBatchSolverSynch (ContingencyAnalysis, TimeSeries,
        // SecurityAnalysis) in tight loops, which never goes through pybind11
        // and so is unaffected by this binding.
        .def("compute_pf",  &Solver::compute_pf_with_input_validation,  py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())
        .def("get_timers",  &Solver::get_timers,  DocSolver::get_timers.c_str())
        .def("solve",       &Solver::compute_pf_with_input_validation,  py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())
        ;
}

// Bind scaling/refactor policy accessors for NRAlgo<> types.
template<typename Solver>
void bind_nr_algo_policies(py::class_<Solver>& cls) {
    cls
        // scaling policy
        .def("get_scaling_policy_type",   &Solver::get_scaling_policy_type,   "Return the current step-scaling policy (ScalingPolicyType)")
        .def("set_scaling_policy",   &Solver::set_scaling_policy,   "Set the step-scaling policy (ScalingPolicyType)",   py::arg("policy"))
        // refactor policy
        .def("get_refactor_policy",  &Solver::get_refactor_policy,  "Return the current Jacobian refactorization policy (RefactorPolicyType)")
        .def("set_refactor_policy",  &Solver::set_refactor_policy,  "Set the Jacobian refactorization policy (RefactorPolicyType)", py::arg("policy"))
        // MaxVoltageChange params
        .def("get_max_dVa",          &Solver::get_max_dVa,          "Max angle step (rad) for MaxVoltageChange policy")
        .def("set_max_dVa",          &Solver::set_max_dVa,          py::arg("value"))
        .def("get_max_dVm",          &Solver::get_max_dVm,          "Max magnitude step (pu) for MaxVoltageChange policy")
        .def("set_max_dVm",          &Solver::set_max_dVm,          py::arg("value"))
        // LineSearch (Armijo) params
        .def("get_ls_c",             &Solver::get_ls_c,             "Armijo sufficient-decrease constant c (LineSearch policy)")
        .def("set_ls_c",             &Solver::set_ls_c,             py::arg("value"))
        .def("get_ls_rho",           &Solver::get_ls_rho,           "Backtracking factor rho in (0,1) (LineSearch policy)")
        .def("set_ls_rho",           &Solver::set_ls_rho,           py::arg("value"))
        .def("get_ls_max_iter",      &Solver::get_ls_max_iter,      "Max backtracking iterations (LineSearch policy)")
        .def("set_ls_max_iter",      &Solver::set_ls_max_iter,      py::arg("value"))
        // Iwamoto params
        .def("get_iw_mu_min",        &Solver::get_iw_mu_min,        "Minimum optimal multiplier (Iwamoto policy)")
        .def("set_iw_mu_min",        &Solver::set_iw_mu_min,        py::arg("value"))
        .def("get_iw_mu_max",        &Solver::get_iw_mu_max,        "Maximum optimal multiplier (Iwamoto policy)")
        .def("set_iw_mu_max",        &Solver::set_iw_mu_max,        py::arg("value"))
        // EveryN param
        .def("get_refactor_every_n", &Solver::get_refactor_every_n, "Refactorize every N-th iteration (EveryN policy)")
        .def("set_refactor_every_n", &Solver::set_refactor_every_n, py::arg("value"))
        // AlgoConfig serialization
        .def("get_config", &Solver::get_config,
            "Return an AlgoConfig capturing all scaling/refactor policy type and parameters.")
        .def("set_config", &Solver::set_config, py::arg("config"),
            "Restore scaling/refactor policy type and parameters from an AlgoConfig.")
        // column (unknown) -> bus-id converters (only valid after a powerflow)
        .def("get_theta_to_J_col", &Solver::get_theta_to_J_col_python,
            "bus_id -> Jacobian column of its voltage-angle (theta) unknown, -1 if none")
        .def("get_vm_to_J_col",    &Solver::get_vm_to_J_col_python,
            "bus_id -> Jacobian column of its voltage-magnitude (vm) unknown, -1 if none")
        .def("get_q_to_J_col",     &Solver::get_q_to_J_col_python,
            "bus_id -> Jacobian column of its reactive (q) unknown, currently always -1")
        ;
}

// Bind get_linear_solver_stats() for solvers holding a single LinearSolver (NR_*/DC_*).
// FDPF_* has two solvers (B'/B'') and is bound separately, see bind_fdpf_linear_solver_stats.
template<typename Solver>
void bind_linear_solver_stats(py::class_<Solver>& cls) {
    cls.def("get_linear_solver_stats", &Solver::get_linear_solver_stats,
            "Per-call counters and timings for the underlying linear solver: "
            "factor/refactor/analyze/solve counts, refactor-failure and "
            "fallback-factor counts, and matching durations (LinearSolverStats). "
            "Counters accumulate over the algorithm's whole lifetime; the timer_* "
            "fields reset every compute_pf call, like get_timers_jacobian().");
}

// Bind the two-solver (B'/B'') equivalent for FDPF_* types.
template<typename Solver>
void bind_fdpf_linear_solver_stats(py::class_<Solver>& cls) {
    cls
        .def("get_linear_solver_stats_bp",  &Solver::get_linear_solver_stats_bp,
            "Per-call counters and timings for the B' linear solver (LinearSolverStats)")
        .def("get_linear_solver_stats_bpp", &Solver::get_linear_solver_stats_bpp,
            "Per-call counters and timings for the B'' linear solver (LinearSolverStats)")
        ;
}

void bind_solvers(py::module_& m) {
    // ---- TimerJac ----
    py::class_<TimerJac>(m, "TimerJac",
        "Named timer record returned by get_timers_jacobian(). "
        "All fields default to -1.0 when not applicable to the active solver. "
        "Supports tuple-style iteration and unpacking.")
        .def_readonly("timer_Fx",         &TimerJac::timer_Fx_)
        .def_readonly("timer_solve",      &TimerJac::timer_solve_)
        .def_readonly("timer_factor",     &TimerJac::timer_factor_)
        .def_readonly("timer_refactor",   &TimerJac::timer_refactor_)
        .def_readonly("timer_initialize", &TimerJac::timer_initialize_)
        .def_readonly("timer_check",      &TimerJac::timer_check_)
        .def_readonly("timer_dSbus",      &TimerJac::timer_dSbus_)
        .def_readonly("timer_fillJ",      &TimerJac::timer_fillJ_)
        .def_readonly("timer_Va_Vm",      &TimerJac::timer_Va_Vm_)
        .def_readonly("timer_pre_proc",   &TimerJac::timer_pre_proc_)
        .def_readonly("timer_scale",      &TimerJac::timer_scale_)
        .def_readonly("timer_mismatch",   &TimerJac::timer_mismatch_)
        .def_readonly("timer_total_nr",   &TimerJac::timer_total_nr_)
        .def("__len__", [](const TimerJac&) { return 13; })
        .def("__getitem__", [](const TimerJac& t, int i) -> double {
            const double vals[13] = {
                t.timer_Fx_, t.timer_solve_, t.timer_factor_,
                t.timer_refactor_, t.timer_initialize_, t.timer_check_,
                t.timer_dSbus_, t.timer_fillJ_, t.timer_Va_Vm_,
                t.timer_pre_proc_, t.timer_scale_, t.timer_mismatch_,
                t.timer_total_nr_
            };
            if (i < 0) i += 13;
            if (i < 0 || i >= 13) throw py::index_error();
            return vals[i];
        })
        .def("__iter__", [](const TimerJac& t) {
            return py::iter(py::make_tuple(
                t.timer_Fx_, t.timer_solve_, t.timer_factor_,
                t.timer_refactor_, t.timer_initialize_, t.timer_check_,
                t.timer_dSbus_, t.timer_fillJ_, t.timer_Va_Vm_,
                t.timer_pre_proc_, t.timer_scale_, t.timer_mismatch_,
                t.timer_total_nr_
            ));
        })
        .def("__repr__", [](const TimerJac& t) {
            return "TimerJac(timer_Fx=" + std::to_string(t.timer_Fx_) +
                   ", timer_solve=" + std::to_string(t.timer_solve_) +
                   ", timer_factor=" + std::to_string(t.timer_factor_) +
                   ", ...)";
        });

    // ---- LinearSolverStats ----
    py::class_<LinearSolverStats>(m, "LinearSolverStats",
        "Per-call counters and timings for a linear solver, as tracked by "
        "LinearSolverPolicy (every built-in solver) and RefactorRetryLinearSolver "
        "(NRRefactorRetry_* solvers). Counters accumulate over the algorithm's whole "
        "lifetime; the timer_* fields reset every compute_pf call.")
        .def_readonly("nb_reset",                      &LinearSolverStats::nb_reset)
        .def_readonly("nb_analyze",                     &LinearSolverStats::nb_analyze)
        .def_readonly("nb_factorize",                   &LinearSolverStats::nb_factorize)
        .def_readonly("nb_refactorize",                 &LinearSolverStats::nb_refactorize)
        .def_readonly("nb_refactorize_failed",          &LinearSolverStats::nb_refactorize_failed)
        .def_readonly("nb_fallback_factorize",          &LinearSolverStats::nb_fallback_factorize)
        .def_readonly("nb_fallback_factorize_failed",   &LinearSolverStats::nb_fallback_factorize_failed)
        .def_readonly("nb_solve",                       &LinearSolverStats::nb_solve)
        .def_readonly("timer_initialize",               &LinearSolverStats::timer_initialize_)
        .def_readonly("timer_factor",                   &LinearSolverStats::timer_factor_)
        .def_readonly("timer_refactor",                 &LinearSolverStats::timer_refactor_)
        .def_readonly("timer_solve",                    &LinearSolverStats::timer_solve_)
        .def("__repr__", [](const LinearSolverStats& s) {
            return "LinearSolverStats(nb_factorize=" + std::to_string(s.nb_factorize) +
                   ", nb_refactorize=" + std::to_string(s.nb_refactorize) +
                   ", nb_refactorize_failed=" + std::to_string(s.nb_refactorize_failed) +
                   ", nb_fallback_factorize=" + std::to_string(s.nb_fallback_factorize) +
                   ", ...)";
        });

    // ---- SparseLU ----
    {
        auto cls = py::class_<NR_SparseLU>(m, "NR_SparseLU", DocSolver::NR_SparseLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_SparseLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRSing_SparseLU>(m, "NRSing_SparseLU", DocSolver::NRSing_SparseLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_SparseLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<DC_SparseLU>(m, "DC_SparseLU", DocSolver::DC_SparseLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_SparseLU>(m, "FDPF_XB_SparseLU", DocSolver::FDPF_XB_SparseLU.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_XB_SparseLU::debug_get_Bp_python,  DocLSGrid::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_XB_SparseLU::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_SparseLU>(m, "FDPF_BX_SparseLU", DocSolver::FDPF_BX_SparseLU.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_BX_SparseLU::debug_get_Bp_python,  DocLSGrid::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_BX_SparseLU::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }

#if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_KLU>(m, "NR_KLU", DocSolver::NR_KLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_KLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRSing_KLU>(m, "NRSing_KLU", DocSolver::NRSing_KLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_KLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<DC_KLU>(m, "DC_KLU", DocSolver::DC_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_KLU>(m, "FDPF_XB_KLU", DocSolver::FDPF_XB_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_KLU>(m, "FDPF_BX_KLU", DocSolver::FDPF_BX_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRRefactorRetry_KLU>(m, "NRRefactorRetry_KLU", DocSolver::NRRefactorRetry_KLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRRefactorRetry_KLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
#endif  // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_NICSLU>(m, "NR_NICSLU", DocSolver::NR_NICSLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_NICSLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRSing_NICSLU>(m, "NRSing_NICSLU", DocSolver::NRSing_NICSLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_NICSLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<DC_NICSLU>(m, "DC_NICSLU", DocSolver::DC_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_NICSLU>(m, "FDPF_XB_NICSLU", DocSolver::FDPF_XB_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_NICSLU>(m, "FDPF_BX_NICSLU", DocSolver::FDPF_BX_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRRefactorRetry_NICSLU>(m, "NRRefactorRetry_NICSLU", DocSolver::NRRefactorRetry_NICSLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRRefactorRetry_NICSLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
#endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_CKTSO>(m, "NR_CKTSO", DocSolver::NR_CKTSO.c_str())
            .def(py::init<>())
            .def("get_J", &NR_CKTSO::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRSing_CKTSO>(m, "NRSing_CKTSO", DocSolver::NRSing_CKTSO.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_CKTSO::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<DC_CKTSO>(m, "DC_CKTSO", DocSolver::DC_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_CKTSO>(m, "FDPF_XB_CKTSO", DocSolver::FDPF_XB_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_CKTSO>(m, "FDPF_BX_CKTSO", DocSolver::FDPF_BX_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
        bind_fdpf_linear_solver_stats(cls);
    }
    {
        auto cls = py::class_<NRRefactorRetry_CKTSO>(m, "NRRefactorRetry_CKTSO", DocSolver::NRRefactorRetry_CKTSO.c_str())
            .def(py::init<>())
            .def("get_J", &NRRefactorRetry_CKTSO::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
        bind_linear_solver_stats(cls);
    }
#endif  // CKTSO_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    {
        auto cls = py::class_<GaussSeidelAlgo>(m, "GaussSeidelAlgo", DocSolver::GaussSeidelAlgo.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<GaussSeidelSynchAlgo>(m, "GaussSeidelSynchAlgo", DocSolver::GaussSeidelSynchAlgo.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }

    // Only "const" methods exported so Python cannot modify the internal solver of a gridmodel
    py::class_<AlgorithmSelector>(m, "AlgorithmSelector", DocSolver::AlgorithmSelector.c_str())
        .def(py::init<>())
        .def("get_type",             &AlgorithmSelector::get_type,             DocSolver::get_type.c_str())
        .def("get_name",             &AlgorithmSelector::get_name,
            "Registry name of the currently active algorithm (eg \"NR_KLU\", or a plugin / "
            "string-only built-in name such as \"NRRefactorRetry_KLU\"). Unlike get_type(), "
            "which reports AlgorithmType.Custom for both plugins and string-only built-ins, "
            "this always identifies the concrete algorithm.")
        .def("get_Va",               &AlgorithmSelector::get_Va,               DocSolver::get_Va.c_str())
        .def("get_Vm",               &AlgorithmSelector::get_Vm,               DocSolver::get_Vm.c_str())
        .def("get_V",                &AlgorithmSelector::get_V,                DocSolver::get_V.c_str())
        .def("get_J",                &AlgorithmSelector::get_J_python,         DocSolver::chooseSolver_get_J_python.c_str())
        .def("get_theta_to_J_col",   &AlgorithmSelector::get_theta_to_J_col_python, "bus_id -> Jacobian column of its voltage-angle (theta) unknown, -1 if none")
        .def("get_vm_to_J_col",      &AlgorithmSelector::get_vm_to_J_col_python,    "bus_id -> Jacobian column of its voltage-magnitude (vm) unknown, -1 if none")
        .def("get_q_to_J_col",       &AlgorithmSelector::get_q_to_J_col_python,     "bus_id -> Jacobian column of its reactive (q) unknown, currently always -1")
        .def("get_error",            &AlgorithmSelector::get_error,            DocSolver::get_V.c_str())
        .def("get_nb_iter",          &AlgorithmSelector::get_nb_iter,          DocSolver::get_nb_iter.c_str())
        .def("converged",            &AlgorithmSelector::converged,            DocSolver::converged.c_str())
        .def("get_computation_time", &AlgorithmSelector::get_computation_time, DocSolver::get_computation_time.c_str())
        .def("get_timers",           &AlgorithmSelector::get_timers,           "TODO")
        .def("get_timers_jacobian",  &AlgorithmSelector::get_timers_jacobian,  "TODO")
        .def("get_timers_ptdf_lodf", &AlgorithmSelector::get_timers_ptdf_lodf, "TODO")
        .def("get_linear_solver_stats", &AlgorithmSelector::get_linear_solver_stats,
            "Per-call counters and timings for the underlying linear solver (LinearSolverStats). "
            "All-zero if the active solver doesn't track them (e.g. GaussSeidel, or the "
            "FDPF family which exposes get_linear_solver_stats_bp/_bpp on its own concrete "
            "Python type instead, since it holds two linear solvers).")
        .def("get_fdpf_xb_lu",       &AlgorithmSelector::get_fdpf_xb_lu,  py::return_value_policy::reference, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_fdpf_bx_lu",       &AlgorithmSelector::get_fdpf_bx_lu,  py::return_value_policy::reference, DocLSGrid::_internal_do_not_use.c_str());
}
