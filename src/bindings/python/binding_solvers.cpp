// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "Solvers.hpp"
#include "ChooseSolver.hpp"
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
        .def("compute_pf",  &Solver::compute_pf,  py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())
        .def("get_timers",  &Solver::get_timers,  DocSolver::get_timers.c_str())
        .def("solve",       &Solver::compute_pf,  py::call_guard<py::gil_scoped_release>(), DocSolver::compute_pf.c_str())
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
        ;
}

void bind_solvers(py::module_& m) {
    // ---- SparseLU ----
    {
        auto cls = py::class_<SparseLUSolver>(m, "SparseLUSolver", DocSolver::SparseLUSolver.c_str())
            .def(py::init<>())
            .def("get_J", &SparseLUSolver::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<SparseLUSolverSingleSlack>(m, "SparseLUSolverSingleSlack", DocSolver::SparseLUSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &SparseLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<DCSolver>(m, "DCSolver", DocSolver::DCSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_SparseLUSolver>(m, "FDPF_XB_SparseLUSolver", DocSolver::FDPF_XB_SparseLUSolver.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_XB_SparseLUSolver::debug_get_Bp_python,  DocGridModel::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_XB_SparseLUSolver::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_SparseLUSolver>(m, "FDPF_BX_SparseLUSolver", DocSolver::FDPF_BX_SparseLUSolver.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_BX_SparseLUSolver::debug_get_Bp_python,  DocGridModel::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_BX_SparseLUSolver::debug_get_Bpp_python, DocGridModel::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
    }

#if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<KLUSolver>(m, "KLUSolver", DocSolver::KLUSolver.c_str())
            .def(py::init<>())
            .def("get_J", &KLUSolver::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<KLUSolverSingleSlack>(m, "KLUSolverSingleSlack", DocSolver::KLUSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &KLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<KLUDCSolver>(m, "KLUDCSolver", DocSolver::KLUDCSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_KLUSolver>(m, "FDPF_XB_KLUSolver", DocSolver::FDPF_XB_KLUSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_KLUSolver>(m, "FDPF_BX_KLUSolver", DocSolver::FDPF_BX_KLUSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
#endif  // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NICSLUSolver>(m, "NICSLUSolver", DocSolver::NICSLUSolver.c_str())
            .def(py::init<>())
            .def("get_J", &NICSLUSolver::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NICSLUSolverSingleSlack>(m, "NICSLUSolverSingleSlack", DocSolver::NICSLUSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &NICSLUSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NICSLUDCSolver>(m, "NICSLUDCSolver", DocSolver::NICSLUDCSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_NICSLUSolver>(m, "FDPF_XB_NICSLUSolver", DocSolver::FDPF_XB_NICSLUSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_NICSLUSolver>(m, "FDPF_BX_NICSLUSolver", DocSolver::FDPF_BX_NICSLUSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
#endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<CKTSOSolver>(m, "CKTSOSolver", DocSolver::CKTSOSolver.c_str())
            .def(py::init<>())
            .def("get_J", &CKTSOSolver::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<CKTSOSolverSingleSlack>(m, "CKTSOSolverSingleSlack", DocSolver::CKTSOSolverSingleSlack.c_str())
            .def(py::init<>())
            .def("get_J", &CKTSOSolverSingleSlack::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<CKTSODCSolver>(m, "CKTSODCSolver", DocSolver::CKTSODCSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_CKTSOSolver>(m, "FDPF_XB_CKTSOSolver", DocSolver::FDPF_XB_CKTSOSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_CKTSOSolver>(m, "FDPF_BX_CKTSOSolver", DocSolver::FDPF_BX_CKTSOSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
#endif  // CKTSO_SOLVER_AVAILABLE (or _READ_THE_DOCS)

    {
        auto cls = py::class_<GaussSeidelAlgo>(m, "GaussSeidelSolver", DocSolver::GaussSeidelSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<GaussSeidelSynchAlgo>(m, "GaussSeidelSynchSolver", DocSolver::GaussSeidelSynchSolver.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }

    // Only "const" methods exported so Python cannot modify the internal solver of a gridmodel
    py::class_<ChooseSolver>(m, "AnySolver", DocSolver::AnySolver.c_str())
        .def(py::init<>())
        .def("get_type",             &ChooseSolver::get_type,             DocSolver::get_type.c_str())
        .def("get_Va",               &ChooseSolver::get_Va,               DocSolver::get_Va.c_str())
        .def("get_Vm",               &ChooseSolver::get_Vm,               DocSolver::get_Vm.c_str())
        .def("get_V",                &ChooseSolver::get_V,                DocSolver::get_V.c_str())
        .def("get_J",                &ChooseSolver::get_J_python,         DocSolver::chooseSolver_get_J_python.c_str())
        .def("get_error",            &ChooseSolver::get_error,            DocSolver::get_V.c_str())
        .def("get_nb_iter",          &ChooseSolver::get_nb_iter,          DocSolver::get_nb_iter.c_str())
        .def("converged",            &ChooseSolver::converged,            DocSolver::converged.c_str())
        .def("get_computation_time", &ChooseSolver::get_computation_time, DocSolver::get_computation_time.c_str())
        .def("get_timers",           &ChooseSolver::get_timers,           "TODO")
        .def("get_timers_jacobian",  &ChooseSolver::get_timers_jacobian,  "TODO")
        .def("get_timers_ptdf_lodf", &ChooseSolver::get_timers_ptdf_lodf, "TODO")
        .def("get_fdpf_xb_lu",       &ChooseSolver::get_fdpf_xb_lu,  py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str())
        .def("get_fdpf_bx_lu",       &ChooseSolver::get_fdpf_bx_lu,  py::return_value_policy::reference, DocGridModel::_internal_do_not_use.c_str());
}
