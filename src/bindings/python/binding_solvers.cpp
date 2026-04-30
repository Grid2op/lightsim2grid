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
        auto cls = py::class_<NR_SparseLU>(m, "NR_SparseLU", DocSolver::NR_SparseLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_SparseLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NRSing_SparseLU>(m, "NRSing_SparseLU", DocSolver::NRSing_SparseLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_SparseLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<DC_SparseLU>(m, "DC_SparseLU", DocSolver::DC_SparseLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_SparseLU>(m, "FDPF_XB_SparseLU", DocSolver::FDPF_XB_SparseLU.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_XB_SparseLU::debug_get_Bp_python,  DocLSGrid::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_XB_SparseLU::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_SparseLU>(m, "FDPF_BX_SparseLU", DocSolver::FDPF_BX_SparseLU.c_str())
            .def(py::init<>())
            .def("debug_get_Bp_python",  &FDPF_BX_SparseLU::debug_get_Bp_python,  DocLSGrid::_internal_do_not_use.c_str())
            .def("debug_get_Bpp_python", &FDPF_BX_SparseLU::debug_get_Bpp_python, DocLSGrid::_internal_do_not_use.c_str());
        bind_algo_methods(cls);
    }

#if defined(KLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_KLU>(m, "NR_KLU", DocSolver::NR_KLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_KLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NRSing_KLU>(m, "NRSing_KLU", DocSolver::NRSing_KLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_KLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<DC_KLU>(m, "DC_KLU", DocSolver::DC_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_KLU>(m, "FDPF_XB_KLU", DocSolver::FDPF_XB_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_KLU>(m, "FDPF_BX_KLU", DocSolver::FDPF_BX_KLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
#endif  // KLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(NICSLU_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_NICSLU>(m, "NR_NICSLU", DocSolver::NR_NICSLU.c_str())
            .def(py::init<>())
            .def("get_J", &NR_NICSLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NRSing_NICSLU>(m, "NRSing_NICSLU", DocSolver::NRSing_NICSLU.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_NICSLU::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<DC_NICSLU>(m, "DC_NICSLU", DocSolver::DC_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_NICSLU>(m, "FDPF_XB_NICSLU", DocSolver::FDPF_XB_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_NICSLU>(m, "FDPF_BX_NICSLU", DocSolver::FDPF_BX_NICSLU.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
#endif  // NICSLU_SOLVER_AVAILABLE (or _READ_THE_DOCS)

#if defined(CKTSO_SOLVER_AVAILABLE) || defined(_READ_THE_DOCS)
    {
        auto cls = py::class_<NR_CKTSO>(m, "NR_CKTSO", DocSolver::NR_CKTSO.c_str())
            .def(py::init<>())
            .def("get_J", &NR_CKTSO::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<NRSing_CKTSO>(m, "NRSing_CKTSO", DocSolver::NRSing_CKTSO.c_str())
            .def(py::init<>())
            .def("get_J", &NRSing_CKTSO::get_J_python, DocSolver::get_J_python.c_str());
        bind_algo_methods(cls);
        bind_nr_algo_policies(cls);
    }
    {
        auto cls = py::class_<DC_CKTSO>(m, "DC_CKTSO", DocSolver::DC_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_XB_CKTSO>(m, "FDPF_XB_CKTSO", DocSolver::FDPF_XB_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
    }
    {
        auto cls = py::class_<FDPF_BX_CKTSO>(m, "FDPF_BX_CKTSO", DocSolver::FDPF_BX_CKTSO.c_str())
            .def(py::init<>());
        bind_algo_methods(cls);
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
        .def("get_Va",               &AlgorithmSelector::get_Va,               DocSolver::get_Va.c_str())
        .def("get_Vm",               &AlgorithmSelector::get_Vm,               DocSolver::get_Vm.c_str())
        .def("get_V",                &AlgorithmSelector::get_V,                DocSolver::get_V.c_str())
        .def("get_J",                &AlgorithmSelector::get_J_python,         DocSolver::chooseSolver_get_J_python.c_str())
        .def("get_error",            &AlgorithmSelector::get_error,            DocSolver::get_V.c_str())
        .def("get_nb_iter",          &AlgorithmSelector::get_nb_iter,          DocSolver::get_nb_iter.c_str())
        .def("converged",            &AlgorithmSelector::converged,            DocSolver::converged.c_str())
        .def("get_computation_time", &AlgorithmSelector::get_computation_time, DocSolver::get_computation_time.c_str())
        .def("get_timers",           &AlgorithmSelector::get_timers,           "TODO")
        .def("get_timers_jacobian",  &AlgorithmSelector::get_timers_jacobian,  "TODO")
        .def("get_timers_ptdf_lodf", &AlgorithmSelector::get_timers_ptdf_lodf, "TODO")
        .def("get_fdpf_xb_lu",       &AlgorithmSelector::get_fdpf_xb_lu,  py::return_value_policy::reference, DocLSGrid::_internal_do_not_use.c_str())
        .def("get_fdpf_bx_lu",       &AlgorithmSelector::get_fdpf_bx_lu,  py::return_value_policy::reference, DocLSGrid::_internal_do_not_use.c_str());
}
