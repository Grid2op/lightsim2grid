// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "binding_declarations.hpp"
#include "BaseConstants.hpp"
#include "AlgorithmSelector.hpp"
#include "Utils.hpp"
#include "powerflow_algorithm/ScalingPolicies.hpp"
#include "powerflow_algorithm/RefactorPolicies.hpp"

using namespace ls2g;

void bind_enums(py::module_& m) {
    py::enum_<FDPFMethod>(m, "FDPFMethod", "This enum controls the type of method you can use for Fast Decoupled Powerflow (XB or BX)")
        .value("XB", FDPFMethod::XB, "denotes the XB method")
        .value("BX", FDPFMethod::BX, "denotes the BX method")
        .export_values();

    py::enum_<AlgorithmType>(m, "AlgorithmType", "This enum controls the powerflow algorithm you want to use.")
        // ---- GaussSeidel (no linear solver choice) ----
        .value("GaussSeidel",      AlgorithmType::GaussSeidel,      "denotes the :class:`lightsim2grid.solver.GaussSeidelAlgo`")
        .value("GaussSeidelSynch", AlgorithmType::GaussSeidelSynch, "denotes the :class:`lightsim2grid.solver.GaussSeidelSynchAlgo`")
        // ---- SparseLU family ----
        .value("NR_SparseLU",     AlgorithmType::NR_SparseLU,     "Newton-Raphson (multi-slack) + SparseLU; see :class:`lightsim2grid.solver.NR_SparseLU`")
        .value("NRSing_SparseLU", AlgorithmType::NRSing_SparseLU, "Newton-Raphson (single-slack) + SparseLU; see :class:`lightsim2grid.solver.NRSing_SparseLU`")
        .value("DC_SparseLU",     AlgorithmType::DC_SparseLU,     "DC approximation + SparseLU; see :class:`lightsim2grid.solver.DC_SparseLU`")
        .value("FDPF_XB_SparseLU", AlgorithmType::FDPF_XB_SparseLU, "Fast-Decoupled PF (XB) + SparseLU; see :class:`lightsim2grid.solver.FDPF_XB_SparseLU`")
        .value("FDPF_BX_SparseLU", AlgorithmType::FDPF_BX_SparseLU, "Fast-Decoupled PF (BX) + SparseLU; see :class:`lightsim2grid.solver.FDPF_BX_SparseLU`")
        // ---- KLU family ----
        .value("NR_KLU",          AlgorithmType::NR_KLU,     "Newton-Raphson (multi-slack) + KLU; see :class:`lightsim2grid.solver.NR_KLU`")
        .value("NRSing_KLU",      AlgorithmType::NRSing_KLU, "Newton-Raphson (single-slack) + KLU; see :class:`lightsim2grid.solver.NRSing_KLU`")
        .value("DC_KLU",          AlgorithmType::DC_KLU,     "DC approximation + KLU; see :class:`lightsim2grid.solver.DC_KLU`")
        .value("FDPF_XB_KLU",     AlgorithmType::FDPF_XB_KLU, "Fast-Decoupled PF (XB) + KLU; see :class:`lightsim2grid.solver.FDPF_XB_KLU`")
        .value("FDPF_BX_KLU",     AlgorithmType::FDPF_BX_KLU, "Fast-Decoupled PF (BX) + KLU; see :class:`lightsim2grid.solver.FDPF_BX_KLU`")
        // ---- NICSLU family ----
        .value("NR_NICSLU",       AlgorithmType::NR_NICSLU,     "Newton-Raphson (multi-slack) + NICSLU; see :class:`lightsim2grid.solver.NR_NICSLU`")
        .value("NRSing_NICSLU",   AlgorithmType::NRSing_NICSLU, "Newton-Raphson (single-slack) + NICSLU; see :class:`lightsim2grid.solver.NRSing_NICSLU`")
        .value("DC_NICSLU",       AlgorithmType::DC_NICSLU,     "DC approximation + NICSLU; see :class:`lightsim2grid.solver.DC_NICSLU`")
        .value("FDPF_XB_NICSLU",  AlgorithmType::FDPF_XB_NICSLU, "Fast-Decoupled PF (XB) + NICSLU; see :class:`lightsim2grid.solver.FDPF_XB_NICSLU`")
        .value("FDPF_BX_NICSLU",  AlgorithmType::FDPF_BX_NICSLU, "Fast-Decoupled PF (BX) + NICSLU; see :class:`lightsim2grid.solver.FDPF_BX_NICSLU`")
        // ---- CKTSO family ----
        .value("NR_CKTSO",        AlgorithmType::NR_CKTSO,     "Newton-Raphson (multi-slack) + CKTSO; see :class:`lightsim2grid.solver.NR_CKTSO`")
        .value("NRSing_CKTSO",    AlgorithmType::NRSing_CKTSO, "Newton-Raphson (single-slack) + CKTSO; see :class:`lightsim2grid.solver.NRSing_CKTSO`")
        .value("DC_CKTSO",        AlgorithmType::DC_CKTSO,     "DC approximation + CKTSO; see :class:`lightsim2grid.solver.DC_CKTSO`")
        .value("FDPF_XB_CKTSO",   AlgorithmType::FDPF_XB_CKTSO, "Fast-Decoupled PF (XB) + CKTSO; see :class:`lightsim2grid.solver.FDPF_XB_CKTSO`")
        .value("FDPF_BX_CKTSO",   AlgorithmType::FDPF_BX_CKTSO, "Fast-Decoupled PF (BX) + CKTSO; see :class:`lightsim2grid.solver.FDPF_BX_CKTSO`")
        // ---- Plugin/external ----
        .value("Custom", AlgorithmType::Custom, "sentinel value for external/plugin solvers loaded via load_solver_plugin()")
        // // ---- Deprecated aliases (old names; kept for backward compat — repr() shows canonical new name) ----
        // .value("SparseLU",          AlgorithmType::NR_SparseLU,     "Deprecated: use NR_SparseLU")
        // .value("SparseLUSingleSlack", AlgorithmType::NRSing_SparseLU, "Deprecated: use NRSing_SparseLU")
        // .value("DC",                AlgorithmType::DC_SparseLU,     "Deprecated: use DC_SparseLU")
        // .value("KLU",               AlgorithmType::NR_KLU,          "Deprecated: use NR_KLU")
        // .value("KLUSingleSlack",    AlgorithmType::NRSing_KLU,      "Deprecated: use NRSing_KLU")
        // .value("KLUDC",             AlgorithmType::DC_KLU,          "Deprecated: use DC_KLU")
        // .value("NICSLU",            AlgorithmType::NR_NICSLU,       "Deprecated: use NR_NICSLU")
        // .value("NICSLUSingleSlack", AlgorithmType::NRSing_NICSLU,   "Deprecated: use NRSing_NICSLU")
        // .value("NICSLUDC",          AlgorithmType::DC_NICSLU,       "Deprecated: use DC_NICSLU")
        // .value("CKTSO",             AlgorithmType::NR_CKTSO,        "Deprecated: use NR_CKTSO")
        // .value("CKTSOSingleSlack",  AlgorithmType::NRSing_CKTSO,    "Deprecated: use NRSing_CKTSO")
        // .value("CKTSODC",           AlgorithmType::DC_CKTSO,        "Deprecated: use DC_CKTSO")
        ;

    py::enum_<ScalingPolicyType>(m, "ScalingPolicyType", "Step-scaling strategy for the Newton-Raphson loop")
        .value("NoScaling",        ScalingPolicyType::NoScaling,        "Full Newton step (alpha = 1), zero overhead")
        .value("MaxVoltageChange", ScalingPolicyType::MaxVoltageChange, "Clamp step so max|dVa| <= max_dVa and max|dVm| <= max_dVm")
        .value("LineSearch",       ScalingPolicyType::LineSearch,       "Armijo backtracking line search")
        .value("Iwamoto",          ScalingPolicyType::Iwamoto,          "Iwamoto optimal multiplier")
        .export_values();

    py::enum_<RefactorPolicyType>(m, "RefactorPolicyType", "Jacobian refactorization strategy for the Newton-Raphson loop")
        .value("AlwaysRefactor", RefactorPolicyType::AlwaysRefactor, "Rebuild and refactorize J every iteration (default)")
        .value("EveryN",         RefactorPolicyType::EveryN,         "Refactorize every N iterations; update values only in between")
        .value("Chord",          RefactorPolicyType::Chord,          "Build J once on the first iteration; reuse factorization throughout")
        .export_values();

    py::enum_<ErrorType>(m, "ErrorType", "This enum controls the error encountered in the solver")
        .value("NoError", ErrorType::NoError, "No error were encountered")
        .value("SingularMatrix", ErrorType::SingularMatrix, "The Jacobian matrix was singular and could not be factorized (most likely, the grid is not connex)")
        .value("TooManyIterations", ErrorType::TooManyIterations, "The solver reached the maximum number of iterations allowed")
        .value("InifiniteValue", ErrorType::InifiniteValue, "Some infinite values were encountered in the update vector (to update Vm or Va)")
        .value("SolverAnalyze", ErrorType::SolverAnalyze, "The linear solver failed at the 'analyze' step (*eg* `analyzePattern` for Eigen, `klu_analyze` for KLU or `Initialize` for NICSLU")
        .value("SolverFactor", ErrorType::SolverFactor, "The linear solver failed to factor the jacobian matrix (*eg* `factorize` for Eigen (first call), `klu_factor` for KLU or `FactorizeMatrix` for NICSLU (first call)")
        .value("SolverReFactor", ErrorType::SolverReFactor, "The linear solver failed to (re)factor the jacobian matrix (*eg* `factorize` for Eigen (later calls), `klu_refactor` for KLU or `FactorizeMatrix` for NICSLU (later calls)")
        .value("SolverSolve", ErrorType::SolverSolve, "The linear solve failed to solve the linear system J.X = b (*eg* `solve` for Eigen, `klu_solve` for KLU or `Solve` for NICSLU")
        .value("NotInitError", ErrorType::NotInitError, "Attempt to perform some powerflow computation when the linear solver is not initiliazed")
        .value("LicenseError", ErrorType::LicenseError, "Impossible to use the linear solver as the license cannot be found (*eg* unable to locate the `nicslu.lic` file")
        .export_values();
}
