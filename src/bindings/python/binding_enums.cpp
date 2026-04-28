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

    py::enum_<AlgorithmType>(m, "AlgorithmType", "This enum controls the solver you want to use.")
        .value("GaussSeidel", AlgorithmType::GaussSeidel, "denotes the :class:`lightsim2grid.solver.GaussSeidelSolver`")
        .value("GaussSeidelSynch", AlgorithmType::GaussSeidelSynch, "denotes the :class:`lightsim2grid.solver.GaussSeidelSynchSolver`")
        .value("SparseLU", AlgorithmType::SparseLU, "denotes the :class:`lightsim2grid.solver.SparseLUSolver`")
        .value("SparseLUSingleSlack", AlgorithmType::SparseLUSingleSlack, "denotes the :class:`lightsim2grid.solver.SparseLUSolverSingleSlack`")
        .value("DC", AlgorithmType::DC, "denotes the :class:`lightsim2grid.solver.DCSolver`")
        .value("KLU", AlgorithmType::KLU, "denotes the :class:`lightsim2grid.solver.KLUSolver`")
        .value("KLUSingleSlack", AlgorithmType::KLUSingleSlack, "denotes the :class:`lightsim2grid.solver.KLUSolverSingleSlack`")
        .value("KLUDC", AlgorithmType::KLUDC, "denotes the :class:`lightsim2grid.solver.KLUDCSolver`")
        .value("NICSLU", AlgorithmType::NICSLU, "denotes the :class:`lightsim2grid.solver.NICSLUSolver`")
        .value("NICSLUSingleSlack", AlgorithmType::NICSLUSingleSlack, "denotes the :class:`lightsim2grid.solver.NICSLUSolverSingleSlack`")
        .value("NICSLUDC", AlgorithmType::NICSLUDC, "denotes the :class:`lightsim2grid.solver.NICSLUDCSolver`")
        .value("CKTSO", AlgorithmType::CKTSO, "denotes the :class:`lightsim2grid.solver.CKTSOSolver`")
        .value("CKTSOSingleSlack", AlgorithmType::CKTSOSingleSlack, "denotes the :class:`lightsim2grid.solver.CKTSOSolverSingleSlack`")
        .value("CKTSODC", AlgorithmType::CKTSODC, "denotes the :class:`lightsim2grid.solver.CKTSODCSolver`")
        .value("FDPF_XB_SparseLU", AlgorithmType::FDPF_XB_SparseLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_SparseLUSolver`")
        .value("FDPF_BX_SparseLU", AlgorithmType::FDPF_BX_SparseLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_SparseLUSolver`")
        .value("FDPF_XB_KLU", AlgorithmType::FDPF_XB_KLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_KLUSolver`")
        .value("FDPF_BX_KLU", AlgorithmType::FDPF_BX_KLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_KLUSolver`")
        .value("FDPF_XB_NICSLU", AlgorithmType::FDPF_XB_NICSLU, "denotes the :class:`lightsim2grid.solver.FDPF_XB_NICSLUSolver`")
        .value("FDPF_BX_NICSLU", AlgorithmType::FDPF_BX_NICSLU, "denotes the :class:`lightsim2grid.solver.FDPF_BX_NICSLUSolver`")
        .value("FDPF_XB_CKTSO", AlgorithmType::FDPF_XB_CKTSO, "denotes the :class:`lightsim2grid.solver.FDPF_XB_CKTSOSolver`")
        .value("FDPF_BX_CKTSO", AlgorithmType::FDPF_BX_CKTSO, "denotes the :class:`lightsim2grid.solver.FDPF_BX_CKTSOSolver`")
        .value("Custom", AlgorithmType::Custom, "sentinel value for external/plugin solvers loaded via load_solver_plugin()")
        .export_values();

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
