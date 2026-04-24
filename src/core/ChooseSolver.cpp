// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ChooseSolver.hpp"
#include "SolverTypeNames.hpp"  // solver_type_to_name / name_to_solver_type

namespace ls2g {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

ChooseSolver::ChooseSolver()
    : _solver_type(SolverType::SparseLU),
      _type_used_for_nr(SolverType::SparseLU),
      _gridmodel_ptr(nullptr)
{
    _solver = SolverRegistry::instance().make("SparseLU");
}

// ---------------------------------------------------------------------------
// change_solver overloads
// ---------------------------------------------------------------------------

void ChooseSolver::change_solver(const SolverType& type)
{
    if (type == SolverType::Custom) {
        throw std::runtime_error(
            "ChooseSolver::change_solver: SolverType::Custom is not a concrete solver; "
            "use the string-based change_solver(name) overload instead.");
    }
    change_solver(solver_type_to_name(type));
}

void ChooseSolver::change_solver(const std::string& name)
{
    SolverType type = name_to_solver_type(name);   // SolverType::Custom if plugin

    if (type == _solver_type && type != SolverType::Custom) return;

    #ifndef KLU_SOLVER_AVAILABLE
        if (type == SolverType::KLU ||
            type == SolverType::KLUDC ||
            type == SolverType::KLUSingleSlack ||
            type == SolverType::FDPF_XB_KLU ||
            type == SolverType::FDPF_BX_KLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using KLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef NICSLU_SOLVER_AVAILABLE
        if (type == SolverType::NICSLU ||
            type == SolverType::NICSLUDC ||
            type == SolverType::NICSLUSingleSlack ||
            type == SolverType::FDPF_XB_NICSLU ||
            type == SolverType::FDPF_BX_NICSLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using NICSLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef CKTSO_SOLVER_AVAILABLE
        if (type == SolverType::CKTSO ||
            type == SolverType::CKTSODC ||
            type == SolverType::CKTSOSingleSlack ||
            type == SolverType::FDPF_XB_CKTSO ||
            type == SolverType::FDPF_BX_CKTSO) {
            throw std::runtime_error(
                "Impossible to change for a solver using CKTSO for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    std::unique_ptr<BaseAlgo> new_solver = SolverRegistry::instance().make(name);

    _solver = std::move(new_solver);
    _solver_type = type;
    _type_used_for_nr = type;

    if (_gridmodel_ptr) _solver->set_gridmodel(_gridmodel_ptr);
    _solver->reset();
}

// ---------------------------------------------------------------------------
// operator<<
// ---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const SolverType& solver_type)
{
    switch (solver_type)
    {
    case SolverType::SparseLU:            out << "SparseLU";            break;
    case SolverType::KLU:                 out << "KLU";                 break;
    case SolverType::GaussSeidel:         out << "GaussSeidel";         break;
    case SolverType::DC:                  out << "DC";                  break;
    case SolverType::GaussSeidelSynch:    out << "GaussSeidelSynch";    break;
    case SolverType::NICSLU:              out << "NICSLU";              break;
    case SolverType::SparseLUSingleSlack: out << "SparseLUSingleSlack"; break;
    case SolverType::KLUSingleSlack:      out << "KLUSingleSlack";      break;
    case SolverType::NICSLUSingleSlack:   out << "NICSLUSingleSlack";   break;
    case SolverType::KLUDC:               out << "KLUDC";               break;
    case SolverType::NICSLUDC:            out << "NICSLUDC";            break;
    case SolverType::CKTSO:               out << "CKTSO";               break;
    case SolverType::CKTSOSingleSlack:    out << "CKTSOSingleSlack";    break;
    case SolverType::CKTSODC:             out << "CKTSODC";             break;
    case SolverType::FDPF_XB_SparseLU:   out << "FDPF_XB_SparseLU";   break;
    case SolverType::FDPF_BX_SparseLU:   out << "FDPF_BX_SparseLU";   break;
    case SolverType::FDPF_XB_KLU:        out << "FDPF_XB_KLU";        break;
    case SolverType::FDPF_BX_KLU:        out << "FDPF_BX_KLU";        break;
    case SolverType::FDPF_XB_NICSLU:     out << "FDPF_XB_NICSLU";     break;
    case SolverType::FDPF_BX_NICSLU:     out << "FDPF_BX_NICSLU";     break;
    case SolverType::FDPF_XB_CKTSO:      out << "FDPF_XB_CKTSO";      break;
    case SolverType::FDPF_BX_CKTSO:      out << "FDPF_BX_CKTSO";      break;
    case SolverType::Custom:              out << "Custom";              break;
    default:                              out << "(unknown)";           break;
    }
    return out;
}

} // namespace ls2g
