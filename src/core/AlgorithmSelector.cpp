// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "AlgorithmSelector.hpp"
#include "AlgorithmTypeNames.hpp"  // algo_type_to_name / name_to_algo_type

namespace ls2g {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

AlgorithmSelector::AlgorithmSelector()
    : _algo_type(AlgorithmType::SparseLU),
      _algo_type_used_for_nr(AlgorithmType::SparseLU),
      _gridmodel_ptr(nullptr)
{
    _algo = AlgorithmRegistry::instance().make("SparseLU");
}

// ---------------------------------------------------------------------------
// change_solver overloads
// ---------------------------------------------------------------------------

void AlgorithmSelector::change_solver(const AlgorithmType& type)
{
    if (type == AlgorithmType::Custom) {
        throw std::runtime_error(
            "AlgorithmSelector::change_solver: AlgorithmType::Custom is not a concrete solver; "
            "use the string-based change_solver(name) overload instead.");
    }
    change_solver(algo_type_to_name(type));
}

void AlgorithmSelector::change_solver(const std::string& name)
{
    AlgorithmType type = name_to_algo_type(name);   // AlgorithmType::Custom if plugin

    if (type == _algo_type && type != AlgorithmType::Custom) return;

    #ifndef KLU_SOLVER_AVAILABLE
        if (type == AlgorithmType::KLU ||
            type == AlgorithmType::KLUDC ||
            type == AlgorithmType::KLUSingleSlack ||
            type == AlgorithmType::FDPF_XB_KLU ||
            type == AlgorithmType::FDPF_BX_KLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using KLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef NICSLU_SOLVER_AVAILABLE
        if (type == AlgorithmType::NICSLU ||
            type == AlgorithmType::NICSLUDC ||
            type == AlgorithmType::NICSLUSingleSlack ||
            type == AlgorithmType::FDPF_XB_NICSLU ||
            type == AlgorithmType::FDPF_BX_NICSLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using NICSLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef CKTSO_SOLVER_AVAILABLE
        if (type == AlgorithmType::CKTSO ||
            type == AlgorithmType::CKTSODC ||
            type == AlgorithmType::CKTSOSingleSlack ||
            type == AlgorithmType::FDPF_XB_CKTSO ||
            type == AlgorithmType::FDPF_BX_CKTSO) {
            throw std::runtime_error(
                "Impossible to change for a solver using CKTSO for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    std::unique_ptr<BaseAlgo> new_algo = AlgorithmRegistry::instance().make(name);

    _algo = std::move(new_algo);
    _algo_type = type;
    _algo_type_used_for_nr = type;

    if (_gridmodel_ptr) _algo->set_gridmodel(_gridmodel_ptr);
    _algo->reset();
}

// ---------------------------------------------------------------------------
// operator<<
// ---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const AlgorithmType& algo_type)
{
    switch (algo_type)
    {
    case AlgorithmType::SparseLU:            out << "SparseLU";            break;
    case AlgorithmType::KLU:                 out << "KLU";                 break;
    case AlgorithmType::GaussSeidel:         out << "GaussSeidel";         break;
    case AlgorithmType::DC:                  out << "DC";                  break;
    case AlgorithmType::GaussSeidelSynch:    out << "GaussSeidelSynch";    break;
    case AlgorithmType::NICSLU:              out << "NICSLU";              break;
    case AlgorithmType::SparseLUSingleSlack: out << "SparseLUSingleSlack"; break;
    case AlgorithmType::KLUSingleSlack:      out << "KLUSingleSlack";      break;
    case AlgorithmType::NICSLUSingleSlack:   out << "NICSLUSingleSlack";   break;
    case AlgorithmType::KLUDC:               out << "KLUDC";               break;
    case AlgorithmType::NICSLUDC:            out << "NICSLUDC";            break;
    case AlgorithmType::CKTSO:               out << "CKTSO";               break;
    case AlgorithmType::CKTSOSingleSlack:    out << "CKTSOSingleSlack";    break;
    case AlgorithmType::CKTSODC:             out << "CKTSODC";             break;
    case AlgorithmType::FDPF_XB_SparseLU:   out << "FDPF_XB_SparseLU";   break;
    case AlgorithmType::FDPF_BX_SparseLU:   out << "FDPF_BX_SparseLU";   break;
    case AlgorithmType::FDPF_XB_KLU:        out << "FDPF_XB_KLU";        break;
    case AlgorithmType::FDPF_BX_KLU:        out << "FDPF_BX_KLU";        break;
    case AlgorithmType::FDPF_XB_NICSLU:     out << "FDPF_XB_NICSLU";     break;
    case AlgorithmType::FDPF_BX_NICSLU:     out << "FDPF_BX_NICSLU";     break;
    case AlgorithmType::FDPF_XB_CKTSO:      out << "FDPF_XB_CKTSO";      break;
    case AlgorithmType::FDPF_BX_CKTSO:      out << "FDPF_BX_CKTSO";      break;
    case AlgorithmType::Custom:              out << "Custom";              break;
    default:                                 out << "(unknown)";           break;
    }
    return out;
}

} // namespace ls2g
