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
// FDPF SparseLU accessors (defined here to avoid implicit BaseFDPFAlgo
// instantiation via dynamic_cast in every TU that includes AlgorithmSelector.hpp)
// ---------------------------------------------------------------------------

FDPF_XB_SparseLU& AlgorithmSelector::get_fdpf_xb_lu() {
    FDPF_XB_SparseLU* p = dynamic_cast<FDPF_XB_SparseLU*>(_algo.get());
    if (!p) throw std::runtime_error("AlgorithmSelector::get_fdpf_xb_lu: current solver is not FDPF_XB_SparseLU");
    return *p;
}

FDPF_BX_SparseLU& AlgorithmSelector::get_fdpf_bx_lu() {
    FDPF_BX_SparseLU* p = dynamic_cast<FDPF_BX_SparseLU*>(_algo.get());
    if (!p) throw std::runtime_error("AlgorithmSelector::get_fdpf_bx_lu: current solver is not FDPF_BX_SparseLU");
    return *p;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

AlgorithmSelector::AlgorithmSelector()
    : _algo_type(AlgorithmType::NR_SparseLU),
      _algo_type_used_for_nr(AlgorithmType::NR_SparseLU),
      _gridmodel_ptr(nullptr)
{
    _algo = AlgorithmRegistry::instance().make("NR_SparseLU");
}

// ---------------------------------------------------------------------------
// change_algorithm overloads
// ---------------------------------------------------------------------------

void AlgorithmSelector::change_algorithm(const AlgorithmType& type)
{
    if (type == AlgorithmType::Custom) {
        throw std::runtime_error(
            "AlgorithmSelector::change_algorithm: AlgorithmType::Custom is not a concrete solver; "
            "use the string-based change_algorithm(name) overload instead.");
    }
    change_algorithm(algo_type_to_name(type));
}

void AlgorithmSelector::change_algorithm(const std::string& name)
{
    AlgorithmType type = name_to_algo_type(name);   // AlgorithmType::Custom if plugin

    if (type == _algo_type && type != AlgorithmType::Custom) return;

    #ifndef KLU_SOLVER_AVAILABLE
        if (type == AlgorithmType::NR_KLU ||
            type == AlgorithmType::DC_KLU ||
            type == AlgorithmType::NRSing_KLU ||
            type == AlgorithmType::FDPF_XB_KLU ||
            type == AlgorithmType::FDPF_BX_KLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using KLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef NICSLU_SOLVER_AVAILABLE
        if (type == AlgorithmType::NR_NICSLU ||
            type == AlgorithmType::DC_NICSLU ||
            type == AlgorithmType::NRSing_NICSLU ||
            type == AlgorithmType::FDPF_XB_NICSLU ||
            type == AlgorithmType::FDPF_BX_NICSLU) {
            throw std::runtime_error(
                "Impossible to change for a solver using NICSLU for linear algebra. "
                "Please compile lightsim2grid from source to benefit from this.");
        }
    #endif

    #ifndef CKTSO_SOLVER_AVAILABLE
        if (type == AlgorithmType::NR_CKTSO ||
            type == AlgorithmType::DC_CKTSO ||
            type == AlgorithmType::NRSing_CKTSO ||
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

    if (_gridmodel_ptr) _algo->set_lsgrid(_gridmodel_ptr);
    _algo->reset();
}

// ---------------------------------------------------------------------------
// operator<<
// ---------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const AlgorithmType& algo_type)
{
    switch (algo_type)
    {
    case AlgorithmType::NR_SparseLU:      out << "NR_SparseLU";      break;
    case AlgorithmType::NR_KLU:           out << "NR_KLU";           break;
    case AlgorithmType::GaussSeidel:      out << "GaussSeidel";      break;
    case AlgorithmType::DC_SparseLU:      out << "DC_SparseLU";      break;
    case AlgorithmType::GaussSeidelSynch: out << "GaussSeidelSynch"; break;
    case AlgorithmType::NR_NICSLU:        out << "NR_NICSLU";        break;
    case AlgorithmType::NRSing_SparseLU:  out << "NRSing_SparseLU";  break;
    case AlgorithmType::NRSing_KLU:       out << "NRSing_KLU";       break;
    case AlgorithmType::NRSing_NICSLU:    out << "NRSing_NICSLU";    break;
    case AlgorithmType::DC_KLU:           out << "DC_KLU";           break;
    case AlgorithmType::DC_NICSLU:        out << "DC_NICSLU";        break;
    case AlgorithmType::NR_CKTSO:         out << "NR_CKTSO";         break;
    case AlgorithmType::NRSing_CKTSO:     out << "NRSing_CKTSO";     break;
    case AlgorithmType::DC_CKTSO:         out << "DC_CKTSO";         break;
    case AlgorithmType::FDPF_XB_SparseLU: out << "FDPF_XB_SparseLU"; break;
    case AlgorithmType::FDPF_BX_SparseLU: out << "FDPF_BX_SparseLU"; break;
    case AlgorithmType::FDPF_XB_KLU:      out << "FDPF_XB_KLU";      break;
    case AlgorithmType::FDPF_BX_KLU:      out << "FDPF_BX_KLU";      break;
    case AlgorithmType::FDPF_XB_NICSLU:   out << "FDPF_XB_NICSLU";   break;
    case AlgorithmType::FDPF_BX_NICSLU:   out << "FDPF_BX_NICSLU";   break;
    case AlgorithmType::FDPF_XB_CKTSO:    out << "FDPF_XB_CKTSO";    break;
    case AlgorithmType::FDPF_BX_CKTSO:    out << "FDPF_BX_CKTSO";    break;
    case AlgorithmType::Custom:            out << "Custom";            break;
    default:                                 out << "(unknown)";           break;
    }
    return out;
}

} // namespace ls2g
