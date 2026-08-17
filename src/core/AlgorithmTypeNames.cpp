// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "AlgorithmSelector.hpp"   // defines AlgorithmType enum
#include "AlgorithmTypeNames.hpp"
#include <stdexcept>

namespace ls2g {

std::string algo_type_to_name(AlgorithmType t) {
    switch (t) {
        case AlgorithmType::NR_SparseLU:      return "NR_SparseLU";
        case AlgorithmType::NRSing_SparseLU:  return "NRSing_SparseLU";
        case AlgorithmType::DC_SparseLU:      return "DC_SparseLU";
        case AlgorithmType::GaussSeidel:      return "GaussSeidel";
        case AlgorithmType::GaussSeidelSynch: return "GaussSeidelSynch";
        case AlgorithmType::FDPF_XB_SparseLU: return "FDPF_XB_SparseLU";
        case AlgorithmType::FDPF_BX_SparseLU: return "FDPF_BX_SparseLU";
#ifdef KLU_SOLVER_AVAILABLE
        case AlgorithmType::NR_KLU:           return "NR_KLU";
        case AlgorithmType::NRSing_KLU:       return "NRSing_KLU";
        case AlgorithmType::DC_KLU:           return "DC_KLU";
        case AlgorithmType::FDPF_XB_KLU:      return "FDPF_XB_KLU";
        case AlgorithmType::FDPF_BX_KLU:      return "FDPF_BX_KLU";
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
        case AlgorithmType::NR_NICSLU:        return "NR_NICSLU";
        case AlgorithmType::NRSing_NICSLU:    return "NRSing_NICSLU";
        case AlgorithmType::DC_NICSLU:        return "DC_NICSLU";
        case AlgorithmType::FDPF_XB_NICSLU:   return "FDPF_XB_NICSLU";
        case AlgorithmType::FDPF_BX_NICSLU:   return "FDPF_BX_NICSLU";
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
        case AlgorithmType::NR_CKTSO:         return "NR_CKTSO";
        case AlgorithmType::NRSing_CKTSO:     return "NRSing_CKTSO";
        case AlgorithmType::DC_CKTSO:         return "DC_CKTSO";
        case AlgorithmType::FDPF_XB_CKTSO:    return "FDPF_XB_CKTSO";
        case AlgorithmType::FDPF_BX_CKTSO:    return "FDPF_BX_CKTSO";
#endif
        case AlgorithmType::Custom:
            throw std::runtime_error(
                "algo_type_to_name: AlgorithmType::Custom has no canonical string name; "
                "use the string-based API instead.");
        default:
            throw std::runtime_error("algo_type_to_name: unknown AlgorithmType enum value");
    }
}

AlgorithmType name_to_algo_type(const std::string& name) {
    if (name == "NR_SparseLU")      return AlgorithmType::NR_SparseLU;
    if (name == "NRSing_SparseLU")  return AlgorithmType::NRSing_SparseLU;
    if (name == "DC_SparseLU")      return AlgorithmType::DC_SparseLU;
    if (name == "GaussSeidel")      return AlgorithmType::GaussSeidel;
    if (name == "GaussSeidelSynch") return AlgorithmType::GaussSeidelSynch;
    if (name == "FDPF_XB_SparseLU") return AlgorithmType::FDPF_XB_SparseLU;
    if (name == "FDPF_BX_SparseLU") return AlgorithmType::FDPF_BX_SparseLU;
#ifdef KLU_SOLVER_AVAILABLE
    if (name == "NR_KLU")           return AlgorithmType::NR_KLU;
    if (name == "NRSing_KLU")       return AlgorithmType::NRSing_KLU;
    if (name == "DC_KLU")           return AlgorithmType::DC_KLU;
    if (name == "FDPF_XB_KLU")      return AlgorithmType::FDPF_XB_KLU;
    if (name == "FDPF_BX_KLU")      return AlgorithmType::FDPF_BX_KLU;
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
    if (name == "NR_NICSLU")        return AlgorithmType::NR_NICSLU;
    if (name == "NRSing_NICSLU")    return AlgorithmType::NRSing_NICSLU;
    if (name == "DC_NICSLU")        return AlgorithmType::DC_NICSLU;
    if (name == "FDPF_XB_NICSLU")   return AlgorithmType::FDPF_XB_NICSLU;
    if (name == "FDPF_BX_NICSLU")   return AlgorithmType::FDPF_BX_NICSLU;
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
    if (name == "NR_CKTSO")         return AlgorithmType::NR_CKTSO;
    if (name == "NRSing_CKTSO")     return AlgorithmType::NRSing_CKTSO;
    if (name == "DC_CKTSO")         return AlgorithmType::DC_CKTSO;
    if (name == "FDPF_XB_CKTSO")    return AlgorithmType::FDPF_XB_CKTSO;
    if (name == "FDPF_BX_CKTSO")    return AlgorithmType::FDPF_BX_CKTSO;
#endif
    // Unknown names are plugin / external solvers
    return AlgorithmType::Custom;
}

} // namespace ls2g
