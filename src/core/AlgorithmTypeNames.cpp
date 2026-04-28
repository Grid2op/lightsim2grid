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
        case AlgorithmType::SparseLU:          return "SparseLU";
        case AlgorithmType::SparseLUSingleSlack: return "SparseLUSingleSlack";
        case AlgorithmType::DC:                return "DC";
        case AlgorithmType::GaussSeidel:       return "GaussSeidel";
        case AlgorithmType::GaussSeidelSynch:  return "GaussSeidelSynch";
        case AlgorithmType::FDPF_XB_SparseLU:  return "FDPF_XB_SparseLU";
        case AlgorithmType::FDPF_BX_SparseLU:  return "FDPF_BX_SparseLU";
#ifdef KLU_SOLVER_AVAILABLE
        case AlgorithmType::KLU:               return "KLU";
        case AlgorithmType::KLUSingleSlack:    return "KLUSingleSlack";
        case AlgorithmType::KLUDC:             return "KLUDC";
        case AlgorithmType::FDPF_XB_KLU:      return "FDPF_XB_KLU";
        case AlgorithmType::FDPF_BX_KLU:      return "FDPF_BX_KLU";
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
        case AlgorithmType::NICSLU:            return "NICSLU";
        case AlgorithmType::NICSLUSingleSlack: return "NICSLUSingleSlack";
        case AlgorithmType::NICSLUDC:          return "NICSLUDC";
        case AlgorithmType::FDPF_XB_NICSLU:   return "FDPF_XB_NICSLU";
        case AlgorithmType::FDPF_BX_NICSLU:   return "FDPF_BX_NICSLU";
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
        case AlgorithmType::CKTSO:             return "CKTSO";
        case AlgorithmType::CKTSOSingleSlack:  return "CKTSOSingleSlack";
        case AlgorithmType::CKTSODC:           return "CKTSODC";
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
    if (name == "SparseLU")           return AlgorithmType::SparseLU;
    if (name == "SparseLUSingleSlack") return AlgorithmType::SparseLUSingleSlack;
    if (name == "DC")                 return AlgorithmType::DC;
    if (name == "GaussSeidel")        return AlgorithmType::GaussSeidel;
    if (name == "GaussSeidelSynch")   return AlgorithmType::GaussSeidelSynch;
    if (name == "FDPF_XB_SparseLU")  return AlgorithmType::FDPF_XB_SparseLU;
    if (name == "FDPF_BX_SparseLU")  return AlgorithmType::FDPF_BX_SparseLU;
#ifdef KLU_SOLVER_AVAILABLE
    if (name == "KLU")               return AlgorithmType::KLU;
    if (name == "KLUSingleSlack")    return AlgorithmType::KLUSingleSlack;
    if (name == "KLUDC")             return AlgorithmType::KLUDC;
    if (name == "FDPF_XB_KLU")      return AlgorithmType::FDPF_XB_KLU;
    if (name == "FDPF_BX_KLU")      return AlgorithmType::FDPF_BX_KLU;
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
    if (name == "NICSLU")            return AlgorithmType::NICSLU;
    if (name == "NICSLUSingleSlack") return AlgorithmType::NICSLUSingleSlack;
    if (name == "NICSLUDC")          return AlgorithmType::NICSLUDC;
    if (name == "FDPF_XB_NICSLU")   return AlgorithmType::FDPF_XB_NICSLU;
    if (name == "FDPF_BX_NICSLU")   return AlgorithmType::FDPF_BX_NICSLU;
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
    if (name == "CKTSO")             return AlgorithmType::CKTSO;
    if (name == "CKTSOSingleSlack")  return AlgorithmType::CKTSOSingleSlack;
    if (name == "CKTSODC")           return AlgorithmType::CKTSODC;
    if (name == "FDPF_XB_CKTSO")    return AlgorithmType::FDPF_XB_CKTSO;
    if (name == "FDPF_BX_CKTSO")    return AlgorithmType::FDPF_BX_CKTSO;
#endif
    // Unknown names are plugin / external solvers
    return AlgorithmType::Custom;
}

} // namespace ls2g
