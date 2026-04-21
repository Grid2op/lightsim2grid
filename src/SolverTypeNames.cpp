// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ChooseSolver.hpp"   // defines SolverType enum
#include "SolverTypeNames.hpp"
#include <stdexcept>

std::string solver_type_to_name(SolverType t) {
    switch (t) {
        case SolverType::SparseLU:          return "SparseLU";
        case SolverType::SparseLUSingleSlack: return "SparseLUSingleSlack";
        case SolverType::DC:                return "DC";
        case SolverType::GaussSeidel:       return "GaussSeidel";
        case SolverType::GaussSeidelSynch:  return "GaussSeidelSynch";
        case SolverType::FDPF_XB_SparseLU:  return "FDPF_XB_SparseLU";
        case SolverType::FDPF_BX_SparseLU:  return "FDPF_BX_SparseLU";
#ifdef KLU_SOLVER_AVAILABLE
        case SolverType::KLU:               return "KLU";
        case SolverType::KLUSingleSlack:    return "KLUSingleSlack";
        case SolverType::KLUDC:             return "KLUDC";
        case SolverType::FDPF_XB_KLU:      return "FDPF_XB_KLU";
        case SolverType::FDPF_BX_KLU:      return "FDPF_BX_KLU";
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
        case SolverType::NICSLU:            return "NICSLU";
        case SolverType::NICSLUSingleSlack: return "NICSLUSingleSlack";
        case SolverType::NICSLUDC:          return "NICSLUDC";
        case SolverType::FDPF_XB_NICSLU:   return "FDPF_XB_NICSLU";
        case SolverType::FDPF_BX_NICSLU:   return "FDPF_BX_NICSLU";
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
        case SolverType::CKTSO:             return "CKTSO";
        case SolverType::CKTSOSingleSlack:  return "CKTSOSingleSlack";
        case SolverType::CKTSODC:           return "CKTSODC";
        case SolverType::FDPF_XB_CKTSO:    return "FDPF_XB_CKTSO";
        case SolverType::FDPF_BX_CKTSO:    return "FDPF_BX_CKTSO";
#endif
        case SolverType::Custom:
            throw std::runtime_error(
                "solver_type_to_name: SolverType::Custom has no canonical string name; "
                "use the string-based API instead.");
        default:
            throw std::runtime_error("solver_type_to_name: unknown SolverType enum value");
    }
}

SolverType name_to_solver_type(const std::string& name) {
    if (name == "SparseLU")           return SolverType::SparseLU;
    if (name == "SparseLUSingleSlack") return SolverType::SparseLUSingleSlack;
    if (name == "DC")                 return SolverType::DC;
    if (name == "GaussSeidel")        return SolverType::GaussSeidel;
    if (name == "GaussSeidelSynch")   return SolverType::GaussSeidelSynch;
    if (name == "FDPF_XB_SparseLU")  return SolverType::FDPF_XB_SparseLU;
    if (name == "FDPF_BX_SparseLU")  return SolverType::FDPF_BX_SparseLU;
#ifdef KLU_SOLVER_AVAILABLE
    if (name == "KLU")               return SolverType::KLU;
    if (name == "KLUSingleSlack")    return SolverType::KLUSingleSlack;
    if (name == "KLUDC")             return SolverType::KLUDC;
    if (name == "FDPF_XB_KLU")      return SolverType::FDPF_XB_KLU;
    if (name == "FDPF_BX_KLU")      return SolverType::FDPF_BX_KLU;
#endif
#ifdef NICSLU_SOLVER_AVAILABLE
    if (name == "NICSLU")            return SolverType::NICSLU;
    if (name == "NICSLUSingleSlack") return SolverType::NICSLUSingleSlack;
    if (name == "NICSLUDC")          return SolverType::NICSLUDC;
    if (name == "FDPF_XB_NICSLU")   return SolverType::FDPF_XB_NICSLU;
    if (name == "FDPF_BX_NICSLU")   return SolverType::FDPF_BX_NICSLU;
#endif
#ifdef CKTSO_SOLVER_AVAILABLE
    if (name == "CKTSO")             return SolverType::CKTSO;
    if (name == "CKTSOSingleSlack")  return SolverType::CKTSOSingleSlack;
    if (name == "CKTSODC")           return SolverType::CKTSODC;
    if (name == "FDPF_XB_CKTSO")    return SolverType::FDPF_XB_CKTSO;
    if (name == "FDPF_BX_CKTSO")    return SolverType::FDPF_BX_CKTSO;
#endif
    // Unknown names are plugin / external solvers
    return SolverType::Custom;
}
