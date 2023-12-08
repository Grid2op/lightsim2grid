// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ChooseSolver.h"

std::ostream& operator<<(std::ostream& out, const SolverType& solver_type)
{
    switch (solver_type)
    {
    case SolverType::SparseLU:
        out << "SparseLU";
        break;
    case SolverType::KLU:
        out << "KLU";
        break;
    case SolverType::GaussSeidel:
        out << "GaussSeidel";
        break;
    case SolverType::DC:
        out << "DC";
        break;
    case SolverType::GaussSeidelSynch:
        out << "GaussSeidelSynch";
        break;
    case SolverType::NICSLU:
        out << "NICSLU";
        break;
    case SolverType::SparseLUSingleSlack:
        out << "SparseLUSingleSlack";
        break;
    case SolverType::KLUSingleSlack:
        out << "KLUSingleSlack";
        break;
    case SolverType::NICSLUSingleSlack:
        out << "NICSLUSingleSlack";
        break;
    case SolverType::KLUDC:
        out << "KLUDC";
        break;
    case SolverType::NICSLUDC:
        out << "NICSLUDC";
        break;
    case SolverType::CKTSO:
        out << "CKTSO";
        break;
    case SolverType::CKTSOSingleSlack:
        out << "CKTSOSingleSlack";
        break;
    case SolverType::CKTSODC:
        out << "CKTSODC";
        break;
    case SolverType::FDPF_XB_SparseLU:
        out << "FDPF_XB_SparseLU";
        break;
    case SolverType::FDPF_BX_SparseLU:
        out << "FDPF_BX_SparseLU";
        break;
    case SolverType::FDPF_XB_KLU:
        out << "FDPF_XB_KLU";
        break;
    case SolverType::FDPF_BX_KLU:
        out << "FDPF_BX_KLU";
        break;
    case SolverType::FDPF_XB_NICSLU:
        out << "FDPF_XB_NICSLU";
        break;
    case SolverType::FDPF_BX_NICSLU:
        out << "FDPF_BX_NICSLU";
        break;
    case SolverType::FDPF_XB_CKTSO:
        out << "FDPF_XB_CKTSO";
        break;
    case SolverType::FDPF_BX_CKTSO:
        out << "FDPF_BX_CKTSO";
        break;
    default:
        out << "(unknown)";
        break;
    }
   return out;
}
