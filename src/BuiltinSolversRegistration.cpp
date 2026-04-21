// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BuiltinSolversRegistration.hpp"
#include "SolverRegistry.hpp"
#include "Solvers.hpp"   // all concrete solver typedefs + #ifdef guards

namespace ls2g {

void register_builtin_solvers(SolverRegistry& reg) {
    reg.register_solver("SparseLU",
        []{ return std::unique_ptr<BaseAlgo>(new SparseLUSolver()); });
    reg.register_solver("SparseLUSingleSlack",
        []{ return std::unique_ptr<BaseAlgo>(new SparseLUSolverSingleSlack()); });
    reg.register_solver("DC",
        []{ return std::unique_ptr<BaseAlgo>(new DCSolver()); });
    reg.register_solver("GaussSeidel",
        []{ return std::unique_ptr<BaseAlgo>(new GaussSeidelAlgo()); });
    reg.register_solver("GaussSeidelSynch",
        []{ return std::unique_ptr<BaseAlgo>(new GaussSeidelSynchAlgo()); });
    reg.register_solver("FDPF_XB_SparseLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_XB_SparseLUSolver()); });
    reg.register_solver("FDPF_BX_SparseLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_BX_SparseLUSolver()); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.register_solver("KLU",
        []{ return std::unique_ptr<BaseAlgo>(new KLUSolver()); });
    reg.register_solver("KLUSingleSlack",
        []{ return std::unique_ptr<BaseAlgo>(new KLUSolverSingleSlack()); });
    reg.register_solver("KLUDC",
        []{ return std::unique_ptr<BaseAlgo>(new KLUDCSolver()); });
    reg.register_solver("FDPF_XB_KLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_XB_KLUSolver()); });
    reg.register_solver("FDPF_BX_KLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_BX_KLUSolver()); });
#endif // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    reg.register_solver("NICSLU",
        []{ return std::unique_ptr<BaseAlgo>(new NICSLUSolver()); });
    reg.register_solver("NICSLUSingleSlack",
        []{ return std::unique_ptr<BaseAlgo>(new NICSLUSolverSingleSlack()); });
    reg.register_solver("NICSLUDC",
        []{ return std::unique_ptr<BaseAlgo>(new NICSLUDCSolver()); });
    reg.register_solver("FDPF_XB_NICSLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_XB_NICSLUSolver()); });
    reg.register_solver("FDPF_BX_NICSLU",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_BX_NICSLUSolver()); });
#endif // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    reg.register_solver("CKTSO",
        []{ return std::unique_ptr<BaseAlgo>(new CKTSOSolver()); });
    reg.register_solver("CKTSOSingleSlack",
        []{ return std::unique_ptr<BaseAlgo>(new CKTSOSolverSingleSlack()); });
    reg.register_solver("CKTSODC",
        []{ return std::unique_ptr<BaseAlgo>(new CKTSODCSolver()); });
    reg.register_solver("FDPF_XB_CKTSO",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_XB_CKTSOSolver()); });
    reg.register_solver("FDPF_BX_CKTSO",
        []{ return std::unique_ptr<BaseAlgo>(new FDPF_BX_CKTSOSolver()); });
#endif // CKTSO_SOLVER_AVAILABLE
}

} // namespace ls2g
