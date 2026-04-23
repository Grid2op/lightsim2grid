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
        []{ return std::make_unique<SparseLUSolver>(); });
    reg.register_solver("SparseLUSingleSlack",
        []{ return std::make_unique<SparseLUSolverSingleSlack>(); });
    reg.register_solver("DC",
        []{ return std::make_unique<DCSolver>(); });
    reg.register_solver("GaussSeidel",
        []{ return std::make_unique<GaussSeidelAlgo>(); });
    reg.register_solver("GaussSeidelSynch",
        []{ return std::make_unique<GaussSeidelSynchAlgo>(); });
    reg.register_solver("FDPF_XB_SparseLU",
        []{ return std::make_unique<FDPF_XB_SparseLUSolver>(); });
    reg.register_solver("FDPF_BX_SparseLU",
        []{ return std::make_unique<FDPF_BX_SparseLUSolver>(); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.register_solver("KLU",
        []{ return std::make_unique<KLUSolver>(); });
    reg.register_solver("KLUSingleSlack",
        []{ return std::make_unique<KLUSolverSingleSlack>(); });
    reg.register_solver("KLUDC",
        []{ return std::make_unique<KLUDCSolver>(); });
    reg.register_solver("FDPF_XB_KLU",
        []{ return std::make_unique<FDPF_XB_KLUSolver>(); });
    reg.register_solver("FDPF_BX_KLU",
        []{ return std::make_unique<FDPF_BX_KLUSolver>(); });
#endif // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    reg.register_solver("NICSLU",
        []{ return std::make_unique<NICSLUSolver>(); });
    reg.register_solver("NICSLUSingleSlack",
        []{ return std::make_unique<NICSLUSolverSingleSlack>(); });
    reg.register_solver("NICSLUDC",
        []{ return std::make_unique<NICSLUDCSolver>(); });
    reg.register_solver("FDPF_XB_NICSLU",
        []{ return std::make_unique<FDPF_XB_NICSLUSolver>(); });
    reg.register_solver("FDPF_BX_NICSLU",
        []{ return std::make_unique<FDPF_BX_NICSLUSolver>(); });
#endif // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    reg.register_solver("CKTSO",
        []{ return std::make_unique<CKTSOSolver>(); });
    reg.register_solver("CKTSOSingleSlack",
        []{ return std::make_unique<CKTSOSolverSingleSlack>(); });
    reg.register_solver("CKTSODC",
        []{ return std::make_unique<CKTSODCSolver>(); });
    reg.register_solver("FDPF_XB_CKTSO",
        []{ return std::make_unique<FDPF_XB_CKTSOSolver>(); });
    reg.register_solver("FDPF_BX_CKTSO",
        []{ return std::make_unique<FDPF_BX_CKTSOSolver>(); });
#endif // CKTSO_SOLVER_AVAILABLE
}

namespace {
    struct _AutoRegister {
        _AutoRegister() { register_builtin_solvers(SolverRegistry::instance()); }
    };
    static const _AutoRegister _auto_reg;
}

} // namespace ls2g
