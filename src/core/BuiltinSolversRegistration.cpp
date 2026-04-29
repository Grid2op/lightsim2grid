// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BuiltinSolversRegistration.hpp"
#include "AlgorithmRegistry.hpp"
#include "Solvers.hpp"   // all concrete solver typedefs + #ifdef guards

namespace ls2g {

void register_builtin_solvers(AlgorithmRegistry& reg) {
    reg.register_solver("NR_SparseLU",
        []{ return std::make_unique<NR_SparseLU>(); });
    reg.register_solver("NRSing_SparseLU",
        []{ return std::make_unique<NRSing_SparseLU>(); });
    reg.register_solver("DC_SparseLU",
        []{ return std::make_unique<DC_SparseLU>(); });
    reg.register_solver("GaussSeidel",
        []{ return std::make_unique<GaussSeidelAlgo>(); });
    reg.register_solver("GaussSeidelSynch",
        []{ return std::make_unique<GaussSeidelSynchAlgo>(); });
    reg.register_solver("FDPF_XB_SparseLU",
        []{ return std::make_unique<FDPF_XB_SparseLU>(); });
    reg.register_solver("FDPF_BX_SparseLU",
        []{ return std::make_unique<FDPF_BX_SparseLU>(); });

#ifdef KLU_SOLVER_AVAILABLE
    reg.register_solver("NR_KLU",
        []{ return std::make_unique<NR_KLU>(); });
    reg.register_solver("NRSing_KLU",
        []{ return std::make_unique<NRSing_KLU>(); });
    reg.register_solver("DC_KLU",
        []{ return std::make_unique<DC_KLU>(); });
    reg.register_solver("FDPF_XB_KLU",
        []{ return std::make_unique<FDPF_XB_KLU>(); });
    reg.register_solver("FDPF_BX_KLU",
        []{ return std::make_unique<FDPF_BX_KLU>(); });
#endif // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    reg.register_solver("NR_NICSLU",
        []{ return std::make_unique<NR_NICSLU>(); });
    reg.register_solver("NRSing_NICSLU",
        []{ return std::make_unique<NRSing_NICSLU>(); });
    reg.register_solver("DC_NICSLU",
        []{ return std::make_unique<DC_NICSLU>(); });
    reg.register_solver("FDPF_XB_NICSLU",
        []{ return std::make_unique<FDPF_XB_NICSLU>(); });
    reg.register_solver("FDPF_BX_NICSLU",
        []{ return std::make_unique<FDPF_BX_NICSLU>(); });
#endif // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    reg.register_solver("NR_CKTSO",
        []{ return std::make_unique<NR_CKTSO>(); });
    reg.register_solver("NRSing_CKTSO",
        []{ return std::make_unique<NRSing_CKTSO>(); });
    reg.register_solver("DC_CKTSO",
        []{ return std::make_unique<DC_CKTSO>(); });
    reg.register_solver("FDPF_XB_CKTSO",
        []{ return std::make_unique<FDPF_XB_CKTSO>(); });
    reg.register_solver("FDPF_BX_CKTSO",
        []{ return std::make_unique<FDPF_BX_CKTSO>(); });
#endif // CKTSO_SOLVER_AVAILABLE
}

namespace {
    struct _AutoRegister {
        _AutoRegister() { register_builtin_solvers(AlgorithmRegistry::instance()); }
    };
    static const _AutoRegister _auto_reg;
}

} // namespace ls2g
