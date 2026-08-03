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
    // Every solver registered here is shipped in-tree, so it is tagged
    // SolverOrigin::Builtin. That tag -- NOT AlgorithmType -- is what tells
    // lightsim2grid apart from a plugin at run time: AlgorithmType is a fixed enum
    // of serialized solver identities and the NRRefactorRetry_* family has no
    // member there, so name_to_algo_type() reports them as AlgorithmType::Custom
    // exactly like a plugin. Going through this helper means a built-in added later
    // cannot silently be treated as external by forgetting the tag.
    const auto add_builtin = [&reg](const std::string& name, AlgorithmRegistry::Factory f) {
        reg.register_solver(name, std::move(f), ls2g_current_abi_tag(), SolverOrigin::Builtin);
    };
    add_builtin("NR_SparseLU",
        []{ return std::make_unique<NR_SparseLU>(); });
    add_builtin("NRSing_SparseLU",
        []{ return std::make_unique<NRSing_SparseLU>(); });
    add_builtin("DC_SparseLU",
        []{ return std::make_unique<DC_SparseLU>(); });
    add_builtin("GaussSeidel",
        []{ return std::make_unique<GaussSeidelAlgo>(); });
    add_builtin("GaussSeidelSynch",
        []{ return std::make_unique<GaussSeidelSynchAlgo>(); });
    add_builtin("FDPF_XB_SparseLU",
        []{ return std::make_unique<FDPF_XB_SparseLU>(); });
    add_builtin("FDPF_BX_SparseLU",
        []{ return std::make_unique<FDPF_BX_SparseLU>(); });

#ifdef KLU_SOLVER_AVAILABLE
    add_builtin("NR_KLU",
        []{ return std::make_unique<NR_KLU>(); });
    add_builtin("NRSing_KLU",
        []{ return std::make_unique<NRSing_KLU>(); });
    add_builtin("DC_KLU",
        []{ return std::make_unique<DC_KLU>(); });
    add_builtin("FDPF_XB_KLU",
        []{ return std::make_unique<FDPF_XB_KLU>(); });
    add_builtin("FDPF_BX_KLU",
        []{ return std::make_unique<FDPF_BX_KLU>(); });
    add_builtin("NRRefactorRetry_KLU",
        []{ return std::make_unique<NRRefactorRetry_KLU>(); });
#endif // KLU_SOLVER_AVAILABLE

#ifdef NICSLU_SOLVER_AVAILABLE
    add_builtin("NR_NICSLU",
        []{ return std::make_unique<NR_NICSLU>(); });
    add_builtin("NRSing_NICSLU",
        []{ return std::make_unique<NRSing_NICSLU>(); });
    add_builtin("DC_NICSLU",
        []{ return std::make_unique<DC_NICSLU>(); });
    add_builtin("FDPF_XB_NICSLU",
        []{ return std::make_unique<FDPF_XB_NICSLU>(); });
    add_builtin("FDPF_BX_NICSLU",
        []{ return std::make_unique<FDPF_BX_NICSLU>(); });
    add_builtin("NRRefactorRetry_NICSLU",
        []{ return std::make_unique<NRRefactorRetry_NICSLU>(); });
#endif // NICSLU_SOLVER_AVAILABLE

#ifdef CKTSO_SOLVER_AVAILABLE
    add_builtin("NR_CKTSO",
        []{ return std::make_unique<NR_CKTSO>(); });
    add_builtin("NRSing_CKTSO",
        []{ return std::make_unique<NRSing_CKTSO>(); });
    add_builtin("DC_CKTSO",
        []{ return std::make_unique<DC_CKTSO>(); });
    add_builtin("FDPF_XB_CKTSO",
        []{ return std::make_unique<FDPF_XB_CKTSO>(); });
    add_builtin("FDPF_BX_CKTSO",
        []{ return std::make_unique<FDPF_BX_CKTSO>(); });
    add_builtin("NRRefactorRetry_CKTSO",
        []{ return std::make_unique<NRRefactorRetry_CKTSO>(); });
#endif // CKTSO_SOLVER_AVAILABLE
}

namespace {
    struct _AutoRegister {
        _AutoRegister() { register_builtin_solvers(AlgorithmRegistry::instance()); }
    };
    static const _AutoRegister _auto_reg;
}

} // namespace ls2g
