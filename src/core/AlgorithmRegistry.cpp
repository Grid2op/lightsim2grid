// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "AlgorithmRegistry.hpp"

namespace ls2g {

AlgorithmRegistry& AlgorithmRegistry::instance() {
    static AlgorithmRegistry reg;
    return reg;
}

void AlgorithmRegistry::register_solver(const std::string& name, Factory f, Ls2gAbiTag caller_tag,
                                        SolverOrigin origin) {
    const Ls2gAbiTag& this_core_tag = core_abi_tag();
    if (caller_tag != this_core_tag) {
        throw std::runtime_error(
            "AlgorithmRegistry: refusing to register solver '" + name + "': it was compiled with "
            "different Eigen SIMD/alignment settings than lightsim2grid_core. Loading it would silently "
            "corrupt the heap the first time an Eigen object crosses the BaseAlgo interface. "
            "lightsim2grid_core was built with [" + this_core_tag.describe() + "], this plugin was built "
            "with [" + caller_tag.describe() + "]. Recompile the plugin with matching compiler/-march "
            "flags -- see docs/solver_plugin.rst, section \"Matching build flags\".");
    }
    if (!is_valid_solver_name(name)) {
        throw std::invalid_argument(
            "AlgorithmRegistry: refusing to register the solver name '" + printable(name) + "': " +
            solver_name_rule() + ". This name is written into every grid saved while the solver is "
            "selected, and is what picks it again on load, so it is restricted to a character set "
            "that cannot be confused with another solver's nor corrupt an error message.");
    }
    if (_factories.count(name)) {
        throw std::invalid_argument("AlgorithmRegistry: a solver named '" + printable(name) + "' is already registered.");
    }
    _factories[name] = std::move(f);
    _origins[name] = origin;
}

void AlgorithmRegistry::register_all(FactoryMap batch, Ls2gAbiTag caller_tag, SolverOrigin origin) {
    const Ls2gAbiTag& this_core_tag = core_abi_tag();
    if (caller_tag != this_core_tag) {
        throw std::runtime_error(
            "AlgorithmRegistry: refusing to register a solver plugin: it was compiled with "
            "different Eigen SIMD/alignment settings than lightsim2grid_core. Loading it would silently "
            "corrupt the heap the first time an Eigen object crosses the BaseAlgo interface. "
            "lightsim2grid_core was built with [" + this_core_tag.describe() + "], this plugin was built "
            "with [" + caller_tag.describe() + "]. Recompile the plugin with matching compiler/-march "
            "flags -- see docs/solver_plugin.rst, section \"Matching build flags\".");
    }
    // Reject the whole batch if any name is invalid or already taken, *before*
    // touching the live registry, so a rejected plugin never half-registers.
    for (const auto& kv : batch) {
        if (!is_valid_solver_name(kv.first)) {
            throw std::invalid_argument(
                "AlgorithmRegistry: refusing to register the solver name '" + printable(kv.first) + "': " +
                solver_name_rule() + ". This name is written into every grid saved while the solver is "
                "selected, and is what picks it again on load, so it is restricted to a character set "
                "that cannot be confused with another solver's nor corrupt an error message.");
        }
        if (_factories.count(kv.first)) {
            throw std::invalid_argument(
                "AlgorithmRegistry: a solver named '" + printable(kv.first) + "' is already registered.");
        }
    }
    // Strong exception guarantee: assemble the merged table off to the side,
    // then commit it with a noexcept swap. If anything above or in the merge
    // throws (e.g. std::bad_alloc), _factories is left exactly as it was.
    FactoryMap merged = _factories;
    merged.reserve(_factories.size() + batch.size());
    std::unordered_map<std::string, SolverOrigin> merged_origins = _origins;
    merged_origins.reserve(_origins.size() + batch.size());
    for (auto& kv : batch) {
        merged_origins.emplace(kv.first, origin);
        merged.emplace(kv.first, std::move(kv.second));
    }
    _factories.swap(merged);
    _origins.swap(merged_origins);
}

std::unique_ptr<BaseAlgo> AlgorithmRegistry::make(const std::string& name) const {
    auto it = _factories.find(name);
    if (it == _factories.end()) {
        throw std::runtime_error("AlgorithmRegistry: unknown solver '" + name + "'. "
                                 "Check available_solver_names() for the list of registered solvers.");
    }
    return it->second();
}

bool AlgorithmRegistry::is_registered(const std::string& name) const {
    return _factories.find(name) != _factories.end();
}

SolverOrigin AlgorithmRegistry::origin_of(const std::string& name) const {
    auto it = _origins.find(name);
    // an unknown name is reported External: the cautious answer (see the header)
    if (it == _origins.end()) return SolverOrigin::External;
    return it->second;
}

std::vector<std::string> AlgorithmRegistry::builtin_algorithm_names() const {
    std::vector<std::string> names;
    for (auto it = _origins.begin(); it != _origins.end(); ++it) {
        if (it->second == SolverOrigin::Builtin) names.push_back(it->first);
    }
    return names;
}

std::vector<std::string> AlgorithmRegistry::available_algorithm_names() const {
    std::vector<std::string> names;
    names.reserve(_factories.size());
    for (auto it = _factories.begin(); it != _factories.end(); ++it) {
        names.push_back(it->first);
    }
    return names;
}

} // namespace ls2g
