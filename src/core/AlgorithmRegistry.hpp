// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ALGORITHM_REGISTRY_H
#define ALGORITHM_REGISTRY_H

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "powerflow_algorithm/BaseAlgo.hpp"

namespace ls2g {

class LS2G_API AlgorithmRegistry {
public:
    using Factory = std::function<std::unique_ptr<BaseAlgo>()>;

    static AlgorithmRegistry& instance();

    void register_solver(const std::string& name, Factory f);
    std::unique_ptr<BaseAlgo> make(const std::string& name) const;
    bool is_registered(const std::string& name) const;
    std::vector<std::string> available_solvers() const;

private:
    AlgorithmRegistry() = default;
    std::unordered_map<std::string, Factory> _factories;
};

// Helper for plugins: declare a static instance of this class in an anonymous
// namespace inside one translation unit to register a solver at load time.
// No macro needed — the static constructor fires when the .so is dlopen'd.
class LS2G_API AlgorithmRegistrar {
public:
    AlgorithmRegistrar(const std::string& name, AlgorithmRegistry::Factory f) {
        AlgorithmRegistry::instance().register_solver(name, std::move(f));
    }
};


} // namespace ls2g

#endif // ALGORITHM_REGISTRY_H
