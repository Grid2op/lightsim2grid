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

void AlgorithmRegistry::register_solver(const std::string& name, Factory f) {
    _factories[name] = std::move(f);
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

std::vector<std::string> AlgorithmRegistry::available_solvers() const {
    std::vector<std::string> names;
    names.reserve(_factories.size());
    for (auto it = _factories.begin(); it != _factories.end(); ++it) {
        names.push_back(it->first);
    }
    return names;
}

} // namespace ls2g
