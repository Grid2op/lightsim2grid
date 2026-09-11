// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LoadContainer.hpp"
#include "BinaryArchive.hpp"

#include <sstream>
#include <iostream>

namespace ls2g {

LoadContainer::StateRes LoadContainer::get_state() const
{
    const auto tmp = get_osc_pq_state();  // osc : one side container
    LoadContainer::StateRes res(tmp);
    return res;
}

void LoadContainer::set_state(LoadContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));  // osc : one side container
    reset_results();
}

void LoadContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

LoadContainer LoadContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<LoadContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
