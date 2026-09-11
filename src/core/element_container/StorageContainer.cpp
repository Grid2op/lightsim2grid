// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "StorageContainer.hpp"
#include "BinaryArchive.hpp"

#include <sstream>
#include <iostream>

namespace ls2g {

StorageContainer::StateRes StorageContainer::get_state() const
{
    const auto tmp = get_osc_pq_state();  // osc : one side container
    StorageContainer::StateRes res(tmp);
    return res;
}

void StorageContainer::set_state(StorageContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<StateResIdx::OSC_PQ_STATE>(my_state));  // osc : one side container
    reset_results();
}

void StorageContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

StorageContainer StorageContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<StorageContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
