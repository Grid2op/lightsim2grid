// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DEFAULT_OUTER_LOOPS_H
#define DEFAULT_OUTER_LOOPS_H

#include <memory>
#include <vector>

#include "BaseOuterLoop.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's default outer-loop list (DefaultAcOuterLoopConfig), in its order, with
 * the default parameters of docs/dev_notes/outer_loops_fixed_sparsity.md: the loops
 * these parameters turn on, among those lightsim2grid implements. A loop that has nothing
 * to do on `grid` is still listed; its is_needed() drops it at solve time.
 */
LS2G_API std::vector<std::shared_ptr<BaseOuterLoop> > make_default_outer_loops(const LSGrid & grid);

}  // namespace ls2g

#endif  // DEFAULT_OUTER_LOOPS_H
