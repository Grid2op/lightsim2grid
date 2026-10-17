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
#include <string>
#include <vector>

#include "BaseOuterLoop.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's default outer-loop list (DefaultAcOuterLoopConfig), in its order, with
 * OpenLoadFlow's own default parameters: the loops these parameters turn on, among those
 * lightsim2grid implements. A loop that has nothing to do on `grid` is still listed; its
 * is_needed() drops it at solve time. A grid loaded from pypowsybl holds instead the list
 * the OpenLoadFlow pypowsybl ships would run (see init_from_pypowsybl).
 */
LS2G_API std::vector<std::shared_ptr<BaseOuterLoop> > make_default_outer_loops(const LSGrid & grid);

/**
 * The loop called `name` (BaseOuterLoop::name) with the parameters `params`, flattened as
 * its get_params() gives them: what a saved loop list is restored from. Throws on an
 * unknown name, a parameter count the loop does not have, or values its constructor refuses.
 */
LS2G_API std::shared_ptr<BaseOuterLoop> make_outer_loop(const std::string & name, const AlgoConfig & params);

/**
 * Detection mode: what each loop of `loops` would do after the solve `ctx` describes, as
 * LimitViolations appended to `out`, in the order of the list. Exactly the loop's own first
 * check, without its action: each loop is initialized -- on a clone, so the list is left as
 * it is -- and detects. `ctx` holds no injection and no control (OuterContext::is_detection).
 */
LS2G_API void outer_detect(const std::vector<std::shared_ptr<BaseOuterLoop> > & loops,
                           OuterContext & ctx, std::vector<LimitViolation> & out);

}  // namespace ls2g

#endif  // DEFAULT_OUTER_LOOPS_H
