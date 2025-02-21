// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TOPO_ACTION_H
#define TOPO_ACTION_H

#include <vector>
#include "GridModel.h"

class TopoAction
{

    public:
        // TODO
        TopoAction() {}

        // TODO
        void apply_to_gridmodel(GridModel & grid) const {

        }

        // TODO check compliance (correct sub, correct bus)
        void check_validity(const GridModel & grid) const {

        }
        
    protected:
        int sub_id;

        // bus id should be global (lightsim2grid / gridmodel) bus and
        // not local bus id
        std::vector<std::pair<int, int> > loads_bus;
        std::vector<std::pair<int, int> > gens_bus;
        std::vector<std::pair<int, int> > static_gens_bus;
        std::vector<std::pair<int, int> > shunts_bus;
        std::vector<std::pair<int, int> > lines_or_bus;
        std::vector<std::pair<int, int> > lines_ex_bus;
        std::vector<std::pair<int, int> > trafo_hv_bus;
        std::vector<std::pair<int, int> > trafo_lv_bus;
        std::vector<std::pair<int, int> > dc_lines_or_bus;
        std::vector<std::pair<int, int> > dc_lines_ex_bus;
        
};

#endif // TOPO_ACTION_H
