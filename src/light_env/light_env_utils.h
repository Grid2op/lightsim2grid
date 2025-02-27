// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LIGHT_ENV_UTIL_H
#define LIGHT_ENV_UTIL_H

enum class ElementType {load, gen, line_or, line_ex, storage, trafo_hv, trafo_lv, 
                        shunt, static_gen, dc_line_or, dc_line_ex};

std::ostream& operator<<(std::ostream& out, const ElementType& solver_type)
{
    switch (solver_type)
    {
    case ElementType::load:
        out << "load";
        break;
    case ElementType::gen:
        out << "gen";
        break;
    case ElementType::line_or:
        out << "line (or side)";
        break;
    case ElementType::line_ex:
        out << "line (ex side)";
        break;
    case ElementType::storage:
        out << "storage";
        break;
    case ElementType::trafo_hv:
        out << "trafo (hv side)";
        break;
    case ElementType::trafo_lv:
        out << "trafo (lv side)";
        break;
    case ElementType::shunt:
        out << "shunt";
        break;
    case ElementType::static_gen:
        out << "static gen";
        break;
    case ElementType::dc_line_or:
        out << "dc line (or side)";
        break;
    case ElementType::dc_line_ex:
        out << "dc line (ex side)";
        break;
    default:
        out << "(unknown)";
        break;
    }
   return out;
}

#endif  // LIGHT_ENV_UTIL_H