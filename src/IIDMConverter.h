// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef IIDM_CONVERTER_AVAILABLE

#ifndef IIDMCONVERTER_H
#define IIDMCONVERTER_H

#include "GridModel.h"

#include <iostream>
#include <libxml/parser.h>
#include "powsybl/PowsyblException.hpp"
#include "powsybl/iidm/Network.hpp"

/**
Supposes that, at the init of the file, every element on the same substation is connected on the same busbar !
**/
class GridModelFromIIDM
{
    GridModelFromIIDM(const std::string & path):network_("not init", "dont use"){
        xmlInitParser();
        powsybl::iidm::Network network_ = powsybl::iidm::Network::readXml("/home/benjamin/Documents/powsybl-iidm4cpp/examples/example1/eurostag-tutorial1.xml");
        xmlCleanupParser();
        // https://javadoc.io/doc/com.powsybl/powsybl-core/latest/com/powsybl/iidm/network/Network.html
        std::cout << "getGeneratorCount: " << network_.getGeneratorCount() << std::endl;
    }

    GridModel get_grid_model() const;

    protected:
        powsybl::iidm::Network network_;

        powsybl::iidm::network::BusView init_bus(GridModel & grid_model) const;
        void init_powerlines(GridModel & grid_model, const powsybl::iidm::network::BusView & bus_view) const;
        void init_trafos(GridModel & grid_model, const powsybl::iidm::network::BusView & bus_view) const;
};

#endif  // IIDMCONVERTER_H
#endif // IIDM_CONVERTER_AVAILABLE
