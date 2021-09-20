// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "IIDMConverter.h"


#include "powsybl/iidm/Substation.hpp"
#include "powsybl/iidm/VoltageLevel.hpp"
#include "powsybl/iidm/Line.hpp"
#include "powsybl/iidm/Terminal.hpp"
#include "powsybl/iidm/Branch.hpp"
#include "powsybl/iidm/TwoWindingsTransformer.hpp"

GridModel GridModelFromIIDM::get_grid_model() const {
    GridModel res;

    init_bus(res);
    res.set_init_vm_pu(1.04);  // by default i will initialize the voltage magnitudes at 1.04 pu
    res.set_sn_mva(1.0);  // TODO for later if i get errors

    init_powerlines(res);
    return res;
}

void GridModelFromIIDM::init_bus(GridModel & grid_model) const{
    // init the buses

    // TODO I need to change that... Trafo are within a substation
    RealVect bus_vn_kv = RealVect(network_.getSubstationCount());
    int sub_id = 0;
    for(const powsybl::iidm::Substation & substation : network_.getSubstations()){
        const auto & levels =  substation.getVoltageLevels();
        auto volt_level_it = levels.begin();
        bus_vn_kv(sub_id) = volt_level_it->getNominalV();
        ++volt_level_it;
        if(volt_level_it != levels.end()){
            std::cout << "WARNING: more that one level found for a substation. This is not taken into account in GridModelFromIIDM at the moment." << std::endl;
        }
        ++sub_id;
    }
    const int nb_line =  network_.getLineCount();
    const int nb_trafo = network_.getTwoWindingsTransformerCount() ;
    grid_model.init_bus(bus_vn_kv, nb_line, nb_trafo);
}

void GridModelFromIIDM::init_powerlines(GridModel & grid_model) const{
    const int nb_line = network_.getLineCount();
    RealVect branch_r = RealVect(nb_line);
    RealVect branch_x = RealVect(nb_line);
    CplxVect branch_h = CplxVect(nb_line);
    Eigen::VectorXi branch_from_id = Eigen::VectorXi(nb_line);
    Eigen::VectorXi branch_to_id = Eigen::VectorXi(nb_line);

    int line_id = 0;
    for(const powsybl::iidm::Line & powerline : network_.getLines()){
        // TODO pair unit !
        branch_r(line_id) = powerline.getR();
        branch_x(line_id) = powerline.getX();
        cplx_type h = {powerline.getG1() + powerline.getG2(), powerline.getB1() + powerline.getB2()};
        if(powerline.getG1() != powerline.getG2()){
            std::cout << "WARNING: some powerline have a different `g` at the origin and extremity side. This is not supported by GridModelFromIIDM and we assume the are equal (to the average)" << std::endl;
        }
        if(powerline.getB1() != powerline.getB2()){
            std::cout << "WARNING: some powerline have a different `b` at the origin and extremity side. This is not supported by GridModelFromIIDM and we assume the are equal (to the average)" << std::endl;
        }
        branch_h(line_id) = h;
//        const auto & sub = powerline.getSubstation1();
        powsybl::iidm::Terminal const * p_terminal1 = &powerline.getTerminal1(); //.NodeBreakerView().getNode();
        powsybl::iidm::Terminal const * p_terminal2 = &powerline.getTerminal1(); //.NodeBreakerView().getNode();
        branch_from_id(line_id) = p_terminal1->getNodeBreakerView().getNode(); // TODO check that !
        branch_to_id(line_id) = p_terminal2->getNodeBreakerView().getNode(); // TODO check that !
        ++line_id;
    }
    grid_model.init_powerlines(branch_r,
                               branch_x,
                               branch_h,
                               branch_from_id,
                               branch_to_id
                               );
}

void GridModelFromIIDM::init_trafos(GridModel & grid_model) const{
    if(network_.getThreeWindingsTransformerCount() > 0){
        std::cout << "WARNING: there are some 3 windings transformer on the grid, they will be discarded by GridModelFromIIDM" << std::endl;
    }
    const int nb_trafo = network_.getTwoWindingsTransformerCount() ;
    RealVect trafo_r = RealVect(nb_trafo);
    RealVect trafo_x = RealVect(nb_trafo);
    CplxVect trafo_b = CplxVect(nb_trafo);
    RealVect trafo_tap_step_pct = RealVect(nb_trafo);
    RealVect trafo_tap_pos = RealVect(nb_trafo);
    RealVect trafo_shift_degree = RealVect(nb_trafo);
    std::vector<bool> trafo_tap_hv = std::vector<bool>(nb_trafo);
    Eigen::VectorXi trafo_hv_id = Eigen::VectorXi(nb_trafo);
    Eigen::VectorXi trafo_lv_id = Eigen::VectorXi(nb_trafo);
    int trafo_id = 0;
    for(const powsybl::iidm::TwoWindingsTransformer & trafo : network_.getTwoWindingsTransformers()){
        // TODO pair unit !
        trafo_r(trafo_id) = trafo.getR();
        trafo_x(trafo_id) = trafo.getX();
        trafo_b(trafo_id) = {trafo.getG(), trafo.getB()};

        // trafo_tap_step_pct() // TODO
        // trafo_tap_pos() // TODO
        // trafo_shift_degree() // TODO
        // trafo_tap_hv() // TODO
        // trafo_hv_id() // TODO
        // trafo_lv_id() // TODO

        ++trafo_id;
    }
    grid_model.init_trafo(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_shift_degree,
                          trafo_tap_hv, trafo_hv_id, trafo_lv_id);
}