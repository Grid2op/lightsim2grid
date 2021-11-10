// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include <unordered_map>
#include "IIDMConverter.h"

#include "powsybl/iidm/Substation.hpp"
#include "powsybl/iidm/Bus.hpp"
#include "powsybl/iidm/VoltageLevel.hpp"
#include "powsybl/iidm/Line.hpp"
#include "powsybl/iidm/Terminal.hpp"
#include "powsybl/iidm/Branch.hpp"
#include "powsybl/iidm/TwoWindingsTransformer.hpp"

GridModel GridModelFromIIDM::get_grid_model() const {
    GridModel res;

    const auto & bus_view = init_bus(res);
    res.set_init_vm_pu(1.04);  // by default i will initialize the voltage magnitudes at 1.04 pu
    res.set_sn_mva(1.0);  // TODO for later if i get errors

    init_powerlines(res, bus_view);
    return res;
}

powsybl::iidm::network::BusView GridModelFromIIDM::init_bus(GridModel & grid_model) const{
    // init the buses

    // TODO I need to change that... Trafo are within a substation
    const auto & bus_view = network_.getBusView();
    const auto & all_buses = bus_view.getBuses();
    unsigned int nb_bus = std::distance(all_buses.begin(), all_buses.end());
    RealVect bus_vn_kv = RealVect(nb_bus);
    int bus_id = 0;
    for(const powsybl::iidm::Bus & bus : all_buses){
        const auto & level =  bus.getVoltageLevel();
        bus_vn_kv(bus_id) = level.getNominalV();
        ++bus_id;
    }
    const int nb_line =  network_.getLineCount();
    const int nb_trafo = network_.getTwoWindingsTransformerCount() ;
    grid_model.init_bus(bus_vn_kv, nb_line, nb_trafo);
    return bus_view;
}

void GridModelFromIIDM::init_powerlines(GridModel & grid_model,
                                        const powsybl::iidm::network::BusView & bus_view) const{
    const int nb_line = network_.getLineCount();
    RealVect branch_r = RealVect(nb_line);
    RealVect branch_x = RealVect(nb_line);
    CplxVect branch_h = CplxVect(nb_line);
    Eigen::VectorXi branch_from_id = Eigen::VectorXi(nb_line);
    Eigen::VectorXi branch_to_id = Eigen::VectorXi(nb_line);

    std::unordered_map<std::string, int> lines_already_there;
    int line_id = 0;
    int bus_id = 0;
    for(const powsybl::iidm::Bus & bus : bus_view.getBuses()){
        const auto & all_lines = bus.getLines();
        for(const powsybl::iidm::Line & powerline : all_lines){
            std::string l_id_str = powerline.getId();

            const auto & it = lines_already_there.find(l_id_str);
            if(it == lines_already_there.end())
            {
                // line has not been added => this means i'm connected to "origin" and that i need to insert it in
                // the "lines" data

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
                branch_from_id(line_id) = bus_id; // TODO check that !
                lines_already_there[l_id_str] = line_id;
                ++line_id;
            }else{
                // lines has already been added in a previous nodes, this means it's the "end" side
                int line_id_tmp = it->second;
                branch_to_id(line_id_tmp) = bus_id; // TODO check that !
            }
        }
        ++bus_id;
    }
    grid_model.init_powerlines(branch_r,
                               branch_x,
                               branch_h,
                               branch_from_id,
                               branch_to_id
                               );
}

void GridModelFromIIDM::init_trafos(GridModel & grid_model,
                                    const powsybl::iidm::network::BusView & bus_view) const{
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

    std::unordered_map<std::string, int> trafos_already_there;
    int trafo_id = 0;
    int bus_id = 0;
    for(const powsybl::iidm::Bus & bus : bus_view.getBuses()){
        const auto & all_trafos = bus.getTwoWindingsTransformers();
        for(const powsybl::iidm::TwoWindingsTransformer & trafo : all_trafos){
            std::string t_id_str = trafo.getId();

            const auto & it = trafos_already_there.find(t_id_str);
            if(it == trafos_already_there.end())
            {
                // line has not been added => this means i'm connected to "origin" and that i need to insert it in
                // the "lines" data

                // TODO pair unit !
                trafo_r(trafo_id) = trafo.getR();
                trafo_x(trafo_id) = trafo.getX();
                trafo_b(trafo_id) = {trafo.getG(), trafo.getB()};

                // trafo_tap_step_pct() // TODO
                // trafo_tap_pos() // TODO
                // trafo_shift_degree() // TODO
                // trafo_tap_hv() // TODO

                // trafo_hv_id(trafo_id) = bus_id; // TODO trafo_hv_id or trafo_lv_id
                trafos_already_there[t_id_str] = trafo_id;
                ++trafo_id;
            }else{
                // lines has already been added in a previous nodes, this means it's the "end" side
                int trafo_id_tmp = it->second;
                // trafo_lv_id(trafo_id_tmp) = trafo_id; // TODO trafo_hv_id or trafo_lv_id
            }
        }
        ++bus_id;
    }
    grid_model.init_trafo(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_shift_degree,
                          trafo_tap_hv, trafo_hv_id, trafo_lv_id);
}