// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ShuntContainer.h"

#include <iostream>

ShuntContainer::StateRes ShuntContainer::get_state() const
{
     ShuntContainer::StateRes res(OneSideContainer::get_osc_state());
     return res;
}

void ShuntContainer::set_state(ShuntContainer::StateRes & my_state )
{
    OneSideContainer::set_osc_state(std::get<0>(my_state));
    reset_results();
}

void ShuntContainer::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                              bool ac,
                              const std::vector<int> & id_grid_to_solver,
                              real_type sn_mva) const
{
    if(!ac) return; // no shunt in DC

    const Eigen::Index nb_shunt = static_cast<int>(q_mvar_.size());
    cplx_type tmp;
    int bus_id_me, bus_id_solver;
    for(Eigen::Index shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        // assign diagonal coefficient
        tmp = {p_mw_(shunt_id), -q_mvar_(shunt_id)};

        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "ShuntContainer::fillYbus: the shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        if(sn_mva != 1.) tmp /= sn_mva;
        res.push_back(Eigen::Triplet<cplx_type> (bus_id_solver, bus_id_solver, tmp));
    }
}

void ShuntContainer::fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<int> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const
{
    const Eigen::Index nb_shunt = static_cast<int>(q_mvar_.size());
    real_type tmp;
    int bus_id_me, bus_id_solver;
    for(Eigen::Index shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "ShuntContainer::fillBp_Bpp: the shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        // assign diagonal coefficient
        tmp = q_mvar_(shunt_id);
        if(sn_mva != 1.) tmp /= sn_mva;
        Bpp.push_back(Eigen::Triplet<real_type> (bus_id_solver, bus_id_solver, tmp));  // -(-tmp) [-tmp for the "correct" value, but then for Bpp i have -(-tmp)]
    }
}

void ShuntContainer::fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const  // in DC i need that
{
    if(ac) return;  // in AC I do not do that
    // std::cout << " ok i use this function" << std::endl;
    // - bus[:, GS] / baseMVA  # in pandapower
    // yish=gish+jbish -> so g is the MW !
    const int nb_shunt = static_cast<int>(q_mvar_.size());
    int bus_id_me, bus_id_solver;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;
        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            throw std::runtime_error("GridModel::fillSbus: A shunt is connected to a disconnected bus.");
        }
        Sbus.coeffRef(bus_id_solver) -= p_mw_(shunt_id);
    }
}

void ShuntContainer::_compute_results(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<int> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac)
{
    const int nb_shunt = static_cast<int>(p_mw_.size());
    for(int shunt_id = 0; shunt_id < nb_shunt; ++shunt_id){
        if(!status_[shunt_id]) {
            res_p_(shunt_id) = my_zero_;
            res_q_(shunt_id) = my_zero_;
            continue;
        }
        int bus_id_me = bus_id_(shunt_id);
        int bus_solver_id = id_grid_to_solver[bus_id_me];
        if(bus_solver_id == _deactivated_bus_id){
            throw std::runtime_error("ShuntContainer::compute_results: A shunt is connected to a disconnected bus.");
        }
        cplx_type E = V(bus_solver_id);
        cplx_type y = -my_one_ * (p_mw_(shunt_id) + my_i * q_mvar_(shunt_id)) / sn_mva;
        cplx_type I = y * E;
        I = std::conj(I);
        cplx_type s = E * I;
        res_p_(shunt_id) = std::real(s) * sn_mva;
        if(ac) res_q_(shunt_id) = std::imag(s) * sn_mva;
        else res_q_(shunt_id) = my_zero_;
    }
}
