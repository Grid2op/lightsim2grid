// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SHUNT_CONTAINER_H
#define SHUNT_CONTAINER_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"

class ShuntContainer;
class ShuntInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        // no members
        inline ShuntInfo(const ShuntContainer & r_data_shunt, int my_id);
};

/**
This class is a container for all shunts on the grid.

The convention used for the shunt is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html#electric-model
**/
class ShuntContainer : public OneSideContainer_PQ, public IteratorAdder<ShuntContainer, ShuntInfo>
{
    friend class ShuntInfo;
    public:
        typedef ShuntInfo DataInfo;

    public:
        typedef std::tuple<OneSideContainer_PQ::StateRes >  StateRes;
        
        ShuntContainer():OneSideContainer_PQ() {};
        
        
        void init(const RealVect & shunt_p_mw,
                  const RealVect & shunt_q_mvar,
                  const Eigen::VectorXi & shunt_bus_id
                  )
        {
            init_osc_pq(shunt_p_mw,
                        shunt_q_mvar,
                        shunt_bus_id,
                        "shunts");
            reset_results();
        }
    
        // pickle (python)
        ShuntContainer::StateRes get_state() const;
        void set_state(ShuntContainer::StateRes & my_state );
        
        virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                              bool ac,
                              const std::vector<SolverBusId> & id_grid_to_solver,
                              real_type sn_mva) const;
        virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<SolverBusId> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const;
        virtual void fillSbus(CplxVect & Sbus, const std::vector<SolverBusId> & id_grid_to_solver, bool ac) const;  // in DC i need that
        
    protected:
        virtual void _change_p(int shunt_id, real_type new_p, bool my_status, SolverControl & solver_control)
        {
            if(abs(target_p_mw_(shunt_id) - new_p) > _tol_equal_float){
                solver_control.tell_recompute_ybus();
                solver_control.tell_recompute_sbus();  // needed for DC
            }
        }

        virtual void _change_q(int shunt_id, real_type new_q, bool my_status, SolverControl & solver_control)
        {
            if(abs(target_q_mvar_(shunt_id) - new_q) > _tol_equal_float){
                solver_control.tell_recompute_ybus();
            }
        }

        virtual void _change_bus(int el_id, GridModelBusId new_bus_id, SolverControl & solver_control, int nb_bus) {
            if(bus_id_(el_id) != new_bus_id){
                solver_control.tell_recompute_ybus();
            }
        };
        virtual void _deactivate(int el_id, SolverControl & solver_control) {
            if(status_[el_id]){
                solver_control.tell_recompute_ybus();
            }
        };
        virtual void _reactivate(int el_id, SolverControl & solver_control) {
            if(!status_[el_id]){
                solver_control.tell_recompute_ybus();
            }
        };
        virtual void _compute_results(
            const Eigen::Ref<const RealVect> & Va,
            const Eigen::Ref<const RealVect> & Vm,
            const Eigen::Ref<const CplxVect> & V,
            const std::vector<SolverBusId> & id_grid_to_solver,
            const RealVect & bus_vn_kv,
            real_type sn_mva,
            bool ac);

    protected:
        // physical properties

        // input data

        //output data

};

inline ShuntInfo::ShuntInfo(const ShuntContainer & r_data_shunt, int my_id):OneSidePQInfo(r_data_shunt, my_id){}

#endif  //SHUNT_CONTAINER_H
