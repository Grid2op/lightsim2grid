// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LOAD_CONTAINER_H
#define LOAD_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"


class LoadContainer;
class LoadInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        inline LoadInfo(const LoadContainer & r_data_load, int my_id) noexcept;
};


/**
This class is a container for all loads on the grid.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/load.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/load.html#electric-model

NOTE: this class is also used for the storage units! So storage units are modeled as load
which entails that negative storage: the unit is discharging, power is injected in the grid,
positive storage: the unit is charging, power is taken from the grid.
**/
class LoadContainer : public OneSideContainer_PQ, public IteratorAdder<LoadContainer, LoadInfo>
{
    friend class LoadInfo;

    public:
        typedef LoadInfo DataInfo;

    // regular implementation
    public:
        typedef std::tuple<
           OneSideContainer_PQ::StateRes  // state of the base class 
           >  StateRes;
        
        LoadContainer() noexcept = default;
        virtual ~LoadContainer() noexcept = default;
        
        // pickle (python)
        LoadContainer::StateRes get_state() const;
        void set_state(LoadContainer::StateRes & my_state);
        
        void init(const RealVect & load_p_mw,
                  const RealVect & load_q_mvar,
                  const Eigen::VectorXi & load_bus_id
                  )
        {
            init_osc_pq(load_p_mw,
                        load_q_mvar,
                        load_bus_id,
                        "loads");
            reset_results();
        }
    
        virtual void fillSbus(CplxVect & Sbus, const std::vector<SolverBusId> & id_grid_to_solver, bool ac) const;

    protected:
        virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                    const Eigen::Ref<const RealVect> & Vm,
                                    const Eigen::Ref<const CplxVect> & V,
                                    const std::vector<SolverBusId> & id_grid_to_solver,
                                    const RealVect & bus_vn_kv,
                                    real_type sn_mva,
                                    bool ac)
                                    {

                                            set_osc_pq_res_p();
                                            set_osc_pq_res_q(ac);
                                    }
};

inline LoadInfo::LoadInfo(const LoadContainer & r_data_load, int my_id) noexcept: 
        OneSidePQInfo(r_data_load, my_id) {}

#endif  //LOAD_CONTAINER_H
