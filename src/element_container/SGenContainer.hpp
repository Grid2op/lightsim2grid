// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SGEN_CONTAINER_H
#define SGEN_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"


class SGenContainer;
class SGenInfo  : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        // members
        real_type min_q_mvar;
        real_type max_q_mvar;
        real_type min_p_mw;
        real_type max_p_mw;

        inline SGenInfo(const SGenContainer & r_data_sgen, int my_id) noexcept;
};

/**
This class is a container for all static generator (PQ generators) on the grid.
They are given in the generator convention: positive sign for P,Q => the power is produced.

The convention used for the static is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html#electric-model
**/
class SGenContainer: public OneSideContainer_PQ, public IteratorAdder<SGenContainer, SGenInfo>
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    friend class SGenInfo;

    public:
        typedef SGenInfo DataInfo;

    public:
        typedef std::tuple<
           OneSideContainer_PQ::StateRes,
           std::vector<real_type>, // p_min
           std::vector<real_type>, //  p_max
           std::vector<real_type>, //  q_min
           std::vector<real_type> //  q_max
           >  StateRes;
        
        SGenContainer() noexcept = default;
        virtual ~SGenContainer() noexcept = default;
        
        // pickle (python)
        SGenContainer::StateRes get_state() const;
        void set_state(SGenContainer::StateRes & my_state );
        
        
        void init(const RealVect & sgen_p,
                  const RealVect & sgen_q,
                  const RealVect & sgen_pmin,
                  const RealVect & sgen_pmax,
                  const RealVect & sgen_qmin,
                  const RealVect & sgen_qmax,
                  const Eigen::VectorXi & sgen_bus_id
                  );
              
        virtual void fillSbus(CplxVect & Sbus, const std::vector<SolverBusId> & id_grid_to_solver, bool ac) const ;

    protected:
        virtual void _compute_results(
            const Eigen::Ref<const RealVect> & Va,
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

    private:
        // physical properties
        RealVect p_min_mw_;
        RealVect p_max_mw_;
        RealVect q_min_mvar_;
        RealVect q_max_mvar_;

        // input data

        //output data
};

inline  SGenInfo::SGenInfo(const SGenContainer & r_data_sgen, int my_id) noexcept:
OneSidePQInfo(r_data_sgen, my_id),
min_q_mvar(0.),
max_q_mvar(0.),
min_p_mw(0.),
max_p_mw(0.)
{
    if((my_id >= 0) && (my_id < r_data_sgen.nb()))
    {
        min_q_mvar = r_data_sgen.q_min_mvar_(my_id);
        max_q_mvar = r_data_sgen.q_max_mvar_(my_id);
        min_p_mw = r_data_sgen.p_min_mw_(my_id);
        max_p_mw = r_data_sgen.p_max_mw_(my_id);
    }
}

#endif  //SGEN_CONTAINER_H
