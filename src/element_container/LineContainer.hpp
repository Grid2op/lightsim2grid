// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LINE_CONTAINER_H
#define LINE_CONTAINER_H

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SubstationContainer.hpp"
#include "OneSideContainer_forBranch.hpp"
#include "TwoSidesContainer_rxh_A.hpp"

class LineContainer;
class LineInfo : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::TwoSidesContainer_rxh_AInfo
{
    public:
        inline LineInfo(const LineContainer & r_data, int my_id) noexcept;
};

/**
This class is a container for all the powerlines on the grid.

**/
class LineContainer : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>, public IteratorAdder<LineContainer, LineInfo>
{
    friend class LineInfo;
    public:
        typedef LineInfo DataInfo;

    public:
        typedef std::tuple<
                   TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::StateRes
                   >  StateRes;
        
        LineContainer() noexcept = default;
        virtual ~LineContainer() noexcept = default;
        
        void init(const RealVect & branch_r,
                  const RealVect & branch_x,
                  const CplxVect & branch_h,
                  const Eigen::VectorXi & branch_from_id,
                  const Eigen::VectorXi & branch_to_id
                  );
              
        void init(const RealVect & branch_r,
                  const RealVect & branch_x,
                  const CplxVect & branch_h_or,
                  const CplxVect & branch_h_ex,
                  const Eigen::VectorXi & branch_from_id,
                  const Eigen::VectorXi & branch_to_id
                  );
              
        // pickle
        StateRes get_state() const
        {
            StateRes res(get_tsc_rxha_state());
            return res;
        }
        void set_state(LineContainer::StateRes & my_state )
        {
            set_tsc_rxha_state(std::get<0>(my_state));
            _update_model_coeffs();
            reset_results();
        }
        
        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const std::vector<SolverBusId> & id_grid_to_solver,
                             const RealVect & bus_vn_kv,
                             real_type sn_mva,
                             bool ac){
            compute_results_tsc_rxha(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        }

        void reset_results() {reset_results_tsc_rxha();}

        // for consistency with trafo, when used for example in BaseMultiplePowerflow...
        Eigen::Ref<const RealVect> dc_x_tau_shift() const {return RealVect();}
    
    protected:
        void _update_model_coeffs();

    protected:
        // physical properties

        // specific grid2op

        // input data

        //output data

        // model coefficients
};

inline LineInfo::LineInfo(const LineContainer & r_data, int my_id) noexcept:
TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::TwoSidesContainer_rxh_AInfo(r_data, my_id) {}

#endif  //LINE_CONTAINER_H
