// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef INJ_ACTION_H
#define INJ_ACTION_H

#include "GridModel.hpp"

#include <vector>


class InjAction
{

    public:
        // TODO
        InjAction(const RealVect & load_p,
                  const RealVect & load_q,
                  const RealVect & gen_p,
                  const RealVect & gen_v,
                  const RealVect & storage_p,
                  const RealVect & shunt_p,
                  const RealVect & shunt_q,
                  const RealVect & sgen_p,
                  const RealVect & sgen_q
                  ){
            load_p_ = load_p;
            load_q_ = load_q;
            gen_p_ = gen_p;
            gen_v_ = gen_v;
            storage_p_ = storage_p;
            shunt_p_ = shunt_p;
            shunt_q_ = shunt_q;
            sgen_p_ = sgen_p;
            sgen_q_ = sgen_q;
        }

        // TODO
        void apply_to_gridmodel(GridModel & grid) const {
            update_grid_values(grid, load_p_,  &GridModel::change_p_load);
            update_grid_values(grid, load_q_,  &GridModel::change_q_load);
            update_grid_values(grid, gen_p_,  &GridModel::change_p_gen);
            update_grid_values(grid, gen_v_,  &GridModel::change_v_gen);
            // TODO rest !
        }

        // TODO check compliance (correct number of elements etc.)
        void check_validity(const GridModel & grid) const {

        }
        
        template<class EigenType, class FunctorType>
        void update_grid_values(GridModel & grid,
                                const EigenType & new_values,
                                FunctorType fun) const
        {
            for(int el_id = 0; el_id < new_values.rows(); ++el_id)
            {
                auto tmp = new_values(el_id);
                if(isfinite(tmp))
                {
                    (grid.*fun)(el_id, static_cast<real_type>(tmp));
                }
            }
        }
    protected:
    
        RealVect load_p_;
        RealVect load_q_;
        RealVect gen_p_;
        RealVect gen_v_;
        RealVect storage_p_;
        RealVect shunt_p_;
        RealVect shunt_q_;
        RealVect sgen_p_;
        RealVect sgen_q_;
        
};

#endif // INJ_ACTION_H
