// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef INJ_ACTION_H
#define INJ_ACTION_H

#include "LSGrid.hpp"

#include <vector>
#include <cmath>

namespace ls2g {


class InjAction
{

    public:
        // TODO
        InjAction(const Eigen::Ref<const RealVect> & load_p,
                  const Eigen::Ref<const RealVect> & load_q,
                  const Eigen::Ref<const RealVect> & gen_p,
                  const Eigen::Ref<const RealVect> & gen_v,
                  const Eigen::Ref<const RealVect> & storage_p,
                  const Eigen::Ref<const RealVect> & shunt_p,
                  const Eigen::Ref<const RealVect> & shunt_q,
                  const Eigen::Ref<const RealVect> & sgen_p,
                  const Eigen::Ref<const RealVect> & sgen_q
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
        void apply_to_gridmodel(LSGrid & grid) const {
            update_grid_values(grid, load_p_,  &LSGrid::change_p_load);
            update_grid_values(grid, load_q_,  &LSGrid::change_q_load);
            update_grid_values(grid, gen_p_,  &LSGrid::change_p_gen);
            update_grid_values(grid, gen_v_,  &LSGrid::change_v_gen);
            // TODO rest !
        }

        // TODO check compliance (correct number of elements etc.)
        void check_validity(const LSGrid & grid) const {

        }
        
        template<class EigenType, class FunctorType>
        void update_grid_values(LSGrid & grid,
                                const EigenType & new_values,
                                FunctorType fun) const
        {
            for(int el_id = 0; el_id < new_values.rows(); ++el_id)
            {
                auto tmp = new_values(el_id);
                if(std::isfinite(tmp))
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

} // namespace ls2g

#endif // INJ_ACTION_H
