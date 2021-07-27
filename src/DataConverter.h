// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATACONVERTER_H
#define DATACONVERTER_H
#include <vector>

#include "Utils.h"

#include "Eigen/Core"
#include "Eigen/Dense"

#include "BaseConstants.h"
/**
This class and all that are in there are provided as examples.

It allows conversion of "higher level" data, such as pandapower into a format that can be
digested by DataModel further used to compute powerflows thanks to KLUSolver.
**/

class PandaPowerConverter : public BaseConstants
{

    public:
        PandaPowerConverter():BaseConstants(),sn_mva_(-1.0),f_hz_(-1.0){};
        void set_f_hz(real_type f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(real_type sn_mva) { sn_mva_ = sn_mva;}

        // data converters
        /**
        This converts the trafo from pandapower to r, x and h (pair unit)
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect>
           get_trafo_param(const RealVect & tap_step_pct,
                           const RealVect & tap_pos,
                           const RealVect & tap_angles,
                           const std::vector<bool> & is_tap_hv_side,
                           const RealVect & vn_hv,  // nominal voltage of hv bus
                           const RealVect & vn_lv,  // nominal voltage of lv bus
                           const RealVect & trafo_vk_percent,
                           const RealVect & trafo_vkr_percent,
                           const RealVect & trafo_sn_trafo_mva,
                           const RealVect & trafo_pfe_kw,
                           const RealVect & trafo_i0_pct);
        /**
        pair unit properly the powerlines
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect>
           get_line_param(const RealVect & branch_r,
                          const RealVect & branch_x,
                          const RealVect & branch_c,
                          const RealVect & branch_g, //TODO g is not supported atm!
                          const RealVect & branch_from_kv,
                          const RealVect & branch_to_kv);

    protected:
        real_type sn_mva_;
        real_type f_hz_;

    private:
        void _check_init();
};

// TODO have a converter from ppc !

#endif // DATACONVERTER_H
