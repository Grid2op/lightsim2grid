// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATACONVERTER_H
#define DATACONVERTER_H

#include "Eigen/Core"
#include "Eigen/Dense"

#include "Utils.h"
/**
This class and all that are in there are provided as examples.

It allows conversion of "higher level" data, such as pandapower into a format that can be
digested by DataModel further used to compute powerflows thanks to KLUSolver.
**/

class PandaPowerConverter{

    public:
        PandaPowerConverter():sn_mva_(-1.0),f_hz_(-1.0){};
        void set_f_hz(double f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(double sn_mva) { sn_mva_ = sn_mva;}

        // data converters
        /**
        This converts the trafo from pandapower to r, x and h (pair unit)
        **/
        std::tuple<Eigen::VectorXd,
                   Eigen::VectorXd,
                   Eigen::VectorXcd>
           get_trafo_param(const Eigen::VectorXd & trafo_vn_hv,
                           const Eigen::VectorXd & trafo_vn_lv,
                           const Eigen::VectorXd & trafo_vk_percent,
                           const Eigen::VectorXd & trafo_vkr_percent,
                           const Eigen::VectorXd & trafo_sn_trafo_mva,
                           const Eigen::VectorXd & trafo_pfe_kw,
                           const Eigen::VectorXd & trafo_i0_pct,
                           const Eigen::VectorXd & trafo_lv_id_vn_kv);
        /**
        pair unit properly the powerlines
        **/
        std::tuple<Eigen::VectorXd,
                   Eigen::VectorXd,
                   Eigen::VectorXcd>
           get_line_param(const Eigen::VectorXd & branch_r,
                          const Eigen::VectorXd & branch_x,
                          const Eigen::VectorXd & branch_c,
                          const Eigen::VectorXd & branch_g, //TODO g is not supported atm!
                          const Eigen::VectorXd & branch_from_kv,
                          const Eigen::VectorXd & branch_to_kv);

    protected:
        double sn_mva_;
        double f_hz_;
        static const cdouble my_i;

    private:
        void _check_init();
};

// TODO have a converter from ppc !

#endif // DATACONVERTER_H
