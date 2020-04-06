// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of PyKLU2Grid, PyKLU2Grid a implements a c++ backend targeting the Grid2Op platform.

#include "DataConverter.h"

void PandaPowerConverter::_check_init(){
    if(sn_mva_ <= 0.){
        throw std::runtime_error("PandaPowerConverter sn_mva has not been initialized");
    }
    if(f_hz_ <= 0.){
        throw std::runtime_error("PandaPowerConverter f_hz has not been initialized");
    }
}


std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXcd>
           PandaPowerConverter::get_trafo_param(const Eigen::VectorXd & trafo_vn_hv,
                                                const Eigen::VectorXd & trafo_vn_lv,
                                                const Eigen::VectorXd & trafo_vk_percent,
                                                const Eigen::VectorXd & trafo_vkr_percent,
                                                const Eigen::VectorXd & trafo_sn_trafo_mva,
                                                const Eigen::VectorXd & trafo_pfe_kw,
                                                const Eigen::VectorXd & trafo_i0_pct,
                                                const Eigen::VectorXd & trafo_lv_id_vn_kv)
{
    //TODO only for "trafo model = t"
    //TODO supposes that the step start at 0 for "no ratio"
    _check_init();

    //TODO consistency: move this class outside of here
    int nb_trafo = trafo_vn_lv.size();

    Eigen::VectorXd vn_trafo_lv = trafo_vn_lv;
    const Eigen::VectorXd & vn_lv = trafo_lv_id_vn_kv;

    // compute r and x
    Eigen::VectorXd tmp = vn_trafo_lv.array() / vn_lv.array();
    tmp = tmp.array() * tmp.array();
    Eigen::VectorXd tap_lv = tmp * sn_mva_;
    Eigen::VectorXd _1_sn_trafo_mva = 1.0 / trafo_sn_trafo_mva.array();
    Eigen::VectorXd z_sc = 0.01 * trafo_vk_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    Eigen::VectorXd r_sc = 0.01 * trafo_vkr_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    Eigen::VectorXd tmp2 = z_sc.array()*z_sc.array() - r_sc.array() * r_sc.array();
    Eigen::VectorXd x_sc = z_sc.cwiseSign().array() * tmp2.cwiseSqrt().array();

    // compute h, the subsceptance
    Eigen::VectorXd baseR = trafo_lv_id_vn_kv.array() * trafo_lv_id_vn_kv.array();
    baseR.array() /= sn_mva_;
    Eigen::VectorXd pfe =  trafo_pfe_kw.array() * 1e-3;

    // Calculate subsceptance ###
    Eigen::VectorXd vnl_squared = trafo_vn_lv.array() * trafo_vn_lv.array();
    Eigen::VectorXd b_real = pfe.array() / vnl_squared.array() * baseR.array();
    tmp2 = (trafo_i0_pct.array() * 0.01 * trafo_sn_trafo_mva.array());
    Eigen::VectorXd b_img =  tmp2.array() * tmp2.array() - pfe.array() * pfe.array();

    for(int i = 0; i<nb_trafo; ++i) {if (b_img(i) < 0.)  b_img(i) = 0.;}
    b_img = b_img.cwiseSqrt();
    b_img.array() *= baseR.array() / vnl_squared.array();
    Eigen::VectorXcd y = - 1.0i * b_real.array().cast<cdouble>() - b_img.array().cast<cdouble>() * trafo_i0_pct.cwiseSign().array();
    Eigen::VectorXcd b_sc = y.array() / tmp.array().cast<cdouble>();

    //transform trafo from t model to pi model, of course...
    // (remove that if trafo model is not t, but directly pi)
    cdouble my_i = 1.0i;
    for(int i = 0; i<nb_trafo; ++i){
        if(b_sc(i) == 0.) continue;
        cdouble za_star = 0.5 * (r_sc(i) + my_i * x_sc(i));
        cdouble zc_star = - my_i / b_sc(i);
        cdouble zSum_triangle = za_star * za_star + 2.0 * za_star * zc_star;
        cdouble zab_triangle = zSum_triangle / zc_star;
        cdouble zbc_triangle = zSum_triangle / za_star;

        r_sc(i) = zab_triangle.real();
        x_sc(i) = zab_triangle.imag();
        b_sc(i) = -2.0 * my_i / zbc_triangle;
    }

    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd> res =
        std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd>(std::move(r_sc), std::move(x_sc), std::move(b_sc));
    return res;
}

std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXcd>
           PandaPowerConverter::get_line_param(const Eigen::VectorXd & branch_r,
                                               const Eigen::VectorXd & branch_x,
                                               const Eigen::VectorXd & branch_c,
                                               const Eigen::VectorXd & branch_g,
                                               const Eigen::VectorXd & branch_from_kv,
                                               const Eigen::VectorXd & branch_to_kv)
{
    //TODO does not use c at the moment!
    _check_init();
    int nb_line = branch_r.size();
    Eigen::VectorXd branch_from_pu = branch_from_kv.array() * branch_from_kv.array() / sn_mva_;

    Eigen::VectorXd powerlines_r = branch_r.array() / branch_from_pu.array();
    Eigen::VectorXd powerlines_x = branch_x.array() / branch_from_pu.array();

    Eigen::VectorXcd powerlines_h = Eigen::VectorXcd::Constant(nb_line, 2.0 * f_hz_ * M_PI * 1e-9);
    powerlines_h.array() *= branch_c.array().cast<cdouble>();
    powerlines_h.array() *=  branch_from_pu.array().cast<cdouble>();
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd> res = std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXcd> (std::move(powerlines_r), std::move(powerlines_x), std::move(powerlines_h));
    return std::move(res);
}