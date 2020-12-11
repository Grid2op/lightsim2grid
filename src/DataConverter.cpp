// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataConverter.h"

void PandaPowerConverter::_check_init(){
    if(sn_mva_ <= 0.){
        throw std::runtime_error("PandaPowerConverter sn_mva has not been initialized");
    }
    if(f_hz_ <= 0.){
        throw std::runtime_error("PandaPowerConverter f_hz has not been initialized");
    }
}

std::tuple<RealVect,
           RealVect,
           CplxVect>
           PandaPowerConverter::get_trafo_param(const RealVect & trafo_vn_hv,
                                                const RealVect & trafo_vn_lv,
                                                const RealVect & trafo_vk_percent,
                                                const RealVect & trafo_vkr_percent,
                                                const RealVect & trafo_sn_trafo_mva,
                                                const RealVect & trafo_pfe_kw,
                                                const RealVect & trafo_i0_pct,
                                                const RealVect & trafo_lv_id_vn_kv)
{
    //TODO only for "trafo model = t"
    //TODO supposes that the step start at 0 for "no ratio"
    _check_init();

    //TODO consistency: move this class outside of here
    int nb_trafo = trafo_vn_lv.size();

    RealVect vn_trafo_lv = trafo_vn_lv;
    const RealVect & vn_lv = trafo_lv_id_vn_kv;

    // compute r and x
    RealVect tmp = vn_trafo_lv.array() / vn_lv.array();
    tmp = tmp.array() * tmp.array();
    RealVect tap_lv = tmp * sn_mva_;
    RealVect _1_sn_trafo_mva = my_one_ / trafo_sn_trafo_mva.array();
    RealVect z_sc = 0.01 * trafo_vk_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    RealVect r_sc = 0.01 * trafo_vkr_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    RealVect tmp2 = z_sc.array()*z_sc.array() - r_sc.array() * r_sc.array();
    RealVect x_sc = z_sc.cwiseSign().array() * tmp2.cwiseSqrt().array();

    // compute h, the subsceptance
    RealVect baseR = trafo_lv_id_vn_kv.array() * trafo_lv_id_vn_kv.array();
    baseR.array() /= sn_mva_;
    RealVect pfe =  trafo_pfe_kw.array() * 1e-3;

    // Calculate subsceptance ###
    RealVect vnl_squared = trafo_vn_lv.array() * trafo_vn_lv.array();
    RealVect b_real = pfe.array() / vnl_squared.array() * baseR.array();
    tmp2 = (trafo_i0_pct.array() * 0.01 * trafo_sn_trafo_mva.array());
    RealVect b_img =  tmp2.array() * tmp2.array() - pfe.array() * pfe.array();

    for(int i = 0; i<nb_trafo; ++i) {if (b_img(i) < 0.)  b_img(i) = 0.;}
    b_img = b_img.cwiseSqrt();
    b_img.array() *= baseR.array() / vnl_squared.array();
    CplxVect y = - my_i * b_real.array().cast<cplx_type>() - b_img.array().cast<cplx_type>() * trafo_i0_pct.cwiseSign().array();
    CplxVect b_sc = y.array() / tmp.array().cast<cplx_type>();

    //transform trafo from t model to pi model, of course...
    // (remove that if trafo model is not t, but directly pi)
    for(int i = 0; i<nb_trafo; ++i){
        if(b_sc(i) == my_zero_) continue;
        cplx_type za_star = my_half_ * (r_sc(i) + my_i * x_sc(i));
        cplx_type zc_star = - my_i / b_sc(i);
        cplx_type zSum_triangle = za_star * za_star + my_two_ * za_star * zc_star;
        cplx_type zab_triangle = zSum_triangle / zc_star;
        cplx_type zbc_triangle = zSum_triangle / za_star;

        r_sc(i) = zab_triangle.real();
        x_sc(i) = zab_triangle.imag();
        b_sc(i) = -my_two_ * my_i / zbc_triangle;
    }

    std::tuple<RealVect, RealVect, CplxVect> res =
        std::tuple<RealVect, RealVect, CplxVect>(std::move(r_sc), std::move(x_sc), std::move(b_sc));
    return res;
}

std::tuple<RealVect,
           RealVect,
           CplxVect>
           PandaPowerConverter::get_line_param(const RealVect & branch_r,
                                               const RealVect & branch_x,
                                               const RealVect & branch_c,
                                               const RealVect & branch_g,
                                               const RealVect & branch_from_kv,
                                               const RealVect & branch_to_kv)
{
    //TODO does not use c at the moment!
    _check_init();
    int nb_line = branch_r.size();
    RealVect branch_from_pu = branch_from_kv.array() * branch_from_kv.array() / sn_mva_;

    RealVect powerlines_r = branch_r.array() / branch_from_pu.array();
    RealVect powerlines_x = branch_x.array() / branch_from_pu.array();

    CplxVect powerlines_h = CplxVect::Constant(nb_line, 2.0 * f_hz_ * M_PI * 1e-9);
    powerlines_h.array() *= branch_c.array().cast<cplx_type>();
    powerlines_h.array() *=  branch_from_pu.array().cast<cplx_type>();
    std::tuple<RealVect, RealVect, CplxVect> res = std::tuple<RealVect,
           RealVect,
           CplxVect> (std::move(powerlines_r), std::move(powerlines_x), std::move(powerlines_h));
    return res;
}
