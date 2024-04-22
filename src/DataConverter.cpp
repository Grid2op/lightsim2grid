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
        throw std::runtime_error("PandaPowerConverter::_check_init: sn_mva has not been initialized");
    }
    if(f_hz_ <= 0.){
        throw std::runtime_error("PandaPowerConverter::_check_init: f_hz has not been initialized");
    }
}

std::tuple<RealVect,
           RealVect,
           CplxVect>
           PandaPowerConverter::get_trafo_param(const RealVect & tap_step_pct,
                                                const RealVect & tap_pos,
                                                const RealVect & tap_angles,
                                                const std::vector<bool> & is_tap_hv_side,
                                                const RealVect & vn_hv,  // nominal voltage of hv bus
                                                const RealVect & vn_lv,  // nominal voltage of lv bus
                                                const RealVect & trafo_vk_percent,
                                                const RealVect & trafo_vkr_percent,
                                                const RealVect & trafo_sn_trafo_mva,
                                                const RealVect & trafo_pfe_kw,
                                                const RealVect & trafo_i0_pct)
{
    //TODO consistency: move this class outside of here
    _check_init();

    const int nb_trafo = static_cast<int>(tap_step_pct.size());
    // TODO check all vectors have the same size

    // compute the adjusted for phase shifter and tap side
    auto tap_steps = 0.01 * tap_step_pct.array() * tap_pos.array();
    auto du_hv = vn_hv.array() * tap_steps.array();
    auto du_lv = vn_lv.array() * tap_steps.array();
    RealVect trafo_vn_hv = vn_hv;
    RealVect trafo_vn_lv = vn_lv;
    for(int i = 0; i < nb_trafo; ++i)
    {
        if(is_tap_hv_side[i]) {
            // adjust the voltage hv side
            double tmp_cos = vn_hv.coeff(i) + du_hv.coeff(i) * cos(tap_angles.coeff(i));
            double tmp_sin = du_hv.coeff(i) * sin(tap_angles.coeff(i));
            trafo_vn_hv.coeffRef(i) = sqrt(tmp_cos * tmp_cos + tmp_sin * tmp_sin);
        }else{
            // adjust the voltage lv side
            double tmp_cos = vn_lv.coeff(i) + du_lv.coeff(i) * cos(tap_angles.coeff(i));
            double tmp_sin = du_lv.coeff(i) * sin(tap_angles.coeff(i));
            trafo_vn_lv.coeffRef(i) = sqrt(tmp_cos * tmp_cos + tmp_sin * tmp_sin);
        }
    }
    const RealVect & vn_trafo_lv = trafo_vn_lv;

    // compute r and x
    // tap_lv = np.square(vn_trafo_lv / vn_lv) * sn_mva  # adjust for low voltage side voltage converter
    RealVect vn_trafo_lv_1_vn_lv_sq = vn_trafo_lv.array() / vn_lv.array();
    vn_trafo_lv_1_vn_lv_sq = vn_trafo_lv_1_vn_lv_sq.array() * vn_trafo_lv_1_vn_lv_sq.array();
    RealVect tap_lv = vn_trafo_lv_1_vn_lv_sq * sn_mva_; // tap_lv = np.square(vn_trafo_lv / vn_lv) * sn_mva

    // z_sc = vk_percent / 100. / sn_trafo_mva * tap_lv
    RealVect _1_sn_trafo_mva = my_one_ / trafo_sn_trafo_mva.array();
    RealVect z_sc = 0.01 * trafo_vk_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    // r_sc = vkr_percent / 100. / sn_trafo_mva * tap_lv
    RealVect r_sc = 0.01 * trafo_vkr_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    // x_sc = np.sign(z_sc) * np.sqrt(z_sc ** 2 - r_sc ** 2)
    RealVect tmp2 = z_sc.array()*z_sc.array() - r_sc.array() * r_sc.array();
    RealVect x_sc = z_sc.cwiseSign().array() * tmp2.cwiseSqrt().array();

    // compute h, the subsceptance
    // baseR = np.square(vn_lv) / sn_mva
    RealVect baseR = vn_lv.array() * vn_lv.array();
    baseR.array() /= sn_mva_;
    // pfe = get_trafo_values(trafo_df, "pfe_kw") * 1e-3
    RealVect pfe =  trafo_pfe_kw.array() * 1e-3;

    // Calculate subsceptance ###
    // vnl_squared = vn_lv_kv ** 2
    RealVect vnl_squared = vn_lv.array() * vn_lv.array();
    // b_real = pfe / vnl_squared * baseR
    RealVect b_real = pfe.array() * baseR.array() / vnl_squared.array();
    // b_img = (i0 / 100. * sn) ** 2 - pfe ** 2
    tmp2 = (trafo_i0_pct.array() * 0.01 * trafo_sn_trafo_mva.array());
    RealVect b_img =  tmp2.array() * tmp2.array() - pfe.array() * pfe.array();

    // b_img[b_img < 0] = 0
    for(int i = 0; i<nb_trafo; ++i) {if (b_img(i) < 0.)  b_img(i) = 0.;}
    //  b_img = np.sqrt(b_img) * baseR / vnl_squared
    b_img = b_img.cwiseSqrt();
    b_img.array() *= baseR.array() / vnl_squared.array();
    // y = - b_real * 1j - b_img * np.sign(i0)
    CplxVect y = - my_i * b_real.array().cast<cplx_type>() - b_img.array().cast<cplx_type>() * trafo_i0_pct.cwiseSign().array();
    // return y / np.square(vn_trafo_lv / vn_lv_kv) * parallel
    CplxVect b_sc = y.array() / vn_trafo_lv_1_vn_lv_sq.array().cast<cplx_type>();


    //transform trafo from t model to pi model, of course...
    // (remove that if trafo model is not t, but directly pi)
    for(int i = 0; i<nb_trafo; ++i){
        if(b_sc(i) == my_zero_) continue;
        cplx_type za_star = my_half_ * (r_sc(i) + my_i * x_sc(i));
        cplx_type zc_star = - my_i / b_sc(i);
        cplx_type zSum_triangle = za_star * za_star + my_two_ * za_star * zc_star;
        cplx_type zab_triangle = zSum_triangle / zc_star;
        cplx_type zbc_triangle = zSum_triangle / za_star;

        r_sc(i) = std::real(zab_triangle);
        x_sc(i) = std::imag(zab_triangle);
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
                                               const RealVect & branch_g,  // TODO
                                               const RealVect & branch_from_kv,
                                               const RealVect & branch_to_kv)
{
    //TODO does not use c at the moment!
    _check_init();
    const int nb_line = static_cast<int>(branch_r.size());
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
