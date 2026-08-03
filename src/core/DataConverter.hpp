// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATACONVERTER_H
#define DATACONVERTER_H
#include <string>
#include <vector>

#include "Utils.hpp"
#include "BaseConstants.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"

namespace ls2g {

/**
This class and all that are in there are provided as examples.

It allows conversion of "higher level" data, such as pandapower into a format that can be
digested by DataModel further used to compute powerflows thanks to KLUSolver.
**/
class LS2G_API PandaPowerConverter final : public BaseConstants
{

    public:
        PandaPowerConverter() noexcept :BaseConstants(),sn_mva_(-1.0),f_hz_(-1.0){};
        ~PandaPowerConverter() noexcept = default;
        
        void set_f_hz(real_type f_hz) { f_hz_ = f_hz;}
        void set_sn_mva(real_type sn_mva) { sn_mva_ = sn_mva;}

        // data converters
        /**
        This converts the trafo from pandapower to r, x and h (pair unit) (for legacy (<= 2.14.somthing) pandapower)
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect>
           get_trafo_param_pp3(const Eigen::Ref<const RealVect> & tap_step_pct,
                               const Eigen::Ref<const RealVect> & tap_pos,
                               const Eigen::Ref<const RealVect> & tap_angles,
                               const std::vector<bool> & is_tap_hv_side,
                               const Eigen::Ref<const RealVect> & vn_hv,  // nominal voltage of hv bus
                               const Eigen::Ref<const RealVect> & vn_lv,  // nominal voltage of lv bus
                               const Eigen::Ref<const RealVect> & trafo_vk_percent,
                               const Eigen::Ref<const RealVect> & trafo_vkr_percent,
                               const Eigen::Ref<const RealVect> & trafo_sn_trafo_mva,
                               const Eigen::Ref<const RealVect> & trafo_pfe_kw,
                               const Eigen::Ref<const RealVect> & trafo_i0_pct,
                               bool trafo_model_is_t) const;

        /**
        This converts the trafo from pandapower to r, x and h (pair unit) (for legacy (<= 2.14.somthing) pandapower)
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect>
           get_trafo_param_pp2(const Eigen::Ref<const RealVect> & tap_step_pct,
                               const Eigen::Ref<const RealVect> & tap_pos,
                               const Eigen::Ref<const RealVect> & tap_angles,
                               const std::vector<bool> & is_tap_hv_side,
                               const Eigen::Ref<const RealVect> & vn_hv,  // nominal voltage of hv bus
                               const Eigen::Ref<const RealVect> & vn_lv,  // nominal voltage of lv bus
                               const Eigen::Ref<const RealVect> & trafo_vk_percent,
                               const Eigen::Ref<const RealVect> & trafo_vkr_percent,
                               const Eigen::Ref<const RealVect> & trafo_sn_trafo_mva,
                               const Eigen::Ref<const RealVect> & trafo_pfe_kw,
                               const Eigen::Ref<const RealVect> & trafo_i0_pct) const;


        /**
        pair unit properly the powerlines (for legacy (<= 2.14.somthing) pandapower)
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect>
           get_line_param_legacy(const Eigen::Ref<const RealVect> & branch_r,
                                 const Eigen::Ref<const RealVect> & branch_x,
                                 const Eigen::Ref<const RealVect> & branch_g,
                                 const Eigen::Ref<const RealVect> & branch_c,
                                 const Eigen::Ref<const RealVect> & branch_from_kv,
                                 const Eigen::Ref<const RealVect> & branch_to_kv) const;
        /**
        pair unit properly the powerlines (for most recent pandapower)
        **/
        std::tuple<RealVect,
                   RealVect,
                   CplxVect,
                   CplxVect>
           get_line_param(const Eigen::Ref<const RealVect> & branch_r,
                          const Eigen::Ref<const RealVect> & branch_x,
                          const Eigen::Ref<const RealVect> & branch_g,
                          const Eigen::Ref<const RealVect> & branch_c,
                          const Eigen::Ref<const RealVect> & branch_from_kv,
                          const Eigen::Ref<const RealVect> & branch_to_kv) const;

    private:
        real_type sn_mva_;
        real_type f_hz_;

    private:
        void _check_init() const;

        // Every method above walks its inputs together, element by element, with
        // unchecked accessors (`vect.coeff(i)`, `is_tap_hv_side[i]`) and combines
        // them in coefficient-wise Eigen expressions whose result is sized from
        // one of them. Eigen's own size assertions are compiled out of release
        // wheels (-O3 -DNDEBUG), so a caller passing arrays of different lengths
        // reads AND writes past the end of the shorter ones (observed as heap
        // corruption, not a clean error). These are public python bindings
        // (`PandaPowerConverter`), so the sizes have to be checked here.
        static void _check_same_size(const std::string & fun_name,
                                     Eigen::Index expected,
                                     const std::string & expected_name,
                                     const std::string & name,
                                     Eigen::Index actual);
        // the shared input-shape check of get_trafo_param_pp2 / get_trafo_param_pp3
        // (identical argument lists)
        static void _check_trafo_inputs(const std::string & fun_name,
                                        const Eigen::Ref<const RealVect> & tap_step_pct,
                                        const Eigen::Ref<const RealVect> & tap_pos,
                                        const Eigen::Ref<const RealVect> & tap_angles,
                                        const std::vector<bool> & is_tap_hv_side,
                                        const Eigen::Ref<const RealVect> & vn_hv,
                                        const Eigen::Ref<const RealVect> & vn_lv,
                                        const Eigen::Ref<const RealVect> & trafo_vk_percent,
                                        const Eigen::Ref<const RealVect> & trafo_vkr_percent,
                                        const Eigen::Ref<const RealVect> & trafo_sn_trafo_mva,
                                        const Eigen::Ref<const RealVect> & trafo_pfe_kw,
                                        const Eigen::Ref<const RealVect> & trafo_i0_pct);
        // ... and of get_line_param / get_line_param_legacy
        static void _check_line_inputs(const std::string & fun_name,
                                       const Eigen::Ref<const RealVect> & branch_r,
                                       const Eigen::Ref<const RealVect> & branch_x,
                                       const Eigen::Ref<const RealVect> & branch_g,
                                       const Eigen::Ref<const RealVect> & branch_c,
                                       const Eigen::Ref<const RealVect> & branch_from_kv,
                                       const Eigen::Ref<const RealVect> & branch_to_kv);
};

// TODO have a converter from ppc !


} // namespace ls2g

#endif // DATACONVERTER_H