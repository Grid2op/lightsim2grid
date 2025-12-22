// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TRAFO_CONTAINER_H
#define TRAFO_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SubstationContainer.hpp"
#include "OneSideContainer_forBranch.hpp"
#include "TwoSidesContainer_rxh_A.hpp"

class TrafoContainer;
class TrafoInfo : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::TwoSidesContainer_rxh_AInfo
{
    public:
        // members
        real_type ratio;
        real_type shift_rad;
        bool is_tap_hv_side;

        inline TrafoInfo(const TrafoContainer & r_data_trafo, int my_id);
};

/**
This class is a container for all transformers on the grid.
Transformers are modeled "in pi" here. If your trafo are given in a "t" model (like in pandapower
for example) use the DataConverter class.

The convention used for the transformer is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html#electric-model
**/
class TrafoContainer : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>, public IteratorAdder<TrafoContainer, TrafoInfo>
{
    //////////////////////////////
    // access data from base class
    public:
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::get_buses_side_1;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::get_buses_side_2;

    protected:
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::side_1_;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::side_2_;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::status_global_;
    //////////////////////////////

    friend class TrafoInfo;

    public:
        typedef TrafoInfo DataInfo;

    public:
        typedef std::tuple<
                   TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::StateRes,
                   std::vector<real_type>, // ratio_
                   std::vector<bool> , // is_tap_hv_side
                   std::vector<real_type> // shift_
               >  StateRes;

        TrafoContainer() {};
        virtual ~TrafoContainer() noexcept = default;

        void init(const RealVect & trafo_r,
                  const RealVect & trafo_x,
                  const CplxVect & trafo_b,
                  const RealVect & trafo_tap_step_pct,
                  const RealVect & trafo_tap_pos,
                  const RealVect & trafo_shift_degree,
                  const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                  const Eigen::VectorXi & trafo_hv_id,
                  const Eigen::VectorXi & trafo_lv_id
                  );

        void init(const RealVect & trafo_r,
                  const RealVect & trafo_x,
                  const CplxVect & trafo_b,
                  const RealVect & trafo_ratio,
                  const RealVect & trafo_shift_degree,
                  const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                  const Eigen::VectorXi & trafo_hv_id,
                  const Eigen::VectorXi & trafo_lv_id
                  );

        //pickle
        StateRes get_state() const;
        void set_state(StateRes & my_state );

        virtual void hack_Sbus_for_dc_phase_shifter(CplxVect & Sbus, bool ac, const std::vector<SolverBusId> & id_grid_to_solver);  // needed for dc mode  

        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const std::vector<SolverBusId> & id_grid_to_solver,
                             const RealVect & bus_vn_kv,
                             real_type sn_mva,
                             bool ac)
        {
            // compute base values
            compute_results_tsc_rxha_no_amps(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
            // adjust for phase shifters
            if(!ac){
                Eigen::Ref<RealVect> res_p_side_1 = get_res_p_side_1();
                Eigen::Ref<RealVect> res_p_side_2 = get_res_p_side_2();

                const int nb_element = nb();
                for(int el_id = 0; el_id < nb_element; ++el_id){
                    res_p_side_1(el_id) -= dc_x_tau_shift_(el_id) * sn_mva;
                    res_p_side_2(el_id) += dc_x_tau_shift_(el_id) * sn_mva;
                }
            }
            // compute amps flow
            compute_amps_after_all_set();
        }
        
        void reset_results(){
            reset_results_tsc_rxha();
        }

        Eigen::Ref<const RealVect> dc_x_tau_shift() const {return dc_x_tau_shift_;}

        void change_ratio(
            int el_id,
            real_type new_ratio,
            SolverControl & solver_control){
                if(std::abs(ratio_(el_id) - new_ratio) >_tol_equal_float){
                    ratio_(el_id) = new_ratio;
                    // TODO speed: only some part needs to be recomputed
                    _update_model_coeffs_one_el(el_id); 
                    solver_control.tell_recompute_ybus();
                }
        }
        
        void change_shift(
            int el_id,
            real_type new_shift,
            SolverControl & solver_control){
                if(std::abs(shift_(el_id) - new_shift) >_tol_equal_float){
                    shift_(el_id) = new_shift;
                    // TODO speed: only some part needs to be recomputed
                    _update_model_coeffs_one_el(el_id); 
                    solver_control.tell_recompute_ybus();
                    solver_control.tell_recompute_sbus();  // only in DC however
                }
        }
        
    protected:
        void _update_model_coeffs();
        void _update_model_coeffs_one_el(int el_id);

    protected:
        // physical properties
        std::vector<bool> is_tap_hv_side_;  // whether the tap is hav side or not

        // input data
        RealVect ratio_;  // transformer ratio
        RealVect shift_;  // phase shifter (in radian !)

        //output data

        // model coefficients
        RealVect dc_x_tau_shift_;

    protected:

        virtual real_type fillBf_for_PTDF_coeff(int el_id) const{
            real_type res = x_(el_id);
            real_type tau = is_tap_hv_side_[el_id] ? ratio_(el_id) : 1. / ratio_(el_id);
            return res * tau;
        }
        virtual FDPFCoeffs get_fdpf_coeffs(int tr_id, FDPFMethod xb_or_bx) const;

        virtual int fillBf_for_PTDF_id(int el_id, int nb_powerline) const{
            return el_id + nb_powerline;
        }
};

inline TrafoInfo::TrafoInfo(const TrafoContainer & r_data_trafo, int my_id):
TwoSidesContainer_rxh_AInfo(r_data_trafo, my_id),
ratio(-1.0),
shift_rad(-1.0),
is_tap_hv_side(true)
{
    if(my_id < 0) return;
    if(my_id >= r_data_trafo.nb()) return;
    is_tap_hv_side = r_data_trafo.is_tap_hv_side_[my_id];
    ratio = r_data_trafo.ratio_.coeff(my_id);
    shift_rad = r_data_trafo.shift_.coeff(my_id);
}

#endif  //TRAFO_CONTAINER_H
