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
#include "BaseSubstation.hpp"
#include "GenericContainer.hpp"
#include "TwoSidesContainer_rxh_A.hpp"

/**
This class is a container for all transformers on the grid.
Transformers are modeled "in pi" here. If your trafo are given in a "t" model (like in pandapower
for example) use the DataConverter class.

The convention used for the transformer is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html#electric-model
**/
class TrafoContainer : public TwoSidesContainer_rxh_A<OneSideContainer>
{
    //////////////////////////////
    // access data from base class
    public:
        using TwoSidesContainer_rxh_A<OneSideContainer>::get_buses_side_1;
        using TwoSidesContainer_rxh_A<OneSideContainer>::get_buses_side_2;

    protected:
        using TwoSidesContainer_rxh_A<OneSideContainer>::side_1_;
        using TwoSidesContainer_rxh_A<OneSideContainer>::side_2_;
        using TwoSidesContainer_rxh_A<OneSideContainer>::status_global_;
    //////////////////////////////

    public:
        class TrafoInfo : public TwoSidesContainer_rxh_A<OneSideContainer>::TwoSidesContainer_rxh_AInfo
        {
            public:
                // members
                // real_type r_pu;
                // real_type x_pu;
                // cplx_type h_pu;
                real_type ratio;
                real_type shift_rad;
                bool is_tap_hv_side;

                // bool has_res;
                // real_type res_a1_ka;
                // real_type res_a2_ka;

                TrafoInfo(const TrafoContainer & r_data_trafo, int my_id):
                TwoSidesContainer_rxh_AInfo(r_data_trafo, my_id),
                // r_pu(-1.0),
                // x_pu(-1.0),
                // h_pu(0., 0.),
                ratio(-1.0),
                shift_rad(-1.0),
                is_tap_hv_side(true)
                // has_res(false),
                // res_a1_ka(0.),
                // res_a2_ka(0.)
                {
                    if(my_id < 0) return;
                    if(my_id >= r_data_trafo.nb()) return;
                    is_tap_hv_side = r_data_trafo.is_tap_hv_side_[my_id];
                    ratio = r_data_trafo.ratio_.coeff(my_id);
                    shift_rad = r_data_trafo.shift_.coeff(my_id);
                }
        };
        typedef TrafoInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<TrafoContainer> TrafoContainerConstIterator;

    // make it iterable
    public:
        typedef TrafoContainerConstIterator const_iterator_type;
        const_iterator_type begin() const {return TrafoContainerConstIterator(this, 0); }
        const_iterator_type end() const {return TrafoContainerConstIterator(this, nb()); }
        TrafoInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a transformer with negative id");
            }
            if(id >= nb())
            {
                throw std::range_error("Trafo out of bound. Not enough transformers on the grid.");
            }
            return TrafoInfo(*this, id);
        }
        /////////////////////////

    public:
        typedef std::tuple<
                   TwoSidesContainer_rxh_A<OneSideContainer>::StateRes,
                //    std::vector<real_type>, // branch_r
                //    std::vector<real_type>, // branch_x
                //    std::vector<cplx_type >, // branch_h
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

        // method used within lightsim
        // void deactivate(int trafo_id, SolverControl & solver_control) {
        //     if(status_global_[trafo_id]){
        //         solver_control.tell_recompute_ybus();
        //         // but sparsity pattern do not change here (possibly one more coeff at 0.)
        //         solver_control.tell_ybus_some_coeffs_zero();
        //         solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        //     }
        //     _generic_deactivate(trafo_id, status_global_);
        // }
        // void reactivate(int trafo_id, SolverControl & solver_control) {
        //     if(!status_global_[trafo_id]){
        //         solver_control.tell_recompute_ybus();
        //         solver_control.tell_ybus_change_sparsity_pattern();  // this might change
        //         solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        //     }
        //     _generic_reactivate(trafo_id, status_global_);
        // }
        // void change_bus_hv(int trafo_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        //     _generic_change_bus(trafo_id, new_bus_id, get_buses_not_const_side_1(), solver_control, nb_bus);
        // }
        // void change_bus_lv(int trafo_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        //     _generic_change_bus(trafo_id, new_bus_id, get_buses_not_const_side_2(), solver_control, nb_bus);
        // }
        // void reconnect_connected_buses(Substation & substation) const;
        // virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component);
        
        // virtual void nb_line_end(std::vector<int> & res) const;
        // virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const;
        
        // virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
        //                       bool ac,
        //                       const std::vector<int> & id_grid_to_solver,
        //                       real_type sn_mva) const;
        virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<int> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const;
        virtual void fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & Bf,
                                     const std::vector<int> & id_grid_to_solver,
                                     real_type sn_mva,
                                     int nb_powerline,
                                     bool transpose) const;
        virtual void hack_Sbus_for_dc_phase_shifter(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);  // needed for dc mode  

        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const std::vector<int> & id_grid_to_solver,
                             const RealVect & bus_vn_kv,
                             real_type sn_mva,
                             bool ac)
        {
            compute_results_tsc_rxha(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
            if(!ac){
                Eigen::Ref<RealVect> res_p_side_1 = get_res_p_side_1();
                Eigen::Ref<RealVect> res_p_side_2 = get_res_p_side_2();

                const int nb_element = nb();
                for(int el_id = 0; el_id < nb_element; ++el_id){
                    res_p_side_1(el_id) -= dc_x_tau_shift_(el_id) * sn_mva;
                    res_p_side_2(el_id) += dc_x_tau_shift_(el_id) * sn_mva;
                }
            }
        }
        void reset_results(){
            reset_results_tsc_rxha();
        }

        tuple4d get_res_hv() const {
            return TwoSidesContainer_rxh_A<OneSideContainer>::get_res_side_1();
        }
        tuple4d get_res_lv() const {
            return TwoSidesContainer_rxh_A<OneSideContainer>::get_res_side_2();
        }
        tuple5d get_res_hv_full() const {
            return TwoSidesContainer_rxh_A<OneSideContainer>::get_res_side_1_full();
        }
        tuple5d get_res_lv_full() const {
            return TwoSidesContainer_rxh_A<OneSideContainer>::get_res_side_2_full();
        }

        // model paramters
        // Eigen::Ref<const CplxVect> yac_ff() const {return yac_ff_;}
        // Eigen::Ref<const CplxVect> yac_ft() const {return yac_ft_;}
        // Eigen::Ref<const CplxVect> yac_tf() const {return yac_tf_;}
        // Eigen::Ref<const CplxVect> yac_tt() const {return yac_tt_;}

        // Eigen::Ref<const CplxVect> ydc_ff() const {return ydc_ff_;}
        // Eigen::Ref<const CplxVect> ydc_ft() const {return ydc_ft_;}
        // Eigen::Ref<const CplxVect> ydc_tf() const {return ydc_tf_;}
        // Eigen::Ref<const CplxVect> ydc_tt() const {return ydc_tt_;}
        Eigen::Ref<const RealVect> dc_x_tau_shift() const {return dc_x_tau_shift_;}

        // TODO !
        // for batched algorithm (need to be removed when powerlines will accept same API)
        const std::vector<bool>& get_status() const {return get_status_global();}
        Eigen::Ref<const Eigen::VectorXi> get_bus_from() const {return get_bus_id_side_1();}
        Eigen::Ref<const Eigen::VectorXi> get_bus_to() const {return get_bus_id_side_2();}
        
    protected:
        void _update_model_coeffs();
        
    protected:
        // physical properties
        // RealVect r_;
        // RealVect x_;
        // CplxVect h_;
        std::vector<bool> is_tap_hv_side_;  // whether the tap is hav side or not

        // input data
        RealVect ratio_;  // transformer ratio
        RealVect shift_;  // phase shifter (in radian !)

        //output data
        // RealVect res_a_hv_;  // in kA
        // RealVect res_a_lv_;  // in kA

        // model coefficients
        // CplxVect yac_ff_;
        // CplxVect yac_ft_;
        // CplxVect yac_tf_;
        // CplxVect yac_tt_;

        // CplxVect ydc_ff_;
        // CplxVect ydc_ft_;
        // CplxVect ydc_tf_;
        // CplxVect ydc_tt_;
        RealVect dc_x_tau_shift_;
};

#endif  //TRAFO_CONTAINER_H
