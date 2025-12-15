// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TWO_SIDES_CONTAINER_RXH_A_H
#define TWO_SIDES_CONTAINER_RXH_A_H

#include "TwoSidesContainer.hpp"
/**
 * Type of container to represent a line or a transformer.
 * 
 * It has results in amps (A), and some physical properties (r, x and h = g+j.b)
 */
template<class OneSideType>
class TwoSidesContainer_rxh_A: public TwoSidesContainer<OneSideType>
{
    //////////////////////////////
    // access data from base class
    public:
        using TwoSidesContainer<OneSideContainer>::_deactivated_bus_id;
        using TwoSidesContainer<OneSideContainer>::v_disco_el_;
        using TwoSidesContainer<OneSideContainer>::theta_disco_el_;
        using TwoSidesContainer<OneSideContainer>::my_180_pi_;
        using TwoSidesContainer<OneSideContainer>::_generic_deactivate;
        using TwoSidesContainer<OneSideContainer>::_generic_reactivate;

        using TwoSidesContainer<OneSideContainer>::nb;
        using TwoSidesContainer<OneSideContainer>::reset_results_tsc;
        using TwoSidesContainer<OneSideContainer>::check_size;
        using TwoSidesContainer<OneSideContainer>::get_bus_side_1;
        using TwoSidesContainer<OneSideContainer>::get_bus_side_2;

    protected:
        using TwoSidesContainer<OneSideContainer>::_get_amps;

        using TwoSidesContainer<OneSideContainer>::get_tsc_state;
        using TwoSidesContainer<OneSideContainer>::set_tsc_state;
        using TwoSidesContainer<OneSideContainer>::status_global_;
        using TwoSidesContainer<OneSideContainer>::side_1_;
        using TwoSidesContainer<OneSideContainer>::side_2_;
        using TwoSidesContainer<OneSideContainer>::get_res_p_side_1;
        using TwoSidesContainer<OneSideContainer>::get_res_p_side_2;
        using TwoSidesContainer<OneSideContainer>::get_res_q_side_1;
        using TwoSidesContainer<OneSideContainer>::get_res_q_side_2;
        using TwoSidesContainer<OneSideContainer>::get_res_v_side_1;
        using TwoSidesContainer<OneSideContainer>::get_res_v_side_2;
        using TwoSidesContainer<OneSideContainer>::get_res_theta_side_1;
        using TwoSidesContainer<OneSideContainer>::get_res_theta_side_2;
    //////////////////////////////

    public:
        class TwoSidesContainer_rxh_AInfo : public TwoSidesContainer<OneSideContainer>::TwoSidesInfo
        {
            public:
                // members
                real_type r_pu;
                real_type x_pu;
                cplx_type h_pu;

                bool has_res;
                real_type res_a1_ka;
                real_type res_a2_ka;

                cplx_type yac_ff;
                cplx_type yac_ft;
                cplx_type yac_tf;
                cplx_type yac_tt;
                cplx_type ydc_ff;
                cplx_type ydc_ft;
                cplx_type ydc_tf;
                cplx_type ydc_tt;

                TwoSidesContainer_rxh_AInfo(const TwoSidesContainer_rxh_A & r_data, int my_id):
                TwoSidesContainer<OneSideContainer>::TwoSidesInfo(r_data, my_id),
                r_pu(-1.0),
                x_pu(-1.0),
                h_pu(0., 0.),
                has_res(false),
                res_a1_ka(0.),
                res_a2_ka(0.),
                yac_ff(0., 0.),
                yac_ft(0., 0.),
                yac_tf(0., 0.),
                yac_tt(0., 0.),
                ydc_ff(0., 0.),
                ydc_ft(0., 0.),
                ydc_tf(0., 0.),
                ydc_tt(0., 0.)
                {
                    if(my_id < 0) return;
                    if(my_id >= r_data.nb()) return;
                    r_pu = r_data.r_.coeff(my_id);
                    x_pu = r_data.x_.coeff(my_id);
                    h_pu = r_data.h_.coeff(my_id);

                    has_res = r_data.side_1_[my_id].has_res;
                    if(has_res)
                    {
                        res_a1_ka = r_data.res_a_side_1_.coeff(my_id);
                        res_a2_ka = r_data.res_a_side_2_.coeff(my_id);
                    }

                    // coeffs
                    yac_ff = r_data.yac_ff_.coeff(my_id);
                    yac_ft = r_data.yac_ft_.coeff(my_id);
                    yac_tf = r_data.yac_tf_.coeff(my_id);
                    yac_tt = r_data.yac_tt_.coeff(my_id);
                    ydc_ff = r_data.ydc_ff_.coeff(my_id);
                    ydc_ft = r_data.ydc_ft_.coeff(my_id);
                    ydc_tf = r_data.ydc_tf_.coeff(my_id);
                    ydc_tt = r_data.ydc_tt_.coeff(my_id);

                }
        };
        typedef TwoSidesContainer_rxh_AInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<TwoSidesContainer_rxh_A> TwoSidesContainer_rxh_AInfoConstIterator;

    // make it iterable
    public:
        TwoSidesContainer_rxh_AInfoConstIterator begin() const {return TwoSidesContainer_rxh_AInfoConstIterator(this, 0); }
        TwoSidesContainer_rxh_AInfoConstIterator end() const {return TwoSidesContainer_rxh_AInfoConstIterator(this, nb()); }
        TwoSidesContainer_rxh_AInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("id should be >= 0");
            }
            if(id >= nb())
            {
                throw std::range_error("Out of bound error.");
            }
            return TwoSidesContainer_rxh_AInfo(*this, id);
        }
        /////////////////////////

    public:
        // pickle
        typedef std::tuple<
                   TwoSidesContainer<OneSideContainer>::StateRes,
                   std::vector<real_type>,  // branch_r
                   std::vector<real_type>,  // branch_x
                   std::vector<cplx_type>   // branch_h
               >  StateRes;

        // setter (states)
        // methods used within lightsim
        void deactivate(int el_id, SolverControl & solver_control) {
            if(status_global_[el_id]){
                solver_control.tell_recompute_ybus();
                // but sparsity pattern do not change here (possibly one more coeff at 0.)
                solver_control.tell_ybus_some_coeffs_zero();
                solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
            }
            _generic_deactivate(el_id, status_global_);
            side_1_.deactivate(el_id, solver_control);
            side_2_.deactivate(el_id, solver_control);
        }
        void reactivate(int el_id, SolverControl & solver_control) {
            if(!status_global_[el_id]){
                solver_control.tell_recompute_ybus();
                solver_control.tell_ybus_change_sparsity_pattern();  // this might change
                solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
            }
            _generic_reactivate(el_id, status_global_);
            side_1_.reactivate(el_id, solver_control);
            side_2_.reactivate(el_id, solver_control);
        }
        void change_bus_side_1(int el_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
            // _generic_change_bus(trafo_id, new_bus_id, get_buses_not_const_side_1(), solver_control, nb_bus);
            if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 1).");
            side_1_.change_bus(el_id, new_bus_id, solver_control, nb_bus);
        }
        void change_bus_side_2(int el_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
            // _generic_change_bus(trafo_id, new_bus_id, get_buses_not_const_side_1(), solver_control, nb_bus);
            if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 2).");
            side_2_.change_bus(el_id, new_bus_id, solver_control, nb_bus);
        }

        // getter (results)
        tuple4d get_res_side_1() const {
            const tuple3d & side_1_res = TwoSidesContainer<OneSideContainer>::get_res_side_1();
            return tuple4d(
                std::get<0>(side_1_res),
                std::get<1>(side_1_res),
                std::get<2>(side_1_res),
                res_a_side_1_);
        }
        tuple4d get_res_side_2() const {
            const tuple3d & side_2_res = TwoSidesContainer<OneSideContainer>::get_res_side_2();
            return tuple4d(
                std::get<0>(side_2_res),
                std::get<1>(side_2_res),
                std::get<2>(side_2_res),
                res_a_side_2_);
        }
        tuple5d get_res_side_1_full() const {
            const tuple4d & side_1_res = TwoSidesContainer<OneSideContainer>::get_res_full_side_1();
            return tuple5d(
                std::get<0>(side_1_res),
                std::get<1>(side_1_res),
                std::get<2>(side_1_res),
                res_a_side_1_,
                std::get<3>(side_1_res));
        }
        tuple5d get_res_side_2_full() const {
            const tuple4d & side_2_res = TwoSidesContainer<OneSideContainer>::get_res_full_side_2();
            return tuple5d(
                std::get<0>(side_2_res),
                std::get<1>(side_2_res),
                std::get<2>(side_2_res),
                res_a_side_2_,
                std::get<3>(side_2_res));
        }

        void compute_results_tsc_rxha(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<int> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac
                                      )
        {
            // it needs to be initialized at 0.
            const int nb_element = nb();
            Eigen::Ref<RealVect> res_p_side_1 = get_res_p_side_1();
            Eigen::Ref<RealVect> res_q_side_1 = get_res_q_side_1();
            Eigen::Ref<RealVect> res_v_side_1 = get_res_v_side_1();
            Eigen::Ref<RealVect> res_theta_side_1 = get_res_theta_side_1();

            Eigen::Ref<RealVect> res_p_lv_ = get_res_p_side_2();
            Eigen::Ref<RealVect> res_q_lv_ = get_res_q_side_2();
            Eigen::Ref<RealVect> res_v_lv_ = get_res_v_side_2();
            Eigen::Ref<RealVect> res_theta_lv_ = get_res_theta_side_2();

            for(int el_id = 0; el_id < nb_element; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) {
                    res_p_side_1(el_id) = 0.0;  // in MW
                    res_q_side_1(el_id) = 0.0;  // in MVar
                    res_v_side_1(el_id) = v_disco_el_;  // in kV
                    res_a_side_1_(el_id) = 0.0;  // in kA
                    res_p_lv_(el_id) = 0.0;  // in MW
                    res_q_lv_(el_id) = 0.0;  // in MVar
                    res_v_lv_(el_id) = v_disco_el_;  // in kV
                    res_a_side_2_(el_id) = 0.0;  // in kA
                    res_theta_side_1(el_id) = theta_disco_el_;  // in degree
                    res_theta_lv_(el_id) = theta_disco_el_;  // in degree            
                    continue;
                }

                // connectivity
                int bus_hv_id_me = get_bus_side_1(el_id);
                int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
                if(bus_hv_solver_id == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TrafoContainer::compute_results: the trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (hv side) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                int bus_lv_id_me = get_bus_side_2(el_id);
                int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
                if(bus_lv_solver_id == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TrafoContainer::compute_results: the trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (lv side) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }

                // retrieve voltages magnitude in kv instead of pu
                real_type v_hv = Vm(bus_hv_solver_id);
                real_type v_lv = Vm(bus_lv_solver_id);
                real_type bus_vn_kv_hv = bus_vn_kv(bus_hv_id_me);
                real_type bus_vn_kv_lv = bus_vn_kv(bus_lv_id_me);           

                // for voltages
                res_v_side_1(el_id) = v_hv * bus_vn_kv_hv;
                res_v_lv_(el_id) = v_lv * bus_vn_kv_lv;

                res_theta_side_1(el_id) = Va(bus_hv_solver_id) * my_180_pi_;
                res_theta_lv_(el_id) = Va(bus_lv_solver_id) * my_180_pi_;

                if(ac){
                    // results of the ac powerflow
                    cplx_type Ehv = V(bus_hv_solver_id);
                    cplx_type Elv = V(bus_lv_solver_id);

                    // TODO for DC with yff, ...
                    // trafo equations
                    cplx_type I_hvlv =  yac_ff_(el_id) * Ehv + yac_ft_(el_id) * Elv;
                    cplx_type I_lvhv =  yac_tt_(el_id) * Elv + yac_tf_(el_id) * Ehv;

                    I_hvlv = std::conj(I_hvlv);
                    I_lvhv = std::conj(I_lvhv);
                    cplx_type s_hvlv = Ehv * I_hvlv;
                    cplx_type s_lvhv = Elv * I_lvhv;

                    res_p_side_1(el_id) = std::real(s_hvlv) * sn_mva;
                    res_q_side_1(el_id) = std::imag(s_hvlv) * sn_mva;
                    res_p_lv_(el_id) = std::real(s_lvhv) * sn_mva;
                    res_q_lv_(el_id) = std::imag(s_lvhv) * sn_mva;
                }else{
                    // result of the dc powerflow
                    res_p_side_1(el_id) = (std::real(ydc_ff_(el_id)) * Va(bus_hv_solver_id) + std::real(ydc_ft_(el_id)) * Va(bus_lv_solver_id)) * sn_mva; // - dc_x_tau_shift_(el_id) ) * sn_mva;
                    res_p_lv_(el_id) = (std::real(ydc_tt_(el_id)) * Va(bus_lv_solver_id) + std::real(ydc_tf_(el_id)) * Va(bus_hv_solver_id)) * sn_mva; // + dc_x_tau_shift_(el_id) ) * sn_mva; 
                
                    res_q_side_1(el_id) = 0.;
                    res_q_lv_(el_id) = 0.;
                    
                    // for voltages, because vm = 1. pu by hypothesis
                    // res_v_hv_(trafo_id) = bus_vn_kv_hv;
                    // res_v_lv_(trafo_id) = bus_vn_kv_lv;
                }

            }
            const auto & res_side1 = side_1_.get_res();
            _get_amps(res_a_side_1_, std::get<0>(res_side1), std::get<1>(res_side1), std::get<2>(res_side1));
            const auto & res_side2 = side_2_.get_res();
            _get_amps(res_a_side_2_, std::get<0>(res_side2), std::get<1>(res_side2), std::get<2>(res_side2));
        }
        // virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component)
        // {
        //     const Eigen::Index nb_el = nb();
        //     SolverControl unused_solver_control;
        //     for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
        //         if(!status_global_[i]) continue;
        //         auto bus_or = get_bus_side_1(el_id);
        //         auto bus_ex = get_bus_side_2(el_id);
        //         if(!busbar_in_main_component[bus_or] || !busbar_in_main_component[bus_ex]){
        //             deactivate(el_id, unused_solver_control);
        //         }
        //     }            
        // }
        
        virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const
        {
            const auto my_size = nb();
            for(Eigen::Index el_id = 0; el_id < my_size; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;
                const auto bus_or = get_bus_side_1(el_id);
                const auto bus_ex = get_bus_side_2(el_id);
                res.push_back(Eigen::Triplet<real_type>(bus_or, bus_ex, 1.));
                res.push_back(Eigen::Triplet<real_type>(bus_ex, bus_or, 1.));
            }
        }

        // model paramters
        Eigen::Ref<const CplxVect> yac_ff() const {return yac_ff_;}
        Eigen::Ref<const CplxVect> yac_ft() const {return yac_ft_;}
        Eigen::Ref<const CplxVect> yac_tf() const {return yac_tf_;}
        Eigen::Ref<const CplxVect> yac_tt() const {return yac_tt_;}

        Eigen::Ref<const CplxVect> ydc_ff() const {return ydc_ff_;}
        Eigen::Ref<const CplxVect> ydc_ft() const {return ydc_ft_;}
        Eigen::Ref<const CplxVect> ydc_tf() const {return ydc_tf_;}
        Eigen::Ref<const CplxVect> ydc_tt() const {return ydc_tt_;}

        // solver interface
        virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                        bool ac,
                        const std::vector<int> & id_grid_to_solver,
                        real_type sn_mva) const
        {
            const Eigen::Index nb_els = nb();
            cplx_type yft, ytf, yff, ytt;
            for(Eigen::Index el_id =0; el_id < nb_els; ++el_id){
                // i don't do anything if the trafo is disconnected
                if(!status_global_[el_id]) continue;

                // compute from / to
                int bus_side1_id_me = get_bus_side_1(el_id);
                int bus_side1_solver_id = id_grid_to_solver[bus_side1_id_me];
                if(bus_side1_solver_id == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillYbus: the trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (side 1) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                int bus_side2_id_me = get_bus_side_2(el_id);
                int bus_side2_solver_id = id_grid_to_solver[bus_side2_id_me];
                if(bus_side2_solver_id == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillYbus: the trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (side 2) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                
                if(ac){
                    // ac mode
                    yft = yac_ft_(el_id);
                    ytf = yac_tf_(el_id);
                    yff = yac_ff_(el_id);
                    ytt = yac_tt_(el_id);
                }else{
                    // dc mode
                    yft = ydc_ft_(el_id);
                    ytf = ydc_tf_(el_id);
                    yff = ydc_ff_(el_id);
                    ytt = ydc_tt_(el_id);
                }
                res.push_back(Eigen::Triplet<cplx_type> (bus_side1_solver_id, bus_side2_solver_id, yft));
                res.push_back(Eigen::Triplet<cplx_type> (bus_side2_solver_id, bus_side1_solver_id, ytf));
                res.push_back(Eigen::Triplet<cplx_type> (bus_side1_solver_id, bus_side1_solver_id, yff));
                res.push_back(Eigen::Triplet<cplx_type> (bus_side2_solver_id, bus_side2_solver_id, ytt));
            }
        }

        // gridmodel utilities
        void reconnect_connected_buses(Substation & substation) const{
            const Eigen::Index nb_els = nb();
            for(Eigen::Index el_id = 0; el_id < nb_els; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;
                
                const auto bus_or_id_me = get_bus_side_1(el_id);        
                if(bus_or_id_me == _deactivated_bus_id){
                    // TODO DEBUG MODE only this in debug mode
                    std::ostringstream exc_;
                    exc_ << "TrafoContainer::reconnect_connected_buses: Trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (hv) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_trafo(...)` ?.";
                    throw std::runtime_error(exc_.str());
                }
                // bus_status[bus_or_id_me] = true;
                substation.reconnect_bus(bus_or_id_me);

                const auto bus_ex_id_me = get_bus_side_2(el_id);        
                if(bus_ex_id_me == _deactivated_bus_id){
                    // TODO DEBUG MODE only this in debug mode
                    std::ostringstream exc_;
                    exc_ << "TrafoContainer::reconnect_connected_buses: Trafo with id ";
                    exc_ << el_id;
                    exc_ << " is connected (lv) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_trafo(...)` ?.";
                    throw std::runtime_error(exc_.str());
                }
                // bus_status[bus_ex_id_me] = true;
                substation.reconnect_bus(bus_ex_id_me);
            }
        }
    protected:

        StateRes get_tsc_rxha_state() const  // tsc: two sides container
        {
            std::vector<real_type> branch_r(r_.begin(), r_.end());
            std::vector<real_type> branch_x(x_.begin(), x_.end());
            std::vector<cplx_type> branch_h(h_.begin(), h_.end());
            StateRes res(
                get_tsc_state(),
                branch_r,
                branch_x,
                branch_h
            );
            return res;
        }

        void set_tsc_rxha_state(StateRes & my_state)  // tsc: two sides container
        {
            set_tsc_state(std::get<0>(my_state));
            const auto size = nb();

            const std::vector<real_type> & branch_r = std::get<1>(my_state);
            const std::vector<real_type> & branch_x = std::get<2>(my_state);
            const std::vector<cplx_type> & branch_h = std::get<3>(my_state);
            check_size(branch_r, size, "branch r");
            check_size(branch_x, size, "branch x");
            check_size(branch_h, size, "branch h (=g+j.b)");

            r_ = RealVect::Map(&branch_r[0], size);
            x_ = RealVect::Map(&branch_x[0], size);
            h_ = CplxVect::Map(&branch_h[0], size);
        }

        void reset_results_tsc_rxha(){
            reset_results_tsc();
            res_a_side_1_ = RealVect(nb());  // in kA
            res_a_side_2_ = RealVect(nb());  // in kA
        }


    protected:
        // physical properties
        RealVect r_;  // in pu
        RealVect x_;  // in pu
        CplxVect h_;  // in pu

        //output data
        RealVect res_a_side_1_;  // in kA
        RealVect res_a_side_2_;  // in kA

        // model coefficients
        CplxVect yac_ff_;
        CplxVect yac_ft_;
        CplxVect yac_tf_;
        CplxVect yac_tt_;

        CplxVect ydc_ff_;
        CplxVect ydc_ft_;
        CplxVect ydc_tf_;
        CplxVect ydc_tt_;
};
#endif  // TWO_SIDES_CONTAINER_RXH_A_H
