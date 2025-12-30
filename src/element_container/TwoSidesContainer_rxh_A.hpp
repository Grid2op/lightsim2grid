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
        using TwoSidesContainer<OneSideType>::_deactivated_bus_id;
        using TwoSidesContainer<OneSideType>::v_disco_el_;
        using TwoSidesContainer<OneSideType>::theta_disco_el_;
        using TwoSidesContainer<OneSideType>::my_180_pi_;
        using TwoSidesContainer<OneSideType>::_generic_deactivate;
        using TwoSidesContainer<OneSideType>::_generic_reactivate;

        using TwoSidesContainer<OneSideType>::nb;
        using TwoSidesContainer<OneSideType>::reset_results_tsc;
        using TwoSidesContainer<OneSideType>::check_size;
        using TwoSidesContainer<OneSideType>::get_bus_side_1;
        using TwoSidesContainer<OneSideType>::get_bus_side_2;

    protected:
        using TwoSidesContainer<OneSideType>::_get_amps;

        using TwoSidesContainer<OneSideType>::get_tsc_state;
        using TwoSidesContainer<OneSideType>::set_tsc_state;
        using TwoSidesContainer<OneSideType>::status_global_;
        using TwoSidesContainer<OneSideType>::get_status_side_1;
        using TwoSidesContainer<OneSideType>::get_status_side_2;
        using TwoSidesContainer<OneSideType>::side_1_;
        using TwoSidesContainer<OneSideType>::side_2_;
        using TwoSidesContainer<OneSideType>::get_res_p_side_1;
        using TwoSidesContainer<OneSideType>::get_res_p_side_2;
        using TwoSidesContainer<OneSideType>::get_res_q_side_1;
        using TwoSidesContainer<OneSideType>::get_res_q_side_2;
        using TwoSidesContainer<OneSideType>::get_res_v_side_1;
        using TwoSidesContainer<OneSideType>::get_res_v_side_2;
        using TwoSidesContainer<OneSideType>::get_res_theta_side_1;
        using TwoSidesContainer<OneSideType>::get_res_theta_side_2;

    typedef typename TwoSidesContainer<OneSideType>::StateRes StateResSuper;
    //////////////////////////////

    public:
        class TwoSidesContainer_rxh_AInfo : public TwoSidesContainer<OneSideType>::TwoSidesInfo
        {
            public:
                // members
                real_type r_pu;
                real_type x_pu;
                cplx_type h1_pu;
                cplx_type h2_pu;

                bool has_res;
                real_type res_a1_ka;
                real_type res_a2_ka;

                cplx_type yac_11;
                cplx_type yac_12;
                cplx_type yac_21;
                cplx_type yac_22;
                cplx_type ydc_11;
                cplx_type ydc_12;
                cplx_type ydc_21;
                cplx_type ydc_22;

                TwoSidesContainer_rxh_AInfo(const TwoSidesContainer_rxh_A & r_data, int my_id) noexcept:
                TwoSidesContainer<OneSideType>::TwoSidesInfo(r_data, my_id),
                r_pu(-1.0),
                x_pu(-1.0),
                h1_pu(0., 0.),
                h2_pu(0., 0.),
                has_res(false),
                res_a1_ka(0.),
                res_a2_ka(0.),
                yac_11(0., 0.),
                yac_12(0., 0.),
                yac_21(0., 0.),
                yac_22(0., 0.),
                ydc_11(0., 0.),
                ydc_12(0., 0.),
                ydc_21(0., 0.),
                ydc_22(0., 0.)
                {
                    if(my_id < 0) return;
                    if(my_id >= r_data.nb()) return;
                    r_pu = r_data.r_.coeff(my_id);
                    x_pu = r_data.x_.coeff(my_id);
                    h1_pu = r_data.h_side_1_.coeff(my_id);
                    h2_pu = r_data.h_side_2_.coeff(my_id);

                    has_res = r_data.side_1_[my_id].has_res;
                    if(has_res)
                    {
                        res_a1_ka = r_data.res_a_side_1_.coeff(my_id);
                        res_a2_ka = r_data.res_a_side_2_.coeff(my_id);
                    }

                    // coeffs
                    yac_11 = r_data.yac_11_.coeff(my_id);
                    yac_12 = r_data.yac_12_.coeff(my_id);
                    yac_21 = r_data.yac_21_.coeff(my_id);
                    yac_22 = r_data.yac_22_.coeff(my_id);
                    ydc_11 = r_data.ydc_11_.coeff(my_id);
                    ydc_12 = r_data.ydc_12_.coeff(my_id);
                    ydc_21 = r_data.ydc_21_.coeff(my_id);
                    ydc_22 = r_data.ydc_22_.coeff(my_id);

                }
        };
        typedef TwoSidesContainer_rxh_AInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<TwoSidesContainer_rxh_A> TwoSidesContainer_rxh_AInfoConstIterator;

    public:
        TwoSidesContainer_rxh_A() noexcept = default;
        virtual ~TwoSidesContainer_rxh_A() noexcept{
            // std::cout << "\tTwoSidesContainer_rxh_A destructor" << std::endl;
        }
        
        // pickle
        typedef std::tuple<
                   StateResSuper,
                   std::vector<real_type>,  // branch_r
                   std::vector<real_type>,  // branch_x
                   std::vector<cplx_type>,   // branch_h1
                   std::vector<cplx_type>   // branch_h2
               >  StateRes;

        // getter (results)
        tuple4d get_res_side_1() const {
            const tuple3d & side_1_res = TwoSidesContainer<OneSideType>::get_res_side_1();
            return tuple4d(
                std::get<0>(side_1_res),
                std::get<1>(side_1_res),
                std::get<2>(side_1_res),
                res_a_side_1_);
        }
        tuple4d get_res_side_2() const {
            const tuple3d & side_2_res = TwoSidesContainer<OneSideType>::get_res_side_2();
            return tuple4d(
                std::get<0>(side_2_res),
                std::get<1>(side_2_res),
                std::get<2>(side_2_res),
                res_a_side_2_);
        }
        tuple5d get_res_full_side_1() const {
            const tuple4d & side_1_res = TwoSidesContainer<OneSideType>::get_res_full_side_1();
            return tuple5d(
                std::get<0>(side_1_res),
                std::get<1>(side_1_res),
                std::get<2>(side_1_res),
                res_a_side_1_,
                std::get<3>(side_1_res));
        }
        tuple5d get_res_full_side_2() const {
            const tuple4d & side_2_res = TwoSidesContainer<OneSideType>::get_res_full_side_2();
            return tuple5d(
                std::get<0>(side_2_res),
                std::get<1>(side_2_res),
                std::get<2>(side_2_res),
                res_a_side_2_,
                std::get<3>(side_2_res));
        }

        void compute_results_tsc_rxha_no_amps(
            const Eigen::Ref<const RealVect> & Va,
            const Eigen::Ref<const RealVect> & Vm,
            const Eigen::Ref<const CplxVect> & V,
            const std::vector<SolverBusId> & id_grid_to_solver,
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

            Eigen::Ref<RealVect> res_p_side_2 = get_res_p_side_2();
            Eigen::Ref<RealVect> res_q_side_2 = get_res_q_side_2();
            Eigen::Ref<RealVect> res_v_side_2 = get_res_v_side_2();
            Eigen::Ref<RealVect> res_theta_side_2 = get_res_theta_side_2();

            const std::vector<bool> & status1 = side_1_.get_status();
            const std::vector<bool> & status2 = side_2_.get_status();
            for(int el_id = 0; el_id < nb_element; ++el_id){

                // don't do anything if the element is disconnected
                if(!status_global_[el_id] || (!status1[el_id] && !status2[el_id])) {
                    res_p_side_1(el_id) = 0.0;  // in MW
                    res_q_side_1(el_id) = 0.0;  // in MVar
                    res_v_side_1(el_id) = v_disco_el_;  // in kV
                    res_a_side_1_(el_id) = 0.0;  // in kA
                    res_p_side_2(el_id) = 0.0;  // in MW
                    res_q_side_2(el_id) = 0.0;  // in MVar
                    res_v_side_2(el_id) = v_disco_el_;  // in kV
                    res_a_side_2_(el_id) = 0.0;  // in kA
                    res_theta_side_1(el_id) = theta_disco_el_;  // in degree
                    res_theta_side_2(el_id) = theta_disco_el_;  // in degree            
                    continue;
                }

                // connectivity
                GridModelBusId bus_hv_id_me, bus_lv_id_me;
                SolverBusId bus_hv_solver_id, bus_lv_solver_id;
                if(status1[el_id]){
                    bus_hv_id_me = get_bus_side_1(el_id);
                    if(bus_hv_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::compute_results: (GlobalBusId) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me.cast_int()];
                    if(bus_hv_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::compute_results: (SolverBusId) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }else{
                    bus_hv_id_me = GridModelBusId(_deactivated_bus_id);
                    bus_hv_solver_id = SolverBusId(_deactivated_bus_id);
                }

                if(status2[el_id]){
                    bus_lv_id_me = get_bus_side_2(el_id);
                    if(bus_lv_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::compute_results: (GlobalBusId) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me.cast_int()];
                    if(bus_lv_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::compute_results: (SolverBusId) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }else{
                    bus_lv_id_me = GridModelBusId(_deactivated_bus_id);
                    bus_lv_solver_id = SolverBusId(_deactivated_bus_id);
                }

                // retrieve voltages magnitude in kv instead of pu
                if(status1[el_id]){
                    real_type v_hv = Vm(bus_hv_solver_id.cast_int());
                    real_type bus_vn_kv_hv = bus_vn_kv(bus_hv_id_me.cast_int());
                    res_v_side_1(el_id) = v_hv * bus_vn_kv_hv;
                    res_theta_side_1(el_id) = Va(bus_hv_solver_id.cast_int()) * my_180_pi_;
                }else{
                    res_v_side_1(el_id) = v_disco_el_;
                    res_theta_side_1(el_id) = theta_disco_el_;
                }

                if(status2[el_id]){
                    real_type v_lv = Vm(bus_lv_solver_id.cast_int());
                    real_type bus_vn_kv_lv = bus_vn_kv(bus_lv_id_me.cast_int());  
                    res_v_side_2(el_id) = v_lv * bus_vn_kv_lv;
                    res_theta_side_2(el_id) = Va(bus_lv_solver_id.cast_int()) * my_180_pi_;
                }else{
                    res_v_side_2(el_id) = v_disco_el_;
                    res_theta_side_2(el_id) = theta_disco_el_;
                }

                if(ac){
                    // results of the ac powerflow
                    cplx_type Ehv = status1[el_id] ? V(bus_hv_solver_id.cast_int()) : - yac_12_(el_id) *  V(bus_lv_solver_id.cast_int()) / yac_11_(el_id);
                    cplx_type Elv = status2[el_id] ? V(bus_lv_solver_id.cast_int()) : - yac_21_(el_id) *  V(bus_hv_solver_id.cast_int()) / yac_22_(el_id);

                    // TODO for DC with yff, ...
                    // trafo equations
                    cplx_type I_hvlv =  (yac_11_(el_id) * Ehv + yac_12_(el_id) * Elv);
                    cplx_type I_lvhv =  (yac_22_(el_id) * Elv + yac_21_(el_id) * Ehv);

                    I_hvlv = std::conj(I_hvlv);
                    I_lvhv = std::conj(I_lvhv);
                    cplx_type s_hvlv = Ehv * I_hvlv;
                    cplx_type s_lvhv = Elv * I_lvhv;

                    res_p_side_1(el_id) = std::real(s_hvlv) * sn_mva;
                    res_q_side_1(el_id) = std::imag(s_hvlv) * sn_mva;
                    res_p_side_2(el_id) = std::real(s_lvhv) * sn_mva;
                    res_q_side_2(el_id) = std::imag(s_lvhv) * sn_mva;
                }else{
                    // result of the dc powerflow
                    if(status1[el_id] && status2[el_id]){
                        real_type Va_hv = Va(bus_hv_solver_id.cast_int());
                        real_type Va_lv = Va(bus_lv_solver_id.cast_int());
                        res_p_side_1(el_id) = (std::real(ydc_11_(el_id)) * Va_hv + std::real(ydc_12_(el_id)) * Va_lv) * sn_mva;
                        res_p_side_2(el_id) = (std::real(ydc_22_(el_id)) * Va_lv + std::real(ydc_21_(el_id)) * Va_hv) * sn_mva;
                    }else{
                        res_p_side_1(el_id) = 0.;
                        res_p_side_2(el_id) = 0.;
                    }
                    res_q_side_1(el_id) = 0.;
                    res_q_side_2(el_id) = 0.;
                    
                    // for voltages, because vm = 1. pu by hypothesis
                    // res_v_hv_(trafo_id) = bus_vn_kv_hv;
                    // res_v_lv_(trafo_id) = bus_vn_kv_lv;
                }

            }
        }
        
        // computes the flows (in A) after p_1, q_1, v_1, p_2, q_2 and v_2 have been computed
        // for example with compute_results_tsc_rxha_no_amps
        void compute_amps_after_all_set(){
            const auto & res_side1 = side_1_.get_res();
            _get_amps(res_a_side_1_, std::get<0>(res_side1), std::get<1>(res_side1), std::get<2>(res_side1));
            const auto & res_side2 = side_2_.get_res();
            _get_amps(res_a_side_2_, std::get<0>(res_side2), std::get<1>(res_side2), std::get<2>(res_side2));
        }
        void compute_results_tsc_rxha(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<SolverBusId> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac
                                      )
        {
            // it needs to be initialized at 0.
            compute_results_tsc_rxha_no_amps(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
            compute_amps_after_all_set();
        }
        
        virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const
        {
            const auto my_size = nb();
            for(Eigen::Index el_id = 0; el_id < my_size; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;
                const GridModelBusId bus_or = get_bus_side_1(el_id);
                const GridModelBusId bus_ex = get_bus_side_2(el_id);
                if((bus_or.cast_int() != _deactivated_bus_id) && 
                   (bus_ex.cast_int() != _deactivated_bus_id)){
                    res.push_back(Eigen::Triplet<real_type>(bus_or.cast_int(), bus_ex.cast_int(), 1.));
                    res.push_back(Eigen::Triplet<real_type>(bus_ex.cast_int(), bus_or.cast_int(), 1.));
                }
            }
        }

        // model paramters
        Eigen::Ref<const CplxVect> yac_11() const {return yac_11_;}
        Eigen::Ref<const CplxVect> yac_12() const {return yac_12_;}
        Eigen::Ref<const CplxVect> yac_21() const {return yac_21_;}
        Eigen::Ref<const CplxVect> yac_22() const {return yac_22_;}

        Eigen::Ref<const CplxVect> ydc_11() const {return ydc_11_;}
        Eigen::Ref<const CplxVect> ydc_12() const {return ydc_12_;}
        Eigen::Ref<const CplxVect> ydc_21() const {return ydc_21_;}
        Eigen::Ref<const CplxVect> ydc_22() const {return ydc_22_;}

        // solver interface
        virtual void fillYbus(
            std::vector<Eigen::Triplet<cplx_type> > & res,
            bool ac,
            const std::vector<SolverBusId> & id_grid_to_solver,
            real_type sn_mva) const
        {
            const Eigen::Index nb_els = nb();
            const std::vector<bool> & status1 = side_1_.get_status();
            const std::vector<bool> & status2 = side_2_.get_status();

            cplx_type yft, ytf, yff, ytt;
            for(Eigen::Index el_id =0; el_id < nb_els; ++el_id){
                // i don't do anything if the trafo is disconnected
                if(!status_global_[el_id]  || (!status1[el_id] && !status2[el_id])) continue;

                // compute from / to
                GlobalBusId bus_side1_id_me, bus_side2_id_me;
                SolverBusId bus_side1_solver_id, bus_side2_solver_id;
                bool status1_me = status1[el_id];
                bool status2_me = status2[el_id];
                if(status1_me){
                    bus_side1_id_me = get_bus_side_1(el_id);
                    if(bus_side1_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillYbus: (GlobalID) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_side1_solver_id = id_grid_to_solver[bus_side1_id_me.cast_int()];
                    if(bus_side1_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillYbus: (SolverID) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }

                if(status2_me){
                    bus_side2_id_me = get_bus_side_2(el_id);
                    if(bus_side2_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillYbus: (GlobalID) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_side2_solver_id = id_grid_to_solver[bus_side2_id_me.cast_int()];
                    if(bus_side2_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillYbus: (SolverID) the branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }
                
                if(ac){
                    // ac mode
                    yft = yac_12_(el_id);
                    ytf = yac_21_(el_id);
                    yff = yac_11_(el_id);
                    ytt = yac_22_(el_id);
                    if(!status1_me){
                        // I know that the powerline is connected on side 2 
                        // otherwise I would not be here
                        // (nothing is done if neither side 1 nor side 2 are connected)
                        ytt -= ytf * yft / yff;
                    }
                    if(!status2_me){
                        // I know that the powerline is connected on side 1
                        // otherwise I would not be here
                        // (nothing is done if neither side 1 nor side 2 are connected)
                        yff -= ytf * yft / ytt;
                    }
                }else{
                    // dc mode
                    yft = ydc_12_(el_id);
                    ytf = ydc_21_(el_id);
                    yff = ydc_11_(el_id);
                    ytt = ydc_22_(el_id);
                    // In DC disconnected on one side == disco on both sides
                    if((!status1_me) || (!status2_me)){
                        status1_me = false;
                        status2_me = false;
                    }
                }
                if(status1_me) res.push_back(Eigen::Triplet<cplx_type> (bus_side1_solver_id.cast_int(), bus_side1_solver_id.cast_int(), yff));
                if(status2_me) res.push_back(Eigen::Triplet<cplx_type> (bus_side2_solver_id.cast_int(), bus_side2_solver_id.cast_int(), ytt));
                if(status1_me && status2_me){
                    res.push_back(Eigen::Triplet<cplx_type> (bus_side1_solver_id.cast_int(), bus_side2_solver_id.cast_int(), yft));
                    res.push_back(Eigen::Triplet<cplx_type> (bus_side2_solver_id.cast_int(), bus_side1_solver_id.cast_int(), ytf));
                }
            }
        }

        virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<SolverBusId> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const
        {

            // For Bp
            // temp_branch[:, BR_B] = zeros(nl)           ## zero out line charging shunts
            // temp_branch[:, TAP] = ones(nl)             ## cancel out taps
            // if alg == 2:                               ## if XB method
            //    temp_branch[:, BR_R] = zeros(nl)       ## zero out line resistance

            // For Bpp
            // temp_branch[:, SHIFT] = zeros(nl)          ## zero out phase shifters
            // if alg == 3:                               ## if BX method
            //     temp_branch[:, BR_R] = zeros(nl)    ## zero out line resistance
            const Eigen::Index nb_trafo = nb();

            for(Eigen::Index el_id=0; el_id < nb_trafo; ++el_id){
                // i only add this if the powerline is connected
                if(!status_global_[el_id]) continue;

                GlobalBusId bus_or_id_me, bus_ex_id_me;
                SolverBusId bus_or_solver_id , bus_ex_solver_id;
                // get the from / to bus id
                if(side_1_.get_status(el_id))
                {
                    bus_or_id_me = get_bus_side_1(el_id);
                    if(bus_or_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillBp_Bpp: (GlobalId) the branch (line or trafo) with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_or_solver_id = id_grid_to_solver[bus_or_id_me.cast_int()];
                    if(bus_or_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillBp_Bpp: (SolverId) the branch (line or trafo) with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }else{
                    throw std::runtime_error("FDPF algorithm does not handle lines / trafos disconnected at only one side at the moment.");
                }
                if(side_2_.get_status(el_id))
                {
                    bus_ex_id_me = get_bus_side_2(el_id);
                    if(bus_ex_id_me.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillBp_Bpp: (GlobalId) the branch (line or trafo) with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                    bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me.cast_int()];
                    if(bus_ex_solver_id.cast_int() == _deactivated_bus_id){
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::fillBp_Bpp: (SolverId) the branch (line or trafo) with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to a disconnected bus while being connected";
                        throw std::runtime_error(exc_.str());
                    }
                }else{
                    throw std::runtime_error("FDPF algorithm does not handle lines / trafos disconnected at only one side at the moment.");
                }

                const FDPFCoeffs & coeffs = this->get_fdpf_coeffs(el_id, xb_or_bx);

                // and now add them
                if(side_1_.get_status(el_id)){
                    Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id.cast_int(), bus_or_solver_id.cast_int(), -coeffs.yff_bp));
                    Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id.cast_int(), bus_or_solver_id.cast_int(), -coeffs.yff_bpp));
                }
                if(side_2_.get_status(el_id)){
                    Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id.cast_int(), bus_ex_solver_id.cast_int(), -coeffs.ytt_bp));
                    Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id.cast_int(), bus_ex_solver_id.cast_int(), -coeffs.ytt_bpp));
                }
                if(side_1_.get_status(el_id) && side_2_.get_status(el_id)){
                    Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id.cast_int(), bus_ex_solver_id.cast_int(), -coeffs.yft_bp));
                    Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id.cast_int(), bus_or_solver_id.cast_int(), -coeffs.ytf_bp));
                    Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id.cast_int(), bus_ex_solver_id.cast_int(), -coeffs.yft_bpp));
                    Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id.cast_int(), bus_or_solver_id.cast_int(), -coeffs.ytf_bpp));
                }
            }
        }

        void fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & Bf,
                             const std::vector<SolverBusId> & id_grid_to_solver,
                             real_type sn_mva,
                             int nb_powerline,
                             bool transpose) const
        {
            const Eigen::Index nb_line = nb();
            const std::vector<bool> & side1_conn = side_1_.get_status();
            const std::vector<bool> & side2_conn = side_2_.get_status();
            for(Eigen::Index line_id=0; line_id < nb_line; ++line_id){
                // i only add this if the powerline is connected
                if(!status_global_[line_id]) continue;
                if(!side1_conn[line_id]) continue;
                if(!side2_conn[line_id]) continue;
                
                // get the from / to bus id
                GlobalBusId bus_or_id_me = get_bus_side_1(line_id);
                if(bus_or_id_me.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillBf_for_PTDF: (GlobalId) the line/trafo with id ";
                    exc_ << line_id;
                    exc_ << " is connected (side 1) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                SolverBusId bus_or_solver_id = id_grid_to_solver[bus_or_id_me.cast_int()];
                if(bus_or_solver_id.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillBf_for_PTDF: (SolverId) the line/trafo with id ";
                    exc_ << line_id;
                    exc_ << " is connected (side 1) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                GlobalBusId bus_ex_id_me = get_bus_side_2(line_id);
                if(bus_ex_id_me.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillBf_for_PTDF: (GlobalId) the line/trafo with id ";
                    exc_ << line_id;
                    exc_ << " is connected (side 2) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                SolverBusId bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me.cast_int()];
                if(bus_ex_solver_id.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << "TwoSidesContainer_rxh_A::fillBf_for_PTDF: (SolverId) the line/trafo with id ";
                    exc_ << line_id;
                    exc_ << " is connected (side 2) to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                real_type x = this->fillBf_for_PTDF_coeff(line_id);
                int id_ = this->fillBf_for_PTDF_id(line_id, nb_powerline);
                
                // TODO
                // Bf (nb_branch, nb_bus) : en dc un truc du genre 1 / x / tap for (1..nb_branch, from_bus)
                // and -1. / x / tap for (1..nb_branch, to_bus) 
                if(transpose){
                    Bf.push_back(Eigen::Triplet<real_type> (bus_or_solver_id.cast_int(), id_, 1. / x));
                    Bf.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id.cast_int(), id_, -1. / x));
                }else{
                    Bf.push_back(Eigen::Triplet<real_type> (id_, bus_or_solver_id.cast_int(), 1. / x));
                    Bf.push_back(Eigen::Triplet<real_type> (id_, bus_ex_solver_id.cast_int(), -1. / x));
                }
            }

        }

        // gridmodel utilities
        void reconnect_connected_buses(SubstationContainer & substation) const{
            const Eigen::Index nb_els = nb();
            const std::vector<bool>& status_side_1_ = get_status_side_1();
            const std::vector<bool>& status_side_2_ = get_status_side_2();
            for(Eigen::Index el_id = 0; el_id < nb_els; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;
                
                if(status_side_1_[el_id])
                {
                    const GlobalBusId bus_or_id_me = get_bus_side_1(el_id);        
                    if(bus_or_id_me.cast_int() == _deactivated_bus_id){
                        // TODO DEBUG MODE only this in debug mode
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::reconnect_connected_buses: branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 1) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_xxx_side_1(...)` ?.";
                        throw std::runtime_error(exc_.str());
                    }
                    substation.reconnect_bus(bus_or_id_me);
                }

                if(status_side_2_[el_id])
                {
                    const GlobalBusId bus_ex_id_me = get_bus_side_2(el_id);        
                    if(bus_ex_id_me.cast_int() == _deactivated_bus_id){
                        // TODO DEBUG MODE only this in debug mode
                        std::ostringstream exc_;
                        exc_ << "TwoSidesContainer_rxh_A::reconnect_connected_buses: branch with id ";
                        exc_ << el_id;
                        exc_ << " is connected (side 2) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_xxx_side_2(...)` ?.";
                        throw std::runtime_error(exc_.str());
                    }
                    // bus_status[bus_ex_id_me] = true;
                    substation.reconnect_bus(bus_ex_id_me);
                }
            }
        }

    protected:

        // TODO save that somewhere (and don't forget to template that)
        // with FDPFXB or FPDFBX if saved !
        struct FDPFCoeffs
        {
            real_type yft_bp;
            real_type ytf_bp;
            real_type yff_bp;
            real_type ytt_bp;
            real_type yft_bpp;
            real_type ytf_bpp;
            real_type yff_bpp;
            real_type ytt_bpp;
        };

        StateRes get_tsc_rxha_state() const  // tsc: two sides container
        {
            std::vector<real_type> branch_r(r_.begin(), r_.end());
            std::vector<real_type> branch_x(x_.begin(), x_.end());
            std::vector<cplx_type> branch_h1(h_side_1_.begin(), h_side_1_.end());
            std::vector<cplx_type> branch_h2(h_side_2_.begin(), h_side_2_.end());
            StateRes res(
                get_tsc_state(),
                branch_r,
                branch_x,
                branch_h1,
                branch_h2
            );
            return res;
        }

        void set_tsc_rxha_state(StateRes & my_state)  // tsc: two sides container
        {
            set_tsc_state(std::get<0>(my_state));
            const auto size = nb();

            const std::vector<real_type> & branch_r = std::get<1>(my_state);
            const std::vector<real_type> & branch_x = std::get<2>(my_state);
            const std::vector<cplx_type> & branch_h1 = std::get<3>(my_state);
            const std::vector<cplx_type> & branch_h2 = std::get<4>(my_state);
            check_size(branch_r, size, "branch r");
            check_size(branch_x, size, "branch x");
            check_size(branch_h1, size, "branch h (=g+j.b), side 1");
            check_size(branch_h2, size, "branch h (=g+j.b), side 2");

            r_ = RealVect::Map(&branch_r[0], size);
            x_ = RealVect::Map(&branch_x[0], size);
            h_side_1_ = CplxVect::Map(&branch_h1[0], size);
            h_side_2_ = CplxVect::Map(&branch_h2[0], size);
        }

        void reset_results_tsc_rxha(){
            reset_results_tsc();
            res_a_side_1_ = RealVect(nb());  // in kA
            res_a_side_2_ = RealVect(nb());  // in kA
        }

        virtual void _deactivate(int el_id, SolverControl & solver_control) {
            if(status_global_[el_id]){
                solver_control.tell_recompute_ybus();
                // but sparsity pattern do not change here (possibly one more coeff at 0.)
                solver_control.tell_ybus_some_coeffs_zero();
                solver_control.tell_one_el_changed_bus();  // if the extremity of the line is alone on a bus, this can happen...
            }
        }
        virtual void _reactivate(int el_id, SolverControl & solver_control) {
            if(!status_global_[el_id]){
                solver_control.tell_recompute_ybus();
                solver_control.tell_ybus_change_sparsity_pattern();  // this might change
                solver_control.tell_one_el_changed_bus();  // if the extremity of the line is alone on a bus, this can happen...
            }
        }

        virtual real_type fillBf_for_PTDF_coeff(int el_id) const{
            return x_(el_id);
        }

        virtual int fillBf_for_PTDF_id(int el_id, int nb_powerline) const{
            return el_id;
        }

        virtual FDPFCoeffs get_fdpf_coeffs(int line_id, FDPFMethod xb_or_bx) const{
            FDPFCoeffs res;        
            cplx_type ys_bp, ys_bpp;
            if(xb_or_bx==FDPFMethod::XB){
                ys_bp = 1. / (cplx_type(0., x_(line_id)));
                ys_bpp = 1. / (cplx_type(r_(line_id), x_(line_id)));
            }else if (xb_or_bx==FDPFMethod::BX){
                ys_bp = 1. / (cplx_type(r_(line_id), x_(line_id)));
                ys_bpp = 1. / (cplx_type(0., x_(line_id)));
            }else{
                std::ostringstream exc_;
                exc_ << "fillBp_Bpp: unknown method for the FDPF powerflow for line id ";
                exc_ << line_id;
                throw std::runtime_error(exc_.str());            
            }
            const real_type ys_bp_r = std::imag(ys_bp); 
            res.yff_bp = ys_bp_r;
            res.ytt_bp = ys_bp_r;
            res.yft_bp = -ys_bp_r;
            res.ytf_bp = -ys_bp_r;
            const real_type ys_bpp_r = std::imag(ys_bpp); 
            res.yff_bpp = ys_bpp_r + std::imag(h_side_1_(line_id));
            res.ytt_bpp = ys_bpp_r + std::imag(h_side_2_(line_id));
            res.yft_bpp = -ys_bpp_r;
            res.ytf_bpp = -ys_bpp_r;
            return res;
        }

    protected:        
        // physical properties
        RealVect r_;  // in pu
        RealVect x_;  // in pu
        CplxVect h_side_1_;  // in pu
        CplxVect h_side_2_;  // in pu

        //output data
        RealVect res_a_side_1_;  // in kA
        RealVect res_a_side_2_;  // in kA

        // model coefficients
        CplxVect yac_11_;
        CplxVect yac_12_;
        CplxVect yac_21_;
        CplxVect yac_22_;

        CplxVect ydc_11_;
        CplxVect ydc_12_;
        CplxVect ydc_21_;
        CplxVect ydc_22_;
};
#endif  // TWO_SIDES_CONTAINER_RXH_A_H
