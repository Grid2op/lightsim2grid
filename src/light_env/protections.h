// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef PROTECTIONS_H
#define PROTECTIONS_H

#include "GridModel.h"
#include "CustTimer.h"
#include <vector>


class Protections
{
    protected:
        typedef Eigen::Array<bool, Eigen::Dynamic, 1 > VectorXb;

    public:
        // TODO
        Protections():
            has_been_checked_(false),
            nb_iter_(0),
            timer_total_(0.),
            timer_check_overflow_(0.),
            timer_powerflow_(0.),
            timer_update_rho_(0.){}

        // input getters
        Eigen::Ref<const RealVect> get_thermal_limit_or() const {return thermal_limit_or_;}
        Eigen::Ref<const RealVect> get_thermal_limit_ex() const {return thermal_limit_ex_;}
        Eigen::Ref<const Eigen::VectorXi> get_max_line_time_step_overflow() const {return max_line_time_step_overflow_;}

        // output getters
        Eigen::Ref<const RealVect> get_rho() const {return rho_;}
        Eigen::Ref<const Eigen::VectorXi> get_line_time_step_overflow() const {return line_time_step_overflow_;}
        Eigen::Ref<const CplxVect> get_v_init_dc() const {return v_init_dc_;}
        Eigen::Ref<const CplxVect> get_v_init_ac() const {return v_init_ac_;}
        Eigen::Ref<const CplxVect> get_v_res_ac() const {return v_res_ac_;}
        int get_nb_iter() const {return nb_iter_;}
        const std::vector< std::vector<int> > & get_line_disconnected() const {return line_disconnected_;}
        double get_total_time() const {return timer_total_;}
        double get_check_overflow_time() const {return timer_check_overflow_;}
        double get_powerflow_time() const {return timer_powerflow_;}
        double get_update_rho_time() const {return timer_update_rho_;}

        // setters
        void set_thermal_limit_or(Eigen::Ref<RealVect> thermal_limit_or){
            if(has_been_checked_){
                // I need to maintain consistency with the grid size
                aux_check_size(thermal_limit_or, "set_thermal_limit_or");
            }
            thermal_limit_or_ = thermal_limit_or;
        } 
        void set_thermal_limit_ex(Eigen::Ref<RealVect> thermal_limit_ex){
            if(has_been_checked_){
                // I need to maintain consistency with the grid size
                aux_check_size(thermal_limit_ex, "set_thermal_limit_ex");
            }
            thermal_limit_ex_ = thermal_limit_ex;
        } 
        void set_max_line_time_step_overflow(Eigen::Ref<Eigen::VectorXi> max_line_time_step_overflow){
            if(has_been_checked_){
                // I need to maintain consistency with the grid size
                aux_check_size(max_line_time_step_overflow, "set_max_line_time_step_overflow");
            }
            max_line_time_step_overflow_ = max_line_time_step_overflow;
        }

        // bellow: use only inside light_env

        // TODO check compliance (correct number of elements etc.)
        void check_validity(const GridModel & grid) {
            // TODO really perform the check
            nb_line_ = grid.nb_powerline();
            nb_line_trafo_ = nb_line_ + grid.nb_trafo();

            aux_check_size(max_line_time_step_overflow_, "check_validity (max_line_time_step_overflow)");
            aux_check_size(thermal_limit_or_, "check_validity (thermal_limit_or)");
            aux_check_size(thermal_limit_ex_, "check_validity (thermal_limit_ex)");

            rho_  = RealVect::Zero(nb_line_trafo_);
            reset_cooldowns();
            has_been_checked_ = true;
        }

        CplxVect next_grid_state(GridModel & grid, int max_iter, float tol) {
            auto timer_total = CustTimer();
            // NB: all of this occurs withing one env step
            VectorXb already_overflow_this = VectorXb::Zero(nb_line_trafo_);
            bool need_another = false;
            nb_iter_ = 0;
            line_disconnected_.clear();
            // std::cout << "after init\n";
            CplxVect res;
            while(true){
                // std::cout << "\tOne more iter\n";
                need_another = false;
                nb_iter_ += 1;
                std::vector<int> line_id_disco_this_iter;
                res = run_powerflow(grid, max_iter, tol);
                // std::cout << "\t\t after powerflow\n";
                if(res.size() == 0){
                    // divergence
                    timer_total_ += timer_total.duration();
                    return res;
                }

                VectorXb is_overflow = check_overflow(grid);
                // std::cout << "\t\t after check_overflow\n";

                // update overflow
                for(int br_id = 0; br_id < nb_line_trafo_; ++br_id){
                    if(is_overflow(br_id) && !already_overflow_this(br_id)){
                        // it's a "new" overflow, so i increase the counter and
                        // tell it's on overflow this step
                        line_time_step_overflow_(br_id) += 1;
                        already_overflow_this(br_id) = true;
                    }
                }
                // std::cout << "\t\t after update overflow\n";

                // disconnect line / trafo if needed
                for(int line_id=0; line_id < nb_line_; ++line_id){
                    if(line_time_step_overflow_(line_id) > max_line_time_step_overflow_(line_id)){
                        // line has been overflow for too much
                        grid.deactivate_powerline(line_id);
                        need_another = true;
                        line_id_disco_this_iter.push_back(line_id);
                    }
                }
                for(int tr_id=nb_line_; tr_id < nb_line_trafo_; ++tr_id){
                    if(line_time_step_overflow_(tr_id) > max_line_time_step_overflow_(tr_id)){
                        // line has been overflow for too much
                        grid.deactivate_trafo(tr_id - nb_line_);
                        need_another = true;
                        line_id_disco_this_iter.push_back(tr_id);
                    }
                }
                // std::cout << "\t\t after disconnect line / trafo \n";

                if(!need_another){
                    // no powerfline disconnection
                    // I can terminate
                    timer_total_ += timer_total.duration();
                    return res;
                }
                line_disconnected_.push_back(line_id_disco_this_iter);
            }
            // TODO update overflow: set them to 0 if not on overflow after the routine
        }

        CplxVect run_powerflow(GridModel & grid, int max_iter, float tol){
            auto timer_powerflows = CustTimer();
            v_init_dc_ = CplxVect::Constant(grid.get_bus_vn_kv().size(), {1.04, 0.});
            grid.deactivate_result_computation();
            v_init_ac_ = grid.dc_pf(v_init_dc_, max_iter, tol);
            grid.reactivate_result_computation();
            if(v_init_ac_.size() == 0) return v_init_ac_;  // TODO divergence
            v_res_ac_ = grid.ac_pf(v_init_ac_, max_iter, tol);
            grid.unset_changes();
            timer_powerflow_ += timer_powerflows.duration();
            return v_res_ac_;
        }
    
    protected:
        VectorXb check_overflow(const GridModel & grid ){
            auto timer_ov = CustTimer();
            update_rho(grid);
            VectorXb res = rho_.array() >= 1.;
            timer_check_overflow_ += timer_ov.duration();
            return res;
        }

        void reset_timers(){
            timer_total_ = 0.;
            timer_check_overflow_ = 0.;
            timer_powerflow_ = 0.;
        }
        
        void reset_cooldowns(){
            line_time_step_overflow_ = Eigen::VectorXi::Zero(nb_line_trafo_);
            nb_iter_ = 0;
            reset_timers();
        }

        template<class EigenType>
        void aux_check_size(const EigenType & new_vect,
                            const std::string & fun_name) const{
            if(new_vect.size() != nb_line_trafo_)
            {
                std::ostringstream exc_;
                exc_ << "Protections::" << fun_name << ": wrong vector size. Please a vector having the size of the number";
                exc_ << "of powerlines in the grid: ";
                exc_ << nb_line_trafo_;
                exc_ << ". You provided a vector with ";
                exc_ << new_vect.size();
                exc_ << ". It is not possible. Please enter new protection data.";
                throw std::runtime_error(exc_.str());
            }
        }

        void update_rho(const GridModel & grid){
            auto timer_update_rho = CustTimer();
            Eigen::Ref<const RealVect> vectl_or = std::get<3>(grid.get_lineor_res());
            Eigen::Ref<const RealVect> vectt_or = std::get<3>(grid.get_trafohv_res());
            // std::cout << "\t\t\t after get or side\n";

            Eigen::Ref<const RealVect> vectl_ex = std::get<3>(grid.get_lineex_res());
            Eigen::Ref<const RealVect> vectt_ex = std::get<3>(grid.get_trafolv_res());
            // std::cout << "\t\t\t after get ex side\n";

            const auto & th_lim_or = get_thermal_limit_or();
            const auto & th_lim_ex = get_thermal_limit_ex();
            // std::cout << "\t\t\t after get_thermal_limit\n";
            RealVect vect_or(nb_line_trafo_);
            vect_or << vectl_or, vectt_or;  // TODO maybe an error if not evaluated on the right order...
            RealVect vect_ex(nb_line_trafo_);
            vect_ex << vectl_ex, vectt_ex;
            // std::cout << "\t\t\t after the << \n";

            rho_.array() = (vect_or.array() / th_lim_or.array()).cwiseMax(vect_ex.array() / th_lim_ex.array());
            // std::cout << "\t\t\t after rho update \n";
            timer_update_rho_ = timer_update_rho.duration();
        }

    protected:
        bool has_been_checked_;
        int nb_line_;
        int nb_line_trafo_;

        // real data
        RealVect thermal_limit_or_;
        RealVect thermal_limit_ex_;
        Eigen::VectorXi max_line_time_step_overflow_;

        // kind of output data
        Eigen::VectorXi line_time_step_overflow_;
        RealVect rho_;
        std::vector< std::vector<int> > line_disconnected_;  // line_disconnected_[iter_id] = [line_id]
        CplxVect v_init_dc_;
        CplxVect v_init_ac_;
        CplxVect v_res_ac_;

        // extra information
        int nb_iter_;

        // timers
        double timer_total_;
        double timer_check_overflow_;
        double timer_powerflow_;
        double timer_update_rho_;
        
};

#endif // PROTECTIONS_H
