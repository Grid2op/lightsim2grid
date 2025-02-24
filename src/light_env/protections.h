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
#include <vector>


class Protections
{
    protected:
        typedef Eigen::Array<bool, Eigen::Dynamic, 1 > VectorXb;

    public:
        // TODO
        Protections():
            has_been_checked_(false),
            nb_iter_(0){};

        // input getters
        Eigen::Ref<const RealVect> get_thermal_limit_or() const {return thermal_limit_or_;}
        Eigen::Ref<const RealVect> get_thermal_limit_ex() const {return thermal_limit_ex_;}
        Eigen::Ref<const Eigen::VectorXi> get_max_line_time_step_overflow() const {return max_line_time_step_overflow_;}

        // output getters
        Eigen::Ref<const RealVect> get_rho() const {return rho_;}
        int get_nb_iter() const {return nb_iter_;}

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

            line_time_step_overflow_ = Eigen::VectorXi::Zero(nb_line_trafo_);
            rho_  = RealVect::Zero(nb_line_trafo_);
            nb_iter_ = 0;
            has_been_checked_ = true;
        }

        void reset_cooldowns(const GridModel & grid ){
            line_time_step_overflow_ = Eigen::VectorXi::Zero(nb_line_trafo_);
            nb_iter_ = 0;
        }

        CplxVect next_grid_state(GridModel & grid, int max_iter, float tol) {
            // NB: all of this occurs withing one env step
            VectorXb already_overflow_this = VectorXb::Zero(nb_line_trafo_);
            bool need_another = false;
            nb_iter_ = 0;
            line_disconnected_.clear();

            CplxVect res;
            while(true){
                need_another = false;
                nb_iter_ += 1;
                std::vector<int> line_id_disco_this_iter;
                res = run_powerflow(grid, max_iter, tol);
                if(res.size() == 0){
                    // divergence
                    return res;
                }

                VectorXb is_overflow = check_overflow(grid);

                // update overflow
                for(int br_id = 0; br_id < nb_line_trafo_; ++br_id){
                    if(is_overflow(br_id) && !already_overflow_this(br_id)){
                        // it's a "new" overflow, so i increase the counter and
                        // tell it's on overflow this step
                        line_time_step_overflow_(br_id) += 1;
                        already_overflow_this(br_id) = true;
                    }
                }

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
                        grid.deactivate_trafo(tr_id);
                        need_another = true;
                        line_id_disco_this_iter.push_back(tr_id);
                    }
                }

                if(!need_another){
                    // no powerfline disconnection
                    // I can terminate
                    return res;
                }
                line_disconnected_.push_back(line_id_disco_this_iter);
            }
            // TODO update overflow: set them to 0 if not on overflow after the routine
        }

        CplxVect run_powerflow(GridModel & grid, int max_iter, float tol) const{
            CplxVect v_init_dc = 1.04 * grid.get_bus_vn_kv();
            grid.deactivate_result_computation();
            CplxVect v_init_ac = grid.dc_pf(v_init_dc, max_iter, tol);
            grid.reactivate_result_computation();
            if(v_init_ac.size() == 0) return v_init_ac;  // TODO divergence
            CplxVect v_res_ac = grid.ac_pf(v_init_ac, max_iter, tol);
            return v_res_ac;
        }
    
    protected:
        VectorXb check_overflow(const GridModel & grid ){
            update_rho(grid);
            return rho_.array() >= 1.;
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

            const auto & vectl_or = std::get<3>(grid.get_lineor_res());
            const auto & vectt_or = std::get<3>(grid.get_trafohv_res());

            const auto & vectl_ex = std::get<3>(grid.get_lineex_res());
            const auto & vectt_ex = std::get<3>(grid.get_trafolv_res());

            const auto & th_lim_or = get_thermal_limit_or();
            const auto & th_lim_ex = get_thermal_limit_ex();

            RealVect vect_or;
            vect_or << vectl_or, vectt_or;  // TODO maybe an error if not evaluated on the right order...
            RealVect vect_ex;
            vect_ex << vectl_ex, vectt_ex;

            rho_.array() = (vect_or.array() / th_lim_or.array()).cwiseMax(vect_ex.array() / th_lim_ex.array());
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

        // extra information
        int nb_iter_;
        
};

#endif // PROTECTIONS_H
