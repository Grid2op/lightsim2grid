// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LIGHT_ENV_H
#define LIGHT_ENV_H

#include "GridModel.h"
#include "topo_action.h"
#include "inj_action.h"
#include "protections.h"

#include <unordered_map>


class LightEnv
{
        typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RealMat;
        // typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CplxMat;

    public:
        typedef std::unordered_map<std::string, std::string> InfoReturnedType;
        typedef std::tuple<RealVect, double, bool, bool,  InfoReturnedType > StepReturnedType;
        typedef std::tuple<RealVect, InfoReturnedType > ResetReturnedType;

        LightEnv(const GridModel & gridmodel):
            has_been_checked_(false),
            grid_(gridmodel)
            {};

        void assign_time_series(const RealMat & load_p,
                                const RealMat & load_q,
                                const RealMat & gen_p,
                                const RealMat & gen_v,
                                const RealMat & storage_p,
                                const RealMat & shunt_p,
                                const RealMat & shunt_q,
                                const RealMat & sgen_p,
                                const RealMat & sgen_q){
            has_been_checked_ = false;
            load_p_ = load_p;
            load_q_ = load_q;
            gen_p_ = gen_p;
            gen_v_ = gen_v;
            storage_p_ = storage_p;
            shunt_p_ = shunt_p;
            shunt_q_ = shunt_q;
            sgen_p_ = sgen_p;
            sgen_q_ = sgen_q;
        }

        void assign_protections(const Protections & protections){
            has_been_checked_ = false;
            protections_ = protections;
        }

        const Protections & get_protections() const {return protections_;}

        ResetReturnedType reset(){
            if(!has_been_checked_){
                perform_internal_checks();
                has_been_checked_ = true;
            }

            current_step_ = 0;
            info_ = InfoReturnedType();

            // reset the cooldown for the action
            // and reset the protections
            reset_cooldowns();

            // TODO assign grid to init topo : init_topo :  topo_action
            // init_topo.apply_to_gridmodel(grid_);

            // set the correct injection to the grid
            const InjAction inj_action(load_p_.row(current_step_),
                                       load_q_.row(current_step_),
                                       gen_p_.row(current_step_),
                                       gen_v_.row(current_step_),
                                       storage_p_.row(current_step_),
                                       shunt_p_.row(current_step_),
                                       shunt_q_.row(current_step_),
                                       sgen_p_.row(current_step_),
                                       sgen_q_.row(current_step_)
                                       );
            inj_action.apply_to_gridmodel(grid_);

            // run the powerflow
            V_ = protections_.run_powerflow(grid_, max_iter_, tol_);
            // TODO if divergence !!!!

            // extract the observation
            extract_observation();

            // return the value
            return ResetReturnedType(obs_, info_);
        }

        StepReturnedType step(int act_id){
            if(!has_been_checked_){
                throw std::runtime_error("Environment cannot be used, you most likely need to call env.reset() ");
            }
            // TODO act_id is not used atm !

            current_step_ += 1;
            info_ = InfoReturnedType();
            if (current_step_ >= max_step_){
                info_["success"] = "true";
                info_["failure"] = "false";
                info_["survival_time"] = std::to_string(current_step_ / max_step_);
                obs_ = RealVect::Zero(obs_.size());
                has_been_checked_ = false;
                return StepReturnedType(obs_, 1., true, false, info_);
            }
            // apply the topology
            // TODO 
            // TODO check cooldowns too !

            // apply the injection
            const InjAction inj_action(load_p_.row(current_step_),
                                       load_q_.row(current_step_),
                                       gen_p_.row(current_step_),
                                       gen_v_.row(current_step_),
                                       storage_p_.row(current_step_),
                                       shunt_p_.row(current_step_),
                                       shunt_q_.row(current_step_),
                                       sgen_p_.row(current_step_),
                                       sgen_q_.row(current_step_)
                                       );
            inj_action.apply_to_gridmodel(grid_);
            // TODO apply redispatching, storage, curtailment and other !

            V_ = protections_.next_grid_state(grid_, max_iter_, tol_);
            if(V_.size() == 0){
                // divergence
                info_["success"] = "false";
                info_["failure"] = "true";
                info_["survival_time"] = std::to_string(current_step_ / max_step_);
                obs_ = RealVect::Zero(obs_.size());
                has_been_checked_ = false;
                return StepReturnedType(obs_, 0., true, true, info_);
            }

            // extract the observation
            extract_observation();

            // update time dependant information
            // TODO decrease cooldown not used
            // TODO increase cooldown if disconnected due to overflow


            // return the value
            return StepReturnedType(obs_, current_step_ / max_iter_, false, false, info_);
        }

    protected:

        void extract_observation(){
            obs_ = protections_.get_rho();
        }

        void reset_cooldowns(){
            time_step_sub_cooldown_ = Eigen::VectorXi::Zero(grid_.get_n_sub());
            time_step_line_cooldown_ = Eigen::VectorXi::Zero(grid_.nb_powerline() + grid_.nb_trafo());
            protections_.reset_cooldowns(grid_);
        }

        void perform_internal_checks(){
            // correct number of rows (steps)
            const int n_ts = load_p_.rows();
            aux_check_row(load_p_, n_ts, "perform_internal_checks (load_p)");
            aux_check_row(load_q_, n_ts, "perform_internal_checks (load_q)");
            aux_check_row(gen_p_, n_ts, "perform_internal_checks (gen_p)");
            aux_check_row(gen_v_, n_ts, "perform_internal_checks (gen_v)");
            aux_check_row(storage_p_, n_ts, "perform_internal_checks (storage_p)");
            aux_check_row(shunt_p_, n_ts, "perform_internal_checks (shunt_p)");
            aux_check_row(shunt_q_, n_ts, "perform_internal_checks (shunt_q)");
            aux_check_row(sgen_p_, n_ts, "perform_internal_checks (sgen_p)");
            aux_check_row(sgen_q_, n_ts, "perform_internal_checks (sgen_q)");
            max_step_ = n_ts;

            // correct number of columns (number of elelments)
            aux_check_col(load_p_, grid_.get_loads().nb(), "perform_internal_checks (load_p)");
            aux_check_col(load_q_, grid_.get_loads().nb(), "perform_internal_checks (load_q)");
            aux_check_col(gen_p_, grid_.get_generators().nb(), "perform_internal_checks (gen_p)");
            aux_check_col(gen_v_, grid_.get_generators().nb(), "perform_internal_checks (gen_v)");
            aux_check_col(storage_p_, grid_.get_storages().nb(), "perform_internal_checks (storage_p)");
            aux_check_col(shunt_p_, grid_.get_shunts().nb(), "perform_internal_checks (shunt_p)");
            aux_check_col(shunt_q_, grid_.get_shunts().nb(), "perform_internal_checks (shunt_q)");
            aux_check_row(sgen_p_, grid_.get_static_generators().nb(), "perform_internal_checks (sgen_p)");
            aux_check_row(sgen_q_, grid_.get_static_generators().nb(), "perform_internal_checks (sgen_q)");

            // TODO topo_actions

            // cooldowns
            reset_cooldowns();

            // protections
            protections_.check_validity(grid_);
        }

        template<class EigenType>
        void aux_check_row(const EigenType & new_vect,
                           int size_th,
                           const std::string & error_detailed) const{
            if(new_vect.rows() != size_th)
            {
                std::ostringstream exc_;
                exc_ << "LightEnv::" << error_detailed << ".";
                exc_ << "Theorical size is: ";
                exc_ << size_th;
                exc_ << " but your provided data with ";
                exc_ << new_vect.rows() << " rows.";
                throw std::runtime_error(exc_.str());
            }
        }
        template<class EigenType>
        void aux_check_col(const EigenType & new_vect,
                           int size_th,
                           const std::string & error_detailed) const{
            if(new_vect.cols() != size_th)
            {
                std::ostringstream exc_;
                exc_ << "LightEnv::" << error_detailed << ".";
                exc_ << "Theorical size is: ";
                exc_ << size_th;
                exc_ << " but your provided data with ";
                exc_ << new_vect.rows() << " columns.";
                throw std::runtime_error(exc_.str());
            }
        }

    protected:
        // consistency
        bool has_been_checked_;

        // time series
        RealMat load_p_;
        RealMat load_q_;
        RealMat gen_p_;
        RealMat gen_v_;
        RealMat storage_p_;
        RealMat shunt_p_;
        RealMat shunt_q_;
        RealMat sgen_p_;
        RealMat sgen_q_;
        int max_step_;  // size of the input data

        // protections and thermal limits
        Protections protections_;

        // opponent

        // powergrid state
        GridModel grid_;
        CplxVect V_;
        int max_iter_;
        float tol_;

        // time related information
        int current_step_;
        Eigen::VectorXi time_step_sub_cooldown_;
        Eigen::VectorXi time_step_line_cooldown_;

        // the observation
        RealVect obs_;
        std::unordered_map<std::string, std::string> info_;

};
#endif // LIGHT_ENV_H
