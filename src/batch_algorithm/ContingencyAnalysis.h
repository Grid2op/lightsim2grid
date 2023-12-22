// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SECURITYANALYSIS_H
#define SECURITYANALYSIS_H

#include "BaseBatchSolverSynch.h"
#include <set>

struct Coeff{
    Eigen::Index row_id;
    Eigen::Index col_id;
    cplx_type value;
};

/**
Class to perform a contingency analysis (security analysis), which consist of performing some powerflow after some powerlines
have been disconnected 
 **/
class ContingencyAnalysis: public BaseBatchSolverSynch
{
    public:
        ContingencyAnalysis(const GridModel & init_grid_model):
                            BaseBatchSolverSynch(init_grid_model),
                            _li_defaults(),
                            _li_coeffs(),
                            _timer_total(0.),
                            _timer_modif_Ybus(0.),
                            _timer_pre_proc(0.)
                            { }

        ContingencyAnalysis(const ContingencyAnalysis&) = delete;

        // utilities to add defaults to simulate
        void add_all_n1(){
            for(int l_id = 0; l_id < n_total_; ++l_id){
                std::set<int> this_default = {l_id};
                _li_defaults.insert(this_default);
            }
        }
        void add_n1(int line_id){
            check_ok_el(line_id);
            std::set<int> this_default = {line_id};
            _li_defaults.insert(this_default);
        }
        void add_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s){
                check_ok_el(line_id);
                std::set<int> this_default = {line_id};
                _li_defaults.insert(this_default);
            }
        }
        void add_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk)
            {
                check_ok_el(line_id);
                this_default.insert(line_id);
            }
            _li_defaults.insert(this_default);
        }

        // utilities to remove defaults to simulate (TODO)
        virtual void clear(){
            BaseBatchSolverSynch::clear();
            _li_defaults.clear();
            _li_coeffs.clear();
            _timer_total = 0.;
            _timer_modif_Ybus = 0.;
            _timer_pre_proc = 0.;
        }
        
        bool remove_n1(int line_id){
            check_ok_el(line_id);
            std::set<int> this_default = {line_id};
            auto nb_removed = _li_defaults.erase(this_default);
            return nb_removed >= 1;
        }
        size_t remove_multiple_n1(const std::vector<int> & vect_n1s){
            size_t nb_removed = 0;
            for(const auto line_id : vect_n1s){
                check_ok_el(line_id);
                std::set<int> this_default = {line_id};
                nb_removed += _li_defaults.erase(this_default);
            }
            return nb_removed;
        }
        bool remove_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk)
            {
                check_ok_el(line_id);
                this_default.insert(line_id);
            }
            auto nb_removed = _li_defaults.erase(this_default);
            return nb_removed >= 1;
        }

        // make the computation
        void compute(const CplxVect & Vinit, int max_iter, real_type tol);
        
        Eigen::Ref<const RealMat > compute_flows() {
            compute_flows_from_Vs();
            clean_flows();
            return _amps_flows;
        }

        Eigen::Ref<const RealMat > compute_power_flows() {
            compute_flows_from_Vs(false);
            clean_flows(false);
            return _active_power_flows;
        }

        // python cannot handle set of sets... (c++ can because set are ordered set, but in python
        // they are unordered one, with hash functions)
        const std::set<std::set<int> > & my_defaults() const {return _li_defaults;}
        std::vector<std::vector<int> > my_defaults_vect() const {
            std::vector<std::vector<int> > res;
            res.reserve(_li_defaults.size());
            for(const auto & set_int: _li_defaults){
                std::vector<int> this_def(set_int.begin(), set_int.end());
                res.push_back(this_def);
            }
            return res;
        }

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}
        double modif_Ybus_time() const {return _timer_modif_Ybus;}

    protected:
        // prevent the insertion of "out of range" elements
        void check_ok_el(Eigen::Index el){
            if(el < 0){
                std::ostringstream exc_;
                exc_ << "SecurityAnalysis: cannot add the contingency with id ";
                exc_ << el << " contingency id should be > 0";
                throw std::runtime_error(exc_.str());
            }
            if(el >= n_total_){
                std::ostringstream exc_;
                exc_ << "SecurityAnalysis: cannot add the contingency with id ";
                exc_ << el << " because the grid counts only " << n_total_ << " powerlines / trafos.";
                throw std::runtime_error(exc_.str());
            }
        }
        void init_li_coeffs(bool ac_solver_used);
        // remove the line parameters from Ybus, this is to emulate its disconnection
        bool remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus, const std::vector<Coeff> & coeffs) const;
        // after the coefficient has been removed with "remove_from_Ybus", add it back to Ybus
        void readd_to_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus, const std::vector<Coeff> & coeffs) const;

        // by default the flows are not 0 when the powerline is connected in the original topology
        // this function sorts this out
        void clean_flows(bool is_amps=true);

        // sometimes, when i perform some disconnection, I make the graph non connexe
        // in this case, well, i don't use the results of the simulation
        bool check_invertible(const Eigen::SparseMatrix<cplx_type> & Ybus) const;
    private:
        // li_default
        std::set<std::set<int> > _li_defaults;  // do not use unordered_set here, we rely on the order for different functions !
        std::vector<std::vector<Coeff> > _li_coeffs;  // for each n-k, stores the coefficients I need to modify in the Ybus

        //timers
        double _timer_total;  // total time spent in "compute"
        double _timer_modif_Ybus;  // time to update the Ybus between the defaults simulation
        double _timer_pre_proc;  // time to compute the coefficients of the Ybus
};
#endif  //COMPUTERS_H
