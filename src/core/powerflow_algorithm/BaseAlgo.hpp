// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEALGO_H
#define BASEALGO_H

#include "ls2g_api.hpp"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI
#include <memory>

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.hpp"
#include "CustTimer.hpp"
#include "BaseConstants.hpp"
#include "AlgoConfig.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

namespace ls2g {

class LSGrid;

struct TimerJac {
    double timer_Fx_         = -1.;
    double timer_solve_      = -1.;
    double timer_factor_     = -1.;
    double timer_refactor_   = -1.;
    double timer_initialize_ = -1.;
    double timer_check_      = -1.;
    double timer_dSbus_      = -1.;
    double timer_fillJ_      = -1.;
    double timer_Va_Vm_      = -1.;
    double timer_pre_proc_   = -1.;
    double timer_scale_      = -1.;
    double timer_mismatch_   = -1.;
    double timer_total_nr_   = -1.;
};
using TimerPTDFLODFType = std::tuple<double, double, double>;

/**
This class represents a algorithm to compute powerflow.

It can be derived for different usecase, for example for DC powerflow, AC powerflow using Newton Raphson method etc.
**/
class LS2G_API BaseAlgo : public BaseConstants
{
    public:
        const bool IS_AC;  // should be static ideally...

    public:
        explicit BaseAlgo(bool is_ac=true) noexcept:
            BaseConstants(),
            IS_AC(is_ac),
            n_(-1),
            err_(ErrorType::NotInitError),
            timer_Fx_(0.),
            timer_solve_(0.),
            timer_check_(0.),
            timer_total_nr_(0.),
            lsgrid_ptr_(nullptr){};

        virtual ~BaseAlgo() noexcept = default;

        // no copy allowed
        BaseAlgo(const BaseAlgo&) = delete;
        BaseAlgo(BaseAlgo&&) = delete;
        BaseAlgo & operator=(BaseAlgo&&) = delete;
        BaseAlgo & operator=(const BaseAlgo&) = delete;

        virtual void set_lsgrid(const LSGrid * gridmodel){
            lsgrid_ptr_ = gridmodel;
        }

        virtual Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J() const {
            throw std::runtime_error("AlgorithmSelector::get_J: There is not Jacobian matrix for this solver type.");
        }

        // bus_id -> Jacobian column of that bus' theta / vm / q unknown (-1 if none).
        // Only Newton-Raphson solvers define a Jacobian, hence these throw by default.
        virtual IntVect get_theta_to_J_col_python() const {
            throw std::runtime_error("get_theta_to_J_col: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_to_J_col_python() const {
            throw std::runtime_error("get_vm_to_J_col: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_to_J_col_python() const {
            throw std::runtime_error("get_q_to_J_col: only available for Newton-Raphson solvers.");
        }

        // bus_id -> Jacobian row of that bus' P / Q mismatch equation (-1 if none).
        // The row counterpart of the *_to_J_col maps. Only Newton-Raphson solvers
        // define a Jacobian, hence these throw by default.
        virtual IntVect get_p_to_J_row_python() const {
            throw std::runtime_error("get_p_to_J_row: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_to_J_row_python() const {
            throw std::runtime_error("get_q_to_J_row: only available for Newton-Raphson solvers.");
        }

        // Compact (bus, row/col) registration pair lists -- the row/col
        // counterpart of the *_to_J_col / *_to_J_row bus-keyed maps, but
        // preserving every registration (a bus may appear more than once, or
        // be absent from the bus-keyed map's CURRENT value if a later
        // registration shadowed it there). Only Newton-Raphson solvers define
        // a Jacobian, hence these throw by default.
        virtual IntVect get_p_buses_python() const {
            throw std::runtime_error("get_p_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_p_rows_python() const {
            throw std::runtime_error("get_p_rows: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_buses_python() const {
            throw std::runtime_error("get_q_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_rows_python() const {
            throw std::runtime_error("get_q_rows: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_theta_buses_python() const {
            throw std::runtime_error("get_theta_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_theta_cols_python() const {
            throw std::runtime_error("get_theta_cols: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_buses_python() const {
            throw std::runtime_error("get_vm_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_cols_python() const {
            throw std::runtime_error("get_vm_cols: only available for Newton-Raphson solvers.");
        }

        // VoltageControl (remote gen + SVC) converged reactive injection per
        // controller (pu) and its (kind, element id) identity, in controller
        // registration order. Empty for algorithms without the extension (DC,
        // Gauss-Seidel, NR without any voltage-mode controller).
        virtual RealVect get_controller_q()       const { return RealVect(); }
        virtual IntVect  get_controller_kind()    const { return IntVect(); }
        virtual IntVect  get_controller_elem_id() const { return IntVect(); }

        // MultiSlack: J column of the slack_absorbed unknown (-1 when the
        // distributed-slack-in-Jacobian extension is not active).
        virtual int      get_slack_col()          const { return -1; }
        // MultiSlack: converged VALUE of the slack_absorbed unknown (pu; 0
        // when the extension is not active). NOT the same as the per-solve
        // initial guess of 0 -- this is the ground truth after convergence.
        virtual real_type get_slack_absorbed()     const { return static_cast<real_type>(0.); }

        Eigen::Ref<const RealVect> get_Va() const{
            return Va_;
        }
        Eigen::Ref<const RealVect> get_Vm() const{
            return Vm_;
        }
        Eigen::Ref<const CplxVect> get_V() const{
            return V_;
        }
        ErrorType get_error() const {
            return err_;
        }
        int get_nb_iter() const {
            return nr_iter_;
        }

        bool converged() const{
            return err_ == ErrorType::NoError;
        }

        std::tuple<double, double, double, double> get_timers() const
        {
            // TODO change the order of the timers here!
            auto res = std::tuple<double, double, double, double>(
              timer_Fx_, timer_solve_, timer_check_, timer_total_nr_);
            return res;
        }
        
        virtual TimerJac get_timers_jacobian() const
        {
            TimerJac res;
            res.timer_Fx_       = timer_Fx_;
            res.timer_solve_    = timer_solve_;
            res.timer_check_    = timer_check_;
            res.timer_total_nr_ = timer_total_nr_;
            return res;
        }
        
        virtual TimerPTDFLODFType get_timers_ptdf_lodf() const
        {
            TimerPTDFLODFType res = {
                -1.,  // not available for non NR solver, so I put -1
                -1.,  // not available for non NR solver, so I put -1
                -1.,  // not available for non NR solver, so I put -1
            };
            return res;
        }

        // Complex AC entry point: solves V . (Ybus . V)* = Sbus.
        // Every AC solver overrides this; the DC solver does not (it uses `compute_pf_dc`) and
        // therefore inherits this throwing default (symmetric with `compute_pf_dc` below).
        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & /*Ybus*/,
                        CplxVect & /*V*/,  // store the results of the powerflow and the Vinit !
                        const CplxVect & /*Sbus*/,
                        Eigen::Ref<const IntVect> /*slack_ids*/,
                        const RealVect & /*slack_weights*/,
                        Eigen::Ref<const IntVect> /*pv*/,
                        Eigen::Ref<const IntVect> /*pq*/,
                        int /*max_iter*/,
                        real_type /*tol*/
                        ){
            throw std::runtime_error("compute_pf (complex AC entry point) is not available for this solver (DC solvers use compute_pf_dc).");
        }

        // Native real-valued DC entry point: DC only needs `Bbus . theta = Pbus` (all real).
        // Only the DC solver overrides this; every other solver type throws.
        // `V` carries the complex initial voltage (slack angle + voltage setpoints) on input
        // and the complex result on output. There is no max_iter / tol: DC is a single linear solve.
        virtual
        bool compute_pf_dc(const Eigen::SparseMatrix<real_type> & /*Bbus*/,
                           CplxVect & /*V*/,
                           const RealVect & /*Pbus*/,
                           Eigen::Ref<const IntVect> /*slack_ids*/,
                           const RealVect & /*slack_weights*/,
                           Eigen::Ref<const IntVect> /*pv*/,
                           Eigen::Ref<const IntVect> /*pq*/){
            throw std::runtime_error("compute_pf_dc is only available for DC solvers.");
        }

        void tell_solver_control(const AlgoControl & solver_control){
            _solver_control = solver_control;
        }
        virtual void reset();
        virtual RealMat get_ptdf(){
            throw std::runtime_error("Impossible to get the PTDF matrix with this solver type.");
        }
        virtual RealMat get_lodf(const IntVect & /*from_bus*/,
                                 const IntVect & /*to_bus*/){  // TODO interface is likely to change
            throw std::runtime_error("Impossible to get the LODF matrix with this solver type.");
        }
        virtual Eigen::SparseMatrix<real_type> get_bsdf(){  // TODO interface is likely to change
            throw std::runtime_error("Impossible to get the BSDF matrix with this solver type.");
        }

        virtual void update_internal_Ybus(const Coeff & /*new_coeffs*/, bool /*add*/){
            throw std::runtime_error("Function update_internal_Ybus not implemented in general.");
        }

        // bus masking (ContingencyAnalysis "handle disconnected grid" mode): forces
        // the given solver buses' equations to identity so an isolated island does
        // not make the system singular, without changing the matrix sparsity. Only
        // the Newton-Raphson family supports it (see supports_bus_masking); the
        // default is a no-op so other algorithms are unaffected.
        virtual bool supports_bus_masking() const { return false; }
        virtual void set_masked_buses(const std::vector<int> & /*solver_bus_ids*/) {}

        virtual AlgoConfig get_config() const { return AlgoConfig{}; }
        virtual void set_config(const AlgoConfig&) {}
        
    protected:
        virtual void reset_timer(){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_check_ = 0.;
            timer_total_nr_ = 0.;
        }

        bool is_linear_solver_valid(){
            // bool res = true;
            // if((err_ == ErrorType::NotInitError) || (err_ == ErrorType::LicenseError)) res = false;  // cannot use a non intialize solver
            // return res;
            return (err_ != ErrorType::LicenseError);
        }
        RealVect _evaluate_Fx(const Eigen::SparseMatrix<cplx_type> &  Ybus,
                              const CplxVect & V,
                              const CplxVect & Sbus,
                              size_t slack_id,  // id of the slack bus
                              real_type slack_absorbed,
                              const RealVect & slack_weights,
                              Eigen::Ref<const IntVect> pv,
                              Eigen::Ref<const IntVect> pq);

        RealVect _evaluate_Fx(const Eigen::SparseMatrix<cplx_type> &  Ybus,
                              const CplxVect & V,
                              const CplxVect & Sbus,
                              Eigen::Ref<const IntVect> pv,
                              Eigen::Ref<const IntVect> pq);

        bool _check_for_convergence(const RealVect & F,
                                    real_type tol);

        bool _check_for_convergence(const RealVect & p,
                                    const RealVect & q,
                                    real_type tol);

        Eigen::VectorXi extract_slack_bus_id(Eigen::Ref<const IntVect> pv,
                                             Eigen::Ref<const IntVect> pq,
                                             unsigned int nb_bus);

        /**
        When there are multiple slacks, add the other "slack buses" in the PV buses indexes
        (behaves as if only the first element is used for the slack !!!, called "ref slack")
        **/
        Eigen::VectorXi retrieve_pv_with_slack(Eigen::Ref<const IntVect> slack_ids,
                                               Eigen::Ref<const IntVect> pv) const {
            if(slack_ids.size() > 1){
                const auto nb_slack_added = slack_ids.size() - 1;
                Eigen::VectorXi my_pv = Eigen::VectorXi(pv.size() + nb_slack_added);
                for(auto i = 0; i < nb_slack_added; ++i){
                    my_pv(i) = slack_ids(i+1);
                }
                for(auto i = 0; i < pv.size(); ++i){
                    my_pv(i + nb_slack_added) = pv(i);
                }
                return my_pv;
            }else{
                return pv;
            }
        }

        /**
        When there are multiple slacks, add the other "slack buses" in the PV buses indexes
        **/
        Eigen::VectorXi add_slack_to_pv(Eigen::Ref<const IntVect> slack_ids,
                                        Eigen::Ref<const IntVect> pv) const {
            Eigen::VectorXi my_pv = Eigen::VectorXi(slack_ids.size() + pv.size());
            my_pv << slack_ids, pv;
            return my_pv;
        }
        
        // terribly inefficient way to know if an element is in a vector
        bool isin(int k, const Eigen::VectorXi vect) const{
            for(auto el : vect){
                if(el == k) return true;
            }
            return false;
        }

        void get_Bf(Eigen::SparseMatrix<real_type> & Bf) const;
        void get_Bf_transpose(Eigen::SparseMatrix<real_type> & Bf_T) const;
        
    protected:
        // solver initialization
        int n_;

        // solution of the problem
        RealVect Vm_;  // voltage magnitude
        RealVect Va_;  // voltage angle
        CplxVect V_;  // complex voltage

        int nr_iter_;  // number of iteration performs by the solver (may vary depending on the solver)
        ErrorType err_; //error message:
        // -1 : the solver has not been initialized (call initialize in this case)
        // 0 everything ok
        // 1: i can't factorize the matrix (klu_factor)
        // 2: i can't refactorize the matrix (klu_refactor)
        // 3: i can't solve the system (klu_solve)
        // 4: end of possible iterations (divergence because nr_iter_ >= max_iter)

        // timers
        double timer_Fx_;
        double timer_solve_;
        double timer_check_;
        double timer_total_nr_;

        const LSGrid * lsgrid_ptr_;  // does not have ownership so that's fine (pointer to the base gridmodel, can be used for some powerflow)
        AlgoControl _solver_control;

};


} // namespace ls2g

#endif // BASEALGO_H