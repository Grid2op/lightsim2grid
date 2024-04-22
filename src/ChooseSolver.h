// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CHOOSESOLVER_H
#define CHOOSESOLVER_H

#include<vector>

// import newton raphson solvers using different linear algebra solvers
#include "Solvers.h"

enum class SolverType {SparseLU, KLU, GaussSeidel, DC, GaussSeidelSynch, NICSLU, 
                       SparseLUSingleSlack, KLUSingleSlack, NICSLUSingleSlack, 
                       KLUDC, NICSLUDC,
                       CKTSO, CKTSOSingleSlack, CKTSODC,
                       FDPF_XB_SparseLU, FDPF_BX_SparseLU, // from 0.7.5
                       FDPF_XB_KLU, FDPF_BX_KLU,  // from 0.7.5
                       FDPF_XB_NICSLU, FDPF_BX_NICSLU,  // from 0.7.5
                       FDPF_XB_CKTSO,  FDPF_BX_CKTSO  // from 0.7.5
                       };


std::ostream& operator<<(std::ostream& out, const SolverType& solver_type);
// TODO define a template class instead of these weird stuff !!!
// TODO export all methods from base class !


// NB: when adding a new solver, you need to specialize the *tmp method (eg get_Va_tmp)
// and also to "forward" the specialisation (adding the if(solvertype==XXX)) in the compute_pf, get_V, get_J, get_Va, get_Vm
// and the "available_solvers" (add it to the list)
// and to add a attribute with the proper class and the reset method
class ChooseSolver
{
    public:
         ChooseSolver():
             _solver_type(SolverType::SparseLU),
             _type_used_for_nr(SolverType::SparseLU)
             {};

        std::vector<SolverType> available_solvers() const
        {
            std::vector<SolverType> res;
            res.reserve(14);

            res.push_back(SolverType::SparseLU);
            res.push_back(SolverType::GaussSeidel);
            res.push_back(SolverType::DC);
            res.push_back(SolverType::GaussSeidelSynch);
            res.push_back(SolverType::SparseLUSingleSlack);
            res.push_back(SolverType::FDPF_XB_SparseLU);
            res.push_back(SolverType::FDPF_BX_SparseLU);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
                res.push_back(SolverType::KLUSingleSlack);
                res.push_back(SolverType::KLUDC);
                res.push_back(SolverType::FDPF_XB_KLU);
                res.push_back(SolverType::FDPF_BX_KLU);
            #endif
            #ifdef NICSLU_SOLVER_AVAILABLE
                res.push_back(SolverType::NICSLU);
                res.push_back(SolverType::NICSLUSingleSlack);
                res.push_back(SolverType::NICSLUDC);
                res.push_back(SolverType::FDPF_XB_NICSLU);
                res.push_back(SolverType::FDPF_BX_NICSLU);
            #endif
            #ifdef CKTSO_SOLVER_AVAILABLE
                res.push_back(SolverType::CKTSO);
                res.push_back(SolverType::CKTSOSingleSlack);
                res.push_back(SolverType::CKTSODC);
                res.push_back(SolverType::FDPF_XB_CKTSO);
                res.push_back(SolverType::FDPF_BX_CKTSO);
            #endif
            return res;
        }
        
        bool is_dc(const SolverType & type){
            bool res;
            res = (type == SolverType::DC) || 
                  (type == SolverType::KLUDC) || 
                  (type == SolverType::NICSLUDC) ||
                  (type == SolverType::CKTSODC);
            return res;
        }
        SolverType get_type() const {return _solver_type;}

        // TODO DO that for all solvers
        FDPF_XB_SparseLUSolver & get_fdpf_xb_lu() { return _solver_fdpf_xb_lu;}
        FDPF_BX_SparseLUSolver & get_fdpf_bx_lu() { return _solver_fdpf_bx_lu;}

        void set_gridmodel(const GridModel * gridmodel){
            _solver_lu.set_gridmodel(gridmodel);
            _solver_lu_single.set_gridmodel(gridmodel);
            _solver_gaussseidel.set_gridmodel(gridmodel);
            _solver_gaussseidelsynch.set_gridmodel(gridmodel);
            _solver_dc.set_gridmodel(gridmodel);
            _solver_fdpf_xb_lu.set_gridmodel(gridmodel);
            _solver_fdpf_bx_lu.set_gridmodel(gridmodel);
            #ifdef KLU_SOLVER_AVAILABLE
                _solver_klu.set_gridmodel(gridmodel);
                _solver_klu_single.set_gridmodel(gridmodel);
                _solver_klu_dc.set_gridmodel(gridmodel);
                _solver_fdpf_xb_klu.set_gridmodel(gridmodel);
                _solver_fdpf_bx_klu.set_gridmodel(gridmodel);
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
                _solver_nicslu.set_gridmodel(gridmodel);
                _solver_nicslu_single.set_gridmodel(gridmodel);
                _solver_nicslu_dc.set_gridmodel(gridmodel);
                _solver_fdpf_xb_nicslu.set_gridmodel(gridmodel);
                _solver_fdpf_bx_nicslu.set_gridmodel(gridmodel);
            #endif  // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
                _solver_cktso.set_gridmodel(gridmodel);
                _solver_cktso_single.set_gridmodel(gridmodel);
                _solver_cktso_dc.set_gridmodel(gridmodel);
                _solver_fdpf_xb_cktso.set_gridmodel(gridmodel);
                _solver_fdpf_bx_cktso.set_gridmodel(gridmodel);
            #endif  // CKTSO_SOLVER_AVAILABLE
        }

        bool ac_solver_used() const{
            auto p_solver = get_prt_solver("ac_solver_used", false);
            return p_solver->IS_AC;
        }

        void change_solver(const SolverType & type)
        {
            if(type == _solver_type) return;

            #ifndef KLU_SOLVER_AVAILABLE
                if((type == SolverType::KLU) || 
                   (type == SolverType::KLUDC) || 
                   (type == SolverType::KLUSingleSlack) ||
                   (type == SolverType::FDPF_XB_KLU) ||
                   (type == SolverType::FDPF_BX_KLU)
                   ){
                    std::string msg;
                    msg = "Impossible to change for a solver using KLU for linear algebra. Please compile lightsim2grid from source to benefit from this.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef NICSLU_SOLVER_AVAILABLE
                if((type == SolverType::NICSLU) || 
                   (type == SolverType::NICSLUDC) || 
                   (type ==  SolverType::NICSLUSingleSlack) || 
                   (type ==  SolverType::FDPF_XB_NICSLU) ||
                   (type ==  SolverType::FDPF_BX_NICSLU)
                   ){
                    std::string msg;
                    msg = "Impossible to change for a solver using NICSLU for linear algebra. Please compile lightsim2grid from source to benefit from this.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef CKTSO_SOLVER_AVAILABLE
                if((type == SolverType::CKTSO) || 
                   (type == SolverType::CKTSODC) || 
                   (type ==  SolverType::CKTSOSingleSlack) ||
                   (type ==  SolverType::FDPF_XB_CKTSO) ||
                   (type ==  SolverType::FDPF_BX_CKTSO)
                   ){
                    std::string msg;
                    msg = "Impossible to change for a solver using CKTSO for linear algebra. Please compile lightsim2grid from source to benefit from this.";
                    throw std::runtime_error(msg);
                }
            #endif

            // now switch the union (see https://en.cppreference.com/w/cpp/language/union)
            // reset the old solver
            reset();
            // and assign the right solver
            _solver_type = type;
            // and now reset the new one
            reset();
        }

        void reset()
        {
            auto p_solver = get_prt_solver("reset", false);  // i should not check if it's the right solver when resetting (used in change_solver)
            return p_solver -> reset();
        }

        // benefit from dynamic stuff and inheritance by having a method that returns a BaseAlgo *
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,  // size (nb_bus, nb_bus)
                        CplxVect & V,  // size nb_bus
                        const CplxVect & Sbus,  // size nb_bus
                        const Eigen::VectorXi & slack_ids,  // bus ids where thare are slack bus
                        const RealVect & slack_weights,  // slack weights (size nb_bus)
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        )
        {
            _type_used_for_nr = _solver_type;
            auto p_solver = get_prt_solver("compute_pf", true);
            return p_solver -> compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);
        }

        Eigen::Ref<const CplxVect> get_V() const
        {
            auto p_solver = get_prt_solver("get_V", true);
            return p_solver -> get_V();
        }

        Eigen::Ref<const RealVect> get_Va() const
        {
            auto p_solver = get_prt_solver( "get_Va", true);
            return p_solver -> get_Va();
        }
        Eigen::Ref<const RealVect> get_Vm() const
        {
            auto p_solver = get_prt_solver("get_Vm", true);
            return p_solver -> get_Vm();
        }

        Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J() const
        {
            check_right_solver("get_J");
            if(_solver_type == SolverType::SparseLU){
                return _solver_lu.get_J();}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                return _solver_klu.get_J();}
            else if(_solver_type == SolverType::KLUSingleSlack){
                return _solver_klu_single.get_J();}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                return _solver_nicslu.get_J();}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                return _solver_nicslu_single.get_J();}
            #endif // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::CKTSO){
                return _solver_cktso.get_J();}
            else if(_solver_type == SolverType::CKTSOSingleSlack){
                return _solver_cktso_single.get_J();}
            #endif // CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::DC || 
                    _solver_type == SolverType::KLUDC || 
                    _solver_type == SolverType::NICSLUDC ||
                    _solver_type == SolverType::CKTSODC){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the DC powerflow.");}
            else if(_solver_type == SolverType::FDPF_XB_SparseLU || 
                    _solver_type == SolverType::FDPF_BX_SparseLU || 
                    _solver_type == SolverType::FDPF_XB_KLU || 
                    _solver_type == SolverType::FDPF_BX_KLU || 
                    _solver_type == SolverType::FDPF_XB_NICSLU ||
                    _solver_type == SolverType::FDPF_BX_NICSLU ||
                    _solver_type == SolverType::FDPF_XB_CKTSO ||
                    _solver_type == SolverType::FDPF_BX_CKTSO){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the FDPF powerflow.");}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidelSynch powerflow.");}
            else if(_solver_type == SolverType::GaussSeidel){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidel powerflow.");}
            else throw std::runtime_error("Unknown solver type encountered (get_J)");
        }

        RealMat get_ptdf(const Eigen::SparseMatrix<cplx_type> & dcYbus){
                if(_solver_type != SolverType::DC && 
                   _solver_type != SolverType::KLUDC && 
                   _solver_type != SolverType::NICSLUDC &&
                   _solver_type != SolverType::CKTSODC){
                throw std::runtime_error("ChooseSolver::get_ptdf: cannot get ptdf for a solver that is not DC.");
                }
            auto p_solver = get_prt_solver("get_ptdf", true);
            const auto & res =  p_solver -> get_ptdf(dcYbus);
            return res;
        }

        void tell_solver_control(const SolverControl & solver_control){
            auto p_solver = get_prt_solver("tell_solver_control", false);
            p_solver -> tell_solver_control(solver_control);
        }

        /** apparently i cannot pass a const ref for a sparse matrix in python**/
        Eigen::SparseMatrix<real_type> get_J_python() const{
            Eigen::SparseMatrix<real_type> res = get_J();
            return res;
        }

        double get_computation_time() const
        {
            auto p_solver = get_prt_solver("get_computation_time", true);
            const auto & res =  p_solver -> get_timers();
            return std::get<3>(res);
        }

        std::tuple<double, double, double, double> get_timers() const
        {
            auto p_solver = get_prt_solver("get_timers", true);
            return p_solver -> get_timers();
        }

        TimerJacType get_timers_jacobian() const
        {
            const BaseAlgo * p_solver = get_prt_solver("get_timers_jacobian", true);
            return p_solver -> get_timers_jacobian();
        }

        ErrorType get_error() const{
            auto p_solver = get_prt_solver("get_error", true);
            return p_solver -> get_error();
        }
        
        int get_nb_iter() const {
            auto p_solver = get_prt_solver("get_nb_iter", true);
            return p_solver -> get_nb_iter();
        }

        bool converged() const{
            auto p_solver = get_prt_solver("converged", true);
            return p_solver -> converged();
        }

    private:
        void check_right_solver(const std::string & error_msg) const
        {
            if(_solver_type != _type_used_for_nr){
                std::ostringstream exc_;
                exc_ << "ChooseSolver: Solver mismatch when calling '";
                exc_ << error_msg;
                exc_ << ": current solver (";
                exc_ << _solver_type;
                exc_ << ") is not the one used to perform a powerflow (";
                exc_ << _type_used_for_nr;
                exc_ << ").";
                throw std::runtime_error(exc_.str());
            }
            
            #ifndef KLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::KLU){
                    std::string msg = "Impossible to use the KLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::KLUSingleSlack){
                    std::string msg = "Impossible to use the KLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::KLUDC){
                    std::string msg = "Impossible to use the KLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::FDPF_XB_KLU){
                    std::string msg = "Impossible to use the KLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                }else if(_solver_type == SolverType::FDPF_BX_KLU){
                    std::string msg = "Impossible to use the KLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                }
            #endif  // KLU_SOLVER_AVAILABLE

            #ifndef NICSLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::NICSLU){
                    std::string msg = "Impossible to use the NICSLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::NICSLUSingleSlack){
                    std::string msg = "Impossible to use the NICSLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::NICSLUDC){
                    std::string msg = "Impossible to use the NICSLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::FDPF_XB_NICSLU){
                    std::string msg = "Impossible to use the NICSLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::FDPF_BX_NICSLU){
                    std::string msg = "Impossible to use the NICSLU linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                }
            #endif  // NICSLU_SOLVER_AVAILABLE

            #ifndef CKTSO_SOLVER_AVAILABLE
                if(_solver_type == SolverType::CKTSO){
                    std::string msg = "Impossible to use the CKTSO linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::CKTSOSingleSlack){
                    std::string msg = "Impossible to use the CKTSO linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::CKTSODC){
                    std::string msg = "Impossible to use the CKTSO linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::FDPF_XB_CKTSO){
                    std::string msg = "Impossible to use the CKTSO linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::FDPF_BX_CKTSO){
                    std::string msg = "Impossible to use the CKTSO linear solver, your version of lightsim2grid has not been compiled to use it.";
                    throw std::runtime_error(msg);
                }
            #endif  // CKTSO_SOLVER_AVAILABLE
        }

    protected:
        /**
        returns a pointer to the current solver used
        **/
        const BaseAlgo * get_prt_solver(const std::string & error_msg, bool check_right_solver_=true) const {
            if (check_right_solver_) check_right_solver(error_msg);
            const BaseAlgo * res;
            if(_solver_type == SolverType::SparseLU){res = &_solver_lu;}
            else if(_solver_type == SolverType::SparseLUSingleSlack){res = &_solver_lu_single;}
            else if(_solver_type == SolverType::DC){res = &_solver_dc;}
            else if(_solver_type == SolverType::FDPF_XB_SparseLU){res = &_solver_fdpf_xb_lu;}
            else if(_solver_type == SolverType::FDPF_BX_SparseLU){res = &_solver_fdpf_bx_lu;}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){res = & _solver_klu;}
            else if(_solver_type == SolverType::KLUSingleSlack){res = &_solver_klu_single;}
            else if(_solver_type == SolverType::KLUDC){res = &_solver_klu_dc;}
            else if(_solver_type == SolverType::FDPF_XB_KLU){res = &_solver_fdpf_xb_klu;}
            else if(_solver_type == SolverType::FDPF_BX_KLU){res = &_solver_fdpf_bx_klu;}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){res = &_solver_nicslu;}
            else if(_solver_type == SolverType::NICSLUSingleSlack){res = &_solver_nicslu_single;}
            else if(_solver_type == SolverType::NICSLUDC){res = &_solver_nicslu_dc;}
            else if(_solver_type == SolverType::FDPF_XB_NICSLU){res = &_solver_fdpf_xb_nicslu;}
            else if(_solver_type == SolverType::FDPF_BX_NICSLU){res = &_solver_fdpf_bx_nicslu;}
            #endif // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::CKTSO){res = &_solver_cktso;}
            else if(_solver_type == SolverType::CKTSOSingleSlack){res = &_solver_cktso_single;}
            else if(_solver_type == SolverType::CKTSODC){res = &_solver_cktso_dc;}
            else if(_solver_type == SolverType::FDPF_XB_CKTSO){res = &_solver_fdpf_xb_cktso;}
            else if(_solver_type == SolverType::FDPF_BX_CKTSO){res = &_solver_fdpf_bx_cktso;}
            #endif // CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){res = &_solver_gaussseidel;}
            else if(_solver_type == SolverType::GaussSeidelSynch){res = &_solver_gaussseidelsynch;}
            else throw std::runtime_error("Unknown solver type encountered (ChooseSolver get_prt_solver const)");
            return res;
        }
        BaseAlgo * get_prt_solver(const std::string & error_msg, bool check_right_solver_=true) {
            if (check_right_solver_) check_right_solver(error_msg);
            BaseAlgo * res;
            if(_solver_type == SolverType::SparseLU){res = &_solver_lu;}
            else if(_solver_type == SolverType::SparseLUSingleSlack){res = &_solver_lu_single;}
            else if(_solver_type == SolverType::DC){res = &_solver_dc;}
            else if(_solver_type == SolverType::FDPF_XB_SparseLU){res = &_solver_fdpf_xb_lu;}
            else if(_solver_type == SolverType::FDPF_BX_SparseLU){res = &_solver_fdpf_bx_lu;}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){res = & _solver_klu;}
            else if(_solver_type == SolverType::KLUSingleSlack){res = &_solver_klu_single;}
            else if(_solver_type == SolverType::KLUDC){res = &_solver_klu_dc;}
            else if(_solver_type == SolverType::FDPF_XB_KLU){res = &_solver_fdpf_xb_klu;}
            else if(_solver_type == SolverType::FDPF_BX_KLU){res = &_solver_fdpf_bx_klu;}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){res = &_solver_nicslu;}
            else if(_solver_type == SolverType::NICSLUSingleSlack){res = &_solver_nicslu_single;}
            else if(_solver_type == SolverType::NICSLUDC){res = &_solver_nicslu_dc;}
            else if(_solver_type == SolverType::FDPF_XB_NICSLU){res = &_solver_fdpf_xb_nicslu;}
            else if(_solver_type == SolverType::FDPF_BX_NICSLU){res = &_solver_fdpf_bx_nicslu;}
            #endif // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::CKTSO){res = &_solver_cktso;}
            else if(_solver_type == SolverType::CKTSOSingleSlack){res = &_solver_cktso_single;}
            else if(_solver_type == SolverType::CKTSODC){res = &_solver_cktso_dc;}
            else if(_solver_type == SolverType::FDPF_XB_CKTSO){res = &_solver_fdpf_xb_cktso;}
            else if(_solver_type == SolverType::FDPF_BX_CKTSO){res = &_solver_fdpf_bx_cktso;}
            #endif // CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){res = &_solver_gaussseidel;}
            else if(_solver_type == SolverType::GaussSeidelSynch){res = &_solver_gaussseidelsynch;}
            else throw std::runtime_error("Unknown solver type encountered (ChooseSolver get_prt_solver non const)");
            return res;
        }
    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;

        // all types
        // TODO have a way to use Union here https://en.cppreference.com/w/cpp/language/union
        SparseLUSolver _solver_lu;
        SparseLUSolverSingleSlack _solver_lu_single;
        GaussSeidelAlgo _solver_gaussseidel;
        GaussSeidelSynchAlgo _solver_gaussseidelsynch;
        DCSolver _solver_dc;
        FDPF_XB_SparseLUSolver _solver_fdpf_xb_lu;
        FDPF_BX_SparseLUSolver _solver_fdpf_bx_lu;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
            KLUSolverSingleSlack _solver_klu_single;
            KLUDCSolver _solver_klu_dc;
            FDPF_XB_KLUSolver _solver_fdpf_xb_klu;
            FDPF_BX_KLUSolver _solver_fdpf_bx_klu;
        #endif  // KLU_SOLVER_AVAILABLE
        #ifdef NICSLU_SOLVER_AVAILABLE
            NICSLUSolver _solver_nicslu;
            NICSLUSolverSingleSlack _solver_nicslu_single;
            NICSLUDCSolver _solver_nicslu_dc;
            FDPF_XB_NICSLUSolver _solver_fdpf_xb_nicslu;
            FDPF_BX_NICSLUSolver _solver_fdpf_bx_nicslu;
        #endif  // NICSLU_SOLVER_AVAILABLE
        #ifdef CKTSO_SOLVER_AVAILABLE
            CKTSOSolver _solver_cktso;
            CKTSOSolverSingleSlack _solver_cktso_single;
            CKTSODCSolver _solver_cktso_dc;
            FDPF_XB_CKTSOSolver _solver_fdpf_xb_cktso;
            FDPF_BX_CKTSOSolver _solver_fdpf_bx_cktso;
        #endif  // CKTSO_SOLVER_AVAILABLE
};

#endif  //CHOOSESOLVER_H
