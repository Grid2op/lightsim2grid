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
#include "GaussSeidelSolver.h"
#include "GaussSeidelSynchSolver.h"
#include "DCSolver.h"

enum class SolverType {SparseLU, KLU, GaussSeidel, DC, GaussSeidelSynch, NICSLU, 
                       SparseLUSingleSlack, KLUSingleSlack, NICSLUSingleSlack, 
                       KLUDC, NICSLUDC,
                       CKTSO, CKTSOSingleSlack, CKTSODC};

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
            res.reserve(11);

            res.push_back(SolverType::SparseLU);
            res.push_back(SolverType::GaussSeidel);
            res.push_back(SolverType::DC);
            res.push_back(SolverType::GaussSeidelSynch);
            res.push_back(SolverType::SparseLUSingleSlack);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
                res.push_back(SolverType::KLUSingleSlack);
                res.push_back(SolverType::KLUDC);
            #endif
            #ifdef NICSLU_SOLVER_AVAILABLE
                res.push_back(SolverType::NICSLU);
                res.push_back(SolverType::NICSLUSingleSlack);
                res.push_back(SolverType::NICSLUDC);
            #endif
            #ifdef CKTSO_SOLVER_AVAILABLE
                res.push_back(SolverType::CKTSO);
                res.push_back(SolverType::CKTSOSingleSlack);
                res.push_back(SolverType::CKTSODC);
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
        
        bool ac_solver_used() const{
            auto p_solver = get_prt_solver("ac_solver_used", false);
            return p_solver->IS_AC;
        }

        void change_solver(const SolverType & type)
        {
            if(type == _solver_type) return;

            #ifndef KLU_SOLVER_AVAILABLE
                if((type == SolverType::KLU) || (type == SolverType::KLUDC) || (type == SolverType::KLUSingleSlack)){
                    std::string msg;
                    msg = "Impossible to change for the KLU solver, that is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef NICSLU_SOLVER_AVAILABLE
                if((type == SolverType::NICSLU) || (type == SolverType::NICSLUDC) || (type ==  SolverType::NICSLUSingleSlack)){
                    std::string msg;
                    msg = "Impossible to change for the NICSLU solver, that is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef CKTSO_SOLVER_AVAILABLE
                if((type == SolverType::CKTSO) || (type == SolverType::CKTSODC) || (type ==  SolverType::CKTSOSingleSlack)){
                    std::string msg;
                    msg = "Impossible to change for the CKTSO solver, that is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif

            // now switch the union (see https://en.cppreference.com/w/cpp/language/union)
            reset();

            // and assign the right solver
            _solver_type = type;
        }

        void reset()
        {
            auto p_solver = get_prt_solver("reset", false);  // i should not check if it's the right solver when resetting (used in change_solver)
            return p_solver -> reset();
        }

        // benefit from dynamic stuff and inheritance by having a method that returns a BaseSolver *
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
            check_right_solver( "get_J");
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
            else if(_solver_type == SolverType::GaussSeidel){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidel powerflow.");}
            else if(_solver_type == SolverType::DC){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the DC powerflow.");}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidelSynch powerflow.");}
            else throw std::runtime_error("Unknown solver type encountered");
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
            if(_solver_type != _type_used_for_nr) throw std::runtime_error("ChooseSolver: Solver mismatch when calling '"+error_msg+"': current solver is not the last solver used to perform a powerflow");
            
            #ifndef KLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::KLU){
                    std::string msg = "Impossible to use the KLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::KLUSingleSlack){
                    std::string msg = "Impossible to use the KLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::KLUDC){
                    std::string msg = "Impossible to use the KLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif  // KLU_SOLVER_AVAILABLE

            #ifndef NICSLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::NICSLU){
                    std::string msg = "Impossible to use the NICSLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::NICSLUSingleSlack){
                    std::string msg = "Impossible to use the NICSLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::NICSLUDC){
                    std::string msg = "Impossible to use the NICSLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif  // NICSLU_SOLVER_AVAILABLE

            #ifndef CKTSO_SOLVER_AVAILABLE
                if(_solver_type == SolverType::CKTSO){
                    std::string msg = "Impossible to use the CKTSO solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::CKTSOSingleSlack){
                    std::string msg = "Impossible to use the CKTSO solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                } else if(_solver_type == SolverType::CKTSODC){
                    std::string msg = "Impossible to use the CKTSO solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif  // CKTSO_SOLVER_AVAILABLE
        }

    protected:
        /**
        returns a pointer to the current solver used
        **/
        const BaseSolver * get_prt_solver(const std::string & error_msg, bool check_right_solver_=true) const {
            if (check_right_solver_) check_right_solver(error_msg);
            const BaseSolver * res;
            if(_solver_type == SolverType::SparseLU){res = &_solver_lu;}
            else if(_solver_type == SolverType::SparseLUSingleSlack){res = &_solver_lu_single;}
            else if(_solver_type == SolverType::DC){res = &_solver_dc;}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){res = & _solver_klu;}
            else if(_solver_type == SolverType::KLUSingleSlack){res = &_solver_klu_single;}
            else if(_solver_type == SolverType::KLUDC){res = &_solver_klu_dc;}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){res = &_solver_nicslu;}
            else if(_solver_type == SolverType::NICSLUSingleSlack){res = &_solver_nicslu_single;}
            else if(_solver_type == SolverType::NICSLUDC){res = &_solver_nicslu_dc;}
            #endif // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::CKTSO){res = &_solver_cktso;}
            else if(_solver_type == SolverType::CKTSOSingleSlack){res = &_solver_cktso_single;}
            else if(_solver_type == SolverType::CKTSODC){res = &_solver_cktso_dc;}
            #endif // CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){res = &_solver_gaussseidel;}
            else if(_solver_type == SolverType::GaussSeidelSynch){res = &_solver_gaussseidelsynch;}
            else throw std::runtime_error("Unknown solver type encountered");
            return res;
        }
        BaseSolver * get_prt_solver(const std::string & error_msg, bool check_right_solver_=true) {
            if (check_right_solver_) check_right_solver(error_msg);
            BaseSolver * res;
            if(_solver_type == SolverType::SparseLU){res = &_solver_lu;}
            else if(_solver_type == SolverType::SparseLUSingleSlack){res = &_solver_lu_single;}
            else if(_solver_type == SolverType::DC){res = &_solver_dc;}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){res = & _solver_klu;}
            else if(_solver_type == SolverType::KLUSingleSlack){res = &_solver_klu_single;}
            else if(_solver_type == SolverType::KLUDC){res = &_solver_klu_dc;}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){res = &_solver_nicslu;}
            else if(_solver_type == SolverType::NICSLUSingleSlack){res = &_solver_nicslu_single;}
            else if(_solver_type == SolverType::NICSLUDC){res = &_solver_nicslu_dc;}
            #endif // NICSLU_SOLVER_AVAILABLE
            #ifdef CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::CKTSO){res = &_solver_cktso;}
            else if(_solver_type == SolverType::CKTSOSingleSlack){res = &_solver_cktso_single;}
            else if(_solver_type == SolverType::CKTSODC){res = &_solver_cktso_dc;}
            #endif // CKTSO_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){res = &_solver_gaussseidel;}
            else if(_solver_type == SolverType::GaussSeidelSynch){res = &_solver_gaussseidelsynch;}
            else throw std::runtime_error("Unknown solver type encountered");
            return res;
        }
    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;

        // all types
        // TODO have a way to use Union here https://en.cppreference.com/w/cpp/language/union
        SparseLUSolver _solver_lu;
        SparseLUSolverSingleSlack _solver_lu_single;
        GaussSeidelSolver _solver_gaussseidel;
        GaussSeidelSynchSolver _solver_gaussseidelsynch;
        DCSolver _solver_dc;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
            KLUSolverSingleSlack _solver_klu_single;
            KLUDCSolver _solver_klu_dc;
        #endif  // KLU_SOLVER_AVAILABLE
        #ifdef NICSLU_SOLVER_AVAILABLE
            NICSLUSolver _solver_nicslu;
            NICSLUSolverSingleSlack _solver_nicslu_single;
            NICSLUDCSolver _solver_nicslu_dc;
        #endif  // NICSLU_SOLVER_AVAILABLE
        #ifdef CKTSO_SOLVER_AVAILABLE
            CKTSOSolver _solver_cktso;
            CKTSOSolverSingleSlack _solver_cktso_single;
            CKTSODCSolver _solver_cktso_dc;
        #endif  // CKTSO_SOLVER_AVAILABLE

};

#endif  //CHOOSESOLVER_H
