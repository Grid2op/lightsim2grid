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
                       SparseLUSingleSlack, KLUSingleSlack, NICSLUSingleSlack};

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

        std::vector<SolverType> available_solvers()
        {
            std::vector<SolverType> res;
            res.push_back(SolverType::SparseLU);
            res.push_back(SolverType::GaussSeidel);
            res.push_back(SolverType::DC);
            res.push_back(SolverType::GaussSeidelSynch);
            res.push_back(SolverType::SparseLUSingleSlack);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
                res.push_back(SolverType::KLUSingleSlack);
            #endif
            #ifdef NICSLU_SOLVER_AVAILABLE
                res.push_back(SolverType::NICSLU);
                res.push_back(SolverType::NICSLUSingleSlack);
            #endif
            return res;
        }
        SolverType get_type() const {return _solver_type;}
        void change_solver(const SolverType & type)
        {
            if(type == _solver_type) return;
            std::string msg;
            #ifndef KLU_SOLVER_AVAILABLE
                // TODO better handling of that :-/
                if(type == SolverType::KLU){
                    msg = "Impossible to change for the KLU solver, that is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef NICSLU_SOLVER_AVAILABLE
                // TODO better handling of that :-/
                if(type == SolverType::NICSLU){
                    msg = "Impossible to change for the NICSLU solver, that is not available on your platform.";
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
            
            if(_solver_type == SolverType::SparseLU){_solver_lu.reset();}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){_solver_klu.reset();}
            else if(_solver_type == SolverType::KLUSingleSlack){_solver_klu_single.reset();}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){_solver_nicslu.reset();}
            else if(_solver_type == SolverType::NICSLUSingleSlack){_solver_nicslu_single.reset();}
            #endif // NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){_solver_gaussseidel.reset();}
            else if(_solver_type == SolverType::DC){_solver_dc.reset();}
            else if(_solver_type == SolverType::GaussSeidelSynch){_solver_gaussseidelsynch.reset();}
            else throw std::runtime_error("Unknown solver type encountered");
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
            check_right_solver();
            if(_solver_type == SolverType::SparseLU){
                return _solver_lu.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                return _solver_klu.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            else if(_solver_type == SolverType::KLUSingleSlack){
                return _solver_klu_single.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                return _solver_nicslu.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                return _solver_nicslu_single.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            #endif // NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){
                return _solver_gaussseidel.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            else if(_solver_type == SolverType::DC){
                return _solver_dc.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                return _solver_gaussseidelsynch.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);}
            else throw std::runtime_error("Unknown solver type encountered");
        }

        Eigen::Ref<const CplxVect> get_V() const
        {
            check_right_solver();
            if(_solver_type == SolverType::SparseLU){
                return _solver_lu.get_V();}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                return _solver_klu.get_V();}
            else if(_solver_type == SolverType::KLUSingleSlack){
                return _solver_klu_single.get_V();}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                return _solver_nicslu.get_V();}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                return _solver_nicslu_single.get_V();}
            #endif // NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){
                return _solver_gaussseidel.get_V();}
            else if(_solver_type == SolverType::DC){
                return _solver_dc.get_V();}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                return _solver_gaussseidelsynch.get_V();}
            else throw std::runtime_error("Unknown solver type encountered");
        }

        Eigen::Ref<const RealVect> get_Va() const
        {
            check_right_solver();
            if(_solver_type == SolverType::SparseLU){
                return _solver_lu.get_Va();}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                return _solver_klu.get_Va();}
            else if(_solver_type == SolverType::KLUSingleSlack){
                return _solver_klu_single.get_Va();}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                return _solver_nicslu.get_Va();}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                return _solver_nicslu_single.get_Va();}
            #endif // NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){
                return _solver_gaussseidel.get_Va();}
            else if(_solver_type == SolverType::DC){
                return _solver_dc.get_Va();}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                return _solver_gaussseidelsynch.get_Va();}
            else throw std::runtime_error("Unknown solver type encountered");
        }
        Eigen::Ref<const RealVect> get_Vm() const
        {
            check_right_solver();
            if(_solver_type == SolverType::SparseLU){
                return _solver_lu.get_Vm();}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                return _solver_klu.get_Vm();}
            else if(_solver_type == SolverType::KLUSingleSlack){
                return _solver_klu_single.get_Vm();}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                return _solver_nicslu.get_Vm();}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                return _solver_nicslu_single.get_Vm();}
            #endif // NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::GaussSeidel){
                return _solver_gaussseidel.get_Vm();}
            else if(_solver_type == SolverType::DC){
                return _solver_dc.get_Vm();}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                return _solver_gaussseidelsynch.get_Vm();}
            else throw std::runtime_error("Unknown solver type encountered");
        }
        Eigen::SparseMatrix<real_type> get_J()
        {
            check_right_solver();
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
            else if(_solver_type == SolverType::GaussSeidel){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidel powerflow.");}
            else if(_solver_type == SolverType::DC){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the DC powerflow.");}
            else if(_solver_type == SolverType::GaussSeidelSynch){
                throw std::runtime_error("ChooseSolver::get_J: There is not Jacobian matrix for the GaussSeidelSynch powerflow.");}
            else throw std::runtime_error("Unknown solver type encountered");
        }

        double get_computation_time()
        {
            check_right_solver();
            if(_solver_type == SolverType::SparseLU){
                const auto & res =  _solver_lu.get_timers();
                return std::get<3>(res);}
            #ifdef KLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::KLU){
                const auto & res =  _solver_klu.get_timers();
                return std::get<3>(res);}
            else if(_solver_type == SolverType::KLUSingleSlack){
                const auto & res =  _solver_klu_single.get_timers();
                return std::get<3>(res);}
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
            else if(_solver_type == SolverType::NICSLU){
                const auto & res =  _solver_nicslu.get_timers();
                return std::get<3>(res);}
            else if(_solver_type == SolverType::NICSLUSingleSlack){
                const auto & res =  _solver_nicslu_single.get_timers();
                return std::get<3>(res);}
            #endif // NICSLU_SOLVER_AVAILABLE
             else if(_solver_type == SolverType::GaussSeidel){
                const auto & res =  _solver_gaussseidel.get_timers();
                return std::get<3>(res);}
             else if(_solver_type == SolverType::DC){
                const auto & res =  _solver_dc.get_timers();
                return std::get<3>(res);}
             else if(_solver_type == SolverType::GaussSeidelSynch){
                const auto & res = _solver_gaussseidelsynch.get_timers();
                return std::get<3>(res);}
            else throw std::runtime_error("Unknown solver type encountered");
        }

    private:
        void check_right_solver() const
        {
            if(_solver_type != _type_used_for_nr) throw std::runtime_error("ChooseSolver: Solver mismatch: current solver is not the last solver used to perform a powerflow");
            
            #ifndef KLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::KLU){
                    std::string msg = "Impossible to use the KLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
                else if(_solver_type == SolverType::KLUSingleSlack){
                    std::string msg = "Impossible to use the KLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif

            #ifndef NICSLU_SOLVER_AVAILABLE
                if(_solver_type == SolverType::NICSLU){
                    std::string msg = "Impossible to use the NICSLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
                if(_solver_type == SolverType::NICSLUSingleSlack){
                    std::string msg = "Impossible to use the NICSLU solver, it is not available on your platform.";
                    throw std::runtime_error(msg);
                }
            #endif
        }

    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;

        // all types
        // TODO have a way to use Union here https://en.cppreference.com/w/cpp/language/union
        SparseLUSolver _solver_lu;
        GaussSeidelSolver _solver_gaussseidel;
        GaussSeidelSynchSolver _solver_gaussseidelsynch;
        DCSolver _solver_dc;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
            KLUSolverSingleSlack _solver_klu_single;
        #endif  // KLU_SOLVER_AVAILABLE
        #ifdef NICSLU_SOLVER_AVAILABLE
            NICSLUSolver _solver_nicslu;
            NICSLUSolverSingleSlack _solver_nicslu_single;
        #endif  // NICSLU_SOLVER_AVAILABLE

};

#endif  //CHOOSESOLVER_H
