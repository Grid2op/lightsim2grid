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

enum class SolverType {SparseLU, KLU, GaussSeidel, DC, GaussSeidelSynch, NICSLU};

// TODO define a template class instead of these weird stuff !!!


// NB: when adding a new solver, you need to specialize the *tmp method (eg get_Va_tmp)
// and also to "forward" the specialisation (adding the if(solvertype==XXX)) in the compute_pf, get_V, get_J, get_Va, get_Vm
// and the "available_solvers" (add it to the list)
// and to add a attribute with the proper class and the reset method
class ChooseSolver
{
    public:
         ChooseSolver():_solver_type(SolverType::SparseLU),_type_used_for_nr(SolverType::SparseLU){};

        std::vector<SolverType> available_solvers()
        {
            std::vector<SolverType> res;
            res.push_back(SolverType::SparseLU);
            res.push_back(SolverType::GaussSeidel);
            res.push_back(SolverType::DC);
            res.push_back(SolverType::GaussSeidelSynch);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
            #endif
            #ifdef NICSLU_SOLVER_AVAILABLE
                res.push_back(SolverType::NICSLU);
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
                msg = "Impossible to change for the KLU solver, that is not available on your platform.";
                if(type == SolverType::KLU) throw std::runtime_error(msg);
            #endif
            #ifndef NICSLU_SOLVER_AVAILABLE
                // TODO better handling of that :-/
                msg = "Impossible to change for the NICSLU solver, that is not available on your platform.";
                if(type == SolverType::NICSLU) throw std::runtime_error(msg);
            #endif
            _solver_type = type;
        }
        void reset()
        {
            // reset all the solvers available
            _solver_lu.reset();
            _solver_gaussseidel.reset();
            _solver_dc.reset();
            #ifdef KLU_SOLVER_AVAILABLE
                _solver_klu.reset();
            #endif  // KLU_SOLVER_AVAILABLE
            #ifdef NICSLU_SOLVER_AVAILABLE
                _solver_nicslu.reset();
            #endif  // NICSLU_SOLVER_AVAILABLE
            _solver_gaussseidelsynch.reset();
        }

        // forward to the right solver used
        //TODO inline all of that
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,  // size (nb_bus, nb_bus)
                        CplxVect & V,  // size nb_bus
                        const CplxVect & Sbus,  // size nb_bus
                        const Eigen::VectorXi & slack_ids,  // bus ids where thare are slack bus
                        const RealVect & slack_weights,  // slack weights (size nb_bus)
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        );
        Eigen::Ref<const CplxVect> get_V() const;
        Eigen::Ref<const RealVect> get_Va() const;
        Eigen::Ref<const RealVect> get_Vm() const;
        Eigen::SparseMatrix<real_type> get_J();
        double get_computation_time();

    private:
        void check_right_solver() const
        {
            if(_solver_type != _type_used_for_nr) throw std::runtime_error("ChooseSolver: Solver mismatch: current solver is not the last solver used to perform a powerflow");
        }

        template<SolverType ST>
        Eigen::Ref<const CplxVect> get_V_tmp() const;

        template<SolverType ST>
        Eigen::Ref<const RealVect> get_Va_tmp() const;

        template<SolverType ST>
        Eigen::Ref<const RealVect> get_Vm_tmp() const;

        template<SolverType ST>
        Eigen::SparseMatrix<real_type> get_J_tmp();

        template<SolverType ST>
        double get_computation_time_tmp();

        template<SolverType ST>
        bool compute_pf_tmp(const Eigen::SparseMatrix<cplx_type> & Ybus,
                            CplxVect & V,
                            const CplxVect & Sbus,
                            const Eigen::VectorXi & slack_ids,
                            const RealVect & slack_weights,
                            const Eigen::VectorXi & pv,
                            const Eigen::VectorXi & pq,
                            int max_iter,
                            real_type tol
                            );

    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;

        // all types
        SparseLUSolver _solver_lu;
        GaussSeidelSolver _solver_gaussseidel;
        GaussSeidelSynchSolver _solver_gaussseidelsynch;
        DCSolver _solver_dc;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
        #endif  // KLU_SOLVER_AVAILABLE
        #ifdef NICSLU_SOLVER_AVAILABLE
            NICSLUSolver _solver_nicslu;
        #endif  // NICSLU_SOLVER_AVAILABLE

};


// template specialization

template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::SparseLU>() const;
template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::KLU>() const;
template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::NICSLU>() const;
template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::GaussSeidel>() const;
template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::GaussSeidelSynch>() const;
template<>
Eigen::Ref<const CplxVect> ChooseSolver::get_V_tmp<SolverType::DC>() const;

template<>
bool ChooseSolver::compute_pf_tmp<SolverType::SparseLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::KLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::NICSLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::GaussSeidel>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::GaussSeidelSynch>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::DC>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & slack_ids,
                       const RealVect & slack_weights,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       );

template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::SparseLU>();
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::KLU>();
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::NICSLU>();
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::GaussSeidel>();
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::GaussSeidelSynch>();
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::DC>();

template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::SparseLU>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::KLU>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::NICSLU>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::GaussSeidel>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::GaussSeidelSynch>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Va_tmp<SolverType::DC>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Vm_tmp<SolverType::SparseLU>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Vm_tmp<SolverType::KLU>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidel>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidelSynch>() const;
template<>
Eigen::Ref<const RealVect> ChooseSolver::get_Vm_tmp<SolverType::DC>() const;

// computation times
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::SparseLU>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::KLU>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::NICSLU>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::GaussSeidel>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::GaussSeidelSynch>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::DC>();

#endif  //CHOOSESOLVER_H
