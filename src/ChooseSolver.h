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
#include "KLUSolver.h"
#include "SparseLUSolver.h"
#include "GaussSeidelSolver.h"
#include "DCSolver.h"

enum class SolverType { SparseLU, KLU, GaussSeidel, DC};


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
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
            #endif
            return res;
        }
        SolverType get_type() const {return _solver_type;}
        void change_solver(const SolverType & type)
        {
            if(type == _solver_type) return;
            #ifndef KLU_SOLVER_AVAILABLE
                // TODO better handling of that :-/
                if(type == SolverType::KLU) throw std::runtime_error("Impossible to change for the KLU solver, that is not available on your platform.");
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
        }

        // forward to the right solver used
        //TODO inline all of that
        bool compute_pf(const Eigen::SparseMatrix<cdouble> & Ybus,
                        Eigen::VectorXcd & V,
                        const Eigen::VectorXcd & Sbus,
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        double tol
                        );
        Eigen::Ref<Eigen::VectorXcd> get_V();
        Eigen::SparseMatrix<double> get_J();
        Eigen::Ref<Eigen::VectorXd> get_Va();
        Eigen::Ref<Eigen::VectorXd> get_Vm();
        double get_computation_time();

    private:
        void check_right_solver()
        {
            if(_solver_type != _type_used_for_nr) throw std::runtime_error("Solver mismatch between the performing of the newton raphson and the retrieval of the result.");
        }

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXcd> get_V_tmp();

        template<SolverType ST>
        Eigen::SparseMatrix<double> get_J_tmp();

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXd> get_Va_tmp();

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXd> get_Vm_tmp();

        template<SolverType ST>
        double get_computation_time_tmp();

        template<SolverType ST>
        bool compute_pf_tmp(const Eigen::SparseMatrix<cdouble> & Ybus,
                            Eigen::VectorXcd & V,
                            const Eigen::VectorXcd & Sbus,
                            const Eigen::VectorXi & pv,
                            const Eigen::VectorXi & pq,
                            int max_iter,
                            double tol
                            );

    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;

        // all types
        SparseLUSolver _solver_lu;
        GaussSeidelSolver _solver_gaussseidel;
        DCSolver _solver_dc;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
        #endif  // KLU_SOLVER_AVAILABLE

};


// template specialization

template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::KLU>();
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::GaussSeidel>();
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::DC>();

template<>
bool ChooseSolver::compute_pf_tmp<SolverType::SparseLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::KLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::GaussSeidel>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::DC>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );

template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::SparseLU>();
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::KLU>();
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::GaussSeidel>();
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::DC>();

template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::KLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::GaussSeidel>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::DC>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::KLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidel>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::DC>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::KLU>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::GaussSeidel>();
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::DC>();

#endif  //CHOOSESOLVER_H
