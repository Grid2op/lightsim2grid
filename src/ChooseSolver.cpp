// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ChooseSolver.h"
#include <iostream>
// template specialization
template<SolverType ST>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp()
{
    throw std::runtime_error("ChooseSolver::get_V_tmp: Unknown solver type.");
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_V();
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_V_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_V();
    #endif
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::NICSLU>()
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_V_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        return _solver_nicslu.get_V();
    #endif
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_V();
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::GaussSeidelSynch>()
{
    return _solver_gaussseidelsynch.get_V();
}
template<>
Eigen::Ref<CplxVect> ChooseSolver::get_V_tmp<SolverType::DC>()
{
    return _solver_dc.get_V();
}


template<SolverType ST>
bool ChooseSolver::compute_pf_tmp(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                  CplxVect & V,
                                  const CplxVect & Sbus,
                                  const Eigen::VectorXi & pv,
                                  const Eigen::VectorXi & pq,
                                  int max_iter,
                                  real_type tol
                                  )
{
    throw std::runtime_error("ChooseSolver::compute_pf_tmp: Unknown solver type.");
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::SparseLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    return _solver_lu.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::GaussSeidel>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    return _solver_gaussseidel.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::GaussSeidelSynch>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    return _solver_gaussseidelsynch.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::DC>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    return _solver_dc.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::KLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::compute_pf_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
    #endif
}
template<>
bool ChooseSolver::compute_pf_tmp<SolverType::NICSLU>(const Eigen::SparseMatrix<cplx_type> & Ybus,
                       CplxVect & V,
                       const CplxVect & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       real_type tol
                       )
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::compute_pf_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        return _solver_nicslu.compute_pf(Ybus, V, Sbus, pv, pq, max_iter, tol);
    #endif
}

template<SolverType ST>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp()
{
    throw std::runtime_error("ChooseSolver::get_J_tmp: Unknown solver type.");
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_J();
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::DC>()
{
    throw std::runtime_error("ChooseSolver::get_J_tmp: There is not Jacobian matrix for a DC powerflow.");
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::GaussSeidel>()
{
    throw std::runtime_error("ChooseSolver::get_J_tmp: There is not Jacobian matrix for the GaussSeidel powerflow.");
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::GaussSeidelSynch>()
{
    throw std::runtime_error("hooseSolver::get_J_tmp: There is not Jacobian matrix for the GaussSeidelSynch powerflow.");
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_J_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_J();
    #endif
}
template<>
Eigen::SparseMatrix<real_type> ChooseSolver::get_J_tmp<SolverType::NICSLU>()
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_J_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        return _solver_nicslu.get_J();
    #endif
}

template<SolverType ST>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp()
{
    throw std::runtime_error("ChooseSolver::get_Va_tmp: Unknown solver type.");
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_Va();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_Va();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::GaussSeidelSynch>()
{
    return _solver_gaussseidelsynch.get_Va();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::DC>()
{
    return _solver_dc.get_Va();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_Va_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_Va();
    #endif
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Va_tmp<SolverType::NICSLU>()
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_Va_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        return _solver_nicslu.get_Va();
    #endif
}
template<SolverType ST>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp()
{
    throw std::runtime_error("ChooseSolver::get_Vm_tmp: Unknown solver type.");
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_Vm();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::DC>()
{
    return _solver_dc.get_Vm();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_Vm();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidelSynch>()
{
    return _solver_gaussseidelsynch.get_Vm();
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_Vm_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_Vm();
    #endif
}
template<>
Eigen::Ref<RealVect> ChooseSolver::get_Vm_tmp<SolverType::NICSLU>()
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_Vm_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        return _solver_nicslu.get_Vm();
    #endif
}

template<SolverType ST>
double ChooseSolver::get_computation_time_tmp()
{
    throw std::runtime_error("ChooseSolver::get_computation_time_tmp: Unknown solver type.");
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::SparseLU>()
{
    const auto & res =  _solver_lu.get_timers();
    return std::get<3>(res);
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::GaussSeidel>()
{
    const auto & res =  _solver_gaussseidel.get_timers();
    return std::get<3>(res);
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::GaussSeidelSynch>()
{
    const auto & res =  _solver_gaussseidelsynch.get_timers();
    return std::get<3>(res);
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_computation_time_tmp: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        const auto & res =  _solver_klu.get_timers();
        return std::get<3>(res);
    #endif
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::NICSLU>()
{
    #ifndef NICSLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("ChooseSolver::get_computation_time_tmp: Impossible to use the NICSLU solver, that is not available on your plaform.");
    #else
        const auto & res =  _solver_nicslu.get_timers();
        return std::get<3>(res);
    #endif
}
template<>
double ChooseSolver::get_computation_time_tmp<SolverType::DC>()
{
   const auto & res =  _solver_dc.get_timers();
   return std::get<3>(res);
}
//TODO refactor all the functions above by making a template function "get_solver"

// function definition
Eigen::Ref<CplxVect> ChooseSolver::get_V(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_V_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_V_tmp<SolverType::KLU>();
    }else if(_solver_type == SolverType::NICSLU){
         return get_V_tmp<SolverType::NICSLU>();
    }else if(_solver_type == SolverType::GaussSeidel){
         return get_V_tmp<SolverType::GaussSeidel>();
    }else if(_solver_type == SolverType::GaussSeidelSynch){
         return get_V_tmp<SolverType::GaussSeidelSynch>();
    }else if(_solver_type == SolverType::DC){
         return get_V_tmp<SolverType::DC>();
    }else{
        throw std::runtime_error("ChooseSolver::get_V: Unknown solver type.");
    }
}

Eigen::Ref<RealVect> ChooseSolver::get_Va(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_Va_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_Va_tmp<SolverType::KLU>();
    }else if(_solver_type == SolverType::NICSLU){
         return get_Va_tmp<SolverType::NICSLU>();
    }else if(_solver_type == SolverType::GaussSeidel){
         return get_Va_tmp<SolverType::GaussSeidel>();
    }else if(_solver_type == SolverType::GaussSeidelSynch){
         return get_Va_tmp<SolverType::GaussSeidelSynch>();
    }else if(_solver_type == SolverType::DC){
         return get_Va_tmp<SolverType::DC>();
    }else{
        throw std::runtime_error("ChooseSolver::get_Va: Unknown solver type.");
    }
}
Eigen::Ref<RealVect> ChooseSolver::get_Vm(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_Vm_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_Vm_tmp<SolverType::KLU>();
    }else if(_solver_type == SolverType::NICSLU){
         return get_Vm_tmp<SolverType::NICSLU>();
    }else if(_solver_type == SolverType::GaussSeidel){
         return get_Vm_tmp<SolverType::GaussSeidel>();
    }else if(_solver_type == SolverType::GaussSeidelSynch){
         return get_Vm_tmp<SolverType::GaussSeidelSynch>();
    }else if(_solver_type == SolverType::DC){
         return get_Vm_tmp<SolverType::DC>();
    }else{
        throw std::runtime_error("ChooseSolver::get_Vm: Unknown solver type.");
    }
}

bool ChooseSolver::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                              CplxVect & V,
                              const CplxVect & Sbus,
                              const Eigen::VectorXi & pv,
                              const Eigen::VectorXi & pq,
                              int max_iter,
                              real_type tol
                              )
{
    _type_used_for_nr = _solver_type;
    if(_solver_type == SolverType::SparseLU)
    {
        return compute_pf_tmp<SolverType::SparseLU>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::KLU){
        return compute_pf_tmp<SolverType::KLU>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::NICSLU){
        return compute_pf_tmp<SolverType::NICSLU>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::GaussSeidel){
        return compute_pf_tmp<SolverType::GaussSeidel>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::GaussSeidelSynch){
        return compute_pf_tmp<SolverType::GaussSeidelSynch>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::DC){
        return compute_pf_tmp<SolverType::DC>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else{
        throw std::runtime_error("ChooseSolver::compute_pf: Unknown solver type.");
    }
}

Eigen::SparseMatrix<real_type> ChooseSolver::get_J(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_J_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_J_tmp<SolverType::KLU>();
    }else if(_solver_type == SolverType::NICSLU){
         return get_J_tmp<SolverType::NICSLU>();
    }else if(_solver_type == SolverType::GaussSeidel){
         return get_J_tmp<SolverType::GaussSeidel>();
    }else if(_solver_type == SolverType::GaussSeidelSynch){
         return get_J_tmp<SolverType::GaussSeidelSynch>();
    }else if(_solver_type == SolverType::DC){
         return get_J_tmp<SolverType::DC>();
    }else{
        throw std::runtime_error("ChooseSolver::get_J: Unknown solver type.");
    }
}

double ChooseSolver::get_computation_time()
{
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_computation_time_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_computation_time_tmp<SolverType::KLU>();
    }else if(_solver_type == SolverType::NICSLU){
         return get_computation_time_tmp<SolverType::NICSLU>();
    }else if(_solver_type == SolverType::GaussSeidel){
         return get_computation_time_tmp<SolverType::GaussSeidel>();
    }else if(_solver_type == SolverType::GaussSeidelSynch){
         return get_computation_time_tmp<SolverType::GaussSeidelSynch>();
    }else if(_solver_type == SolverType::DC){
         return get_computation_time_tmp<SolverType::DC>();
    }else{
        throw std::runtime_error("ChooseSolver::get_computation_time: Unknown solver type.");
    }
}
