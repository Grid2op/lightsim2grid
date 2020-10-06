#include "ChooseSolver.h"
#include <iostream>
// template specialization
template<SolverType ST>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp()
{
    throw std::runtime_error("Unknown solver type.");
}
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_V();
}
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("get_V: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_V();
    #endif
}
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_V();
}


template<SolverType ST>
bool ChooseSolver::do_newton_tmp(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       )
{
    throw std::runtime_error("Unknown solver type.");
}
template<>
bool ChooseSolver::do_newton_tmp<SolverType::SparseLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       )
{
    return _solver_lu.do_newton(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::do_newton_tmp<SolverType::GaussSeidel>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       )
{
    return _solver_gaussseidel.do_newton(Ybus, V, Sbus, pv, pq, max_iter, tol);
}
template<>
bool ChooseSolver::do_newton_tmp<SolverType::KLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       )
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("do_newton: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.do_newton(Ybus, V, Sbus, pv, pq, max_iter, tol);
    #endif
}

template<SolverType ST>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp()
{
    throw std::runtime_error("Unknown solver type.");
}
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_J();
}
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_J();
}
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("get_J: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_J();
    #endif
}

template<SolverType ST>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp()
{
    throw std::runtime_error("Unknown solver type.");
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_Va();
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_Va();
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("get_Va: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_Va();
    #endif
}
template<SolverType ST>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp()
{
    throw std::runtime_error("Unknown solver type.");
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::SparseLU>()
{
    return _solver_lu.get_Vm();
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::GaussSeidel>()
{
    return _solver_gaussseidel.get_Vm();
}
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::KLU>()
{
    #ifndef KLU_SOLVER_AVAILABLE
        // I asked result of KLU solver without the required libraries
        throw std::runtime_error("get_Va: Impossible to use the KLU solver, that is not available on your plaform.");
    #else
        return _solver_klu.get_Vm();
    #endif
}

// function definition
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_V_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_V_tmp<SolverType::KLU>();
    }else{
        throw std::runtime_error("Unknown solver type.");
    }
}

Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_Va_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_Va_tmp<SolverType::KLU>();
    }else{
        throw std::runtime_error("Unknown solver type.");
    }
}
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_Vm_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_Vm_tmp<SolverType::KLU>();
    }else{
        throw std::runtime_error("Unknown solver type.");
    }
}

bool ChooseSolver::do_newton(const Eigen::SparseMatrix<cdouble> & Ybus,
                             Eigen::VectorXcd & V,
                             const Eigen::VectorXcd & Sbus,
                             const Eigen::VectorXi & pv,
                             const Eigen::VectorXi & pq,
                             int max_iter,
                             double tol
                             )
{
    _type_used_for_nr = _solver_type;
    if(_solver_type == SolverType::SparseLU)
    {
        return do_newton_tmp<SolverType::SparseLU>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else if(_solver_type == SolverType::KLU){
        return do_newton_tmp<SolverType::KLU>(Ybus, V, Sbus, pv, pq, max_iter, tol);
    }else{
        throw std::runtime_error("Unknown solver type.");
    }
}

Eigen::SparseMatrix<double> ChooseSolver::get_J(){
    check_right_solver();
    if(_solver_type == SolverType::SparseLU)
    {
         return get_J_tmp<SolverType::SparseLU>();
    }else if(_solver_type == SolverType::KLU){
         return get_J_tmp<SolverType::KLU>();
    }else{
        throw std::runtime_error("Unknown solver type.");
    }
}
