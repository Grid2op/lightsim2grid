// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CHOOSESOLVER_H
#define CHOOSESOLVER_H

#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Concrete solver typedefs (SparseLUSolver, KLUSolver, …) + Solvers.hpp guards
#include "Solvers.hpp"
// Registry for solver creation by name
#include "SolverRegistry.hpp"
#include "AlgoConfig.hpp"

namespace ls2g {

enum class LS2G_API SolverType {SparseLU, KLU, GaussSeidel, DC, GaussSeidelSynch, NICSLU,
                       SparseLUSingleSlack, KLUSingleSlack, NICSLUSingleSlack,
                       KLUDC, NICSLUDC,
                       CKTSO, CKTSOSingleSlack, CKTSODC,
                       FDPF_XB_SparseLU, FDPF_BX_SparseLU, // from 0.7.5
                       FDPF_XB_KLU, FDPF_BX_KLU,           // from 0.7.5
                       FDPF_XB_NICSLU, FDPF_BX_NICSLU,     // from 0.7.5
                       FDPF_XB_CKTSO, FDPF_BX_CKTSO,       // from 0.7.5
                       Custom                               // external/plugin solvers
                       };

std::ostream& operator<<(std::ostream& out, const SolverType& solver_type);

class LS2G_API ChooseSolver final
{
    public:
        // Default-constructs a SparseLU solver via the registry.
        // The registry must be populated (call register_builtin_solvers) before
        // any ChooseSolver instance is created — the PYBIND11_MODULE init does this.
        ChooseSolver();

        ~ChooseSolver() noexcept = default;

        ChooseSolver(const ChooseSolver&) = delete;
        ChooseSolver(ChooseSolver&&) = delete;
        ChooseSolver& operator=(ChooseSolver&&) = delete;
        ChooseSolver& operator=(const ChooseSolver&) = delete;

        // Returns the enum-typed solvers available in this build.
        // Does not include plugin (Custom) solvers; use SolverRegistry::available_solvers()
        // for the full list.
        std::vector<SolverType> available_solvers() const
        {
            std::vector<SolverType> res;
            res.reserve(22);
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

        bool is_dc(const SolverType& type) const noexcept {
            return (type == SolverType::DC) ||
                   (type == SolverType::KLUDC) ||
                   (type == SolverType::NICSLUDC) ||
                   (type == SolverType::CKTSODC);
        }
        bool is_fdpf(const SolverType& type) const noexcept {
            return (type == SolverType::FDPF_XB_SparseLU) ||
                   (type == SolverType::FDPF_BX_SparseLU) ||
                   (type == SolverType::FDPF_XB_KLU) ||
                   (type == SolverType::FDPF_BX_KLU) ||
                   (type == SolverType::FDPF_XB_NICSLU) ||
                   (type == SolverType::FDPF_BX_NICSLU) ||
                   (type == SolverType::FDPF_XB_CKTSO) ||
                   (type == SolverType::FDPF_BX_CKTSO);
        }

        SolverType get_type() const { return _solver_type; }

        // Convenience accessors for the two FDPF SparseLU variants.
        // These exist mainly for internal diagnostic use (exposed as AnySolver.get_fdpf_*).
        FDPF_XB_SparseLUSolver& get_fdpf_xb_lu() {
            FDPF_XB_SparseLUSolver* p = dynamic_cast<FDPF_XB_SparseLUSolver*>(_solver.get());
            if (!p) throw std::runtime_error("ChooseSolver::get_fdpf_xb_lu: current solver is not FDPF_XB_SparseLU");
            return *p;
        }
        FDPF_BX_SparseLUSolver& get_fdpf_bx_lu() {
            FDPF_BX_SparseLUSolver* p = dynamic_cast<FDPF_BX_SparseLUSolver*>(_solver.get());
            if (!p) throw std::runtime_error("ChooseSolver::get_fdpf_bx_lu: current solver is not FDPF_BX_SparseLU");
            return *p;
        }

        void set_gridmodel(const GridModel* gridmodel) {
            _gridmodel_ptr = gridmodel;
            if (_solver) _solver->set_gridmodel(gridmodel);
        }

        bool ac_solver_used() const {
            return get_prt_solver("ac_solver_used", false)->IS_AC;
        }

        // Change solver by enum value (backward-compatible API).
        // Converts to the registry name via SolverTypeNames and calls the string overload.
        void change_solver(const SolverType& type);

        // Change solver by registry name (primary API; supports plugin solvers).
        void change_solver(const std::string& name);

        void reset() {
            get_prt_solver("reset", false)->reset();
        }

        bool compute_pf(const Eigen::SparseMatrix<cplx_type>& Ybus,
                        CplxVect& V,
                        const CplxVect& Sbus,
                        Eigen::Ref<const IntVect> slack_ids,
                        const RealVect& slack_weights,
                        Eigen::Ref<const IntVect> pv,
                        Eigen::Ref<const IntVect> pq,
                        int max_iter,
                        real_type tol)
        {
            _type_used_for_nr = _solver_type;
            return get_prt_solver("compute_pf", true)->compute_pf(
                Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, max_iter, tol);
        }

        Eigen::Ref<const CplxVect> get_V() const {
            return get_prt_solver("get_V", true)->get_V();
        }
        Eigen::Ref<const RealVect> get_Va() const {
            return get_prt_solver("get_Va", true)->get_Va();
        }
        Eigen::Ref<const RealVect> get_Vm() const {
            return get_prt_solver("get_Vm", true)->get_Vm();
        }

        Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const {
            check_right_solver("get_J");
            return get_prt_solver("get_J", false)->get_J();
        }

        RealMat get_ptdf() {
            if (!is_dc(_solver_type)) {
                throw std::runtime_error("ChooseSolver::get_ptdf: cannot get ptdf for a solver that is not DC.");
            }
            return get_prt_solver("get_ptdf", true)->get_ptdf();
        }

        RealMat get_lodf(const IntVect& from_bus, const IntVect& to_bus) {
            if (!is_dc(_solver_type)) {
                throw std::runtime_error("ChooseSolver::get_lodf: cannot get lodf for a solver that is not DC.");
            }
            return get_prt_solver("get_lodf", true)->get_lodf(from_bus, to_bus);
        }

        void update_internal_Ybus(const Coeff& new_coeffs, bool add) {
            if (!is_dc(_solver_type)) {
                throw std::runtime_error("ChooseSolver::update_internal_Ybus: cannot update internal Ybus for a solver that is not DC.");
            }
            get_prt_solver("update_internal_Ybus", true)->update_internal_Ybus(new_coeffs, add);
        }

        void tell_solver_control(const SolverControl& solver_control) {
            get_prt_solver("tell_solver_control", false)->tell_solver_control(solver_control);
        }

        Eigen::SparseMatrix<real_type> get_J_python() const {
            Eigen::SparseMatrix<real_type> res = get_J();
            return res;
        }

        double get_computation_time() const {
            return std::get<3>(get_prt_solver("get_computation_time", true)->get_timers());
        }

        std::tuple<double, double, double, double> get_timers() const {
            return get_prt_solver("get_timers", true)->get_timers();
        }

        TimerJacType get_timers_jacobian() const {
            return get_prt_solver("get_timers_jacobian", true)->get_timers_jacobian();
        }

        TimerPTDFLODFType get_timers_ptdf_lodf() const {
            return get_prt_solver("get_timers_ptdf_lodf", true)->get_timers_ptdf_lodf();
        }

        ErrorType get_error() const {
            return get_prt_solver("get_error", true)->get_error();
        }

        int get_nb_iter() const {
            return get_prt_solver("get_nb_iter", true)->get_nb_iter();
        }

        bool converged() const {
            return get_prt_solver("converged", true)->converged();
        }

        AlgoConfig get_config() const {
            return get_prt_solver("get_config", false)->get_config();
        }
        void set_config(const AlgoConfig& cfg) {
            get_prt_solver("set_config", false)->set_config(cfg);
        }

    protected:
        const BaseAlgo* get_prt_solver(const std::string& error_msg, bool check_right_solver_ = true) const {
            if (check_right_solver_) check_right_solver(error_msg);
            if (!_solver) throw std::runtime_error("ChooseSolver: no solver is active (not initialized?)");
            return _solver.get();
        }
        BaseAlgo* get_prt_solver(const std::string& error_msg, bool check_right_solver_ = true) {
            if (check_right_solver_) check_right_solver(error_msg);
            if (!_solver) throw std::runtime_error("ChooseSolver: no solver is active (not initialized?)");
            return _solver.get();
        }

    private:
        void check_right_solver(const std::string& error_msg) const {
            if (_solver_type != _type_used_for_nr) {
                std::ostringstream exc_;
                exc_ << "ChooseSolver: Solver mismatch when calling '";
                exc_ << error_msg;
                exc_ << "': current solver (";
                exc_ << _solver_type;
                exc_ << ") is not the one used to perform a powerflow (";
                exc_ << _type_used_for_nr;
                exc_ << ").";
                throw std::runtime_error(exc_.str());
            }
        }

        std::unique_ptr<BaseAlgo> _solver;
        SolverType _solver_type;
        SolverType _type_used_for_nr;
        const GridModel* _gridmodel_ptr;
};


} // namespace ls2g

#endif  // CHOOSESOLVER_H