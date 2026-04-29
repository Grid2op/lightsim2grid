// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ALGORITHMSELECTOR_H
#define ALGORITHMSELECTOR_H

#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Concrete algorithm typedefs (NR_SparseLU, NR_KLU, …) + Solvers.hpp guards
#include "Solvers.hpp"
// Registry for solver creation by name
#include "AlgorithmRegistry.hpp"
#include "AlgoConfig.hpp"

namespace ls2g {

enum class LS2G_API AlgorithmType {NR_SparseLU, NR_KLU, GaussSeidel, DC_SparseLU, GaussSeidelSynch, NR_NICSLU,
                       NRSing_SparseLU, NRSing_KLU, NRSing_NICSLU,
                       DC_KLU, DC_NICSLU,
                       NR_CKTSO, NRSing_CKTSO, DC_CKTSO,
                       FDPF_XB_SparseLU, FDPF_BX_SparseLU, // from 0.7.5
                       FDPF_XB_KLU, FDPF_BX_KLU,           // from 0.7.5
                       FDPF_XB_NICSLU, FDPF_BX_NICSLU,     // from 0.7.5
                       FDPF_XB_CKTSO, FDPF_BX_CKTSO,       // from 0.7.5
                       Custom                               // external/plugin solvers
                       };

std::ostream& operator<<(std::ostream& out, const AlgorithmType& algo_type);

class LS2G_API AlgorithmSelector final
{
    public:
        // Default-constructs a SparseLU solver via the registry.
        // The registry must be populated (call register_builtin_solvers) before
        // any AlgorithmSelector instance is created — the PYBIND11_MODULE init does this.
        AlgorithmSelector();

        ~AlgorithmSelector() noexcept = default;

        AlgorithmSelector(const AlgorithmSelector&) = delete;
        AlgorithmSelector(AlgorithmSelector&&) = delete;
        AlgorithmSelector& operator=(AlgorithmSelector&&) = delete;
        AlgorithmSelector& operator=(const AlgorithmSelector&) = delete;

        // Returns the enum-typed solvers available in this build.
        // Does not include plugin (Custom) solvers; use AlgorithmRegistry::available_default_algorithms()
        // for the full list.
        std::vector<AlgorithmType> available_default_algorithms() const
        {
            std::vector<AlgorithmType> res;
            res.reserve(22);
            res.push_back(AlgorithmType::NR_SparseLU);
            res.push_back(AlgorithmType::GaussSeidel);
            res.push_back(AlgorithmType::DC_SparseLU);
            res.push_back(AlgorithmType::GaussSeidelSynch);
            res.push_back(AlgorithmType::NRSing_SparseLU);
            res.push_back(AlgorithmType::FDPF_XB_SparseLU);
            res.push_back(AlgorithmType::FDPF_BX_SparseLU);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(AlgorithmType::NR_KLU);
                res.push_back(AlgorithmType::NRSing_KLU);
                res.push_back(AlgorithmType::DC_KLU);
                res.push_back(AlgorithmType::FDPF_XB_KLU);
                res.push_back(AlgorithmType::FDPF_BX_KLU);
            #endif
            #ifdef NICSLU_SOLVER_AVAILABLE
                res.push_back(AlgorithmType::NR_NICSLU);
                res.push_back(AlgorithmType::NRSing_NICSLU);
                res.push_back(AlgorithmType::DC_NICSLU);
                res.push_back(AlgorithmType::FDPF_XB_NICSLU);
                res.push_back(AlgorithmType::FDPF_BX_NICSLU);
            #endif
            #ifdef CKTSO_SOLVER_AVAILABLE
                res.push_back(AlgorithmType::NR_CKTSO);
                res.push_back(AlgorithmType::NRSing_CKTSO);
                res.push_back(AlgorithmType::DC_CKTSO);
                res.push_back(AlgorithmType::FDPF_XB_CKTSO);
                res.push_back(AlgorithmType::FDPF_BX_CKTSO);
            #endif
            return res;
        }

        bool is_dc(const AlgorithmType& type) const noexcept {
            return (type == AlgorithmType::DC_SparseLU) ||
                   (type == AlgorithmType::DC_KLU) ||
                   (type == AlgorithmType::DC_NICSLU) ||
                   (type == AlgorithmType::DC_CKTSO);
        }
        bool is_fdpf(const AlgorithmType& type) const noexcept {
            return (type == AlgorithmType::FDPF_XB_SparseLU) ||
                   (type == AlgorithmType::FDPF_BX_SparseLU) ||
                   (type == AlgorithmType::FDPF_XB_KLU) ||
                   (type == AlgorithmType::FDPF_BX_KLU) ||
                   (type == AlgorithmType::FDPF_XB_NICSLU) ||
                   (type == AlgorithmType::FDPF_BX_NICSLU) ||
                   (type == AlgorithmType::FDPF_XB_CKTSO) ||
                   (type == AlgorithmType::FDPF_BX_CKTSO);
        }

        AlgorithmType get_type() const { return _algo_type; }

        // Convenience accessors for the two FDPF SparseLU variants.
        // These exist mainly for internal diagnostic use (exposed as AlgorithmSelector.get_fdpf_*).
        // Defined in AlgorithmSelector.cpp to avoid implicit instantiation of BaseFDPFAlgo in every TU.
        FDPF_XB_SparseLU& get_fdpf_xb_lu();
        FDPF_BX_SparseLU& get_fdpf_bx_lu();

        void set_gridmodel(const GridModel* gridmodel) {
            _gridmodel_ptr = gridmodel;
            if (_algo) _algo->set_gridmodel(gridmodel);
        }

        bool ac_solver_used() const {
            return get_prt_solver("ac_solver_used", false)->IS_AC;
        }

        // Change solver by enum value (backward-compatible API).
        // Converts to the registry name via AlgorithmTypeNames and calls the string overload.
        void change_algorithm(const AlgorithmType& type);

        // Change solver by registry name (primary API; supports plugin solvers).
        void change_algorithm(const std::string& name);

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
            _algo_type_used_for_nr = _algo_type;
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
            if (!is_dc(_algo_type)) {
                throw std::runtime_error("AlgorithmSelector::get_ptdf: cannot get ptdf for a solver that is not DC.");
            }
            return get_prt_solver("get_ptdf", true)->get_ptdf();
        }

        RealMat get_lodf(const IntVect& from_bus, const IntVect& to_bus) {
            if (!is_dc(_algo_type)) {
                throw std::runtime_error("AlgorithmSelector::get_lodf: cannot get lodf for a solver that is not DC.");
            }
            return get_prt_solver("get_lodf", true)->get_lodf(from_bus, to_bus);
        }

        void update_internal_Ybus(const Coeff& new_coeffs, bool add) {
            if (!is_dc(_algo_type)) {
                throw std::runtime_error("AlgorithmSelector::update_internal_Ybus: cannot update internal Ybus for a solver that is not DC.");
            }
            get_prt_solver("update_internal_Ybus", true)->update_internal_Ybus(new_coeffs, add);
        }

        void tell_solver_control(const AlgoControl& solver_control) {
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
            if (!_algo) throw std::runtime_error("AlgorithmSelector: no solver is active (not initialized?)");
            return _algo.get();
        }
        BaseAlgo* get_prt_solver(const std::string& error_msg, bool check_right_solver_ = true) {
            if (check_right_solver_) check_right_solver(error_msg);
            if (!_algo) throw std::runtime_error("AlgorithmSelector: no solver is active (not initialized?)");
            return _algo.get();
        }

    private:
        void check_right_solver(const std::string& error_msg) const {
            if (_algo_type != _algo_type_used_for_nr) {
                std::ostringstream exc_;
                exc_ << "AlgorithmSelector: Solver mismatch when calling '";
                exc_ << error_msg;
                exc_ << "': current solver (";
                exc_ << _algo_type;
                exc_ << ") is not the one used to perform a powerflow (";
                exc_ << _algo_type_used_for_nr;
                exc_ << ").";
                throw std::runtime_error(exc_.str());
            }
        }

        std::unique_ptr<BaseAlgo> _algo;
        AlgorithmType _algo_type;
        AlgorithmType _algo_type_used_for_nr;
        const GridModel* _gridmodel_ptr;
};


} // namespace ls2g

#endif  // ALGORITHMSELECTOR_H
