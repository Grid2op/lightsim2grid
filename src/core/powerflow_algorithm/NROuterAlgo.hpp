// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_OUTER_ALGO_H
#define NR_OUTER_ALGO_H

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "NRAlgo.hpp"
#include "outer_loop/BaseOuterLoop.hpp"
#include "outer_loop/OuterLoopDriver.hpp"

namespace ls2g {

/**
 * Single-slack Newton-Raphson with OpenLoadFlow's outer loops around it (NROuter_*).
 *
 * The driver is OpenLoadFlow's AcloadFlowEngine.run (v2.3.0):
 *   1. keep the loops whose is_needed holds and initialize them, in order;
 *   2. first Newton solve; the loops only run if it converged;
 *   3. passes over the loops: each one repeats check(), re-solving after every UNSTABLE
 *      result, until it is stable, a solve fails or the iteration cap is reached. A pass
 *      stops before the loop that was the last one to be unstable (every other loop has
 *      then been checked since its change); passes repeat while one added Newton
 *      iterations;
 *   4. FAILED stops everything;
 *   5. the unrealistic-voltage check runs on every solve, or -- in robust mode -- only
 *      from the end of the last loop able to fix an unrealistic state;
 *   6. cleanup() in reverse order.
 *
 * Between two solves the state is kept: voltages, controllers, the sparsity of J and its
 * factorization; a loop edits values only (the private Sbus copy, ...), and the next solve
 * refactorizes on its first iteration. The sparsity holds the union of what every loop
 * declared, so a whole solve is one symbolic analysis. The refactorization fallback is
 * always on: a value edit may shrink a pivot the fixed pivot sequence then finds at zero.
 *
 * The loop list is the grid's (LSGrid::add_outer_loop), handed over before each solve
 * through set_outer_loops; the algorithm runs clones of it.
 */
template<class LinearSolver>
class NROuterAlgo final : public NRAlgo<LinearSolver, SingleSlackNRSystem>
{
    using Parent = NRAlgo<LinearSolver, SingleSlackNRSystem>;

public:
    NROuterAlgo() noexcept : Parent() {
        this->_linear_solver.set_refactor_fallback(true);
    }
    ~NROuterAlgo() noexcept override = default;

    // ----- outer loops -----------------------------------------------------------

    bool supports_outer_loops() const noexcept override { return true; }

    void set_outer_loops(const std::vector<std::shared_ptr<const BaseOuterLoop> > & loops) override {
        std::vector<std::string> signature;
        signature.reserve(loops.size());
        for(const auto & loop : loops) signature.push_back(_signature(*loop));
        if(signature == signature_) return;  // same list: keep the loops (and the sparsity)
        loops_.clear();
        loops_.reserve(loops.size());
        for(const auto & loop : loops) loops_.push_back(loop->clone());
        signature_ = signature;
        // what the loops reserve is part of the sparsity: rebuild it on the next solve
        this->need_factorize_ = true;
    }

    OuterLoopStats get_outer_loop_stats() const override { return stats_; }

    void get_outer_target_p(std::vector<real_type> & gen_p_mw,
                            std::vector<real_type> & storage_p_mw) const override {
        gen_p_mw = state_.gen_target_p;
        storage_p_mw = state_.storage_target_p;
    }
    void get_outer_hvdc_status(std::vector<int> & status) const override {
        status = state_.hvdc_status;
    }
    void get_outer_phase_tap(std::vector<int> & positions) const override {
        positions.clear();
        const BranchControl * phase = this->_system.branch_control();
        if (phase == nullptr) return;
        phase->positions(positions);  // TAP_KEEP for the transformers it does not handle
    }
    void get_outer_ratio_tap(std::vector<int> & positions) const override {
        positions.clear();
        const BranchControl * branch = this->_system.branch_control();
        if (branch != nullptr) branch->ratio_positions(positions);
    }

    const OuterLoopDriverParams & get_driver_params() const { return params_; }
    void set_driver_params(const OuterLoopDriverParams & params) {
        _check_driver_params(params);
        params_ = params;
    }

    // ----- PV / PQ relabelling at constant sparsity --------------------------------
    // A caller (a batch's generator contingencies) and the loops may both reserve
    // switchable buses: the system gets the union of the two, see _before_init_topology.
    void set_switchable_vm_buses(const std::vector<int> & solver_bus_ids) override {
        caller_switchable_ = solver_bus_ids;
        this->_system.set_switchable_vm_buses(solver_bus_ids);
    }
    // a caller's pinned buses (a batch's PV buses among its switchable ones) are pinned on
    // top of the loops' own, see _solve
    void set_pv_pinned_buses(const std::vector<int> & solver_bus_ids) override {
        caller_pinned_ = solver_bus_ids;
    }

    // ----- AlgoConfig: the Newton's parameters, then the driver's ------------------
    // int_params:  [the Newton's 4], max_outer_iterations, robust_mode
    // real_params: [the Newton's 6], min_realistic_voltage, max_realistic_voltage,
    //              min_nominal_voltage_realistic_check
    AlgoConfig get_config() const override {
        AlgoConfig cfg = Parent::get_config();
        cfg.int_params.push_back(params_.max_outer_iterations);
        cfg.int_params.push_back(params_.voltage_remote_control_robust_mode ? 1 : 0);
        cfg.real_params.push_back(static_cast<double>(params_.min_realistic_voltage));
        cfg.real_params.push_back(static_cast<double>(params_.max_realistic_voltage));
        cfg.real_params.push_back(static_cast<double>(params_.min_nominal_voltage_realistic_check));
        return cfg;
    }

    void set_config(const AlgoConfig & cfg) override {
        // a config of a plain NR_* (shorter) leaves the driver's parameters as they are
        OuterLoopDriverParams params = params_;
        if(cfg.int_params.size() >= 6){
            params.max_outer_iterations = cfg.int_params[4];
            params.voltage_remote_control_robust_mode = cfg.int_params[5] != 0;
        }
        if(cfg.real_params.size() >= 9){
            params.min_realistic_voltage = static_cast<real_type>(cfg.real_params[6]);
            params.max_realistic_voltage = static_cast<real_type>(cfg.real_params[7]);
            params.min_nominal_voltage_realistic_check = static_cast<real_type>(cfg.real_params[8]);
        }
        // all or nothing, like NRAlgo::set_config: validate the driver's part, apply the
        // Newton's (which validates itself before writing anything), then the driver's
        _check_driver_params(params);
        Parent::set_config(cfg);
        params_ = params;
    }

    // ----- powerflow -------------------------------------------------------------

    bool compute_pf(
        const EigenRefConstCplxSpMat     & Ybus,
        const Eigen::Ref<const CplxVect> & V,
        const Eigen::Ref<const CplxVect> & Sbus,
        const Eigen::Ref<const IntVect>  & slack_ids,
        const Eigen::Ref<const RealVect> & slack_weights,
        const Eigen::Ref<const IntVect>  & pv,
        const Eigen::Ref<const IntVect>  & pq,
        int                              max_iter,
        real_type                        tol
    ) override;

    void reset() override {
        Parent::reset();
        stats_.clear();
    }

protected:
    // the loops claim their slots before the system registers its rows / columns
    void _before_init_topology() override;

private:
    // the Newton from the current state, then the voltages published and, if asked, the
    // unrealistic-voltage check. Returns whether the solve is usable by the loops.
    bool _solve(int max_iter, real_type tol, bool & need_init, bool check_unrealistic);
    // OuterContext on the algorithm's current state
    OuterContext _context(OuterState * state);
    std::vector<bool> _vm_unknown() const;

    static void _check_driver_params(const OuterLoopDriverParams & params) {
        if(params.max_outer_iterations < 0){
            throw std::runtime_error("NROuterAlgo: max_outer_iterations must be >= 0.");
        }
        if(!(params.min_realistic_voltage < params.max_realistic_voltage)){
            throw std::runtime_error("NROuterAlgo: min_realistic_voltage must be lower than "
                                     "max_realistic_voltage.");
        }
    }

    static std::string _signature(const BaseOuterLoop & loop) {
        std::string res = loop.name();
        const AlgoConfig params = loop.get_params();
        for(int v : params.int_params) res += "|" + std::to_string(v);
        for(double v : params.real_params) res += "|" + std::to_string(v);
        return res;
    }

    std::vector<std::unique_ptr<BaseOuterLoop> > loops_;
    std::vector<int> caller_switchable_;  // see set_switchable_vm_buses
    std::vector<int> caller_pinned_;      // see set_pv_pinned_buses
    // the buses the loops declared switchable (they are PV unless a loop made them PQ,
    // OuterState::pq_buses), and the ones pinned in the last solve
    std::vector<int> declared_switchable_;
    std::vector<int> pinned_;
    std::vector<std::string> signature_;
    OuterLoopDriverParams params_;
    OuterLoopStats stats_;
    int total_nr_iterations_ = 0;
    bool unrealistic_ = false;

    // private copies the loops edit; the system reads Sbus_ through a pointer, so it must
    // not move during a solve (it is only reassigned before _setup)
    CplxVect Sbus_;
    CplxVect Sbus_init_;
    CplxVect Sbus_target_;  // see OuterState::Sbus_target
    Eigen::SparseMatrix<cplx_type> Ybus_;  // the grid's, its values patched by the phase shifters
    RealVect controller_q_;  // the context's copy, refreshed with it
    int slack_bus_ = -1;     // solver id of the slack bus of the current solve
    OuterState state_;       // what the loops edit; kept after the solve for its results
};

template<class LinearSolver>
OuterContext NROuterAlgo<LinearSolver>::_context(OuterState * state)
{
    OuterContext ctx;
    ctx.grid = this->lsgrid_ptr_;
    ctx.V = &this->V_;
    ctx.Va = &this->Va_;
    ctx.bus_mismatch = &this->mis_bus_;
    controller_q_ = this->_system.controller_q();
    ctx.controller_q = &controller_q_;
    ctx.slack_bus = slack_bus_;
    ctx.slack_absorbed = this->_system.slack_absorbed();
    ctx.state = state;
    ctx.branch_control = this->_system.branch_control();
    return ctx;
}

#include "NROuterAlgo.tpp"

}  // namespace ls2g

#endif  // NR_OUTER_ALGO_H
