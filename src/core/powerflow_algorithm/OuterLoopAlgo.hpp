// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_LOOP_ALGO_H
#define OUTER_LOOP_ALGO_H

#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "BaseAlgo.hpp"
#include "outer_loop/BaseOuterLoop.hpp"
#include "outer_loop/OuterLoopDriver.hpp"

namespace ls2g {

/**
 * OpenLoadFlow's outer loops around an algorithm it wraps (the NROuter_* family wraps a
 * single-slack Newton-Raphson, see NROuterInner).
 *
 * The driver is OpenLoadFlow's AcloadFlowEngine.run (v2.3.0):
 *   1. keep the loops whose is_needed holds and initialize them, in order;
 *   2. first solve; the loops only run if it converged;
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
 * Between two solves the inner algorithm keeps its state: voltages, controllers, the
 * sparsity of J and its factorization; a loop edits values only (the private Sbus copy,
 * ...), and the next solve refactorizes on its first iteration. The sparsity holds the
 * union of what every loop declared, so a whole solve is one symbolic analysis.
 *
 * `Inner` is what the driver needs from the algorithm beyond BaseAlgo's interface (see
 * NROuterInner, the one there is): algo() (the wrapped BaseAlgo), controls() (the
 * OuterControls the loops reserve and act through), begin(), setup(..., declare), newton(),
 * finalize(), request_rebuild(), apply_state(), fill_context(),
 * vm_unknown(), and the phase_tap / ratio_tap / shunt_sections results. Everything else
 * BaseAlgo declares is forwarded to algo(), and the results of each solve are copied back
 * into this object, which is what LSGrid reads.
 *
 * The loop list is the grid's (LSGrid::add_outer_loop), handed over before each solve
 * through set_outer_loops; the algorithm runs clones of it.
 */
template<class Inner>
class OuterLoopAlgo final : public BaseAlgo
{
public:
    OuterLoopAlgo() noexcept : BaseAlgo(true) {}
    ~OuterLoopAlgo() noexcept override = default;

    // the wrapped algorithm (the bindings reach the Newton's own parameters through it)
    Inner & inner() { return inner_; }
    const Inner & inner() const { return inner_; }

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
        inner_.request_rebuild();
    }

    OuterLoopStats get_outer_loop_stats() const override { return stats_; }

    void get_outer_target_p(std::vector<real_type> & gen_p_mw,
                            std::vector<real_type> & storage_p_mw) const override {
        gen_p_mw = state_.gen_target_p;
        storage_p_mw = state_.storage_target_p;
    }
    void get_outer_hvdc_status(std::vector<int> & status) const override {
        status = inner_.controls().hvdc_regimes();
    }
    void get_outer_phase_tap(std::vector<int> & positions) const override { inner_.phase_tap(positions); }
    void get_outer_shunt_sections(std::vector<int> & counts) const override { inner_.shunt_sections(counts); }
    void get_outer_ratio_tap(std::vector<int> & positions) const override { inner_.ratio_tap(positions); }

    const OuterLoopDriverParams & get_driver_params() const { return params_; }
    void set_driver_params(const OuterLoopDriverParams & params) {
        _check_driver_params(params);
        params_ = params;
    }

    // ----- PV / PQ relabelling at constant sparsity --------------------------------
    // A caller (a batch's generator contingencies) and the loops may both reserve
    // switchable buses: the controls hold the caller's next to the loops' own.
    void set_switchable_vm_buses(const std::vector<int> & solver_bus_ids) override {
        inner_.controls().set_caller_switchable(solver_bus_ids);
        inner_.algo().set_switchable_vm_buses(solver_bus_ids);
    }
    void set_pv_pinned_buses(const std::vector<int> & solver_bus_ids) override {
        inner_.controls().set_caller_pinned(solver_bus_ids);
    }
    // the fallback stays on whatever a caller asks: a value edit between two solves may
    // shrink a pivot the fixed pivot sequence then finds at zero (see NROuterInner)
    void set_refactor_fallback(bool /*val*/) override {}

    // ----- AlgoConfig: the inner algorithm's parameters, then the driver's ------------
    // int_params:  [the Newton's 4], max_outer_iterations, robust_mode
    // real_params: [the Newton's 6], min_realistic_voltage, max_realistic_voltage,
    //              min_nominal_voltage_realistic_check
    AlgoConfig get_config() const override {
        AlgoConfig cfg = inner_.algo().get_config();
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
        // inner algorithm's (which validates itself before writing anything), then the driver's
        _check_driver_params(params);
        inner_.algo().set_config(cfg);
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
        BaseAlgo::reset();
        inner_.algo().reset();
        err_ = inner_.algo().get_error();
        stats_.clear();
    }

    void set_lsgrid(const LSGrid * gridmodel) override {
        BaseAlgo::set_lsgrid(gridmodel);
        inner_.algo().set_lsgrid(gridmodel);
    }

    // ----- forwarded to the inner algorithm ------------------------------------------

    bool supports_hvdc_droop() const noexcept override { return inner_.algo().supports_hvdc_droop(); }
    bool supports_remote_voltage_control() const noexcept override { return inner_.algo().supports_remote_voltage_control(); }
    bool fills_bus_mismatch() const noexcept override { return inner_.algo().fills_bus_mismatch(); }

    EigenRefConstRealSpMat get_J() const override { return inner_.algo().get_J(); }
    IntVect get_theta_to_J_col_python() const override { return inner_.algo().get_theta_to_J_col_python(); }
    IntVect get_vm_to_J_col_python()    const override { return inner_.algo().get_vm_to_J_col_python(); }
    IntVect get_q_to_J_col_python()     const override { return inner_.algo().get_q_to_J_col_python(); }
    IntVect get_p_to_J_row_python()     const override { return inner_.algo().get_p_to_J_row_python(); }
    IntVect get_q_to_J_row_python()     const override { return inner_.algo().get_q_to_J_row_python(); }
    IntVect get_p_buses_python()        const override { return inner_.algo().get_p_buses_python(); }
    IntVect get_p_rows_python()         const override { return inner_.algo().get_p_rows_python(); }
    IntVect get_q_buses_python()        const override { return inner_.algo().get_q_buses_python(); }
    IntVect get_q_rows_python()         const override { return inner_.algo().get_q_rows_python(); }
    IntVect get_theta_buses_python()    const override { return inner_.algo().get_theta_buses_python(); }
    IntVect get_theta_cols_python()     const override { return inner_.algo().get_theta_cols_python(); }
    IntVect get_vm_buses_python()       const override { return inner_.algo().get_vm_buses_python(); }
    IntVect get_vm_cols_python()        const override { return inner_.algo().get_vm_cols_python(); }

    RealVect  get_controller_q()       const override { return inner_.algo().get_controller_q(); }
    IntVect   get_controller_kind()    const override { return inner_.algo().get_controller_kind(); }
    IntVect   get_controller_elem_id() const override { return inner_.algo().get_controller_elem_id(); }
    IntVect   get_controller_q_col()   const override { return inner_.algo().get_controller_q_col(); }
    IntVect   get_group_v_row()        const override { return inner_.algo().get_group_v_row(); }
    int       get_slack_col()          const override { return inner_.algo().get_slack_col(); }
    real_type get_slack_absorbed()     const override { return inner_.algo().get_slack_absorbed(); }

    // the inner algorithm's timers, over every solve of the last compute_pf, with the
    // driver's total
    TimerJac get_timers_jacobian() const override {
        TimerJac res = inner_.algo().get_timers_jacobian();
        res.timer_total_nr_ = timer_total_nr_;
        return res;
    }
    LinearSolverStats get_linear_solver_stats() const override { return inner_.algo().get_linear_solver_stats(); }

    bool supports_bus_masking() const override { return inner_.algo().supports_bus_masking(); }
    bool supports_jacobian() const override { return inner_.algo().supports_jacobian(); }
    void refresh_J_at_solution() override { inner_.algo().refresh_J_at_solution(); }
    void set_masked_buses(const std::vector<int> & solver_bus_ids) override { inner_.algo().set_masked_buses(solver_bus_ids); }
    void set_may_mask_voltage_control(bool val) override { inner_.algo().set_may_mask_voltage_control(val); }
    void set_voltage_control_v_set(const Eigen::Ref<const RealVect> & v_set) override { inner_.algo().set_voltage_control_v_set(v_set); }
    bool supports_pv_pinning() const override { return inner_.algo().supports_pv_pinning(); }
    void set_start_polar_cache(bool val) override { inner_.algo().set_start_polar_cache(val); }

    bool supports_cpf() const noexcept override { return inner_.algo().supports_cpf(); }
    bool cpf_tangent(const Eigen::Ref<const CplxVect> & dir_solver, RealVect & z) override {
        const bool res = inner_.algo().cpf_tangent(dir_solver, z);
        if(!res) err_ = inner_.algo().get_error();
        return res;
    }
    void cpf_predict(const Eigen::Ref<const RealVect> & z, real_type coeff, CplxVect & V_pred) const override {
        inner_.algo().cpf_predict(z, coeff, V_pred);
    }
    bool cpf_refactorize_at_current() override {
        const bool res = inner_.algo().cpf_refactorize_at_current();
        if(!res) err_ = inner_.algo().get_error();
        return res;
    }

private:
    // the inner solve from the current state, then its results copied here and, if asked,
    // the unrealistic-voltage check. Returns whether the solve is usable by the loops.
    bool _solve(int max_iter, real_type tol, bool & need_init, bool check_unrealistic);
    // the loops' declarations, on a rebuild of the topology (see Inner::setup)
    OuterDeclaration _declare();
    // OuterContext on the algorithm's current state
    OuterContext _context(OuterState * state);

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

    Inner inner_;

    std::vector<std::unique_ptr<BaseOuterLoop> > loops_;
    std::vector<std::string> signature_;
    OuterLoopDriverParams params_;
    OuterLoopStats stats_;
    int total_nr_iterations_ = 0;
    bool unrealistic_ = false;

    // private copies the loops edit; the inner algorithm reads Sbus_ and Ybus_ through a
    // pointer, so they must not move during a solve (they are only reassigned before setup)
    CplxVect Sbus_;
    CplxVect Sbus_init_;
    CplxVect Sbus_target_;  // see OuterState::Sbus_target
    Eigen::SparseMatrix<cplx_type> Ybus_;  // the grid's, its values patched by the phase shifters
    RealVect controller_q_;  // the context's copy, refreshed with it
    int slack_bus_ = -1;     // solver id of the slack bus of the current solve
    OuterState state_;       // what the loops edit; kept after the solve for its results
};

#include "OuterLoopAlgo.tpp"

}  // namespace ls2g

#endif  // OUTER_LOOP_ALGO_H
