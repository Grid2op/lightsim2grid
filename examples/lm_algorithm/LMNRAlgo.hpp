// Copyright (c) 2026
// Experimental plugin, not part of the lightsim2grid core distribution.
//
// LMNRAlgo: Newton-Raphson with Levenberg-Marquardt diagonal damping,
// (J + lambda*D) dx = -F instead of J dx = -F, lambda adapted by an
// accept/reject loop on the actual (nonlinear) mismatch reduction of each
// trial step.
//
// Motivation: on some real grids, clusters that are high-resistance,
// heavily loaded and short of local generation make the local Jacobian
// sub-block ill-conditioned/near-singular. Plain NR can diverge there, and
// so can the Iwamoto optimal multiplier (ScalingPolicyType::Iwamoto) or the
// MaxVoltageChange clamp (ScalingPolicyType::MaxVoltageChange): both of
// those only ever rescale the *same* Newton direction J^-1 F, which is
// already the problem when J itself is ill-conditioned. LM instead
// regularizes the direction, by damping J's diagonal before solving.
//
// This is deliberately a SIBLING of NRAlgo, not built on top of it: the
// accept/reject inner loop (potentially several factorize+solve+mismatch-eval
// cycles per outer NR iteration) does not fit NRAlgo's single
// solve-per-iteration control flow, nor its RefactorPolicy (which assumes J's
// numeric values are stable enough across iterations to reuse a
// factorization -- false here, since lambda legitimately changes) or
// ScalingPolicy (whose scale() hook only rescales an already-solved Newton
// step post-hoc; it cannot touch the matrix being solved).
//
// The diagonal damping itself is implemented as a NRSystem extension
// (DiagonalEntry, see DiagonalDamping.hpp) plus one small, generic addition
// to core (NRSystem::update_trailing_feature_values) that directly pokes the
// n reserved diagonal J_.valuePtr() positions -- bypassing
// fill_internal_variables()/fill_J() entirely, so a rejected trial costs
// only a numeric refactor + solve + (allocation-free) mismatch evaluation,
// not a full O(nnz) Jacobian rebuild.
//
// Build/load, once compiled (see CMakeLists.txt in this directory):
//   import lightsim2grid
//   lightsim2grid.load_algorithm_plugin("examples/lm_algorithm/build/liblm_algorithm.so")
//   grid.change_algorithm("NR_LM_KLU")

#ifndef LM_NR_ALGO_H
#define LM_NR_ALGO_H

#include <Solvers.hpp>
#include <powerflow_algorithm/BaseAlgo.hpp>
#include <powerflow_algorithm/NRSystem.hpp>
#include <linear_solvers/LinearSolverStats.hpp>

#include "DiagonalDamping.hpp"

namespace ls2g {

template<class LinearSolver, class NRSystemT>
class LMNRAlgo final : public BaseAlgo
{
public:
    using DampedSystem = typename with_diagonal_entry<NRSystemT>::type;

    LMNRAlgo() noexcept :
        BaseAlgo(true),
        need_factorize_(true),
        lambda_(static_cast<real_type>(0.)),
        lambda_init_(static_cast<real_type>(1e-6)),
        nu_(static_cast<real_type>(2.)),
        lambda_min_(static_cast<real_type>(0.)),
        lambda_max_(static_cast<real_type>(1e8)),
        lambda_eps_(static_cast<real_type>(1e-12)),
        max_trials_(8),
        timer_dSbus_(0.),
        timer_fillJ_(0.) {}

    ~LMNRAlgo() noexcept override = default;

    // The underlying system still carries Hvdc + VoltageControl (same
    // reasoning as NRAlgo: both SingleSlackNRSystem and MultiSlackNRSystem
    // include them unconditionally).
    static constexpr bool SUPPORTS_HVDC_DROOP = true;
    bool supports_hvdc_droop() const noexcept override { return SUPPORTS_HVDC_DROOP; }
    static constexpr bool SUPPORTS_REMOTE_VOLTAGE_CONTROL = true;
    bool supports_remote_voltage_control() const noexcept override { return SUPPORTS_REMOTE_VOLTAGE_CONTROL; }

    // ----- Jacobian / VoltageControl introspection: forward to _system, same as NRAlgo -----
    Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const override { return _system.J(); }
    IntVect get_theta_to_J_col_python() const override { return _to_intvect(_system.theta_to_J_col()); }
    IntVect get_vm_to_J_col_python()    const override { return _to_intvect(_system.vm_to_J_col()); }
    IntVect get_q_to_J_col_python()     const override { return _to_intvect(_system.q_to_J_col()); }
    IntVect get_p_to_J_row_python()     const override { return _to_intvect(_system.p_to_J_row()); }
    IntVect get_q_to_J_row_python()     const override { return _to_intvect(_system.q_to_J_row()); }
    IntVect get_p_buses_python()        const override { return _to_intvect(_system.p_buses()); }
    IntVect get_p_rows_python()         const override { return _to_intvect(_system.p_rows()); }
    IntVect get_q_buses_python()        const override { return _to_intvect(_system.q_buses()); }
    IntVect get_q_rows_python()         const override { return _to_intvect(_system.q_rows()); }
    IntVect get_theta_buses_python()    const override { return _to_intvect(_system.theta_buses()); }
    IntVect get_theta_cols_python()     const override { return _to_intvect(_system.theta_cols()); }
    IntVect get_vm_buses_python()       const override { return _to_intvect(_system.vm_buses()); }
    IntVect get_vm_cols_python()        const override { return _to_intvect(_system.vm_cols()); }
    RealVect get_controller_q()         const override { return _system.controller_q(); }
    IntVect  get_controller_kind()      const override { return _system.controller_kind(); }
    IntVect  get_controller_elem_id()   const override { return _system.controller_elem_id(); }
    IntVect  get_controller_q_col()     const override { return _system.controller_q_col(); }
    int      get_slack_col()            const override { return _system.slack_col(); }
    real_type get_slack_absorbed()      const override { return _system.slack_absorbed(); }

    TimerJac get_timers_jacobian() const override
    {
        const LinearSolverStats lsstats = detail::get_stats_impl(_linear_solver, 0);
        TimerJac res;
        res.timer_Fx_         = timer_Fx_;
        res.timer_solve_      = lsstats.timer_solve_;
        res.timer_factor_     = lsstats.timer_factor_;
        res.timer_refactor_   = lsstats.timer_refactor_;
        res.timer_initialize_ = lsstats.timer_initialize_;
        res.timer_check_      = timer_check_;
        res.timer_dSbus_      = timer_dSbus_;
        res.timer_fillJ_      = timer_fillJ_;
        res.timer_total_nr_   = timer_total_nr_;
        return res;
    }
    LinearSolverStats get_linear_solver_stats() const override {
        return detail::get_stats_impl(_linear_solver, 0);
    }

    bool supports_bus_masking() const override { return true; }
    void set_masked_buses(const std::vector<int>& solver_bus_ids) override {
        _system.set_masked_buses(solver_bus_ids);
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
        _system.clear_jacobian();
        need_factorize_ = true;
        n_ = -1;
        ErrorType reset_status = _linear_solver.reset();
        if (reset_status != ErrorType::NoError) err_ = reset_status;
    }

    // ----- Levenberg-Marquardt tunables -------------------------------------------
    real_type get_lambda_init() const { return lambda_init_; }
    void      set_lambda_init(real_type v) { lambda_init_ = v; }
    real_type get_nu()          const { return nu_; }
    void      set_nu(real_type v) { nu_ = v; }
    real_type get_lambda_min()  const { return lambda_min_; }
    void      set_lambda_min(real_type v) { lambda_min_ = v; }
    real_type get_lambda_max()  const { return lambda_max_; }
    void      set_lambda_max(real_type v) { lambda_max_ = v; }
    real_type get_lambda_eps()  const { return lambda_eps_; }
    void      set_lambda_eps(real_type v) { lambda_eps_ = v; }
    int       get_max_trials()  const { return max_trials_; }
    void      set_max_trials(int v) { max_trials_ = v; }
    // Value lambda was left at when the last compute_pf call returned
    // (informational only -- inspect it to see how hard the last solve was).
    real_type get_last_lambda() const { return lambda_; }

    // int_params  = [max_trials_]
    // real_params = [lambda_init_, nu_, lambda_min_, lambda_max_, lambda_eps_]
    AlgoConfig get_config() const override {
        AlgoConfig cfg;
        cfg.int_params  = { max_trials_ };
        cfg.real_params = { static_cast<double>(lambda_init_),
                            static_cast<double>(nu_),
                            static_cast<double>(lambda_min_),
                            static_cast<double>(lambda_max_),
                            static_cast<double>(lambda_eps_) };
        return cfg;
    }
    void set_config(const AlgoConfig& cfg) override {
        if (cfg.int_params.size()  < 1) throw std::runtime_error("LMNRAlgo::set_config: int_params must have at least 1 element");
        if (cfg.real_params.size() < 5) throw std::runtime_error("LMNRAlgo::set_config: real_params must have at least 5 elements");
        max_trials_ = cfg.int_params[0];
        lambda_init_ = static_cast<real_type>(cfg.real_params[0]);
        nu_          = static_cast<real_type>(cfg.real_params[1]);
        lambda_min_  = static_cast<real_type>(cfg.real_params[2]);
        lambda_max_  = static_cast<real_type>(cfg.real_params[3]);
        lambda_eps_  = static_cast<real_type>(cfg.real_params[4]);
    }

protected:
    void reset_timer() override {
        BaseAlgo::reset_timer();
        detail::reset_stats_timers_impl(_linear_solver, 0);
        timer_dSbus_ = 0.;
        timer_fillJ_ = 0.;
        _system.reset_timers();
    }

private:
    static IntVect _to_intvect(const std::vector<int>& v) {
        return Eigen::Map<const IntVect>(v.data(), static_cast<Eigen::Index>(v.size()));
    }

    LinearSolver  _linear_solver;
    DampedSystem  _system;

    bool need_factorize_;

    // Levenberg-Marquardt state/tunables
    real_type lambda_;       // current damping factor, persists across outer NR iterations
    real_type lambda_init_;  // value lambda_ is reset to at the start of each compute_pf call
    real_type nu_;           // growth factor on rejection
    real_type lambda_min_;   // floor applied when easing off after an accepted step
    real_type lambda_max_;   // ceiling
    real_type lambda_eps_;   // additive floor on rejection, so lambda_==0 is not a fixed point
    int       max_trials_;   // per outer-iteration cap on reject-and-retry

    double timer_dSbus_;
    double timer_fillJ_;

    // No copy
    LMNRAlgo(const LMNRAlgo&) = delete;
    LMNRAlgo(LMNRAlgo&&) = delete;
    LMNRAlgo& operator=(LMNRAlgo&&) = delete;
    LMNRAlgo& operator=(const LMNRAlgo&) = delete;
};

#include "LMNRAlgo.tpp"

} // namespace ls2g

#endif // LM_NR_ALGO_H
