// Copyright (c) 2026
// Experimental plugin, not part of the lightsim2grid core distribution.
//
// NRAlgoDistSlack: distributed slack implemented as an OUTER loop around a
// single-slack `NRAlgo` (SingleSlackNRSystem = NRSystem<Base, VoltageControl,
// Hvdc>), instead of lightsim2grid's usual in-Jacobian MultiSlack extension.
//
// Rationale (see lightsim2grid session notes / memory
// project_orell_multislack_voltagecontrol_nr_divergence): PowSyBl OpenLoadFlow
// implements distributed slack as a classic outer loop
// (DistributedSlackOuterLoop.java) around a clean single-reference-slack inner
// Newton-Raphson solve -- the MAX_VOLTAGE_CHANGE step-damping policy never sees
// a distributed-slack unknown in the same linearization as VoltageControl's
// Q-sharing/voltage-pinning equations. lightsim2grid's MultiSlackNRSystem bakes
// both into the same Jacobian simultaneously, which was found to diverge under
// MaxVoltageChange damping specifically when a >=2-controller VoltageControl
// sharing group is also entirely composed of distributed-slack participants.
// This class tests whether mimicking OLF's architecture avoids the issue.
//
// Only ACTIVE power is redistributed across outer rounds (mirroring OLF's
// DistributedSlackOuterLoop, which only ever touches targetP) -- reactive
// power / voltage control is handled entirely by the (unchanged) VoltageControl
// extension inside each inner solve.

#ifndef NR_ALGO_DIST_SLACK_H
#define NR_ALGO_DIST_SLACK_H

// Include Solvers.hpp (not powerflow_algorithm/NRAlgo.hpp directly): it pulls
// in SparseLULinearSolver/KLULinearSolver plus, when LS2G_BUILDING_CORE is
// NOT defined (true for any out-of-tree plugin build), `extern template`
// declarations for NRAlgo<SparseLULinearSolver/KLULinearSolver,
// SingleSlackNRSystem> -- so `inner_` below links against the symbols
// already compiled into the main lightsim2grid_cpp library instead of being
// re-instantiated (and re-emitted) inside this plugin's .so.
#include <Solvers.hpp>
#include <powerflow_algorithm/BaseAlgo.hpp>
#include <LSGrid.hpp>

namespace ls2g {

template<class LinearSolver>
class NRAlgoDistSlack final : public BaseAlgo {
public:
    NRAlgoDistSlack() noexcept :
        BaseAlgo(/*is_ac=*/true),
        slack_absorbed_accum_(static_cast<real_type>(0.)),
        outer_iter_(0),
        max_outer_iter_(20),
        outer_tol_(static_cast<real_type>(1e-4)) {}

    ~NRAlgoDistSlack() noexcept override = default;

    void set_lsgrid(const LSGrid* gridmodel) override {
        BaseAlgo::set_lsgrid(gridmodel);
        inner_.set_lsgrid(gridmodel);
    }

    // ----- Jacobian / VoltageControl introspection: forward to the last inner round -----
    Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const override { return inner_.get_J(); }
    IntVect get_theta_to_J_col_python() const override { return inner_.get_theta_to_J_col_python(); }
    IntVect get_vm_to_J_col_python()    const override { return inner_.get_vm_to_J_col_python(); }
    IntVect get_q_to_J_col_python()     const override { return inner_.get_q_to_J_col_python(); }
    IntVect get_p_to_J_row_python()     const override { return inner_.get_p_to_J_row_python(); }
    IntVect get_q_to_J_row_python()     const override { return inner_.get_q_to_J_row_python(); }
    IntVect get_p_buses_python()        const override { return inner_.get_p_buses_python(); }
    IntVect get_p_rows_python()         const override { return inner_.get_p_rows_python(); }
    IntVect get_q_buses_python()        const override { return inner_.get_q_buses_python(); }
    IntVect get_q_rows_python()         const override { return inner_.get_q_rows_python(); }
    IntVect get_theta_buses_python()    const override { return inner_.get_theta_buses_python(); }
    IntVect get_theta_cols_python()     const override { return inner_.get_theta_cols_python(); }
    IntVect get_vm_buses_python()       const override { return inner_.get_vm_buses_python(); }
    IntVect get_vm_cols_python()        const override { return inner_.get_vm_cols_python(); }
    RealVect get_controller_q()         const override { return inner_.get_controller_q(); }
    IntVect  get_controller_kind()      const override { return inner_.get_controller_kind(); }
    IntVect  get_controller_elem_id()   const override { return inner_.get_controller_elem_id(); }

    // No in-Jacobian slack column (that's the whole point) -- but there IS a
    // meaningful total: the accumulated active power redistributed by the outer
    // loop across all its rounds, for THIS compute_pf call.
    int       get_slack_col()      const override { return -1; }
    real_type get_slack_absorbed() const override { return slack_absorbed_accum_; }

    // The inner NRAlgo<LinearSolver, SingleSlackNRSystem> includes both the
    // Hvdc and VoltageControl extensions (see NRAlgo.hpp), forwarded to
    // unconditionally on every outer round -- so this wrapper genuinely
    // supports both, unlike a generic plugin (BaseAlgo's conservative
    // default of `false` for both).
    static constexpr bool SUPPORTS_HVDC_DROOP = true;
    bool supports_hvdc_droop() const noexcept override { return SUPPORTS_HVDC_DROOP; }
    static constexpr bool SUPPORTS_REMOTE_VOLTAGE_CONTROL = true;
    bool supports_remote_voltage_control() const noexcept override { return SUPPORTS_REMOTE_VOLTAGE_CONTROL; }

    // Number of outer redistribution rounds performed by the last compute_pf call.
    int get_outer_iter() const { return outer_iter_; }

    TimerJac get_timers_jacobian() const override { return inner_.get_timers_jacobian(); }

    bool supports_bus_masking() const override { return inner_.supports_bus_masking(); }
    void set_masked_buses(const std::vector<int>& solver_bus_ids) override {
        inner_.set_masked_buses(solver_bus_ids);
    }

    // int_params[0] = max outer rounds; real_params[0] = outer active-power
    // mismatch tolerance (same units as Sbus, i.e. per-unit on Sn). Delegates
    // everything else to the inner NRAlgo's own config (scaling/refactor policy).
    AlgoConfig get_config() const override {
        AlgoConfig cfg = inner_.get_config();
        cfg.int_params.insert(cfg.int_params.begin(), max_outer_iter_);
        cfg.real_params.insert(cfg.real_params.begin(), static_cast<double>(outer_tol_));
        return cfg;
    }
    void set_config(const AlgoConfig& cfg) override {
        if (cfg.int_params.empty() || cfg.real_params.empty())
            throw std::runtime_error("NRAlgoDistSlack::set_config: expects int_params[0]=max_outer_iter, "
                                      "real_params[0]=outer_tol, followed by the inner NRAlgo's own config.");
        max_outer_iter_ = cfg.int_params[0];
        outer_tol_ = static_cast<real_type>(cfg.real_params[0]);
        AlgoConfig inner_cfg;
        inner_cfg.int_params.assign(cfg.int_params.begin() + 1, cfg.int_params.end());
        inner_cfg.real_params.assign(cfg.real_params.begin() + 1, cfg.real_params.end());
        if (!inner_cfg.int_params.empty() || !inner_cfg.real_params.empty()) inner_.set_config(inner_cfg);
    }

    void reset() override {
        BaseAlgo::reset();
        inner_.reset();
        slack_absorbed_accum_ = static_cast<real_type>(0.);
        outer_iter_ = 0;
    }

    bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                    CplxVect & V,
                    const Eigen::Ref<const CplxVect> & Sbus,
                    const Eigen::Ref<const IntVect> & slack_ids,
                    const Eigen::Ref<const RealVect> & slack_weights,
                    const Eigen::Ref<const IntVect> & pv,
                    const Eigen::Ref<const IntVect> & pq,
                    int max_iter,
                    real_type tol) override
    {
        err_ = ErrorType::NoError;
        slack_absorbed_accum_ = static_cast<real_type>(0.);
        outer_iter_ = 0;
        nr_iter_ = 0;

        // Every distributed-slack participant except the reference (slack_ids(0))
        // becomes an ordinary PV bus for the inner single-slack solve; Base's own
        // `free_vm_slack_buses_` mechanism (driven by LSGrid's persistent slack
        // participant set, independent of the slack_ids we pass here) still
        // correctly gives any non-locally-pinned bus among them a free Vm unknown
        // + Q equation, so remote voltage control (incl. multi-controller sharing
        // groups) keeps working exactly as it does for the multi-slack case.
        IntVect my_pv = retrieve_pv_with_slack(slack_ids, pv);
        IntVect ref_slack_id(1);
        ref_slack_id(0) = slack_ids(0);

        CplxVect Sbus_working = Sbus;
        CplxVect V_round = V;

        AlgoControl round_control = _solver_control;

        const int n_participants = static_cast<int>(slack_ids.size());
        const real_type w_ref = slack_weights(slack_ids(0));
        const real_type one_minus_w_ref = static_cast<real_type>(1.) - w_ref;

        bool conv = false;
        for (outer_iter_ = 0; outer_iter_ < max_outer_iter_; ++outer_iter_) {
            inner_.tell_solver_control(round_control);
            conv = inner_.compute_pf(Ybus, V_round, Sbus_working, ref_slack_id, slack_weights,
                                      my_pv, pq, max_iter, tol);
            nr_iter_ += inner_.get_nb_iter();
            if (!conv) { err_ = inner_.get_error(); break; }

            // Subsequent inner rounds only ever change Sbus (redistributed P):
            // reuse the symbolic AND numeric factorization / Jacobian sparsity.
            round_control.tell_none_changed();
            V_round = inner_.get_V();

            // Full-vector mismatch (same formula as LSGrid::compute_results),
            // we only need the reference bus' component.
            const CplxVect mismatch = V_round.array() * (Ybus * V_round).conjugate().array() - Sbus_working.array();
            const real_type mis = std::real(mismatch(slack_ids(0)));

            if (std::abs(mis) < outer_tol_) break;

            if (n_participants > 1 && one_minus_w_ref > BaseConstants::_tol_equal_float) {
                for (int k = 1; k < n_participants; ++k) {
                    const int bus = slack_ids(k);
                    const real_type w = slack_weights(bus);
                    if (w <= BaseConstants::my_zero_) continue;
                    Sbus_working(bus) += cplx_type(w / one_minus_w_ref * mis, static_cast<real_type>(0.));
                }
            }
            slack_absorbed_accum_ += mis;
        }

        V_ = inner_.get_V();
        Va_ = inner_.get_Va();
        Vm_ = inner_.get_Vm();
        n_ = static_cast<int>(V_.size());
        return conv;
    }

private:
    NRAlgo<LinearSolver, SingleSlackNRSystem> inner_;
    real_type slack_absorbed_accum_;
    int outer_iter_;
    int max_outer_iter_;
    real_type outer_tol_;

    NRAlgoDistSlack(const NRAlgoDistSlack&) = delete;
    NRAlgoDistSlack(NRAlgoDistSlack&&) = delete;
    NRAlgoDistSlack& operator=(NRAlgoDistSlack&&) = delete;
    NRAlgoDistSlack& operator=(const NRAlgoDistSlack&) = delete;
};

} // namespace ls2g

#endif // NR_ALGO_DIST_SLACK_H
