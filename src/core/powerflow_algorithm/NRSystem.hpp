// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_SYSTEM_H
#define NR_SYSTEM_H

#include "Utils.hpp"
#include "BaseConstants.hpp"
#include "CustTimer.hpp"
#include "NRLedger.hpp"
#include "HvdcDroopData.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <stdexcept>

// Public API (used by NRAlgo, in order)
// init_topology : if Ybus changed
// update_state : unconditionnally
// build_J_sparsity : if init_topology was called
// mismatch (before NR loop)
// ---------------- ENTERING NR LOOP
// if need factorize:
//     fill_internal_variables
//     fill_J
// apply_step
// mismatch
// ---------------- END NR LOOP
// V(), Vm(), Va()
// timer_dSbus()
// timer_fillJ()
// J()  (public accessor of NRAlgo)

// Public API (bus id -> Jacobian column converters, built in init_topology)
// theta_to_J_col() : bus id -> column of that bus' theta (ΔVa) unknown, -1 if none
// vm_to_J_col()    : bus id -> column of that bus' vm (ΔVm) unknown, -1 if none
// q_to_J_col()     : bus id -> column of that bus' q unknown, -1 if none

// Public API (used by scaling policy)
// mismatch_sq_norm_at()
// max_abs_dtheta()
// max_abs_dvm()


namespace ls2g {

class LSGrid;  // only a pointer travels through the component protocol

// ---- Primary template declaration (no definition) -----------------------------

template <typename... Extensions>
class NRSystem;

enum LS2G_API dSdX { dSdVa_r, dSdVa_i, dSdVm_r, dSdVm_i };

// One dS-derived Jacobian entry: J(jrow, jcol) takes the real or imaginary
// part (whichmat) of the dS_dVa / dS_dVm value stored at Ybus nnz position
// ybus_k. Produced by the generic dS pass of NRSystem::build_J_sparsity.
class LS2G_API Contrib final
{
    private:
        int jrow_;
        int jcol_;
        int ybus_k_;
        dSdX whichmat_;

    public:
       Contrib(int jrow, int jcol, int ybus_k, dSdX whichmat) noexcept:
           jrow_(jrow),
           jcol_(jcol),
           ybus_k_(ybus_k),
           whichmat_(whichmat) {}

       int jrow() const {return jrow_;}
       int jcol() const {return jcol_;}
       int ybus_k() const {return ybus_k_;}
       dSdX whichmat() const {return whichmat_;}

};

// ---- Component protocol --------------------------------------------------------
//
// The Base block and every extension implement the same interface; NRSystem
// composes them (Base is the first, special-cased member; extensions live in a
// tuple). A component NEVER manipulates global J offsets: it claims its rows /
// columns in the NRLedger and keeps the returned indices.
//
//  - update_state(...)            : per-solve state refresh (cheap, always).
//  - init_topology(...)           : computes the component's own index sets only.
//  - register_in(NRLedger&)       : claims bus-owned rows/cols and allocates
//                                   custom rows/cols; stores the returned indices.
//  - declare_feature_entries(FeatureSink&) : declares the component's non-dS
//                                   Jacobian entries; stores the handles.
//  - fill_feature_values(FeatureWriter&, Va) : writes the feature values
//                                   (called at every fill_J; J was zeroed first).
//                                   Va is the CURRENT voltage-angle vector (all
//                                   solver buses), for state-dependent slopes.
//  - adjust_mismatch(V_t, dx, mis): adds the component's physical injections to
//                                   the per-bus complex mismatch, evaluated at
//                                   the trial voltages V_t. dx is the FULL
//                                   step vector (zero vector of full size when
//                                   evaluating at the current state); the
//                                   component indexes it with its stored column ids.
//  - fill_custom_rows(res, Va, Vm, dx) : accumulates (+=) residuals ONLY into
//                                   custom rows the component allocated. Bus-owned
//                                   P/Q rows are filled generically by NRSystem.
//  - apply_step(dx)               : updates ONLY the component's own non-voltage
//                                   state. Voltage updates are generic.
//  - clear()                      : reset.

/**
 * Base Newton-Raphson system block (single-slack core).
 *
 * Jacobian orientation (applies to the whole augmented J):
 *   - each ROW is an EQUATION   : a power-mismatch residual
 *       (ΔP at a pv/pq bus, ΔQ at a pq bus)
 *   - each COLUMN is an UNKNOWN  : a state-variable correction
 *       (ΔVa at a pv/pq bus, ΔVm at a pq bus, ... extensions add more)
 *
 * i.e. J[row, col] = d(equation row) / d(unknown col), used to solve
 *   J . dx = -F(x). So there is one row per equation and one column per
 *   unknown (NOT one row per unknown).
 *
 * The base block registers, in INCREASING BUS INDEX order (the pv / pq input
 * arrays are NOT assumed sorted; they are sorted internally):
 *   - theta unknowns for sorted(pv ∪ pq),
 *   - vm unknowns for sorted(pq),
 *   - P equations for sorted(pv ∪ pq),
 *   - Q equations for sorted(pq).
 */
class LS2G_API Base
{
    public:
        Base():
            nb_pv_(0),
            nb_pq_(0),
            nb_pvpq_(0)
            {}

        static int find_J_pos (
            Eigen::Ref<const Eigen::SparseMatrix<real_type, Eigen::ColMajor> > J_csc,
            int row,
            int col){
            int start = J_csc.outerIndexPtr()[col];
            int end   = J_csc.outerIndexPtr()[col + 1];
            const int* inner = J_csc.innerIndexPtr();
            auto it = std::lower_bound(inner + start, inner + end, row);
            if (it == inner + end || *it != row) return -1;
            return (int) (it - inner);
        };

        // call at the beginning of each solve
        void update_state(
            const LSGrid *                         /*lsgrid_ptr*/,
            const Eigen::SparseMatrix<cplx_type>&  /*Ybus*/,
            const CplxVect&                        /*Sbus*/,
            const RealVect&                        /*slack_weights*/
        ){}

        // call after update_state
        // at the beginning of each solve
        // only if the topology has changed
        void init_topology(
            Eigen::Ref<const IntVect>              /*slack_ids*/,
            const RealVect&                        /*slack_weights*/,
            Eigen::Ref<const IntVect>              pv,
            Eigen::Ref<const IntVect>              pq
        ) {
            pv_ = IntVect(pv);
            pq_ = IntVect(pq);

            nb_pv_ = static_cast<int>(pv_.size());
            nb_pq_ = static_cast<int>(pq_.size());

            pvpq_.resize(nb_pv_ + nb_pq_);
            pvpq_ << pv_, pq_;

            nb_pvpq_ = static_cast<int>(pvpq_.size());
        }

        void register_in(NRLedger& ledger) const
        {
            std::vector<int> pvpq_sorted(pvpq_.data(), pvpq_.data() + nb_pvpq_);
            std::sort(pvpq_sorted.begin(), pvpq_sorted.end());
            std::vector<int> pq_sorted(pq_.data(), pq_.data() + nb_pq_);
            std::sort(pq_sorted.begin(), pq_sorted.end());

            for (int bus : pvpq_sorted) ledger.add_theta_unknown(bus);
            for (int bus : pq_sorted)   ledger.add_vm_unknown(bus);
            for (int bus : pvpq_sorted) ledger.add_p_equation(bus);
            for (int bus : pq_sorted)   ledger.add_q_equation(bus);
        }

        void declare_feature_entries(FeatureSink& /*sink*/) {}
        void fill_feature_values(FeatureWriter& /*writer*/, const RealVect& /*Va*/) const {}
        void adjust_mismatch(const CplxVect& /*V_t*/, const RealVect& /*dx*/, CplxVect& /*mis*/) const {}
        void fill_custom_rows(RealVect& /*res*/,
                              const RealVect& /*Va*/,
                              const RealVect& /*Vm*/,
                              const RealVect& /*dx*/) const {}
        void apply_step(const RealVect& /*dx*/) {}

        void clear() {
            pv_ = IntVect();
            pq_ = IntVect();
            pvpq_ = IntVect();
            nb_pv_ = 0;
            nb_pq_ = 0;
            nb_pvpq_ = 0;
        }

    private:
        IntVect pv_, pq_, pvpq_;
        int     nb_pv_, nb_pq_, nb_pvpq_;

    public:
        Eigen::Ref<const IntVect> pv() const { return pv_; }
        Eigen::Ref<const IntVect> pq() const { return pq_; }
        Eigen::Ref<const IntVect> pvpq() const { return pvpq_; }
        int                       nb_pv() const {return nb_pv_;}
        int                       nb_pq() const {return nb_pq_;}
        int                       nb_pvpq() const {return nb_pvpq_;}
};

// ---- Extension tag types ------------------------------------------------------

/**
 * This extension adds the ability to have a distributed slack directly in the jacobian.
 *
 * It adds exactly "nb slack" extra rows and columns (nb_slack being the
 * number of slacks in the grid). It registers, after the base block:
 *   - for each non-ref slack bus (sorted by bus id): a theta unknown and a
 *     P equation;
 *   - a P equation for the ref slack bus (slack_ids[0]);
 *   - one custom column: the slack_absorbed unknown.
 *
 * All the dS-derived Jacobian entries of those rows / columns are generated by
 * the generic dS pass of NRSystem. The only feature-specific entries are the
 * slack_absorbed column coefficients: slack_weights(bus) at each slack bus'
 * P row.
 */
class LS2G_API MultiSlack   // distributed-slack extension
{
    public:
        MultiSlack():
            my_size_(0),
            ref_slack_id_(0),
            slack_col_(-1),
            slack_absorbed_(static_cast<real_type>(0.)) {}

        // call at the beginning of each NR solve
        // once. It is not called for
        // NR iterations
        void update_state(
            const Base *                           /*nr_system_base_ptr*/,
            const LSGrid *                         /*lsgrid_ptr*/,
            const Eigen::SparseMatrix<cplx_type>&  /*Ybus*/,
            const CplxVect&                        Sbus,
            const RealVect&                        slack_weights
        ){
            slack_weights_ = slack_weights;
            // initial slack absorbed (see MultiSlackPolicy::initial_slack_absorbed)
            slack_absorbed_ = std::real(Sbus.sum());
        }

        // call after update_state
        // at the beginning of each solve
        // only if the topology has changed
        void init_topology(
            Eigen::Ref<const IntVect>              slack_ids,
            const RealVect&                        /*slack_weights*/,
            Eigen::Ref<const IntVect>              /*pv*/,
            Eigen::Ref<const IntVect>              /*pq*/
        ) {
            my_size_ = static_cast<int>(slack_ids.size());
            ref_slack_id_ = slack_ids[0];
            // slack buses in registration order: non-ref slacks sorted by bus
            // id, then the ref slack last
            slack_buses_.assign(slack_ids.data() + 1, slack_ids.data() + my_size_);
            std::sort(slack_buses_.begin(), slack_buses_.end());
            slack_buses_.push_back(ref_slack_id_);
        }

        void register_in(NRLedger& ledger)
        {
            slack_p_rows_.clear();
            slack_p_rows_.reserve(my_size_);
            for (int k = 0; k + 1 < my_size_; ++k) {
                ledger.add_theta_unknown(slack_buses_[k]);
                slack_p_rows_.push_back(ledger.add_p_equation(slack_buses_[k]));
            }
            slack_p_rows_.push_back(ledger.add_p_equation(ref_slack_id_));
            slack_col_ = ledger.add_custom_col();  // slack_absorbed unknown
        }

        // one entry per slack bus: (that bus' P row, slack_absorbed column)
        void declare_feature_entries(FeatureSink& sink)
        {
            feature_handles_.clear();
            feature_handles_.reserve(my_size_);
            for (int k = 0; k < my_size_; ++k)
                feature_handles_.push_back(sink.add(slack_p_rows_[k], slack_col_));
        }

        void fill_feature_values(FeatureWriter& writer, const RealVect& /*Va*/) const
        {
            for (int k = 0; k < my_size_; ++k)
                writer.add(feature_handles_[k], slack_weights_(slack_buses_[k]));
        }

        // adjust the per-bus complex power mismatch by (slack_absorbed + dx step) * slack_weights
        void adjust_mismatch(const CplxVect& /*V_t*/, const RealVect& dx, CplxVect& mis) const
        {
            const real_type sa = slack_absorbed_ + dx(slack_col_);
            mis.array() += (sa * slack_weights_.array()).cast<cplx_type>();
        }

        // all this extension's rows are bus-owned P equations,
        // filled generically by NRSystem
        void fill_custom_rows(RealVect& /*res*/,
                              const RealVect& /*Va*/,
                              const RealVect& /*Vm*/,
                              const RealVect& /*dx*/) const {}

        // voltage updates (the non-ref slack thetas) are generic; only the
        // slack_absorbed state is owned by this extension
        void apply_step(const RealVect& dx)
        {
            slack_absorbed_ += dx(slack_col_);
        }

        void clear(){
            my_size_ = 0;
            ref_slack_id_ = 0;
            slack_col_ = -1;
            slack_buses_.clear();
            slack_p_rows_.clear();
            feature_handles_.clear();
            slack_weights_ = RealVect();
            slack_absorbed_ = static_cast<real_type>(0.);
        }

    private:
        int                my_size_;
        int                ref_slack_id_;
        int                slack_col_;        // J column of the slack_absorbed unknown
        std::vector<int>   slack_buses_;      // non-ref slacks sorted, then ref slack
        std::vector<int>   slack_p_rows_;     // P row of each slack bus (same order)
        std::vector<int>   feature_handles_;  // FeatureSink handles (same order)
        RealVect           slack_weights_;    // size: nb_bus
        real_type          slack_absorbed_;

};

/**
 * Hvdc angle-droop ("AC emulation") extension.
 *
 * For every CONNECTED droop-enabled hvdc line, the active power follows
 *   raw = p0 + k * (theta1 - theta2)
 * instead of a fixed setpoint (the HvdcLineContainer does NOT stamp those
 * lines in Sbus in AC). The flows / derivatives mirror pyloadflow
 * `equations/ac_emulation.py` (itself a port of open-loadflow's
 * HvdcAcEmulationSide{1,2}ActiveFlowEquationTerm).
 *
 * The extension claims NO row / column: its equations are the existing bus
 * P mismatches and its unknowns the existing bus thetas. It only:
 *   - declares up to 4 dP/dtheta Jacobian entries per line,
 *     (p_row(b1) | p_row(b2)) x (theta_col(b1) | theta_col(b2)), where the
 *     row / column exists (a slack end simply drops its entries: its theta is
 *     not an unknown / its P balance not an equation of the base block).
 *     They are declared whatever the droop regime, so a `status_droop` flip
 *     between two solves does NOT change the J sparsity pattern and the
 *     symbolic factorization of the linear solver is reused;
 *   - adds the theta-dependent flows to the per-bus mismatch;
 *   - writes the (piecewise constant) droop slopes: on the controller side
 *     +/- k, on the non-controller side +/- k * (1 - lf1) * (1 - lf2), the
 *     derivative of the resistive dc-line loss being neglected (OLF parity).
 *     Saturated lines are pure constant injections: zero slopes.
 *
 * When the grid has no droop hvdc line, every loop below is empty: the
 * behaviour (and the J pattern) is strictly identical to a build without
 * this extension.
 */
class LS2G_API Hvdc
{
    public:
        Hvdc(): my_size_(0) {}

        // pulls the droop data (solver bus ids, pu) from the grid;
        // defined in NRSystemHvdc.cpp (needs the full LSGrid type)
        void update_state(
            const Base *                           nr_system_base_ptr,
            const LSGrid *                         lsgrid_ptr,
            const Eigen::SparseMatrix<cplx_type>&  Ybus,
            const CplxVect&                        Sbus,
            const RealVect&                        slack_weights
        );

        void init_topology(
            Eigen::Ref<const IntVect>              /*slack_ids*/,
            const RealVect&                        /*slack_weights*/,
            Eigen::Ref<const IntVect>              /*pv*/,
            Eigen::Ref<const IntVect>              /*pq*/
        ) {}

        // claims nothing: only caches the rows / columns of the two end buses
        // (the ledger is fully populated by the previous components: Hvdc must
        // be the LAST extension of the tuple)
        void register_in(NRLedger& ledger)
        {
            p_row_1_.resize(my_size_);
            p_row_2_.resize(my_size_);
            theta_col_1_.resize(my_size_);
            theta_col_2_.resize(my_size_);
            for (int k = 0; k < my_size_; ++k) {
                p_row_1_[k] = ledger.p_row(data_.bus1(k));
                p_row_2_[k] = ledger.p_row(data_.bus2(k));
                theta_col_1_[k] = ledger.theta_col(data_.bus1(k));
                theta_col_2_[k] = ledger.theta_col(data_.bus2(k));
            }
        }

        // up to 4 dP/dtheta entries per line, declared whatever the regime
        void declare_feature_entries(FeatureSink& sink)
        {
            h11_.assign(my_size_, -1);
            h12_.assign(my_size_, -1);
            h21_.assign(my_size_, -1);
            h22_.assign(my_size_, -1);
            for (int k = 0; k < my_size_; ++k) {
                if ((p_row_1_[k] >= 0) && (theta_col_1_[k] >= 0)) h11_[k] = sink.add(p_row_1_[k], theta_col_1_[k]);
                if ((p_row_1_[k] >= 0) && (theta_col_2_[k] >= 0)) h12_[k] = sink.add(p_row_1_[k], theta_col_2_[k]);
                if ((p_row_2_[k] >= 0) && (theta_col_1_[k] >= 0)) h21_[k] = sink.add(p_row_2_[k], theta_col_1_[k]);
                if ((p_row_2_[k] >= 0) && (theta_col_2_[k] >= 0)) h22_[k] = sink.add(p_row_2_[k], theta_col_2_[k]);
            }
        }

        void fill_feature_values(FeatureWriter& writer, const RealVect& Va) const
        {
            for (int k = 0; k < my_size_; ++k) {
                if (data_.status(k) != 0) continue;  // saturated: constant injection, zero slopes
                const real_type raw = data_.p0(k) + data_.k(k) * (Va(data_.bus1(k)) - Va(data_.bus2(k)));
                const real_type loss_mult = (1. - data_.lf1(k)) * (1. - data_.lf2(k));
                // dp1 = dp1/dtheta1, dp2 = dp2/dtheta1; d/dtheta2 = -d/dtheta1
                const real_type dp1 = (raw >= 0.) ? data_.k(k) : data_.k(k) * loss_mult;
                const real_type dp2 = (raw < 0.) ? -data_.k(k) : -data_.k(k) * loss_mult;
                if (h11_[k] >= 0) writer.add(h11_[k], dp1);
                if (h12_[k] >= 0) writer.add(h12_[k], -dp1);
                if (h21_[k] >= 0) writer.add(h21_[k], dp2);
                if (h22_[k] >= 0) writer.add(h22_[k], -dp2);
            }
        }

        // the flows LEAVE the buses into the hvdc: they ADD to the computed
        // power, ie to the mismatch (mis = V conj(Ybus V) - Sbus + ...)
        void adjust_mismatch(const CplxVect& V_t, const RealVect& /*dx*/, CplxVect& mis) const
        {
            real_type p1_flow, p2_flow;
            for (int k = 0; k < my_size_; ++k) {
                const real_type theta_1 = std::arg(V_t(data_.bus1(k)));
                const real_type theta_2 = std::arg(V_t(data_.bus2(k)));
                const real_type raw = data_.p0(k) + data_.k(k) * (theta_1 - theta_2);
                data_.flows_pu(k, raw, p1_flow, p2_flow);
                mis(data_.bus1(k)) += p1_flow;
                mis(data_.bus2(k)) += p2_flow;
            }
        }

        void fill_custom_rows(RealVect& /*res*/,
                              const RealVect& /*Va*/,
                              const RealVect& /*Vm*/,
                              const RealVect& /*dx*/) const {}

        void apply_step(const RealVect& /*dx*/) {}

        void clear() {
            my_size_ = 0;
            data_.clear();
            p_row_1_.clear();
            p_row_2_.clear();
            theta_col_1_.clear();
            theta_col_2_.clear();
            h11_.clear();
            h12_.clear();
            h21_.clear();
            h22_.clear();
        }

    private:
        int                  my_size_;
        HvdcDroopSolverData  data_;
        std::vector<int>     p_row_1_, p_row_2_;          // P row of each end bus (-1 if none)
        std::vector<int>     theta_col_1_, theta_col_2_;  // theta column of each end bus (-1 if none)
        std::vector<int>     h11_, h12_, h21_, h22_;      // FeatureSink handles (-1 if entry dropped)
};


// ---- NRSystem<Base, Rest...> ----------------------------------------------------

/**
 * Composable Newton-Raphson system: a Base block (single-slack core) plus any
 * number of extensions (e.g. MultiSlack), all registered in a central NRLedger.
 *
 * 3-phase interface:
 *   Phase 1  — init_topology(): components compute their index sets, then claim
 *              their equations (J rows) / unknowns (J columns) in the ledger.
 *              Call when pv/pq/slack topology changes.
 *   Phase 1.5— update_state():  updates V/Sbus pointers. Call every compute_pf.
 *   Phase 2  — build_J_sparsity(): symbolic J build + value maps. Call when
 *              topology changes.
 *   Phase 3  — fill_J(): fast numerical fill via the value maps. Call each
 *              factorisation.
 *
 * All the dS-derived Jacobian entries (∂(P or Q mismatch at bus i)/∂(θ or Vm
 * at bus j), which exist iff Ybus(i, j) != 0) are generated by ONE generic
 * pass over the Ybus nonzeros, for every component at once: the ledger says
 * which buses own a P/Q equation and a θ/Vm unknown. Components only declare
 * their feature-specific (non dS-derived) entries through the FeatureSink.
 *
 * Likewise the residual P/Q rows, the voltage updates (apply_step /
 * _compute_trial_V) and the scaling reductions (max_abs_dtheta / max_abs_dvm)
 * are generic loops over the ledger's (bus, row) / (bus, col) pair lists;
 * components only handle their own custom rows and non-voltage state.
 *
 * Layout: each component registers its bus-owned rows / columns in INCREASING
 * BUS INDEX order, components in declaration order (Base first). Several
 * contributions may resolve to the same J position, so fill_J zeroes J then
 * accumulates (+=).
 */
template <typename... Rest>
class NRSystem<Base, Rest...>
{
public:
    NRSystem() noexcept:
        timer_dSbus_(0.),
        timer_fillJ_(0.),
        lsgrid_ptr_(nullptr),
        Ybus_ptr_(nullptr),
        Sbus_ptr_(nullptr) {}

    virtual ~NRSystem() = default;

    // ----- Phase 1: topology init (call when pv/pq/slack topology changes) -------

    void init_topology(
        Eigen::Ref<const IntVect>              slack_ids,
        const RealVect&                        slack_weights,
        Eigen::Ref<const IntVect>              pv,
        Eigen::Ref<const IntVect>              pq);

    // ----- Phase 1.5: per-compute_pf state update (cheap) -----------------------

    void update_state(
        const LSGrid *                         lsgrid_ptr,
        const Eigen::SparseMatrix<cplx_type>&  Ybus,
        const CplxVect&                        V_init,
        const CplxVect&                        Sbus,
        Eigen::Ref<const RealVect>             slack_weights);

    // ----- Phase 2: build J sparsity + value maps -------------------------------

    void build_J_sparsity();

    // ----- Phase 3: fill J numerically -------------------------------------------

    void fill_J();
    void fill_internal_variables();

    // ----- NR iteration primitives -----------------------------------------------

    virtual RealVect   mismatch()                           const;
    virtual void       apply_step(const RealVect& dx);
    virtual real_type  mismatch_sq_norm_at(const RealVect& dx) const;

    // ----- Housekeeping ----------------------------------------------------------

    void clear_jacobian() {
        J_         = Eigen::SparseMatrix<real_type, Eigen::ColMajor>();
        base_.clear();
        _clear_extensions(std::make_index_sequence<sizeof...(Rest)>{});

        dS_dVm_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();
        dS_dVa_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();

        map_dsdva_r_.clear();
        map_dsdva_i_.clear();
        map_dsdvm_r_.clear();
        map_dsdvm_i_.clear();
        sink_.clear();
        feature_pos_.clear();
        ledger_.reset(0);
    }

    Eigen::Ref<const Eigen::SparseMatrix<real_type> > J()  const { return J_; }
    Eigen::Ref<const CplxVect> V()  const { return V_; }
    Eigen::Ref<const RealVect> Va() const { return Va_; }
    Eigen::Ref<const RealVect> Vm() const { return Vm_; }

    // bus_id -> Jacobian column of that bus' theta / vm / q unknown (-1 if none).
    // Each vector has size n_bus and spans the full augmented J (base + extensions).
    const std::vector<int>& theta_to_J_col() const { return ledger_.theta_col_of_bus(); }
    const std::vector<int>& vm_to_J_col()    const { return ledger_.vm_col_of_bus(); }
    const std::vector<int>& q_to_J_col()     const { return ledger_.q_col_of_bus(); }

    size_t total_state_variables() const { return static_cast<size_t>(ledger_.size()); }

    // ----- Scaling reductions ----------------------------------------------------
    // max |angle step| / max |voltage-magnitude step| across all state variables.
    real_type max_abs_dtheta(const RealVect& dx) const {
        real_type m = static_cast<real_type>(0.);
        for (int col : ledger_.theta_cols()) m = std::max(m, std::abs(dx(col)));
        return m;
    }
    real_type max_abs_dvm(const RealVect& dx) const {
        real_type m = static_cast<real_type>(0.);
        for (int col : ledger_.vm_cols()) m = std::max(m, std::abs(dx(col)));
        return m;
    }

    // ----- Timers ----------------------------------------------------------------

    double timer_dSbus() const { return timer_dSbus_; }
    double timer_fillJ() const { return timer_fillJ_; }
    void   reset_timers()      { timer_dSbus_ = 0.; timer_fillJ_ = 0.; }

private:
    // ---- Shared data (one copy, shared by all components) -----------------------
    RealVect                               Va_, Vm_;
    CplxVect                               V_;
    Eigen::SparseMatrix<real_type, Eigen::ColMajor>         J_;
    double                                 timer_dSbus_, timer_fillJ_;
    Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>         dS_dVm_, dS_dVa_;

    // dS value maps: Ybus nnz position -> position in J_.valuePtr() (-1 if unused)
    std::vector<int>                       map_dsdva_r_;
    std::vector<int>                       map_dsdva_i_;
    std::vector<int>                       map_dsdvm_r_;
    std::vector<int>                       map_dsdvm_i_;

    // central registry of equations / unknowns (defines the J layout)
    NRLedger                               ledger_;
    // feature (non dS-derived) entries: (row, col) list + resolved valuePtr positions
    FeatureSink                            sink_;
    std::vector<int>                       feature_pos_;

    // Holds the base things
    Base                                   base_;

    // Holds the state for HVDC, DistSlack, etc.
    std::tuple<Rest...>                    extensions_;

protected:
    // visible attribute for derived class (non owning ptr)
    const LSGrid *                                         lsgrid_ptr_;
    const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>* Ybus_ptr_;
    const CplxVect*                                        Sbus_ptr_;

    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
    CplxVect _compute_trial_V(const RealVect& dx) const;
    // assemble the (negated) residual at trial voltages V_t; dx is the step that
    // produced V_t (used by components that carry extra state, e.g. slack absorbed).
    RealVect _residual(const CplxVect& V_t, const RealVect& dx) const;

private:
    // ---- component hook fold helpers (C++14 index_sequence/dummy-array idiom) ---
    template <std::size_t... Is>
    void _init_topology_extensions(
        Eigen::Ref<const IntVect>              slack_ids,
        const RealVect&                        slack_weights,
        Eigen::Ref<const IntVect>              pv,
        Eigen::Ref<const IntVect>              pq,
        std::index_sequence<Is...>) {
        int dummy[] = { 0, (std::get<Is>(extensions_).init_topology(
            slack_ids,
            slack_weights,
            pv,
            pq
            ), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _update_state_extensions(
        const LSGrid *                         lsgrid_ptr,
        const Eigen::SparseMatrix<cplx_type>&  Ybus,
        const CplxVect&                        Sbus,
        const RealVect&                        slack_weights,
        std::index_sequence<Is...>){
        int dummy[] = { 0, (std::get<Is>(extensions_).update_state(
            &base_,
            lsgrid_ptr,
            Ybus,
            Sbus,
            slack_weights
            ), 0)... };
        (void)dummy;
        // silence unused-parameter warnings when the extension pack is empty
        (void)lsgrid_ptr; (void)Ybus; (void)Sbus; (void)slack_weights;
    }

    template <std::size_t... Is>
    void _register_in_extensions(NRLedger& ledger, std::index_sequence<Is...>) {
        int dummy[] = { 0, (std::get<Is>(extensions_).register_in(ledger), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _declare_feature_entries_extensions(FeatureSink& sink, std::index_sequence<Is...>) {
        int dummy[] = { 0, (std::get<Is>(extensions_).declare_feature_entries(sink), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _fill_feature_values_extensions(FeatureWriter& writer, std::index_sequence<Is...>) const {
        int dummy[] = { 0, (std::get<Is>(extensions_).fill_feature_values(writer, Va_), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _adjust_mismatch_extensions(const CplxVect& V_t, const RealVect& dx, CplxVect& mis,
                                     std::index_sequence<Is...>) const {
        int dummy[] = { 0, (std::get<Is>(extensions_).adjust_mismatch(V_t, dx, mis), 0)... };
        (void)dummy;
        (void)V_t;
    }

    template <std::size_t... Is>
    void _fill_custom_rows_extensions(RealVect& res,
                                      const RealVect& Va, const RealVect& Vm,
                                      const RealVect& dx,
                                      std::index_sequence<Is...>) const {
        int dummy[] = { 0, (std::get<Is>(extensions_).fill_custom_rows(res, Va, Vm, dx), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _apply_step_extensions(const RealVect& dx, std::index_sequence<Is...>) {
        int dummy[] = { 0, (std::get<Is>(extensions_).apply_step(dx), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _clear_extensions(std::index_sequence<Is...>) {
        int dummy[] = { 0, (std::get<Is>(extensions_).clear(), 0)... };
        (void)dummy;
    }

private:
    NRSystem(const NRSystem&)            = delete;
    NRSystem(NRSystem&&)                 = delete;
    NRSystem& operator=(const NRSystem&) = delete;
    NRSystem& operator=(NRSystem&&)      = delete;
};

// ---- Type aliases (keep existing names working) --------------------------------
// Hvdc must be the LAST extension (it reads the ledger populated by the others)

using SingleSlackNRSystem = NRSystem<Base, Hvdc>;
using MultiSlackNRSystem  = NRSystem<Base, MultiSlack, Hvdc>;

} // namespace ls2g

#include "NRSystem.tpp"

#endif // NR_SYSTEM_H
