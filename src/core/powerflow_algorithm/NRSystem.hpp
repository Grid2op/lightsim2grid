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
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <vector>
#include <stdexcept>

namespace ls2g {

// ---- Extension tag types ------------------------------------------------------

struct MultiSlack {};   // distributed-slack extension (+1 row/col always)
// struct HVDC {};      // future (+n_hvdc_lines rows/cols at runtime)

// ---- Primary template declaration (no definition) -----------------------------

template <typename... Extensions>
class NRSystem;

// ---- Base specialisation: NRSystem<> (single-slack, lag = 0) ------------------

/**
 * Base Newton-Raphson system: single-slack, (n_pvpq + n_pq) x (n_pvpq + n_pq) Jacobian.
 *
 * 3-phase interface:
 *   Phase 1  — init_topology(): builds pvpq maps, sets lag_. Call when pv/pq change.
 *   Phase 1.5— update_state():  updates V/Sbus pointers. Call every compute_pf.
 *   Phase 2  — build_J_sparsity(): symbolic J build + value_map. Call when topology changes.
 *   Phase 3  — fill_J(): fast numerical fill via value_map_. Call each factorisation.
 *
 * Extension variables (slack, future HVDC) go at the LAST rows/cols so the core
 * (n_pvpq + n_pq) x (n_pvpq + n_pq) block is index-compatible for all depths.
 */
template <>
class NRSystem<>
{
public:
    NRSystem() noexcept
        : Ybus_ptr_(nullptr), Sbus_ptr_(nullptr),
          nb_pv_(0), nb_pq_(0), lag_(0),
          need_full_rebuild_(true),
          timer_dSbus_(0.), timer_fillJ_(0.) {}

    virtual ~NRSystem() = default;

    // ----- Phase 1: topology init (call when pv/pq/slack topology changes) -------

    virtual void init_topology(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        const CplxVect&                        Sbus,
        Eigen::Ref<const IntVect>              slack_ids,
        const RealVect&                        slack_weights,
        Eigen::Ref<const IntVect>              pv,
        Eigen::Ref<const IntVect>              pq);

    // ----- Phase 1.5: per-compute_pf state update (cheap) -----------------------

    virtual void update_state(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        const CplxVect&                        V_init,
        const CplxVect&                        Sbus);

    // ----- Phase 2: build J sparsity + value_map (non-virtual) ------------------

    void build_J_sparsity();

    // ----- Phase 3: fill J numerically (non-virtual) ----------------------------

    void fill_J();

    // ----- NR iteration primitives -----------------------------------------------

    virtual RealVect   mismatch()                           const;
    virtual void       apply_step(const RealVect& dx);
    virtual real_type  mismatch_sq_norm_at(const RealVect& dx) const;

    // ----- Housekeeping ----------------------------------------------------------

    void clear_jacobian() {
        J_         = Eigen::SparseMatrix<real_type>();
        dS_dVm_    = Eigen::SparseMatrix<cplx_type>();
        dS_dVa_    = Eigen::SparseMatrix<cplx_type>();
        value_map_.clear();
        need_full_rebuild_ = true;
    }

    const Eigen::SparseMatrix<real_type>& J()  const { return J_; }
    const CplxVect& V()  const { return V_; }
    const RealVect& Va() const { return Va_; }
    const RealVect& Vm() const { return Vm_; }

    // ----- Size / segment accessors (replace NRLayout) ---------------------------

    int nb_pv()      const { return nb_pv_; }
    int nb_pq()      const { return nb_pq_; }
    int lag()        const { return lag_; }
    int theta_size() const { return nb_pv_ + nb_pq_; }
    int vm_size()    const { return nb_pq_; }
    int total()      const { return lag_ + nb_pv_ + 2 * nb_pq_; }

    // Core state vector segments: theta first, then vm, then extension variables.
    RealVect::SegmentReturnType      theta(RealVect& x)       const { return x.segment(0, theta_size()); }
    RealVect::ConstSegmentReturnType theta(const RealVect& x) const { return x.segment(0, theta_size()); }
    RealVect::SegmentReturnType      vm   (RealVect& x)       const { return x.segment(theta_size(), vm_size()); }
    RealVect::ConstSegmentReturnType vm   (const RealVect& x) const { return x.segment(theta_size(), vm_size()); }

    // ----- Timers ----------------------------------------------------------------

    double timer_dSbus() const { return timer_dSbus_; }
    double timer_fillJ() const { return timer_fillJ_; }
    void   reset_timers()      { timer_dSbus_ = 0.; timer_fillJ_ = 0.; }

protected:
    // ---- Shared data (one copy, inherited by all extensions) --------------------
    const Eigen::SparseMatrix<cplx_type>* Ybus_ptr_;
    const CplxVect*                        Sbus_ptr_;
    Eigen::VectorXi                        pv_, pq_, pvpq_;
    std::vector<int>                       pvpq_inv_, pq_inv_;
    int                                    nb_pv_, nb_pq_, lag_;
    RealVect                               Va_, Vm_;
    CplxVect                               V_;
    Eigen::SparseMatrix<real_type>         J_;
    Eigen::SparseMatrix<cplx_type>         dS_dVm_, dS_dVa_;
    std::vector<cplx_type*>                value_map_;
    bool                                   need_full_rebuild_;
    double                                 timer_dSbus_, timer_fillJ_;

    // ---- Virtual hooks (overridden by extensions) --------------------------------

    // Collect ALL structural triplets for J.  Base fills the core block;
    // each extension override calls Base::_collect_J_triplets(coeffs) first then appends.
    virtual void _collect_J_triplets(std::vector<Eigen::Triplet<double>>& coeffs) const;

    // Per-entry value-pointer hook: return pointer into dS_dVa_/dS_dVm_ for
    // variable entries, or nullptr for constant entries (e.g. slack row).
    // Called once per J nonzero during build_J_sparsity — virtual overhead OK here.
    // Base handles the core block; extensions override to handle their own rows.
    virtual cplx_type* _get_entry_ptr(int row, int col);

    // ---- Shared helpers (non-virtual, never overridden) --------------------------

    void _collect_value_map();   // called by build_J_sparsity; iterates J col-major

    void _dSbus_dV(const Eigen::SparseMatrix<cplx_type>& Ybus, const CplxVect& V);
    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
    CplxVect _compute_trial_V(const RealVect& dx) const;
    RealVect _mismatch_core(const CplxVect& V_trial) const;

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       const Eigen::VectorXi& index_col,
                       size_t col_id,
                       size_t row_lag,
                       size_t col_lag) const;

    void _get_values_J(int& nb_obj_this_col,
                       std::vector<Eigen::Index>& inner_index,
                       std::vector<real_type>& values,
                       const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
                       const std::vector<int>& index_row_inv,
                       size_t col_id_mat,
                       size_t row_lag,
                       size_t col_lag) const;

private:
    NRSystem(const NRSystem&)            = delete;
    NRSystem(NRSystem&&)                 = delete;
    NRSystem& operator=(const NRSystem&) = delete;
    NRSystem& operator=(NRSystem&&)      = delete;
};

// ---- Extension specialisation: NRSystem<MultiSlack, Rest...> ------------------

/**
 * Distributed-slack extension.  Inherits all data and logic from NRSystem<Rest...>
 * and adds:
 *   - one extra row/col at the end of J (lag_ += 1)
 *   - slack_absorbed_ tracking in apply_step / mismatch
 *
 * The slack row in J is set once in build_J_sparsity and kept CONSTANT across
 * NR iterations (a valid quasi-Newton approximation for this equation).
 */
template <typename... Rest>
class NRSystem<MultiSlack, Rest...> : public NRSystem<Rest...>
{
    using Base = NRSystem<Rest...>;
public:
    NRSystem() noexcept
        : Base(), slack_bus_id_(0), slack_absorbed_(static_cast<real_type>(0.)) {}

    virtual ~NRSystem() = default;

    // ----- Phase 1 ---------------------------------------------------------------
    virtual void init_topology(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        const CplxVect&                        Sbus,
        Eigen::Ref<const IntVect>              slack_ids,
        const RealVect&                        slack_weights,
        Eigen::Ref<const IntVect>              pv,
        Eigen::Ref<const IntVect>              pq) override;

    // ----- Phase 1.5 -------------------------------------------------------------
    virtual void update_state(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        const CplxVect&                        V_init,
        const CplxVect&                        Sbus) override;

    // ----- NR primitives ---------------------------------------------------------
    virtual RealVect  mismatch()                              const override;
    virtual void      apply_step(const RealVect& dx)                override;
    virtual real_type mismatch_sq_norm_at(const RealVect& dx) const override;

protected:
    // Appends slack row + slack column triplets after the core block
    virtual void _collect_J_triplets(std::vector<Eigen::Triplet<double>>& coeffs) const override;

    // Returns nullptr for the slack row (constant); delegates to Base otherwise
    virtual cplx_type* _get_entry_ptr(int row, int col) override;

private:
    size_t    slack_bus_id_;
    RealVect  slack_weights_;
    real_type slack_absorbed_;

    int _J_slack_row() const { return this->theta_size() + this->vm_size(); }

    void   _append_slack_triplets(std::vector<Eigen::Triplet<double>>& coeffs) const;
    RealVect _mismatch_with_slack(const CplxVect& V_trial, real_type sa) const;
};

// ---- Type aliases (keep existing names working) --------------------------------

using SingleSlackNRSystem = NRSystem<>;
using MultiSlackNRSystem  = NRSystem<MultiSlack>;

} // namespace ls2g

#include "NRSystem.tpp"

#endif // NR_SYSTEM_H
