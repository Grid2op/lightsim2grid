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

// TODO need to be cleaned up and simplified

// Public API (used by NRAlgo, in order)
// init_topology : if Ybus changed
// update_state : unconditionnally
// build_J_sparsity : if init_topology was called
// mismatch (before NR loop)
// ---------------- ENTERING NR LOOP
// if need factorize:
//     fill_J
// apply_step
// mismatch
// ---------------- END NRR LOOP 
// V(), Vm(), Va()
// timer_dSbus()
// timer_fillJ()
// J()  (public accessor of NRAlgo)

// Public API (used by scaling policy)
// mismatch_sq_norm_at()
// theta()
// vm()
// theta_size()
// vm_size()


namespace ls2g {

// ---- Primary template declaration (no definition) -----------------------------

template <typename... Extensions>
class NRSystem;


// ---- Extension tag types ------------------------------------------------------

class MultiSlack   // distributed-slack extension
{
    public:
        template<class NRSystemCLS>
        void update_state(
            const NRSystemCLS *                    nr_system_ptr,
            const LSGrid *                         lsgrid_ptr,
            const Eigen::SparseMatrix<cplx_type>&  Ybus,
            const CplxVect&                        Sbus,
            const RealVect&                        slack_weights
        ){
            // TODO remember slack weights !
        }

        void init_topology(
            Eigen::Ref<const IntVect>              slack_ids,
            const RealVect&                        slack_weights,
            Eigen::Ref<const IntVect>              pv,
            Eigen::Ref<const IntVect>              pq
        ) {
            my_size_ = slack_ids.size();
        }

    private:
        size_t my_size_;

};

// struct HVDC {};      // future (+n_hvdc_lines rows/cols at runtime)


// ---- Base specialisation: NRSystem<> (single-slack) ------------------

/**
 * Base Newton-Raphson system: single-slack, (n_pvpq + n_pq) x (n_pvpq + n_pq) Jacobian.
 *
 * 3-phase interface : TODO refacto NR description
 *   Phase 1  — init_topology(): builds pvpq maps, sets lag_. Call when pv/pq change.
 *   Phase 1.5— update_state():  updates V/Sbus pointers. Call every compute_pf.
 *   Phase 2  — build_J_sparsity(): symbolic J build + value_map. Call when topology changes.
 *   Phase 3  — fill_J(): fast numerical fill via value_map_. Call each factorisation.
 *
 * J has the structure:
 * 
 *   | J11 | J12 |               | (pvpq, pvpq) | (pvpq, pq) |
 *   | --------- | = dimensions: | ------------------------- |
 *   | J21 | J22 |               |  (pq, pvpq)  | (pq, pq)   |
 * 
 * Python implementation:
 * 
 *   J11 = dS_dVa[array([pvpq]).T, pvpq].real
 *   J12 = dS_dVm[array([pvpq]).T, pq].real
 *   J21 = dS_dVa[array([pq]).T, pvpq].imag
 *   J22 = dS_dVm[array([pq]).T, pq].imag
 *
 * Where pvpq is the concatenation of pv then pq bus (might be sorted in the future)
 * And pv is the pv indexes.
 * 
 * The "state variables" for this base class are then:
 *    - all the theta (complex voltage angle) at each pv or pq nodes 
 *      (theta at slack buses are 0 by definition)
 *    - all the V (complex voltage magnitude) at each pq nodes
 * 
 * This means there are exactly 2 * n_pq + n_pv "state variables" for this base class.
 * 
 * This structure will not be changed by future additions.
 * 
 * And now, for each extension you have the following structure:
 * 
 *               |J_base      | J_coupling |
 * J_augmented = |------------|------------|
 *               |J_interface | J_features |
 * 
 * And they can "recursively" defined J_coupling, J_interface and J_feature
 * depending on the structure of the previous jacobians build.
 * 
 * They should NEVER rewrite J_base. They can add as many row / columns
 * as they want: each row / colum represents an added "state variable".
 * 
 * Rewriting J_base would mean modifying the state variables of the base
 * class which would break the "inheritance" / "composition" pattern that
 * we defined here.
 */
template <typename... Rest>
class NRSystem
{
protected:
    struct Contrib { int jrow, jcol, ybus_k; };

public:
    NRSystem() noexcept:
        Ybus_ptr_(nullptr),
        Sbus_ptr_(nullptr),
        nb_pv_(0),
        nb_pq_(0),
        nb_pvpq_(0),
        need_full_rebuild_(true),
        timer_dSbus_(0.),
        timer_fillJ_(0.) {}

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

    // ----- Phase 2: build J sparsity + value_map (non-virtual) ------------------

    void build_J_sparsity();  // calls build_value_maps

    // ----- Phase 3: fill J numerically (non-virtual) ----------------------------

    void fill_J();
    void fill_internal_variables();

    // ----- NR iteration primitives -----------------------------------------------

    virtual RealVect   mismatch()                           const;
    virtual void       apply_step(const RealVect& dx);
    virtual real_type  mismatch_sq_norm_at(const RealVect& dx) const;

    // ----- Housekeeping ----------------------------------------------------------

    void clear_jacobian() {
        J_         = Eigen::SparseMatrix<real_type, Eigen::ColMajor>();
        dS_dVm_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();
        dS_dVa_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();

        map_j11_.clear();
        map_j21_.clear();
        map_j12_.clear();
        map_j22_.clear();
        need_full_rebuild_ = true;
    }

    Eigen::Ref<const Eigen::SparseMatrix<real_type> > J()  const { return J_; }
    Eigen::Ref<const CplxVect> V()  const { return V_; }
    Eigen::Ref<const RealVect> Va() const { return Va_; }
    Eigen::Ref<const RealVect> Vm() const { return Vm_; }

    // to be implemented by all other classes
    virtual size_t total_state_variables() const {return nb_pvpq_ + nb_pq_;}

    // ----- Size / segment accessors ---------------------------
    int theta_size() const { return nb_pvpq_; }
    int vm_size()    const { return nb_pq_; }

    // Core state vector segments: theta first, then vm, then extension variables.
    RealVect::SegmentReturnType      theta(RealVect& x)       const { return x.segment(0, theta_size()); }
    RealVect::ConstSegmentReturnType theta(const RealVect& x) const { return x.segment(0, theta_size()); }
    RealVect::SegmentReturnType      vm   (RealVect& x)       const { return x.segment(theta_size(), vm_size()); }
    RealVect::ConstSegmentReturnType vm   (const RealVect& x) const { return x.segment(theta_size(), vm_size()); }

    // ----- Timers ----------------------------------------------------------------

    double timer_dSbus() const { return timer_dSbus_; }
    double timer_fillJ() const { return timer_fillJ_; }
    void   reset_timers()      { timer_dSbus_ = 0.; timer_fillJ_ = 0.; }

private:
    // ---- Shared data (one copy, inherited by all extensions) --------------------
    Eigen::VectorXi                        pv_, pq_, pvpq_;
    std::vector<int>                       pvpq_inv_, pq_inv_;
    int                                    nb_pv_, nb_pq_, nb_pvpq_;
    RealVect                               Va_, Vm_;
    CplxVect                               V_;
    Eigen::SparseMatrix<real_type, Eigen::ColMajor>         J_;
    Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>         dS_dVm_, dS_dVa_;
    // std::vector<cplx_type*>                value_map_;  // todo remove
    std::vector<int>                       map_j11_;
    std::vector<int>                       map_j21_;
    std::vector<int>                       map_j12_;
    std::vector<int>                       map_j22_;
    bool                                   need_full_rebuild_;
    double                                 timer_dSbus_, timer_fillJ_;

    std::tuple<Rest...> extensions_; // Holds the state for HVDC, DistSlack, etc.

protected:
    // visible attribute for derived class (non owning ptr)
    const LSGrid *                                         lsgrid_ptr_;
    const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>* Ybus_ptr_;
    const CplxVect*                                        Sbus_ptr_;

    // protected getters (const)
    Eigen::Ref<const Eigen::VectorXi> pv() const { return pv_; } 
    Eigen::Ref<const Eigen::VectorXi> pq() const { return pq_; }
    Eigen::Ref<const Eigen::VectorXi> pvpq() const { return pvpq_; }
    const std::vector<int> &          pvpq_inv() const {return pvpq_inv_; }
    const std::vector<int> &          pq_inv() const {return pq_inv_; }
    int                               nb_pv() const {return nb_pv_;}
    int                               nb_pq() const {return nb_pq_;}
    int                               nb_pvpq() const {return nb_pvpq_;}
    
    Eigen::Ref<const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor> > dS_dVm()  const { return dS_dVm_; }
    Eigen::Ref<const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor> > dS_dVa()  const { return dS_dVa_; }
    const std::vector<int> &          map_j11() const {return map_j11_; }
    const std::vector<int> &          map_j21() const {return map_j21_; }
    const std::vector<int> &          map_j12() const {return map_j12_; }
    const std::vector<int> &          map_j22() const {return map_j22_; }

    // TODO NR refacto
    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
    CplxVect _compute_trial_V(const RealVect& dx) const;
    RealVect _mismatch_core(const CplxVect& V_trial) const;


    // protected, but might be private, I don't really know.
    void _build_value_map(
        const std::vector< std::vector<Contrib> > & cijs
    );

    // definitely protected for this one
    int find_J_pos(int row, int col) const {
        int start = J_.outerIndexPtr()[col];
        int end   = J_.outerIndexPtr()[col + 1];
        const int* inner = J_.innerIndexPtr();
        auto it = std::lower_bound(inner + start, inner + end, row);
        if (it == inner + end || *it != row) return -1;
        return (int)(it - inner);
    };

private:
    // private members to combine  the extension features
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
            this,
            lsgrid_ptr,
            Ybus,
            Sbus,
            slack_weights
            ), 0)... };
        (void)dummy;
    }

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
// template <typename... Rest>
// class NRSystem<MultiSlack, Rest...>
// {
// public:
//     NRSystem() noexcept
//         : Base(), slack_bus_id_(0), slack_absorbed_(static_cast<real_type>(0.)) {}

//     virtual ~NRSystem() = default;

//     // ----- Phase 1 ---------------------------------------------------------------
//     void init_topology(
//         const Eigen::SparseMatrix<cplx_type>& Ybus,
//         const CplxVect&                        Sbus,
//         Eigen::Ref<const IntVect>              slack_ids,
//         const RealVect&                        slack_weights,
//         Eigen::Ref<const IntVect>              pv,
//         Eigen::Ref<const IntVect>              pq) override;

//     // ----- Phase 1.5 -------------------------------------------------------------
//     virtual void update_state(
//         const Eigen::SparseMatrix<cplx_type>& Ybus,
//         const CplxVect&                        V_init,
//         const CplxVect&                        Sbus) override;

//     // ----- NR primitives ---------------------------------------------------------
//     virtual RealVect  mismatch()                              const override;
//     virtual void      apply_step(const RealVect& dx)                override;
//     virtual real_type mismatch_sq_norm_at(const RealVect& dx) const override;

// protected:
//     // Appends slack row + slack column triplets after the core block
//     // virtual void _collect_J_triplets(std::vector<Eigen::Triplet<double>>& coeffs) const override;

//     // Returns nullptr for the slack row (constant); delegates to Base otherwise
//     // virtual cplx_type* _get_entry_ptr(int row, int col) override;

// private:
//     size_t    slack_bus_id_;
//     RealVect  slack_weights_;
//     real_type slack_absorbed_;

//     int _J_slack_row() const { return this->theta_size() + this->vm_size(); }

//     void   _append_slack_triplets(std::vector<Eigen::Triplet<double>>& coeffs) const;
//     RealVect _mismatch_with_slack(const CplxVect& V_trial, real_type sa) const;
// };

// ---- Type aliases (keep existing names working) --------------------------------

using SingleSlackNRSystem = NRSystem<>;
using MultiSlackNRSystem  = NRSystem<MultiSlack>;

} // namespace ls2g

#include "NRSystem.tpp"

#endif // NR_SYSTEM_H
