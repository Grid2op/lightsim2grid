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

enum LS2G_API dSdX { dSdVa_r, dSdVa_i, dSdVm_r, dSdVm_i, NotIndSdV};

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

// ---- Base class definition ----------------------------------------------------
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

        static int find_J_pos_csr(
            Eigen::Ref<const Eigen::SparseMatrix<real_type, Eigen::RowMajor> > J_csr,
            int row,
            int col){
            int start = J_csr.outerIndexPtr()[row];
            int end   = J_csr.outerIndexPtr()[row + 1];
            const int* inner = J_csr.innerIndexPtr();
            auto it = std::lower_bound(inner + start, inner + end, col);
            if (it == inner + end || *it != col) return -1;
            return (int) (it - inner);
        };

        // call at the beginning of each solve
        void update_state(
            const LSGrid *                         lsgrid_ptr,
            const Eigen::SparseMatrix<cplx_type>&  Ybus,
            const CplxVect&                        Sbus,
            const RealVect&                        slack_weights
        ){
            lsgrid_ptr_ = lsgrid_ptr;
            Ybus_ptr_ = &Ybus;
            Sbus_ptr_ = &Sbus;
        }

        // call after update_state
        // at the beginning of each solve
        // only if the topology has changed
        void init_topology(
            Eigen::Ref<const IntVect>              slack_ids,
            const RealVect&                        slack_weights,
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
            const int n_bus  = static_cast<int>(Ybus_ptr_->rows());

            pvpq_inv_.assign(n_bus, -1);
            for (int i = 0; i < nb_pvpq_; ++i) pvpq_inv_[pvpq_(i)] = i;
            pq_inv_.assign(n_bus, -1);
            for (int i = 0; i < nb_pq_; ++i) pq_inv_[pq_(i)] = i;
        }

        size_t get_size() const { return nb_pvpq_ + nb_pq_;}

        std::vector<Contrib> build_J_contrib() noexcept
        {
            // compute its sparsity pattern
            const Eigen::SparseMatrix<cplx_type> & Ybus = *Ybus_ptr_;
            const int nnz_Y  = Ybus.nonZeros();

            // get the triplets (will be in a virtual function later)
            std::vector<Contrib> cntrb;
            cntrb.reserve(4 * nnz_Y);
            // TODO this might be overly pessimistic (only pvpq comp are kept for c11 and c21)
            // TODO this might be overly pessimistic (only pq comp are kept for c12 and c22) 

            int k = 0;
            for (int outer = 0; outer < Ybus.outerSize(); ++outer) {
                for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator
                    it(Ybus, outer); it; ++it, ++k)
                {
                    int i = (int)it.row(), j = (int)it.col();
                    int ri = pvpq_inv_[i], rq = pq_inv_[i];
                    int ci = pvpq_inv_[j], cq = pq_inv_[j];
                    if (ri >= 0 && ci >= 0) cntrb.push_back({ri,          ci,          k, dSdVa_r});
                    if (ri >= 0 && cq >= 0) cntrb.push_back({ri,          nb_pvpq_ + cq, k, dSdVm_r});
                    if (rq >= 0 && ci >= 0) cntrb.push_back({nb_pvpq_ + rq, ci,          k, dSdVa_i});
                    if (rq >= 0 && cq >= 0) cntrb.push_back({nb_pvpq_ + rq, nb_pvpq_ + cq, k, dSdVm_i});
                }
            }

            return cntrb;
        }

        void clear_jacobian() {}

        void fill_J() {}

        const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor> * Ybus_ptr() const {return Ybus_ptr_;}

    private:
        Eigen::VectorXi                        pv_, pq_, pvpq_;
        std::vector<int>                       pvpq_inv_, pq_inv_;
        int                                    nb_pv_, nb_pq_, nb_pvpq_;


        // visible attribute for derived class (non owning ptr)
        const LSGrid *                                         lsgrid_ptr_;
        const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>* Ybus_ptr_;
        const CplxVect*                                        Sbus_ptr_;


    public:
        Eigen::Ref<const Eigen::VectorXi> pv() const { return pv_; } 
        Eigen::Ref<const Eigen::VectorXi> pq() const { return pq_; }
        Eigen::Ref<const Eigen::VectorXi> pvpq() const { return pvpq_; }
        const std::vector<int> &          pvpq_inv() const {return pvpq_inv_; }
        const std::vector<int> &          pq_inv() const {return pq_inv_; }
        int                               nb_pv() const {return nb_pv_;}
        int                               nb_pq() const {return nb_pq_;}
        int                               nb_pvpq() const {return nb_pvpq_;}

        int theta_size() const { return nb_pvpq_; }
        int vm_size()    const { return nb_pq_; }
};

// ---- Extension tag types ------------------------------------------------------

/**
 * This extension adds the ability to have a distributed slack directly in the jacobian.
 * 
 * It adds exactly "nb slack" extra variables (nb_slack being the
 * number of slacks in the grid)
 * 
 * The first "nb_slack" - 1 are pv nodes: theta is the extra state variables.
 * 
 * The last state variable represents the p mismatch at each slack bus.
 */
class LS2G_API MultiSlack   // distributed-slack extension
{
    public:
        // call at the beginning of each NR solve
        // once. It is not called for
        // NR iterations
        void update_state(
            const Base *                           nr_system_base_ptr,
            const LSGrid *                         lsgrid_ptr,
            const Eigen::SparseMatrix<cplx_type>&  Ybus,
            const CplxVect&                        Sbus,
            const RealVect&                        slack_weights
        ){
            // TODO remember slack weights !
            nr_system_base_ptr_ = nr_system_base_ptr;
            slack_weights_ = slack_weights;
        }

        // call after update_state
        // at the beginning of each solve
        // only if the topology has changed
        void init_topology(
            Eigen::Ref<const IntVect>              slack_ids,
            const RealVect&                        slack_weights,
            Eigen::Ref<const IntVect>              pv,
            Eigen::Ref<const IntVect>              pq
        ) {
            my_size_ = slack_ids.size();
            // TODO: sort the slack_ids
            // TODO take the last slack (easier to fill the J matrix !)
            ref_slack_id_ = static_cast<size_t>(slack_ids[0]);
            slack_ids_ = slack_ids;
            pv_slack_inv_.assign(nr_system_base_ptr_->Ybus_ptr()->cols(), -1);
            if(my_size_ > 1){
                // real distributed slack
                const int nb_slack_pv = my_size_ - 1;
                pv_slack_ = slack_ids.segment(1, nb_slack_pv);  // TODO sort this
                for (int i = 0; i < nb_slack_pv; ++i) pv_slack_inv_[pv_slack_(i)] = i;
            }

        }

        // called after update_state
        // (and init_topology)
        // at the beginning of each solve
        size_t get_size() const {
            return my_size_;
        }

        // called when topology has changed
        // after update_state, init_topology (and get_size)

        /**
        J has the shape:
        
        | J11 | J12 | S13 |   |               | (pvpq_base, pvpq_base) | (pvpq_base, pq) | (pvpq_base, pvslack) | (pvpq_base, 1) |
        | --------------- |   |               | --------------------------------------------------------------- |                |
        | J21 | J22 | S23 | s |               |    (pq, pvpq_base)     |    (pq, pq)     |    (pq, pvslack)     |     (pq, 1)    |
        | --------------- | l | = dimensions: | --------------------------------------------------------------- |                |
        | S31 | S32 | S33 | a |               |  (pvslack, pvpq_base)  |  (pvslack, pq)  |  (pvslack, pvslack)  |  (pvslack, 1)  |
        |-----------------| c |               | --------------------------------------------------------------- |                |
        |    slack_bus    | k |                

        python implementation (base, not modified) :

        `J11` = dS_dVa[array([pvpq]).T, pvpq].real
        `J12` = dS_dVm[array([pvpq]).T, pq].real
        `J21` = dS_dVa[array([pq]).T, pvpq].imag
        `J22` = dS_dVm[array([pq]).T, pq].imag

        The J_coupling is :

        - `S13` = dS_dVa[array([pvpq]).T, pv_slack].real
        - `S23` = dS_dVa[array([pq]).T, pv_slack].imag

        The J_interface is :

        - `S31` = dS_dVa[array([pv_slack]).T, pvpq_base].real
        - `S32` = dS_dVm[array([pv_slack]).T, pvpq_base].real

        Some part of J feature are:
        
        - `S33` = dS_dVa[pvslack, pvslack].real

        `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real
        and dS_dVm[slack_bus_id, pq].real

        `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
        the remaining pq components are all 0.

        ┌───────┬────────────────┬─────────────────┬───────────────────────────────────┐
        │ Block │ Row equations  │  Col unknowns   │              Formula              │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ J11   │ ΔP / pvpq_base │ ΔVa / pvpq_base │ dS_dVa[pvpq_base, pvpq_base].real │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ J12   │ ΔP / pvpq_base │ ΔVm / pq        │ dS_dVm[pvpq_base, pq].real        │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ S13   │ ΔP / pvpq_base │ ΔVa / pvslack   │ dS_dVa[pvpq_base, pvslack].real   │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ J21   │ ΔQ / pq        │ ΔVa / pvpq_base │ dS_dVa[pq, pvpq_base].imag        │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ J22   │ ΔQ / pq        │ ΔVm / pq        │ dS_dVm[pq, pq].imag               │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ S23   │ ΔQ / pq        │ ΔVa / pvslack   │ dS_dVa[pq, pvslack].imag          │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ S31   │ ΔP / pvslack   │ ΔVa / pvpq_base │ dS_dVa[pvslack, pvpq_base].real   │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ S32   │ ΔP / pvslack   │ ΔVm / pq        │ dS_dVm[pvslack, pq].real          │
        ├───────┼────────────────┼─────────────────┼───────────────────────────────────┤
        │ S33   │ ΔP / pvslack   │ ΔVa / pvslack   │ dS_dVa[pvslack, pvslack].real     │
        └───────┴────────────────┴─────────────────┴───────────────────────────────────┘

        **/
        std::vector<Contrib> build_J_contrib(
            int offset
        ) noexcept {
            const auto & Ybus = * nr_system_base_ptr_->Ybus_ptr();
            const size_t nnz_Y = Ybus.nonZeros();
            my_offset_ = static_cast<size_t>(offset);
            int my_size = static_cast<int>(my_size_);

            std::vector<Contrib> cntrb; 
            cntrb.reserve(5 * nnz_Y);  // 5 * nnz_Y is overly pessimistic (too much memory allocated)

            const std::vector<int> & pvpq_inv = nr_system_base_ptr_ ->pvpq_inv();
            const std::vector<int> & pq_inv = nr_system_base_ptr_ ->pq_inv();
            int nb_pvpq_base = nr_system_base_ptr_->nb_pvpq();

            int k = 0;
            for (int outer = 0; outer < Ybus.outerSize(); ++outer) {
                for (Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>::InnerIterator
                    it(Ybus, outer); it; ++it, ++k)
                {
                    int i = (int)it.row(), j = (int)it.col();
                    int ri = pvpq_inv[i], rq = pq_inv[i];
                    int ci = pvpq_inv[j], cq = pq_inv[j];
                    int rs = pv_slack_inv_[i];
                    int cs = pv_slack_inv_[j];

                    // adapt the base triplets (J_interface)
                    // S13
                    if (ri >= 0 && cs >= 0) cntrb.push_back({ri, offset + cs, k, dSdVa_r});
                    // S23
                    if (rq >= 0 && cs >= 0) cntrb.push_back({nb_pvpq_base + rq, offset + cs, k, dSdVa_i});
                    
                    // add the new columns (J_coupling)
                    // S31
                    if (rs >= 0 && ci >= 0) cntrb.push_back({offset + rs, ci, k, dSdVa_r});

                    // S32
                    if (rs >= 0 && cq >= 0) cntrb.push_back({offset + rs, nb_pvpq_base + cq, k, dSdVm_r});

                    // add the new  J feature
                    // S33
                    if (rs >= 0 && cs >= 0) cntrb.push_back({offset + rs,     offset + cs, k, dSdVa_r});

                    // add ref slack equations
                    if (i == ref_slack_id_){
                        // ref slack P wrt theta pvpq
                        if(ci >= 0) {
                            cntrb.push_back({offset + my_size - 1, ci, k, dSdVa_r});
                        }
                        // ref slack P wrt Vm pq
                        if(cq >= 0){
                            cntrb.push_back({offset + my_size - 1, nb_pvpq_base + cq, k, dSdVm_r});
                        }
                        // ref slack P wrt theta pvslack
                        if(cs >= 0)
                        {
                            cntrb.push_back({offset + my_size - 1, offset + cs, k, dSdVa_r});
                        }
                    }
                }

            }
            // and now add the slack equations (last column), it is 
            // 0 for pvpq, pq and the slack weights for pvslack
            for(size_t i = 0; i < my_size_; ++i)
            {
                cntrb.push_back({offset + static_cast<int>(i), offset + my_size - 1, -1, NotIndSdV});
            }

            // now store them in a defined order
            return cntrb;
        }

        // this function is called to update the "value map"
        // It is part of possible "computation speed optimization"
        // so we advise to ignore it at first.
        // It is called each time a new J is defined.
        // It is called after update_state, init_topology and build_J_contrib.
        // The idea behind this is that it is faster to update the 
        // jacobian at each iteration of the Newton Raphson,
        // you can use a "value_map" that maps the new values you want 
        // to compute to their index in Jacobian.valuePtr()
        // So intead of computing Jacobian.coeffRef(row, col) = XXX
        // you will be able to call Jacobian.valuePtr()[value_map[i]] = XXX[i]
        void build_value_map(
            // Jacobian is the inialized (correct sparsity pattern) of the Jacobian
            Eigen::Ref<const Eigen::SparseMatrix<real_type, Eigen::ColMajor> > Jacobian,
            // my_contrib is the result of the call to build_J_contrib
            const std::vector<Contrib> & my_contrib
        )
        {
            // when filling the values I want to fill by 
            // looking to slack_ids directly in order
            // as the ref slack is the first one (TODO we might change that...)
            // there is a "lag" of 1
            map_dist_slack_ = std::vector<int>(my_size_); 
            for(size_t i = 0; i < my_size_ - 1; ++i)
            {
                // other slack (non ref)
                map_dist_slack_[i + 1] = Base::find_J_pos(Jacobian, my_offset_ + i, my_offset_ + (my_size_ - 1));
                // the [i+1] = XXX(my_offset_ + i, ...) (and not i+1) is volontary
                // because the first slack bus is taken as the ref slack.
            }
            // ref slack (first one)
            map_dist_slack_[0] = Base::find_J_pos(Jacobian, my_offset_ + (my_size_ - 1), my_offset_ + (my_size_ - 1));
        }


        // can be explicitly called any time.
        // if called, then the update_state, init_topology, add_triplets etc.
        // will be called again if another powerflow needs
        // to be run.
        void clear_jacobian(){
            my_size_ = 0;
            ref_slack_id_ = 0;
            my_offset_ = 0;
            slack_ids_ = IntVect::Zero(0);
            pv_slack_ = IntVect::Zero(0);
            pv_slack_inv_.clear();
            nr_system_base_ptr_ = nullptr;
            slack_weights_ = RealVect::Zero(0);
            map_dist_slack_.clear();
        }

        // called at each iteration of the Newton Raphson
        // when this is called, it assumes that update_state, init_topology 
        // etc. have been properly called.
        void fill_J(Eigen::Ref<Eigen::SparseMatrix<real_type, Eigen::ColMajor> > Jacobian) const
        {
            size_t i = 0;
            for(auto & c : map_dist_slack_){
                if(c == -1){
                    i++;
                    continue;
                }
                Jacobian.valuePtr()[c] = slack_weights_(slack_ids_(i));
                i++;
            }

        }

    private:
        size_t             my_size_;
        size_t             ref_slack_id_;
        size_t             my_offset_;
        Eigen::VectorXi    slack_ids_;
        Eigen::VectorXi    pv_slack_;
        std::vector<int>   pv_slack_inv_;
        const Base *       nr_system_base_ptr_;
        RealVect           slack_weights_;  // size: nb_bus
        std::vector<int>   map_dist_slack_;

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
class LS2G_API NRSystem<Base, Rest...>
{
public:
    NRSystem() noexcept:
        Ybus_ptr_(nullptr),
        Sbus_ptr_(nullptr),
        need_full_rebuild_(true),
        timer_dSbus_(0.),
        timer_fillJ_(0.),
        nb_total_state_variables_(0) {}

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
        base_.clear_jacobian();
        _clear_jacobian_extensions(std::make_index_sequence<sizeof...(Rest)>{});

        dS_dVm_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();
        dS_dVa_    = Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>();

        map_dsdva_r_.clear();
        map_dsdva_i_.clear();
        map_dsdvm_r_.clear();
        map_dsdvm_i_.clear();
        need_full_rebuild_ = true;
    }

    Eigen::Ref<const Eigen::SparseMatrix<real_type> > J()  const { return J_; }
    Eigen::Ref<const CplxVect> V()  const { return V_; }
    Eigen::Ref<const RealVect> Va() const { return Va_; }
    Eigen::Ref<const RealVect> Vm() const { return Vm_; }

    // to be implemented by all other classes
    size_t total_state_variables() const {return nb_total_state_variables_;}

    // ----- Size / segment accessors ---------------------------
    int theta_size() const { return base_.theta_size(); }
    int vm_size()    const { return base_.vm_size(); }

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
    RealVect                               Va_, Vm_;
    CplxVect                               V_;
    Eigen::SparseMatrix<real_type, Eigen::ColMajor>         J_;
    bool                                   need_full_rebuild_;
    double                                 timer_dSbus_, timer_fillJ_;
    size_t                                 nb_total_state_variables_;
    Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>         dS_dVm_, dS_dVa_;
    std::vector<int>                       map_dsdva_r_;
    std::vector<int>                       map_dsdva_i_;
    std::vector<int>                       map_dsdvm_r_;
    std::vector<int>                       map_dsdvm_i_;

    // Holds the base things
    Base                                   base_;

    // Holds the state for HVDC, DistSlack, etc.
    std::tuple<Rest...>                    extensions_; 

protected:
    // visible attribute for derived class (non owning ptr)
    const LSGrid *                                         lsgrid_ptr_;
    const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor>* Ybus_ptr_;
    const CplxVect*                                        Sbus_ptr_;

    // protected getters (const)
    Eigen::Ref<const Eigen::VectorXi> pv() const { return base_.pv(); } // TODO and also the slack "pv" for dist slack
    Eigen::Ref<const Eigen::VectorXi> pq() const { return base_.pq(); }
    Eigen::Ref<const Eigen::VectorXi> pvpq() const { return base_.pvpq(); }  // TODO and also the slack "pv" for dist slack
    const std::vector<int> &          pvpq_inv() const {return base_.pvpq_inv(); }  // TODO and also the slack "pv" for dist slack
    const std::vector<int> &          pq_inv() const {return base_.pq_inv(); }
    int                               nb_pv() const {return base_.nb_pv();}
    int                               nb_pq() const {return base_.nb_pq();}
    int                               nb_pvpq() const {return base_.nb_pvpq();}
    
    Eigen::Ref<const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor> > dS_dVm()  const { return dS_dVm_; }
    Eigen::Ref<const Eigen::SparseMatrix<cplx_type, Eigen::ColMajor> > dS_dVa()  const { return dS_dVa_; }

    // TODO NR refacto
    static CplxVect _reconstruct_V(const RealVect& Va, const RealVect& Vm);
    CplxVect _compute_trial_V(const RealVect& dx) const;
    RealVect _mismatch_core(const CplxVect& V_trial) const;


    void _build_value_map(
        const std::vector< std::vector<Contrib> > & cijs
    )
    {
        const auto & Ybus = * Ybus_ptr_;
        const size_t nnz_Y = Ybus.nonZeros();

        // now synch the map_j{1,2}{1,2} with the new J_ matrix
        map_dsdva_r_.assign(nnz_Y, -1);
        map_dsdva_i_.assign(nnz_Y, -1);
        map_dsdvm_r_.assign(nnz_Y, -1);
        map_dsdvm_i_.assign(nnz_Y, -1);

        // same order as build_J_contrib
        // int c11 = 0, c12 = 1, c21 = 2, c22 = 3;  
        // this would also avoid copy
        for(const auto& cij : cijs)
        {
            for (auto& c : cij)
            {
                switch (c.whichmat())
                {
                case dSdVa_r:
                    map_dsdva_r_[c.ybus_k()] = Base::find_J_pos(J_, c.jrow(), c.jcol());
                    break;
                case dSdVa_i:
                    map_dsdva_i_[c.ybus_k()] = Base::find_J_pos(J_, c.jrow(), c.jcol());
                    break;
                case dSdVm_r:
                    map_dsdvm_r_[c.ybus_k()] = Base::find_J_pos(J_, c.jrow(), c.jcol());
                    break;
                case dSdVm_i:
                    map_dsdvm_i_[c.ybus_k()] = Base::find_J_pos(J_, c.jrow(), c.jcol());
                    break;

                default:
                    break;
                }
                
            }
        }

        // build value map of the extensions
        _build_value_map_extensions(cijs, std::make_index_sequence<sizeof...(Rest)>{});

        // indicate it is fully initialized
        need_full_rebuild_ = false;
    }
    
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
            &base_,
            lsgrid_ptr,
            Ybus,
            Sbus,
            slack_weights
            ), 0)... };
        (void)dummy;
    }

    template <std::size_t... Is>
    void _update_total_state_variables(
        std::index_sequence<Is...>){
        nb_total_state_variables_ = base_.get_size();
        int dummy[] = { 0, (nb_total_state_variables_ += std::get<Is>(extensions_).get_size(), 0)... };
        (void) dummy;
    }

    template <std::size_t... Is>
    void _build_J_contrib_extensions(
        std::vector<std::vector<Contrib> >& contribs,
        std::index_sequence<Is...>) {
        int current_offset = base_.get_size();
        
        // We can't use the simple dummy array if we need to update 'current_offset' 
        // inside the expansion. We use a small recursive call or a braced-init loop.
        
        int dummy[] = { 0, (
            contribs.push_back(std::get<Is>(extensions_).build_J_contrib(current_offset)), 
            current_offset += std::get<Is>(extensions_).get_size(), 
            0
        )... };
        (void) dummy;
    }

    template <std::size_t... Is>
    void _clear_jacobian_extensions(
        std::index_sequence<Is...>) {                
        int dummy[] = { 0, (
            std::get<Is>(extensions_).clear_jacobian(),
            0
        )... };
    }

    template <std::size_t... Is>
    void _build_value_map_extensions(
        const std::vector< std::vector<Contrib> > & cijs,
        std::index_sequence<Is...>
    ) {
        size_t current_id = 1;  // nothing for base, so assigned to 1
        int dummy[] = { 0, (
            std::get<Is>(extensions_).build_value_map(J_, cijs[current_id]), 
            current_id += 1, 
            0
        )... };
    }

    template <std::size_t... Is>
    void _fill_J_extensions(
        std::index_sequence<Is...>) {                
        int dummy[] = { 0, (
            std::get<Is>(extensions_).fill_J(J_),
            0
        )... };
        (void) dummy;
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

using SingleSlackNRSystem = NRSystem<Base>;
using MultiSlackNRSystem  = NRSystem<Base, MultiSlack>;

} // namespace ls2g

#include "NRSystem.tpp"

#endif // NR_SYSTEM_H
