// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef REBUILD_OUTER_INNER_H
#define REBUILD_OUTER_INNER_H

// Solvers.hpp (not NRAlgo.hpp): out of the core build it declares the NRAlgo instantiations
// `extern template`, so the Newton wrapped below links against the core's.
#include <Solvers.hpp>
#include <powerflow_algorithm/NROuterInner.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <vector>

namespace ls2g {

/**
 * The outer loops' Newton the way OpenLoadFlow (v2.3.0) runs it: the Jacobian holds only the
 * equations active now, and a change of that set rebuilds it -- a new sparsity, a new symbolic
 * analysis -- where NROuterInner reserves every state up front and switches by value.
 *
 * The same driver (OuterLoopAlgo), the same loops, the same controls and the same Newton
 * formulation as NROuter_*: only what is handed to the Newton differs. The loops still reserve
 * every control they may use (they hold pointers to them); this class hands the Newton the
 * state the controls are in:
 *  - a reserved PV bus a loop (or a batch row) made PQ is a pq bus, the others pv buses: no
 *    switchable bus, no pinned row;
 *  - the slots holding a voltage controller by value only while one is held;
 *  - a solved column for a transformer's shift or ratio, or a shunt's susceptance, only while
 *    its control is on.
 * The masked buses (a batch row's islands) stay identity rows.
 *
 * What it costs is OpenLoadFlow's (equations/JacobianMatrix.java):
 *  - another set of active equations (that state, the held controllers and the masked buses
 *    included): the pattern is rebuilt and analyzed again (STRUCTURE_INVALID);
 *  - other Ybus values (a batch row opening a branch) or another hvdc regime: a numeric
 *    factorization on the analysis kept (VALUES_AND_ZEROS_INVALID);
 *  - anything else (a tap or a section moved, a target): a refactorization, redone as a
 *    factorization when it fails or its pivot growth is too small (VALUES_INVALID).
 *
 * A rebuild in the middle of a solve sets the Newton up again (NRAlgo::setup), which starts its
 * state again from the grid: the voltages, the controllers' reactive outputs and the ratio taps
 * and shunt sections a loop moved are carried over; Ybus is put back as the grid stamped it
 * first, since the transformers and shunts patch it from there. The continuous shift / ratio /
 * susceptance of a control that stays on restarts from its tap.
 */
template<class LinearSolver>
class RebuildOuterInner
{
public:
    using NRSys = SingleSlackNRSystem;
    using Algo = NRAlgo<LinearSolver, NRSys>;

    RebuildOuterInner() : controls_(algo_.system()) {
        // OpenLoadFlow's incremental update: a refactorization that fails, or whose pivot
        // growth is too small, is redone as a factorization
        algo_.set_refactor_fallback(true);
    }

    Algo & algo() { return algo_; }
    const Algo & algo() const { return algo_; }
    OuterControls & controls() { return controls_; }
    const OuterControls & controls() const { return controls_; }

    bool begin() {
        in_solve_ = nullptr;  // the caller's, for this solve only (see set_in_solve)
        return algo_.begin_solve();
    }

    // The driver's setup, on its private Ybus. The structure handed to the Newton is the one
    // its J was last built for: the loops have not initialized yet, and the first
    // apply_controls rebuilds it if they want another.
    template<class Declare>
    bool setup(Eigen::SparseMatrix<cplx_type>   & Ybus,
               const Eigen::Ref<const CplxVect> & V,
               const Eigen::Ref<const CplxVect> & Sbus,
               const Eigen::Ref<const IntVect>  & slack_ids,
               const Eigen::Ref<const RealVect> & slack_weights,
               const Eigen::Ref<const IntVect>  & pv,
               const Eigen::Ref<const IntVect>  & pq,
               bool                             & need_init,
               Declare                          && declare)
    {
        algo_.system().set_branch_mutable_ybus(&Ybus);
        nb_bus_ = V.size();
        ref_bus_ = slack_ids.size() > 0 ? slack_ids(0) : -1;
        // what a rebuild in the middle of the solve sets up again: Sbus is the driver's
        // private copy, which does not move during a solve
        ybus_ = &Ybus;
        sbus_ = Sbus.data();
        slack_ids_ = slack_ids;
        slack_weights_ = slack_weights;
        pv_ = pv;
        pq_ = pq;
        // the grid's stamp, before the transformers and shunts patch it; other values than
        // the last solve's are terms switched on or off (a row's contingency)
        std::vector<cplx_type> values(Ybus.valuePtr(), Ybus.valuePtr() + Ybus.nonZeros());
        const bool terms_changed = has_built_ && values != ybus_values_;
        ybus_values_.swap(values);

        Structure s = built_;
        s.masked = _masked();
        if (!has_built_) algo_.request_rebuild();
        IntVect pv_s, pq_s;
        _split(s, pv_s, pq_s);
        // what NROuterInner rebuilds on (the first solve, the driver's request, a change the
        // caller told): the loops reserve their controls again, as they do there
        bool ok = algo_.setup(Ybus, V, Sbus, slack_ids, slack_weights, pv_s, pq_s, need_init,
                              [&](){ controls_.clear_reservations(); declare(); _reserve(s); });
        // other masked buses alone (a batch row's islands): J is rebuilt for them on the
        // reservations kept, as _restructure does. The loops declare for the row's view, so
        // declaring again here would take out what this row opens for every row after it.
        if (ok && !need_init && has_built_ && !(s == built_)) {
            algo_.request_rebuild();
            ok = algo_.setup(Ybus, V, Sbus, slack_ids, slack_weights, pv_s, pq_s, need_init,
                             [&](){ _reserve(s); });
        }
        if (!ok) return false;
        built_ = s;
        has_built_ = true;
        analyze_pending_ = analyze_pending_ || need_init;
        if (terms_changed && !analyze_pending_) algo_.request_factorize();
        return true;
    }

    // The transformers the next setup's Ybus holds. Kept for the rebuilds of this solve, which
    // set the Newton up again after the driver has handed it back (nullptr).
    void set_in_solve(std::function<bool(int, int, int)> in_solve) {
        if (in_solve) in_solve_ = in_solve;
        BranchControl * branch = algo_.system().branch_control();
        if (branch != nullptr) branch->set_in_solve(std::move(in_solve));
    }

    bool newton(int max_iter, real_type tol, bool need_init) {
        const bool res = algo_.newton(max_iter, tol, need_init || analyze_pending_);
        // a solve that converged in zero iterations factorized nothing
        if (algo_.get_nb_iter() > 0) analyze_pending_ = false;
        return res;
    }

    // a single-slack Newton shares nothing
    real_type take_slack_absorbed() { return 0.; }
    real_type published_slack_absorbed() const { return 0.; }
    int reference_bus() const { return ref_bus_; }
    void finalize() { algo_.finalize(); }
    void request_rebuild() { algo_.request_rebuild(); }

    // Before each solve: rebuild the Newton if the active equations changed, then push what
    // the loops asked of their controls, as NROuterInner does (no pinning: a PQ bus is one).
    void apply_controls();

    void fill_context(OuterContext & ctx, RealVect & controller_q) const {
        controller_q = algo_.system().controller_q();
        ctx.controller_q = &controller_q;
        ctx.slack_absorbed = 0.;
    }

    // a PV bus has no magnitude unknown
    std::vector<bool> vm_unknown() const {
        const std::vector<int> & vm_col = algo_.system().vm_to_J_col();
        std::vector<bool> res(vm_col.size(), false);
        for (std::size_t bus = 0; bus < vm_col.size(); ++bus) res[bus] = vm_col[bus] >= 0;
        return res;
    }

    void control_layout(OuterControlLayout & out) const {
        out = OuterControlLayout();
        const BranchControl * branch = algo_.system().branch_control();
        if (branch != nullptr) branch->describe(out);
        const ShuntControl * shunt = algo_.system().shunt_control();
        if (shunt != nullptr) shunt->describe(out);
    }
    void phase_tap(std::vector<int> & positions) const {
        positions.clear();
        const BranchControl * branch = algo_.system().branch_control();
        if (branch != nullptr) branch->positions(positions);
    }
    void ratio_tap(std::vector<int> & positions) const {
        positions.clear();
        const BranchControl * branch = algo_.system().branch_control();
        if (branch != nullptr) branch->ratio_positions(positions);
    }
    void shunt_sections(std::vector<int> & counts) const {
        counts.clear();
        const ShuntControl * shunt = algo_.system().shunt_control();
        if (shunt != nullptr) shunt->section_counts(counts);
    }

    // how many times the active equations changed (a rebuild each) since the algorithm was made
    int nb_restructure() const { return nb_restructure_; }

private:
    // the set of active equations, OpenLoadFlow's equation index
    struct Structure {
        std::vector<int> pq_switched;  // reserved PV buses that are PQ (solver ids)
        std::vector<int> held;         // voltage controllers held (plan positions)
        std::vector<int> shift_on;     // transformers whose shift is solved
        std::vector<int> ratio_on;     // transformers whose ratio control is on
        std::vector<int> shunt_on;     // shunt controller buses whose control is on
        std::vector<int> masked;       // buses out of the solve (identity rows)
        bool operator==(const Structure & o) const {
            return pq_switched == o.pq_switched && held == o.held && shift_on == o.shift_on &&
                   ratio_on == o.ratio_on && shunt_on == o.shunt_on && masked == o.masked;
        }
    };

    // the state the controls are in now
    Structure _current() const;
    std::vector<int> _masked() const {
        std::vector<int> res = algo_.system().masked_buses();
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        return res;
    }
    // the labelling's pv / pq, the switched buses moved to pq
    void _split(const Structure & s, IntVect & pv, IntVect & pq) const;
    // the structure, in the system, right before it claims its rows and columns
    void _reserve(const Structure & s);
    // the Newton set up again for `s`, its state carried over
    bool _restructure(const Structure & s);

    static bool _contains(const std::vector<int> & sorted, int v) {
        return std::binary_search(sorted.begin(), sorted.end(), v);
    }

    Algo algo_;
    NROuterControls<NRSys> controls_;
    Eigen::Index nb_bus_ = 0;
    int ref_bus_ = -1;

    // the last setup's input, for a rebuild in the middle of the solve
    Eigen::SparseMatrix<cplx_type> * ybus_ = nullptr;
    std::vector<cplx_type> ybus_values_;  // its values as the grid stamped them
    const cplx_type * sbus_ = nullptr;
    IntVect slack_ids_;
    RealVect slack_weights_;
    IntVect pv_;
    IntVect pq_;
    std::function<bool(int, int, int)> in_solve_;

    Structure built_;            // what J was last built for
    bool has_built_ = false;
    bool analyze_pending_ = false;
    std::vector<int> hvdc_applied_;  // the regimes the last solve ran with
    int nb_restructure_ = 0;
};

template<class LinearSolver>
typename RebuildOuterInner<LinearSolver>::Structure RebuildOuterInner<LinearSolver>::_current() const
{
    Structure s;
    const std::vector<int> all = controls_.switchable_buses();  // sorted
    const std::vector<int> pv = controls_.pinned_buses();       // sorted
    std::set_difference(all.begin(), all.end(), pv.begin(), pv.end(), std::back_inserter(s.pq_switched));
    const std::vector<real_type> held = controls_.held_q();
    for (std::size_t c = 0; c < held.size(); ++c) {
        if (std::isfinite(held[c])) s.held.push_back(static_cast<int>(c));
    }
    for (const auto & ts : controls_.phase_shifters()) {
        if (ts.second->solves_shift() && ts.second->requested_control() == 1) s.shift_on.push_back(ts.first);
    }
    for (const auto & tr : controls_.ratio_taps()) {
        if (tr.second->requested_control() == 1) s.ratio_on.push_back(tr.first);
    }
    for (const auto & bc : controls_.shunt_controllers()) {
        if (bc.second->requested_control() == 1) s.shunt_on.push_back(bc.first);
    }
    s.masked = _masked();
    return s;
}

template<class LinearSolver>
void RebuildOuterInner<LinearSolver>::_split(const Structure & s, IntVect & pv, IntVect & pq) const
{
    std::vector<int> pv_s, pq_s(pq_.data(), pq_.data() + pq_.size());
    for (Eigen::Index i = 0; i < pv_.size(); ++i) {
        if (_contains(s.pq_switched, pv_(i))) pq_s.push_back(pv_(i));
        else pv_s.push_back(pv_(i));
    }
    std::sort(pq_s.begin(), pq_s.end());
    pv = Eigen::Map<const IntVect>(pv_s.data(), static_cast<Eigen::Index>(pv_s.size()));
    pq = Eigen::Map<const IntVect>(pq_s.data(), static_cast<Eigen::Index>(pq_s.size()));
}

template<class LinearSolver>
void RebuildOuterInner<LinearSolver>::_reserve(const Structure & s)
{
    NRSys & system = algo_.system();
    // a PQ bus is a pq bus: nothing switches by value. The reference bus is in neither pv nor
    // pq: switched PQ (its generators frozen), it keeps its angle and gets a magnitude unknown
    // and its Q equation, the free magnitude of a switchable bus left unpinned.
    std::vector<int> free_vm;
    for (Eigen::Index i = 0; i < slack_ids_.size(); ++i) {
        if (_contains(s.pq_switched, slack_ids_(i))) free_vm.push_back(slack_ids_(i));
    }
    algo_.set_switchable_vm_buses(free_vm);
    system.set_may_hold_voltage_controllers(!s.held.empty());
    {
        std::vector<int> trafos;
        std::vector<char> with_column;
        for (const PhaseShifterControl * shifter : controls_.phase_shifters_reserved()) {
            trafos.push_back(shifter->trafo());
            with_column.push_back(_contains(s.shift_on, shifter->trafo()) ? 1 : 0);
        }
        system.set_phase_controllers(trafos, with_column);
    }
    // a group keeps its transformers (their taps patch Ybus), solved while one of them is on
    std::vector<BranchControl::RatioGroupDecl> groups;
    for (const auto & g : controls_.ratio_groups()) {
        BranchControl::RatioGroupDecl d;
        d.bus_solver = g.bus;
        d.target_vm = g.target_vm;
        d.trafos = g.trafos;
        d.solved = false;
        for (int t : g.trafos) d.solved = d.solved || (g.solved && _contains(s.ratio_on, t));
        groups.push_back(d);
    }
    system.set_ratio_groups(groups);
    std::vector<ShuntControl::GroupDecl> shunt_groups;
    for (const auto & g : controls_.shunt_groups()) {
        ShuntControl::GroupDecl d;
        d.bus_solver = g.bus;
        d.target_vm = g.target_vm;
        d.controller_buses = g.controller_buses;
        d.shunts = g.shunts;
        d.solved = false;
        for (int bus : g.controller_buses) d.solved = d.solved || (g.solved && _contains(s.shunt_on, bus));
        shunt_groups.push_back(d);
    }
    system.set_shunt_groups(shunt_groups);
}

template<class LinearSolver>
bool RebuildOuterInner<LinearSolver>::_restructure(const Structure & s)
{
    NRSys & system = algo_.system();
    // what setting the Newton up again starts from the grid
    const CplxVect V = system.V();
    const RealVect q = system.controller_q();
    std::vector<int> ratio_pos, sections;
    BranchControl * branch = system.branch_control();
    ShuntControl * shunt = system.shunt_control();
    if (branch != nullptr) branch->ratio_positions(ratio_pos);
    if (shunt != nullptr) shunt->section_counts(sections);

    // Ybus as the grid stamped it: the transformers and shunts patch it from there
    std::copy(ybus_values_.begin(), ybus_values_.end(), ybus_->valuePtr());
    IntVect pv, pq;
    _split(s, pv, pq);
    algo_.request_rebuild();
    if (branch != nullptr) branch->set_in_solve(in_solve_);
    bool need_init = false;
    const Eigen::Map<const CplxVect> Sbus(sbus_, nb_bus_);
    const bool ok = algo_.setup(*ybus_, V, Sbus, slack_ids_, slack_weights_, pv, pq, need_init,
                                [&](){ _reserve(s); });
    if (branch != nullptr) branch->set_in_solve(nullptr);
    if (!ok) return false;

    system.set_controller_q(q);
    const int none = std::numeric_limits<int>::min();
    if (branch != nullptr) {
        for (std::size_t t = 0; t < ratio_pos.size(); ++t) {
            if (ratio_pos[t] != none) branch->set_ratio_tap(static_cast<int>(t), ratio_pos[t]);
        }
    }
    if (shunt != nullptr && !sections.empty()) {
        for (const auto & g : controls_.shunt_groups()) {
            for (std::size_t c = 0; c < g.controller_buses.size() && c < g.shunts.size(); ++c) {
                std::vector<int> counts;
                for (int sh : g.shunts[c]) {
                    if (sh >= 0 && static_cast<std::size_t>(sh) < sections.size()) counts.push_back(sections[static_cast<std::size_t>(sh)]);
                }
                if (counts.size() == g.shunts[c].size()) shunt->set_sections(g.controller_buses[c], counts);
            }
        }
    }
    built_ = s;
    analyze_pending_ = true;
    ++nb_restructure_;
    return true;
}

template<class LinearSolver>
void RebuildOuterInner<LinearSolver>::apply_controls()
{
    // OpenLoadFlow's STRUCTURE_INVALID: another set of active equations
    const Structure s = _current();
    if (!(s == built_)) _restructure(s);

    NRSys & system = algo_.system();
    // OpenLoadFlow's VALUES_AND_ZEROS_INVALID: the hvdc lines' terms switched
    const std::vector<int> regimes = controls_.hvdc_regimes();
    if (!regimes.empty()) system.set_hvdc_status_override(regimes, HvdcRegimeControl::KEEP);
    if (regimes != hvdc_applied_ && !hvdc_applied_.empty() && !analyze_pending_) algo_.request_factorize();
    hvdc_applied_ = regimes;
    // what follows is NROuterInner::apply_controls, without the pinning
    const std::vector<real_type> svc_target_vm = controls_.svc_target_vm();
    if (!svc_target_vm.empty()) system.release_held_svcs(svc_target_vm);
    system.set_held_voltage_controllers(controls_.held_q());
    BranchControl * phase = system.branch_control();
    if (phase != nullptr) {
        for (const auto & ts : controls_.phase_shifters()) {
            if (ts.second->tap_moved()) phase->set_tap(ts.first, ts.second->requested_tap());
        }
        for (const auto & ts : controls_.phase_shifters()) {
            if (ts.second->requested_control() >= 0) phase->set_control_on(ts.first, ts.second->requested_control() == 1);
        }
        for (const auto & tr : controls_.ratio_taps()) {
            int position;
            if (tr.second->take_tap(position)) phase->set_ratio_tap(tr.first, position);
        }
        for (const auto & tr : controls_.ratio_taps()) {
            if (tr.second->requested_control() >= 0) phase->set_ratio_control_on(tr.first, tr.second->requested_control() == 1);
        }
    }
    ShuntControl * shunt = system.shunt_control();
    if (shunt != nullptr) {
        std::vector<int> counts;
        for (const auto & bc : controls_.shunt_controllers()) {
            if (bc.second->take_sections(counts)) shunt->set_sections(bc.first, counts);
        }
        for (const auto & bc : controls_.shunt_controllers()) {
            if (bc.second->requested_control() >= 0) shunt->set_control_on(bc.first, bc.second->requested_control() == 1);
        }
    }
    std::vector<int> buses;
    std::vector<real_type> vm;
    controls_.take_pending_vm(buses, vm);
    if (!buses.empty()) system.set_vm_at(buses, vm);
}

}  // namespace ls2g

#endif  // REBUILD_OUTER_INNER_H
