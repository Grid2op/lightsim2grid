// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_OUTER_INNER_H
#define NR_OUTER_INNER_H

#include <vector>

#include "NRAlgo.hpp"
#include "outer_loop/BaseOuterLoop.hpp"
#include "outer_loop/OuterControls.hpp"

namespace ls2g {

/**
 * The Newton-Raphson as the outer-loop driver (OuterLoopAlgo) sees it: what the driver
 * needs from the algorithm it wraps beyond BaseAlgo's interface. Everything here goes
 * through the Newton's system -- the controls the loops reserve in it (phase shifters,
 * ratio and shunt groups, held voltage controllers), the values they edit between two
 * solves, and what the loops read back -- so that the driver itself never sees an NRSystem.
 *
 * NRSystem must carry the extensions the loops act through (VoltageControl, BranchControl,
 * ShuntControl, Hvdc); a missing BranchControl or ShuntControl only disables what needs it.
 */
template<class LinearSolver, class NRSystem = SingleSlackNRSystem>
class NROuterInner
{
public:
    using Algo = NRAlgo<LinearSolver, NRSystem>;

    NROuterInner() {
        // always on: a value edit between two solves may shrink a pivot the fixed pivot
        // sequence then finds at zero
        algo_.set_refactor_fallback(true);
    }

    // the wrapped algorithm, for everything BaseAlgo already says (results, J, timers, the
    // batch setters) and the Newton's own parameters
    Algo & algo() { return algo_; }
    const Algo & algo() const { return algo_; }

    // what the loops reserve and act through
    OuterControls & controls() { return controls_; }
    const OuterControls & controls() const { return controls_; }

    // a new solve: false if the algorithm cannot run
    bool begin() { return algo_.begin_solve(); }

    // The Newton's setup, on the driver's private Ybus (whose values the phase shifters and
    // the shunts patch). On a rebuild of the topology only, the reservations are cleared and
    // `declare()` is called right before the system claims its rows / columns: the loops
    // reserve their controls (controls()) and return the rest of what they need.
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
        return algo_.setup(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, need_init,
                           [&](){ controls_.clear_reservations(); _reserve(declare()); });
    }

    bool newton(int max_iter, real_type tol, bool need_init) { return algo_.newton(max_iter, tol, need_init); }
    void finalize() { algo_.finalize(); }
    // what the loops reserve changed: rebuild the sparsity on the next setup
    void request_rebuild() { algo_.request_rebuild(); }

    // Push what the loops changed to the system, before a solve.
    void apply_state(OuterState & state);

    // the solver's side of a loop's context; `controller_q` is the driver's buffer the
    // context points to
    void fill_context(OuterContext & ctx, RealVect & controller_q) const {
        controller_q = algo_.system().controller_q();
        ctx.controller_q = &controller_q;
        ctx.slack_absorbed = algo_.system().slack_absorbed();
        ctx.branch_control = algo_.system().branch_control();
        ctx.shunt_control = algo_.system().shunt_control();
    }

    // per solver bus, whether its magnitude is an unknown of the last solve; a pinned bus
    // (its Q row pinned, so PV) is not
    std::vector<bool> vm_unknown() const {
        const std::vector<int> & vm_col = algo_.system().vm_to_J_col();
        std::vector<bool> res(vm_col.size(), false);
        for(std::size_t bus = 0; bus < vm_col.size(); ++bus) res[bus] = vm_col[bus] >= 0;
        for(int bus : pinned_) {
            if(bus >= 0 && static_cast<std::size_t>(bus) < res.size()) res[static_cast<std::size_t>(bus)] = false;
        }
        return res;
    }

    // where the loops left the phase / ratio taps and the shunt sections, see
    // BaseAlgo::get_outer_phase_tap and friends
    void phase_tap(std::vector<int> & positions) const {
        positions.clear();
        const BranchControl * phase = algo_.system().branch_control();
        if(phase != nullptr) phase->positions(positions);  // TAP_KEEP for the transformers it does not handle
    }
    void ratio_tap(std::vector<int> & positions) const {
        positions.clear();
        const BranchControl * branch = algo_.system().branch_control();
        if(branch != nullptr) branch->ratio_positions(positions);
    }
    void shunt_sections(std::vector<int> & counts) const {
        counts.clear();
        const ShuntControl * shunt = algo_.system().shunt_control();
        if(shunt != nullptr) shunt->section_counts(counts);
    }

private:
    // what the loops reserved, in the system
    void _reserve(const OuterDeclaration & decl);

    Algo algo_;
    OuterControls controls_;
    std::vector<int> pinned_;  // the buses pinned PV in the last solve
};

template<class LinearSolver, class NRSystem>
void NROuterInner<LinearSolver, NRSystem>::_reserve(const OuterDeclaration & decl)
{
    NRSystem & system = algo_.system();
    // the switchable buses: a caller's (a batch), then the loops'
    algo_.set_switchable_vm_buses(controls_.switchable_buses());  // a set there: duplicates are fine
    system.set_may_hold_voltage_controllers(decl.holds_voltage_controllers());
    system.set_phase_controllers(decl.phase_shifters(), decl.phase_shifter_column());
    std::vector<BranchControl::RatioGroupDecl> groups;
    for(const auto & g : decl.ratio_groups()) {
        BranchControl::RatioGroupDecl d;
        d.bus_solver = g.bus_solver;
        d.target_vm = g.target_vm;
        d.trafos = g.trafos;
        d.solved = g.solved;
        groups.push_back(d);
    }
    system.set_ratio_groups(groups);
    std::vector<ShuntControl::GroupDecl> shunt_groups;
    for(const auto & g : decl.shunt_groups()) {
        ShuntControl::GroupDecl d;
        d.bus_solver = g.bus_solver;
        d.target_vm = g.target_vm;
        d.controller_buses = g.controller_buses;
        d.shunts = g.shunts;
        d.solved = g.solved;
        shunt_groups.push_back(d);
    }
    system.set_shunt_groups(shunt_groups);
}

template<class LinearSolver, class NRSystem>
void NROuterInner<LinearSolver, NRSystem>::apply_state(OuterState & state)
{
    NRSystem & system = algo_.system();
    // what the loops changed outside the injection: the hvdc lines' regimes
    if(!state.hvdc_status.empty()) {
        system.set_hvdc_status_override(state.hvdc_status, OuterState::HVDC_KEEP);
    }
    // ... and the idle standby SVCs switched on (a release is never undone in a solve)
    if(!state.svc_target_vm.empty()) {
        system.release_held_svcs(state.svc_target_vm);
    }
    // the voltage controllers a loop holds (none: the plan's own state)
    system.set_held_voltage_controllers(state.controller_hold_q);
    // the phase taps a loop moved, then the shifts it lets the Newton solve for
    BranchControl * phase = system.branch_control();
    if(phase != nullptr) {
        for(std::size_t t = 0; t < state.phase_tap.size(); ++t) {
            if(state.phase_tap[t] != OuterState::TAP_KEEP) phase->set_tap(static_cast<int>(t), state.phase_tap[t]);
        }
        for(std::size_t t = 0; t < state.phase_control.size(); ++t) {
            if(state.phase_control[t] >= 0) phase->set_control_on(static_cast<int>(t), state.phase_control[t] == 1);
        }
        // the ratio taps (a move once), then the voltage controls
        for(std::size_t t = 0; t < state.ratio_tap.size(); ++t) {
            if(state.ratio_tap[t] == OuterState::TAP_KEEP) continue;
            phase->set_ratio_tap(static_cast<int>(t), state.ratio_tap[t]);
            state.ratio_tap[t] = OuterState::TAP_KEEP;
        }
        for(std::size_t t = 0; t < state.ratio_control.size(); ++t) {
            if(state.ratio_control[t] >= 0) phase->set_ratio_control_on(static_cast<int>(t), state.ratio_control[t] == 1);
        }
    }
    // the shunt sections (a switch once), then the voltage controls
    ShuntControl * shunt = system.shunt_control();
    if(shunt != nullptr) {
        for(const auto & bs : state.shunt_sections) shunt->set_sections(bs.first, bs.second);
        state.shunt_sections.clear();
        for(std::size_t b = 0; b < state.shunt_control.size(); ++b) {
            if(state.shunt_control[b] >= 0) shunt->set_control_on(static_cast<int>(b), state.shunt_control[b] == 1);
        }
    }
    // the switchable buses: PV (pinned) unless a loop made them PQ, a caller's on top
    pinned_ = controls_.pinned_buses();
    system.set_pv_pinned_buses(pinned_);
    // the magnitudes a loop reset (a bus back to PV at its set-point, robust mode)
    std::vector<int> buses;
    std::vector<real_type> vm;
    controls_.take_pending_vm(buses, vm);
    if(!buses.empty()) system.set_vm_at(buses, vm);
}

}  // namespace ls2g

#endif  // NR_OUTER_INNER_H
