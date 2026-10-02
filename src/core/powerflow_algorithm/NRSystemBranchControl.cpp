// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

#include "LSGrid.hpp"

#include <limits>

namespace ls2g {

namespace {

// position of (row, col) in a column-major sparse matrix' values, -1 if not stored
int ybus_position(const Eigen::SparseMatrix<cplx_type> & ybus, int row, int col)
{
    return Base::find_J_pos(ybus.outerIndexPtr(), ybus.innerIndexPtr(), row, col);
}

const cplx_type J_(0., 1.);

}  // namespace

void BranchControl::update_state(const Base * /*nr_system_base_ptr*/, const LSGrid * lsgrid_ptr,
                                 const EigenRefConstCplxSpMat & /*Ybus*/, const Eigen::Ref<const CplxVect> & /*Sbus*/,
                                 const Eigen::Ref<const RealVect> & /*slack_weights*/)
{
    lsgrid_ = lsgrid_ptr;
    // a new solve: the grid's taps, shifts and ratios again, the block Ybus was built with
    for (Entry & e : entries_) _reset(e);
}

int BranchControl::_add_entry(int t)
{
    const int known = _index(t);
    if (known >= 0) return known;
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    if (t < 0 || t >= trafos.nb()) return -1;
    // connected at both ends, on two solver buses
    if (!trafos.get_status_global()[t] || !trafos.get_status_side_1()[t] || !trafos.get_status_side_2()[t]) return -1;
    const SolverBusIdVect & to_solver = lsgrid_->id_me_to_ac_solver();
    const int b1 = to_solver[trafos.get_bus_side_1(t).cast_int()].cast_int();
    const int b2 = to_solver[trafos.get_bus_side_2(t).cast_int()].cast_int();
    if (b1 < 0 || b2 < 0 || b1 == b2) return -1;
    Entry e;
    e.trafo = t;
    e.b1 = b1;
    e.b2 = b2;
    e.k11 = ybus_position(*ybus_, b1, b1);
    e.k12 = ybus_position(*ybus_, b1, b2);
    e.k21 = ybus_position(*ybus_, b2, b1);
    e.k22 = ybus_position(*ybus_, b2, b2);
    if (e.k11 < 0 || e.k12 < 0 || e.k21 < 0 || e.k22 < 0) return -1;
    const int k = static_cast<int>(entries_.size());
    index_of_trafo_[static_cast<std::size_t>(t)] = k;
    entries_.push_back(e);
    return k;
}

void BranchControl::init_topology(const Eigen::Ref<const IntVect> & /*slack_ids*/,
                                  const Eigen::Ref<const RealVect> & /*slack_weights*/,
                                  const Eigen::Ref<const IntVect> & /*pv*/,
                                  const Eigen::Ref<const IntVect> & /*pq*/)
{
    entries_.clear();
    groups_.clear();
    index_of_trafo_.clear();
    if (lsgrid_ == nullptr || ybus_ == nullptr || (declared_.empty() && declared_groups_.empty())) return;
    index_of_trafo_.assign(static_cast<std::size_t>(lsgrid_->get_trafos().nb()), -1);
    for (std::size_t k = 0; k < declared_.size(); ++k) {
        const int i = _add_entry(declared_[k]);
        if (i < 0) continue;
        Entry & e = entries_[static_cast<std::size_t>(i)];
        e.phase = true;
        e.column = e.column || (k < declared_column_.size() && declared_column_[k] != 0);
    }
    for (const RatioGroupDecl & decl : declared_groups_) {
        RatioGroup g;
        g.bus = decl.bus_solver;
        g.target = decl.target_vm;
        for (int t : decl.trafos) {
            const int i = _add_entry(t);
            if (i < 0) continue;
            Entry & e = entries_[static_cast<std::size_t>(i)];
            if (e.ratio) continue;  // already in a group
            e.ratio = true;
            // a phase controller keeps its column alone (see the class comment)
            if (decl.solved && !e.column && g.bus >= 0) {
                e.r_column = true;
                g.members.push_back(i);
            }
        }
        if (!g.members.empty()) groups_.push_back(g);
    }
    for (Entry & e : entries_) _reset(e);
}

void BranchControl::_reset(Entry & e)
{
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    const TapChangers & rtc = trafos.get_tap_changers(false);
    const TrafoInfo info = trafos[e.trafo];
    e.sign = trafos.shift_sign(e.trafo);
    e.rexp = trafos.ratio_exponent(e.trafo);
    e.control_on = false;
    e.r_control_on = false;
    e.side = ptc.has(e.trafo) && ptc.regulated(e.trafo) == 2 ? 2 : 1;
    e.p_target = ptc.has(e.trafo) ? ptc.target(e.trafo) / lsgrid_->get_sn_mva() : 0.;
    e.pos = ptc.has(e.trafo) ? ptc.position(e.trafo) : 0;
    e.rpos = rtc.has(e.trafo) ? rtc.position(e.trafo) : 0;
    e.a_tap = info.shift_rad;
    e.rho_tap = (ptc.has(e.trafo) || rtc.has(e.trafo)) ? trafos.ratio_at(e.trafo, e.rpos, e.pos) : info.ratio;
    e.a = e.a_tap;
    e.a_target = e.a_tap;
    e.rho = e.rho_tap;
    e.rho_target = e.rho_tap;
    e.tap = trafos.pi_block_at(e.trafo, e.rpos, e.pos, e.a_tap);
    // Ybus holds what the container stamped (a fresh copy of the grid's, see NROuterAlgo)
    e.applied = {info.yac_eff_11, info.yac_eff_12, info.yac_eff_21, info.yac_eff_22};
    _patch(e);
}

void BranchControl::_retap(Entry & e, int rpos, int pos)
{
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    e.pos = pos;
    e.rpos = rpos;
    if (ptc.has(e.trafo)) e.a_tap = ptc.alpha_at(e.trafo, pos);
    e.rho_tap = trafos.ratio_at(e.trafo, rpos, pos);
    e.tap = trafos.pi_block_at(e.trafo, rpos, pos, e.a_tap);
    e.a = e.a_tap;
    e.a_target = e.a_tap;
    e.rho = e.rho_tap;
    e.rho_target = e.rho_tap;
    _patch(e);
}

void BranchControl::_patch(Entry & e)
{
    const std::array<cplx_type, 4> block = _block_at(e, e.a, e.rho);
    if (ybus_ != nullptr) {
        cplx_type * values = ybus_->valuePtr();
        values[e.k11] += block[0] - e.applied[0];
        values[e.k12] += block[1] - e.applied[1];
        values[e.k21] += block[2] - e.applied[2];
        values[e.k22] += block[3] - e.applied[3];
    }
    e.applied = block;
}

void BranchControl::register_in(NRLedger & ledger)
{
    for (Entry & e : entries_) {
        if (!e.column && !e.r_column) continue;
        if (e.column) {
            e.col = ledger.add_custom_col();
            e.row = ledger.add_custom_row();
        }
        if (e.r_column) e.rcol = ledger.add_custom_col();
        e.p1 = ledger.p_row(e.b1);
        e.q1 = ledger.q_row(e.b1);
        e.p2 = ledger.p_row(e.b2);
        e.q2 = ledger.q_row(e.b2);
        e.th1 = ledger.theta_col(e.b1);
        e.th2 = ledger.theta_col(e.b2);
        e.vm1 = ledger.vm_col(e.b1);
        e.vm2 = ledger.vm_col(e.b2);
    }
    for (RatioGroup & g : groups_) {
        g.vm_col = ledger.vm_col(g.bus);
        g.rows.clear();
        for (std::size_t j = 0; j < g.members.size(); ++j) g.rows.push_back(ledger.add_custom_row());
    }
}

void BranchControl::declare_feature_entries(FeatureSink & sink)
{
    for (Entry & e : entries_) {
        if (e.column) {
            // dS / da on the two buses
            e.h_p1 = e.p1 >= 0 ? sink.add(e.p1, e.col) : -1;
            e.h_q1 = e.q1 >= 0 ? sink.add(e.q1, e.col) : -1;
            e.h_p2 = e.p2 >= 0 ? sink.add(e.p2, e.col) : -1;
            e.h_q2 = e.q2 >= 0 ? sink.add(e.q2, e.col) : -1;
            // the row, both of its forms
            e.h_a = sink.add(e.row, e.col);
            e.h_th1 = e.th1 >= 0 ? sink.add(e.row, e.th1) : -1;
            e.h_th2 = e.th2 >= 0 ? sink.add(e.row, e.th2) : -1;
            e.h_vm1 = e.vm1 >= 0 ? sink.add(e.row, e.vm1) : -1;
            e.h_vm2 = e.vm2 >= 0 ? sink.add(e.row, e.vm2) : -1;
        }
        if (e.r_column) {
            // dS / drho on the two buses
            e.h_rp1 = e.p1 >= 0 ? sink.add(e.p1, e.rcol) : -1;
            e.h_rq1 = e.q1 >= 0 ? sink.add(e.q1, e.rcol) : -1;
            e.h_rp2 = e.p2 >= 0 ? sink.add(e.p2, e.rcol) : -1;
            e.h_rq2 = e.q2 >= 0 ? sink.add(e.q2, e.rcol) : -1;
        }
    }
    // every slot of a group, every one of its forms: the group's columns and the bus' Vm
    for (RatioGroup & g : groups_) {
        const std::size_t n = g.members.size();
        g.h_cols.assign(n, std::vector<int>(n, -1));
        g.h_vm.assign(n, -1);
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t m = 0; m < n; ++m) {
                g.h_cols[j][m] = sink.add(g.rows[j], entries_[static_cast<std::size_t>(g.members[m])].rcol);
            }
            if (g.vm_col >= 0) g.h_vm[j] = sink.add(g.rows[j], g.vm_col);
        }
    }
}

void BranchControl::fill_feature_values(FeatureWriter & writer, const Eigen::Ref<const RealVect> & /*Va*/) const
{
    if (V_ == nullptr) return;
    for (const Entry & e : entries_) {
        if (!e.column && !e.r_column) continue;
        const cplx_type V1 = (*V_)(e.b1);
        const cplx_type V2 = (*V_)(e.b2);
        const std::array<cplx_type, 4> & y = e.applied;
        const cplx_type v1_y12_v2 = V1 * std::conj(y[1] * V2);
        const cplx_type v2_y21_v1 = V2 * std::conj(y[2] * V1);
        if (e.r_column) {
            // dy11/drho = 2 e y11 / rho, dy12/drho = e y12 / rho, dy21/drho = e y21 / rho
            const real_type k = e.rexp / e.rho;
            const cplx_type dS1_dr = k * (2. * V1 * std::conj(y[0] * V1) + v1_y12_v2);
            const cplx_type dS2_dr = k * v2_y21_v1;
            if (e.h_rp1 >= 0) writer.add(e.h_rp1, std::real(dS1_dr));
            if (e.h_rq1 >= 0) writer.add(e.h_rq1, std::imag(dS1_dr));
            if (e.h_rp2 >= 0) writer.add(e.h_rp2, std::real(dS2_dr));
            if (e.h_rq2 >= 0) writer.add(e.h_rq2, std::imag(dS2_dr));
        }
        if (!e.column) continue;
        // S1 = V1 conj(y11 V1 + y12 V2), S2 = V2 conj(y21 V1 + y22 V2); dy12/da = j s y12,
        // dy21/da = -j s y21
        const cplx_type dS1_da = -J_ * e.sign * v1_y12_v2;
        const cplx_type dS2_da = J_ * e.sign * v2_y21_v1;
        if (e.h_p1 >= 0) writer.add(e.h_p1, std::real(dS1_da));
        if (e.h_q1 >= 0) writer.add(e.h_q1, std::imag(dS1_da));
        if (e.h_p2 >= 0) writer.add(e.h_p2, std::real(dS2_da));
        if (e.h_q2 >= 0) writer.add(e.h_q2, std::imag(dS2_da));
        if (!e.control_on) {
            writer.add(e.h_a, 1.);
            continue;
        }
        const real_type vm1 = std::abs(V1);
        const real_type vm2 = std::abs(V2);
        real_type d_th1, d_th2, d_vm1, d_vm2, d_a;
        if (e.side == 1) {
            const cplx_type S1 = V1 * std::conj(y[0] * V1 + y[1] * V2);
            d_th1 = std::real(J_ * v1_y12_v2);
            d_th2 = -d_th1;
            d_vm1 = std::real((S1 + vm1 * vm1 * std::conj(y[0])) / vm1);
            d_vm2 = std::real(v1_y12_v2) / vm2;
            d_a = std::real(dS1_da);
        } else {
            const cplx_type S2 = V2 * std::conj(y[2] * V1 + y[3] * V2);
            d_th2 = std::real(J_ * v2_y21_v1);
            d_th1 = -d_th2;
            d_vm2 = std::real((S2 + vm2 * vm2 * std::conj(y[3])) / vm2);
            d_vm1 = std::real(v2_y21_v1) / vm1;
            d_a = std::real(dS2_da);
        }
        writer.add(e.h_a, d_a);
        if (e.h_th1 >= 0) writer.add(e.h_th1, d_th1);
        if (e.h_th2 >= 0) writer.add(e.h_th2, d_th2);
        if (e.h_vm1 >= 0) writer.add(e.h_vm1, d_vm1);
        if (e.h_vm2 >= 0) writer.add(e.h_vm2, d_vm2);
    }
    for (const RatioGroup & g : groups_) {
        const std::size_t n = g.members.size();
        const real_type inv_n = 1. / static_cast<real_type>(n);
        bool first_on = true;
        for (std::size_t j = 0; j < n; ++j) {
            const Entry & e = entries_[static_cast<std::size_t>(g.members[j])];
            if (!e.r_control_on) {
                writer.add(g.h_cols[j][j], 1.);
            } else if (first_on) {
                first_on = false;
                if (g.h_vm[j] >= 0) writer.add(g.h_vm[j], 1.);
            } else {
                for (std::size_t m = 0; m < n; ++m) writer.add(g.h_cols[j][m], m == j ? inv_n - 1. : inv_n);
            }
        }
    }
}

void BranchControl::adjust_mismatch(const Eigen::Ref<const CplxVect> & V_t, const Eigen::Ref<const RealVect> & dx,
                                    Eigen::Ref<CplxVect> mis) const
{
    // Ybus holds the block at the current shift and ratio: a trial step on them adds the difference
    for (const Entry & e : entries_) {
        const real_type da = e.column && e.col < dx.size() ? dx(e.col) : 0.;
        const real_type dr = e.r_column && e.rcol < dx.size() ? dx(e.rcol) : 0.;
        if (da == 0. && dr == 0.) continue;
        const std::array<cplx_type, 4> block = _block_at(e, e.a + da, e.rho + dr);
        const cplx_type V1 = V_t(e.b1);
        const cplx_type V2 = V_t(e.b2);
        mis(e.b1) += V1 * std::conj((block[0] - e.applied[0]) * V1 + (block[1] - e.applied[1]) * V2);
        mis(e.b2) += V2 * std::conj((block[2] - e.applied[2]) * V1);
    }
}

void BranchControl::fill_custom_rows(Eigen::Ref<RealVect> res, const Eigen::Ref<const RealVect> & Va,
                                     const Eigen::Ref<const RealVect> & Vm, const Eigen::Ref<const RealVect> & dx) const
{
    for (const Entry & e : entries_) {
        if (!e.column) continue;
        const real_type a = e.a + dx(e.col);
        if (!e.control_on) {
            res(e.row) -= a - e.a_target;
            continue;
        }
        const real_type va1 = Va(e.b1) + (e.th1 >= 0 ? dx(e.th1) : 0.);
        const real_type va2 = Va(e.b2) + (e.th2 >= 0 ? dx(e.th2) : 0.);
        const real_type vm1 = Vm(e.b1) + (e.vm1 >= 0 ? dx(e.vm1) : 0.);
        const real_type vm2 = Vm(e.b2) + (e.vm2 >= 0 ? dx(e.vm2) : 0.);
        const cplx_type V1 = std::polar(vm1, va1);
        const cplx_type V2 = std::polar(vm2, va2);
        const std::array<cplx_type, 4> y = _block_at(e, a, e.rho);
        const real_type p = e.side == 1 ? std::real(V1 * std::conj(y[0] * V1 + y[1] * V2))
                                        : std::real(V2 * std::conj(y[2] * V1 + y[3] * V2));
        res(e.row) -= p - e.p_target;
    }
    for (const RatioGroup & g : groups_) {
        const std::size_t n = g.members.size();
        real_type mean = 0.;
        for (int i : g.members) {
            const Entry & e = entries_[static_cast<std::size_t>(i)];
            mean += e.rho + dx(e.rcol);
        }
        mean /= static_cast<real_type>(n);
        bool first_on = true;
        for (std::size_t j = 0; j < n; ++j) {
            const Entry & e = entries_[static_cast<std::size_t>(g.members[j])];
            const real_type rho = e.rho + dx(e.rcol);
            if (!e.r_control_on) {
                res(g.rows[j]) -= rho - e.rho_target;
            } else if (first_on) {
                first_on = false;
                const real_type vm = Vm(g.bus) + (g.vm_col >= 0 ? dx(g.vm_col) : 0.);
                res(g.rows[j]) -= vm - g.target;
            } else {
                res(g.rows[j]) -= mean - rho;
            }
        }
    }
}

void BranchControl::apply_step(const Eigen::Ref<const RealVect> & dx)
{
    for (Entry & e : entries_) {
        const real_type da = e.column ? dx(e.col) : 0.;
        const real_type dr = e.r_column ? dx(e.rcol) : 0.;
        if (da == 0. && dr == 0.) continue;
        e.a += da;
        e.rho += dr;
        _patch(e);
    }
}

void BranchControl::clear()
{
    entries_.clear();
    groups_.clear();
    index_of_trafo_.clear();
}

void BranchControl::set_control_on(int trafo_id, bool on)
{
    const int k = _index(trafo_id);
    if (k < 0) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    if (!e.column) return;
    e.control_on = on;
    if (!on) e.a_target = e.a;
}

void BranchControl::set_ratio_control_on(int trafo_id, bool on)
{
    const int k = _index(trafo_id);
    if (k < 0) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    if (!e.r_column) return;
    e.r_control_on = on;
    if (!on) e.rho_target = e.rho;
}

void BranchControl::set_tap(int trafo_id, int position)
{
    const int k = _index(trafo_id);
    if (k < 0 || lsgrid_ == nullptr) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    const TapChangers & ptc = lsgrid_->get_trafos().get_tap_changers(true);
    if (!e.phase || !ptc.has(trafo_id) || position < ptc.low_tap(trafo_id) || position > ptc.high_tap(trafo_id)) return;
    if (position == e.pos && e.a == e.a_tap && e.rho == e.rho_tap) return;
    _retap(e, e.rpos, position);
}

void BranchControl::set_ratio_tap(int trafo_id, int position)
{
    const int k = _index(trafo_id);
    if (k < 0 || lsgrid_ == nullptr) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    const TapChangers & rtc = lsgrid_->get_trafos().get_tap_changers(false);
    if (!e.ratio || !rtc.has(trafo_id) || position < rtc.low_tap(trafo_id) || position > rtc.high_tap(trafo_id)) return;
    if (position == e.rpos && e.a == e.a_tap && e.rho == e.rho_tap) return;
    _retap(e, position, e.pos);
}

real_type BranchControl::shift(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? std::numeric_limits<real_type>::quiet_NaN() : entries_[static_cast<std::size_t>(k)].a;
}

real_type BranchControl::ratio(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? std::numeric_limits<real_type>::quiet_NaN() : entries_[static_cast<std::size_t>(k)].rho;
}

int BranchControl::position(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? std::numeric_limits<int>::min() : entries_[static_cast<std::size_t>(k)].pos;
}

int BranchControl::ratio_position(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? std::numeric_limits<int>::min() : entries_[static_cast<std::size_t>(k)].rpos;
}

void BranchControl::positions(std::vector<int> & out) const
{
    out.clear();
    bool any = false;
    for (const Entry & e : entries_) any = any || e.phase;
    if (!any) return;
    out.assign(index_of_trafo_.size(), std::numeric_limits<int>::min());
    for (const Entry & e : entries_) {
        if (e.phase) out[static_cast<std::size_t>(e.trafo)] = e.pos;
    }
}

void BranchControl::ratio_positions(std::vector<int> & out) const
{
    out.clear();
    bool any = false;
    for (const Entry & e : entries_) any = any || e.ratio;
    if (!any) return;
    out.assign(index_of_trafo_.size(), std::numeric_limits<int>::min());
    for (const Entry & e : entries_) {
        if (e.ratio) out[static_cast<std::size_t>(e.trafo)] = e.rpos;
    }
}

void BranchControl::current(int trafo_id, int side, real_type & i_pu, real_type & di_da) const
{
    i_pu = 0.;
    di_da = 0.;
    const int k = _index(trafo_id);
    if (k < 0 || V_ == nullptr) return;
    const Entry & e = entries_[static_cast<std::size_t>(k)];
    const cplx_type V1 = (*V_)(e.b1);
    const cplx_type V2 = (*V_)(e.b2);
    const std::array<cplx_type, 4> & y = e.applied;
    const cplx_type I = side == 1 ? y[0] * V1 + y[1] * V2 : y[2] * V1 + y[3] * V2;
    const cplx_type dI = side == 1 ? J_ * e.sign * y[1] * V2 : -J_ * e.sign * y[2] * V1;
    i_pu = std::abs(I);
    if (i_pu > 0.) di_da = std::real(std::conj(I) * dI) / i_pu;
}

}  // namespace ls2g
