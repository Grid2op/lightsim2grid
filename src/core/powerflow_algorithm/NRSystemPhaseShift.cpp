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

void PhaseShift::update_state(const Base * /*nr_system_base_ptr*/, const LSGrid * lsgrid_ptr,
                              const EigenRefConstCplxSpMat & /*Ybus*/, const Eigen::Ref<const CplxVect> & /*Sbus*/,
                              const Eigen::Ref<const RealVect> & /*slack_weights*/)
{
    lsgrid_ = lsgrid_ptr;
    // a new solve: the grid's taps and shifts again, the block Ybus was built with
    for (Entry & e : entries_) _reset(e);
}

void PhaseShift::init_topology(const Eigen::Ref<const IntVect> & /*slack_ids*/,
                               const Eigen::Ref<const RealVect> & /*slack_weights*/,
                               const Eigen::Ref<const IntVect> & /*pv*/,
                               const Eigen::Ref<const IntVect> & /*pq*/)
{
    entries_.clear();
    index_of_trafo_.clear();
    if (lsgrid_ == nullptr || ybus_ == nullptr || declared_.empty()) return;
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    const SolverBusIdVect & to_solver = lsgrid_->id_me_to_ac_solver();
    index_of_trafo_.assign(static_cast<std::size_t>(trafos.nb()), -1);
    const std::vector<bool> & status = trafos.get_status_global();
    const std::vector<bool> & status1 = trafos.get_status_side_1();
    const std::vector<bool> & status2 = trafos.get_status_side_2();
    for (std::size_t k = 0; k < declared_.size(); ++k) {
        const int t = declared_[k];
        if (t < 0 || t >= trafos.nb() || index_of_trafo_[static_cast<std::size_t>(t)] >= 0) continue;
        // connected at both ends, on two solver buses
        if (!status[t] || !status1[t] || !status2[t]) continue;
        const int b1 = to_solver[trafos.get_bus_side_1(t).cast_int()].cast_int();
        const int b2 = to_solver[trafos.get_bus_side_2(t).cast_int()].cast_int();
        if (b1 < 0 || b2 < 0 || b1 == b2) continue;
        Entry e;
        e.trafo = t;
        e.column = k < declared_column_.size() && declared_column_[k] != 0;
        e.b1 = b1;
        e.b2 = b2;
        e.k11 = ybus_position(*ybus_, b1, b1);
        e.k12 = ybus_position(*ybus_, b1, b2);
        e.k21 = ybus_position(*ybus_, b2, b1);
        e.k22 = ybus_position(*ybus_, b2, b2);
        if (e.k11 < 0 || e.k12 < 0 || e.k21 < 0 || e.k22 < 0) continue;
        index_of_trafo_[static_cast<std::size_t>(t)] = static_cast<int>(entries_.size());
        entries_.push_back(e);
    }
    for (Entry & e : entries_) _reset(e);
}

void PhaseShift::_reset(Entry & e)
{
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    const TrafoInfo info = trafos[e.trafo];
    e.sign = trafos.shift_sign(e.trafo);
    e.pos = ptc.has(e.trafo) ? ptc.position(e.trafo) : 0;
    e.a_tap = info.shift_rad;
    e.a = e.a_tap;
    e.a_target = e.a_tap;
    e.control_on = false;
    e.side = ptc.has(e.trafo) && ptc.regulated(e.trafo) == 2 ? 2 : 1;
    e.p_target = ptc.has(e.trafo) ? ptc.target(e.trafo) / lsgrid_->get_sn_mva() : 0.;
    e.tap = trafos.pi_block_at(e.trafo, e.pos, e.a_tap);
    // Ybus holds what the container stamped (a fresh copy of the grid's, see NROuterAlgo)
    e.applied = {info.yac_eff_11, info.yac_eff_12, info.yac_eff_21, info.yac_eff_22};
    _patch(e);
}

void PhaseShift::_patch(Entry & e)
{
    const std::array<cplx_type, 4> block = _block_at(e, e.a);
    if (ybus_ != nullptr) {
        cplx_type * values = ybus_->valuePtr();
        values[e.k11] += block[0] - e.applied[0];
        values[e.k12] += block[1] - e.applied[1];
        values[e.k21] += block[2] - e.applied[2];
        values[e.k22] += block[3] - e.applied[3];
    }
    e.applied = block;
}

void PhaseShift::register_in(NRLedger & ledger)
{
    for (Entry & e : entries_) {
        if (!e.column) continue;
        e.col = ledger.add_custom_col();
        e.row = ledger.add_custom_row();
        e.p1 = ledger.p_row(e.b1);
        e.q1 = ledger.q_row(e.b1);
        e.p2 = ledger.p_row(e.b2);
        e.q2 = ledger.q_row(e.b2);
        e.th1 = ledger.theta_col(e.b1);
        e.th2 = ledger.theta_col(e.b2);
        e.vm1 = ledger.vm_col(e.b1);
        e.vm2 = ledger.vm_col(e.b2);
    }
}

void PhaseShift::declare_feature_entries(FeatureSink & sink)
{
    for (Entry & e : entries_) {
        if (!e.column) continue;
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
}

void PhaseShift::fill_feature_values(FeatureWriter & writer, const Eigen::Ref<const RealVect> & /*Va*/) const
{
    if (V_ == nullptr) return;
    for (const Entry & e : entries_) {
        if (!e.column) continue;
        const cplx_type V1 = (*V_)(e.b1);
        const cplx_type V2 = (*V_)(e.b2);
        const std::array<cplx_type, 4> & y = e.applied;
        // S1 = V1 conj(y11 V1 + y12 V2), S2 = V2 conj(y21 V1 + y22 V2); dy12/da = j s y12,
        // dy21/da = -j s y21
        const cplx_type v1_y12_v2 = V1 * std::conj(y[1] * V2);
        const cplx_type v2_y21_v1 = V2 * std::conj(y[2] * V1);
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
}

void PhaseShift::adjust_mismatch(const Eigen::Ref<const CplxVect> & V_t, const Eigen::Ref<const RealVect> & dx,
                                 Eigen::Ref<CplxVect> mis) const
{
    // Ybus holds the block at the current shift: a trial step on the shift adds the difference
    for (const Entry & e : entries_) {
        if (!e.column || e.col >= dx.size() || dx(e.col) == 0.) continue;
        const std::array<cplx_type, 4> block = _block_at(e, e.a + dx(e.col));
        mis(e.b1) += V_t(e.b1) * std::conj((block[1] - e.applied[1]) * V_t(e.b2));
        mis(e.b2) += V_t(e.b2) * std::conj((block[2] - e.applied[2]) * V_t(e.b1));
    }
}

void PhaseShift::fill_custom_rows(Eigen::Ref<RealVect> res, const Eigen::Ref<const RealVect> & Va,
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
        const std::array<cplx_type, 4> y = _block_at(e, a);
        const real_type p = e.side == 1 ? std::real(V1 * std::conj(y[0] * V1 + y[1] * V2))
                                        : std::real(V2 * std::conj(y[2] * V1 + y[3] * V2));
        res(e.row) -= p - e.p_target;
    }
}

void PhaseShift::apply_step(const Eigen::Ref<const RealVect> & dx)
{
    for (Entry & e : entries_) {
        if (!e.column || dx(e.col) == 0.) continue;
        e.a += dx(e.col);
        _patch(e);
    }
}

void PhaseShift::clear()
{
    entries_.clear();
    index_of_trafo_.clear();
}

void PhaseShift::set_control_on(int trafo_id, bool on)
{
    const int k = _index(trafo_id);
    if (k < 0) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    if (!e.column) return;
    e.control_on = on;
    if (!on) e.a_target = e.a;
}

void PhaseShift::set_tap(int trafo_id, int position)
{
    const int k = _index(trafo_id);
    if (k < 0 || lsgrid_ == nullptr) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    const TrafoContainer & trafos = lsgrid_->get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    if (!ptc.has(trafo_id) || position < ptc.low_tap(trafo_id) || position > ptc.high_tap(trafo_id)) return;
    if (position == e.pos && e.a == e.a_tap) return;
    e.pos = position;
    e.a_tap = ptc.alpha_at(trafo_id, position);
    e.tap = trafos.pi_block_at(trafo_id, position, e.a_tap);
    e.a = e.a_tap;
    e.a_target = e.a_tap;
    _patch(e);
}

real_type PhaseShift::shift(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? std::numeric_limits<real_type>::quiet_NaN() : entries_[static_cast<std::size_t>(k)].a;
}

int PhaseShift::position(int trafo_id) const
{
    const int k = _index(trafo_id);
    return k < 0 ? -1 : entries_[static_cast<std::size_t>(k)].pos;
}

void PhaseShift::positions(std::vector<int> & out) const
{
    out.clear();
    if (entries_.empty()) return;
    out.assign(index_of_trafo_.size(), std::numeric_limits<int>::min());
    for (const Entry & e : entries_) out[static_cast<std::size_t>(e.trafo)] = e.pos;
}

void PhaseShift::current(int trafo_id, int side, real_type & i_pu, real_type & di_da) const
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
