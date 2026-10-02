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

void ShuntControl::update_state(const Base * /*nr_system_base_ptr*/, const LSGrid * lsgrid_ptr,
                                const EigenRefConstCplxSpMat & /*Ybus*/, const Eigen::Ref<const CplxVect> & /*Sbus*/,
                                const Eigen::Ref<const RealVect> & /*slack_weights*/)
{
    lsgrid_ = lsgrid_ptr;
    // a new solve: the grid's sections again, the susceptance Ybus was built with
    for (Entry & e : entries_) _reset(e);
}

void ShuntControl::init_topology(const Eigen::Ref<const IntVect> & /*slack_ids*/,
                                 const Eigen::Ref<const RealVect> & /*slack_weights*/,
                                 const Eigen::Ref<const IntVect> & /*pv*/,
                                 const Eigen::Ref<const IntVect> & /*pq*/)
{
    entries_.clear();
    groups_.clear();
    index_of_bus_.clear();
    if (lsgrid_ == nullptr || ybus_ == nullptr || declared_.empty()) return;
    const ShuntContainer & shunts = lsgrid_->get_shunts();
    nb_shunts_ = shunts.nb();
    index_of_bus_.assign(static_cast<std::size_t>(ybus_->rows()), -1);
    for (const GroupDecl & decl : declared_) {
        Group g;
        g.bus = decl.bus_solver;
        g.target = decl.target_vm;
        for (std::size_t c = 0; c < decl.controller_buses.size(); ++c) {
            const int bus = decl.controller_buses[c];
            if (bus < 0 || bus >= ybus_->rows() || index_of_bus_[static_cast<std::size_t>(bus)] >= 0) continue;
            Entry e;
            e.bus = bus;
            e.k = Base::find_J_pos(ybus_->outerIndexPtr(), ybus_->innerIndexPtr(), bus, bus);
            if (e.k < 0) continue;
            for (int s : decl.shunts[c]) {
                if (s >= 0 && s < shunts.nb() && shunts.has_sections(s)) e.shunts.push_back(s);
            }
            if (e.shunts.empty()) continue;
            e.column = decl.solved && g.bus >= 0;
            index_of_bus_[static_cast<std::size_t>(bus)] = static_cast<int>(entries_.size());
            if (e.column) g.members.push_back(static_cast<int>(entries_.size()));
            entries_.push_back(e);
        }
        if (!g.members.empty()) groups_.push_back(g);
    }
    for (Entry & e : entries_) _reset(e);
}

void ShuntControl::_reset(Entry & e)
{
    const ShuntContainer & shunts = lsgrid_->get_shunts();
    const real_type sn = lsgrid_->get_sn_mva();
    e.counts.clear();
    e.b = 0.;
    for (int s : e.shunts) {
        e.counts.push_back(shunts.get_section_count(s));
        // the shunt's own q (MVar at 1 pu, positive absorbing): what Ybus holds now
        if (shunts.get_status()[static_cast<std::size_t>(s)]) e.b -= shunts.get_target_q()(s) / sn;
    }
    e.b_applied = e.b;
    e.b_target = e.b;
    e.on = false;
}

void ShuntControl::_patch(Entry & e)
{
    if (ybus_ != nullptr) ybus_->valuePtr()[e.k] += cplx_type(0., e.b - e.b_applied);
    e.b_applied = e.b;
}

void ShuntControl::register_in(NRLedger & ledger)
{
    for (Entry & e : entries_) {
        if (!e.column) continue;
        e.col = ledger.add_custom_col();
        e.q_row = ledger.q_row(e.bus);
    }
    for (Group & g : groups_) {
        g.vm_col = ledger.vm_col(g.bus);
        g.rows.clear();
        for (std::size_t j = 0; j < g.members.size(); ++j) g.rows.push_back(ledger.add_custom_row());
    }
}

void ShuntControl::declare_feature_entries(FeatureSink & sink)
{
    for (Entry & e : entries_) {
        if (e.column && e.q_row >= 0) e.h_q = sink.add(e.q_row, e.col);
    }
    for (Group & g : groups_) {
        const std::size_t n = g.members.size();
        g.h_cols.assign(n, std::vector<int>(n, -1));
        g.h_vm.assign(n, -1);
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t m = 0; m < n; ++m) {
                g.h_cols[j][m] = sink.add(g.rows[j], entries_[static_cast<std::size_t>(g.members[m])].col);
            }
            if (g.vm_col >= 0) g.h_vm[j] = sink.add(g.rows[j], g.vm_col);
        }
    }
}

void ShuntControl::fill_feature_values(FeatureWriter & writer, const Eigen::Ref<const RealVect> & /*Va*/) const
{
    if (V_ == nullptr) return;
    for (const Entry & e : entries_) {
        // S = V conj(j B V) = -j B |V|^2
        if (e.h_q >= 0) writer.add(e.h_q, -std::norm((*V_)(e.bus)));
    }
    for (const Group & g : groups_) {
        const std::size_t n = g.members.size();
        const real_type inv_n = 1. / static_cast<real_type>(n);
        bool first_on = true;
        for (std::size_t j = 0; j < n; ++j) {
            const Entry & e = entries_[static_cast<std::size_t>(g.members[j])];
            if (!e.on) {
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

void ShuntControl::adjust_mismatch(const Eigen::Ref<const CplxVect> & V_t, const Eigen::Ref<const RealVect> & dx,
                                   Eigen::Ref<CplxVect> mis) const
{
    for (const Entry & e : entries_) {
        if (!e.column || e.col >= dx.size() || dx(e.col) == 0.) continue;
        mis(e.bus) += cplx_type(0., -dx(e.col) * std::norm(V_t(e.bus)));
    }
}

void ShuntControl::fill_custom_rows(Eigen::Ref<RealVect> res, const Eigen::Ref<const RealVect> & /*Va*/,
                                    const Eigen::Ref<const RealVect> & Vm, const Eigen::Ref<const RealVect> & dx) const
{
    for (const Group & g : groups_) {
        const std::size_t n = g.members.size();
        real_type mean = 0.;
        for (int i : g.members) {
            const Entry & e = entries_[static_cast<std::size_t>(i)];
            mean += e.b + dx(e.col);
        }
        mean /= static_cast<real_type>(n);
        bool first_on = true;
        for (std::size_t j = 0; j < n; ++j) {
            const Entry & e = entries_[static_cast<std::size_t>(g.members[j])];
            const real_type b = e.b + dx(e.col);
            if (!e.on) {
                res(g.rows[j]) -= b - e.b_target;
            } else if (first_on) {
                first_on = false;
                const real_type vm = Vm(g.bus) + (g.vm_col >= 0 ? dx(g.vm_col) : 0.);
                res(g.rows[j]) -= vm - g.target;
            } else {
                res(g.rows[j]) -= mean - b;
            }
        }
    }
}

void ShuntControl::apply_step(const Eigen::Ref<const RealVect> & dx)
{
    for (Entry & e : entries_) {
        if (!e.column || dx(e.col) == 0.) continue;
        e.b += dx(e.col);
        _patch(e);
    }
}

void ShuntControl::clear()
{
    entries_.clear();
    groups_.clear();
    index_of_bus_.clear();
}

void ShuntControl::set_control_on(int bus, bool on)
{
    const int k = _index(bus);
    if (k < 0) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    if (!e.column) return;
    e.on = on;
    if (!on) e.b_target = e.b;
}

void ShuntControl::set_sections(int bus, const std::vector<int> & counts)
{
    const int k = _index(bus);
    if (k < 0 || lsgrid_ == nullptr) return;
    Entry & e = entries_[static_cast<std::size_t>(k)];
    if (counts.size() != e.shunts.size()) return;
    const ShuntContainer & shunts = lsgrid_->get_shunts();
    const real_type sn = lsgrid_->get_sn_mva();
    real_type b = 0.;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        const int s = e.shunts[i];
        const int count = std::max(0, std::min(counts[i], shunts.get_max_section_count(s)));
        e.counts[i] = count;
        if (shunts.get_status()[static_cast<std::size_t>(s)]) b -= shunts.section_q(s, count) / sn;
    }
    // the conductance stays the one the sections had, as OpenLoadFlow's dispatchB
    e.b = b;
    e.b_target = b;
    _patch(e);
}

real_type ShuntControl::b(int bus) const
{
    const int k = _index(bus);
    return k < 0 ? std::numeric_limits<real_type>::quiet_NaN() : entries_[static_cast<std::size_t>(k)].b;
}

void ShuntControl::section_counts(std::vector<int> & out) const
{
    out.clear();
    if (entries_.empty()) return;
    out.assign(static_cast<std::size_t>(nb_shunts_), std::numeric_limits<int>::min());
    for (const Entry & e : entries_) {
        for (std::size_t i = 0; i < e.shunts.size(); ++i) out[static_cast<std::size_t>(e.shunts[i])] = e.counts[i];
    }
}

}  // namespace ls2g
