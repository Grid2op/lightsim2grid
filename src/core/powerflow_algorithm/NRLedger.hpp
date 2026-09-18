// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_LEDGER_H
#define NR_LEDGER_H

#include "Utils.hpp"
#include "ls2g_api.hpp"

#include <cassert>
#include <vector>

namespace ls2g {

/**
 * Central registry of the equations (Jacobian rows) and unknowns (Jacobian
 * columns) of a Newton-Raphson system.
 *
 * Each component of an NRSystem (the Base block and every extension) claims
 * its rows / columns here during topology initialisation, in registration
 * order: the allocation order defines the layout of the augmented Jacobian.
 *
 * Rows / columns can be:
 *   - bus-owned: a P or Q power-mismatch equation at a bus (row), or a
 *     theta / vm voltage unknown at a bus (column). These are recorded in
 *     bus-keyed maps (size n_bus, -1 when absent) AND in compact
 *     (bus, index) pair lists for fast iteration.
 *   - custom: an index with no bus attached (e.g. the slack_absorbed
 *     unknown, a voltage-setpoint equation).
 *
 * Multiplicity rules (no injectivity assumption):
 *   - a bus has at most one entry in each bus-keyed map (last registration
 *     wins; current use never re-registers);
 *   - several buses MAY claim the same row, and several buses MAY share the
 *     same column;
 *   - several Jacobian contributions MAY resolve to the same J position, so
 *     ALL numeric fills must use accumulate (+=) semantics.
 */
class LS2G_API NRLedger
{
    public:
        NRLedger(): n_rows_(0), n_cols_(0) {}

        void reset(int n_bus) {
            p_row_of_bus_.assign(n_bus, -1);
            q_row_of_bus_.assign(n_bus, -1);
            theta_col_of_bus_.assign(n_bus, -1);
            vm_col_of_bus_.assign(n_bus, -1);
            q_col_of_bus_.assign(n_bus, -1);
            p_buses_.clear();     p_rows_.clear();
            q_buses_.clear();     q_rows_.clear();
            theta_buses_.clear(); theta_cols_.clear();
            vm_buses_.clear();    vm_cols_.clear();
            n_rows_ = 0;
            n_cols_ = 0;
        }

        // ---- registration (allocation order defines the J layout) -----------
        int add_p_equation(int bus) {
            const int row = n_rows_++;
            p_row_of_bus_[bus] = row;
            p_buses_.push_back(bus);
            p_rows_.push_back(row);
            return row;
        }
        int add_q_equation(int bus) {
            const int row = n_rows_++;
            q_row_of_bus_[bus] = row;
            q_buses_.push_back(bus);
            q_rows_.push_back(row);
            return row;
        }
        int add_theta_unknown(int bus) {
            const int col = n_cols_++;
            theta_col_of_bus_[bus] = col;
            theta_buses_.push_back(bus);
            theta_cols_.push_back(col);
            return col;
        }
        int add_vm_unknown(int bus) {
            const int col = n_cols_++;
            vm_col_of_bus_[bus] = col;
            vm_buses_.push_back(bus);
            vm_cols_.push_back(col);
            return col;
        }
        int add_custom_row() { return n_rows_++; }
        int add_custom_col() { return n_cols_++; }
        // sugar: a custom column recorded in q_col_of_bus_ for Python-facing
        // introspection (q_to_J_col). Nothing else.
        int add_q_unknown(int bus) {
            const int col = add_custom_col();
            q_col_of_bus_[bus] = col;
            return col;
        }

        // ---- lookups ---------------------------------------------------------
        int p_row(int bus)     const { return p_row_of_bus_[bus]; }
        int q_row(int bus)     const { return q_row_of_bus_[bus]; }
        int theta_col(int bus) const { return theta_col_of_bus_[bus]; }
        int vm_col(int bus)    const { return vm_col_of_bus_[bus]; }

        const std::vector<int>& p_buses()     const { return p_buses_; }
        const std::vector<int>& p_rows()      const { return p_rows_; }
        const std::vector<int>& q_buses()     const { return q_buses_; }
        const std::vector<int>& q_rows()      const { return q_rows_; }
        const std::vector<int>& theta_buses() const { return theta_buses_; }
        const std::vector<int>& theta_cols()  const { return theta_cols_; }
        const std::vector<int>& vm_buses()    const { return vm_buses_; }
        const std::vector<int>& vm_cols()     const { return vm_cols_; }

        // bus_id -> Jacobian column converters (-1 when the bus owns no such
        // unknown); exposed as-is through NRSystem's public API.
        const std::vector<int>& theta_col_of_bus() const { return theta_col_of_bus_; }
        const std::vector<int>& vm_col_of_bus()    const { return vm_col_of_bus_; }
        const std::vector<int>& q_col_of_bus()     const { return q_col_of_bus_; }

        // bus_id -> Jacobian row converters (-1 when the bus owns no such
        // equation); the row counterpart of the *_col_of_bus maps. Consumed by
        // external solvers (e.g. gpusim2grid) that re-derive the dS scatter and
        // residual layout of the augmented J from the ledger.
        const std::vector<int>& p_row_of_bus() const { return p_row_of_bus_; }
        const std::vector<int>& q_row_of_bus() const { return q_row_of_bus_; }

        int size() const {
            assert(n_rows_ == n_cols_);
            return n_rows_;
        }

    private:
        // bus-keyed maps (size n_bus, -1 when absent)
        std::vector<int> p_row_of_bus_, q_row_of_bus_;
        std::vector<int> theta_col_of_bus_, vm_col_of_bus_, q_col_of_bus_;

        // compact pair lists for fast iteration
        std::vector<int> p_buses_, p_rows_;
        std::vector<int> q_buses_, q_rows_;
        std::vector<int> theta_buses_, theta_cols_;
        std::vector<int> vm_buses_, vm_cols_;

        int n_rows_, n_cols_;
};

/**
 * Collects the feature-specific (non dS-derived) Jacobian entries declared by
 * the components during the sparsity build. Each add(row, col) returns a
 * handle the component keeps; once the J pattern is compressed, NRSystem
 * resolves each (row, col) to a position in J.valuePtr().
 */
class LS2G_API FeatureSink
{
    public:
        void clear() { rows_.clear(); cols_.clear(); }
        int add(int row, int col) {
            const int h = static_cast<int>(rows_.size());
            rows_.push_back(row);
            cols_.push_back(col);
            return h;
        }
        int size()      const { return static_cast<int>(rows_.size()); }
        int row(int h)  const { return rows_[h]; }
        int col(int h)  const { return cols_[h]; }

    private:
        std::vector<int> rows_, cols_;
};

/**
 * Write access to the feature entries during fill_J. J has been zeroed
 * beforehand, so writes accumulate (several contributions may share a
 * position).
 */
class LS2G_API FeatureWriter
{
    public:
        FeatureWriter(real_type* J_values, const std::vector<int>& positions):
            J_values_(J_values),
            positions_(positions) {}

        void add(int h, real_type value) const { J_values_[positions_[h]] += value; }

    private:
        real_type*              J_values_;
        const std::vector<int>& positions_;
};

} // namespace ls2g

#endif // NR_LEDGER_H
