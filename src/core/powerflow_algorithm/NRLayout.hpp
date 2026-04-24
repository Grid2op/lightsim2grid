// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_LAYOUT_H
#define NR_LAYOUT_H

#include "Utils.hpp"
#include "Eigen/Core"

namespace ls2g {

/**
 * Encodes everything about how the NR state vector and Jacobian blocks are indexed.
 *
 *   lag == 1  → multi-slack  (extra row/col for the distributed-slack variable)
 *   lag == 0  → single-slack (no extra row/col)
 *
 * All attributes are private; the entire public interface is const.
 * The attributes are non-const so that the owning class can replace the layout with
 *   _layout = NRLayout(nb_pv, nb_pq, lag);
 * at the start of each compute_pf call.  Making them truly const would delete the
 * copy-assignment operator and break that reassignment pattern, which only becomes
 * formally well-defined in C++17 (via placement-new on trivially-destructible objects).
 */
class NRLayout final {
public:
    NRLayout() noexcept : nb_pv_(0), nb_pq_(0), lag_(0) {}
    NRLayout(int nb_pv, int nb_pq, int lag) noexcept
        : nb_pv_(nb_pv), nb_pq_(nb_pq), lag_(lag) {}

    // --- scalar accessors ---------------------------------------------------
    int nb_pv()      const { return nb_pv_; }
    int nb_pq()      const { return nb_pq_; }
    int lag()        const { return lag_; }
    int theta_size() const { return nb_pv_ + nb_pq_; }
    int vm_size()    const { return nb_pq_; }
    int total()      const { return lag_ + nb_pv_ + 2 * nb_pq_; }

    // --- state-vector segment accessors -------------------------------------
    // Both segments start at offset lag.  "vm" not "V" to avoid confusion with
    // the complex voltage vector.
    RealVect::SegmentReturnType theta(RealVect& x) const {
        return x.segment(lag_, theta_size());
    }
    RealVect::SegmentReturnType vm(RealVect& x) const {
        return x.segment(lag_ + theta_size(), vm_size());
    }
    RealVect::ConstSegmentReturnType theta(const RealVect& x) const {
        return x.segment(lag_, theta_size());
    }
    RealVect::ConstSegmentReturnType vm(const RealVect& x) const {
        return x.segment(lag_ + theta_size(), vm_size());
    }

    // --- Jacobian block offsets (rows and cols both shifted by lag) ----------
    int J_dP_row()        const { return lag_; }
    int J_dQ_row()        const { return lag_ + theta_size(); }
    int J_dP_dTheta_col() const { return lag_; }
    int J_dP_dV_col()     const { return lag_ + theta_size(); }
    int J_dQ_dTheta_col() const { return lag_; }
    int J_dQ_dV_col()     const { return lag_ + theta_size(); }

    // --- distributed-slack row/col (only meaningful when lag == 1) ----------
    int J_slack_row() const { return 0; }
    int J_slack_col() const { return 0; }

private:
    int nb_pv_;
    int nb_pq_;
    int lag_;
};

} // namespace ls2g

#endif // NR_LAYOUT_H
