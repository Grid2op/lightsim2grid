// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef HVDC_DROOP_DATA_H
#define HVDC_DROOP_DATA_H

#include "Utils.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

/**
 * Per-solve data of the CONNECTED angle-droop ("AC emulation") hvdc lines, in
 * SOLVER bus labelling and per-unit, as consumed by the solvers (the Hvdc
 * extension of the Newton-Raphson system, the DC algorithm).
 * Built by LSGrid::fill_hvdc_droop_solver_data.
 *
 * The flow equations mirror pyloadflow `equations/ac_emulation.py`:
 *   raw = p0 + k * (theta1 - theta2)
 * with the converter / dc-line losses applied on the received side, see
 * `recv_pu`. `status` is an INPUT (0 linear, +1 saturated 1->2,
 * -1 saturated 2->1), constant across one solve.
 */
struct LS2G_API HvdcDroopSolverData
{
    Eigen::VectorXi bus1;  // solver bus id of the side 1 station
    Eigen::VectorXi bus2;  // solver bus id of the side 2 station
    Eigen::VectorXi status;  // droop regime (INPUT)
    RealVect p0;      // pu
    RealVect k;       // pu / radian
    RealVect lf1;     // side 1 station loss factor (fraction)
    RealVect lf2;     // side 2 station loss factor (fraction)
    RealVect r;       // dc line resistance, pu (r_ohm * sn_mva / nominal_v^2)
    RealVect pmax12;  // pu
    RealVect pmax21;  // pu

    int size() const {return static_cast<int>(bus1.size());}

    void clear() {
        bus1 = Eigen::VectorXi();
        bus2 = Eigen::VectorXi();
        status = Eigen::VectorXi();
        p0 = RealVect();
        k = RealVect();
        lf1 = RealVect();
        lf2 = RealVect();
        r = RealVect();
        pmax12 = RealVect();
        pmax21 = RealVect();
    }

    /**
     * Active power received by the non-controller side for `p_ctrl_abs` (pu)
     * entering the controller side (pyloadflow `_ac_emul_abs_p_with_losses`).
     */
    real_type recv_pu(int i, real_type p_ctrl_abs, bool side_1_is_controller) const {
        const real_type lf_ctrl = side_1_is_controller ? lf1(i) : lf2(i);
        const real_type lf_recv = side_1_is_controller ? lf2(i) : lf1(i);
        const real_type line_in = (1. - lf_ctrl) * p_ctrl_abs;
        return (1. - lf_recv) * (line_in - r(i) * line_in * line_in);
    }

    /**
     * The two active flows LEAVING the AC buses into the hvdc (pu), given the
     * already computed `raw = p0 + k * (theta1 - theta2)` and the regime
     * (pyloadflow `ac_emulation_flows`).
     */
    void flows_pu(int i, real_type raw, real_type & p1_flow, real_type & p2_flow) const {
        const int st = status(i);
        if(st == 0){
            if(raw >= 0.){
                p1_flow = raw;
                p2_flow = -recv_pu(i, raw, true);
            }else{
                p1_flow = -recv_pu(i, -raw, false);
                p2_flow = -raw;
            }
        }else if(st > 0){
            p1_flow = pmax12(i);
            p2_flow = -recv_pu(i, pmax12(i), true);
        }else{
            p1_flow = -recv_pu(i, pmax21(i), false);
            p2_flow = pmax21(i);
        }
    }
};

class LSGrid;

/**
 * Free-function wrapper around LSGrid::fill_hvdc_droop_solver_data, defined in
 * LSGrid.cpp: lets template solver code (eg BaseDCAlgo.tpp) pull the droop
 * data through a forward-declared LSGrid pointer. A null grid yields empty data.
 */
LS2G_API void fill_hvdc_droop_data_from_grid(const LSGrid * lsgrid_ptr, HvdcDroopSolverData & data, bool ac);

} // namespace ls2g

#endif // HVDC_DROOP_DATA_H
