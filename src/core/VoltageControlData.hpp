// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef VOLTAGE_CONTROL_DATA_H
#define VOLTAGE_CONTROL_DATA_H

#include "Utils.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

/**
 * Per-solve data of the ACTIVE voltage-mode controllers (remote-regulating
 * generators, voltage-mode SVCs, and voltage-regulating hvdc converter stations
 * whose bus a group claims), in SOLVER bus labelling and per-unit, as
 * consumed by the `VoltageControl` extension of the Newton-Raphson system.
 * Built by LSGrid::fill_voltage_control_solver_data; AC only (DC ignores it).
 *
 * A "control group" g = { controllers c1..cN regulating the same solver bus
 * reg_bus(g) at setpoint v_set(g) }. The bordered formulation (see the comment
 * block of class VoltageControl in NRSystem.hpp) adds, per group: N reactive
 * injection columns Q_c (one per controller), 1 voltage-constraint row and N-1
 * reactive-sharing rows. Every regulated bus and every controller bus REMAINS an
 * ordinary PQ bus (it keeps its theta/vm unknowns and P/Q equations in Base).
 *
 * Controllers are laid out grouped contiguously: group g owns the controllers
 * [grp_start(g), grp_start(g) + grp_count(g)). The FIRST controller of each
 * group (index grp_start(g)) is the reference of the N-1 sharing rows.
 *
 * Conventions (pinned against pypowsybl OpenLoadFlow, Phase 0 probes):
 *   - q (filled by the extension, not here) is the reactive injection in
 *     GENERATOR convention (= -q_pypowsybl_terminal);
 *   - the voltage constraint is  Vm(reg) + s.Q_c - v_set = 0  with s = slope(c)
 *     (0 for generators and non-sloped SVCs);
 *   - the sharing keys w_i give pure proportional sharing Q_i / w_i = const, with
 *     w_i = qmax_i - qmin_i (OLF Q_EQUAL_PROPORTION).
 */
struct LS2G_API VoltageControlSolverData
{
    // Controller kind tag (stored per controller, used for result write-back).
    // HVDC_SIDE_1 / HVDC_SIDE_2 are the two converter stations of an hvdc line;
    // `elem_id` carries the hvdc LINE id for both, the kind saying which end. A
    // station is a controller only when a group regulates its bus -- on its own it
    // pins that bus through the ordinary PV path, like a local generator.
    enum Kind { GEN = 0, SVC = 1, HVDC_SIDE_1 = 2, HVDC_SIDE_2 = 3 };

    // ---- per controller (flat, grouped contiguously) ------------------------
    Eigen::VectorXi bus;       // controller solver bus (must own a Q equation)
    Eigen::VectorXi kind;      // VoltageControlSolverData::Kind
    Eigen::VectorXi elem_id;   // id of the controller in its own container
    RealVect        slope;     // pu (0 except for a sloped SVC)
    RealVect        weight;    // sharing key w_i (= qmax_i - qmin_i, pu cancels)
    Eigen::VectorXi group;     // group index this controller belongs to

    // ---- per group ----------------------------------------------------------
    Eigen::VectorXi reg_bus;   // regulated solver bus (must own a Vm unknown)
    RealVect        v_set;     // voltage setpoint, pu
    Eigen::VectorXi grp_start; // index of the first controller of the group
    Eigen::VectorXi grp_count; // number of controllers of the group

    int n_controllers() const { return static_cast<int>(bus.size()); }
    int n_groups()      const { return static_cast<int>(reg_bus.size()); }

    void clear() {
        bus = Eigen::VectorXi();
        kind = Eigen::VectorXi();
        elem_id = Eigen::VectorXi();
        slope = RealVect();
        weight = RealVect();
        group = Eigen::VectorXi();
        reg_bus = Eigen::VectorXi();
        v_set = RealVect();
        grp_start = Eigen::VectorXi();
        grp_count = Eigen::VectorXi();
    }
};

} // namespace ls2g

#endif // VOLTAGE_CONTROL_DATA_H
