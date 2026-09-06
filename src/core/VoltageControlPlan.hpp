// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef VOLTAGE_CONTROL_PLAN_H
#define VOLTAGE_CONTROL_PLAN_H

#include <set>

#include "TaggedIdVec.hpp"
#include "Utils.hpp"
#include "VoltageControlData.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

class GeneratorContainer;
class SvcContainer;
class HvdcLineContainer;

/**
 * Everything one AC solve needs to know about the "fancy" voltage controllers --
 * remote-regulating generators, several machines regulating the same bus, and
 * voltage-mode SVCs -- derived once, in one place, out of the element containers
 * and one bus labelling.
 *
 * ---------------------------------------------------------------------------
 * WHY THIS IS ONE OBJECT
 *
 * Three answers used to be re-derived independently, by three functions on
 * LSGrid, each of them called from a different place in a powerflow:
 *
 *   get_group_controlled_buses()      -- fillpv_pq, while it is still building
 *                                        the very labelling the other two need;
 *   get_free_vm_slack_solver_buses()  -- Base::update_state, and again inside
 *                                        fill_voltage_control_solver_data;
 *   fill_voltage_control_solver_data()-- VoltageControl::update_state.
 *
 * They are not three answers, they are three LAYERS of one: the controller list
 * is built out of the group layout and the free-Vm slack set, and the free-Vm
 * slack set is built out of the group layout. Re-deriving each of them where it
 * happened to be needed meant walking every generator of the grid three to four
 * times per powerflow -- and, worse, it meant three separate opportunities for
 * the layers to disagree, because each walk read the containers again and
 * nothing said they had to agree.
 *
 * Built as one object they cannot: layer 2 is built from layer 1's own result,
 * layer 3 from layer 2's, and the whole thing is built once and then read.
 *
 * ---------------------------------------------------------------------------
 * THE THREE LAYERS, AND WHAT EACH NEEDS
 *
 *  1. `group_controlled_buses()` -- GRID bus ids. Reads only input data
 *     (generators, SVCs), never the labelling, which is what lets `fillpv_pq`
 *     ask for it while it is still deciding the pv-pq split.
 *  2. `free_vm_slack_buses()`    -- SOLVER bus ids. Needs layer 1, the labelling
 *     and the slack.
 *  3. `controllers()`            -- the full VoltageControlSolverData, SOLVER bus
 *     ids and per unit. Needs layers 1 and 2, the labelling and the pv-pq split.
 *
 * Hence one build entry point per layer: `build_groups` is layer 1 alone (all
 * `fillpv_pq` can have and all it needs), `build_free_vm_slack` and
 * `build_controllers` are the two that need the labelling, with
 * `build_solver_side` as the pair of them in the only order that works.
 * `build_controllers` throws, with the messages the configuration deserves, on
 * what the v1 bordered formulation cannot express; the other two cannot fail.
 *
 * ---------------------------------------------------------------------------
 * LIFETIME
 *
 * A plan is a member of SolverBusLayout, i.e. of a solver-side cache: it is one
 * picture of the grid in ONE bus labelling, exactly like the labelling, the
 * slack and the pv-pq split it sits next to, and it is built, reused and retired
 * with them. Both families carry layer 1 (`fillpv_pq` is one function, and its
 * rule does not change between them); only the AC family fills layers 2 and 3,
 * because a DC solve has no voltage control at all.
 *
 * A plan built against another labelling is not stale data, it is a different
 * grid: every bus id in it would silently mean another bus. That is why it does
 * not live in the NR extensions that consume it and is never carried across a
 * rebuild of the labelling -- see AlgoControl::need_recompute_voltage_control()
 * for what makes the cached one reusable, and LSGrid::_build_into_cache for
 * where the layers are built, around the pv-pq split they straddle.
 */
class LS2G_API VoltageControlPlan
{
    public:
        VoltageControlPlan() = default;

        /// back to "nothing known": every layer empty
        void clear() noexcept;

        // ---- layer 1 -----------------------------------------------------------
        /**
         * GRID-bus ids whose voltage magnitude is set by a control GROUP rather than
         * by the classical PV treatment. A bus lands here as soon as an ACTIVE
         * remote-regulating generator or an active voltage-mode SVC aims at it,
         * because the group's bordered voltage row needs that bus to keep a Vm
         * unknown (and hence a Q equation).
         *
         * Buses regulated ONLY by things standing on them -- local generators, or
         * voltage-regulating hvdc converter stations -- are deliberately NOT
         * included: several machines sharing a bus and all regulating it locally is
         * the ordinary PV case, handled as before by the per-bus reactive
         * redistribution. It is the arrival of a remote controller (or an SVC, which
         * is always a group controller) that switches the bus over to the bordered
         * formulation -- and then every local regulator on that bus, generator or
         * converter station, joins the group as a co-controller instead of pinning it.
         */
        void build_groups(const GeneratorContainer & generators,
                          const SvcContainer & svcs);

        // ---- layer 2 -------------------------------------------------------------
        /**
         * The free-Vm slack buses, in the labelling passed in. `build_groups` must have
         * run against the same containers first -- this reads its result rather than
         * re-deriving it. Never throws: a slack bus either is locally pinned or is not.
         */
        void build_free_vm_slack(const GeneratorContainer & generators,
                                 const SolverBusIdVect & id_me_to_solver,
                                 const GlobalBusIdVect & id_solver_to_me,
                                 const SolverBusIdVect & slack_bus_id_solver);

        // ---- layer 3 -------------------------------------------------------------
        /**
         * The controller list. Both layers above must have been built against the same
         * containers and the same labelling first. THROWS, with the messages the
         * configuration deserves, on what the v1 bordered formulation cannot express:
         * a controller whose own bus owns no Q equation, a regulated bus with no Vm
         * unknown, conflicting setpoints inside one group, an SVC sharing its group.
         */
        void build_controllers(const GeneratorContainer & generators,
                               const SvcContainer & svcs,
                               const HvdcLineContainer & hvdc_lines,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const SolverBusIdVect & bus_pq);

        // ---- layers 2 and 3, which is what a solve wants -------------------------
        /**
         * `build_free_vm_slack` then `build_controllers`, in the only order that works.
         * AC only, and the caller says so by not calling this at all for a DC cache.
         * `build_groups` must have run first, against the same containers.
         */
        void build_solver_side(const GeneratorContainer & generators,
                               const SvcContainer & svcs,
                               const HvdcLineContainer & hvdc_lines,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const SolverBusIdVect & slack_bus_id_solver,
                               const SolverBusIdVect & bus_pq);

        // ---- what a solve reads --------------------------------------------------
        [[nodiscard]] const std::set<int> & group_controlled_buses() const noexcept {
            return group_reg_buses_;
        }
        /**
         * Solver-bus ids of the slack buses that need a free Vm unknown and a Q
         * equation (added by the Base block of the NR system), i.e. every slack bus
         * whose magnitude is NOT pinned by a local voltage-regulating generator.
         * This covers distributed-slack participants whose generator is PQ
         * (voltage_regulator_on == false), slack buses hosting a remote-voltage
         * controller, and SVC-regulated slack buses. A slack bus that DOES host a
         * local PV generator stays Vm-fixed (PV-like) with no Q equation.
         */
        [[nodiscard]] const std::set<int> & free_vm_slack_buses() const noexcept {
            return free_vm_slack_buses_;
        }
        /// the bordered block's own data, as the VoltageControl extension consumes it
        [[nodiscard]] const VoltageControlSolverData & controllers() const noexcept {
            return controllers_;
        }

    private:
        /// one candidate controller, before the grouping pass
        struct Raw {
            int bus;          ///< controller solver bus
            int reg_bus;      ///< regulated solver bus
            real_type v_set;  ///< pu
            real_type slope;  ///< pu (0 except for a sloped SVC)
            real_type weight; ///< sharing key
            int kind;         ///< VoltageControlSolverData::Kind
            int elem_id;
        };

        /// the three per-container passes that fill `raws`
        void _collect_gen_controllers(const GeneratorContainer & generators,
                                      const SolverBusIdVect & id_me_to_solver,
                                      const std::vector<bool> & is_pq,
                                      const std::vector<bool> & has_free_q,
                                      std::vector<Raw> & raws) const;
        void _collect_svc_controllers(const SvcContainer & svcs,
                                      const SolverBusIdVect & id_me_to_solver,
                                      const std::vector<bool> & is_pq,
                                      const std::vector<bool> & has_free_q,
                                      std::vector<Raw> & raws) const;
        void _collect_station_controllers(const HvdcLineContainer & hvdc_lines,
                                          const SolverBusIdVect & id_me_to_solver,
                                          const std::vector<bool> & is_pq,
                                          const std::vector<bool> & has_free_q,
                                          std::vector<Raw> & raws) const;
        /// group by regulated bus, check the setpoints agree, emit `controllers_`
        void _group_and_emit(const std::vector<Raw> & raws);

        std::set<int> group_reg_buses_;      ///< layer 1, GRID bus ids
        std::set<int> free_vm_slack_buses_;  ///< layer 2, SOLVER bus ids
        VoltageControlSolverData controllers_;  ///< layer 3
};

} // namespace ls2g

#endif // VOLTAGE_CONTROL_PLAN_H
