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
#include <utility>
#include <vector>

#include "TaggedIdVec.hpp"
#include "Utils.hpp"
#include "VoltageControlData.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

class GeneratorContainer;
class GenericContainer;
class HvdcLineContainer;
class SvcContainer;

/**
 * How every bus of one solved system has its voltage magnitude pinned, and by what.
 *
 * The classical part -- the pv/pq split -- and the "fancy" part -- remote-regulating
 * generators, several machines regulating the same bus, voltage-mode SVCs -- derived
 * once, in one place, out of the element containers and one bus labelling.
 *
 * ---------------------------------------------------------------------------
 * WHY THIS IS ONE OBJECT
 *
 * Four answers used to be derived in four different places, by four functions on
 * LSGrid, each called from somewhere else in a powerflow:
 *
 *   get_group_controlled_buses()      -- fillpv_pq, while it is still building
 *                                        the very split the others are keyed on;
 *   fillpv_pq()                       -- _build_into_cache;
 *   get_free_vm_slack_solver_buses()  -- Base::update_state, and again inside
 *                                        fill_voltage_control_solver_data;
 *   fill_voltage_control_solver_data()-- VoltageControl::update_state.
 *
 * They are not four answers, they are four LAYERS of one, each built out of the
 * one before it. Deriving each where it happened to be needed meant walking every
 * generator of the grid three to four times per powerflow -- and, worse, it meant
 * as many separate opportunities for the layers to disagree, because each walk
 * read the containers again and nothing said they had to agree.
 *
 * Built as one object they cannot: each layer is built from the previous one's own
 * result, and the whole thing is built once and then read.
 *
 * ---------------------------------------------------------------------------
 * THE FOUR LAYERS, AND WHAT EACH NEEDS
 *
 *  1. `build_groups` -> `group_controlled_buses()` -- GRID bus ids: the buses a
 *     control GROUP regulates. Reads only input data (generators, SVCs), never the
 *     labelling, because layer 2 -- which produces part of that labelling -- needs
 *     it as an input.
 *  2. `build_pv_pq` -- the pv/pq split, SOLVER bus ids, written into the caller's
 *     vectors. Every container says which buses it pins (the classical part), then
 *     layer 1 takes back the ones a group regulates. Needs layer 1 and the slack.
 *  3. `build_free_vm_slack` -> `free_vm_slack_buses()` -- SOLVER bus ids of the
 *     slack buses that keep a free Vm unknown. Needs layer 1 and the labelling.
 *  4. `build_controllers` -> `controllers()` -- the full VoltageControlSolverData,
 *     SOLVER bus ids and per unit. Needs layers 1 to 3.
 *
 * Only layer 4 can fail: it throws, with the message the configuration deserves,
 * on what the v1 bordered formulation cannot express. The other three cannot.
 *
 * ---------------------------------------------------------------------------
 * ALGORITHMS THAT CANNOT DO ANY OF THIS
 *
 * Layers 3 and 4 are consumed by the `Base` and `VoltageControl` components of
 * NRSystem, i.e. by the Newton-Raphson algorithms and nothing else. Fast-decoupled
 * and Gauss-Seidel hold no NRSystem at all, so for them those two layers are work
 * whose result nobody reads -- and layer 2's group step is worse than useless: it
 * would take a bus out of PV and leave NOTHING to pin its magnitude, which is a
 * converged, plausible and wrong answer rather than an error.
 *
 * So the whole thing is conditional on `BaseAlgo::supports_remote_voltage_control()`,
 * passed to `build_groups` as `supports_voltage_control`. When it is false the plan
 * stays empty and layer 2 produces the classical split, and it is `LSGrid::ac_pf`
 * that refuses -- before any of this runs, and naming every offending element -- to
 * solve a grid that actually has controllers with an algorithm that cannot honour
 * them. See `list_unsupported`.
 *
 * The DC family is the one exception, and it is deliberate: its layer 2 keeps the
 * group step whatever the AC algorithm can do, because that is what it has always
 * done and BaseDCAlgo does read `pv`. See LSGrid::_build_into_cache.
 *
 * ---------------------------------------------------------------------------
 * LIFETIME
 *
 * A plan is a member of SolverBusLayout, i.e. of a solver-side cache: it is one
 * picture of the grid in ONE bus labelling, exactly like the labelling, the
 * slack and the pv-pq split it sits next to, and it is built, reused and retired
 * with them. Both families carry layers 1 and 2 (the pv/pq split is one rule, and
 * it does not change between them); only the AC family fills layers 3 and 4,
 * because a DC solve has no voltage control at all.
 *
 * A plan built against another labelling is not stale data, it is a different
 * grid: every bus id in it would silently mean another bus. That is why it does
 * not live in the NR extensions that consume it and is never carried across a
 * rebuild of the labelling -- see AlgoControl::need_recompute_voltage_control()
 * for what makes the cached one reusable, and LSGrid::_build_into_cache for
 * where the layers are built.
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
         *
         * `supports_voltage_control` false (a fast-decoupled or Gauss-Seidel solve)
         * leaves the set EMPTY, which is what makes layer 2 produce the classical
         * split. It is not a shortcut past a real configuration: `LSGrid::ac_pf`
         * refuses such a grid outright, see `list_unsupported`.
         */
        void build_groups(const GeneratorContainer & generators,
                          const SvcContainer & svcs,
                          bool supports_voltage_control = true);

        // ---- layer 2 -------------------------------------------------------------
        /**
         * The pv/pq split of one solved system, written into the caller's vectors
         * (they live in the solver-side cache, next to the labelling they are
         * expressed in -- this class owns the RULE, not the storage).
         *
         * Two steps. First every container says which buses it pins, through the
         * `fillpv` each of them overrides -- that is the classical part, and it is
         * why the containers arrive as a flat list: the rule is "ask all of them",
         * not "ask these eight in this order".
         *
         * Then layer 1 takes back the buses a control GROUP regulates. The per-container
         * pass only knows how to keep a controller's OWN bus out of PV (a generator
         * regulating remotely does not claim its own bus); nothing there stops a LOCAL
         * regulator sitting on somebody else's regulated bus from claiming it. When that
         * happened the group's bordered voltage row found no Vm column to write its +1
         * into -- a structurally empty row, hence a singular Jacobian -- so layer 4 had
         * to reject the whole configuration ("regulates bus X which has no voltage (Vm)
         * unknown") even though it is perfectly well posed: the local regulator simply
         * belongs in the group, and the sharing row then supplies the equation that
         * fixes the reactive split. Dropped from PV here, picked up as PQ, and enrolled
         * by `build_controllers`.
         *
         * `build_groups` must have run first -- the controller list built afterwards is
         * keyed on this very split, and both must come from one group layout.
         */
        void build_pv_pq(const std::vector<const GenericContainer *> & pv_sources,
                         const SolverBusIdVect & id_me_to_solver,
                         const GlobalBusIdVect & id_solver_to_me,
                         const SolverBusIdVect & slack_bus_id_solver,
                         SolverBusIdVect & bus_pv_out,
                         SolverBusIdVect & bus_pq_out) const;

        // ---- layer 3 -------------------------------------------------------------
        /**
         * The free-Vm slack buses, in the labelling passed in. `build_groups` must have
         * run against the same containers first -- this reads its result rather than
         * re-deriving it. Never throws: a slack bus either is locally pinned or is not.
         */
        void build_free_vm_slack(const GeneratorContainer & generators,
                                 const SolverBusIdVect & id_me_to_solver,
                                 const GlobalBusIdVect & id_solver_to_me,
                                 const SolverBusIdVect & slack_bus_id_solver);

        // ---- layer 4 -------------------------------------------------------------
        /**
         * The controller list. Every layer above must have been built against the same
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

        // ---- layers 3 and 4, which is what an AC solve wants ----------------------
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

        // ---- the guard an algorithm without the bordered block needs ---------------
        /**
         * Everything on this grid that would need the bordered block, so that an
         * algorithm which does not implement it can refuse the grid by name rather
         * than solve a system with a bus nothing pins.
         *
         * `build_groups` must have run first (with `supports_voltage_control` true,
         * or this answers "nothing" by construction). A converter station is never
         * an offender on its own -- it pins its bus the classical way unless a group
         * already claims it -- but it is listed when a group does, because the user
         * has to know which elements are affected, not only which ones caused it.
         */
        struct Unsupported {
            std::vector<int> gen_ids;   ///< generators regulating a bus that is not their own
            std::vector<int> svc_ids;   ///< voltage-mode SVCs (always group controllers)
            std::vector<std::pair<int, int> > station_ids;  ///< (hvdc line id, side) enrolled in a group
            [[nodiscard]] bool empty() const noexcept {
                return gen_ids.empty() && svc_ids.empty() && station_ids.empty();
            }
        };
        [[nodiscard]] Unsupported list_unsupported(const GeneratorContainer & generators,
                                                   const SvcContainer & svcs,
                                                   const HvdcLineContainer & hvdc_lines) const;

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
        std::set<int> free_vm_slack_buses_;  ///< layer 3, SOLVER bus ids
        VoltageControlSolverData controllers_;  ///< layer 4
        // layer 2 has no member: the split it produces belongs to the solver-side
        // cache (SolverBusLayout::bus_pv / bus_pq), where everything downstream
        // already reads it. This class owns the rule, not a second copy of the answer.
};

} // namespace ls2g

#endif // VOLTAGE_CONTROL_PLAN_H
