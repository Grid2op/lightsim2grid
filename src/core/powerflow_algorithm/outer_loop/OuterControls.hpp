// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_CONTROLS_H
#define OUTER_CONTROLS_H

#include <deque>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "Utils.hpp"

namespace ls2g {

class OuterControls;

/**
 * A bus whose units hold its voltage (PV) and may stop doing so (PQ). Its Vm unknown and its
 * Q row are reserved either way, so switching is a value edit (the Q row pinned while PV).
 * One per bus, shared by every loop that reserved it.
 */
class LS2G_API BusVoltageControl final
{
    public:
        int bus() const { return bus_; }
        bool is_pq() const { return pq_; }
        void set_pq() { pq_ = true; }
        void set_pv() { pq_ = false; }

    private:
        friend class OuterControls;
        explicit BusVoltageControl(int bus) : bus_(bus) {}
        void _reset() { pq_ = false; }

        int bus_;
        bool pq_ = false;
};

/**
 * Collects what the outer loops reserve in the solve and hands out the controls they act
 * through. A loop reserves, when the solver input is (re)built (BaseOuterLoop::declare),
 * every slot any of its states may need; the union is what lets a whole solve run on one
 * symbolic analysis. Afterwards it looks its controls up and edits them by value; the inner
 * algorithm applies them before each solve, in a fixed order. What nobody reserved costs
 * nothing.
 *
 * Algorithm-agnostic: the inner algorithm derives from it (NROuterInner's) and says what it
 * supports; a control it does not support is never handed out (nullptr).
 *
 * Solver bus ids throughout. The controls live until the next reservation (a rebuild of the
 * solver input), and are reset to "as the grid has it" at the start of every solve.
 */
class LS2G_API OuterControls
{
    public:
        virtual ~OuterControls() = default;

        // ----- reservation (BaseOuterLoop::declare) ---------------------------------------
        /// a bus that is PV in the labelling but may become PQ (or back); the same control for
        /// every loop asking
        BusVoltageControl * reserve_bus_voltage(int bus);

        // ----- lookup (after the reservation) -----------------------------------------------
        /// the control of `bus`, nullptr when no loop reserved it
        BusVoltageControl * bus_voltage(int bus);
        const BusVoltageControl * bus_voltage(int bus) const;

        // ----- edits that need no reservation --------------------------------------------
        /// restart the magnitude of `bus` from `vm_pu` at the next solve, once (a bus back to
        /// PV at its set-point, the robust mode's 1 pu); the last one asked for a bus wins
        void reset_vm(int bus, real_type vm_pu) { pending_vm_.emplace_back(bus, vm_pu); }
        /// a controller bus whose voltage control a loop took over for a while
        /// (TransformerVoltageControl): the other loops leave it alone meanwhile
        void set_suspended(int bus, bool val) { if (val) suspended_.insert(bus); else suspended_.erase(bus); }
        bool suspended(int bus) const { return suspended_.count(bus) > 0; }

        // ----- the driver's side -------------------------------------------------------------
        /// forget every reservation, before the loops declare again
        void clear_reservations();
        /// every control back to "as the grid has it", at the start of a solve
        void reset_states();

        /// the PV buses a caller (a batch) may switch, and those it keeps PV in the next solve,
        /// see BaseAlgo::set_switchable_vm_buses / set_pv_pinned_buses
        void set_caller_switchable(const std::vector<int> & buses) { caller_switchable_ = buses; }
        void set_caller_pinned(const std::vector<int> & buses) { caller_pinned_ = buses; }
        /// the caller's switchable buses, then the reserved ones (sorted)
        std::vector<int> switchable_buses() const;
        /// the caller's pinned buses, then the reserved ones that are PV (sorted)
        std::vector<int> pinned_buses() const;

        // ----- the inner algorithm's side ------------------------------------------------
        /// the magnitudes to reset (bus, pu), in the order they were asked for, then forgotten
        void take_pending_vm(std::vector<int> & buses, std::vector<real_type> & vm);

    private:
        // stable addresses: the loops keep pointers to these
        std::deque<BusVoltageControl> bus_voltage_store_;
        std::map<int, BusVoltageControl *> bus_voltage_of_;
        std::vector<int> caller_switchable_;
        std::vector<int> caller_pinned_;
        std::vector<std::pair<int, real_type> > pending_vm_;
        std::set<int> suspended_;
};

}  // namespace ls2g

#endif  // OUTER_CONTROLS_H
