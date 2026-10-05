// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_CONTROLS_H
#define OUTER_CONTROLS_H

#include <cmath>
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
 * A voltage controller of a group (remote regulation, an SVC, a station) that may be held at
 * a reactive output, by value: its share of the group's rows pins it there, the group's other
 * controllers regulating on. Reserving one writes every group's rows in the form where any
 * controller may be held (VoltageControl::set_may_hold_controllers).
 */
class LS2G_API VoltageControllerHold final
{
    public:
        /// its position in the voltage-control plan's controller list
        int controller() const { return controller_; }
        bool is_held() const { return std::isfinite(q_); }
        /// the reactive output it is held at, pu (NaN when it regulates)
        real_type q() const { return q_; }
        void hold(real_type q_pu) { q_ = q_pu; }
        void release() { q_ = std::numeric_limits<real_type>::quiet_NaN(); }

    private:
        friend class OuterControls;
        explicit VoltageControllerHold(int controller) : controller_(controller) {}

        int controller_;
        real_type q_ = std::numeric_limits<real_type>::quiet_NaN();
};

/**
 * An idle standby SVC (held at Q = 0 by the plan, alone in its group) that may be switched on
 * to regulate its bus: the group's "Q = 0" row turns back into its voltage row, entries
 * declared either way (VoltageControl::release_held_svcs).
 */
class LS2G_API StandbySvcControl final
{
    public:
        int svc() const { return svc_; }
        bool is_on() const { return std::isfinite(target_vm_); }
        /// the set-point it regulates at once switched on, pu (NaN while idle)
        real_type target_vm() const { return target_vm_; }
        void switch_on(real_type target_vm_pu) { target_vm_ = target_vm_pu; }

    private:
        friend class OuterControls;
        explicit StandbySvcControl(int svc) : svc_(svc) {}
        void _reset() { target_vm_ = std::numeric_limits<real_type>::quiet_NaN(); }

        int svc_;
        real_type target_vm_ = std::numeric_limits<real_type>::quiet_NaN();
};

/**
 * The regime of an hvdc line in AC emulation: 0 linear, +1 saturated 1 -> 2, -1 saturated
 * 2 -> 1, KEEP as the grid has it. Every regime's entries are declared, so a change is a
 * value edit (Hvdc::set_status_override).
 */
class LS2G_API HvdcRegimeControl final
{
    public:
        static constexpr int KEEP = 2;
        int line() const { return line_; }
        int regime() const { return regime_; }
        void set(int regime) { regime_ = regime; }

    private:
        friend class OuterControls;
        explicit HvdcRegimeControl(int line) : line_(line) {}
        void _reset() { regime_ = KEEP; }

        int line_;
        int regime_ = KEEP;
};

/**
 * A phase shifter (a transformer, grid id): its tap may move and, when the solve solves for
 * its shift (`solves_shift`, a column and a row), its active power control may be switched
 * on or off. Its block of Ybus is patched by value (BranchControl). The reads are the last
 * solve's, and mean nothing unless the solve handles it (handled()).
 */
class LS2G_API PhaseShifterControl final
{
    public:
        int trafo() const { return trafo_; }
        bool solves_shift() const { return solves_shift_; }
        /// whether the solve handles it (a transformer connected at both ends, on two buses)
        bool handled() const;
        /// its shift now, rad
        real_type shift() const;
        /// its tap position now
        int position() const;
        /// the current through its `side` (1 or 2), pu of that side's base, and its derivative
        /// with respect to the shift
        void current(int side, real_type & i_pu, real_type & di_da) const;

        /// switch its active power control on or off (off: its shift stays where it is)
        void set_control(bool on) { control_ = on ? 1 : 0; }
        /// move its tap to `position`: its shift is that tap's from the next solve on
        void move_tap(int position) { tap_ = position; tap_moved_ = true; }

        /// what a loop asked, for the inner algorithm: the control (1 on, 0 off, -1 as it is)
        /// and the tap (when moved)
        int requested_control() const { return control_; }
        bool tap_moved() const { return tap_moved_; }
        int requested_tap() const { return tap_; }

    private:
        friend class OuterControls;
        PhaseShifterControl(const OuterControls * owner, int trafo, bool solves_shift)
            : owner_(owner), trafo_(trafo), solves_shift_(solves_shift) {}
        void _reset() { control_ = -1; tap_moved_ = false; tap_ = 0; }

        const OuterControls * owner_;
        int trafo_;
        bool solves_shift_;
        int control_ = -1;
        bool tap_moved_ = false;
        int tap_ = 0;
};

/**
 * A transformer of a ratio group (grid id): transformers regulating the voltage of one bus.
 * Its ratio tap may move and, when the group is solved (`solved`, a column per transformer and
 * the group's rows), its voltage control may be switched on or off. The reads are the last
 * solve's, and mean nothing unless the solve handles it (handled()).
 */
class LS2G_API RatioTapControl final
{
    public:
        int trafo() const { return trafo_; }
        /// whether the solve handles it (a transformer connected at both ends, on two buses)
        bool handled() const;
        /// its ratio now
        real_type ratio() const;
        /// its ratio tap position now
        int position() const;

        /// switch its voltage control on or off (off: its ratio stays where it is)
        void set_control(bool on) { control_ = on ? 1 : 0; }
        /// move its ratio tap to `position`, once, at the next solve
        void move_tap(int position) { tap_ = position; tap_moved_ = true; }

        /// what a loop asked, for the inner algorithm: the control (1 on, 0 off, -1 as it is)
        int requested_control() const { return control_; }
        /// the tap a loop moved it to, once: false when it did not
        bool take_tap(int & position) {
            if (!tap_moved_) return false;
            position = tap_;
            tap_moved_ = false;
            return true;
        }

    private:
        friend class OuterControls;
        RatioTapControl(const OuterControls * owner, int trafo) : owner_(owner), trafo_(trafo) {}
        void _reset() { control_ = -1; tap_moved_ = false; tap_ = 0; }

        const OuterControls * owner_;
        int trafo_;
        int control_ = -1;
        bool tap_moved_ = false;
        int tap_ = 0;
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
        OuterControls() = default;
        virtual ~OuterControls() = default;
        // the controls point back to it
        OuterControls(const OuterControls &) = delete;
        OuterControls & operator=(const OuterControls &) = delete;

        // ----- reservation (BaseOuterLoop::declare) ---------------------------------------
        /// a bus that is PV in the labelling but may become PQ (or back); the same control for
        /// every loop asking
        BusVoltageControl * reserve_bus_voltage(int bus);
        /// a voltage controller (its position in the plan's list) that may be held
        VoltageControllerHold * reserve_controller_hold(int controller);
        /// an idle standby SVC (grid id) that may be switched on
        StandbySvcControl * reserve_standby_svc(int svc);
        /// an hvdc line (grid id) in AC emulation whose regime may change
        HvdcRegimeControl * reserve_hvdc_regime(int line);
        /// a phase shifter (transformer grid id) whose tap may move and, with `solves_shift`,
        /// whose shift the solve solves for; nullptr when the inner algorithm has none
        PhaseShifterControl * reserve_phase_shifter(int trafo, bool solves_shift);
        /// transformers (grid ids, in order: the first one on holds the voltage) regulating
        /// the voltage of `bus` at `target_vm` pu, their ratios solved for when `solved`; a
        /// RatioTapControl per transformer. False when the inner algorithm has none.
        bool reserve_ratio_group(int bus, real_type target_vm, const std::vector<int> & trafos, bool solved);

        // ----- lookup (after the reservation) -----------------------------------------------
        /// the control of `bus`, nullptr when no loop reserved it
        BusVoltageControl * bus_voltage(int bus);
        const BusVoltageControl * bus_voltage(int bus) const;
        VoltageControllerHold * controller_hold(int controller);
        const VoltageControllerHold * controller_hold(int controller) const;
        StandbySvcControl * standby_svc(int svc);
        const StandbySvcControl * standby_svc(int svc) const;
        HvdcRegimeControl * hvdc_regime(int line);
        const HvdcRegimeControl * hvdc_regime(int line) const;
        PhaseShifterControl * phase_shifter(int trafo);
        const PhaseShifterControl * phase_shifter(int trafo) const;
        RatioTapControl * ratio_tap(int trafo);
        const RatioTapControl * ratio_tap(int trafo) const;
        /// the buses a reserved ratio group regulates (a control of a higher priority than a
        /// shunt's, see ShuntVoltageControlLoop)
        std::set<int> ratio_group_buses() const;
        /// whether a loop holds that controller now
        bool is_held(int controller) const {
            const VoltageControllerHold * hold = controller_hold(controller);
            return hold != nullptr && hold->is_held();
        }

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
        /// whether any controller may be held
        bool holds_voltage_controllers() const { return !controller_hold_of_.empty(); }
        /// the output each controller is held at, by position in the plan (NaN: not held)
        std::vector<real_type> held_q() const;
        /// the set-point each standby SVC was switched on at, by grid id (NaN: idle); empty
        /// when none is reserved
        std::vector<real_type> svc_target_vm() const;
        /// the regime of each hvdc line, by grid id (HvdcRegimeControl::KEEP: the grid's);
        /// empty when none is reserved
        std::vector<int> hvdc_regimes() const;
        /// the phase shifters, in the order they were reserved / by transformer id
        const std::vector<PhaseShifterControl *> & phase_shifters_reserved() const { return phase_shifter_order_; }
        const std::map<int, PhaseShifterControl *> & phase_shifters() const { return phase_shifter_of_; }
        /// the ratio groups, in the order they were reserved, and their transformers by id
        struct RatioGroup {
            int bus;
            real_type target_vm;
            std::vector<int> trafos;
            bool solved;
        };
        const std::vector<RatioGroup> & ratio_groups() const { return ratio_groups_; }
        const std::map<int, RatioTapControl *> & ratio_taps() const { return ratio_tap_of_; }

        /// the magnitudes to reset (bus, pu), in the order they were asked for, then forgotten
        void take_pending_vm(std::vector<int> & buses, std::vector<real_type> & vm);

    protected:
        // what the solve reads back, for the controls above (the inner algorithm's)
        friend class PhaseShifterControl;
        friend class RatioTapControl;
        virtual bool _supports_phase_shifters() const { return false; }
        virtual bool _supports_ratio_groups() const { return false; }
        virtual bool _ratio_handled(int /*trafo*/) const { return false; }
        virtual real_type _ratio(int /*trafo*/) const { return std::numeric_limits<real_type>::quiet_NaN(); }
        virtual int _ratio_position(int /*trafo*/) const { return 0; }
        virtual bool _phase_handled(int /*trafo*/) const { return false; }
        virtual real_type _phase_shift(int /*trafo*/) const { return std::numeric_limits<real_type>::quiet_NaN(); }
        virtual int _phase_position(int /*trafo*/) const { return 0; }
        virtual void _phase_current(int /*trafo*/, int /*side*/, real_type & i_pu, real_type & di_da) const {
            i_pu = std::numeric_limits<real_type>::quiet_NaN();
            di_da = std::numeric_limits<real_type>::quiet_NaN();
        }

    private:
        // stable addresses: the loops keep pointers to these
        std::deque<BusVoltageControl> bus_voltage_store_;
        std::map<int, BusVoltageControl *> bus_voltage_of_;
        std::deque<VoltageControllerHold> controller_hold_store_;
        std::map<int, VoltageControllerHold *> controller_hold_of_;
        std::deque<StandbySvcControl> standby_svc_store_;
        std::map<int, StandbySvcControl *> standby_svc_of_;
        std::deque<HvdcRegimeControl> hvdc_regime_store_;
        std::map<int, HvdcRegimeControl *> hvdc_regime_of_;
        std::deque<PhaseShifterControl> phase_shifter_store_;
        std::map<int, PhaseShifterControl *> phase_shifter_of_;
        std::vector<PhaseShifterControl *> phase_shifter_order_;
        std::vector<RatioGroup> ratio_groups_;
        std::deque<RatioTapControl> ratio_tap_store_;
        std::map<int, RatioTapControl *> ratio_tap_of_;
        std::vector<int> caller_switchable_;
        std::vector<int> caller_pinned_;
        std::vector<std::pair<int, real_type> > pending_vm_;
        std::set<int> suspended_;
};

}  // namespace ls2g

#endif  // OUTER_CONTROLS_H
