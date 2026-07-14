// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef HVDCLINECONTAINER_H
#define HVDCLINECONTAINER_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "Utils.hpp"
#include "SubstationContainer.hpp"
#include "TwoSidesContainer.hpp"
#include "ConverterStationContainer.hpp"

namespace ls2g {

class HvdcLineContainer;
class LS2G_API HvdcLineInfo final : public TwoSidesContainer<ConverterStationContainer>::TwoSidesInfo
{
    public:
        // members (the legacy DCLine names are kept: target_p_1_mw is the
        // side-1 station active power, generator sign convention)
        real_type target_p_1_mw;
        real_type p_2_mw;
        real_type target_vm_1_pu;
        real_type target_vm_2_pu;
        real_type loss_pct;
        real_type loss_mw;

        // IIDM model
        int converters_mode;       // 0 = side 1 rectifier, 1 = side 2 rectifier
        real_type p_setpoint_mw;   // active power drawn at the rectifier (>= 0)
        real_type r_ohm;
        real_type nominal_v_kv;

        // angle-droop (AC emulation)
        bool droop_enabled;
        real_type droop_p0_mw;
        real_type droop_k_mw_per_rad;
        real_type pmax_1to2_mw;
        real_type pmax_2to1_mw;
        int status_droop;          // 0 linear, +1 saturated 1->2, -1 saturated 2->1

        ConverterStationInfo station_side_1;
        ConverterStationInfo station_side_2;

        inline HvdcLineInfo(const HvdcLineContainer & r_data_hvdc, int my_id) noexcept;
};


/**
This class represents the list of HVDC lines (it replaces the legacy
pandapower-shaped "DC lines", which are now a special case of it).

The model follows powsybl IIDM / open-loadflow: the line owns the active power
(`p_setpoint_mw_` drawn at the rectifier, `converters_mode_` saying which side
rectifies), the two embedded converter stations (side_1_ / side_2_, see
`ConverterStationContainer`) own the reactive / voltage behaviour.

Loss model (mirrors open-loadflow `HvdcUtils.getConverterStationTargetP`,
extended with the pandapower fixed loss for the exact `loss_mw` round trip):

    line_in  = (1 - lf_rect) * (1 - loss_percent/100) * p_setpoint
    line_loss = r_ohm * line_in^2 / nominal_v^2        (0 when nominal_v == 0)
    received = (1 - lf_inv) * (line_in - line_loss) - loss_mw

where lf_rect / lf_inv are the rectifier / inverter station loss factors
(fractions). The legacy pandapower dcline maps onto this exactly with
station loss factors = 0, r_ohm = 0 (see `init_legacy`).

Angle-droop (IIDM `HvdcAngleDroopActivePowerControl`, OLF "AC emulation"):
when `droop_enabled_`, the AC active power follows
`raw = p0 + k * (theta1 - theta2)` instead of the fixed setpoint. The flows /
derivatives mirror pyloadflow `equations/ac_emulation.py`. `status_droop_` is
an INPUT (0 linear, +1 saturated 1->2, -1 saturated 2->1), constant across one
solve: saturation is decided by an outer loop (in Python) between solves.
For droop lines, this container does NOT stamp the AC active power in Sbus:
the NR system Hvdc extension owns it (all regimes). In DC, linear-mode lines
are handled by the DC algorithm (theta-dependent term) and saturated lines are
stamped here as fixed injections.
**/
class LS2G_API HvdcLineContainer final : public TwoSidesContainer<ConverterStationContainer>, public IteratorAdder<HvdcLineContainer, HvdcLineInfo>
{
    friend class HvdcLineInfo;

    public:
        using DataInfo = HvdcLineInfo;

        // /!\ these integer values are serialized verbatim (binary files and
        // pickles, as converters_mode_): renumbering requires bumping
        // BINARY_FORMAT_VERSION (BinaryArchive.hpp). Guarded python side by
        // TestSerializedEnumValues in test_binary_serialization.py.
        enum ConvertersMode {
            SIDE_1_RECTIFIER = 0,
            SIDE_2_RECTIFIER = 1
        };

        // Unlike a normal branch, the two converter stations of an HVDC line are NOT
        // synchronous partners: disconnecting one (its DC partner switched off) does
        // not electrically justify tripping the other, which OpenLoadFlow keeps
        // active as a local reactive/voltage-support device. So, unlike
        // TwoSidesContainer's generic default (mirror both sides), HVDC lines default
        // to independent sides.
        HvdcLineContainer() noexcept { synch_status_both_side_ = false; }
        virtual ~HvdcLineContainer() noexcept = default;

        // pickle
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
                TwoSidesContainer<ConverterStationContainer>::StateRes,
                std::vector<real_type>, // loss_percent_
                std::vector<real_type>, // loss_mw_
                std::vector<int>,       // converters_mode_
                std::vector<real_type>, // p_setpoint_mw_
                std::vector<real_type>, // r_ohm_
                std::vector<real_type>, // nominal_v_kv_
                std::vector<bool>,      // droop_enabled_
                std::vector<real_type>, // p0_mw_
                std::vector<real_type>, // k_mw_per_rad_
                std::vector<real_type>, // pmax_1to2_mw_
                std::vector<real_type>, // pmax_2to1_mw_
                std::vector<int>        // status_droop_
                >;
        enum StateResIdx {
            TSC_STATE = 0,
            LOSS_PERCENT,
            LOSS_MW,
            CONVERTERS_MODE,
            P_SETPOINT_MW,
            R_OHM,
            NOMINAL_V_KV,
            DROOP_ENABLED,
            P0_MW,
            K_MW_PER_RAD,
            PMAX_1TO2_MW,
            PMAX_2TO1_MW,
            STATUS_DROOP,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "HvdcLineContainer::StateRes and StateResIdx do not match");
        HvdcLineContainer::StateRes get_state() const;
        void set_state(HvdcLineContainer::StateRes & my_state);

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static HvdcLineContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "HvdcLineContainer"; }  // written into / checked against the binary file header

        /**
         * Full IIDM-style initialization (used by the pypowsybl converter and
         * the direct API). Loss factors are fractions (0 - 1), the droop is
         * given in MW per degree (converted to MW per radian internally).
         */
        void init(const Eigen::Ref<const Eigen::VectorXi> & bus1_id,
                  const Eigen::Ref<const Eigen::VectorXi> & bus2_id,
                  const std::vector<int> & type1,
                  const std::vector<int> & type2,
                  const Eigen::Ref<const RealVect> & loss_factor1,
                  const Eigen::Ref<const RealVect> & loss_factor2,
                  const std::vector<bool> & vreg1_on,
                  const std::vector<bool> & vreg2_on,
                  const Eigen::Ref<const RealVect> & vm1_pu,
                  const Eigen::Ref<const RealVect> & vm2_pu,
                  const Eigen::Ref<const RealVect> & q1_setpoint_mvar,
                  const Eigen::Ref<const RealVect> & q2_setpoint_mvar,
                  const Eigen::Ref<const RealVect> & min_q1,
                  const Eigen::Ref<const RealVect> & max_q1,
                  const Eigen::Ref<const RealVect> & min_q2,
                  const Eigen::Ref<const RealVect> & max_q2,
                  const Eigen::Ref<const RealVect> & power_factor1,
                  const Eigen::Ref<const RealVect> & power_factor2,
                  const std::vector<int> & converters_mode,
                  const Eigen::Ref<const RealVect> & p_setpoint_mw,
                  const Eigen::Ref<const RealVect> & r_ohm,
                  const Eigen::Ref<const RealVect> & nominal_v_kv,
                  const Eigen::Ref<const RealVect> & loss_percent,
                  const Eigen::Ref<const RealVect> & loss_mw,
                  const std::vector<bool> & droop_enabled,
                  const Eigen::Ref<const RealVect> & droop_p0_mw,
                  const Eigen::Ref<const RealVect> & droop_mw_per_deg,
                  const Eigen::Ref<const RealVect> & pmax_1to2_mw,
                  const Eigen::Ref<const RealVect> & pmax_2to1_mw
                  );

        /**
         * Legacy pandapower-shaped initialization (signature of the old
         * DCLineContainer::init): two regulating VSC stations, `p_mw` is the
         * side-1 station active power in generator convention (the python
         * layer negates pandapower's `p_mw` before calling this).
         */
        void init_legacy(const Eigen::Ref<const Eigen::VectorXi> & branch_from_id,
                         const Eigen::Ref<const Eigen::VectorXi> & branch_to_id,
                         const Eigen::Ref<const RealVect> & p_mw,
                         const Eigen::Ref<const RealVect> & loss_percent,
                         const Eigen::Ref<const RealVect> & loss_mw,
                         const Eigen::Ref<const RealVect> & vm_or_pu,
                         const Eigen::Ref<const RealVect> & vm_ex_pu,
                         const Eigen::Ref<const RealVect> & min_q_or,
                         const Eigen::Ref<const RealVect> & max_q_or,
                         const Eigen::Ref<const RealVect> & min_q_ex,
                         const Eigen::Ref<const RealVect> & max_q_ex
                         );

        // accessor / modifiers
        virtual void reconnect_connected_buses(SubstationContainer & substation) const override {
            side_1_.reconnect_connected_buses(substation);
            side_2_.reconnect_connected_buses(substation);
        }

        virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & /*res*/) const override {
            // for buses only connected through a hvdc line, i don't add them
            // they are not in the same "connected component"
        }

        // An HVDC line bridges two AC grids that are NOT synchronous (no edge is
        // added in `get_graph`), so the connectivity BFS of
        // `LSGrid::consider_only_main_component` is really splitting the grid into
        // *synchronous* components. When the two converters land in different
        // components, only one of them can be in the main (solved) synchronous
        // component. OpenLoadFlow keeps that in-main converter as a fixed boundary
        // injection (it still imports / exports its scheduled HVDC power); dropping
        // it would silently lose hundreds of MW. So, unlike a normal branch, we do
        // NOT deactivate the whole line when a single side is outside the main
        // component: we keep the in-main converter injecting and open only the
        // out-of-main one. A line with BOTH sides outside is still fully dropped.
        virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component) override {
            const int nb_el = nb();
            DualAlgoControl unused_solver_control;
            const GlobalBusIdVect & bus_side_1_id_ = get_buses_side_1();
            const GlobalBusIdVect & bus_side_2_id_ = get_buses_side_2();
            for(int i = 0; i < nb_el; ++i){
                if(!status_global_[i]){
                    side_1_.deactivate(i, unused_solver_control);
                    side_2_.deactivate(i, unused_solver_control);
                    continue;
                }
                const int b1 = bus_side_1_id_(i).cast_int();
                const int b2 = bus_side_2_id_(i).cast_int();
                const bool s1_outside = (b1 != _deactivated_bus_id) && !busbar_in_main_component[b1];
                const bool s2_outside = (b2 != _deactivated_bus_id) && !busbar_in_main_component[b2];
                if(s1_outside && s2_outside){
                    // both converters in (an)other synchronous component: drop it all
                    side_1_.deactivate(i, unused_solver_control);
                    side_2_.deactivate(i, unused_solver_control);
                    if(!ignore_status_global_) status_global_[i] = false;
                } else if(s1_outside || s2_outside){
                    // exactly one converter in the main synchronous component: keep
                    // the HVDC line active and that converter injecting its scheduled
                    // power, open only the converter that is out of the main component.
                    // Angle-droop ("AC emulation") cannot run across the cut (the
                    // remote angle is gone) -> fall back to the fixed power setpoint.
                    droop_enabled_[i] = false;
                    if(s1_outside) side_1_.deactivate(i, unused_solver_control);
                    else           side_2_.deactivate(i, unused_solver_control);
                    // status_global_[i] stays true: the line still injects in-main
                }
            }
        }

        real_type get_qmin_or(int hvdc_id) const {return side_1_.get_qmin(hvdc_id);}
        real_type get_qmax_or(int hvdc_id) const {return side_1_.get_qmax(hvdc_id);}
        real_type get_qmin_ex(int hvdc_id) const {return side_2_.get_qmin(hvdc_id);}
        real_type get_qmax_ex(int hvdc_id) const {return side_2_.get_qmax(hvdc_id);}

        /**
         * Change the active power of the line, legacy convention: `new_p` is
         * the side-1 station active power, generator convention (it is the
         * negated pandapower `p_mw`).
         */
        void change_p(int hvdc_id, real_type new_p, DualAlgoControl & solver_control);
        void change_v_side_1(int hvdc_id, real_type new_v_pu, DualAlgoControl & solver_control){
            side_1_.change_v(hvdc_id, new_v_pu, solver_control);
        }
        void change_v_side_2(int hvdc_id, real_type new_v_pu, DualAlgoControl & solver_control){
            side_2_.change_v(hvdc_id, new_v_pu, solver_control);
        }

        // droop API
        bool is_droop_active(int hvdc_id) const {
            return droop_enabled_[hvdc_id] && status_global_[hvdc_id];
        }
        // Angle-droop ("AC emulation") needs both remote angles: cannot run once either
        // converter is individually opened while the line stays connected_global (a
        // "half-open" HVDC, see deactivate_side_1/2 and
        // disconnect_if_not_in_main_component, which does the same for the
        // main-component case) -- fall back to the fixed power setpoint. Called by
        // LSGrid::deactivate_dcline_side1/2; without it, fill_hvdc_droop_solver_data
        // would try to resolve the now-disconnected (-1) side's bus and read out of
        // bounds.
        void disable_droop(int hvdc_id) { droop_enabled_[hvdc_id] = false; }
        bool has_droop_active() const {
            const int nb_hvdc = static_cast<int>(nb());
            for(int i = 0; i < nb_hvdc; ++i) if(is_droop_active(i)) return true;
            return false;
        }
        int get_status_droop(int hvdc_id) const {
            // hvdc_id indexes status_droop_ with an unchecked Eigen operator() (OOB read);
            // the matching setter set_status_droop already guards the same way.
            _check_in_range(hvdc_id, status_droop_, "get_status_droop");
            return status_droop_(hvdc_id);
        }
        std::vector<int> get_status_droop_vect() const {
            return std::vector<int>(status_droop_.begin(), status_droop_.end());
        }
        const std::vector<bool> & get_droop_enabled() const {return droop_enabled_;}
        real_type get_droop_p0_mw(int hvdc_id) const {return p0_mw_(hvdc_id);}
        real_type get_droop_k_mw_per_rad(int hvdc_id) const {return k_mw_per_rad_(hvdc_id);}
        real_type get_pmax_1to2_mw(int hvdc_id) const {return pmax_1to2_mw_(hvdc_id);}
        real_type get_pmax_2to1_mw(int hvdc_id) const {return pmax_2to1_mw_(hvdc_id);}
        real_type get_loss_factor_side_1(int hvdc_id) const {return side_1_.get_loss_factor(hvdc_id);}
        real_type get_loss_factor_side_2(int hvdc_id) const {return side_2_.get_loss_factor(hvdc_id);}
        real_type get_r_ohm(int hvdc_id) const {return r_ohm_(hvdc_id);}
        real_type get_nominal_v_kv(int hvdc_id) const {return nominal_v_kv_(hvdc_id);}
        /**
         * Set the droop regime (INPUT of the solver): 0 = linear,
         * +1 = saturated 1 -> 2, -1 = saturated 2 -> 1.
         * Only `tell_recompute_sbus` is signalled: the Jacobian sparsity
         * pattern is unchanged (the droop entries are declared whatever the
         * regime), so the symbolic factorization of the linear solver is reused.
         */
        void set_status_droop(int hvdc_id, int status, DualAlgoControl & solver_control);

        /**
         * Active power received by the non-controller side (MW), for
         * `p_ctrl_abs_mw` (>= 0) entering the controller side.
         * Mirror of pyloadflow `_ac_emul_abs_p_with_losses` (in SI units),
         * minus the pandapower fixed loss.
         */
        real_type recv_mw(int hvdc_id, real_type p_ctrl_abs_mw, bool side_1_is_controller) const;

        /**
         * The two active flows LEAVING the AC buses into the HVDC (MW), per
         * the droop equations (`raw_mw = p0 + k * (theta1 - theta2)`), regime
         * given by `status_droop_`. Mirror of pyloadflow `ac_emulation_flows`.
         */
        void droop_flows_mw(int hvdc_id, real_type raw_mw, real_type & p1_flow_mw, real_type & p2_flow_mw) const;

        // solver stuff
        virtual void fillSbus(CplxVect & Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;

        virtual void fillpv(std::vector<int>& bus_pv,
                            std::vector<bool> & has_bus_been_added,
                            const SolverBusIdVect & slack_bus_id_solver,
                            const SolverBusIdVect & id_grid_to_solver) const override {
            side_1_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_grid_to_solver);
            side_2_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_grid_to_solver);
        }

        virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & /*Bp*/,
                                std::vector<Eigen::Triplet<real_type> > & /*Bpp*/,
                                const SolverBusIdVect & /*id_grid_to_solver*/,
                                real_type /*sn_mva*/,
                                FDPFMethod /*xb_or_bx*/) const override {
                                    // no Bp coeffs for hvdc lines
                                }

        void init_q_vector(int nb_bus,
                           Eigen::VectorXi & total_gen_per_bus,
                           RealVect & total_q_min_per_bus,
                           RealVect & total_q_max_per_bus) const {
            side_1_.init_q_vector(nb_bus, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
            side_2_.init_q_vector(nb_bus, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
        }

        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const SolverBusIdVect & id_grid_to_solver,
                             const Eigen::Ref<const RealVect> & bus_vn_kv,
                             real_type sn_mva,
                             bool ac);

        void reset_results(){
            reset_results_tsc();
        }

        void set_q(const RealVect & reactive_mismatch,
                   const SolverBusIdVect & id_grid_to_solver,
                   bool ac,
                   const Eigen::VectorXi & total_gen_per_bus,
                   const RealVect & total_q_min_per_bus,
                   const RealVect & total_q_max_per_bus){
            side_1_.set_q(reactive_mismatch, id_grid_to_solver, ac, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
            side_2_.set_q(reactive_mismatch, id_grid_to_solver, ac, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
        }
        void get_vm_for_dc(RealVect & Vm){
            side_1_.get_vm_for_dc(Vm);
            side_2_.get_vm_for_dc(Vm);
        }
        void set_vm_or(CplxVect & V, const SolverBusIdVect & id_grid_to_solver) const{
            side_1_.set_vm(V, id_grid_to_solver);
        }
        void set_vm_ex(CplxVect & V, const SolverBusIdVect & id_grid_to_solver) const{
            side_2_.set_vm(V, id_grid_to_solver);
        }

        /**
        this functions makes sure that the voltage magnitude of every connected bus is properly used to initialize
        the ac powerflow
        **/
        void set_vm(CplxVect & V, const SolverBusIdVect & id_grid_to_solver) const{
            side_1_.set_vm(V, id_grid_to_solver);
            side_2_.set_vm(V, id_grid_to_solver);
        }

        const ConverterStationContainer & get_stations_side_1() const {return side_1_;}
        const ConverterStationContainer & get_stations_side_2() const {return side_2_;}

    private:
        /**
         * Recompute the two station active powers (and the LCC reactive
         * consumptions) from `converters_mode_` / `p_setpoint_mw_`.
         */
        void update_station_targets(int hvdc_id, DualAlgoControl & solver_control);
        /**
         * Legacy parametrization: given the side-1 station active power `p1`
         * (generator convention), derive `converters_mode_` / `p_setpoint_mw_`
         * (matches the legacy DCLineContainer::get_to_mw in both branches)
         * then update the station targets.
         */
        void update_targets_from_p1(int hvdc_id, real_type p1_mw, DualAlgoControl & solver_control);

    private:
        // line data (legacy pandapower loss model: multiplicative + fixed)
        RealVect loss_percent_;
        RealVect loss_mw_;

        // IIDM model
        IntVect converters_mode_;   // ConvertersMode
        RealVect p_setpoint_mw_;    // active power drawn at the rectifier (>= 0)
        RealVect r_ohm_;            // DC line resistance
        RealVect nominal_v_kv_;     // DC nominal voltage (for the r loss term)

        // angle-droop (AC emulation)
        std::vector<bool> droop_enabled_;
        RealVect p0_mw_;
        RealVect k_mw_per_rad_;
        RealVect pmax_1to2_mw_;
        RealVect pmax_2to1_mw_;
        IntVect status_droop_;      // INPUT: 0 linear, +1 sat 1->2, -1 sat 2->1
};

inline HvdcLineInfo::HvdcLineInfo(const HvdcLineContainer & r_data_hvdc, int my_id) noexcept:
    TwoSidesContainer<ConverterStationContainer>::TwoSidesInfo(r_data_hvdc, my_id),
    target_p_1_mw(0.),
    p_2_mw(0.),
    target_vm_1_pu(0.),
    target_vm_2_pu(0.),
    loss_pct(0.),
    loss_mw(0.),
    converters_mode(HvdcLineContainer::ConvertersMode::SIDE_1_RECTIFIER),
    p_setpoint_mw(0.),
    r_ohm(0.),
    nominal_v_kv(0.),
    droop_enabled(false),
    droop_p0_mw(0.),
    droop_k_mw_per_rad(0.),
    pmax_1to2_mw(0.),
    pmax_2to1_mw(0.),
    status_droop(0),
    station_side_1(r_data_hvdc.side_1_, my_id),
    station_side_2(r_data_hvdc.side_2_, my_id)
{
    if (my_id < 0) return;
    if (my_id >= static_cast<int>(r_data_hvdc.nb())) return;
    loss_pct = r_data_hvdc.loss_percent_(my_id);
    loss_mw = r_data_hvdc.loss_mw_(my_id);

    target_p_1_mw = station_side_1.target_p_mw;
    p_2_mw = station_side_2.target_p_mw;

    target_vm_1_pu = station_side_1.target_vm_pu;
    target_vm_2_pu = station_side_2.target_vm_pu;

    converters_mode = r_data_hvdc.converters_mode_(my_id);
    p_setpoint_mw = r_data_hvdc.p_setpoint_mw_(my_id);
    r_ohm = r_data_hvdc.r_ohm_(my_id);
    nominal_v_kv = r_data_hvdc.nominal_v_kv_(my_id);

    droop_enabled = r_data_hvdc.droop_enabled_[my_id];
    droop_p0_mw = r_data_hvdc.p0_mw_(my_id);
    droop_k_mw_per_rad = r_data_hvdc.k_mw_per_rad_(my_id);
    pmax_1to2_mw = r_data_hvdc.pmax_1to2_mw_(my_id);
    pmax_2to1_mw = r_data_hvdc.pmax_2to1_mw_(my_id);
    status_droop = r_data_hvdc.status_droop_(my_id);
}

} // namespace ls2g

#endif  //HVDCLINECONTAINER_H
