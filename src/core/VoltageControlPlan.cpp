// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "VoltageControlPlan.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "BaseConstants.hpp"
#include "element_container/GenericContainer.hpp"
#include "element_container/GeneratorContainer.hpp"
#include "element_container/HvdcLineContainer.hpp"
#include "element_container/StorageContainer.hpp"
#include "element_container/SvcContainer.hpp"

namespace ls2g {

void VoltageControlPlan::clear() noexcept
{
    group_reg_buses_.clear();
    free_vm_slack_buses_.clear();
    controllers_.clear();
}

// ---------------------------------------------------------------------------
// layer 1: which buses a GROUP regulates (grid ids, no labelling needed)
// ---------------------------------------------------------------------------
void VoltageControlPlan::build_groups(const GeneratorContainer & generators,
                                      const SvcContainer & svcs,
                                      bool supports_voltage_control,
                                      bool hold_frozen)
{
    group_reg_buses_.clear();
    // An algorithm with no bordered block cannot honour a group, and taking a bus
    // out of PV for one it will not build leaves that bus' magnitude pinned by
    // nothing at all. Leaving the set empty is what gives layer 2 the classical
    // split; refusing the grid outright, when it really does have controllers, is
    // LSGrid::ac_pf's job (see list_unsupported) and happens before this runs.
    if(!supports_voltage_control) return;
    // an ACTIVE remote-regulating generator: is_remote_voltage_controller() already
    // means "connected, regulator on, not pseudo-off, and regulating a bus that is
    // NOT its own". A purely local regulator therefore never lands here, which is
    // what keeps the ordinary (possibly multi-generator) PV bus untouched.
    const int nb_gen = static_cast<int>(generators.nb());
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators.is_remote_voltage_controller(gen_id)) continue;
        const int reg = generators.get_regulated_bus_id(gen_id);
        if(reg >= 0) group_reg_buses_.insert(reg);
    }
    // a frozen remote regulator kept, held, in the group it would join (see
    // LSGrid::set_hold_frozen_regulators): that group's voltage row needs the bus'
    // Vm unknown exactly as an active remote regulator's would
    if(hold_frozen){
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
            if(!generators.is_frozen_remote_regulator(gen_id)) continue;
            group_reg_buses_.insert(generators.get_regulated_bus_id(gen_id));
        }
    }
    // a voltage-mode SVC is ALWAYS a group controller (even local and non-sloped),
    // so the bus it regulates always needs the bordered treatment
    const int nb_svc = static_cast<int>(svcs.nb());
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svcs.is_voltage_controller(svc_id)) continue;
        const int reg = svcs.get_regulated_bus_id(svc_id);
        if(reg >= 0) group_reg_buses_.insert(reg);
    }
    // and NOT the hvdc converter stations, which is not an oversight. A station has no
    // regulated bus of its own to read -- ConverterStationContainer stores none, a
    // regulating station always regulates the bus it stands on -- so it can never be
    // the REMOTE controller that makes a bus need the bordered treatment. It pins its
    // own bus through the ordinary PV path (ConverterStationContainer::fillpv), exactly
    // like a local generator, and becomes a group member only when a group formed by
    // one of the two loops above already claims that bus. That enrolment is
    // build_controllers' job (_collect_station_controllers), and it is keyed on the set
    // this function returns -- which is why adding stations here would be circular as
    // well as wrong. Remote regulation BY a station is simply not modelled; see the
    // changelog TODO.
}

// ---------------------------------------------------------------------------
// layer 2: the pv/pq split
// ---------------------------------------------------------------------------
void VoltageControlPlan::build_pv_pq(const std::vector<const GenericContainer *> & pv_sources,
                                     const SolverBusIdVect & id_me_to_solver,
                                     const GlobalBusIdVect & id_solver_to_me,
                                     const SolverBusIdVect & slack_bus_id_solver,
                                     SolverBusIdVect & bus_pv_out,
                                     SolverBusIdVect & bus_pq_out) const
{
    const int nb_bus = static_cast<int>(id_solver_to_me.size());  // number of bus in the solver!
    std::vector<int> bus_pq;
    bus_pq.reserve(nb_bus);
    std::vector<int> bus_pv;
    bus_pv.reserve(nb_bus);
    std::vector<bool> has_bus_been_added(nb_bus, false);

    bus_pv_out = SolverBusIdVect();
    bus_pq_out = SolverBusIdVect();

    // the classical part: every container says which buses it pins
    for(const GenericContainer * container : pv_sources){
        if(container == nullptr) continue;
        container->fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    }

    // ... and layer 1 takes back the ones a control GROUP regulates: their magnitude
    // is set by the group's bordered voltage row, so they must keep their own Vm
    // unknown (and hence their Q equation). See the doc on this function for the
    // configuration this fixes. Whatever pinned such a bus through the PV path -- a
    // local generator, or a voltage-regulating hvdc converter station -- is enrolled
    // as a member of the group by build_controllers instead.
    if(!group_reg_buses_.empty()){
        std::vector<int> bus_pv_kept;
        bus_pv_kept.reserve(bus_pv.size());
        for(int bus_id_solver : bus_pv){
            const int bus_id_me = id_solver_to_me[bus_id_solver].cast_int();
            if(bus_id_me < 0 || !group_reg_buses_.count(bus_id_me)){
                bus_pv_kept.push_back(bus_id_solver);
                continue;
            }
            has_bus_been_added[bus_id_solver] = false;  // let the PQ loop take it
        }
        bus_pv.swap(bus_pv_kept);
    }

    // TODO remove the order here..., i could be faster in this piece of code
    // (looping once through the buses)
    for(int bus_id = 0; bus_id < nb_bus; ++bus_id){
        if(GenericContainer::is_in_vect(bus_id, slack_bus_id_solver.to_int_vector())) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_out = SolverBusIdVect(bus_pv.size(), SolverBusId(0));
    for(int i = 0; i < static_cast<int>(bus_pv.size()); ++i){
        bus_pv_out(i) = SolverBusId(bus_pv[i]);
    }
    bus_pq_out = SolverBusIdVect(bus_pq.size(), SolverBusId(0));
    for(int i = 0; i < static_cast<int>(bus_pq.size()); ++i){
        bus_pq_out(i) = SolverBusId(bus_pq[i]);
    }
}

// ---------------------------------------------------------------------------
// the guard: what an algorithm without the bordered block must refuse
// ---------------------------------------------------------------------------
VoltageControlPlan::Unsupported
VoltageControlPlan::list_unsupported(const GeneratorContainer & generators,
                                     const SvcContainer & svcs,
                                     const HvdcLineContainer & hvdc_lines) const
{
    Unsupported res;
    if(group_reg_buses_.empty()) return res;  // nothing needs the bordered block

    const int nb_gen = static_cast<int>(generators.nb());
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(generators.is_remote_voltage_controller(gen_id)) res.gen_ids.push_back(gen_id);
    }
    const int nb_svc = static_cast<int>(svcs.nb());
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(svcs.is_voltage_controller(svc_id)) res.svc_ids.push_back(svc_id);
    }
    // a station is never an offender on its own: it pins its own bus the classical
    // way unless a group already claims that bus, in which case it is enrolled and
    // the user needs to know it is affected too
    const int nb_hvdc = static_cast<int>(hvdc_lines.nb());
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            if(!hvdc_lines.station_is_voltage_controller(hvdc_id, side)) continue;
            const int bus = hvdc_lines.get_station_bus(hvdc_id, side).cast_int();
            if(group_reg_buses_.count(bus)) res.station_ids.push_back(std::make_pair(hvdc_id, side));
        }
    }
    return res;
}

// ---------------------------------------------------------------------------
// layers 3 and 4
// ---------------------------------------------------------------------------
void VoltageControlPlan::build_solver_side(const GeneratorContainer & generators,
                                           const StorageContainer & storages,
                                           const SvcContainer & svcs,
                                           const HvdcLineContainer & hvdc_lines,
                                           const SolverBusIdVect & id_me_to_solver,
                                           const GlobalBusIdVect & id_solver_to_me,
                                           const SolverBusIdVect & slack_bus_id_solver,
                                           const SolverBusIdVect & bus_pq,
                                           real_type sn_mva,
                                           bool hold_frozen)
{
    build_free_vm_slack(generators, storages, id_me_to_solver, id_solver_to_me, slack_bus_id_solver);
    build_controllers(generators, storages, svcs, hvdc_lines, id_me_to_solver, id_solver_to_me, bus_pq,
                      sn_mva, hold_frozen);
}

std::vector<int> VoltageControlPlan::group_controlled_solver_buses(const SolverBusIdVect & id_me_to_solver) const
{
    std::vector<int> res;
    res.reserve(group_reg_buses_.size());
    const int nb_bus_me = static_cast<int>(id_me_to_solver.size());
    for(const int bus_me : group_reg_buses_){
        if(bus_me < 0 || bus_me >= nb_bus_me) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver < 0) continue;  // not solved in this labelling (deactivated, other component)
        res.push_back(bus_solver);
    }
    return res;
}

void VoltageControlPlan::build_free_vm_slack(const GeneratorContainer & generators,
                                             const StorageContainer & storages,
                                             const SolverBusIdVect & id_me_to_solver,
                                             const GlobalBusIdVect & id_solver_to_me,
                                             const SolverBusIdVect & slack_bus_id_solver)
{
    free_vm_slack_buses_.clear();
    // Nothing has been labelled yet (a grid that never solved, or a cache just
    // retired): there is no solver bus to express any of this in, and every id below
    // would index an empty labelling. Layer 1 stays as it was -- it is expressed in
    // grid ids and does not depend on any of this.
    if(id_solver_to_me.size() == 0) return;

    // solver-bus ids of the slack buses
    std::set<int> slack;
    for(int k = 0; k < static_cast<int>(slack_bus_id_solver.size()); ++k){
        slack.insert(slack_bus_id_solver(k).cast_int());
    }
    if(slack.empty()) return;

    // A slack bus is Vm-fixed (PV-like, no Q equation) only when a LOCAL
    // voltage-regulating generator pins its magnitude. Collect those buses.
    std::set<int> locally_vfixed;
    // ... except that a local regulator on a bus a control GROUP regulates does not
    // pin it: it is enrolled as a member of that group instead (see build_groups and
    // the reclassification in build_pv_pq), and the group's voltage row needs
    // the free Vm this grants.
    const int nb_gen = static_cast<int>(generators.nb());
    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators.is_local_voltage_controller(gen_id)) continue;
        const int ctrl_grid = gen_buses(gen_id).cast_int();
        if(group_reg_buses_.count(ctrl_grid)) continue;
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        if(ctrl_solver == GenericContainer::_deactivated_bus_id) continue;
        locally_vfixed.insert(ctrl_solver);
    }
    // ... and so does a voltage-regulating storage unit, which pins its own bus through
    // the same PV path (it only regulates locally, see StorageContainer::_check_valid).
    // Same exception: build_pv_pq takes a group-regulated bus back from PV whatever
    // pinned it. Without this loop a slack bus held by a regulating battery -- a
    // storage participant of the distributed slack, or a PQ slack generator sharing
    // the bus -- would get a free Vm and the battery's setpoint would be ignored.
    const int nb_storage = static_cast<int>(storages.nb());
    const GlobalBusIdVect & storage_buses = storages.get_bus_id();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!storages.is_local_voltage_controller(storage_id)) continue;
        const int ctrl_grid = storage_buses(storage_id).cast_int();
        if(group_reg_buses_.count(ctrl_grid)) continue;
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        if(ctrl_solver == GenericContainer::_deactivated_bus_id) continue;
        locally_vfixed.insert(ctrl_solver);
    }

    // Every slack bus whose magnitude is NOT pinned locally needs a free Vm
    // unknown + Q equation: distributed-slack PQ participants (the common case),
    // remote-voltage controllers, and SVC-regulated slack buses all fall here.
    for(int b : slack){
        if(!locally_vfixed.count(b)) free_vm_slack_buses_.insert(b);
    }
}

void VoltageControlPlan::build_controllers(const GeneratorContainer & generators,
                                           const StorageContainer & storages,
                                           const SvcContainer & svcs,
                                           const HvdcLineContainer & hvdc_lines,
                                           const SolverBusIdVect & id_me_to_solver,
                                           const GlobalBusIdVect & id_solver_to_me,
                                           const SolverBusIdVect & bus_pq,
                                           real_type sn_mva,
                                           bool hold_frozen)
{
    controllers_.clear();
    const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
    if(nb_bus_solver == 0) return;  // see build_free_vm_slack

    // PQ membership: a bus owns a Q equation AND a Vm unknown iff it is a PQ bus
    // (PV buses have only theta/P, the slack none). `bus_pq` is layer 2's own output.
    std::vector<bool> is_pq(nb_bus_solver, false);
    for(int k = 0; k < static_cast<int>(bus_pq.size()); ++k){
        const int b = bus_pq(k).cast_int();
        if(b >= 0 && b < nb_bus_solver) is_pq[b] = true;
    }
    // Slack buses are not PQ in the base block, but a slack bus that is not pinned
    // by a local PV generator is given a Q equation + free Vm by the Base block (see
    // build_free_vm_slack), so a controller on such a slack bus is supported even
    // though `is_pq` is false there. A slack bus that IS locally pinned (another
    // generator regulates it directly) gets no such Q equation at all -- checking
    // membership of the whole slack list here (as opposed to just this "free"
    // subset) would wrongly accept that case: its Q equation lookup then resolves to
    // -1, the controller's own reactive-injection column ends up with no Jacobian
    // entry anywhere, and the factorization fails with ErrorType::SolverFactor
    // instead of this function's own clear error.
    std::vector<bool> has_free_q(nb_bus_solver, false);
    for(int b : free_vm_slack_buses_){
        if(b >= 0 && b < nb_bus_solver) has_free_q[b] = true;
    }

    // 1. collect the active voltage-mode controllers. Per controller: solver bus,
    //    regulated solver bus, v_set (pu), sharing key, kind, elem id.
    std::vector<Raw> raws;
    _collect_gen_controllers(generators, id_me_to_solver, is_pq, has_free_q, raws);
    _collect_svc_controllers(svcs, id_me_to_solver, is_pq, has_free_q, raws);
    _collect_station_controllers(hvdc_lines, id_me_to_solver, is_pq, has_free_q, raws);
    // last, so that within a group every held controller comes after the active ones
    if(hold_frozen) _collect_held_gen_controllers(generators, id_me_to_solver, is_pq, has_free_q, raws);
    if(raws.empty()) return;

    _group_and_emit(raws, _collect_passive_gens(generators, storages, svcs, hvdc_lines,
                                                id_me_to_solver, raws, sn_mva));
}

std::map<int, VoltageControlPlan::PassiveBus> VoltageControlPlan::_collect_passive_gens(
    const GeneratorContainer & generators,
    const StorageContainer & storages,
    const SvcContainer & svcs,
    const HvdcLineContainer & hvdc_lines,
    const SolverBusIdVect & id_me_to_solver,
    const std::vector<Raw> & raws,
    real_type sn_mva) const
{
    // OpenLoadFlow counts every LfGenerator of a controller bus: the generators, but also
    // the batteries, the VSC converter stations and the SVCs (an LCC station is a load
    // there). Only a generator can carry a reactive key, so any other unit makes its bus
    // unkeyed -- and the whole group falls back on the reactive ranges.
    std::map<int, PassiveBus> out;
    std::set<int> ctrl_buses;
    std::set<std::pair<int, int> > ctrl_units;  // (kind, elem_id) of the controllers
    for(const Raw & r : raws){
        ctrl_buses.insert(r.bus);
        ctrl_units.insert(std::make_pair(r.kind, r.elem_id));
    }
    const int nb_bus_solver = static_cast<int>(id_me_to_solver.size());
    // the controller bus a connected unit sits on, -1 if none
    const auto ctrl_bus_of = [&](int bus_me){
        if(bus_me < 0 || bus_me >= nb_bus_solver) return -1;
        const int bus = id_me_to_solver[bus_me].cast_int();
        return ctrl_buses.count(bus) ? bus : -1;
    };
    // a range nothing can share on (unset limits) adds nothing, but still unkeys the bus
    const auto add = [&](int bus, real_type range, real_type key){
        PassiveBus & p = out[bus];
        if(std::isfinite(range)) p.range += range;
        if(std::isfinite(key) && key > 0.) p.keys += key; else p.keyed = false;
    };
    const real_type no_key = std::numeric_limits<real_type>::quiet_NaN();

    const int nb_gen = static_cast<int>(generators.nb());
    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators.get_status(gen_id)) continue;
        if(ctrl_units.count(std::make_pair(static_cast<int>(VoltageControlSolverData::GEN), gen_id))) continue;
        const int bus = ctrl_bus_of(gen_buses(gen_id).cast_int());
        if(bus < 0) continue;
        add(bus, generators.get_max_q(gen_id) - generators.get_min_q(gen_id), generators.get_reactive_key(gen_id));
    }
    const int nb_storage = static_cast<int>(storages.nb());
    const GlobalBusIdVect & storage_buses = storages.get_bus_id();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!storages.get_status(storage_id)) continue;
        const int bus = ctrl_bus_of(storage_buses(storage_id).cast_int());
        if(bus < 0) continue;
        add(bus, storages.get_max_q(storage_id) - storages.get_min_q(storage_id), no_key);
    }
    const int nb_svc = static_cast<int>(svcs.nb());
    const GlobalBusIdVect & svc_buses = svcs.get_bus_id();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svcs.get_status(svc_id)) continue;
        if(ctrl_units.count(std::make_pair(static_cast<int>(VoltageControlSolverData::SVC), svc_id))) continue;
        const int bus = ctrl_bus_of(svc_buses(svc_id).cast_int());
        if(bus < 0) continue;
        // its susceptance range (pu) at 1 pu, in MVAr like the others
        add(bus, (svcs.get_b_max(svc_id) - svcs.get_b_min(svc_id)) * sn_mva, no_key);
    }
    const int nb_hvdc = static_cast<int>(hvdc_lines.nb());
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            const ConverterStationContainer & stations = (side == 1) ? hvdc_lines.get_stations_side_1()
                                                                     : hvdc_lines.get_stations_side_2();
            if(stations.is_lcc(hvdc_id) || !stations.get_status(hvdc_id)) continue;
            const int kind = (side == 1) ? VoltageControlSolverData::HVDC_SIDE_1
                                         : VoltageControlSolverData::HVDC_SIDE_2;
            if(ctrl_units.count(std::make_pair(kind, hvdc_id))) continue;
            const int bus = ctrl_bus_of(hvdc_lines.get_station_bus(hvdc_id, side).cast_int());
            if(bus < 0) continue;
            add(bus, hvdc_lines.get_station_q_range_mvar(hvdc_id, side), no_key);
        }
    }
    return out;
}

void VoltageControlPlan::_collect_gen_controllers(const GeneratorContainer & generators,
                                                  const SolverBusIdVect & id_me_to_solver,
                                                  const std::vector<bool> & is_pq,
                                                  const std::vector<bool> & has_free_q,
                                                  std::vector<Raw> & raws) const
{
    const int nb_gen = static_cast<int>(generators.nb());
    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        // remote regulators are always controllers; a local one only when the bus it
        // regulates is group-controlled (something else remote/an SVC aims at it too)
        if(!generators.is_remote_voltage_controller(gen_id)){
            if(!generators.is_local_voltage_controller(gen_id)) continue;
            if(!group_reg_buses_.count(generators.get_regulated_bus_id(gen_id))) continue;
        }
        const int ctrl_grid = gen_buses(gen_id).cast_int();
        const int reg_grid  = generators.get_regulated_bus_id(gen_id);
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        const int reg_solver  = (reg_grid >= 0) ? id_me_to_solver[reg_grid].cast_int()
                                                : GenericContainer::_deactivated_bus_id;
        if(ctrl_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " is a voltage controller but its bus is disconnected.";
            throw std::runtime_error(exc_.str());
        }
        if(reg_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates a disconnected bus.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[ctrl_solver] && !has_free_q[ctrl_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates a remote bus but its OWN bus has no reactive (Q) equation"
                    " (it is a PV bus that is not a slack, or a slack bus already locally"
                    " pinned by another voltage-regulating generator). This is not supported"
                    " in v1.";
            throw std::runtime_error(exc_.str());
        }
        // The regulated bus needs a Vm unknown for the bordered row to act on. An
        // ordinary PQ bus has one; so does a slack bus that nothing pins locally (same
        // escape hatch as the controller bus just above -- `has_free_q`). A bus that is
        // PV despite being group-controlled cannot occur any more (layer 2
        // reclassifies it), so what is left here is a bus pinned by something that
        // cannot be enrolled.
        if(!is_pq[reg_solver] && !has_free_q[reg_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates bus " << reg_grid << " which has no voltage (Vm) unknown"
                    " (its magnitude is pinned by something that cannot join a control"
                    " group). This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        const real_type w = generators.get_max_q(gen_id) - generators.get_min_q(gen_id);
        raws.push_back({ctrl_solver, reg_solver, generators.get_target_vm_pu(gen_id),
                        static_cast<real_type>(0.), w, VoltageControlSolverData::GEN, gen_id,
                        generators.get_reactive_key(gen_id), false});
    }
}

void VoltageControlPlan::_collect_held_gen_controllers(const GeneratorContainer & generators,
                                                       const SolverBusIdVect & id_me_to_solver,
                                                       const std::vector<bool> & is_pq,
                                                       const std::vector<bool> & has_free_q,
                                                       std::vector<Raw> & raws) const
{
    // The same rules as an active remote regulator, with "left out" in place of every
    // error: a held controller only exists to be released later, and a grid that solves
    // without it must keep solving with the option on.
    const int nb_gen = static_cast<int>(generators.nb());
    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    const int nb_bus_solver = static_cast<int>(is_pq.size());
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators.is_frozen_remote_regulator(gen_id)) continue;
        const int ctrl_grid = gen_buses(gen_id).cast_int();
        const int reg_grid  = generators.get_regulated_bus_id(gen_id);
        if(ctrl_grid < 0 || reg_grid < 0) continue;
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        const int reg_solver  = id_me_to_solver[reg_grid].cast_int();
        if(ctrl_solver < 0 || ctrl_solver >= nb_bus_solver) continue;
        if(reg_solver < 0 || reg_solver >= nb_bus_solver) continue;
        if(!is_pq[ctrl_solver] && !has_free_q[ctrl_solver]) continue;   // no Q equation of its own
        if(!is_pq[reg_solver] && !has_free_q[reg_solver]) continue;     // nothing to regulate
        const real_type w = generators.get_max_q(gen_id) - generators.get_min_q(gen_id);
        raws.push_back({ctrl_solver, reg_solver, generators.get_target_vm_pu(gen_id),
                        static_cast<real_type>(0.), w, VoltageControlSolverData::GEN, gen_id,
                        generators.get_reactive_key(gen_id), true});
    }
}

void VoltageControlPlan::_collect_svc_controllers(const SvcContainer & svcs,
                                                  const SolverBusIdVect & id_me_to_solver,
                                                  const std::vector<bool> & is_pq,
                                                  const std::vector<bool> & has_free_q,
                                                  std::vector<Raw> & raws) const
{
    // the active VOLTAGE-mode SVCs (local or remote, with or without slope)
    const int nb_svc = static_cast<int>(svcs.nb());
    const GlobalBusIdVect & svc_buses = svcs.get_bus_id();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svcs.is_voltage_controller(svc_id)) continue;
        const int ctrl_grid = svc_buses(svc_id).cast_int();
        const int reg_grid  = svcs.get_regulated_bus_id(svc_id);
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        const int reg_solver  = (reg_grid >= 0) ? id_me_to_solver[reg_grid].cast_int()
                                                : GenericContainer::_deactivated_bus_id;
        if(ctrl_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " is a voltage controller but its bus is disconnected.";
            throw std::runtime_error(exc_.str());
        }
        if(reg_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " regulates a disconnected bus.";
            throw std::runtime_error(exc_.str());
        }
        // same two escape hatches as the generator branch above: a slack bus that
        // nothing pins locally owns a Q equation and a free Vm
        if(!is_pq[ctrl_solver] && !has_free_q[ctrl_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " is at a bus with no reactive (Q) equation (it is a PV bus, or a slack"
                    " bus pinned by a local voltage-regulating generator)."
                    " This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[reg_solver] && !has_free_q[reg_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " regulates bus " << reg_grid << " which has no voltage (Vm) unknown"
                    " (its magnitude is pinned by something that cannot join a control"
                    " group). This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        const real_type w = svcs.get_b_max(svc_id) - svcs.get_b_min(svc_id);
        raws.push_back({ctrl_solver, reg_solver, svcs.get_target_vm_pu(svc_id),
                        svcs.get_slope_pu(svc_id), w, VoltageControlSolverData::SVC, svc_id});
    }
}

void VoltageControlPlan::_collect_station_controllers(const HvdcLineContainer & hvdc_lines,
                                                      const SolverBusIdVect & id_me_to_solver,
                                                      const std::vector<bool> & is_pq,
                                                      const std::vector<bool> & has_free_q,
                                                      std::vector<Raw> & raws) const
{
    // the voltage-regulating hvdc converter stations. A VSC station with
    // voltage_regulator_on pins its own bus through the PV path, exactly like a local
    // generator, so it is a controller only when a GROUP regulates that bus instead
    // (layer 2 then kept the bus out of PV and it needs its members). Its sharing
    // key is a reactive range in MVAr -- the same currency as a generator's -- so a
    // mixed generator/station group shares reactive power correctly, unlike an SVC
    // (whose key is a susceptance range; hence the SVC-alone restriction below).
    const int nb_hvdc = static_cast<int>(hvdc_lines.nb());
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            if(!hvdc_lines.station_is_voltage_controller(hvdc_id, side)) continue;
            const int ctrl_grid = hvdc_lines.get_station_bus(hvdc_id, side).cast_int();
            // not group-controlled: the station keeps pinning its bus the classical way
            if(!group_reg_buses_.count(ctrl_grid)) continue;
            if(ctrl_grid == GenericContainer::_deactivated_bus_id) continue;
            const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
            if(ctrl_solver == GenericContainer::_deactivated_bus_id){
                std::ostringstream exc_;
                exc_ << "LSGrid::fill_voltage_control_solver_data: hvdc line " << hvdc_id
                     << " side " << side << " regulates voltage but its bus is disconnected.";
                throw std::runtime_error(exc_.str());
            }
            if(!is_pq[ctrl_solver] && !has_free_q[ctrl_solver]){
                std::ostringstream exc_;
                exc_ << "LSGrid::fill_voltage_control_solver_data: hvdc line " << hvdc_id
                     << " side " << side << " takes part in a voltage-control group but its"
                        " bus has no reactive (Q) equation. This is not supported in v1.";
                throw std::runtime_error(exc_.str());
            }
            const int kind = (side == 1) ? VoltageControlSolverData::HVDC_SIDE_1
                                         : VoltageControlSolverData::HVDC_SIDE_2;
            raws.push_back({ctrl_solver, ctrl_solver,   // a station regulates its OWN bus
                            hvdc_lines.get_station_target_vm_pu(hvdc_id, side),
                            static_cast<real_type>(0.),
                            hvdc_lines.get_station_q_range_mvar(hvdc_id, side),
                            kind, hvdc_id});
        }
    }
}

void VoltageControlPlan::_group_and_emit(const std::vector<Raw> & raws,
                                         const std::map<int, PassiveBus> & passive)
{
    // 2. group by regulated solver bus (merge gens that share a regulated bus),
    //    checking the v_set agree within tolerance.
    std::vector<int> grp_reg;
    std::vector<real_type> grp_vset;
    std::vector<std::vector<int> > grp_members;  // indices into raws
    // the held controllers come last in `raws` (build_controllers): a group made by an
    // active controller is never joined by one whose set-point differs, nor one that
    // would share an SVC's group -- both left out, never an error
    auto has_svc = [&](int g){
        for(int idx : grp_members[g])
            if(raws[idx].kind == VoltageControlSolverData::SVC) return true;
        return false;
    };
    for(int i = 0; i < static_cast<int>(raws.size()); ++i){
        int g = -1;
        for(int gg = 0; gg < static_cast<int>(grp_reg.size()); ++gg)
            if(grp_reg[gg] == raws[i].reg_bus){ g = gg; break; }
        if(g == -1){
            g = static_cast<int>(grp_reg.size());
            grp_reg.push_back(raws[i].reg_bus);
            grp_vset.push_back(raws[i].v_set);
            grp_members.push_back(std::vector<int>());
        } else if(std::abs(grp_vset[g] - raws[i].v_set) > BaseConstants::_tol_equal_float){
            if(raws[i].held) continue;
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: several controllers regulate the"
                    " same bus with conflicting voltage setpoints (" << grp_vset[g] << " vs "
                 << raws[i].v_set << " pu).";
            throw std::runtime_error(exc_.str());
        } else if(raws[i].held && has_svc(g)){
            continue;
        }
        grp_members[g].push_back(i);
    }

    // 2b. v1 restriction: an SVC may only be ALONE in its control group. The
    //     cross-weight sharing of an SVC with other controllers (and any sloped
    //     SVC sharing a regulated bus, cf Phase 0 probe #3) is not supported yet.
    for(int g = 0; g < static_cast<int>(grp_members.size()); ++g){
        if(grp_members[g].size() <= 1) continue;
        for(int idx : grp_members[g]){
            if(raws[idx].kind == VoltageControlSolverData::SVC){
                std::ostringstream exc_;
                exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << raws[idx].elem_id
                     << " shares a regulated bus with other controllers, which is not"
                        " supported in v1 (an SVC must be the only controller of its bus).";
                throw std::runtime_error(exc_.str());
            }
        }
    }

    // 3. emit, controllers grouped contiguously
    const int ng = static_cast<int>(grp_reg.size());
    int nc = 0;
    bool any_held = false;
    for(const auto & m : grp_members){
        nc += static_cast<int>(m.size());
        for(int idx : m) any_held = any_held || raws[idx].held;
    }
    VoltageControlSolverData & data = controllers_;
    data.bus = Eigen::VectorXi(nc);
    data.kind = Eigen::VectorXi(nc);
    data.elem_id = Eigen::VectorXi(nc);
    data.slope = RealVect(nc);
    data.weight = RealVect(nc);
    data.group = Eigen::VectorXi(nc);
    if(any_held) data.held = Eigen::VectorXi::Zero(nc);
    data.reg_bus = Eigen::VectorXi(ng);
    data.v_set = RealVect(ng);
    data.grp_start = Eigen::VectorXi(ng);
    data.grp_count = Eigen::VectorXi(ng);
    // The sharing key, OpenLoadFlow's rule, which works bus by bus. The buses of the
    // group share by the sum of the keys of ALL their generators (a connected one that
    // controls nothing, or is held, counts too, and so do the batteries, VSC stations and
    // SVCs there, which have no key: `passive`) when every one of them has a key, by the
    // sum of their reactive ranges otherwise. Inside one controller bus
    // the controllers share that by their keys when they all have one, by their
    // reactive ranges otherwise. The sharing rows hold Q_i / w_i equal across the
    // group, so w_i = (share of its bus) * (its share inside the bus).
    // the sums over some controllers of one bus, accumulated in the order they are added
    struct BusSum {
        real_type keys = 0.;
        real_type range = 0.;
        bool keyed = true;
        void add(const Raw & o){
            range += o.weight;
            if(std::isfinite(o.key) && o.key > 0.) keys += o.key; else keyed = false;
        }
    };

    int cursor = 0;
    for(int g = 0; g < ng; ++g){
        data.reg_bus(g) = grp_reg[g];
        data.v_set(g) = grp_vset[g];
        data.grp_start(g) = cursor;
        data.grp_count(g) = static_cast<int>(grp_members[g].size());
        // The keys of the ACTIVE members are shared among themselves (a held one counts
        // in its bus' share, as it would were it not held), so that the held ones change
        // nothing of the system the grid poses; a held one gets the key it would have
        // among the active ones once released.
        // Everything but the controller's own membership is the same for the whole group,
        // so it is summed once here: the share of each bus (every member there, and its
        // passive generators), whether every bus of the active members is keyed, and the
        // sums inside each bus over its active members. A held controller's peers are the
        // active members plus itself, so it adds its own bus (to the keyed test) and its own
        // terms (to the sums inside its bus) on top of them.
        const std::vector<int> & members = grp_members[g];
        std::map<int, BusSum> bus_share;    // per controller bus, over the members
        for(int idx : members) bus_share[raws[idx].bus].add(raws[idx]);
        for(auto & bs : bus_share){
            const auto it = passive.find(bs.first);
            if(it != passive.end()){
                bs.second.keys += it->second.keys;
                bs.second.range += it->second.range;
                bs.second.keyed = bs.second.keyed && it->second.keyed;
            }
        }
        bool active_buses_keyed = true;
        std::map<int, BusSum> inside;       // per controller bus, over the active members
        for(int idx : members){
            if(raws[idx].held) continue;
            active_buses_keyed = bus_share[raws[idx].bus].keyed && active_buses_keyed;
            inside[raws[idx].bus].add(raws[idx]);
        }
        for(int idx : members){
            const Raw & r = raws[idx];
            data.bus(cursor) = r.bus;
            data.kind(cursor) = r.kind;
            data.elem_id(cursor) = r.elem_id;
            data.slope(cursor) = r.slope;
            if(r.held) data.held(cursor) = 1;
            // the share of its bus: by keys when every bus of its peers is keyed
            const BusSum & own_bus = bus_share[r.bus];
            const bool all_buses_keyed = active_buses_keyed && own_bus.keyed;
            const real_type share = all_buses_keyed ? own_bus.keys : own_bus.range;
            // its share inside the bus, among its peers there
            BusSum in_bus;
            const auto it = inside.find(r.bus);
            if(it != inside.end()) in_bus = it->second;
            if(r.held) in_bus.add(r);
            const real_type w = in_bus.keyed ? share * (r.key / in_bus.keys)
                                             : share * (r.weight / in_bus.range);
            // floor the sharing key to keep the N>1 sharing rows non-singular
            data.weight(cursor) = (std::abs(w) > BaseConstants::_tol_equal_float) ? w : BaseConstants::_tol_equal_float;
            data.group(cursor) = g;
            ++cursor;
        }
    }
}

} // namespace ls2g
