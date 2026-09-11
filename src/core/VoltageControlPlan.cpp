// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "VoltageControlPlan.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "BaseConstants.hpp"
#include "element_container/GenericContainer.hpp"
#include "element_container/GeneratorContainer.hpp"
#include "element_container/HvdcLineContainer.hpp"
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
                                      bool supports_voltage_control)
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
                                           const SvcContainer & svcs,
                                           const HvdcLineContainer & hvdc_lines,
                                           const SolverBusIdVect & id_me_to_solver,
                                           const GlobalBusIdVect & id_solver_to_me,
                                           const SolverBusIdVect & slack_bus_id_solver,
                                           const SolverBusIdVect & bus_pq)
{
    build_free_vm_slack(generators, id_me_to_solver, id_solver_to_me, slack_bus_id_solver);
    build_controllers(generators, svcs, hvdc_lines, id_me_to_solver, id_solver_to_me, bus_pq);
}

void VoltageControlPlan::build_free_vm_slack(const GeneratorContainer & generators,
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

    // Every slack bus whose magnitude is NOT pinned locally needs a free Vm
    // unknown + Q equation: distributed-slack PQ participants (the common case),
    // remote-voltage controllers, and SVC-regulated slack buses all fall here.
    for(int b : slack){
        if(!locally_vfixed.count(b)) free_vm_slack_buses_.insert(b);
    }
}

void VoltageControlPlan::build_controllers(const GeneratorContainer & generators,
                                           const SvcContainer & svcs,
                                           const HvdcLineContainer & hvdc_lines,
                                           const SolverBusIdVect & id_me_to_solver,
                                           const GlobalBusIdVect & id_solver_to_me,
                                           const SolverBusIdVect & bus_pq)
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
    if(raws.empty()) return;

    _group_and_emit(raws);
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
                        static_cast<real_type>(0.), w, VoltageControlSolverData::GEN, gen_id});
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

void VoltageControlPlan::_group_and_emit(const std::vector<Raw> & raws)
{
    // 2. group by regulated solver bus (merge gens that share a regulated bus),
    //    checking the v_set agree within tolerance.
    std::vector<int> grp_reg;
    std::vector<real_type> grp_vset;
    std::vector<std::vector<int> > grp_members;  // indices into raws
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
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: several controllers regulate the"
                    " same bus with conflicting voltage setpoints (" << grp_vset[g] << " vs "
                 << raws[i].v_set << " pu).";
            throw std::runtime_error(exc_.str());
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
    const int nc = static_cast<int>(raws.size());
    VoltageControlSolverData & data = controllers_;
    data.bus = Eigen::VectorXi(nc);
    data.kind = Eigen::VectorXi(nc);
    data.elem_id = Eigen::VectorXi(nc);
    data.slope = RealVect(nc);
    data.weight = RealVect(nc);
    data.group = Eigen::VectorXi(nc);
    data.reg_bus = Eigen::VectorXi(ng);
    data.v_set = RealVect(ng);
    data.grp_start = Eigen::VectorXi(ng);
    data.grp_count = Eigen::VectorXi(ng);
    int cursor = 0;
    for(int g = 0; g < ng; ++g){
        data.reg_bus(g) = grp_reg[g];
        data.v_set(g) = grp_vset[g];
        data.grp_start(g) = cursor;
        data.grp_count(g) = static_cast<int>(grp_members[g].size());
        for(int idx : grp_members[g]){
            const Raw & r = raws[idx];
            data.bus(cursor) = r.bus;
            data.kind(cursor) = r.kind;
            data.elem_id(cursor) = r.elem_id;
            data.slope(cursor) = r.slope;
            // floor the sharing key to keep the N>1 sharing rows non-singular
            data.weight(cursor) = (std::abs(r.weight) > BaseConstants::_tol_equal_float) ? r.weight : BaseConstants::_tol_equal_float;
            data.group(cursor) = g;
            ++cursor;
        }
    }
}

} // namespace ls2g
