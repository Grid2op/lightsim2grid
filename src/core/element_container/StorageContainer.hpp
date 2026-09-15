// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef STORAGE_CONTAINER_H
#define STORAGE_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SlackParticipation.hpp"
#include "VoltageSourceContainer.hpp"

namespace ls2g {


class StorageContainer;
class LS2G_API StorageInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        bool voltage_regulator_on;
        real_type target_vm_pu;
        real_type min_q_mvar;
        real_type max_q_mvar;
        int regulated_bus_id;   // grid bus id whose voltage is regulated (== bus_id: local control, the only kind supported)
        bool is_slack;
        real_type slack_weight;

        inline StorageInfo(const StorageContainer & r_data_storage, int my_id) noexcept;
};


/**
This class is a container for all storage units (batteries) on the grid.

Storage units are modeled as PQ injections using the **load convention** (same as
in pandapower and grid2op): positive `target_p` means the unit is charging (power is
taken from the grid), negative `target_p` means the unit is discharging (power is
injected in the grid).

A storage unit can also **regulate the voltage of its own bus** (an IIDM battery
carrying a ``voltageRegulation`` extension, which OpenLoadFlow runs exactly like a
PV generator): it then holds `target_vm_pu` there and its reactive output is
solved for -- within `min_q_` / `max_q_` for the reactive-residual split -- instead
of being the `target_q_mvar` setpoint. The voltage side (regulating or not, the
setpoint, `set_vm`, the PV path, `set_q`) is VoltageSourceContainer's, shared with
the generators, converter stations and SVCs; what this leaf adds is the load
convention of its injection and of the reactive output written back to it. Only
LOCAL regulation is supported: a storage unit regulating another bus is refused
by `check_grid` (a remote controller would have to be enrolled in the
VoltageControl plan, which knows generators, SVCs and converter stations).

A storage unit can also **take part in the distributed slack**, like a generator
(OpenLoadFlow distributes the slack on batteries with the same rule): see
SlackParticipation. The share it absorbs is written back to its active result in
the load convention.

NOTE: PowSyBl / IIDM batteries use the opposite (generator) convention for their
`target_p` / `target_q` setpoints; the pypowsybl converter
(`init_from_pypowsybl`) negates them before feeding this container.

This is a dedicated container (rather than reusing :class:`LoadContainer`) so that
storage-specific behaviour can diverge from loads -- as the voltage regulation does.
**/
class LS2G_API StorageContainer final: public VoltageSourceContainer<StorageContainer>, public IteratorAdder<StorageContainer, StorageInfo>
{
    friend class StorageInfo;
    friend class VoltageSourceContainer<StorageContainer>;

    public:
        using DataInfo = StorageInfo;

    // regular implementation
    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
           OneSideContainer_PQ::StateRes,  // state of the base class
           std::vector<bool>,              // voltage_regulator_on_
           std::vector<real_type>,         // target_vm_pu_
           std::vector<real_type>,         // min_q_
           std::vector<real_type>,         // max_q_
           std::vector<int>,               // regulated_bus_id_ (== own bus)
           std::vector<bool>,              // slack participation flag
           std::vector<real_type>          // slack weight
           > ;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
            VREG_ON,
            TARGET_VM_PU,
            MIN_Q,
            MAX_Q,
            REGULATED_BUS_ID,
            SLACKBUS,
            SLACK_WEIGHT,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "StorageContainer::StateRes and StateResIdx do not match");

        StorageContainer() noexcept = default;
        ~StorageContainer() noexcept override = default;

        // pickle (python)
        StorageContainer::StateRes get_state() const;
        void set_state(StorageContainer::StateRes & my_state);

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static StorageContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "StorageContainer"; }  // written into / checked against the binary file header

        /// plain PQ storage units (no voltage regulation), the historical entry point
        void init(const Eigen::Ref<const RealVect> & storage_p_mw,
                  const Eigen::Ref<const RealVect> & storage_q_mvar,
                  const Eigen::Ref<const Eigen::VectorXi> & storage_bus_id
                  );

        /// same, plus the voltage side: which units regulate their bus, at what
        /// magnitude (pu), and their reactive range (MVAr, generator convention:
        /// `min_q <= max_q`, what the unit can inject) used to split the bus'
        /// reactive residual between the machines holding it
        void init_full(const Eigen::Ref<const RealVect> & storage_p_mw,
                       const Eigen::Ref<const RealVect> & storage_q_mvar,
                       const std::vector<bool> & voltage_regulator_on,
                       const Eigen::Ref<const RealVect> & storage_target_vm_pu,
                       const Eigen::Ref<const RealVect> & storage_min_q,
                       const Eigen::Ref<const RealVect> & storage_max_q,
                       const Eigen::Ref<const Eigen::VectorXi> & storage_bus_id
                       );

        real_type get_min_q(int storage_id) const {return min_q_.coeff(storage_id);}
        real_type get_max_q(int storage_id) const {return max_q_.coeff(storage_id);}

        /// the reactive output (MVAr, GENERATOR convention, what the split hands every
        /// machine) of a regulating unit, stored in this container's load convention.
        /// Hides VoltageSourceContainer's, which LSGrid reaches statically.
        void set_voltage_control_q(int storage_id, real_type q_mvar) {res_q_(storage_id) = -q_mvar;}

        // ---- distributed slack (see SlackParticipation; LSGrid validates the ids) ----
        void add_slackbus(int storage_id, real_type weight, DualAlgoControl & solver_control){
            slack_.add(storage_id, weight, solver_control, "StorageContainer::add_slackbus");
        }
        void remove_slackbus(int storage_id, DualAlgoControl & solver_control){
            slack_.remove(storage_id, solver_control);
        }
        void remove_all_slackbus(){ slack_.remove_all(); }
        bool is_slack(int storage_id) const {return slack_.is_slack(storage_id);}
        /// the unit's own (un-normalised) share of the distributed slack
        real_type get_slack_weight(int storage_id) const {return slack_.weight(storage_id);}
        /// add every participating unit's raw weight to its solver bus
        void accumulate_slack_weights_solver(RealVect & res, const SolverBusIdVect & id_grid_to_solver) const {
            slack_.accumulate_raw(res, status_, bus_id_, id_grid_to_solver, nullptr, _element_name());
        }
        void append_slack_bus_id(std::vector<int> & buses) const {slack_.append_slack_buses(buses, bus_id_);}
        void slack_summary(bool & any_flagged, bool & any_connected) const {slack_.summary(status_, any_flagged, any_connected);}
        /// write the share of the slack each participating unit absorbed (load convention)
        void set_p_slack(const Eigen::Ref<const RealVect> & node_mismatch,
                         const SolverBusIdVect & id_grid_to_solver,
                         const Eigen::Ref<const RealVect> & bus_raw_total){
            slack_.split(res_p_, -1., node_mismatch, bus_raw_total, status_, bus_id_, id_grid_to_solver,
                         "StorageContainer::set_p_slack");
        }

    protected:
        // ---- what VoltageSourceContainer asks of its leaf -------------------------
        static real_type _vm_scale(real_type target_vm, real_type current_vm) { return (1.0 / current_vm) * target_vm; }
        static const char * _element_name() { return "storage"; }
        bool _treated_as_off(int /*storage_id*/) const { return false; }
        bool _set_vm_skips(int /*storage_id*/) const { return false; }
        static constexpr bool set_vm_throws_on_unresolved = true;

        // load convention: the active setpoint is drawn from the grid; the reactive
        // one too, unless the unit regulates its bus (its Q is then solved for)
        void _fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool /*ac*/) const override;

        // the voltage-source checks, plus "local regulation only" and the slack weights
        void _check_valid(int nb_bus,
                          int nb_sub,
                          const SubstationContainer & substations,
                          std::vector<int> & all_pos_topo_vect) const override;

        // the voltage-source flags plus the slack role
        void _on_deactivate(int storage_id, DualAlgoControl & solver_control) override final;
        void _on_reactivate(int storage_id, DualAlgoControl & solver_control) override final;
        void _on_change_bus(int storage_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control) override final;

    protected:
        bool _in_topo_vect() const override { return true; }

    private:
        // reactive range (MVAr, generator convention), read by the reactive-residual split
        RealVect min_q_;
        RealVect max_q_;

        // distributed slack participation
        SlackParticipation slack_;
};

inline StorageInfo::StorageInfo(const StorageContainer & r_data_storage, int my_id) noexcept:
        OneSidePQInfo(r_data_storage, my_id),
        voltage_regulator_on(false),
        target_vm_pu(0.),
        min_q_mvar(0.),
        max_q_mvar(0.),
        regulated_bus_id(-1),
        is_slack(false),
        slack_weight(-1.0)
{
    if((my_id >= 0) && (my_id < r_data_storage.nb()))
    {
        voltage_regulator_on = r_data_storage.voltage_regulator_on_[my_id];
        target_vm_pu = r_data_storage.target_vm_pu_.coeff(my_id);
        min_q_mvar = r_data_storage.min_q_.coeff(my_id);
        max_q_mvar = r_data_storage.max_q_.coeff(my_id);
        regulated_bus_id = r_data_storage.regulated_bus_id_(my_id);
        is_slack = r_data_storage.slack_.is_slack(my_id);
        slack_weight = r_data_storage.slack_.weight(my_id);
    }
}


} // namespace ls2g

#endif  //STORAGE_CONTAINER_H
