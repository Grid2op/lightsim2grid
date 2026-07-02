// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SVCCONTAINER_H
#define SVCCONTAINER_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"

namespace ls2g {

class SvcContainer;

class LS2G_API SvcInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        int regulation_mode;       // SvcContainer::RegulationMode
        real_type target_vm_pu;    // VOLTAGE mode
        real_type slope_pu;        // VOLTAGE mode (0 = no slope / droop)
        real_type b_min;           // stored, NEVER enforced
        real_type b_max;           // stored, NEVER enforced
        int regulated_bus_id;      // grid bus id whose voltage is regulated (VOLTAGE mode)

        inline SvcInfo(const SvcContainer & r_data_svc, int my_id) noexcept;
};

/**
This class represents the list of Static Var Compensators (SVC).

An SVC injects reactive power only (its active power is ALWAYS 0). It follows the
IIDM model of powsybl: three regulation modes
  - VOLTAGE        : regulates the voltage of a bus (local or remote), optionally
                     with a voltage/reactive slope ("droop"). It stamps NOTHING in
                     Sbus, is NEVER PV, and is ALWAYS a controller of a
                     VoltageControl group (the bordered formulation), even for the
                     local non-sloped case;
  - REACTIVE_POWER : a fixed reactive injection (behaves like a non-regulating
                     generator / a load): stamps Q (and P = 0) into Sbus;
  - OFF            : behaves as if disconnected (nothing anywhere).

`b_min_` / `b_max_` are stored for introspection but NEVER enforced (no outer
loop, no limit check), mirroring the generator Qmin/Qmax handling.
**/
class LS2G_API SvcContainer : public OneSideContainer_PQ, public IteratorAdder<SvcContainer, SvcInfo>
{
    friend class SvcInfo;

    public:
        using DataInfo = SvcInfo;

        enum RegulationMode {
            OFF = 0,
            VOLTAGE = 1,
            REACTIVE_POWER = 2
        };

        using StateRes = std::tuple<
           OneSideContainer_PQ::StateRes,
           std::vector<int>,        // regulation_mode_
           std::vector<real_type>,  // target_vm_pu_
           std::vector<real_type>,  // slope_pu_
           std::vector<real_type>,  // b_min_
           std::vector<real_type>,  // b_max_
           std::vector<int>         // regulated_bus_id_
        >;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
            REGULATION_MODE,
            TARGET_VM_PU,
            SLOPE_PU,
            B_MIN,
            B_MAX,
            REGULATED_BUS_ID,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "SvcContainer::StateRes and StateResIdx do not match");

        SvcContainer() noexcept = default;
        virtual ~SvcContainer() noexcept = default;

        void init(const std::vector<int> & regulation_mode,
                  const RealVect & target_vm_pu,
                  const RealVect & q_setpoint_mvar,
                  const RealVect & slope_pu,
                  const RealVect & b_min,
                  const RealVect & b_max,
                  const Eigen::VectorXi & regulated_bus_id,
                  const Eigen::VectorXi & bus_id);

        // pickle
        SvcContainer::StateRes get_state() const;
        void set_state(SvcContainer::StateRes & my_state);

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path) const;
        static SvcContainer load_binary(const std::string & path);

        // accessors
        int get_regulation_mode(int svc_id) const {return regulation_mode_(svc_id);}
        real_type get_target_vm_pu(int svc_id) const {return target_vm_pu_(svc_id);}
        real_type get_slope_pu(int svc_id) const {return slope_pu_.coeff(svc_id);}
        real_type get_b_min(int svc_id) const {return b_min_.coeff(svc_id);}
        real_type get_b_max(int svc_id) const {return b_max_.coeff(svc_id);}
        int get_regulated_bus_id(int svc_id) const {return regulated_bus_id_(svc_id);}
        bool regulates_remote(int svc_id) const {
            return regulated_bus_id_(svc_id) != bus_id_(svc_id).cast_int();
        }
        // true iff this SVC is an ACTIVE voltage-mode controller (joins a group)
        bool svc_is_voltage_controller(int svc_id) const {
            return status_[svc_id] && (regulation_mode_(svc_id) == RegulationMode::VOLTAGE);
        }
        bool has_slope(int svc_id) const {
            return std::abs(slope_pu_.coeff(svc_id)) > _tol_equal_float;
        }
        // write the converged reactive output (MVAr) of a voltage-mode SVC,
        // supplied by the VoltageControl extension (LSGrid::compute_results)
        void set_voltage_control_q(int svc_id, real_type q_mvar) {res_q_(svc_id) = q_mvar;}

        // solver interface
        virtual void fillSbus(CplxVect & Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;
        void get_vm_for_dc(RealVect & Vm);
        void set_vm(CplxVect & V, const SolverBusIdVect & id_grid_to_solver) const;

    protected:
        virtual void _compute_results(
            const Eigen::Ref<const RealVect> & Va,
            const Eigen::Ref<const RealVect> & Vm,
            const Eigen::Ref<const CplxVect> & V,
            const SolverBusIdVect & id_grid_to_solver,
            const RealVect & bus_vn_kv,
            real_type sn_mva,
            bool ac) override;
        bool _deactivate(int svc_id, DualAlgoControl & solver_control) final;
        bool _reactivate(int svc_id, DualAlgoControl & solver_control) final;
        bool _change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int nb_bus) final;

    private:
        IntVect regulation_mode_;             // RegulationMode, per SVC
        RealVect target_vm_pu_;               // VOLTAGE mode
        RealVect slope_pu_;                   // VOLTAGE mode (0 = no slope)
        RealVect b_min_;                      // stored, NEVER enforced
        RealVect b_max_;                      // stored, NEVER enforced
        Eigen::VectorXi regulated_bus_id_;    // grid bus id (defaults to own bus = local)
};

inline SvcInfo::SvcInfo(const SvcContainer & r_data_svc, int my_id) noexcept:
OneSidePQInfo(r_data_svc, my_id),
regulation_mode(SvcContainer::RegulationMode::OFF),
target_vm_pu(0.),
slope_pu(0.),
b_min(0.),
b_max(0.),
regulated_bus_id(-1)
{
    if((my_id >= 0) && (my_id < r_data_svc.nb()))
    {
        regulation_mode = r_data_svc.regulation_mode_(my_id);
        target_vm_pu = r_data_svc.target_vm_pu_.coeff(my_id);
        slope_pu = r_data_svc.slope_pu_.coeff(my_id);
        b_min = r_data_svc.b_min_.coeff(my_id);
        b_max = r_data_svc.b_max_.coeff(my_id);
        regulated_bus_id = r_data_svc.regulated_bus_id_(my_id);
    }
}

} // namespace ls2g

#endif  //SVCCONTAINER_H
