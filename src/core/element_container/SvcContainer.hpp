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
#include "VoltageSourceContainer.hpp"

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

The voltage side is VoltageSourceContainer's, with `voltage_regulator_on_`
derived from the mode (VOLTAGE) -- the mode has no setter, so the two cannot
drift -- and the one thing an SVC does differently: it never takes the PV path.

`b_min_` / `b_max_` are stored for introspection but NEVER enforced (no outer
loop, no limit check), mirroring the generator Qmin/Qmax handling.
**/
class LS2G_API SvcContainer final : public VoltageSourceContainer<SvcContainer>, public IteratorAdder<SvcContainer, SvcInfo>
{
    friend class SvcInfo;
    friend class VoltageSourceContainer<SvcContainer>;

    public:
        using DataInfo = SvcInfo;

        // /!\ these integer values are serialized verbatim (binary files and
        // pickles, as regulation_mode_): renumbering requires bumping
        // BINARY_FORMAT_VERSION (BinaryArchive.hpp). Guarded python side by
        // TestSerializedEnumValues in test_binary_serialization.py.
        enum RegulationMode {
            OFF = 0,
            VOLTAGE = 1,
            REACTIVE_POWER = 2
        };

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

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
        ~SvcContainer() noexcept override = default;

        void init(const std::vector<int> & regulation_mode,
                  const Eigen::Ref<const RealVect> & target_vm_pu,
                  const Eigen::Ref<const RealVect> & q_setpoint_mvar,
                  const Eigen::Ref<const RealVect> & slope_pu,
                  const Eigen::Ref<const RealVect> & b_min,
                  const Eigen::Ref<const RealVect> & b_max,
                  const Eigen::Ref<const Eigen::VectorXi> & regulated_bus_id,
                  const Eigen::Ref<const Eigen::VectorXi> & bus_id);

        // pickle
        SvcContainer::StateRes get_state() const;
        void set_state(SvcContainer::StateRes & my_state);

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static SvcContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "SvcContainer"; }  // written into / checked against the binary file header

        // accessors
        real_type get_slope_pu(int svc_id) const {return slope_pu_.coeff(svc_id);}
        real_type get_b_min(int svc_id) const {return b_min_.coeff(svc_id);}
        real_type get_b_max(int svc_id) const {return b_max_.coeff(svc_id);}
        bool has_slope(int svc_id) const {
            return std::abs(slope_pu_.coeff(svc_id)) > _tol_equal_float;
        }

    protected:
        // ---- what VoltageSourceContainer asks of its leaf -------------------------
        static const char * _element_name() { return "svc"; }
        bool _treated_as_off(int /*svc_id*/) const { return false; }
        // a sloped SVC does NOT hold its regulated bus exactly at the setpoint
        // (Vm = v_set - s.Q): forcing it there would corrupt an already-solved V
        // (e.g. in check_solution). The flat init is good enough for the NR.
        bool _set_vm_skips(int svc_id) const { return has_slope(svc_id); }
        // an SVC may point at a bus outside the solved system: skipped, not an error
        static constexpr bool set_vm_throws_on_unresolved = false;

        // an SVC is NEVER PV: in VOLTAGE mode it is always a group controller,
        // local and non-sloped included (the bordered formulation)
        void _fillpv(std::vector<int> & /*bus_pv*/,
                     std::vector<bool> & /*has_bus_been_added*/,
                     const SolverBusIdVect & /*slack_bus_id_solver*/,
                     const SolverBusIdVect & /*id_grid_to_solver*/) const override final {}

        void _fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;
        // P is always 0; Q is the setpoint in REACTIVE_POWER mode, the write-back
        // in VOLTAGE mode, nothing when OFF
        void _compute_res_pq(
            const Eigen::Ref<const RealVect> & Va,
            const Eigen::Ref<const RealVect> & Vm,
            const Eigen::Ref<const CplxVect> & V,
            const SolverBusIdVect & id_grid_to_solver,
            const Eigen::Ref<const RealVect> & bus_vn_kv,
            real_type sn_mva,
            bool ac) override;
        // the voltage-source flags, plus `one_el_changed_bus` (kept from before the
        // shared base: an SVC's status change raised it, a generator's did not)
        void _on_deactivate(int svc_id, DualAlgoControl & solver_control) override final;
        void _on_reactivate(int svc_id, DualAlgoControl & solver_control) override final;

    private:
        // voltage_regulator_on_ (the base's) is `regulation_mode_ == VOLTAGE`
        void _derive_voltage_regulator_on();

        IntVect regulation_mode_;             // RegulationMode, per SVC
        RealVect slope_pu_;                   // VOLTAGE mode (0 = no slope)
        RealVect b_min_;                      // stored, NEVER enforced
        RealVect b_max_;                      // stored, NEVER enforced
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
