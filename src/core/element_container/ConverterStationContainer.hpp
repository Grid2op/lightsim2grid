// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CONVERTERSTATIONCONTAINER_H
#define CONVERTERSTATIONCONTAINER_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "Utils.hpp"
#include "VoltageSourceContainer.hpp"

namespace ls2g {

class ConverterStationContainer;

class LS2G_API ConverterStationInfo final : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        int converter_type;  // 0 = VSC, 1 = LCC (ConverterStationContainer::ConverterType)
        real_type loss_factor;  // fraction (0 - 1)

        bool voltage_regulator_on;
        real_type target_vm_pu;
        real_type min_q_mvar;
        real_type max_q_mvar;
        real_type power_factor;  // LCC only

        inline ConverterStationInfo(const ConverterStationContainer & r_data_station, int my_id) noexcept;
};

/**
This class represents the list of HVDC converter stations (one container per HVDC line side).

It follows the IIDM model of powsybl: a station is either a VSC (voltage source
converter: behaves like a generator, either regulating voltage - PV - or with a
fixed reactive setpoint - PQ) or a LCC (line commutated converter: behaves like
a load, always consuming Q = |P| * tan(acos(power_factor))).

The active power of a station (`target_p_mw_` of the PQ base, generator sign
convention) is NOT an input: it is derived and written by the owning
`HvdcLineContainer` (from the line active power setpoint and the loss model).
Same for the LCC reactive consumption (stored in `target_q_mvar_`).

Reactive "personality" (input data): `voltage_regulator_on_` + `target_vm_pu_`
+ [`min_q_`, `max_q_`] when regulating (VSC only), `target_q_mvar_` otherwise.
The voltage side is VoltageSourceContainer's; a station always regulates its
own bus (`regulated_bus_id_` is not an input and follows the bus), which is
what leaves remote regulation by a station for later.

A station with a (derived) active power of exactly 0 MW is considered
"pseudo off": it is not PV and does not contribute to Sbus when regulating
(this mirrors the behaviour of the legacy DC lines, whose embedded generators
were configured with `turnedoff_no_pv`).
**/
class LS2G_API ConverterStationContainer final : public VoltageSourceContainer<ConverterStationContainer>, public IteratorAdder<ConverterStationContainer, ConverterStationInfo>
{
    friend class ConverterStationInfo;
    friend class VoltageSourceContainer<ConverterStationContainer>;

    public:
        using DataInfo = ConverterStationInfo;

        // /!\ these integer values are serialized verbatim (binary files and
        // pickles, as type_): renumbering requires bumping
        // BINARY_FORMAT_VERSION (BinaryArchive.hpp). Guarded python side by
        // TestSerializedEnumValues in test_binary_serialization.py.
        enum ConverterType {
            VSC = 0,
            LCC = 1
        };

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
           OneSideContainer_PQ::StateRes,
           std::vector<int>,        // type_
           std::vector<real_type>,  // loss_factor_
           std::vector<bool>,       // voltage_regulator_on_
           std::vector<real_type>,  // target_vm_pu_
           std::vector<real_type>,  // min_q_
           std::vector<real_type>,  // max_q_
           std::vector<real_type>   // power_factor_
        >;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
            TYPE,
            LOSS_FACTOR,
            VREG_ON,
            TARGET_VM_PU,
            MIN_Q,
            MAX_Q,
            POWER_FACTOR,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "ConverterStationContainer::StateRes and StateResIdx do not match");

        ConverterStationContainer() noexcept = default;
        ~ConverterStationContainer() noexcept override = default;

        void init(const std::vector<int> & type,
                  const Eigen::Ref<const RealVect> & loss_factor,
                  const std::vector<bool> & voltage_regulator_on,
                  const Eigen::Ref<const RealVect> & target_vm_pu,
                  const Eigen::Ref<const RealVect> & q_setpoint_mvar,
                  const Eigen::Ref<const RealVect> & min_q,
                  const Eigen::Ref<const RealVect> & max_q,
                  const Eigen::Ref<const RealVect> & power_factor,
                  const Eigen::Ref<const Eigen::VectorXi> & bus_id);

        // pickle
        ConverterStationContainer::StateRes get_state() const;
        void set_state(ConverterStationContainer::StateRes & my_state);

        // accessors (used by HvdcLineContainer, which owns the active power)
        bool is_lcc(int station_id) const {return type_(station_id) == ConverterType::LCC;}
        real_type get_loss_factor(int station_id) const {return loss_factor_(station_id);}
        real_type get_min_q(int station_id) const {return min_q_.coeff(station_id);}
        real_type get_max_q(int station_id) const {return max_q_.coeff(station_id);}
        real_type get_target_p(int station_id) const {return target_p_mw_(station_id);}

        /**
         * Set the (derived) station active power, generator convention.
         * For LCC stations the reactive consumption is updated accordingly.
         * Only called by the owning HvdcLineContainer.
         */
        void set_station_p(int station_id, real_type p_mw, DualAlgoControl & solver_control);

        /**
         * Add the station injection to Sbus.
         * `skip_p` says, per station, that the active power is NOT handled here
         * (it is handled by the HVDC droop extension of the NR system, or by
         * the DC algorithm): only the reactive personality is stamped then.
         * This is why the station's own `_fillSbus` hook stays empty: the owning
         * HvdcLineContainer, which knows the droop regime, stamps its stations.
         */
        void fillSbus_station(Eigen::Ref<CplxVect> Sbus,
                              const SolverBusIdVect & id_grid_to_solver,
                              bool ac,
                              const std::vector<bool> & skip_p) const;

    protected:
        // ---- what VoltageSourceContainer asks of its leaf -------------------------
        static const char * _element_name() { return "converter station"; }
        bool _treated_as_off(int /*station_id*/) const { return false; }
        bool _set_vm_skips(int /*station_id*/) const { return false; }
        static constexpr bool set_vm_throws_on_unresolved = true;

        void _on_change_p(int station_id, real_type new_p, DualAlgoControl & solver_control) override final;

    private:
        // input data
        IntVect type_;                           // ConverterType, per station
        RealVect loss_factor_;                   // fraction (0 - 1)
        RealVect min_q_;                         // when regulating
        RealVect max_q_;                         // when regulating
        RealVect power_factor_;                  // LCC only
};

inline ConverterStationInfo::ConverterStationInfo(const ConverterStationContainer & r_data_station, int my_id) noexcept:
OneSidePQInfo(r_data_station, my_id),
converter_type(ConverterStationContainer::ConverterType::VSC),
loss_factor(0.),
voltage_regulator_on(false),
target_vm_pu(0.),
min_q_mvar(0.),
max_q_mvar(0.),
power_factor(1.)
{
    if((my_id >= 0) && (my_id < r_data_station.nb()))
    {
        converter_type = r_data_station.type_(my_id);
        loss_factor = r_data_station.loss_factor_(my_id);
        voltage_regulator_on = r_data_station.voltage_regulator_on_[my_id];
        target_vm_pu = r_data_station.target_vm_pu_.coeff(my_id);
        min_q_mvar = r_data_station.min_q_.coeff(my_id);
        max_q_mvar = r_data_station.max_q_.coeff(my_id);
        power_factor = r_data_station.power_factor_.coeff(my_id);
    }
}

} // namespace ls2g

#endif  //CONVERTERSTATIONCONTAINER_H
