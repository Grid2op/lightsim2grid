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
#include "OneSideContainer_PQ.hpp"

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

A station with a (derived) active power of exactly 0 MW is considered
"pseudo off": it is not PV and does not contribute to Sbus when regulating
(this mirrors the behaviour of the legacy DC lines, whose embedded generators
were configured with `turnedoff_no_pv`).
**/
class LS2G_API ConverterStationContainer : public OneSideContainer_PQ, public IteratorAdder<ConverterStationContainer, ConverterStationInfo>
{
    friend class ConverterStationInfo;

    public:
        using DataInfo = ConverterStationInfo;

        enum ConverterType {
            VSC = 0,
            LCC = 1
        };

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
        virtual ~ConverterStationContainer() noexcept = default;

        void init(const std::vector<int> & type,
                  const RealVect & loss_factor,
                  const std::vector<bool> & voltage_regulator_on,
                  const RealVect & target_vm_pu,
                  const RealVect & q_setpoint_mvar,
                  const RealVect & min_q,
                  const RealVect & max_q,
                  const RealVect & power_factor,
                  const Eigen::VectorXi & bus_id);

        // pickle
        ConverterStationContainer::StateRes get_state() const;
        void set_state(ConverterStationContainer::StateRes & my_state);

        // accessors (used by HvdcLineContainer, which owns the active power)
        bool is_vsc(int station_id) const {return type_(station_id) == ConverterType::VSC;}
        bool is_lcc(int station_id) const {return type_(station_id) == ConverterType::LCC;}
        real_type get_loss_factor(int station_id) const {return loss_factor_(station_id);}
        bool get_voltage_regulator_on(int station_id) const {return voltage_regulator_on_[station_id];}
        real_type get_power_factor(int station_id) const {return power_factor_(station_id);}
        real_type get_qmin(int station_id) const {return min_q_.coeff(station_id);}
        real_type get_qmax(int station_id) const {return max_q_.coeff(station_id);}
        real_type get_target_p(int station_id) const {return target_p_mw_(station_id);}
        real_type get_target_vm_pu(int station_id) const {return target_vm_pu_(station_id);}

        /**
         * Set the (derived) station active power, generator convention.
         * For LCC stations the reactive consumption is updated accordingly.
         * Only called by the owning HvdcLineContainer.
         */
        void set_station_p(int station_id, real_type p_mw, AlgoControl & solver_control);

        void change_v(int station_id, real_type new_v_pu, AlgoControl & solver_control);

        // solver interface (mirrors GeneratorContainer, with the converter personality)
        /**
         * Add the station injection to Sbus.
         * `skip_p` says, per station, that the active power is NOT handled here
         * (it is handled by the HVDC droop extension of the NR system, or by
         * the DC algorithm): only the reactive personality is stamped then.
         */
        void fillSbus_station(CplxVect & Sbus,
                              const SolverBusIdVect & id_grid_to_solver,
                              bool ac,
                              const std::vector<bool> & skip_p) const;
        virtual void fillpv(std::vector<int>& bus_pv,
                            std::vector<bool> & has_bus_been_added,
                            const SolverBusIdVect & slack_bus_id_solver,
                            const SolverBusIdVect & id_grid_to_solver) const;
        void init_q_vector(int nb_bus,
                           Eigen::VectorXi & total_gen_per_bus,
                           RealVect & total_q_min_per_bus,
                           RealVect & total_q_max_per_bus) const;
        void set_q(const RealVect & reactive_mismatch,
                   const SolverBusIdVect & id_grid_to_solver,
                   bool ac,
                   const Eigen::VectorXi & total_gen_per_bus,
                   const RealVect & total_q_min_per_bus,
                   const RealVect & total_q_max_per_bus);
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

        void _change_p(int station_id, real_type new_p, bool my_status, DualAlgoControl & solver_control) final;
        bool _deactivate(int station_id, DualAlgoControl & solver_control) final;
        bool _reactivate(int station_id, DualAlgoControl & solver_control) final;
        bool _change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control, int nb_bus) final;

        /**
         * Station with a derived active power of 0: not PV, no Sbus contribution
         * when regulating (legacy dcline behaviour, cf `turnedoff_no_pv`).
         */
        bool is_pseudo_off(int station_id) const {
            return (abs(target_p_mw_(station_id)) < _tol_equal_float);
        }
        bool is_regulating(int station_id) const {
            return voltage_regulator_on_[station_id];
        }

    private:
        // input data
        IntVect type_;                           // ConverterType, per station
        RealVect loss_factor_;                   // fraction (0 - 1)
        std::vector<bool> voltage_regulator_on_; // VSC only, always false for LCC
        RealVect target_vm_pu_;                  // when regulating
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
