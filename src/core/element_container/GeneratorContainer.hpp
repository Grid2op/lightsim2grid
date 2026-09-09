// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENERATORCONTAINER_H
#define GENERATORCONTAINER_H

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "VoltageSourceContainer.hpp"

namespace ls2g {

class GeneratorContainer;
class LS2G_API GenInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        bool is_slack;
        real_type slack_weight;

        bool voltage_regulator_on;
        real_type target_vm_pu;
        real_type min_q_mvar;
        real_type max_q_mvar;
        int regulated_bus_id;   // grid bus id whose voltage is regulated (== bus_id for local control)

        inline GenInfo(const GeneratorContainer & r_data_gen, int my_id) noexcept;
};

/**
This class represents the list of all generators.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/gen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/gen.html#electric-model

The voltage side (regulating or not, the setpoint, the regulated bus, the PV
path, the voltage initialisation) is VoltageSourceContainer's; what is specific
here is the active side -- the distributed slack -- and the `turnedoff_gen_pv_`
rule, which decides whether a generator at 0 MW still holds its voltage.
**/
class LS2G_API GeneratorContainer final: public VoltageSourceContainer<GeneratorContainer>, public IteratorAdder<GeneratorContainer, GenInfo>
{
    friend class GenInfo;
    friend class VoltageSourceContainer<GeneratorContainer>;

    public:
        using DataInfo = GenInfo;

    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
           OneSideContainer_PQ::StateRes,
           bool,                    // turnedoff_gen_pv_
           std::vector<bool>,       // voltage_regulator_on
           std::vector<real_type>,  // target_vm_pu_
           std::vector<real_type>,  // min_q_
           std::vector<real_type>,  // max_q_
           std::vector<bool>,       // gen_slackbus
           std::vector<real_type>,  // gen_slack_weight_
           std::vector<int>         // regulated_bus_id_ (appended; defaults to own bus)
        > ;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
            TURNEDOFF_GEN_PV,
            VREG_ON,
            TARGET_VM_PU,
            MIN_Q,
            MAX_Q,
            GEN_SLACKBUS,
            GEN_SLACK_WEIGHT,
            REGULATED_BUS_ID,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "GeneratorContainer::StateRes and StateResIdx do not match");

        GeneratorContainer() noexcept :VoltageSourceContainer<GeneratorContainer>(), turnedoff_gen_pv_(true){};
        explicit GeneratorContainer(bool turnedoff_gen_pv) noexcept :VoltageSourceContainer<GeneratorContainer>(), turnedoff_gen_pv_(turnedoff_gen_pv) {};
        ~GeneratorContainer() noexcept override = default;

        // TODO add pmin and pmax here !
        void init(const Eigen::Ref<const RealVect> & generators_p,
                  const Eigen::Ref<const RealVect> & generators_v,
                  const Eigen::Ref<const RealVect> & generators_min_q,
                  const Eigen::Ref<const RealVect> & generators_max_q,
                  const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id
                  );

        void init_full(const Eigen::Ref<const RealVect> & generators_p,
                       const Eigen::Ref<const RealVect> & generators_v,
                       const Eigen::Ref<const RealVect> & generators_q,
                       const std::vector<bool> & voltage_regulator_on,
                       const Eigen::Ref<const RealVect> & generators_min_q,
                       const Eigen::Ref<const RealVect> & generators_max_q,
                       const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id
                       );

        // pickle
        GeneratorContainer::StateRes get_state() const;
        void set_state(GeneratorContainer::StateRes & my_state );

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static GeneratorContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "GeneratorContainer"; }  // written into / checked against the binary file header

        // slack handling
        /**
        we suppose that the data are correct (ie gen_id in the proper range, and weight > 0.)
        This is checked in GridModel, and not at this stage
        **/
        void add_slackbus(int gen_id, real_type weight, DualAlgoControl & solver_control){
            // TODO DEBUG MODE
            if(weight <= 0.) throw std::runtime_error("GeneratorContainer::add_slackbus Cannot assign a negative (<=0) weight to the slack bus.");
            // Why `tell_slack_participate_changed()` and nothing about voltage: taking
            // the slack role is an ACTIVE-power role, but the pv/pq split reads the
            // slack set directly -- `fillpv` skips a bus that is in
            // `slack_bus_id_solver` ("slack bus is not PV"), and the PQ loop right
            // after it skips it too -- so a bus joining or leaving the slack set moves
            // the split, whatever it does to voltage. That flag is what says so, and
            // it is measured, not assumed: see the `[pv_pq]` cases in
            // test_cache_reuse.cpp, which reach a state where it is the only term
            // raised and fail if it is dropped.
            //
            // The voltage side needs no flag of its own for the same reason. It is
            // real but secondary -- `is_pseudo_off()` answers false for a slack
            // generator whatever its active power, so a zero-P slack generator IS a
            // voltage controller where an ordinary one would not be -- and the
            // voltage-control plan is built around the split this same flag rebuilds.
            // See AlgoControl::need_recompute_voltage_control.
            if(!gen_slackbus_[gen_id]){ solver_control.tell_slack_participate_changed(); }
            gen_slackbus_[gen_id] = true;
            if(abs(gen_slack_weight_[gen_id] - weight) > _tol_equal_float){
                solver_control.tell_slack_weight_changed();
                gen_slack_weight_[gen_id] = weight;
            }
        }
        void remove_slackbus(int gen_id, DualAlgoControl & solver_control){
            if(gen_slackbus_[gen_id]){ solver_control.tell_slack_participate_changed(); }
            if(abs(gen_slack_weight_[gen_id]) > _tol_equal_float){ solver_control.tell_slack_weight_changed(); }
            gen_slackbus_[gen_id] = false;
            gen_slack_weight_[gen_id] = 0.;
        }
        void remove_all_slackbus(){
            const int nb_gen = nb();
            DualAlgoControl unused_solver_control;
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
            {
                remove_slackbus(gen_id, unused_solver_control);
            }
        }

        // returns only the gen_id with the highest p that is connected to this bus !
        int assign_slack_bus(int slack_bus_id,
                             const std::vector<real_type> & gen_p_per_bus,
                             DualAlgoControl & solver_control){
            const int nb_gen = nb();
            int res_gen_id = -1;
            real_type max_p = -1.;
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
            {
                if(!status_[gen_id]) continue;
                if(bus_id_(gen_id).cast_int() != slack_bus_id) continue;
                const real_type p_mw = target_p_mw_(gen_id);
                if (p_mw > 0.) add_slackbus(gen_id, p_mw / gen_p_per_bus[slack_bus_id], solver_control);
                if((p_mw > max_p) || (res_gen_id == -1) ){
                    res_gen_id = gen_id;
                    max_p = p_mw;
                }
            }
            // TODO DEBUG MODE
            if(res_gen_id == -1) throw std::runtime_error("GeneratorContainer::assign_slack_bus No generator connected to the desired buses");
            return res_gen_id;
        }

        /**
        Retrieve the normalized (=sum to 1.000) slack weights for all the buses
        **/
        RealVect get_slack_weights_solver(size_t nb_bus_solver, const SolverBusIdVect & id_grid_to_solver);
        /**
         * Same normalised per-solver-bus slack weights as get_slack_weights_solver
         * above, but evaluated as if the generators flagged in `gen_off` (sized nb())
         * were disconnected -- what a batch sweep needs to re-weight the distributed
         * slack for a row whose contingency takes a participating machine out.
         *
         * Shares get_slack_weights_solver's rules through _raw_slack_weights_solver,
         * so the two can never drift apart. Unlike it, this one is const and does NOT
         * refresh the cached bus_slack_weight_ (that cache belongs to the grid's own
         * solve, not to a hypothetical row). Returns an ALL-ZERO vector when every
         * participating generator is off -- there is no meaningful normalisation
         * then, and it is the caller's business to decide what to do about it.
         */
        RealVect get_slack_weights_solver_without(size_t nb_bus_solver,
                                                  const SolverBusIdVect & id_grid_to_solver,
                                                  const std::vector<bool> & gen_off) const;

        GlobalBusIdVect get_slack_bus_id() const;
        /** distribute the active mismatch of the slack buses onto the participating generators **/
        void set_p_slack(const Eigen::Ref<const RealVect>& node_mismatch, const SolverBusIdVect & id_grid_to_solver);

        // modification
        void turnedoff_no_pv(DualAlgoControl & solver_control){
            solver_control.tell_slack_participate_changed();
            solver_control.tell_slack_weight_changed();
            turnedoff_gen_pv_=false;  // turned off generators are not pv. This is NOT the default.
            }
        void turnedoff_pv(DualAlgoControl & solver_control){
            solver_control.tell_slack_participate_changed();
            solver_control.tell_slack_weight_changed();
            turnedoff_gen_pv_=true;  // turned off generators are pv. This is the default.
            }
        bool get_turnedoff_gen_pv() const {return turnedoff_gen_pv_;}
        void update_slack_weights(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & could_be_slack,
                                  DualAlgoControl & solver_control);
        void update_slack_weights_by_id(const Eigen::Ref<const IntVect> & gen_slack_id, DualAlgoControl & solver_control);

        real_type get_min_q(int gen_id) const {return min_q_.coeff(gen_id);}
        real_type get_max_q(int gen_id) const {return max_q_.coeff(gen_id);}
        // the reactive setpoint a NON voltage-regulating generator injects (a
        // regulating one's reactive output is solved for, not set -- see fillSbus,
        // which only stamps this when voltage_regulator_on_ is false)
        real_type get_target_q_mvar(int gen_id) const {return target_q_mvar_(gen_id);}
        // the generator's own (un-normalised) share of the distributed slack, as
        // aggregated per bus by get_slack_weights_solver
        real_type get_gen_slack_weight(int gen_id) const {return gen_slack_weight_[gen_id];}

        /**
         * pseudo off generator (with p == 0) and with no contribution to the the slack bus
         */
        bool is_pseudo_off(int gen_id) const{
            if (gen_slackbus_[gen_id]) return false;  // slack is not pseudo off
            if ((abs(gen_slack_weight_[gen_id]) >= _tol_equal_float)) return false;  // slack is not pseudo off
            // pseudo-off <=> target_p == 0.
            return (abs(target_p_mw_(gen_id)) < _tol_equal_float);
        }

    protected:
        // ---- what VoltageSourceContainer asks of its leaf -------------------------
        static const char * _element_name() { return "generator"; }
        // a generator at 0 MW is switched off, voltage included, when `turnedoff_no_pv`
        // is set (except a slack generator, which is never pseudo off)
        bool _treated_as_off(int gen_id) const { return (!turnedoff_gen_pv_) && is_pseudo_off(gen_id); }
        bool _set_vm_skips(int /*gen_id*/) const { return false; }
        static constexpr bool set_vm_throws_on_unresolved = true;

        bool _in_topo_vect() const override { return true; }

        // the voltage-source checks plus the generator-specific ones -- slack weights
        void _check_valid(int nb_bus,
                          int nb_sub,
                          const SubstationContainer & substations,
                          std::vector<int> & all_pos_topo_vect) const override;

        void _fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;

        // the voltage-source flags plus the slack role and the `turnedoff_no_pv` rule
        void _on_change_p(int gen_id, real_type new_p, DualAlgoControl & solver_control) override final;
        void _on_deactivate(int gen_id, DualAlgoControl & solver_control) override final;
        void _on_reactivate(int gen_id, DualAlgoControl & solver_control) override final;
        void _on_change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & solver_control) override final;

        /**
         * The un-normalised per-solver-bus slack weights: one place holding the rule
         * for who participates (connected, flagged a slack bus, non-zero weight) and
         * the disconnected-bus checks. `gen_off`, when non-null, is a nb()-sized mask
         * of generators to leave out on top of that.
         */
        RealVect _raw_slack_weights_solver(size_t nb_bus_solver,
                                           const SolverBusIdVect & id_grid_to_solver,
                                           const std::vector<bool> * gen_off) const;

    private:
        // physical properties
        RealVect min_q_;
        RealVect max_q_;

        // remember which generators are "slack bus"
        std::vector<bool> gen_slackbus_;  // say for each generator if it's a slack or not
        std::vector<real_type> gen_slack_weight_;

        // intermediate data
        RealVect bus_slack_weight_;  // do not sum to 1., for each node of the grid, say the raw contribution for the generator

        // different parameter of the behaviour of the class
        bool turnedoff_gen_pv_;  // are turned off generators (including one with p=0) pv ?
};

inline GenInfo::GenInfo(const GeneratorContainer & r_data_gen, int my_id) noexcept:
OneSidePQInfo(r_data_gen, my_id),
is_slack(false),
slack_weight(-1.0),
voltage_regulator_on(false),
target_vm_pu(0.),
min_q_mvar(0.),
max_q_mvar(0.),
regulated_bus_id(-1)
{
    if((my_id >= 0) && (my_id < r_data_gen.nb()))
    {
        is_slack = r_data_gen.gen_slackbus_[my_id];
        slack_weight = r_data_gen.gen_slack_weight_[my_id];

        voltage_regulator_on = r_data_gen.voltage_regulator_on_[my_id];
        target_vm_pu = r_data_gen.target_vm_pu_.coeff(my_id);
        min_q_mvar = r_data_gen.min_q_.coeff(my_id);
        max_q_mvar = r_data_gen.max_q_.coeff(my_id);
        regulated_bus_id = r_data_gen.regulated_bus_id_(my_id);
    }
}


} // namespace ls2g

#endif  //GENERATORCONTAINER_H
