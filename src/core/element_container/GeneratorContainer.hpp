// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENERATORCONTAINER_H
#define GENERATORCONTAINER_H

#include <cmath>
#include <limits>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SlackParticipation.hpp"
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
        // active power limits, in MW -- OPTIONAL: NaN when the grid was never given any
        // (see LSGrid::set_gen_p_limits). Nothing enforces them; they are what says whether
        // the active power a distributed slack ended up asking of this machine is one it
        // could actually deliver (see batch_algorithm/GenPCheck.hpp).
        real_type min_p_mw;
        real_type max_p_mw;
        // a hint from the caller that this PQ machine is one an outer loop pinned at a
        // reactive limit, so that its PQ -> PV release is worth checking (see
        // LSGrid::set_gen_can_be_pv); false by default, never read by a powerflow
        bool can_be_pv;
        int regulated_bus_id;   // grid bus id whose voltage is regulated (== bus_id for local control)
        real_type reactive_key; // reactive sharing key, NaN when there is none

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
here is the active side -- the distributed slack, whose rules are
SlackParticipation's, shared with the storage units -- and the `turnedoff_gen_pv_`
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
           std::vector<int>,        // regulated_bus_id_ (appended; defaults to own bus)
           std::vector<real_type>,  // p_min_mw_ (appended, optional: empty if unset)
           std::vector<real_type>,  // p_max_mw_ (appended, optional: empty if unset)
           std::vector<real_type>,  // reactive_key_ (appended; NaN: no key)
           std::vector<bool>        // can_be_pv_ (appended; all false by default)
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
            P_MIN_MW,
            P_MAX_MW,
            REACTIVE_KEY,
            CAN_BE_PV,
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
        This is checked in GridModel, and not at this stage.
        See SlackParticipation::add for the flags this raises, and why.
        **/
        void add_slackbus(int gen_id, real_type weight, DualAlgoControl & solver_control){
            slack_.add(gen_id, weight, solver_control, "GeneratorContainer::add_slackbus");
        }
        void remove_slackbus(int gen_id, DualAlgoControl & solver_control){
            slack_.remove(gen_id, solver_control);
        }
        void remove_all_slackbus(){ slack_.remove_all(); }
        /// the participants outside the main component leave the slack (LSGrid::consider_only_main_component)
        void remove_slackbus_not_in_main_component(const std::vector<bool> & busbar_in_main_component,
                                                   DualAlgoControl & solver_control){
            slack_.remove_if_bus_not_in(busbar_in_main_component, bus_id_, solver_control);
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
         * Add every participating generator's raw (un-normalised) slack weight to its
         * solver bus in `res`. `gen_off`, when non-null, is a nb()-sized mask of
         * generators to evaluate as if they were disconnected -- what a batch sweep
         * needs to re-weight the distributed slack for a row whose contingency takes a
         * participating machine out. LSGrid adds the storage units' on top and
         * normalises (LSGrid::get_slack_weights_solver_without).
         */
        void accumulate_slack_weights_solver(RealVect & res,
                                             const SolverBusIdVect & id_grid_to_solver,
                                             const std::vector<bool> * gen_off) const {
            slack_.accumulate_raw(res, status_, bus_id_, id_grid_to_solver, gen_off, _element_name());
        }
        /// append the grid buses of the flagged generators not in `buses` yet
        void append_slack_bus_id(std::vector<int> & buses) const {slack_.append_slack_buses(buses, bus_id_);}
        void slack_summary(bool & any_flagged, bool & any_connected) const {slack_.summary(status_, any_flagged, any_connected);}
        /** distribute the active mismatch of the slack buses onto the participating generators **/
        void set_p_slack(const Eigen::Ref<const RealVect> & node_mismatch,
                         const SolverBusIdVect & id_grid_to_solver,
                         const Eigen::Ref<const RealVect> & bus_raw_total){
            slack_.split(res_p_, 1., node_mismatch, bus_raw_total, status_, bus_id_, id_grid_to_solver,
                         "GeneratorContainer::set_p_slack");
        }

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

        /**
         * Active power limits (MW), OPTIONAL -- exactly like a branch's thermal rating
         * (BranchContainer::set_limit_a1_ka): nothing in the powerflow reads them, they are
         * what a limit check compares against, and a grid that was never given any simply
         * has none (the two vectors stay empty, and `get_min_p` / `get_max_p` answer NaN).
         *
         * They matter because the distributed slack is solved INSIDE the Newton system
         * (`MultiSlack`), with fixed participation factors and no notion of a limit: a
         * participating machine's converged active power is `target_p + its share of the
         * imbalance`, which can land anywhere. See batch_algorithm/GenPCheck.hpp.
         *
         * Pass two empty vectors to drop them again.
         */
        void set_p_limits(const Eigen::Ref<const RealVect> & p_min_mw,
                          const Eigen::Ref<const RealVect> & p_max_mw){
            if((p_min_mw.size() == 0) && (p_max_mw.size() == 0)){
                p_min_mw_ = RealVect();
                p_max_mw_ = RealVect();
                return;
            }
            check_size(p_min_mw, nb(), "GeneratorContainer::set_p_limits (p_min_mw)");
            check_size(p_max_mw, nb(), "GeneratorContainer::set_p_limits (p_max_mw)");
            p_min_mw_ = p_min_mw;
            p_max_mw_ = p_max_mw;
        }
        Eigen::Ref<const RealVect> get_p_min_mw() const {return p_min_mw_;}
        Eigen::Ref<const RealVect> get_p_max_mw() const {return p_max_mw_;}
        /// NaN where no limit was given -- for the whole grid (never set) or for that one
        /// machine (a NaN in the vector handed to set_p_limits)
        real_type get_min_p(int gen_id) const {
            return p_min_mw_.size() > 0 ? p_min_mw_.coeff(gen_id)
                                        : std::numeric_limits<real_type>::quiet_NaN();
        }
        real_type get_max_p(int gen_id) const {
            return p_max_mw_.size() > 0 ? p_max_mw_.coeff(gen_id)
                                        : std::numeric_limits<real_type>::quiet_NaN();
        }
        /**
         * Flag, per generator, the machines an outer loop pinned at a reactive limit as
         * PQ (see LSGrid::set_gen_can_be_pv). Nothing in a powerflow reads it -- no
         * solver flag to raise -- it only says which PQ machines the physical checks may
         * report as "would regulate again". One entry per generator.
         */
        void set_can_be_pv(const std::vector<bool> & can_be_pv){
            check_size(can_be_pv, nb(), "GeneratorContainer::set_can_be_pv");
            can_be_pv_ = can_be_pv;
        }
        bool get_can_be_pv(int gen_id) const {return can_be_pv_[gen_id];}
        const std::vector<bool> & get_can_be_pv() const {return can_be_pv_;}
        // reactive sharing key among the generators holding one bus together (see
        // VoltageControlPlan::build_controllers); no key -- NaN, 0 or negative -- lets
        // the reactive range decide
        real_type get_reactive_key(int gen_id) const {return reactive_key_.coeff(gen_id);}
        void set_reactive_key(int gen_id, real_type key, DualAlgoControl & solver_control){
            _check_in_range(gen_id, reactive_key_, "set_reactive_key");
            const real_type old_key = reactive_key_(gen_id);
            if((std::isnan(old_key) && std::isnan(key)) || old_key == key) return;
            reactive_key_(gen_id) = key;
            // the sharing weights of the voltage-control plan (AC only); the pattern of
            // the sharing rows does not depend on them
            if(voltage_regulator_on_[gen_id]) solver_control.ac_algo_controler().tell_voltage_control_changed();
        }
        // the reactive setpoint a NON voltage-regulating generator injects (a
        // regulating one's reactive output is solved for, not set -- see fillSbus,
        // which only stamps this when voltage_regulator_on_ is false)
        real_type get_target_q_mvar(int gen_id) const {return target_q_mvar_(gen_id);}
        // the generator's own (un-normalised) share of the distributed slack, as
        // aggregated per bus by accumulate_slack_weights_solver
        real_type get_gen_slack_weight(int gen_id) const {return slack_.weight(gen_id);}
        /// the same, under the name every slack-participating container answers to (the
        /// storage units have one too), so that code checking both families can be written
        /// once -- see batch_algorithm/GenPCheck.hpp
        real_type get_slack_weight(int gen_id) const {return slack_.weight(gen_id);}
        bool is_slack(int gen_id) const {return slack_.is_slack(gen_id);}

        /**
         * pseudo off generator (with p == 0) and with no contribution to the the slack bus
         */
        bool is_pseudo_off(int gen_id) const{
            if (slack_.is_slack(gen_id)) return false;  // slack is not pseudo off
            if (slack_.has_weight(gen_id)) return false;  // slack is not pseudo off
            // pseudo-off <=> target_p == 0.
            return (abs(target_p_mw_(gen_id)) < _tol_equal_float);
        }

    protected:
        // ---- what VoltageSourceContainer asks of its leaf -------------------------
        static real_type _vm_scale(real_type target_vm, real_type current_vm) { return (1.0 / current_vm) * target_vm; }
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

    private:
        // physical properties
        RealVect min_q_;
        RealVect max_q_;
        // active power limits (MW), optional: empty when the grid was never given any
        RealVect p_min_mw_;
        RealVect p_max_mw_;
        RealVect reactive_key_;  // reactive sharing key, NaN when there is none
        // the PQ machines a caller knows an outer loop pinned at a reactive limit (see
        // set_can_be_pv); all false unless told otherwise, never read by a powerflow
        std::vector<bool> can_be_pv_;

        // which generators take part in the distributed slack, and with what weight
        SlackParticipation slack_;

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
min_p_mw(std::numeric_limits<real_type>::quiet_NaN()),
max_p_mw(std::numeric_limits<real_type>::quiet_NaN()),
can_be_pv(false),
regulated_bus_id(-1),
reactive_key(std::numeric_limits<real_type>::quiet_NaN())
{
    if((my_id >= 0) && (my_id < r_data_gen.nb()))
    {
        is_slack = r_data_gen.slack_.is_slack(my_id);
        slack_weight = r_data_gen.slack_.weight(my_id);

        voltage_regulator_on = r_data_gen.voltage_regulator_on_[my_id];
        target_vm_pu = r_data_gen.target_vm_pu_.coeff(my_id);
        min_q_mvar = r_data_gen.min_q_.coeff(my_id);
        max_q_mvar = r_data_gen.max_q_.coeff(my_id);
        min_p_mw = r_data_gen.get_min_p(my_id);
        max_p_mw = r_data_gen.get_max_p(my_id);
        can_be_pv = r_data_gen.can_be_pv_[my_id];
        regulated_bus_id = r_data_gen.regulated_bus_id_(my_id);
        reactive_key = r_data_gen.reactive_key_.coeff(my_id);
    }
}


} // namespace ls2g

#endif  //GENERATORCONTAINER_H
