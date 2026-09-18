// Copyright (c) 2024-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ONE_SIDE_CONTAINER_PQ_H
#define ONE_SIDE_CONTAINER_PQ_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer.hpp"

namespace ls2g {


/**
 * A one-sided element with an active and a reactive setpoint: loads, static
 * generators, storage units, shunts, and (through VoltageSourceContainer) the
 * generators, SVCs and HVDC converter stations.
 *
 * On top of OneSideContainer it owns `target_p_mw_` / `target_q_mvar_` and their
 * setters, and provides what most such elements share: the default results
 * (`_compute_res_pq` publishes the setpoints), the injection stamp
 * (`_stamp_pq`, used by the leaves' `_fillSbus` with their own sign), and the
 * flags a status or bus change raises for an element that lives in Sbus only.
 * A leaf that is also in Ybus (a shunt) or that pins a bus (a generator)
 * overrides the corresponding `_on_xxx` hook.
 */
class OneSideContainer_PQ : public OneSideContainer
{

    public:
        class OneSidePQInfo: public OneSideContainer::OneSideInfo
        {
            public:
                real_type target_p_mw;
                real_type target_q_mvar;

                OneSidePQInfo(const OneSideContainer_PQ & r_data_pq, int my_id) noexcept:
                OneSideInfo(r_data_pq, my_id),
                target_p_mw(0.),
                target_q_mvar(0.)
                {
                    if (my_id < 0) return;
                    if (my_id >= r_data_pq.nb()) return;

                    target_p_mw = r_data_pq.target_p_mw_.coeff(my_id);
                    target_q_mvar = r_data_pq.target_q_mvar_.coeff(my_id);
                }
        };
    
    // regular implementation
    public:
        OneSideContainer_PQ() noexcept = default;
        ~OneSideContainer_PQ() noexcept override = default;

        // public generic API

        Eigen::Ref<const RealVect> get_target_p() const {return target_p_mw_;}
        Eigen::Ref<const RealVect> get_target_q() const {return target_q_mvar_;}

    protected:
        void _gen_p_per_bus(std::vector<real_type> & res) const override
        {
            const int nb_gen = nb();
            for(int sgen_id = 0; sgen_id < nb_gen; ++sgen_id)
            {
                if(!status_[sgen_id]) continue;
                const GlobalBusId my_bus = bus_id_(sgen_id);
                res[my_bus.cast_int()] += target_p_mw_(sgen_id);
            }
        }

    public:
        void change_p_nothrow(int el_id, real_type new_p, DualAlgoControl & solver_control)
        {
            _check_in_range(el_id, status_, "change_p");
            // notified BEFORE the write: the hooks compare the old and the new value
            _on_change_p(el_id, new_p, solver_control);
            if (abs(target_p_mw_(el_id) - new_p) > _tol_equal_float) {
                target_p_mw_(el_id) = new_p;
            }
        }
        void change_q_nothrow(int el_id, real_type new_q, DualAlgoControl & solver_control)
        {
            _check_in_range(el_id, status_, "change_q");
            _on_change_q(el_id, new_q, solver_control);
            if (abs(target_q_mvar_(el_id) - new_q) > _tol_equal_float) {
                target_q_mvar_(el_id) = new_q;
            }
        }

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            OneSideContainer::StateRes,
            std::vector<real_type>, // p_mw
            std::vector<real_type> // q_mvar
            >;
        enum StateResIdx {
            OSC_STATE = 0,
            TARGET_P_MW,
            TARGET_Q_MVAR,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "OneSideContainer_PQ::StateRes and StateResIdx do not match");

    protected:
        OneSideContainer_PQ::StateRes get_osc_pq_state() const  // osc: one side element
        {
            std::vector<real_type> target_p_mw(target_p_mw_.begin(), target_p_mw_.end());
            std::vector<real_type> target_q_mvar(target_q_mvar_.begin(), target_q_mvar_.end());
            OneSideContainer_PQ::StateRes res(
                get_osc_state(),
                target_p_mw,
                target_q_mvar);
            return res;
        }

        void set_osc_pq_state(OneSideContainer_PQ::StateRes & my_state)  // osc: one side element
        {
            // read data from my_state
            set_osc_state(std::get<StateResIdx::OSC_STATE>(my_state));

            // init target_p and target_q
            std::vector<real_type> & p_mw = std::get<StateResIdx::TARGET_P_MW>(my_state);
            std::vector<real_type> & q_mvar = std::get<StateResIdx::TARGET_Q_MVAR>(my_state);

            // check sizes
            const auto size = nb();
            check_size(p_mw, size, "p_mw");
            check_size(q_mvar, size, "q_mvar");

            // input data
            target_p_mw_ = RealVect::Map(p_mw.data(), p_mw.size());
            target_q_mvar_ = RealVect::Map(q_mvar.data(), q_mvar.size());

            // initialize properly the right "results" vectors (ie res_XXX RealVect)
            this->reset_results();
        }
        
        void init_osc_pq(const Eigen::Ref<const RealVect> & els_p,
                         const Eigen::Ref<const RealVect> & els_q,
                         const Eigen::Ref<const Eigen::VectorXi> & els_bus_id,
                         const std::string & name_el
                         )  // osc: one side element
        {
            init_osc(els_bus_id);
            int size = nb();
            check_size(els_p, size, name_el + "_p");
            check_size(els_q, size, name_el + "_q");

            target_p_mw_ = els_p;
            target_q_mvar_ = els_q;
        }

        void set_osc_pq_res_p(){
            const int nb_els = nb();
            for(int el_id = 0; el_id < nb_els; ++el_id){
                if(!status_[el_id]) res_p_[el_id] = 0.;
                else res_p_[el_id] = target_p_mw_(el_id);
            }
        }

        void set_osc_pq_res_q(bool ac){
            if(!ac){
                set_osc_res_q(ac);  // no q in DC mode
                return;
            }
            const int nb_els = nb();
            for(int el_id = 0; el_id < nb_els; ++el_id){
                if(!status_[el_id]) res_q_[el_id] = 0.;
                else res_q_[el_id] = target_q_mvar_(el_id);
            }
        }

        /**
         * Stamp `sign * (target_p + j.target_q)` of every active element into Sbus.
         * `sign` is +1 for an element in the generator convention (static
         * generators), -1 for one in the load convention (loads, storage units).
         * A leaf with a richer personality (a generator whose reactive output is
         * solved for, a shunt that only stamps in DC) writes its own `_fillSbus`.
         */
        void _stamp_pq(Eigen::Ref<CplxVect> Sbus,
                       const SolverBusIdVect & id_grid_to_solver,
                       real_type sign,
                       const char * fun_name) const
        {
            const int nb_els = nb();
            for(int el_id = 0; el_id < nb_els; ++el_id){
                if(!status_[el_id]) continue;
                const SolverBusId bus_id_solver = _solver_bus(el_id, bus_id_(el_id), id_grid_to_solver, fun_name);
                const cplx_type tmp = {target_p_mw_(el_id), target_q_mvar_(el_id)};
                Sbus.coeffRef(bus_id_solver.cast_int()) += sign * tmp;
            }
        }

    protected:
        // ---- the leaf hooks --------------------------------------------------------
        // an element that is in Sbus only: its setpoints are its results
        void _compute_res_pq(const Eigen::Ref<const RealVect> & /*Va*/,
                             const Eigen::Ref<const RealVect> & /*Vm*/,
                             const Eigen::Ref<const CplxVect> & /*V*/,
                             const SolverBusIdVect & /*id_grid_to_solver*/,
                             const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
                             real_type /*sn_mva*/,
                             bool ac) override
        {
            set_osc_pq_res_p();
            set_osc_pq_res_q(ac);
        }
        // ... and moving it, or switching it, only moves the injections
        void _on_deactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_sbus();
            solver_control.tell_one_el_changed_bus();
        }
        void _on_reactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_sbus();
            solver_control.tell_one_el_changed_bus();
        }
        void _on_change_bus(int /*el_id*/, GridModelBusId /*new_bus_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_sbus();
            solver_control.tell_one_el_changed_bus();
        }
        // Setpoint changes, notified BEFORE the write (both values are at hand).
        virtual void _on_change_p(int el_id, real_type new_p, DualAlgoControl & solver_control) {
            if (abs(target_p_mw_(el_id) - new_p) > _tol_equal_float) {
                solver_control.tell_recompute_sbus();
            }
        }
        virtual void _on_change_q(int el_id, real_type new_q, DualAlgoControl & solver_control) {
            if (abs(target_q_mvar_(el_id) - new_q) > _tol_equal_float) {
                solver_control.tell_recompute_sbus();
            }
        }

    protected:
        // physical properties

        // data for grid2op compat

        // input data
        RealVect target_p_mw_;
        RealVect target_q_mvar_;

};


} // namespace ls2g

#endif  //ONE_SIDE_CONTAINER_PQ_H
