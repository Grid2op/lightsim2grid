// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TWO_SIDES_CONTAINER_H
#define TWO_SIDES_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "GenericContainer.hpp"

namespace ls2g {

/**
 * An element with two terminals: a line, a transformer (BranchContainer) or an
 * HVDC line (HvdcLineContainer). Each terminal is a full OneSideContainer of the
 * `OneSideType` given, and the element adds a global status on top.
 *
 * What lives here is everything that is about the pair and not about either
 * end: the per-bus element counts (the branch counts ONCE around whatever its
 * two ends do -- see OneSideContainer), the topology-vector update, the mutators
 * of one side or of the whole element, and the state.
 *
 * What a leaf writes:
 *   - `_on_connectivity_changed(el_id, ctrl)`: called after ANY change of the
 *     element's connectivity (either end's status or bus, the global status),
 *     once per element and per mutation, with the new state in place. A branch
 *     redoes its Kron-reduced coefficients there, a phase shifter raises the DC
 *     Sbus flag;
 *   - `_on_deactivate` / `_on_reactivate`: flags for the GLOBAL status flipping;
 *   - `_contribute_to_buses`: whether the global status gates the ends (a line)
 *     or not (an HVDC line, whose stations stand alone);
 *   - `_fill*`, `_compute_results`, `_get_graph`, ...: per GenericContainer.
 */
template<class OneSideType>
class TwoSidesContainer : public GenericContainer
{
    public:
        class TwoSidesInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                std::string name;
                int sub_1_id;
                int sub_2_id;
                int pos_1_topo_vect;
                int pos_2_topo_vect;

                bool connected_global;
                bool connected_1;
                bool connected_2;

                int bus_1_id;
                int bus_2_id;

                bool has_res;
                real_type res_p1_mw;
                real_type res_q1_mvar;
                real_type res_v1_kv;
                real_type res_theta1_deg;
                real_type res_p2_mw;
                real_type res_q2_mvar;
                real_type res_v2_kv;
                real_type res_theta2_deg;

                TwoSidesInfo(const TwoSidesContainer & r_data_two_sides, int my_id) noexcept:
                id(-1),
                name(""),
                sub_1_id(-1),
                sub_2_id(-1),
                pos_1_topo_vect(-1),
                pos_2_topo_vect(-1),
                connected_global(false),
                connected_1(false),
                connected_2(false),
                bus_1_id(_deactivated_bus_id),
                bus_2_id(_deactivated_bus_id),
                has_res(false),
                res_p1_mw(0.),
                res_q1_mvar(0.),
                res_v1_kv(0.),
                res_theta1_deg(0.),
                res_p2_mw(0.),
                res_q2_mvar(0.),
                res_v2_kv(0.),
                res_theta2_deg(0.)
                {
                    if (my_id < 0) return;
                    if (my_id >= r_data_two_sides.nb()) return;
                    id = my_id;

                    if(r_data_two_sides.names_.size()){
                        name = r_data_two_sides.names_[my_id];
                    }

                    connected_global = r_data_two_sides.status_global_[my_id];
                    
                    // (a side's Info, built directly: the side is not iterable by itself)
                    const typename OneSideType::DataInfo side_1_info(r_data_two_sides.side_1_, my_id);
                    const typename OneSideType::DataInfo side_2_info(r_data_two_sides.side_2_, my_id);
                    sub_1_id = side_1_info.sub_id;
                    sub_2_id = side_2_info.sub_id;
                    pos_1_topo_vect = side_1_info.pos_topo_vect;
                    pos_2_topo_vect = side_2_info.pos_topo_vect;
                    connected_1 = side_1_info.connected;
                    connected_2 = side_2_info.connected;
                    bus_1_id = side_1_info.bus_id;
                    bus_2_id = side_2_info.bus_id;

                    if(side_1_info.has_res)
                    {
                        res_p1_mw = side_1_info.res_p_mw;
                        res_q1_mvar = side_1_info.res_q_mvar;
                        res_v1_kv = side_1_info.res_v_kv;
                        res_theta1_deg = side_1_info.res_theta_deg;
                    }
                    if(side_2_info.has_res)
                    {
                        res_p2_mw = side_2_info.res_p_mw;
                        res_q2_mvar = side_2_info.res_q_mvar;
                        res_v2_kv = side_2_info.res_v_kv;
                        res_theta2_deg = side_2_info.res_theta_deg;
                    }
                }
        };

    public:
        TwoSidesContainer() noexcept :ignore_status_global_(false), synch_status_both_side_(true){}
        ~TwoSidesContainer() noexcept override = default;

        // public generic API
        int nb() const { return side_1_.nb(); }

    protected:
        // Whole-grid semantic validation (see GenericContainer::check_valid):
        // each side is a full one-side container, so validate both. Derived
        // classes (eg BranchContainer) add the branch electrical checks.
        void _check_valid(int nb_bus,
                          int nb_sub,
                          const SubstationContainer & substations,
                          std::vector<int> & all_pos_topo_vect) const override
        {
            side_1_.check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);
            side_2_.check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);
        }

        void _compute_results(const Eigen::Ref<const RealVect> & Va,
                              const Eigen::Ref<const RealVect> & Vm,
                              const Eigen::Ref<const CplxVect> & V,
                              const SolverBusIdVect & id_grid_to_solver,
                              const Eigen::Ref<const RealVect> & bus_vn_kv,
                              real_type sn_mva,
                              bool ac) override
        {
            side_1_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
            side_2_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        }
        void _reset_results() override {
            side_1_.reset_results();
            side_2_.reset_results();
        }

    public:

        GridModelBusId get_bus_side_1(int el_id) const {return side_1_.get_bus(el_id);}
        GridModelBusId get_bus_side_2(int el_id) const {return side_2_.get_bus(el_id);}
        // Per-side connectivity: an element can be `connected_global` (the
        // TwoSidesContainer-level status, e.g. an HVDC line kept active
        // because at least one converter is in the main synchronous
        // component -- see disconnect_if_not_in_main_component) while ONE
        // side is individually open (real RTE grids: a half-open HVDC line
        // with its remote converter in another synchronous island). Callers
        // that pull per-side data (e.g. droop flows, Q-limit masking) MUST
        // check these, not just nb()/connected_global, or they will silently
        // treat an open side as a normal, both-ends-connected element.
        bool get_connected_side_1(int el_id) const {return side_1_.get_status(el_id);}
        bool get_connected_side_2(int el_id) const {return side_2_.get_status(el_id);}

        void init_tsc(
            const Eigen::Ref<const Eigen::VectorXi> & els_bus1_id,
            const Eigen::Ref<const Eigen::VectorXi> & els_bus2_id,
            const std::string & name_elements
        )  // tsc: two sides container
        {
            auto size = els_bus1_id.size();
            check_size(els_bus2_id, size, name_elements);
            side_1_.init_osc(els_bus1_id);
            side_2_.init_osc(els_bus2_id);
            status_global_ = std::vector<bool>(els_bus1_id.size(), true);
        }

        tuple3d get_res_side_1() const {return side_1_.get_res();}
        tuple3d get_res_side_2() const {return side_2_.get_res();}

        tuple4d get_res_full_side_1() const {return side_1_.get_res_full();}
        tuple4d get_res_full_side_2() const {return side_2_.get_res_full();}

        Eigen::Ref<const RealVect> get_theta_side_1() const {return side_1_.get_theta();}
        Eigen::Ref<const RealVect> get_theta_side_2() const {return side_2_.get_theta();}

        const std::vector<bool>& get_status_global() const {return status_global_;}
        const std::vector<bool>& get_status_side_1() const {return side_1_.get_status();}
        const std::vector<bool>& get_status_side_2() const {return side_2_.get_status();}

        const GlobalBusIdVect & get_bus_id_side_1() const {return side_1_.get_bus_id();}
        const GlobalBusIdVect & get_bus_id_side_2() const {return side_2_.get_bus_id();}

        Eigen::Ref<const IntVect> get_bus_id_side_1_numpy() const {return side_1_.get_bus_id_numpy();}
        Eigen::Ref<const IntVect> get_bus_id_side_2_numpy() const {return side_2_.get_bus_id_numpy();}

    protected:
        // same, for an el_id one of our own loops produced (see _get_bus_internal)
        GridModelBusId get_bus_side_1_internal(int el_id) const {return side_1_.get_bus_internal(el_id);}
        GridModelBusId get_bus_side_2_internal(int el_id) const {return side_2_.get_bus_internal(el_id);}

        /**
         * The default has NO global gate: the two ends stand alone (an HVDC line,
         * whose converter stations can legitimately be one on, one off). A line or a
         * transformer overrides this to gate on `status_global_` first: see
         * BranchContainer.
         */
        void _contribute_to_buses(int el_id, SubstationContainer & substation,
                                  int sign, bool & crossed) const override {
            side_1_.contribute_to_buses(el_id, substation, sign, crossed);
            side_2_.contribute_to_buses(el_id, substation, sign, crossed);
        }

        /**
         * Open the ends of element `el_id` that are asked for (`side_1` / `side_2`),
         * and the whole element too (`whole`, honouring `ignore_status_global_`),
         * without going through the public mutators. ONE per-bus count bracket
         * around all of it, through the element's own contribution rule, which
         * is what the main-component clean-up needs: a side's own `deactivate()`
         * would count with the one-sided rule and decrement a bus the global gate
         * says the branch never held. Returns whether anything changed.
         */
        bool _open_sides(int el_id, bool side_1, bool side_2, bool whole,
                         DualAlgoControl & solver_control, SubstationContainer & substation)
        {
            bool changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                if(side_1) changed = side_1_.deactivate_no_bus_tracking(el_id, solver_control) || changed;
                if(side_2) changed = side_2_.deactivate_no_bus_tracking(el_id, solver_control) || changed;
                if(whole && !ignore_status_global_ && status_global_[el_id]){
                    status_global_[el_id] = false;
                    changed = true;
                }
            });
            if(changed) _on_connectivity_changed(el_id, solver_control);
            return changed;
        }

        void _disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component, SubstationContainer & substation, DualAlgoControl & solver_control) override {
            const int nb_el = nb();
            const GlobalBusIdVect & bus_side_1_id_ = get_bus_id_side_1();
            const GlobalBusIdVect & bus_side_2_id_ = get_bus_id_side_2();
            for(int i = 0; i < nb_el; ++i){
                if(!status_global_[i]){
                    // the branch is globally off, so by contribute_to_buses it already
                    // holds nothing: the bracket takes nothing away and puts nothing
                    // back. Applied rather than assumed, and through the branch's rule
                    // -- a side's own deactivate() would use the one-sided one and
                    // decrement a bus the gate says the branch never held.
                    _open_sides(i, true, true, false, solver_control, substation);
                    continue;
                }
                // A side is "outside the main component" only if it is CONNECTED
                // (its bus is a real bus, not the deactivated/open marker) AND that
                // bus is not flagged in the main component. An open side (bus ==
                // _deactivated_bus_id, e.g. a half-open line) imposes no constraint:
                // such a branch stays as long as its connected side(s) are in main.
                const int b1 = bus_side_1_id_(i).cast_int();
                const int b2 = bus_side_2_id_(i).cast_int();
                const bool s1_outside = (b1 != _deactivated_bus_id) && !busbar_in_main_component[b1];
                const bool s2_outside = (b2 != _deactivated_bus_id) && !busbar_in_main_component[b2];
                if(s1_outside || s2_outside){
                    // island, boundary, or (defensively) a branch straddling two
                    // components: drop the whole element rather than throw. Keeping
                    // the main component well-posed is the goal of this function.
                    // One bracket around both sides AND the `status_global_` flip, as
                    // everywhere else: the gate is what contribute_to_buses reads first.
                    _open_sides(i, true, true, true, solver_control, substation);
                }
            }
        }
        void _nb_line_end(std::vector<int> & res) const override final {
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;

                const GlobalBusId bus_or = get_bus_side_1_internal(el_id);
                if(bus_or.cast_int() != _deactivated_bus_id) res[bus_or.cast_int()] += 1;
                const GlobalBusId bus_ex = get_bus_side_2_internal(el_id);
                if(bus_ex.cast_int() != _deactivated_bus_id) res[bus_ex.cast_int()] += 1;
            }
        }

    public:
        void set_pos_topo_vect_side_1(const Eigen::Ref<const IntVect> & pos_topo_vect)
        {
            side_1_.set_pos_topo_vect(pos_topo_vect);
        }
        void set_pos_topo_vect_side_2(const Eigen::Ref<const IntVect> & pos_topo_vect)
        {
            side_2_.set_pos_topo_vect(pos_topo_vect);
        }

        void set_subid_side_1(const Eigen::Ref<const IntVect> & subid)
        {
            side_1_.set_subid(subid);
        }
        void set_subid_side_2(const Eigen::Ref<const IntVect> & subid)
        {
            side_2_.set_subid(subid);
        }

    protected:
        std::vector<bool> _update_topo(
            const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            const Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            DualAlgoControl & solver_control,
            SubstationContainer & substations
        ) override final
        {
            side_1_._check_pos_topo_vect_filled();
            side_2_._check_pos_topo_vect_filled();

            const int nb_el = nb();
            const int nb_topo = static_cast<int>(has_changed.rows());
            std::vector<bool> side1_changed(nb_el, false);
            std::vector<bool> side2_changed(nb_el, false);
            std::vector<bool> real_changed(nb_el, false);

            for(int el_id=0; el_id<nb_el; ++el_id)
            {
                const int pos1 = side_1_.checked_pos_topo_vect(el_id, nb_topo);
                const int pos2 = side_2_.checked_pos_topo_vect(el_id, nb_topo);
                const bool touched_1 = has_changed(pos1);
                const bool touched_2 = has_changed(pos2);
                // an element the caller did not touch mutates nothing, so it must not
                // be bracketed either -- see OneSideContainer::update_topo.
                if(!touched_1 && !touched_2) continue;

                // The branch counts ONCE, around everything this element's entry does:
                // both side moves AND `resolve_status`, which sets `status_global_` --
                // the gate contribute_to_buses reads FIRST. Leaving resolve_status
                // outside the bracket is how a bus whose only element was the ex side
                // of a line stayed counted after the line was disconnected: the system
                // lost a bus, and nothing told the solver its dimension had changed.
                // The sides use the *_no_bus_tracking entry point for the same reason
                // deactivate() does: a line end does not own its contribution.
                bool real_change = false;
                _apply_and_track_buses(el_id, substations, solver_control, [&]{
                    const bool s1 = side_1_.update_topo_one_el_no_bus_tracking(el_id, has_changed, new_values, solver_control, substations);
                    const bool s2 = side_2_.update_topo_one_el_no_bus_tracking(el_id, has_changed, new_values, solver_control, substations);
                    side1_changed[el_id] = s1;
                    side2_changed[el_id] = s2;
                    real_change = s1 || s2;
                    // set the global status
                    if(touched_1){
                        real_change = resolve_status(el_id, true, solver_control) || real_change;
                    }
                    if(touched_2){
                        real_change = resolve_status(el_id, false, solver_control) || real_change;
                    }
                });
                real_changed[el_id] = real_change;
            }
            for(int el_id=0; el_id<nb_el; ++el_id)
            {
                if(real_changed[el_id]) _on_connectivity_changed(el_id, solver_control);
            }
            return real_changed;
        }

    public:
        // setter (states)
        // methods used within lightsim
        // The branch as a whole does the bus counting, ONCE, around whatever the two
        // sides do -- see OneSideContainer::deactivate_no_bus_tracking for why the
        // sides must not count for themselves.
        void deactivate(int el_id, DualAlgoControl & solver_control,
                        SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "deactivate");  // before _apply_and_track_buses reads status_global_[el_id]
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_1_.deactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = side_2_.deactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                if(status_global_[el_id]){
                    _on_deactivate(el_id, solver_control);
                    one_changed = true;
                }
                status_global_[el_id] = ignore_status_global_;  // off, unless the gate is ignored
            });
            if(one_changed) _on_connectivity_changed(el_id, solver_control);
        }
        void reactivate(int el_id, DualAlgoControl & solver_control,
                        SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "reactivate");  // before _apply_and_track_buses reads status_global_[el_id]
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_1_.reactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = side_2_.reactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                if(!status_global_[el_id]){
                    _on_reactivate(el_id, solver_control);
                    one_changed = true;
                }
                status_global_[el_id] = true;
            });
            if(one_changed) _on_connectivity_changed(el_id, solver_control);
        }
        void deactivate_side_1(int el_id, DualAlgoControl & solver_control,
                               SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "deactivate_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            _open_sides(el_id, true, false, false, solver_control, substation);
        }
        void deactivate_side_2(int el_id, DualAlgoControl & solver_control,
                               SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "deactivate_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            _open_sides(el_id, false, true, false, solver_control, substation);
        }
        void reactivate_side_1(int el_id, DualAlgoControl & solver_control,
                               SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "reactivate_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            bool changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                changed = side_1_.reactivate_no_bus_tracking(el_id, solver_control);
            });
            if(changed) _on_connectivity_changed(el_id, solver_control);
        }
        void reactivate_side_2(int el_id, DualAlgoControl & solver_control,
                               SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "reactivate_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            bool changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                changed = side_2_.reactivate_no_bus_tracking(el_id, solver_control);
            });
            if(changed) _on_connectivity_changed(el_id, solver_control);
        }

        /**
         * Change the bus on "side 1" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */        
        void change_bus_side_1(int el_id, GridModelBusId new_gridmodel_bus_id, DualAlgoControl & solver_control, SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "change_bus_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            // and the BUS id, before the bracket takes this element's contribution
            // away -- a call the grid will refuse must not touch the counts at all.
            _check_new_bus_id(new_gridmodel_bus_id, substation.nb_bus());
            // the branch counts once, around both the side move and resolve_status --
            // a line END must not count for itself, status_global_ gates it
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 1).");
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_1_.change_bus_no_bus_tracking(el_id, new_gridmodel_bus_id, solver_control);
                one_changed = resolve_status(el_id, true, solver_control) || one_changed;
            });
            if(one_changed) _on_connectivity_changed(el_id, solver_control);
        }
        /**
         * Change the bus on "side 2" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */  
        void change_bus_side_2(int el_id, GridModelBusId new_gridmodel_bus_id, DualAlgoControl & solver_control, SubstationContainer & substation) {
            _check_in_range(el_id, status_global_, "change_bus_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            // and the BUS id, before the bracket takes this element's contribution
            // away -- a call the grid will refuse must not touch the counts at all.
            _check_new_bus_id(new_gridmodel_bus_id, substation.nb_bus());
            // the branch counts once, around both the side move and resolve_status --
            // a line END must not count for itself, status_global_ gates it
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 2).");
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_2_.change_bus_no_bus_tracking(el_id, new_gridmodel_bus_id, solver_control);
                one_changed = resolve_status(el_id, false, solver_control) || one_changed;
            });
            if(one_changed) _on_connectivity_changed(el_id, solver_control);
        }

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            bool,  // ignore_status_global_
            bool,  // synch_status_both_side_
            std::vector<std::string>,
            std::vector<bool>,          // status_global
            typename OneSideType::StateRes, // side_1
            typename OneSideType::StateRes  // side_2
            >;
        enum StateResIdx {
            IGNORE_STATUS_GLOBAL = 0,
            SYNCH_STATUS_BOTH_SIDE,
            NAMES,
            STATUS_GLOBAL,
            SIDE_1,
            SIDE_2,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "TwoSidesContainer::StateRes and StateResIdx do not match");

        void set_ignore_status_global(bool ignore_status_global){
            ignore_status_global_ = ignore_status_global;
        }
        bool get_ignore_status_global() const{
            return ignore_status_global_;
        }
        void set_synch_status_both_side(bool synch_status_both_side){
            synch_status_both_side_=synch_status_both_side;
        }
        bool get_synch_status_both_side() const{
            return synch_status_both_side_;
        }

    protected:
        bool ignore_status_global_;
        bool synch_status_both_side_;
        
        OneSideType side_1_;
        OneSideType side_2_;

        std::vector<bool> status_global_;

    protected:
        StateRes get_tsc_state() const  // tsc: two sides container
        {
            StateRes res(
                ignore_status_global_,
                synch_status_both_side_,
                names_,
                status_global_,
                side_1_.get_state(),
                side_2_.get_state()
            );
            return res;
        }

        void set_tsc_state(TwoSidesContainer::StateRes & my_state)  // tsc: two sides container
        {
            ignore_status_global_ = std::get<StateResIdx::IGNORE_STATUS_GLOBAL>(my_state);
            synch_status_both_side_ = std::get<StateResIdx::SYNCH_STATUS_BOTH_SIDE>(my_state);
            names_ = std::get<StateResIdx::NAMES>(my_state);
            status_global_ = std::get<StateResIdx::STATUS_GLOBAL>(my_state);
            side_1_.set_state(std::get<StateResIdx::SIDE_1>(my_state));
            side_2_.set_state(std::get<StateResIdx::SIDE_2>(my_state));
            const int size = nb();
            if(names_.size() > 0) check_size(names_, size, "names");  // names are optional
            if(side_1_.nb() != size) throw std::runtime_error("Side_1 do not have the proper size");
            if(side_2_.nb() != size) throw std::runtime_error("Side_2 do not have the proper size");
            // `nb()` is side_1_.nb(), NOT status_global_.size(): nothing above ties the
            // two together, yet status_global_ is indexed with element ids bounded by
            // nb() all over this class (resolve_status, _deactivate, fillYbus, the batch
            // solvers...) with an unchecked operator[]. A pickle / binary file declaring
            // a shorter (in particular empty) status_global_ therefore reads and writes
            // past its end -- and check_grid() never sees it, it runs later and only
            // looks at the per-side data. Demand the exact length here.
            check_size(status_global_, size, "status_global");
        }

        bool resolve_status(int el_id, bool side_1_modif, DualAlgoControl & solver_control){
            OneSideType & side_modified = side_1_modif ? side_1_: side_2_;
            OneSideType & side_to_update = side_1_modif ? side_2_: side_1_;
            bool res = false;
            if(synch_status_both_side_){
                if(side_modified.get_status(el_id)){
                    // element has been reconnected
                    // I need to reconnect other side
                    res = res || side_to_update.reactivate_no_bus_tracking(el_id, solver_control);
                    status_global_[el_id] = true;
                    res = true;
                }else{
                    res = res || side_to_update.deactivate_no_bus_tracking(el_id, solver_control);
                    status_global_[el_id] = false;
                }
            }
            if(ignore_status_global_) status_global_[el_id] = true;  // always true in this case
            else{
                if(side_modified.get_status(el_id) == side_to_update.get_status(el_id)){
                    res = res || (status_global_[el_id] != side_modified.get_status(el_id));
                    status_global_[el_id] = side_modified.get_status(el_id);
                }
            }
            return res;
        }

        // ---- the leaf hooks --------------------------------------------------------
        /**
         * The connectivity of element `el_id` -- either end's status or bus, or the
         * global status -- has just changed, and the new state is in place. Called
         * once per element and per mutation, from every mutator above (whole
         * element, one side, topology vector, main-component clean-up). Derived
         * state that depends on which ends are connected is refreshed here (the
         * Kron-reduced coefficients of a branch), and so are the flags a leaf owes
         * for that change on top of what the ends already raised (the DC Sbus term
         * of a phase-shifting transformer).
         */
        virtual void _on_connectivity_changed(int /*el_id*/, DualAlgoControl & /*solver_control*/) {}
        /// the GLOBAL status of `el_id` is about to flip (the ends have their own hooks)
        virtual void _on_deactivate(int /*el_id*/, DualAlgoControl & /*solver_control*/) {}
        virtual void _on_reactivate(int /*el_id*/, DualAlgoControl & /*solver_control*/) {}

        // writable results, for a leaf that computes both ends' flows itself
        Eigen::Ref<RealVect> get_res_theta_side_1() {return side_1_.get_res_theta();}
        Eigen::Ref<RealVect> get_res_p_side_1() {return side_1_.get_res_p();}
        Eigen::Ref<RealVect> get_res_q_side_1() {return side_1_.get_res_q();}
        Eigen::Ref<RealVect> get_res_v_side_1() {return side_1_.get_res_v();}
        Eigen::Ref<RealVect> get_res_theta_side_2() {return side_2_.get_res_theta();}
        Eigen::Ref<RealVect> get_res_p_side_2() {return side_2_.get_res_p();}
        Eigen::Ref<RealVect> get_res_q_side_2() {return side_2_.get_res_q();}
        Eigen::Ref<RealVect> get_res_v_side_2() {return side_2_.get_res_v();}

};



} // namespace ls2g

#endif  // TWO_SIDES_CONTAINER_H
