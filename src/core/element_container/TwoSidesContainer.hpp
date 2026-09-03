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

// TODO other part of the API, like deactivate, reactivate etc.
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
                    if (static_cast<size_t>(my_id) >= r_data_two_sides.nb()) return;
                    id = my_id;

                    if(r_data_two_sides.names_.size()){
                        name = r_data_two_sides.names_[my_id];
                    }

                    connected_global = r_data_two_sides.status_global_[my_id];
                    
                    const auto & side_1_info = r_data_two_sides.side_1_[my_id];
                    const auto & side_2_info = r_data_two_sides.side_2_[my_id];
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
        size_t nb() const { return side_1_.nb(); }

        // Whole-grid semantic validation (see GenericContainer::check_valid):
        // each side is a full one-side container, so validate both. Derived
        // classes (eg TwoSidesContainer_rxh_A) add the branch electrical checks.
        void check_valid(int nb_bus,
                         int nb_sub,
                         const SubstationContainer & substations,
                         std::vector<int> & all_pos_topo_vect) const override
        {
            side_1_.check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);
            side_2_.check_valid(nb_bus, nb_sub, substations, all_pos_topo_vect);
        }

        GridModelBusId get_bus_side_1(int el_id) const {return side_1_.get_bus(el_id);}
        GridModelBusId get_bus_side_2(int el_id) const {return side_2_.get_bus(el_id);}
        // same, for an el_id one of our own loops produced (see _get_bus_internal)
        GridModelBusId get_bus_side_1_internal(int el_id) const {return side_1_.get_bus_internal(el_id);}
        GridModelBusId get_bus_side_2_internal(int el_id) const {return side_2_.get_bus_internal(el_id);}
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

        const GlobalBusIdVect & get_buses_side_1() const {return side_1_.get_buses();}
        const GlobalBusIdVect & get_buses_side_2() const {return side_2_.get_buses();}

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

        /**
         * An HVDC line has NO global gate: the two converter stations stand alone,
         * one can be on while the other is off, which is legitimate. So this simply
         * delegates -- deliberately unlike TwoSidesContainer_rxh_A above, which gates
         * on `status_global_` first.
         */
        void contribute_to_buses(int el_id, SubstationContainer & substation,
                                 int sign, bool & crossed) const override {
            side_1_.contribute_to_buses(el_id, substation, sign, crossed);
            side_2_.contribute_to_buses(el_id, substation, sign, crossed);
        }

        void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component, SubstationContainer & substation, DualAlgoControl & solver_control) override {
            const int nb_el = nb();
            const GlobalBusIdVect & bus_side_1_id_ = get_buses_side_1();
            const GlobalBusIdVect & bus_side_2_id_ = get_buses_side_2();
            for(int i = 0; i < nb_el; ++i){
                if(!status_global_[i]){
                    // the branch is globally off, so by contribute_to_buses it already
                    // holds nothing: the bracket takes nothing away and puts nothing
                    // back. Applied rather than assumed, and through the branch's rule
                    // -- a side's own deactivate() would use the one-sided one and
                    // decrement a bus the gate says the branch never held.
                    _apply_and_track_buses(i, substation, solver_control, [&]{
                        side_1_.deactivate_no_bus_tracking(i, solver_control);
                        side_2_.deactivate_no_bus_tracking(i, solver_control);
                    });
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
                    _apply_and_track_buses(i, substation, solver_control, [&]{
                        side_1_.deactivate_no_bus_tracking(i, solver_control);
                        side_2_.deactivate_no_bus_tracking(i, solver_control);
                        if(!ignore_status_global_) status_global_[i] = false;
                    });
                }
            }
        }
        void nb_line_end(std::vector<int> & res) const override final {
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

        virtual void update_topo(
            const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            const Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            DualAlgoControl & solver_control,
            SubstationContainer & substations
        ) final
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
            // used for updating derived class for example.
            this->_update_topo(solver_control, substations, side1_changed, side2_changed);

            for(int el_id=0; el_id<nb_el; ++el_id)
            {
                if(real_changed[el_id]) this->_update_effective_coeffs_one_el(el_id);
            }
        }

        // setter (states)
        // methods used within lightsim
        // The branch as a whole does the bus counting, ONCE, around whatever the two
        // sides do -- see OneSideContainer::deactivate_no_bus_tracking for why the
        // sides must not count for themselves.
        virtual void deactivate(int el_id, DualAlgoControl & solver_control,
                                SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "deactivate");  // before _apply_and_track_buses reads status_global_[el_id]
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_1_.deactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = side_2_.deactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = this->_deactivate(el_id, solver_control) || one_changed;
                _generic_deactivate(el_id, status_global_);
                if(ignore_status_global_) status_global_[el_id] = true;
            });
            if(one_changed){
                // update coefficient for Ybus
                this->_update_effective_coeffs_one_el(el_id);
            }
        }
        virtual void reactivate(int el_id, DualAlgoControl & solver_control,
                                SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "reactivate");  // before _apply_and_track_buses reads status_global_[el_id]
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                one_changed = side_1_.reactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = side_2_.reactivate_no_bus_tracking(el_id, solver_control) || one_changed;
                one_changed = this->_reactivate(el_id, solver_control) || one_changed;
                _generic_reactivate(el_id, status_global_);
                if(ignore_status_global_) status_global_[el_id] = true;
            });
            if(one_changed){
                this->_update_effective_coeffs_one_el(el_id);
            }
        }
        virtual void deactivate_side_1(int el_id, DualAlgoControl & solver_control,
                                       SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "deactivate_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                if(side_1_.deactivate_no_bus_tracking(el_id, solver_control)) this->_update_effective_coeffs_one_el(el_id);
            });
        }
        virtual void deactivate_side_2(int el_id, DualAlgoControl & solver_control,
                                       SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "deactivate_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                if(side_2_.deactivate_no_bus_tracking(el_id, solver_control)) this->_update_effective_coeffs_one_el(el_id);
            });
        }
        virtual void reactivate_side_1(int el_id, DualAlgoControl & solver_control,
                                       SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "reactivate_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                if(side_1_.reactivate_no_bus_tracking(el_id, solver_control)) this->_update_effective_coeffs_one_el(el_id);
            });
        }
        virtual void reactivate_side_2(int el_id, DualAlgoControl & solver_control,
                                       SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "reactivate_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
                if(side_2_.reactivate_no_bus_tracking(el_id, solver_control)) this->_update_effective_coeffs_one_el(el_id);
            });
        }

        void reset_results_tsc(){
            side_1_.reset_results();
            side_2_.reset_results();
        };

        /**
         * Change the bus on "side 1" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */        
        virtual void change_bus_side_1(int el_id, GridModelBusId new_gridmodel_bus_id, DualAlgoControl & solver_control, SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "change_bus_side_1");  // before _apply_and_track_buses reads status_global_[el_id]
            // the branch counts once, around both the side move and resolve_status --
            // a line END must not count for itself, status_global_ gates it
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 1).");
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
            one_changed = side_1_.change_bus_no_bus_tracking(el_id, new_gridmodel_bus_id, solver_control, substation);
            one_changed = resolve_status(el_id, true, solver_control) || one_changed;
            });
            this-> _change_bus_side_1(el_id, new_gridmodel_bus_id, solver_control, substation, one_changed);
            if(one_changed){
                // update coefficient for Ybus
                this->_update_effective_coeffs_one_el(el_id);
            }
        }
        /**
         * Change the bus on "side 2" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */  
        virtual void change_bus_side_2(int el_id, GridModelBusId new_gridmodel_bus_id, DualAlgoControl & solver_control, SubstationContainer & substation) final {
            _check_in_range(el_id, status_global_, "change_bus_side_2");  // before _apply_and_track_buses reads status_global_[el_id]
            // the branch counts once, around both the side move and resolve_status --
            // a line END must not count for itself, status_global_ gates it
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 2).");
            bool one_changed = false;
            _apply_and_track_buses(el_id, substation, solver_control, [&]{
            one_changed = side_2_.change_bus_no_bus_tracking(el_id, new_gridmodel_bus_id, solver_control, substation);
            one_changed = resolve_status(el_id, false, solver_control) || one_changed;
            });
            this-> _change_bus_side_2(el_id, new_gridmodel_bus_id, solver_control, substation, one_changed);
            if(one_changed){
                // update coefficient for Ybus
                this->_update_effective_coeffs_one_el(el_id);
            }
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
            ignore_status_global_ = std::get<0>(my_state);
            synch_status_both_side_ = std::get<1>(my_state);
            names_ = std::get<2>(my_state);
            status_global_ = std::get<3>(my_state);
            side_1_.set_state(std::get<4>(my_state));
            side_2_.set_state(std::get<5>(my_state));
            auto size = nb();
            if(names_.size() > 0) check_size(names_, size, "names");  // names are optional
            if(static_cast<size_t>(side_1_.nb()) != size) throw std::runtime_error("Side_1 do not have the proper size");
            if(static_cast<size_t>(side_2_.nb()) != size) throw std::runtime_error("Side_2 do not have the proper size");
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

        // hook when disconnecting or changing the bus a given element
        // for example used when disconnecting a powerline on only one side
        // to update yac_eff_12_, yac_eff_21_ etc. in TwoSidesContainer_rxh_A
        virtual void _update_effective_coeffs_one_el(int /*el_id*/) {
            // nothing to do by default
        }

        virtual bool _deactivate(int el_id, DualAlgoControl & /*solver_control*/) {
            // nothing to do by default: handled in derived class
            if(status_global_[el_id]) return true;
            return false;
        }
        virtual bool _reactivate(int el_id, DualAlgoControl & /*solver_control*/) {
            // nothing to do by default: handled in derived class
            if(!status_global_[el_id]) return true;
            return false;
        }

        virtual void _change_bus_side_1(
            int /*el_id*/, 
            GridModelBusId /*new_gridmodel_bus_id*/, 
            DualAlgoControl & /*solver_control*/, 
            const SubstationContainer & /*substation*/,
            bool /*has_effectively_changed*/
        ) {
            // nothing to do by default: handled in derived class
        }
        virtual void _change_bus_side_2(
            int /*el_id*/,
            GridModelBusId /*new_gridmodel_bus_id*/,
            DualAlgoControl & /*solver_control*/,
            const SubstationContainer & /*substation*/,
            bool /*has_effectively_changed*/
        ) {
            // nothing to do by default: handled in derived class
        }

        virtual void _update_topo(
            DualAlgoControl & /*solver_control*/,
            SubstationContainer & /*substations*/,
            const std::vector<bool> & /*side1_changed*/,
            const std::vector<bool> & /*side2_changed*/
        )
        {
            // nothing to do by default: handled in derived class
        }

        // DANGER ZONE, for modifiers
        // when it will be fully refactorize, it should disappear
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
