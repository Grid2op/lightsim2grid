// Copyright (c) 2024-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ONE_SIDE_CONTAINER_H
#define ONE_SIDE_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "GenericContainer.hpp"
#include "SubstationContainer.hpp"

namespace ls2g {

// same for all
// - X nb 
// - X get_bus
// - get_buses
// - get_res
// - get_res_full
// - get_theta
// - get_status
// - get_bus_id
// - gen_p_per_bus

// same public api but need overriden in private api
// - deactivate
// - reactivate
// - change_bus
// - change_p
// - change_q
// - reset_results
// - compute_results

// need to modify in overriden class
// - get_state
// - set_state
// - init

template<class OneSideType>
class TwoSidesContainer;

/**
 * This is the most generic part of the "one side container".
 * 
 * It can be used to represent side of "multi sided elements" 
 * (such as Lines or Transformers) or element directly connected
 * to one bus (such as Loads or Generators).
 *
 * It is not meant to be used directly.
 */
class OneSideContainer : public GenericContainer
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    // provide access to all instanciation of "TwoSidesContainer" class 
    // to protected members of "OneSideContainer" (eg set_osc_state)
    template<class T>
    friend class TwoSidesContainer;

    // regular implementation
    public:

        class OneSideInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                std::string name;
                int sub_id;
                int pos_topo_vect;

                bool connected;
                int bus_id;

                bool has_res;
                real_type res_p_mw;
                real_type res_q_mvar;
                real_type res_v_kv;
                real_type res_theta_deg;

                OneSideInfo(const OneSideContainer & r_data_one_side, int my_id) noexcept:
                id(-1),
                name(""),
                sub_id(-1),
                pos_topo_vect(-1),
                connected(false),
                bus_id(_deactivated_bus_id),
                has_res(false),
                res_p_mw(0.),
                res_q_mvar(0.),
                res_v_kv(0.),
                res_theta_deg(0.)
                {
                    if((my_id >= 0) && (my_id < r_data_one_side.nb()))
                    {
                        id = my_id;
                        if(r_data_one_side.names_.size()){
                            name = r_data_one_side.names_[my_id];
                        }
                        if(r_data_one_side.subid_.size()){
                            sub_id = r_data_one_side.subid_(my_id);
                        }
                        if(r_data_one_side.pos_topo_vect_.size()){
                            pos_topo_vect = r_data_one_side.pos_topo_vect_(my_id);
                        }
                        connected = r_data_one_side.status_[my_id];
                        if(connected) bus_id = r_data_one_side.bus_id_[my_id].cast_int();

                        has_res = r_data_one_side.res_p_.size() > 0;
                        if(has_res)
                        {
                            res_p_mw = r_data_one_side.res_p_.coeff(my_id);
                            res_q_mvar = r_data_one_side.res_q_.coeff(my_id);
                            res_v_kv = r_data_one_side.res_v_.coeff(my_id);
                            res_theta_deg = r_data_one_side.res_theta_.coeff(my_id);
                        }
                    }
                }
        };
        using DataInfo = OneSideInfo;

    /////////////////////////////////////
    // iterator
    // private:
    //     typedef GenericContainerConstIterator<OneSideContainer> OSCConstIterator;

    // public:
    //     OSCConstIterator begin() const {return OSCConstIterator(this, 0); }
    //     OSCConstIterator end() const {return OSCConstIterator(this, nb()); }
    //     OneSideInfo operator[](int id) const
    //     {
    //         if(id < 0)
    //         {
    //             throw std::range_error("You cannot ask for a negative load id.");
    //         }
    //         if(id >= nb())
    //         {
    //             throw std::range_error("Load out of bound. Not enough loads on the grid.");
    //         }
    //         return OneSideInfo(*this, id);
    //     }
    /////////////////////////////////////

    public:
        OneSideContainer() noexcept = default;
        ~OneSideContainer() noexcept override = default;
        // OneSideInfo get_osc_info(int id_) {return OneSideInfo(*this, id_);}

        // public generic API
        int nb() const { return static_cast<int>(bus_id_.size()); }
        GridModelBusId get_bus(int el_id) const {return _get_bus(el_id, status_, bus_id_);}
        const GlobalBusIdVect & get_buses() const {return bus_id_;}

        tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
        tuple4d get_res_full() const {return tuple4d(res_p_, res_q_, res_v_, res_theta_);}
        
        Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
        const std::vector<bool>& get_status() const {return status_;}
        bool get_status(int el_id) const {return status_.at(el_id);}
        const GlobalBusIdVect & get_bus_id() const {return bus_id_;}
        Eigen::Ref<const IntVect> get_bus_id_numpy() const {
            return bus_id_.as_eigen();
        }

        /// one-sided: this element holds its own bus, and only while it is active
        void contribute_to_buses(int el_id, SubstationContainer & substation,
                                 int sign, bool & crossed) const override {
            if(!status_[el_id]) return;                 // inactive: holds nothing
            const GlobalBusId my_bus = bus_id_(el_id);
            if(my_bus.cast_int() == _deactivated_bus_id) return;
            crossed |= (sign > 0) ? substation.bus_gained_element(my_bus)
                                  : substation.bus_lost_element(my_bus);
        }

        void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component, SubstationContainer & substation, DualAlgoControl & solver_control) override final {
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                if(!status_[el_id]) continue;
                const GlobalBusId my_bus = bus_id_(el_id);
                if(!busbar_in_main_component[my_bus.cast_int()]){
                    deactivate(el_id, solver_control, substation);
                }
            }    
        }

        virtual bool deactivate(int el_id, DualAlgoControl & solver_control,
                                SubstationContainer & substation) final {
            // validate BEFORE _apply_and_track_buses: it asks contribute_to_buses for
            // the element's current contribution first, which indexes status_[el_id]
            // with an unchecked operator[] (a negative id wraps to a huge size_t).
            // The check inside *_no_bus_tracking is too late to stop that.
            _check_in_range(el_id, status_, "deactivate");
            bool res = false;
            _apply_and_track_buses(el_id, substation, solver_control,
                                   [&]{ res = deactivate_no_bus_tracking(el_id, solver_control); });
            return res;
        }
        virtual bool reactivate(int el_id, DualAlgoControl & solver_control,
                                SubstationContainer & substation) final {
            // validate BEFORE _apply_and_track_buses: it asks contribute_to_buses for
            // the element's current contribution first, which indexes status_[el_id]
            // with an unchecked operator[] (a negative id wraps to a huge size_t).
            // The check inside *_no_bus_tracking is too late to stop that.
            _check_in_range(el_id, status_, "reactivate");
            bool res = false;
            _apply_and_track_buses(el_id, substation, solver_control,
                                   [&]{ res = reactivate_no_bus_tracking(el_id, solver_control); });
            return res;
        }

        /**
         * The same mutation WITHOUT touching the per-bus element counts.
         *
         * For a container that owns its own contribution (a load, a generator, an
         * HVDC converter station) this is only ever called through `deactivate`
         * above, which brackets it with the counting. It is public because
         * TwoSidesContainer must call it directly: a line's two ends do NOT own
         * their contribution -- `status_global_` gates it, and a side knows nothing
         * about that -- so the branch as a whole does the counting, once, around
         * both sides. Letting each side count here would decrement a bus the gate
         * says the branch never held.
         */
        bool deactivate_no_bus_tracking(int el_id, DualAlgoControl & solver_control) {
            // validate el_id *before* dispatching: `_deactivate` indexes status_[el_id]
            // with an unchecked operator[] (a negative id would wrap to a huge size_t),
            // and `_generic_deactivate` only checks afterwards.
            _check_in_range(el_id, status_, "deactivate");
            bool res = this->_deactivate(el_id, solver_control);
            _generic_deactivate(el_id, status_);
            return res;
        }
        /// change_bus WITHOUT touching the per-bus counts; see
        /// deactivate_no_bus_tracking for why TwoSidesContainer needs this.
        bool change_bus_no_bus_tracking(int el_id, GridModelBusId new_gridmodel_bus_id,
                                        DualAlgoControl & solver_control,
                                        const SubstationContainer & substation) {
            _check_in_range(el_id, bus_id_, "change_bus");
            if(bus_id_(el_id) == new_gridmodel_bus_id) return false;
            bool res = this->_change_bus(el_id, new_gridmodel_bus_id, solver_control, substation.nb_bus());
            _generic_change_bus(el_id, new_gridmodel_bus_id, bus_id_, solver_control, substation.nb_bus());
            return res;
        }
        bool reactivate_no_bus_tracking(int el_id, DualAlgoControl & solver_control) {
            _check_in_range(el_id, status_, "reactivate");
            bool res = this->_reactivate(el_id, solver_control);
            _generic_reactivate(el_id, status_);
            return res;
        }

        /**
         * This function changes the bus. The bus_id is here given in the
         * "gridmodel" bus.
         * 
         * Not the "solver" bus, nor the "substation" / "local" bus.
         */
        virtual bool change_bus(
            int load_id,
            GridModelBusId new_gridmodel_bus_id,
            DualAlgoControl & solver_control,
            SubstationContainer & substation) final {
                // validate load_id *before* dispatching: `_change_bus` reads bus_id_(load_id)
                // with an unchecked Eigen operator(); `_generic_change_bus` only checks afterwards.
                _check_in_range(load_id, bus_id_, "change_bus");
                // a move to the bus it is already on holds exactly the same bus
                // afterwards. Tracking it would take the contribution away and put
                // it straight back -- correct counts, but a bus that is alone would
                // transiently hit 0 and report a crossing that never happened,
                // costing a full rebuild. grid2op sends this every step.
                if(bus_id_(load_id) == new_gridmodel_bus_id) return false;
                bool res = false;
                _apply_and_track_buses(load_id, substation, solver_control, [&]{
                    res = change_bus_no_bus_tracking(load_id, new_gridmodel_bus_id, solver_control, substation);
                });
                return res;
        }

        virtual void compute_results(const Eigen::Ref<const RealVect> & Va,
                                     const Eigen::Ref<const RealVect> & Vm,
                                     const Eigen::Ref<const CplxVect> & V,
                                     const SolverBusIdVect & id_grid_to_solver,
                                     const Eigen::Ref<const RealVect> & bus_vn_kv,
                                     real_type sn_mva,
                                     bool ac) final
        {
            const int nb_els = nb();
            v_kv_from_vpu(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
            v_deg_from_va(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
            this->_compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        }

        virtual void reset_results() final {
            reset_osc_results();
        }

        void set_pos_topo_vect(const Eigen::Ref<const IntVect> & pos_topo_vect)
        {
            check_size(pos_topo_vect, nb(), "pos_topo_vect");
            // The upper bound (dim_topo) is the sum of every topology-participating
            // container's element count, not just this one's -- it is not known here
            // (see update_topo(), which validates the full [0, dim_topo) range with
            // that context). A negative position, though, is never valid regardless of
            // dim_topo, so reject it immediately rather than only at the point of use.
            for(Eigen::Index i = 0; i < pos_topo_vect.size(); ++i){
                if(pos_topo_vect(i) < 0){
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::set_pos_topo_vect: element " << i
                         << " has a negative position in the topology vector ("
                         << pos_topo_vect(i) << ").";
                    throw std::out_of_range(exc_.str());
                }
            }
            pos_topo_vect_.array() = pos_topo_vect;
        }

        void set_subid(const Eigen::Ref<const IntVect> & subid)
        {
            check_size(subid, nb(), "subid");
            // The upper bound (n_sub) is owned by SubstationContainer, not available
            // here (see update_topo(), which validates the full [0, n_sub) range with
            // that context). A negative substation id, though, is never valid
            // regardless of n_sub, so reject it immediately.
            for(Eigen::Index i = 0; i < subid.size(); ++i){
                if(subid(i) < 0){
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::set_subid: element " << i
                         << " has a negative substation id (" << subid(i) << ").";
                    throw std::out_of_range(exc_.str());
                }
            }
            subid_.array() = subid;
        }

        /**
         * The position this element occupies in the grid2op topology vector, checked
         * against the length of the arrays that position is about to index.
         *
         * `pos_topo_vect_` indexes has_changed / new_values with an unchecked Eigen
         * operator(). check_grid() proves the stored positions form a permutation of
         * [0, dim_topo), but the set_pos_topo_vect() setters only check the vector's
         * *length*, not its values -- so a position written straight through a setter
         * (bypassing check_grid) would read past the caller arrays. Release wheels are
         * -O3 -DNDEBUG, so neither Eigen nor the STL catches it: validate here, which
         * every update_topo path goes through before indexing.
         */
        int checked_pos_topo_vect(int el_id, int nb_topo) const {
            const int el_pos = pos_topo_vect_(el_id);
            if((el_pos < 0) || (el_pos >= nb_topo)){
                std::ostringstream exc_;
                exc_ << "OneSideContainer::update_topo: element " << el_id << " has position "
                     << el_pos << " in the topology vector, out of range [0, " << nb_topo
                     << "). The stored pos_topo_vect is inconsistent (run check_grid()).";
                throw std::out_of_range(exc_.str());
            }
            return el_pos;
        }

        /**
         * Apply this element's entry of the grid2op topology vector, WITHOUT touching
         * the per-bus element counts. Returns whether anything actually changed.
         *
         * Same contract, and the same reason, as deactivate_no_bus_tracking: a line
         * END does not own its contribution -- `status_global_` gates it -- so when
         * this side belongs to a branch it is the BRANCH that brackets the whole
         * per-element update with the counting, once, around both sides *and*
         * `resolve_status` (which flips that very gate). A standalone container owns
         * its contribution and brackets this itself, in update_topo below.
         */
        bool update_topo_one_el_no_bus_tracking(
            int el_id,
            const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            const Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            DualAlgoControl & solver_control,
            SubstationContainer & substations
        )
        {
            const int el_pos = checked_pos_topo_vect(el_id, static_cast<int>(has_changed.rows()));
            if(!has_changed(el_pos)) return false;
            LocalBusId new_bus = LocalBusId(new_values(el_pos));  // it is a LocalBusId
            if(new_bus.cast_int() < _deactivated_bus_id){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "OneSideContainer::update_topo: bus id should be between -1 and ";
                exc_ << substations.nmax_busbar_per_sub();
                exc_ << " you provided ";
                exc_ << new_bus.cast_int();
                exc_ << ".";
                throw std::out_of_range(exc_.str());
            }
            if(new_bus.cast_int() > substations.nmax_busbar_per_sub()){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "OneSideContainer::update_topo: bus id should be between -1 and ";
                exc_ << substations.nmax_busbar_per_sub();
                exc_ << " you provided ";
                exc_ << new_bus.cast_int();
                exc_ << ".";
                throw std::out_of_range(exc_.str());
            }

            if(new_bus.cast_int() > 0){
                // new bus is a real bus, so i need to make sure to have it turned on, and then change the bus
                if(subid_.size() == 0){
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::update_topo: cannot reconnect element " << el_id
                         << " to a bus: no substation id was ever set for this container "
                         << "(set_subid was never called).";
                    throw std::runtime_error(exc_.str());
                }
                // subid_ is only ever assigned by set_subid() (checked against nb() at the
                // time of the call) or set_osc_state(); neither is re-run when the container
                // is re-initialized with a different element count (init() does not touch
                // subid_), so a container whose element count grew after set_subid() was last
                // called leaves subid_ shorter than the CURRENT nb() -- indexing el_id below
                // would read past its end. Same class of bug _check_pos_topo_vect_filled()
                // already guards against for pos_topo_vect_ (its size() != nb() check); mirror
                // it here.
                if(subid_.size() != nb()){
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::update_topo: cannot reconnect element " << el_id
                         << " to a bus: subid_ has " << subid_.size() << " entries but this "
                         << "container currently has " << nb() << " elements (set_subid was "
                         << "called for a different element count -- call it again after "
                         << "re-initializing this container).";
                    throw std::runtime_error(exc_.str());
                }
                int sub_id = subid_(el_id);
                // `sub_id` feeds local_to_gridmodel's arithmetic (sub_id + (busbar-1)*n_sub),
                // whose OUTPUT is bounds-checked before being stored as this element's bus id
                // -- but an out-of-range `sub_id` can still combine with a valid busbar to land
                // BY COINCIDENCE on another substation's legitimate bus id, silently reconnecting
                // this element to the WRONG bus instead of raising. set_subid() only rejects
                // negative ids (it has no access to n_sub); validate the full range here, where
                // `substations` gives us that context.
                if((sub_id < 0) || (sub_id >= substations.nb_sub())){
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::update_topo: element " << el_id
                         << " has substation id " << sub_id << ", out of range [0, "
                         << substations.nb_sub() << "). The stored subid is inconsistent "
                         << "(run check_grid()).";
                    throw std::out_of_range(exc_.str());
                }
                GridModelBusId new_bus_backend = substations.local_to_gridmodel(sub_id, new_bus);
                bool change_effective = reactivate_no_bus_tracking(el_id, solver_control); // eg reactivate_load(load_id);
                change_effective = change_bus_no_bus_tracking(el_id, new_bus_backend, solver_control, substations) || change_effective; // eg change_bus_load(load_id, new_bus_backend);
                return change_effective;
            } else if (new_bus.cast_int() == _deactivated_bus_id){
                // new bus is negative, we deactivate it
                // the bus is taken out of the system in GridModel.update_topo
                // and a bus is activated if (and only if) one element is connected to it.
                // I must not take `new_bus_backend` out of the system in this case !
                return deactivate_no_bus_tracking(el_id, solver_control);// eg deactivate_load(load_id);
            }
            return false;
        }

        /**
         * Only the values of "new_values" corresponding to "has_changed" == true are used.
         *
         * The bus labelling in "new_values" are local bus (between 1 and n_max_busbar_per_sub).
         */
        virtual std::vector<bool> update_topo(
            const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            const Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            DualAlgoControl & solver_control,
            SubstationContainer & substations
        ) final
        {
            std::vector<bool> res(nb(), false);
            _check_pos_topo_vect_filled();
            const int nb_topo = static_cast<int>(has_changed.rows());
            for(int el_id = 0; el_id < nb(); ++el_id)
            {
                // an entry the caller did not touch mutates nothing, so it must not be
                // bracketed either: taking a contribution away and putting it straight
                // back leaves the counts right, but a bus held by this element alone
                // would transiently hit 0 and report a crossing that never happened.
                if(!has_changed(checked_pos_topo_vect(el_id, nb_topo))) continue;
                // ONE bracket for the whole entry: reactivating an element and then
                // moving it is a single change of which bus it holds, not two.
                _apply_and_track_buses(el_id, substations, solver_control, [&]{
                    res[el_id] = update_topo_one_el_no_bus_tracking(el_id, has_changed,
                                                                   new_values, solver_control,
                                                                   substations);
                });
            }
            return res;
        }

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            std::vector<std::string>,
            std::vector<int>, // bus_id
            std::vector<bool>, // status
            bool,  // has subid info
            std::vector<int>,  // sub_id
            bool,  // has pos_topo_vect info
            std::vector<int>  // pos_topo_vect
            >;

    protected:

        OneSideContainer::StateRes get_osc_state() const  // osc: one side element
        {
            std::vector<int> bus_id(bus_id_.to_int_vector());
            std::vector<bool> status = status_;
            bool has_subid_info = subid_.size();
            std::vector<int> subid(subid_.begin(), subid_.end());
            bool has_topo_vect_info = pos_topo_vect_.size() > 0;
            std::vector<int> pos_topo_vect(pos_topo_vect_.begin(), pos_topo_vect_.end());
            OneSideContainer::StateRes res(
                names_,
                bus_id,
                status,
                has_subid_info,
                subid,
                has_topo_vect_info,
                pos_topo_vect);
            return res;
        }

        void set_osc_state(OneSideContainer::StateRes & my_state)  // osc: one side element
        {
            // read data
            names_ = std::get<0>(my_state);
            std::vector<int> & bus_id = std::get<1>(my_state);
            std::vector<bool> & status = std::get<2>(my_state);
            bool has_subid_info = std::get<3>(my_state);
            bool has_topo_vect_info = std::get<5>(my_state);

            // check sizes
            size_t size = bus_id.size();
            if(names_.size() > 0) check_size(names_, size, "names");  // names are optional
            check_size(bus_id, size, "bus_id");
            check_size(status, size, "status");
            if(has_subid_info)
            {
                const std::vector<int> & subid = std::get<4>(my_state);
                check_size(subid, size, "subid");
                subid_ = IntVect::Map(subid.data(), subid.size());
            }
            if(has_topo_vect_info)
            {
                const std::vector<int> & topo_vect = std::get<6>(my_state);
                check_size(topo_vect, size, "topo_vect");
                pos_topo_vect_ = IntVect::Map(topo_vect.data(), topo_vect.size());
            }

            // input data
            bus_id_ = GlobalBusIdVect(bus_id);
            status_ = status;
        }
        
        void init_osc(
            const Eigen::Ref<const Eigen::VectorXi> & els_bus_id
        )  // osc: one side container
        {
            bus_id_ = GlobalBusIdVect(els_bus_id);
            status_ = std::vector<bool>(els_bus_id.size(), true);
        }

        void set_osc_res_p(){
            const int nb_els = nb();
            for(int el_id = 0; el_id < nb_els; ++el_id){
                if(!status_[el_id]) res_p_[el_id] = 0.;
            }
        }

        void set_osc_res_q(bool ac){
            const int nb_els = nb();
            if(ac){
                for(int el_id = 0; el_id < nb_els; ++el_id){
                    if(!status_[el_id]) res_q_[el_id] = 0.;
                }
            }
            else{
                // no q in DC mode
                for(int el_id = 0; el_id < nb_els; ++el_id) res_q_(el_id) = 0.;
            }
        }

        void reset_osc_results()
        {
            res_p_ = RealVect(nb());  // in MW
            res_q_ =  RealVect(nb());  // in MVar
            res_v_ = RealVect(nb());  // in kV
            res_theta_ = RealVect(nb());  // in deg
            this->_reset_results();
        }

    protected:
        virtual void _reset_results() {
            // nothing to do by default
        };
        virtual void _compute_results(const Eigen::Ref<const RealVect> & /*Va*/,
                                      const Eigen::Ref<const RealVect> & /*Vm*/,
                                      const Eigen::Ref<const CplxVect> & /*V*/,
                                      const SolverBusIdVect & /*id_grid_to_solver*/,
                                      const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
                                      real_type /*sn_mva*/,
                                      bool /*ac*/) {
                                        // nothing to do by default
                                      };
        virtual bool _deactivate(int el_id, DualAlgoControl & /*solver_control*/) {
            // nothing do to by default
            if(status_[el_id]) return true;
            return false;
        };
        virtual bool _reactivate(int el_id, DualAlgoControl & /*solver_control*/) {
            // nothing to do by default
            if(!status_[el_id]) return false;
            return true;
        };
        virtual bool _change_bus(int el_id, GridModelBusId new_bus_id, DualAlgoControl & /*solver_control*/, int /*nb_bus*/) {
            // nothing to do by default
            if(bus_id_(el_id) == new_bus_id) return false;  // nothing to do if the bus did not changed
            return true;
        };
        virtual void _change_p(int /*el_id*/, real_type /*new_p*/, bool /*my_status*/, DualAlgoControl & /*solver_control*/) {
            // nothing to do by default
            };
        virtual void _change_q(int /*el_id*/, real_type /*new_p*/, bool /*my_status*/,DualAlgoControl & /*solver_control*/) {
            // nothing to do by default
        };

    public:
        // Whole-grid semantic validation (see GenericContainer::check_valid / LSGrid::check_grid).
        void check_valid(int nb_bus,
                         int nb_sub,
                         const SubstationContainer & substations,
                         std::vector<int> & all_pos_topo_vect) const override
        {
            check_valid_osc(nb_bus, nb_sub, substations, all_pos_topo_vect, "element");
        }

    protected:
        void _check_pos_topo_vect_filled(){
            if((nb() > 0) && (pos_topo_vect_.size() == 0)){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "update_topo: can only be used if the pos_topo_vect has been set.";
                throw std::runtime_error(exc_.str());
            }
            if(pos_topo_vect_.size() != nb()){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "update_topo: pos_topo_vect_ has not the size of the number of elements: ";
                exc_ << pos_topo_vect_.size() << " vs " << nb() << " elements.";
                throw std::runtime_error(exc_.str());
            }
        }
    protected:
        // used for example when trafo.change_bus_hv need to access 
        GlobalBusIdVect & get_buses_not_const() {return bus_id_;}

        // DANGER zone, neede for trafoContainer and lineContainer
        // because TwoSidesContainer is not fully made
        Eigen::Ref<RealVect> get_res_theta() {return res_theta_;}
        Eigen::Ref<RealVect> get_res_p() {return res_p_;}
        Eigen::Ref<RealVect> get_res_q() {return res_q_;}
        Eigen::Ref<RealVect> get_res_v() {return res_v_;}

        // Range-checks bus_id_ / subid_ / pos_topo_vect_ for a grid with `nb_bus`
        // buses and `nb_sub` substations, and appends the (optional) pos_topo_vect
        // entries to `all_pos_topo_vect`. `el_name` is used only in error messages.
        // NB subid_ and pos_topo_vect_ are optional (may be empty) -- they are only
        // checked when populated.
        void check_valid_osc(int nb_bus,
                             int nb_sub,
                             const SubstationContainer & substations,
                             std::vector<int> & all_pos_topo_vect,
                             const std::string & el_name) const
        {
            const int nb_el = nb();
            const bool has_subid = subid_.size() > 0;          // optional
            const bool has_topo  = pos_topo_vect_.size() > 0;  // optional
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                const int bus = bus_id_(el_id).cast_int();
                const bool connected = status_[el_id];
                if((bus != _deactivated_bus_id) && ((bus < 0) || (bus >= nb_bus)))
                {
                    std::ostringstream exc_;
                    exc_ << "LSGrid::check_grid: " << el_name << " id " << el_id;
                    if(el_id < static_cast<int>(names_.size())) exc_ << " ('" << names_[el_id] << "')";
                    exc_ << " is on bus id " << bus << " which is out of range [0, " << nb_bus
                         << ") (only " << _deactivated_bus_id << " is allowed, meaning disconnected).";
                    throw std::out_of_range(exc_.str());
                }
                if(connected)
                {
                    if(bus == _deactivated_bus_id)
                    {
                        std::ostringstream exc_;
                        exc_ << "LSGrid::check_grid: " << el_name << " id " << el_id
                             << " is marked as connected (its status is true) but its bus id is "
                             << _deactivated_bus_id << " (meaning disconnected).";
                        throw std::runtime_error(exc_.str());
                    }
                    // NB there is no "... and its bus must be active" check here any
                    // more. A bus is active iff an element holds it, and an active
                    // element on bus `bus` IS one such element, so the condition was
                    // a tautology: it could only ever fire when the per-bus counts
                    // had not been established yet, which says nothing about the
                    // grid. What it used to catch -- a saved bus-status vector
                    // contradicting the elements -- cannot be expressed any more.
                }
                if(has_subid)
                {
                    const int sub = subid_(el_id);
                    if((sub < 0) || (sub >= nb_sub))
                    {
                        std::ostringstream exc_;
                        exc_ << "LSGrid::check_grid: " << el_name << " id " << el_id
                             << " has substation id " << sub << " out of range [0, " << nb_sub << ").";
                        throw std::out_of_range(exc_.str());
                    }
                }
                if(has_topo)
                {
                    const int pos = pos_topo_vect_(el_id);
                    if(pos < 0)
                    {
                        std::ostringstream exc_;
                        exc_ << "LSGrid::check_grid: " << el_name << " id " << el_id
                             << " has a negative position in the topology vector (" << pos << ").";
                        throw std::out_of_range(exc_.str());
                    }
                    all_pos_topo_vect.push_back(pos);
                }
            }
        }

    protected:
        // physical properties

        // data for grid2op compat
        IntVect subid_;
        IntVect pos_topo_vect_;

        // input data
        GlobalBusIdVect bus_id_;
        std::vector<bool> status_;

        //output data
        RealVect res_p_;  // in MW
        RealVect res_q_;  // in MVar
        RealVect res_v_;  // in kV
        RealVect res_theta_;  // in degree
};


} // namespace ls2g

#endif  //ONE_SIDE_CONTAINER_H
