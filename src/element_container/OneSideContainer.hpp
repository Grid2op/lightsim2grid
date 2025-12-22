// Copyright (c) 2024, RTE (https://www.rte-france.com)
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

// same for all
// - X nb 
// - X get_bus
// - get_buses
// - get_res
// - get_res_full
// - get_theta
// - get_status
// - get_bus_id
// - reconnect_connected_buses
// - update_bus_status
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

                OneSideInfo(const OneSideContainer & r_data_one_side, int my_id):
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
                    if((my_id >= 0) & (my_id < r_data_one_side.nb()))
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
                        if(connected) bus_id = r_data_one_side.bus_id_[my_id];

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
        typedef OneSideInfo DataInfo;

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
        OneSideContainer() {};
        // OneSideInfo get_osc_info(int id_) {return OneSideInfo(*this, id_);}

        // public generic API
        int nb() const { return static_cast<int>(bus_id_.size()); }
        GridModelBusId get_bus(int el_id) const {return _get_bus(el_id, status_, bus_id_);}
        Eigen::Ref<const GlobalBusIdVect> get_buses() const {return bus_id_;}

        tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
        tuple4d get_res_full() const {return tuple4d(res_p_, res_q_, res_v_, res_theta_);}
        
        Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
        const std::vector<bool>& get_status() const {return status_;}
        bool get_status(int el_id) const {return status_.at(el_id);}
        Eigen::Ref<const GlobalBusIdVect> get_bus_id() const {return bus_id_;}
        Eigen::Ref<const IntVect> get_bus_id_numpy() const {
            return IntVect::Map(reinterpret_cast<const int *>(&bus_id_(0)), bus_id_.size());
        }

        void reconnect_connected_buses(SubstationContainer & substation) const{
            const int nb_els = nb();
            for(int el_id = 0; el_id < nb_els; ++el_id)
            {
                if(!status_[el_id]) continue;
                const GlobalBusId my_bus = bus_id_(el_id);
                if(my_bus == _deactivated_bus_id){
                    // TODO DEBUG MODE only this in debug mode
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::reconnect_connected_buses: element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_xxx(...)` ?.";
                    throw std::runtime_error(exc_.str());
                }
                substation.reconnect_bus(my_bus);  // this bus is connected
            }
        }

        void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
            const int nb_el = nb();
            SolverControl unused_solver_control;
            for(int el_id = 0; el_id < nb_el; ++el_id)
            {
                if(!status_[el_id]) continue;
                const GlobalBusId my_bus = bus_id_(el_id);
                if(!busbar_in_main_component[my_bus.cast_int()]){
                    deactivate(el_id, unused_solver_control);
                }
            }    
        }
        void update_bus_status(SubstationContainer & substation) const {
            const int nb_ = nb();
            for(int el_id = 0; el_id < nb_; ++el_id)
            {
                if(!status_[el_id]) continue;
                substation.reconnect_bus(bus_id_[el_id]);
            }
        }    

        void deactivate(int el_id, SolverControl & solver_control) {
            this->_deactivate(el_id, solver_control);
            _generic_deactivate(el_id, status_);
        }
        void reactivate(int el_id, SolverControl & solver_control) {
            this->_reactivate(el_id, solver_control);
            _generic_reactivate(el_id, status_);
        }

        /**
         * This function changes the bus. The bus_id is here given in the
         * "gridmodel" bus.
         * 
         * Not the "solver" bus, nor the "substation" / "local" bus.
         */
        void change_bus(int load_id, GridModelBusId new_gridmodel_bus_id, SolverControl & solver_control, const SubstationContainer & substation) {
            this->_change_bus(load_id, new_gridmodel_bus_id, solver_control, substation.nb_bus());
            _generic_change_bus(load_id, new_gridmodel_bus_id, bus_id_, solver_control, substation.nb_bus());
        }

        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const std::vector<SolverBusId> & id_grid_to_solver,
                             const RealVect & bus_vn_kv,
                             real_type sn_mva,
                             bool ac)
        {
            const int nb_els = nb();
            v_kv_from_vpu(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
            v_deg_from_va(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
            this->_compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        }

        // can be overriden, but has a default behaviour
        virtual void reset_results(){
            reset_osc_results();
        }

        void set_pos_topo_vect(Eigen::Ref<const IntVect> pos_topo_vect)
        {
            pos_topo_vect_.array() = pos_topo_vect;
        }

        void set_subid(Eigen::Ref<const IntVect> subid)
        {
            subid_.array() = subid;
        }

        /**
         * Only the values of "new_values" corresponding to "has_changed" == true are used.
         * 
         * The bus labelling in "new_values" are local bus (between 1 and n_max_busbar_per_sub).
         */
        void update_topo(
            Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            SolverControl & solver_control,
            SubstationContainer & substations
        )
        {
            _check_pos_topo_vect_filled();
            for(int el_id = 0; el_id < nb(); ++el_id)
            {
                int el_pos = pos_topo_vect_(el_id);
                if(!has_changed(el_pos)) continue;
                LocalBusId new_bus = new_values(el_pos);  // it is a LocalBusId
                if(new_bus < _deactivated_bus_id){
                    // TODO DEBUG MODE: only check in debug mode
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::update_topo: bus id should be between -1 and ";
                    exc_ << substations.nmax_busbar_per_sub();
                    exc_ << " you provided ";
                    exc_ << new_bus;
                    exc_ << ".";
                    throw std::out_of_range(exc_.str());
                }
                if(new_bus > substations.nmax_busbar_per_sub()){
                    // TODO DEBUG MODE: only check in debug mode
                    std::ostringstream exc_;
                    exc_ << "OneSideContainer::update_topo: bus id should be between -1 and ";
                    exc_ << substations.nmax_busbar_per_sub();
                    exc_ << " you provided ";
                    exc_ << new_bus;
                    exc_ << ".";
                    throw std::out_of_range(exc_.str());
                }

                if(new_bus > 0){
                    // new bus is a real bus, so i need to make sure to have it turned on, and then change the bus
                    int sub_id = subid_(el_id);
                    GridModelBusId new_bus_backend = substations.local_to_gridmodel(sub_id, new_bus);
                    substations.reconnect_bus(new_bus_backend);
                    reactivate(el_id, solver_control); // eg reactivate_load(load_id);
                    change_bus(el_id, new_bus_backend, solver_control, substations); // eg change_bus_load(load_id, new_bus_backend);
                } else if (new_bus == _deactivated_bus_id){
                    // new bus is negative, we deactivate it
                    deactivate(el_id, solver_control);// eg deactivate_load(load_id);
                    // bus_status_ is set to "false" in GridModel.update_topo
                    // and a bus is activated if (and only if) one element is connected to it.
                    // I must not set `bus_status_[new_bus_backend] = false;` in this case !
                }
            }
        }

        typedef std::tuple<
            std::vector<std::string>,
            std::vector<int>, // bus_id
            std::vector<bool>, // status
            bool,  // has subid info
            std::vector<int>,  // sub_id
            bool,  // has pos_topo_vect info
            std::vector<int>  // pos_topo_vect
            >  StateRes;

    protected:

        OneSideContainer::StateRes get_osc_state() const  // osc: one side element
        {
            std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
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
                subid_ = IntVect::Map(&subid[0], subid.size());
            }
            if(has_topo_vect_info)
            {
                const std::vector<int> & topo_vect = std::get<6>(my_state);
                check_size(topo_vect, size, "topo_vect");
                pos_topo_vect_ = IntVect::Map(&topo_vect[0], topo_vect.size());
            }

            // input data
            bus_id_ = GlobalBusIdVect::Map(reinterpret_cast<GlobalBusId *>(&bus_id[0]), bus_id.size());
            status_ = status;
        }
        
        void init_osc(
            const Eigen::VectorXi & els_bus_id
        )  // osc: one side container
        {
            bus_id_ = els_bus_id.cast<GridModelBusId>();
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
        virtual void _reset_results() {};
        virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<SolverBusId> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac) {};
        virtual void _deactivate(int el_id, SolverControl & solver_control) {};
        virtual void _reactivate(int el_id, SolverControl & solver_control) {};
        virtual void _change_bus(int load_id, GridModelBusId new_bus_id, SolverControl & solver_control, int nb_bus) {};
        virtual void _change_p(int el_id, real_type new_p, bool my_status, SolverControl & solver_control) {};
        virtual void _change_q(int el_id, real_type new_p, bool my_status,SolverControl & solver_control) {};
        // virtual void _change_v(int el_id, real_type new_p, SolverControl & solver_control) {};

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
        Eigen::Ref<GlobalBusIdVect> get_buses_not_const() {return bus_id_;}

        // DANGER zone, neede for trafoContainer and lineContainer
        // because TwoSidesContainer is not fully made
        Eigen::Ref<RealVect> get_res_theta() {return res_theta_;}
        Eigen::Ref<RealVect> get_res_p() {return res_p_;}
        Eigen::Ref<RealVect> get_res_q() {return res_q_;}
        Eigen::Ref<RealVect> get_res_v() {return res_v_;}

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

#endif  //ONE_SIDE_CONTAINER_H
