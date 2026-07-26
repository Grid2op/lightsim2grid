// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SUBSTATIONCONTAINER_H
#define SUBSTATIONCONTAINER_H

#include <iostream>
#include <vector>
// #include <set>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI, std::abs
#include <sstream>
#include <stdexcept>

#include "Utils.hpp"
#include "BaseConstants.hpp"
#include "element_container/Container_IteratorUtils.hpp"

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

namespace ls2g {


class SubstationContainer;

class SubstationInfo
{
    public:
        int id;
        std::string name;
        int nb_max_busbars;
        real_type vn_kv;

        inline SubstationInfo(const SubstationContainer & r_data, int my_id);
};

class LS2G_API SubstationContainer final : public IteratorAdder<SubstationContainer, SubstationInfo>
{
    friend class SubstationInfo;

    public:
        using DataInfo = SubstationInfo;
        
    public:

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            int,  // n_sub_
            int,  // nmax_busbar_per_sub
            // sub_vn_kv_: OPTIONAL, and empty in practice on every grid produced
            // today. It is only ever filled by init_sub() / the two-argument
            // constructor, neither of which is called anywhere (the live path is
            // the default constructor + init_bus(), which fills bus_vn_kv_ only),
            // and nothing reads it back. It is therefore NOT valid to require it
            // whenever bus_vn_kv_ is set: doing so rejects every grid built by
            // init_from_pandapower / _pypowsybl / _matpower / _powermodels, and
            // every already-saved binary file -- including this repo's own
            // format-4 compatibility fixture (tests/binary_format_fixture).
            std::vector<real_type>, // sub_vn_kv_ (optional, see above)
            std::vector<bool>,  // bus_status_;
            std::vector<real_type>,  // bus_vn_kv_;
            std::vector<std::string>,  // sub_names_
            std::vector<real_type>,  // bus_vmin_kv_ (optional, empty if unset)
            std::vector<real_type>  // bus_vmax_kv_ (optional, empty if unset)
            >;
        
        int nb() const {return n_sub_;}

        SubstationContainer::StateRes get_state() const;
        void set_state(SubstationContainer::StateRes & my_state);

        /**
         * Whole-container consistency validation, the substation half of
         * LSGrid::check_grid() (see there). set_state() already enforces this on
         * anything coming from a pickle / a binary file; this is what covers a
         * container reached the other way -- built element by element through the
         * python API, or mutated afterwards. Throws std::runtime_error, returns
         * normally for a consistent container.
         */
        void check_valid() const;

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static SubstationContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "SubstationContainer"; }  // written into / checked against the binary file header

        SubstationContainer() noexcept:
            n_sub_(-1), 
            nmax_busbar_per_sub_(-1),
            n_bus_max_(-1){}
            
        SubstationContainer(int n_sub, int nmax_busbar_per_sub) noexcept:
            n_sub_(n_sub),
            nmax_busbar_per_sub_(nmax_busbar_per_sub),
            n_bus_max_(n_sub * nmax_busbar_per_sub),
            sub_vn_kv_(n_sub),
            bus_status_(n_bus_max_, false),
            bus_vn_kv_(n_bus_max_){}
            
        ~SubstationContainer() noexcept = default;

        void reset_bus_status(){
            // iterate the vector actually being written, not n_bus_max_ (a separate
            // field): set_state() now guarantees they agree, but not depending on
            // that invariant costs nothing and cannot go wrong if it is ever broken.
            const std::size_t nb_bus_total = bus_status_.size();
            for(std::size_t i = 0; i < nb_bus_total; ++ i) bus_status_[i] = false;
        }

        void init_sub(const Eigen::Ref<const RealVect> & sub_vn_kv){
            if((n_sub_ <= 0) || (nmax_busbar_per_sub_ <= 0)){
                throw std::runtime_error("SubstationContainer::init_sub: the substations must be declared "
                                         "first (see init_bus / the two-argument constructor).");
            }
            if(sub_vn_kv.size() != n_sub_){
                throw std::range_error("SubstationContainer::init_sub: sub_vn_kv should have the size of the number of substations on the grid.");
            }
            // Both destinations are written with an unchecked operator[] below, and
            // neither is guaranteed to be sized yet: only the two-argument
            // constructor sizes sub_vn_kv_, and the live construction path is the
            // default constructor followed by init_bus(), which leaves sub_vn_kv_
            // EMPTY. Writing n_sub_ doubles into it would be an out-of-bounds heap
            // write. Size both here instead of assuming.
            const Eigen::Index nb_bus_total =
                static_cast<Eigen::Index>(n_sub_) * static_cast<Eigen::Index>(nmax_busbar_per_sub_);
            if(sub_vn_kv_.size() != n_sub_) sub_vn_kv_ = RealVect::Zero(n_sub_);
            if(bus_vn_kv_.size() != nb_bus_total) bus_vn_kv_ = RealVect::Zero(nb_bus_total);
            for(int i = 0; i < n_sub_; ++i){
                // store substation vn kv
                sub_vn_kv_[i] = sub_vn_kv[i];
                for(int j = 0; j < nmax_busbar_per_sub_; ++j){
                    // store the buses vn kv, for all buses
                    bus_vn_kv_[i + j * n_sub_] = sub_vn_kv[i];
                }
            }

        }

        void init_sub_names(const std::vector<std::string> & sub_names){
            if(sub_names.size() != static_cast<size_t>(n_sub_)){
                throw std::runtime_error("Wrong number of substation when setting their names.");
            }
            sub_names_ = sub_names;
        }
        const std::vector<std::string> & get_sub_names() const {
            return sub_names_;
        }

        // return the number of substations on the grid
        int nb_sub() const {return n_sub_;}

        // return the maximum number of possible buses on the grid
        size_t nb_bus() const {return bus_vn_kv_.size();}
        
        /**
        Retrieve the number of connected buses
        **/
        size_t nb_connected_bus() const
        {
            size_t res = 0;
            for(const auto & el : bus_status_)
            {
                if(el) ++res;
            }
            return res;
        }

        
        int nmax_busbar_per_sub() const {return nmax_busbar_per_sub_;}

        Eigen::Ref<const RealVect> get_bus_vn_kv() const {return bus_vn_kv_;}

        // per-bus min/max operating voltage (kV), optional: empty if never set (e.g. pandapower-origin grids)
        void init_bus_voltage_limits(const Eigen::Ref<const RealVect> & bus_vmin_kv, const Eigen::Ref<const RealVect> & bus_vmax_kv){
            if(static_cast<size_t>(bus_vmin_kv.size()) != nb_bus()){
                throw std::runtime_error("SubstationContainer::init_bus_voltage_limits: bus_vmin_kv does not have the proper size.");
            }
            if(static_cast<size_t>(bus_vmax_kv.size()) != nb_bus()){
                throw std::runtime_error("SubstationContainer::init_bus_voltage_limits: bus_vmax_kv does not have the proper size.");
            }
            bus_vmin_kv_ = bus_vmin_kv;
            bus_vmax_kv_ = bus_vmax_kv;
        }
        Eigen::Ref<const RealVect> get_bus_vmin_kv() const {return bus_vmin_kv_;}
        Eigen::Ref<const RealVect> get_bus_vmax_kv() const {return bus_vmax_kv_;}

        bool is_bus_connected(const GridModelBusId & global_bus_id) const {return bus_status_[global_bus_id.cast_int()];}
        bool is_bus_connected(int sub_id, const LocalBusId & local_bus_id) const {
            return bus_status_[local_to_gridmodel(sub_id, local_bus_id).cast_int()];
        }

        void init_bus(int n_sub, int nmax_busbar_per_sub, const Eigen::Ref<const RealVect> & bus_vn_kv)
        {
            if(bus_vn_kv.size() != n_sub * nmax_busbar_per_sub){
                std::ostringstream exc_;
                exc_ << "Substation::init_bus: ";
                exc_ << "your model counts ";
                exc_ << n_sub_  << " substations and ";
                exc_ << nmax_busbar_per_sub_  << " maximum busbars per substation. But you provided a bus_vn_kv with ";
                exc_ << bus_vn_kv.size() << " elements.";
                exc_ << "Both should match.";
                throw std::runtime_error(exc_.str());
            }

            n_sub_ = n_sub;
            nmax_busbar_per_sub_ = nmax_busbar_per_sub;
            bus_vn_kv_ = bus_vn_kv;  // base_kv

            bus_status_ = std::vector<bool>(nb_bus(), true); // by default everything is connected

            // check that a "substation" always has the same vn_kv for all of its buses
            check_bus_vn_kv_uniform_per_sub();
        }

        /**
         * Every busbar of a substation must carry the same nominal voltage.
         *
         * A gridmodel bus id is laid out as `sub_id + (local_bus_id - 1) * n_sub_`
         * (see local_to_gridmodel), so the buses of one substation are strided by
         * **n_sub_**, not by nmax_busbar_per_sub_, and LocalBusId is 1-BASED:
         * local id 1 is the first busbar (gridmodel id == sub_id), local id
         * nmax_busbar_per_sub_ is the last one.
         *
         * That 1-based convention is what made the original loop here a no-op: it
         * ran local ids `[1, nmax_busbar_per_sub_)`, so it started by comparing the
         * reference bus with itself and stopped one busbar short -- with the usual
         * nmax_busbar_per_sub_ == 2 it compared bus `sub_id` to bus `sub_id` and
         * checked nothing at all. Hence the range below: local ids 2 .. nmax
         * inclusive, i.e. every busbar *after* the reference one.
         */
        void check_bus_vn_kv_uniform_per_sub() const
        {
            for(int sub_id = 0; sub_id < n_sub_; sub_id++){
                const real_type ref_vn_kv = bus_vn_kv_(sub_id);
                for(int local_bus_id = 2; local_bus_id <= nmax_busbar_per_sub_; local_bus_id++)
                {
                    const int bus_id = local_to_gridmodel(sub_id, LocalBusId(local_bus_id)).cast_int();
                    const real_type this_bus_vn_kv = bus_vn_kv_(bus_id);
                    if(std::abs(this_bus_vn_kv - ref_vn_kv) > BaseConstants::_tol_equal_float){
                        std::ostringstream exc_;
                        exc_ << "SubstationContainer: each bus of a substation must have the same nominal "
                             << "voltage. Substation " << sub_id << " has " << ref_vn_kv
                             << " kV on its busbar 1 (bus id " << sub_id << ") but " << this_bus_vn_kv
                             << " kV on its busbar " << local_bus_id << " (bus id " << bus_id << ").";
                        throw std::runtime_error(exc_.str());
                    }
                }
            }
        }
        const std::vector<bool> & get_bus_status() const {return bus_status_;}
        
        void disconnect_all_buses(){
            const size_t nb_bus_total = nb_bus();
            for(size_t i = 0; i < nb_bus_total; ++i) bus_status_[i] = false;
        }
        void reconnect_bus(const GridModelBusId& global_bus_id){
            bus_status_[global_bus_id.cast_int()] = true;
        }
        void reconnect_bus(int sub_id, const LocalBusId & local_bus_id){
            reconnect_bus(local_to_gridmodel(sub_id, local_bus_id));
        }
        void disconnect_bus(const GridModelBusId & global_bus_id){
            bus_status_[global_bus_id.cast_int()] = false;
        }
        void disconnect_bus(int sub_id, const LocalBusId & local_bus_id){
            disconnect_bus(local_to_gridmodel(sub_id, local_bus_id));
        }
        GridModelBusId local_to_gridmodel(int sub_id, const LocalBusId & local_bus_id) const{
            if(local_bus_id.cast_int() == BaseConstants::_deactivated_bus_id){
                return GlobalBusId(BaseConstants::_deactivated_bus_id);
            }
            if(local_bus_id.cast_int() == 0){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "SubstationContainer::local_to_gridmodel at this stage, local_bus_id should not be 0.";
                exc_ << " A local_bus_id should be either -1 (for disconnected) or between 1 and n_bus_max_per_sub";
                throw std::runtime_error(exc_.str());
            }
            if(local_bus_id.cast_int() > nmax_busbar_per_sub_){
                // TODO DEBUG MODE: only check in debug mode
                std::ostringstream exc_;
                exc_ << "SubstationContainer::local_to_gridmodel at this stage, ";
                exc_ << "local_bus_id should be < ";
                exc_ << nmax_busbar_per_sub_ << "(nmax_busbar_per_sub, max number of buses per substations) ";
                exc_ << "But you provided " << local_bus_id.cast_int();
                throw std::runtime_error(exc_.str());
            }

            return GridModelBusId(sub_id + (local_bus_id.cast_int() - 1) * n_sub_);
        }

        // reverse of local_to_gridmodel's substation part: a gridmodel bus id is laid out as
        // `sub_id + (local_bus_id - 1) * n_sub_`, so `% n_sub_` recovers the substation id
        // regardless of which local busbar the bus is.
        int sub_id_of_bus(int gridmodel_bus_id) const {
            return gridmodel_bus_id % n_sub_;
        }

    private:
        int n_sub_;
        int nmax_busbar_per_sub_;
        int n_bus_max_;
        RealVect sub_vn_kv_;
        std::vector<bool> bus_status_;
        RealVect bus_vn_kv_;
        std::vector<std::string> sub_names_;
        RealVect bus_vmin_kv_;  // optional, empty if unset
        RealVect bus_vmax_kv_;  // optional, empty if unset

};

inline SubstationInfo::SubstationInfo(const SubstationContainer & r_data, int my_id):
    id(my_id),
    name(""),
    nb_max_busbars(-1),
    vn_kv(-1.)
{
    if(my_id < 0) return;
    if(my_id >= r_data.nb()) return;
    // sub_names_ is OPTIONAL (empty unless set_substation_names() was called, which
    // eg the pandapower / matpower / powermodels loaders never do) and bus_vn_kv_ is
    // empty on a substation container that was declared but never given its buses:
    // both must be bounds-checked here, `my_id < nb()` only bounds them against
    // n_sub_. Same reason as everywhere else: -O3 -DNDEBUG release wheels have no
    // std::vector / Eigen bounds check to fall back on.
    if(my_id < static_cast<int>(r_data.sub_names_.size())) name = r_data.sub_names_[my_id];
    nb_max_busbars = r_data.nmax_busbar_per_sub_;
    if(my_id < static_cast<int>(r_data.bus_vn_kv_.size())) vn_kv = r_data.bus_vn_kv_[my_id];
}


} // namespace ls2g

#endif // SUBSTATIONCONTAINER_H