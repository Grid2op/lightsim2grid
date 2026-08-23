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
#include <cassert>
#include <vector>
// #include <set>
#include <stdio.h>
#include <cstdint> // for int32, int64
#include <limits>
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
            nb_elements_per_bus_(n_bus_max_, 0),
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
            // positivity + no int overflow, same guarantee as init_bus: `i + j * n_sub_`
            // below is int arithmetic
            const int nb_bus_checked = checked_nb_bus(n_sub_, nmax_busbar_per_sub_, "SubstationContainer::init_sub");
            // Both destinations are written with an unchecked operator[] below, and
            // neither is guaranteed to be sized yet: only the two-argument
            // constructor sizes sub_vn_kv_, and the live construction path is the
            // default constructor followed by init_bus(), which leaves sub_vn_kv_
            // EMPTY. Writing n_sub_ doubles into it would be an out-of-bounds heap
            // write. Size both here instead of assuming.
            if(sub_vn_kv_.size() != n_sub_) sub_vn_kv_ = RealVect::Zero(n_sub_);
            if(bus_vn_kv_.size() != nb_bus_checked) bus_vn_kv_ = RealVect::Zero(nb_bus_checked);
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
            // NB: a band is NOT required to bracket the bus nominal voltage. Real IIDM data
            // does violate it (a 380 kV level declared with an operating range of
            // [390, 450] kV is ordinary practice on the European 400 kV network, and shows
            // up on real data), so rejecting it here would refuse to load perfectly usable
            // grids. ContingencyAnalysis's `violation_threshold` copes by clamping its
            // interpolation anchor into [vmin, vmax] instead -- see check_bus_voltage_violations.
            bus_vmin_kv_ = bus_vmin_kv;
            bus_vmax_kv_ = bus_vmax_kv;
        }
        Eigen::Ref<const RealVect> get_bus_vmin_kv() const {return bus_vmin_kv_;}
        Eigen::Ref<const RealVect> get_bus_vmax_kv() const {return bus_vmax_kv_;}

        bool is_bus_connected(const GridModelBusId & global_bus_id) const {return bus_status_[global_bus_id.cast_int()];}
        bool is_bus_connected(int sub_id, const LocalBusId & local_bus_id) const {
            return bus_status_[local_to_gridmodel(sub_id, local_bus_id).cast_int()];
        }

        /**
         * `n_sub * nmax_busbar_per_sub` -- the total bus count -- computed safely.
         *
         * Both factors reach us from outside the C++ core (init_bus / init_sub
         * arguments, or a pickle / binary file via set_state), and the product is
         * stored in an `int` (n_bus_max_) and used in `int` index arithmetic
         * (`i + j * n_sub`). Computing it directly overflows for large inputs, which
         * is UB and yields a negative bus count. Do the multiplication in 64 bits and
         * refuse anything that does not fit, so no caller has to remember to.
         * `caller` only names the function in the error message.
         */
        static int checked_nb_bus(int n_sub, int nmax_busbar_per_sub, const char * caller)
        {
            if(n_sub <= 0){
                std::ostringstream exc_;
                exc_ << caller << ": the number of substations must be strictly positive, got "
                     << n_sub << " (it is also used as a modulo divisor in sub_id_of_bus, so 0 "
                     << "would be a division by zero).";
                throw std::runtime_error(exc_.str());
            }
            if(nmax_busbar_per_sub <= 0){
                std::ostringstream exc_;
                exc_ << caller << ": the maximum number of busbars per substation must be strictly "
                     << "positive, got " << nmax_busbar_per_sub << ".";
                throw std::runtime_error(exc_.str());
            }
            const std::int64_t total = static_cast<std::int64_t>(n_sub) *
                                       static_cast<std::int64_t>(nmax_busbar_per_sub);
            if(total > static_cast<std::int64_t>(std::numeric_limits<int>::max())){
                std::ostringstream exc_;
                exc_ << caller << ": n_sub (" << n_sub << ") * nmax_busbar_per_sub ("
                     << nmax_busbar_per_sub << ") = " << total << " buses, which does not fit in the "
                     << "int used to index them.";
                throw std::runtime_error(exc_.str());
            }
            return static_cast<int>(total);
        }

        /**
         * Every bus nominal voltage must be a finite, strictly positive number.
         *
         * NB this is deliberately only about `bus_vn_kv_` / `sub_vn_kv_`, never about
         * bus_vmin_kv_ / bus_vmax_kv_: those use NaN as the documented "no limit set
         * for this bus" sentinel. Same reason the element-level finiteness checks were
         * dropped earlier (non-finite thermal limits / reactive bounds are legitimate
         * sentinels on real grids) -- a NOMINAL voltage has no such convention, 0 or a
         * negative value is simply meaningless and would silently produce nonsense
         * per-unit conversions (v_kv_from_vpu multiplies by it).
         */
        void check_vn_kv_positive() const
        {
            const auto check_one = [](real_type v, int id, const char * what){
                if(!std::isfinite(v) || (v <= 0.)){
                    std::ostringstream exc_;
                    exc_ << "SubstationContainer: " << what << " " << id << " has a nominal voltage of "
                         << v << " kV; it must be a finite, strictly positive number.";
                    throw std::runtime_error(exc_.str());
                }
            };
            for(Eigen::Index i = 0; i < bus_vn_kv_.size(); ++i) check_one(bus_vn_kv_(i), static_cast<int>(i), "bus");
            for(Eigen::Index i = 0; i < sub_vn_kv_.size(); ++i) check_one(sub_vn_kv_(i), static_cast<int>(i), "substation");
        }

        void init_bus(int n_sub, int nmax_busbar_per_sub, const Eigen::Ref<const RealVect> & bus_vn_kv)
        {
            // positivity + no int overflow on n_sub * nmax_busbar_per_sub
            const int nb_bus_total = checked_nb_bus(n_sub, nmax_busbar_per_sub, "Substation::init_bus");
            if(bus_vn_kv.size() != nb_bus_total){
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
            reset_bus_element_counts();  // recounted from the elements by LSGrid::init_bus_status()

            // check that a "substation" always has the same vn_kv for all of its buses
            check_bus_vn_kv_uniform_per_sub();
            check_vn_kv_positive();
        }

        /**
         * Every bus of a substation carries that substation's nominal voltage.
         *
         * With `n` = n_sub_ and `k` = nmax_busbar_per_sub_, the buses of substation
         * `I` are the ids `I, I + n, I + 2n, ..., I + (k-1)n` -- so the invariant is
         *
         *     for all 0 <= i < n, for all 0 <= j < k:  bus_vn_kv[i] == bus_vn_kv[i + j*n]
         *
         * and that is exactly the loop below, written as plain index arithmetic on
         * bus_vn_kv_ (the stride is `n`, NOT 1 and not `k`: consecutive entries of
         * bus_vn_kv_ are different SUBSTATIONS, which is why eg the case14 fixture
         * legitimately has bus_vn_kv_[4] == 138 next to bus_vn_kv_[5] == 20).
         *
         * NB this used to be expressed via local_to_gridmodel(sub_id, LocalBusId(..)),
         * which is where it went wrong: LocalBusId is 1-based (bus id ==
         * `sub_id + (local_bus_id - 1) * n`), and the loop ran local ids `[1, k)`, so
         * it opened by comparing bus `i` with itself and stopped one busbar short.
         * At the usual k == 2 that compared bus `i` to bus `i` and checked nothing at
         * all -- the invariant was silently unenforced on every path.
         */
        void check_bus_vn_kv_uniform_per_sub() const
        {
            const int n = n_sub_;
            const int k = nmax_busbar_per_sub_;
            for(int i = 0; i < n; ++i){
                const real_type ref_vn_kv = bus_vn_kv_(i);
                for(int j = 1; j < k; ++j)
                {
                    const int bus_id = i + j * n;
                    const real_type this_bus_vn_kv = bus_vn_kv_(bus_id);
                    if(std::abs(this_bus_vn_kv - ref_vn_kv) > BaseConstants::_tol_equal_float){
                        std::ostringstream exc_;
                        exc_ << "SubstationContainer: each bus of a substation must have the same nominal "
                             << "voltage. Substation " << i << " has " << ref_vn_kv
                             << " kV on bus id " << i << " but " << this_bus_vn_kv
                             << " kV on bus id " << bus_id << " (its busbar " << (j + 1) << ").";
                        throw std::runtime_error(exc_.str());
                    }
                }
            }
        }
        // ---- per-bus element count ---------------------------------------------
        /**
         * How many elements currently hold each bus alive.
         *
         * DERIVED state, deliberately absent from StateRes: it is a function of the
         * elements, so it is recomputed (recompute_bus_element_counts, driven by
         * LSGrid) rather than serialized. That is also the recovery path -- if the
         * incremental bookkeeping ever drifts, any rebuild restores it from the
         * elements themselves.
         *
         * A bus is in the solved system iff its count is non-zero, so the two
         * transitions that matter are 0 -> 1 and 1 -> 0: those, and only those,
         * change the DIMENSION of the system and renumber every bus after them.
         * `bus_gained_element` / `bus_lost_element` report exactly that; what the
         * caller does with it (raising `tell_dimension_changed`) is the caller's
         * business, and nothing else about the change flags is decided here.
         */
        const std::vector<std::size_t> & get_nb_elements_per_bus() const {return nb_elements_per_bus_;}

        /**
         * Have the counts been established from the elements at least once?
         *
         * They start at zero on a grid that has only just been built, which is
         * indistinguishable from "every bus is empty" -- so until something has
         * counted, tracking a decrement would underflow on the first
         * `deactivate_load(...)`. Tracking is therefore inert until then, and
         * `recompute_bus_element_counts()` is what turns it on. That is not a
         * special case bolted on: it is the same "recount whenever there is no
         * cache to protect" rule that makes drift unable to outlive an
         * invalidation, applied to the moment before the first powerflow.
         */
        bool bus_counts_ready() const {return bus_counts_ready_;}
        void mark_bus_counts_ready(){bus_counts_ready_ = true;}

        /// one more element holds `global_bus_id`; true iff it just crossed 0 -> 1
        bool bus_gained_element(const GlobalBusId & global_bus_id){
            if(!bus_counts_ready_) return false;
            return ++nb_elements_per_bus_[global_bus_id.cast_int()] == 1;
        }
        /// one fewer element holds `global_bus_id`; true iff it just crossed 1 -> 0
        bool bus_lost_element(const GlobalBusId & global_bus_id){
            if(!bus_counts_ready_) return false;
            // a count that would go negative means an element was removed twice, ie
            // the bookkeeping has drifted. Not something to paper over silently.
            assert(nb_elements_per_bus_[global_bus_id.cast_int()] > 0 &&
                   "SubstationContainer::bus_lost_element: count already 0 -- an element "
                   "was removed from this bus twice, the incremental counts have drifted");
            return --nb_elements_per_bus_[global_bus_id.cast_int()] == 0;
        }
        /// back to "nothing holds any bus", ready to be recounted from the elements
        void reset_bus_element_counts(){
            nb_elements_per_bus_.assign(nb_bus(), 0);
            bus_counts_ready_ = false;   // re-established by recompute_bus_element_counts()
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
        /**
         * TODO sub_vn_kv_ IS DEAD STATE: nothing ever sets it, and nothing reads it.
         *
         * The only two writers are `init_sub()` and the two-argument constructor
         * `SubstationContainer(n_sub, nmax_busbar_per_sub)` -- and NEITHER is called
         * anywhere in the codebase (nor exposed to python). Every grid is built by
         * the default constructor followed by `init_bus()`, which fills bus_vn_kv_
         * only, so this vector is EMPTY on every grid produced by
         * init_from_pandapower / _pypowsybl / _matpower / _powermodels, and empty in
         * every binary file saved so far (including the format-4 fixture in
         * lightsim2grid/tests/binary_format_fixture/). It is nevertheless part of
         * StateRes, so it is serialized -- as an empty vector -- into every pickle
         * and every binary file.
         *
         * It is also fully REDUNDANT: a substation's nominal voltage is by
         * definition the common vn_kv of its buses (see
         * check_bus_vn_kv_uniform_per_sub), i.e. sub_vn_kv_[I] == bus_vn_kv_[I] for
         * 0 <= I < n_sub_. check_valid() enforces exactly that whenever it is
         * non-empty.
         *
         * So it must stay OPTIONAL: requiring it whenever bus_vn_kv_ is set would
         * reject every grid and every saved file that exists today.
         *
         * Decide before a release which way to go, and do not leave it as is:
         *   - either DROP it from StateRes (it is derivable from bus_vn_kv_ and read
         *     by nobody) -- needs a BINARY_FORMAT_VERSION bump; or
         *   - have `init_bus()` populate it (`sub_vn_kv_ = bus_vn_kv.head(n_sub)`)
         *     and make check_valid() require it -- this changes what newly-saved
         *     files contain, and the requirement can only be turned on once no
         *     already-saved file needs to load.
         */
        RealVect sub_vn_kv_;
        std::vector<bool> bus_status_;
        // derived from the elements, never serialized: see get_nb_elements_per_bus()
        std::vector<std::size_t> nb_elements_per_bus_;
        bool bus_counts_ready_ = false;  // see bus_counts_ready()
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