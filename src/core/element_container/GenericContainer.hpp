// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENERIC_CONTAINER_H
#define GENERIC_CONTAINER_H

#include <algorithm>  // for std::find

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "BaseConstants.hpp"
#include "Container_IteratorUtils.hpp"
#include "SubstationContainer.hpp"

namespace ls2g {

/**
Base class for every object that can be manipulated
**/
class LS2G_API GenericContainer : public BaseConstants
{
    public:

        virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & /*res*/,
                              bool /*ac*/,
                              const SolverBusIdVect & /*id_grid_to_solver*/,
                              real_type /*sn_mva*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                              };

        // Real-valued DC admittance matrix (Bbus) contribution. DC only needs `Bbus . theta = Pbus`
        // (all real), so this fills real triplets directly (no complex temporary).
        // Only "branches" (lines, transformers) contribute to the DC Bbus.
        virtual void fillBdc(std::vector<Eigen::Triplet<real_type> > & /*res*/,
                             const SolverBusIdVect & /*id_grid_to_solver*/,
                             real_type /*sn_mva*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                             };

        virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & /*Bp*/,
                                std::vector<Eigen::Triplet<real_type> > & /*Bpp*/,
                                const SolverBusIdVect & /*id_grid_to_solver*/,
                                real_type /*sn_mva*/,
                                FDPFMethod /*xb_or_bx*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                                };
                                
        virtual void fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & /*Bf*/,
                                     const SolverBusIdVect & /*id_grid_to_solver*/,
                                     real_type /*sn_mva*/,
                                     int /*nb_line*/,
                                     bool /*transpose*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                                };

        virtual void fillSbus(Eigen::Ref<CplxVect> /*Sbus*/, const SolverBusIdVect & /*id_grid_to_solver*/, bool /*ac*/) const {
                                // nothing to do by default
                                // is overriden mainly for "one side elements" (loads, generators etc.)
                                };
        virtual void fillpv(std::vector<int>& /*bus_pv*/,
                            std::vector<bool> & /*has_bus_been_added*/,
                            const SolverBusIdVect& /*slack_bus_id_solver*/,
                            const SolverBusIdVect & /*id_grid_to_solver*/) const {
                                // nothing to do by default
                                // is overriden mainly for "generators"
                            };
        
        virtual void get_q(std::vector<real_type>& /*q_by_bus*/) {
                                // nothing to do by default
                                // is overriden mainly for "generators"
                                };
        
        virtual void set_p_slack(const Eigen::Ref<const RealVect>& /*node_mismatch*/, const SolverBusIdVect & /*id_grid_to_solver*/) {
                                // nothing to do by default
                                // is overriden mainly for "generators"
                                };
    
        static const int _deactivated_bus_id;
        // ---- which buses does an element hold alive? ---------------------------
        /**
         * Apply this element's contribution to `substation`'s per-bus element count:
         * `sign > 0` adds what it holds NOW, `sign < 0` removes it. `crossed` is
         * OR-ed with true for every bus that thereby crossed 0 -- ie every bus that
         * entered or left the solved system.
         *
         * THIS IS THE ONLY STATEMENT OF "WHICH BUSES DOES THIS ELEMENT HOLD".
         * The from-scratch recount is built from it and every mutator tracks its own
         * effect through it (see `_apply_and_track_buses`); there is no longer a
         * `reconnect_connected_buses` stating the same rule a second time. Anything
         * that states it twice can drift, and the failure mode is a bus labelling
         * that no longer matches the matrices -- a converged, plausible, wrong
         * answer.
         *
         * The rule is genuinely different per container, which is why this is
         * virtual and not a helper:
         *   - one-sided elements hold their bus iff they are active;
         *   - a line or transformer is gated by `status_global_` FIRST, and only
         *     then does each side count -- a line end cannot know that, which is
         *     why this cannot live in OneSideContainer alone;
         *   - an HVDC line has no such gate: each converter station stands alone.
         */
        virtual void contribute_to_buses(int /*el_id*/,
                                         SubstationContainer & /*substation*/,
                                         int /*sign*/,
                                         bool & /*crossed*/) const {
                                // nothing to do by default
                                };

        /**
         * Run `mutate`, and tell `substation` which buses this element stopped or
         * started holding as a result.
         *
         * Written once, used by every mutator that can change bus membership. It
         * never restates the contribution rule: it asks `contribute_to_buses` to
         * take the element's current contribution away, lets the mutation happen,
         * then asks it to put the new one back. Whatever the container's rule is,
         * and whatever the mutation did, the counts end up right.
         *
         * `crossed` is what the caller turns into `tell_dimension_changed()`. Every
         * OTHER change flag -- ybus sparsity, ybus values, sbus -- stays where it
         * has always been decided, inside the container's own `_deactivate` /
         * `_reactivate` / `_change_bus`, which `mutate` calls.
         *
         * Note a caller that hands us a no-op (deactivating an already-inactive
         * element, changing a bus to itself) gets the right COUNTS either way -- it
         * removes and re-adds the same contribution -- but a bus that is alone would
         * transiently hit 0 and report a crossing that did not happen, costing a
         * needless rebuild. Mutators guard the no-op cases before calling in.
         */
        template<class Mutation>
        void _apply_and_track_buses(int el_id,
                                    SubstationContainer & substation,
                                    DualAlgoControl & solver_control,
                                    Mutation && mutate){
            bool crossed = false;
            contribute_to_buses(el_id, substation, -1, crossed);
            // Everything a caller can get wrong -- the element id, the bus id -- is
            // checked by the mutators BEFORE they call in here (_check_in_range and
            // _check_new_bus_id), so a call the grid refuses never reaches this
            // bracket and the counts it would have left short are never touched.
            //
            // This used to be a try / catch that put the contribution back on the way
            // out of an exception. It was not free: an unwind edge through this header
            // made GCC keep every std::vector<bool> access in fillYbus live across it,
            // for 4.9M instructions per rebuild solve of case9241pegase, in a function
            // that never calls any of this. An exception from deeper inside a mutation
            // can still leave the counts short; that is a grid which must be rebuilt
            // and its caches dropped, not one to carry on solving with.
            mutate();
            contribute_to_buses(el_id, substation, +1, crossed);
            if(crossed){
                solver_control.ac_algo_controler().tell_dimension_changed();
                solver_control.dc_algo_controler().tell_dimension_changed();
            }
        }

        /**computes the total amount of power for each bus (for generator only)**/
        virtual void gen_p_per_bus(std::vector<real_type> & /*res*/) const {
                                // nothing to do by default
                                // is overriden mainly for "one side elements" (loads, generators etc.)
                                };
        virtual void nb_line_end(std::vector<int> & /*res*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                                };
        virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & /*res*/) const {
                                // nothing to do by default
                                // is overriden mainly for "branches" (lines, transformers etc.)
                                };
        virtual void disconnect_if_not_in_main_component(std::vector<bool> & /*busbar_in_main_component*/, SubstationContainer & /*substation*/, DualAlgoControl & /*solver_control*/) {
                                // nothing to do by default
                                };

        /**
        Whole-grid semantic validation (see LSGrid::check_grid).

        Checks that every index this container carries (bus ids, substation ids,
        position in the topology vector, slack references...) is in range for a
        grid with `nb_bus` buses and `nb_sub` substations. Throws std::out_of_range
        on a bad index and std::runtime_error on a structural error.

        `all_pos_topo_vect` is an accumulator: each container appends the
        `pos_topo_vect` values it actually carries (the field is optional and may
        be empty), so LSGrid can afterwards check they form a valid permutation.

        The default does nothing; element containers override it.
        **/
        virtual void check_valid(int /*nb_bus*/,
                                 int /*nb_sub*/,
                                 const SubstationContainer & /*substations*/,
                                 std::vector<int> & /*all_pos_topo_vect*/) const {
                                // nothing to validate by default
                                };

        void set_names(const std::vector<std::string> & names){
            names_ = names;
        }
        // empty if set_names() was never called on this container
        const std::vector<std::string> & get_names() const {
            return names_;
        }

        /**"define" the destructor for compliance with clang (otherwise lots of warnings)**/
        GenericContainer() noexcept = default;
        virtual ~GenericContainer() noexcept = default;
        
        /**
        check the size of the elements
        **/
        template<class T, class intType>
        static void check_size(const T & container, intType size, const std::string & container_name)
        {
            if(static_cast<intType>(container.size()) != size) throw std::runtime_error(container_name + " do not have the proper size");
        }

        /**
        activation / deactivation of elements
        **/
        static void _generic_reactivate(int el_id, std::vector<bool> & status);
        static void _generic_deactivate(int el_id, std::vector<bool> & status);


        /**
        check if an element is in a vector or an Eigen Vector, do not use for other types of containers (might not be efficient at all)
        **/
        template<class ScalarCLS, class VectCLS>  // a std::vector, or an Eigen::Vector                                                 
        static bool is_in_vect(const ScalarCLS & val, const VectCLS & cont) {
            return std::find(
                cont.begin(),
                cont.end(),
                static_cast<typename VectCLS::value_type>(val)) != cont.end();}

        // Bounds check for an element id that came from OUTSIDE this library --
        // a python call, a grid2op action, anything a user chose. Always
        // compiled in: such an id is never trusted, and the alternative to
        // throwing is an out-of-bounds read on a bit-packed std::vector<bool>,
        // which in a release wheel is a segfault or silent garbage rather than
        // an error.
        template<typename Cont, typename FunName, typename IntType>
        // todo automatically "unwrap" IntType to be either cont::size_type for stl container and
        // Eigen::Index for Eigen containers
        static void _check_in_range(IntType el_id, const Cont & cont, FunName fun_name="")
        {
            if(el_id >= static_cast<IntType>(cont.size()))
            {
                std::ostringstream exc_;
                exc_ << "GenericContainer::"<<fun_name<<": Cannot access element with id";
                exc_ << el_id;
                exc_ << " while the grid counts ";
                exc_ << cont.size();
                exc_ << " such elements (id too high)";
                throw std::out_of_range(exc_.str());
            }
            if(el_id < 0)
            {
                std::ostringstream exc_;
                exc_ << "GenericContainer::"<< fun_name <<" Cannot change the bus of element with id ";
                exc_ << el_id;
                exc_ << " (id should be >= 0)";
                throw std::out_of_range(exc_.str());
            }
        }

        // The same check for an id THIS library generated: the counter of a loop
        // over a container's own elements, which cannot be out of range unless
        // there is a bug here rather than in the caller. Compiled out with
        // NDEBUG; the assertion and sanitizer CI builds keep it, since
        // USE_DEBUG_ASSERTS clears NDEBUG (see src/core/CMakeLists.txt).
        //
        // It was not free. On a 9241-bus grid the internal callers alone reached
        // it 128k times per powerflow -- four passes over both ends of every
        // branch, in fillYbus, compute_results and reconnect_connected_buses --
        // and it plus the _get_bus it guards were ~1.9% of a solve, spent
        // re-deriving a bound the loop already guarantees. Note the check is
        // also a function call there, not an inlined compare: the error paths
        // build an ostringstream, so the compiler keeps the whole thing
        // out of line.
        template<typename Cont, typename FunName, typename IntType>
        static void _check_in_range_internal(IntType el_id, const Cont & cont, FunName fun_name="")
        {
#ifndef NDEBUG
            _check_in_range(el_id, cont, fun_name);
#else
            (void) el_id;
            (void) cont;
            (void) fun_name;
#endif
        }

    protected:
        std::vector<std::string> names_;

    protected:
        
        /**
         * Change the bus of the element "el_id" and performs some basic check that the new bus is valid.
         * 
         * The new_gridmodel_bus_id is given in the gridmodel convention, between 0 and `n_sub * n_busbar_per_sub` 
         */
        /**
         * Validate a caller-supplied bus id, BEFORE anything is mutated.
         *
         * This used to live inside `_generic_change_bus`, which runs from inside the
         * `_apply_and_track_buses` bracket -- so a bus id the grid was going to refuse
         * was rejected only after the element's contribution had been taken away.
         * `_apply_and_track_buses` compensated with a try / catch that put the
         * contribution back on the way out.
         *
         * Checking first is both simpler and cheaper. Simpler because a refused call
         * never touches the counts at all, rather than touching them and undoing it
         * with a restore that has to reason about a half-applied mutation. Cheaper
         * because the catch cost far more than the exceptional path it protected: an
         * unwind edge through this header made GCC keep every `std::vector<bool>`
         * access in `fillYbus` live across it, for 4.9M instructions per rebuild solve
         * of case9241pegase -- in a function that never calls any of this.
         *
         * Always active: the id comes from the caller, so this is a user-facing check.
         */
        static void _check_new_bus_id(const GridModelBusId & new_gridmodel_bus_id, int nb_max_bus)
        {
            if(new_gridmodel_bus_id.cast_int() >= nb_max_bus)
            {
                std::ostringstream exc_;
                exc_ << "GenericContainer::_change_bus: Cannot change an element to bus ";
                exc_ << new_gridmodel_bus_id.cast_int();
                exc_ << " There are only ";
                exc_ << nb_max_bus;
                exc_ << " distinct buses on this grid.";
                throw std::out_of_range(exc_.str());
            }
            if(new_gridmodel_bus_id.cast_int() < 0)
            {
                std::ostringstream exc_;
                exc_ << "GenericContainer::_change_bus: new bus id should be >=0 and not ";
                exc_ << new_gridmodel_bus_id.cast_int();
                throw std::out_of_range(exc_.str());
            }
        }

        void _generic_change_bus(
            int el_id,
            const GridModelBusId & new_gridmodel_bus_id,
            GlobalBusIdVect &  el_bus_ids,
            DualAlgoControl & solver_control,
            int nb_bus) const;
        // el_id supplied by a caller: bounds-checked, throws out_of_range.
        GridModelBusId _get_bus(int el_id, const std::vector<bool> & status_, const GlobalBusIdVect & bus_id_) const;
        // el_id produced by one of our own loops over this container: the bound
        // is a property of the loop, so it is only asserted (see
        // _check_in_range_internal). Defined here rather than in the .cpp so the
        // loops that call it per element can actually inline it.
        GridModelBusId _get_bus_internal(int el_id, const std::vector<bool> & status_, const GlobalBusIdVect & bus_id_) const
        {
            _check_in_range_internal(static_cast<std::vector<bool>::size_type>(el_id),
                                     status_,
                                     "_get_bus_internal");
            if(!status_[el_id]) return GridModelBusId(_deactivated_bus_id);
            return bus_id_(el_id);
        }

        /**
        compute the amps from the p, the q and the v (v should NOT be pair unit)
        **/
        void _get_amps(Eigen::Ref<RealVect> a,
                       const Eigen::Ref<const RealVect> & p,
                       const Eigen::Ref<const RealVect> & q,
                       const Eigen::Ref<const RealVect> & v) const;

        /**
        Convert this container's bus voltages to the per-element results: v from
        pu to kV, theta from rad to degrees.

        This was two functions, v_kv_from_vpu and v_deg_from_va, called back to
        back on the same elements. Everything before the last line of each was
        the same work -- read the element's bus, map it to a solver bus, check
        both are connected -- so the walk is done once and both results are
        written from it.
        **/
        void v_kv_theta_from_vpu(const Eigen::Ref<const RealVect> & Va,
                                 const Eigen::Ref<const RealVect> & Vm,
                                 const std::vector<bool> & status,
                                 int nb_element,
                                 const GlobalBusIdVect & bus_me_id,
                                 const SolverBusIdVect & id_grid_to_solver,
                                 const Eigen::Ref<const RealVect> & bus_vn_kv,
                                 Eigen::Ref<RealVect> v,
                                 Eigen::Ref<RealVect> theta) const;
};


} // namespace ls2g

#endif // GENERIC_CONTAINER_H
