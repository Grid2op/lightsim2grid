// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef UTILS_H
#define UTILS_H

/**
Some typedef and other structures define here and used everywhere else
**/
#include <iostream>
#include <complex>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "TaggedIdVec.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

// using real_type = double;  // type for real numbers: can be changed if installed from source
using real_type = double;  // type for real numbers: can be changed if installed from source

using cplx_type = std::complex<real_type>;  // type for complex number

using EigenPythonNumType = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;  // Eigen::VectorXd
using IntVect = Eigen::Matrix<int, Eigen::Dynamic, 1>;
using RealVect = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
using CplxVect = Eigen::Matrix<cplx_type, Eigen::Dynamic, 1>;

using tuple3d = std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> >;
using tuple4d = std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> >;
using tuple5d = std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> >;

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;
using CplxMat = Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic> ;

using EigenRefConstCplxSpMat = Eigen::Ref<const Eigen::SparseMatrix<cplx_type> >;
using EigenRefConstRealSpMat = Eigen::Ref<const Eigen::SparseMatrix<real_type> >;

// type of error in the different solvers
enum class ErrorType {NoError,
                      SingularMatrix,
                      TooManyIterations,
                      InifiniteValue,
                      SolverAnalyze,
                      SolverFactor,
                      SolverReFactor,
                      SolverSolve,
                      NotInitError,
                      LicenseError};
std::ostream& operator<<(std::ostream& out, const ErrorType & error_type);

// Escape (and truncate to 64 chars) a string of untrusted origin -- read from a
// possibly-corrupted file, or supplied by a plugin -- before embedding it in an
// exception message. pybind11 converts what() to a python str as UTF-8, so raw
// garbage bytes would turn the intended RuntimeError into a UnicodeDecodeError;
// control characters would also let a name inject newlines or terminal escape
// sequences into logs. Truncation keeps a corrupted length from producing a
// message megabytes long.
LS2G_API std::string printable(const std::string & s);


struct Coeff{
    Eigen::Index row_id;
    Eigen::Index col_id;
    cplx_type value;
};

// define some constant for compilation outside of "setup.py"
#ifndef VERSION_MAJOR
#define VERSION_MAJOR "-1"
#endif

#ifndef VERSION_MEDIUM
#define VERSION_MEDIUM "-1"
#endif

#ifndef VERSION_MINOR
#define VERSION_MINOR "-1"
#endif

class AlgoControl final
{
    public:
        AlgoControl() noexcept:
            change_dimension_(true),
            pv_changed_(true),
            pq_changed_(true),
            slack_participate_changed_(true),
            need_reset_solver_(true), 
            need_recompute_sbus_(true),
            need_recompute_ybus_(true),
            v_changed_(true),
            slack_weight_changed_(true),
            ybus_some_coeffs_zero_(true),
            ybus_change_sparsity_pattern_(true),
            one_el_change_bus_(true),
            cache_maybe_poisoned_(true),
            voltage_control_changed_(true)
            {};

        ~AlgoControl() noexcept = default;

        void tell_all_changed(){
            change_dimension_ = true;
            pv_changed_ = true;
            pq_changed_ = true;
            slack_participate_changed_ = true;
            need_reset_solver_ = true;
            need_recompute_sbus_ = true;
            need_recompute_ybus_ = true;
            v_changed_ = true;
            slack_weight_changed_ = true;
            ybus_some_coeffs_zero_ = true;
            ybus_change_sparsity_pattern_ = true;
            one_el_change_bus_ = true;
            cache_maybe_poisoned_ = true;
            voltage_control_changed_ = true;
        }

        /**
         * Is there nothing outstanding -- has every change this control tracks
         * already been consumed by a solve?
         *
         * The exact negation of tell_all_changed(), and the only honest way to ask
         * "is a cache built against this grid still valid?". A cache cannot answer
         * that by looking at itself: changing a line's impedance, a tap, or an
         * injection leaves every vector size and every bus status exactly as it was.
         * Only these flags know. See LSGrid::unset_changes().
         */
        [[nodiscard]] bool nothing_changed() const noexcept {
            return !change_dimension_ && !pv_changed_ && !pq_changed_ &&
                   !slack_participate_changed_ && !need_reset_solver_ &&
                   !need_recompute_sbus_ && !need_recompute_ybus_ && !v_changed_ &&
                   !slack_weight_changed_ && !ybus_some_coeffs_zero_ &&
                   !ybus_change_sparsity_pattern_ && !one_el_change_bus_ &&
                   !cache_maybe_poisoned_ && !voltage_control_changed_;
        }

        void tell_none_changed(){
            change_dimension_ = false;
            pv_changed_ = false;
            pq_changed_ = false;
            slack_participate_changed_ = false;
            need_reset_solver_ = false;
            need_recompute_sbus_ = false;
            need_recompute_ybus_ = false;
            v_changed_ = false;
            slack_weight_changed_ = false;
            ybus_some_coeffs_zero_ = false;
            ybus_change_sparsity_pattern_ = false;
            one_el_change_bus_ = false;
            cache_maybe_poisoned_ = false;
            voltage_control_changed_ = false;
        }

        // the dimension of the Ybus matrix / Sbus vector has changed (eg. topology changes)
        void tell_dimension_changed(){change_dimension_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some pv generators are now pq or the opposite
        void tell_pv_changed(){pv_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some pq nodes are now pv or the opposite
        void tell_pq_changed(){pq_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some generators that participated to the slack bus now do not, or the opposite
        void tell_slack_participate_changed(){slack_participate_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // ybus need to be recomputed for some reason
        void tell_recompute_ybus(){need_recompute_ybus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // sbus need to be recomputed for some reason
        void tell_recompute_sbus(){need_recompute_sbus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // solver needs to be reset from scratch for some reason
        void tell_solver_need_reset(){need_reset_solver_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // the sparsity pattern of ybus changed
        void tell_ybus_change_sparsity_pattern(){ybus_change_sparsity_pattern_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // tell at least one generator changed its v setpoint
        void tell_v_changed(){v_changed_ = true;}
        // at least one generator has changed its slack participation
        void tell_slack_weight_changed(){slack_weight_changed_ = true;}
        // tells that some coeff of ybus might have been set to 0.
        // (and ybus compressed again, so these coeffs are really completely hidden)
        // might need to trigger some recomputation of some solvers (eg NR based ones)
        void tell_ybus_some_coeffs_zero(){ybus_some_coeffs_zero_ = true;}
        void tell_one_el_changed_bus(){one_el_change_bus_ = true;}
        /**
         * A voltage SETPOINT a control group's bordered rows carry has moved: a
         * regulating generator's or hvdc converter station's target magnitude.
         *
         * Deliberately narrow, and deliberately NOT raised for everything the
         * voltage-control plan is made of. Who is in a group, and which bus a group
         * regulates, are the same inputs the pv/pq split reads -- and since that split
         * IS a layer of the plan (VoltageControlPlan::build_pv_pq), a change to one is
         * a change to the other, already carried by `pv_changed_` and friends. See
         * `need_recompute_pv_pq()`. What those flags do NOT carry is a setpoint: moving
         * a remote regulator's target changes no bus' pv/pq class at all -- that is the
         * point of the bordered formulation, the regulated bus stays PQ -- so it needs
         * a flag of its own, and this is it.
         *
         * If you add an input to the plan that the pv/pq split does not read, raise
         * this from the modifier that moves it.
         */
        void tell_voltage_control_changed(){voltage_control_changed_ = true;}
        /**
         * The per-bus element counts may no longer be what the elements say.
         *
         * Deliberately NOT the same question as `need_reset_solver()`, which asks
         * whether the SOLVER-SIDE data has to be rebuilt. Everything that data is
         * made of is derived from the elements when it is rebuilt, so a reset costs
         * time and nothing else. The counts are different: they are maintained
         * incrementally (+1 / -1 in GenericContainer::_apply_and_track_buses), they
         * are never recomputed on the ordinary path, and since bus connectivity IS
         * the counts a wrong one is not a slow path but a different grid. Rebuilding
         * them is O(all elements), so it must be asked for by the thing that can
         * actually make them wrong -- not by every caller who merely wants a fresh
         * solve.
         *
         * Raised by: construction, `tell_all_changed()` (so: reset, a copy, set_state,
         * a powerflow that threw part way through), and LSGrid's public
         * `tell_bus_counts_maybe_poisoned()`, which is what a caller who mutated the
         * containers behind LSGrid's back -- or who caught an exception out of a
         * mutator -- uses to say so. NOT raised by `prevent_cache_reuse()` /
         * `tell_solver_need_reset()`: those say the solver data is stale, which says
         * nothing about the counts.
         */
        void tell_cache_maybe_poisoned(){
            cache_maybe_poisoned_ = true;
            // Poisoned counts imply a solver reset, and the implication only runs this
            // way. The counts decide which buses exist, so EVERYTHING the solver side
            // is made of is built on them -- the bus labelling, the dimension of Ybus,
            // the pv-pq split, the slack weights. Recounting while re-stamping the rest
            // would repair the counts and then solve the old bus set anyway, which is
            // the same wrong grid with a tidier bookkeeping. The converse does not
            // hold, and that is the whole point of having two flags: a solver reset
            // says nothing about the counts.
            need_reset_solver_ = true;
        }

        bool has_dimension_changed() const {return change_dimension_;}
        bool has_pv_changed() const {return pv_changed_;}
        bool has_pq_changed() const {return pq_changed_;}
        bool has_slack_participate_changed() const {return slack_participate_changed_;}
        bool need_reset_solver() const {return need_reset_solver_;}
        bool need_recompute_sbus() const {return need_recompute_sbus_;}
        bool need_recompute_ybus() const {return need_recompute_ybus_;}
        bool ybus_change_sparsity_pattern() const {return ybus_change_sparsity_pattern_;}
        bool has_slack_weight_changed() const {return slack_weight_changed_;}
        bool has_v_changed() const {return v_changed_;}
        bool has_ybus_some_coeffs_zero() const {return ybus_some_coeffs_zero_;}
        bool has_one_el_changed_bus() const {return one_el_change_bus_;}
        // see tell_cache_maybe_poisoned()
        bool cache_maybe_poisoned() const {return cache_maybe_poisoned_;}
        // see tell_voltage_control_changed()
        bool has_voltage_control_changed() const {return voltage_control_changed_;}

        /**
         * Must the pv/pq split be rebuilt?
         *
         * Every term is a reason the split itself changes: the system was rebuilt from
         * scratch (`need_reset_solver_`, which `tell_cache_maybe_poisoned()` implies),
         * the bus set changed (`change_dimension_`), the slack set moved
         * (`slack_participate_changed_`, and the slack is not PV), or a bus changed
         * class (`pv_changed_` / `pq_changed_`).
         *
         * `ybus_change_sparsity_pattern_` is deliberately NOT a term, though it looks
         * like one: it is the flag for "the bus LABELLING may have moved". It is only
         * ever raised by a BRANCH-side mutation (reconnecting a line or a trafo,
         * moving one of its ends), and every one of those goes through
         * `GenericContainer::_apply_and_track_buses`, which raises `change_dimension_`
         * exactly when such a mutation empties or fills a bus -- which is exactly when
         * the labelling moves. If no bus crossed, `id_me_to_solver` is unchanged, and a
         * branch is neither a voltage controller nor a slack, so the split it produces
         * is identical. The term was there, and dropping it was measured, not argued:
         * see the `[pv_pq]` cases in `test_cache_reuse.cpp`, which put the grid in a
         * state where this flag is the ONLY term raised and check the reused split
         * against a cold one. With the term dropped they still pass; the same
         * experiment run on `slack_participate_changed_` fails, which is why that one
         * stays.
         *
         * Read by `LSGrid::_build_into_cache`, which is the only thing that rebuilds it.
         */
        [[nodiscard]] bool need_recompute_pv_pq() const noexcept {
            return need_reset_solver_ || change_dimension_ ||
                   slack_participate_changed_ || pv_changed_ || pq_changed_;
        }

        /**
         * Must the voltage-control plan be rebuilt, or does the one the previous solve
         * of this family left in its cache still describe the grid?
         *
         * The pv/pq split is a LAYER of that plan (VoltageControlPlan::build_pv_pq), so
         * whatever rebuilds the split rebuilds the plan -- that is the first half, and
         * it is not a grab-bag of loosely-related flags but literally the same question
         * asked one layer down. Whoever is in a control group, and which bus it
         * regulates, are exactly the inputs the split reads; there is no way to change
         * one without changing the other.
         *
         * The second half is what the split does NOT read: a setpoint. See
         * `tell_voltage_control_changed()`.
         *
         * What is deliberately in NEITHER half is the set of changes an ordinary
         * grid2op step makes: moving a load's P and Q raises `need_recompute_sbus_` and
         * nothing else, so the plan survives such a step untouched -- which is the
         * whole point of caching it.
         *
         * All-or-nothing on purpose: the plan's layers are rebuilt together or not at
         * all. Refreshing only the setpoints of an otherwise unchanged group layout is
         * a finer question this does not try to answer yet.
         */
        [[nodiscard]] bool need_recompute_voltage_control() const noexcept {
            return need_recompute_pv_pq() || voltage_control_changed_;
        }

    private:    
        bool change_dimension_;
        bool pv_changed_;
        bool pq_changed_;
        bool slack_participate_changed_;
        bool need_reset_solver_;  // some matrices change size, needs to be computed
        bool need_recompute_sbus_;  // some coeff of sbus changed, need to recompute it
        bool need_recompute_ybus_;  // some coeff of ybus changed, but not its sparsity pattern
        bool v_changed_;
        bool slack_weight_changed_;
        bool ybus_some_coeffs_zero_;  // tells that some coeff of ybus might have been set to 0. (and ybus compressed again, so these coeffs are really completely hidden)
        bool ybus_change_sparsity_pattern_;  // sparsity pattern of ybus changed (and so are its coeff), or ybus change of dimension
        bool one_el_change_bus_;  // whether one element has change of bus (or being reconnected / disconnected)
        bool cache_maybe_poisoned_;  // the per-bus element counts may have drifted: see tell_cache_maybe_poisoned()
        bool voltage_control_changed_;  // an input of the voltage-control plan moved: see tell_voltage_control_changed()
};

/**
Change-tracking control for both solver families (AC and DC).

A grid modification (eg. disconnecting a line, changing a setpoint) invalidates the cached
matrices of *both* the AC and the DC solver. `DualAlgoControl` simply holds one independent
`AlgoControl` per family so each solver keeps its own change tracking (the AC solver consumes
and resets `ac_algo_controler()` on an AC powerflow, the DC solver consumes and resets
`dc_algo_controler()` on a DC powerflow, without clobbering each other).

It is a plain composition (no inheritance, no virtual dispatch). Almost every change concerns
both families the same way, so the `tell_xxx()` forwarders below raise a flag on both at once,
which is what the element containers use. A change that only one family reads (the DC Sbus
term of a phase shifter, the AC-only voltage-control plan) names its family explicitly:
    dual.dc_algo_controler().tell_recompute_sbus();
**/
class DualAlgoControl final
{
    public:
        DualAlgoControl() noexcept = default;
        ~DualAlgoControl() noexcept = default;

        AlgoControl & ac_algo_controler() noexcept {return ac_algo_controler_;}
        AlgoControl & dc_algo_controler() noexcept {return dc_algo_controler_;}
        const AlgoControl & ac_algo_controler() const noexcept {return ac_algo_controler_;}
        const AlgoControl & dc_algo_controler() const noexcept {return dc_algo_controler_;}

        // raise the same flag on both families (see the AlgoControl method of the same name)
        void tell_dimension_changed() noexcept {ac_algo_controler_.tell_dimension_changed(); dc_algo_controler_.tell_dimension_changed();}
        void tell_pv_changed() noexcept {ac_algo_controler_.tell_pv_changed(); dc_algo_controler_.tell_pv_changed();}
        void tell_pq_changed() noexcept {ac_algo_controler_.tell_pq_changed(); dc_algo_controler_.tell_pq_changed();}
        void tell_slack_participate_changed() noexcept {ac_algo_controler_.tell_slack_participate_changed(); dc_algo_controler_.tell_slack_participate_changed();}
        void tell_recompute_ybus() noexcept {ac_algo_controler_.tell_recompute_ybus(); dc_algo_controler_.tell_recompute_ybus();}
        void tell_recompute_sbus() noexcept {ac_algo_controler_.tell_recompute_sbus(); dc_algo_controler_.tell_recompute_sbus();}
        void tell_solver_need_reset() noexcept {ac_algo_controler_.tell_solver_need_reset(); dc_algo_controler_.tell_solver_need_reset();}
        void tell_ybus_change_sparsity_pattern() noexcept {ac_algo_controler_.tell_ybus_change_sparsity_pattern(); dc_algo_controler_.tell_ybus_change_sparsity_pattern();}
        void tell_v_changed() noexcept {ac_algo_controler_.tell_v_changed(); dc_algo_controler_.tell_v_changed();}
        void tell_slack_weight_changed() noexcept {ac_algo_controler_.tell_slack_weight_changed(); dc_algo_controler_.tell_slack_weight_changed();}
        void tell_ybus_some_coeffs_zero() noexcept {ac_algo_controler_.tell_ybus_some_coeffs_zero(); dc_algo_controler_.tell_ybus_some_coeffs_zero();}
        void tell_one_el_changed_bus() noexcept {ac_algo_controler_.tell_one_el_changed_bus(); dc_algo_controler_.tell_one_el_changed_bus();}

    private:
        AlgoControl ac_algo_controler_;  // change tracking consumed by the AC solver
        AlgoControl dc_algo_controler_;  // change tracking consumed by the DC solver
};

template<int U>
class IntClass
{
    private:
        int m_bus_id;

    public:
        explicit IntClass(const int bus_id) noexcept : m_bus_id(bus_id) {} 
        IntClass() noexcept : m_bus_id(0) {}

        // rule of 5 for the same IntClass
        IntClass(const IntClass & oth) noexcept = default;
        IntClass(IntClass && oth) noexcept = default;
        IntClass& operator=(const IntClass & oth) noexcept = default;
        IntClass& operator=(IntClass && oth) noexcept = default;
        IntClass& operator=(const int bus_id) noexcept{
            return IntClass(bus_id);
        }
        ~ IntClass() noexcept = default;

        // cast to int (if needed)
        const int & cast_int() const & noexcept {return m_bus_id;}
        const int & cast_int() const && noexcept {return m_bus_id;}
        int cast_int() & noexcept {return m_bus_id;}  // do a copy in this case

        // automatic conversion to int
        explicit operator int() const noexcept {return m_bus_id;}
        bool operator==(const IntClass & oth) const noexcept{
            return m_bus_id == oth.m_bus_id;
        }
        bool operator!=(const IntClass & oth) const noexcept{
            return !(m_bus_id == oth.m_bus_id);
        }

        // prevent conversion between the different IntClass
        template<int V>
        IntClass(const IntClass<V> &) = delete;
        template<int V>
        IntClass(IntClass<V> &&) = delete;
        template<int V>
        IntClass& operator=(const IntClass<V> &) = delete;
        template<int V>
        IntClass& operator=(IntClass<V> &&) = delete;

};
        
/**
 * This class is used to store the id of a bus in the local "convention".
 * 
 * This means that this number is dependant of the substation.
 * 
 * It is typically between -1 (for disconnected) or between 1 and "n_max_busbar_per_sub"
 */
using LocalBusId = IntClass<LOCAL_BUS>;
static_assert(sizeof(LocalBusId)==sizeof(int));  // make sure I can safely "reinterpret_cast" LocalBusId to int and vice versa

/**
 * The classes GridModelBusId and GlobalBusId (both are the same type in cpp) denotes
 * buses number expressed in "GridModel bus" or in "GlobalBus".
 * 
 * GlobalBus / GridModelBus are typically between 0 and `n_sub * n_max_busbar_per_sub`
 */
using GridModelBusId = IntClass<GRIDMODEL_BUS>;
using GlobalBusId = IntClass<GLOBAL_BUS>;
static_assert(sizeof(GlobalBusId)==sizeof(int));  // make sure I can safely "reinterpret_cast" GlobalBusId to int and vice versa


/**
 * These are "solver bus id". They depends on the actuall topolgy of the grid.
 *
 * They are used to pass information from the gridmodel to the solver and should not be
 * used anywhere else.
 */
using SolverBusId = IntClass<SOLVER_BUS>;
static_assert(sizeof(SolverBusId)==sizeof(int));  // make sure I can safely "reinterpret_cast" SolverBusId to int and vice versa

} // namespace ls2g

#endif // UTILS_H