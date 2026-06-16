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
            one_el_change_bus_(true)
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
};

/**
DC-specific change-tracking control.

DC power flow only solves the real linear system `Bbus . theta = Pbus`, so it does
not need the full AC change tracking (no pv/pq distinction, no reactive power, etc.).
This holds exactly the subset of flags the DC solver reacts to, expressed in DC terms
(Bbus / Pbus / Vm). It is built from an `AlgoControl` so that all the existing call
sites (which set an `AlgoControl`) keep working unchanged.
**/
class DCControl final
{
    public:
        DCControl() noexcept:
            need_reset_(true),
            dim_changed_(true),
            slack_changed_(true),
            bbus_sparsity_changed_(true),
            bbus_some_zero_(true),
            recompute_bbus_(true),
            recompute_pbus_(true),
            bus_class_changed_(true),
            vm_changed_(true),
            slack_weight_changed_(true)
            {};

        // build the DC control from the (AC oriented) AlgoControl.
        // pv / pq changes are both mapped to "bus classification changed" (DC does not
        // distinguish pv from pq, but both can affect which bus is the reference slack).
        explicit DCControl(const AlgoControl & ac) noexcept:
            need_reset_(ac.need_reset_solver()),
            dim_changed_(ac.has_dimension_changed()),
            slack_changed_(ac.has_slack_participate_changed()),
            bbus_sparsity_changed_(ac.ybus_change_sparsity_pattern()),
            bbus_some_zero_(ac.has_ybus_some_coeffs_zero()),
            recompute_bbus_(ac.need_recompute_ybus()),
            recompute_pbus_(ac.need_recompute_sbus()),
            bus_class_changed_(ac.has_pv_changed() || ac.has_pq_changed()),
            vm_changed_(ac.has_v_changed()),
            slack_weight_changed_(ac.has_slack_weight_changed())
            {};

        ~DCControl() noexcept = default;

        bool need_reset() const {return need_reset_;}
        bool has_dim_changed() const {return dim_changed_;}
        bool has_slack_changed() const {return slack_changed_;}
        bool has_bbus_sparsity_changed() const {return bbus_sparsity_changed_;}
        bool has_bbus_some_zero() const {return bbus_some_zero_;}
        bool need_recompute_bbus() const {return recompute_bbus_;}
        bool need_recompute_pbus() const {return recompute_pbus_;}
        bool has_bus_class_changed() const {return bus_class_changed_;}
        bool has_vm_changed() const {return vm_changed_;}
        bool has_slack_weight_changed() const {return slack_weight_changed_;}

    private:
        bool need_reset_;               // solver needs to be reset from scratch
        bool dim_changed_;              // dimension of Bbus / Pbus changed (eg topology change)
        bool slack_changed_;           // slack participation changed
        bool bbus_sparsity_changed_;   // sparsity pattern of Bbus changed
        bool bbus_some_zero_;          // some coeff of Bbus might have been set to 0
        bool recompute_bbus_;          // some coeff of Bbus changed (same sparsity)
        bool recompute_pbus_;          // some coeff of Pbus changed
        bool bus_class_changed_;       // pv / pq classification changed (affects ref slack)
        bool vm_changed_;              // a generator changed its voltage setpoint (affects Vm)
        bool slack_weight_changed_;    // slack participation weights changed (distributed slack)
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