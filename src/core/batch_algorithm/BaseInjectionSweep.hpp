// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEINJECTIONSWEEP_H
#define BASEINJECTIONSWEEP_H

#include "BaseBatchSolverSynch.hpp"

#include <exception>

namespace ls2g {

/**
How each element of a batch is initialized, ie which complex voltage the solver starts
from when it solves step `i`.

This is a compile time property of the algorithm, not a user setting: it is what
distinguishes the two aliases declared at the bottom of this file. What the *seed* itself
is (the caller's `Vinit` or the converged "n" powerflow) remains a runtime setting, shared
by every batch algorithm: see BaseBatchSolverSynch::set_init_from_n_powerflow.
 **/
enum class LS2G_API BatchInitKind
{
    /** warm start: step `i` starts from the solution of step `i-1` (see TimeSeries) **/
    FromPreviousStep,
    /** every step restarts from the same seed, so the steps are independent
        of one another and of their ordering (see InjectionSweep) **/
    FromSeed
};

/**
Batch of powerflows over a *fixed grid topology* and a varying injection: one step per row
of the (generator / static generator / load) matrices handed to `compute_Vs`.

Nothing here assumes the steps are consecutive instants -- that assumption belongs to the
`TimeSeries` alias alone, and is carried entirely by `INIT`, which selects how each step is
initialized. Everything else -- the inputs, the results, the whole public interface -- is
identical between the two instantiations aliased at the bottom of this file.
 **/
template<BatchInitKind INIT>
class LS2G_API BaseInjectionSweep: public BaseBatchSolverSynch
{
    public:
        explicit BaseInjectionSweep(const LSGrid & init_grid_model):
            BaseBatchSolverSynch(init_grid_model),
            _Sbuses(),
            _status(1), // 1: success, 0: failure
            _compute_flows(true)
            {}
        ~BaseInjectionSweep() noexcept override = default;

        BaseInjectionSweep(const BaseInjectionSweep&) = delete;
        BaseInjectionSweep(BaseInjectionSweep&&) = delete;
        BaseInjectionSweep & operator=(BaseInjectionSweep&&) = delete;
        BaseInjectionSweep & operator=(const BaseInjectionSweep&) = delete;

        // Name of this instantiation, used in the error messages below. The code is
        // shared by both aliases, so hardcoding "TimeSeries" there would mislead an
        // InjectionSweep user (and vice versa).
        static const char * algo_name() {
            return INIT == BatchInitKind::FromPreviousStep ? "TimeSeries" : "InjectionSweep";
        }

        // see BaseBatchSolverSynch::supports_multithread. Derived from INIT rather than
        // set independently, so the two can never disagree: chaining the steps is
        // exactly what makes them impossible to split over threads.
        bool supports_multithread() const override {
            return INIT != BatchInitKind::FromPreviousStep;
        }

        // control on whether I compute the flows or not
        void deactivate_flow_computations() {_compute_flows = false;}
        void activate_flow_computations() {_compute_flows = true;}

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}

        // status
        int get_status() const {return _status;}

        /**
        This function computes the results of running as many powerflow when varying the 
        injection (Sbus). 

        Each line of `Sbuses` will be a time step, and each column with 
        **/
        int compute_Vs(const Eigen::Ref<const RealMat> & gen_p,
                       const Eigen::Ref<const RealMat> & sgen_p,
                       const Eigen::Ref<const RealMat> & load_p,
                       const Eigen::Ref<const RealMat> & load_q,
                       const Eigen::Ref<const CplxVect> & Vinit,
                       const int max_iter,
                       const real_type tol);
                       
        void clear() override {
            BaseBatchSolverSynch::clear();
            _Sbuses = CplxMat();
            _status = 1;
            _compute_flows = true;
            _timer_total = 0.;
            _timer_pre_proc = 0.;
        }

        Eigen::Ref<const CplxMat > get_sbuses() const {return _Sbuses;}
        Eigen::Ref<RealMat > compute_flows() {
            compute_flows_from_Vs();
            return _amps_flows;
        }
        Eigen::Ref<RealMat > compute_power_flows() {
            compute_flows_from_Vs(false);
            return _active_power_flows;
        }

    protected:
        /**
         * Per-unit injection that the four per-step matrices of `compute_Vs` do NOT
         * account for, and which is therefore constant across the steps: storage
         * units, REACTIVE_POWER-mode SVCs, the reactive setpoint of a generator whose
         * voltage regulation is off, a static generator's reactive setpoint, the hvdc
         * injections and (dc only) the phase-shifter term.
         *
         * Computed as `complete - accounted_for`: the complete injection the
         * gridmodel built in `Sbus_` (ac) / `Pbus_` (dc), minus the same
         * reconstruction `compute_Vs` performs, evaluated at the gridmodel's own
         * target values. Deriving it this way -- instead of listing the elements --
         * is what keeps it from falling out of date when a new element type starts
         * contributing to LSGrid::fillSbus_me.
         *
         * Only valid after prepare_solver_input_base() has run.
         */
        CplxVect _constant_sbus_pu(bool ac_solver_used) const;

        /**
         * Solve the steps in [step_begin, step_end) with the passed solver (one per
         * thread) and its own AlgoControl, writing directly to the (disjoint) rows of
         * _voltages; book-keeping goes to the passed accumulators, merged by the caller
         * once every thread has joined. This is the body shared by the single- and
         * multi-threaded paths.
         *
         * Unlike ContingencyAnalysis::run_contingency_range there is no per-thread Ybus:
         * the topology is fixed here, so every thread simply reads the member Ybus_
         * (ac) / Bbus_ (dc), which no one modifies for the whole duration of the call.
         *
         * A step that does not converge stops this range (the remaining rows keep their
         * zeros) and is reported in `first_diverging_step`; any thrown exception is
         * captured into `err` instead of crossing the thread boundary.
         */
        void run_step_range(size_t step_begin,
                            size_t step_end,
                            AlgorithmSelector & algo,
                            AlgoControl & control,
                            const Eigen::Ref<const CplxVect> & Vinit_solver,
                            bool ac_solver_used,
                            int max_iter,
                            real_type tol_solver,  // already divided by sn_mva
                            int & nb_solved,
                            double & timer_solver,
                            int & first_diverging_step,
                            std::exception_ptr & err,
                            bool needs_solver_init);

        template<class T>
        void fill_SBus_real(Eigen::Ref<CplxMat> Sbuses,
                            const T & structure_data,
                            const Eigen::Ref<const RealMat> & temporal_data,
                            const SolverBusIdVect & id_me_to_ac_solver,
                            bool add  // if true call += else calls -=
                            ) const 
        {
            size_t nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            SolverBusId bus_id_solver;
            GlobalBusId bus_id_me;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                if(bus_id_me.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::fill_SBus_real: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                bus_id_solver = id_me_to_ac_solver[bus_id_me.cast_int()];
                if(bus_id_solver.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::fill_SBus_real: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver.cast_int()) += tmp;
                else Sbuses.col(bus_id_solver.cast_int()) -= tmp;
            }
        }

        template<class T>
        void fill_SBus_imag(Eigen::Ref<CplxMat> Sbuses,
                            const T & structure_data,
                            const Eigen::Ref<const RealMat> & temporal_data,
                            const SolverBusIdVect & id_me_to_ac_solver,
                            bool add  // if true call += else calls -=
                            ) const 
        {
            size_t nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            SolverBusId  bus_id_solver;
            GlobalBusId bus_id_me;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                if(bus_id_me.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::fill_SBus_imag: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                bus_id_solver = id_me_to_ac_solver[bus_id_me.cast_int()];
                if(bus_id_solver.cast_int() == _deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::fill_SBus_imag: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver.cast_int()) += BaseConstants::my_i * tmp;
                else Sbuses.col(bus_id_solver.cast_int()) -= BaseConstants::my_i * tmp;
            }
        }

    private:
        // inputs
        CplxMat _Sbuses;

        // outputs
        int _status;

        // parameters
        bool _compute_flows;

};

/**
Time series: the same grid topology is used along with time series of injections
(productions and loads) to compute powerflows. Each step is warm started with the solution
of the previous one, which is what makes a *time* series: consecutive steps are assumed to
be close to one another, and the steps are therefore solved in order, on a single thread.
 **/
using TimeSeries = BaseInjectionSweep<BatchInitKind::FromPreviousStep>;

/**
Injection sweep: same inputs, same outputs and same interface as TimeSeries, but every step
restarts from the same seed (the caller's `Vinit`, or the converged "n" powerflow when
`init_from_n_powerflow` is set -- exactly like ContingencyAnalysis does for each of its
contingencies).

Use it when the "steps" are independent scenarios rather than consecutive instants: the
result of a step then does not depend on the steps computed before it, nor on the order in
which they were given, and the sweep can be split over several threads (see
`set_nb_thread`).
 **/
using InjectionSweep = BaseInjectionSweep<BatchInitKind::FromSeed>;

#ifndef LS2G_BUILDING_CORE
    // both instantiations are compiled once, into the core library (see the bottom of
    // BaseInjectionSweep.cpp); every other translation unit links against them instead of
    // instantiating the template again.
    extern template class LS2G_API BaseInjectionSweep<BatchInitKind::FromPreviousStep>;
    extern template class LS2G_API BaseInjectionSweep<BatchInitKind::FromSeed>;
#endif  // LS2G_BUILDING_CORE

} // namespace ls2g

#endif  //BASEINJECTIONSWEEP_H