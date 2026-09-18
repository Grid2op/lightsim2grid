// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SBUSPOLICY_H
#define SBUSPOLICY_H

#include "LSGrid.hpp"

#include <sstream>

namespace ls2g {

/**
Policy controlling how (if at all) a `BaseBatchSweep` instantiation's injection (Sbus)
varies from step to step. The two nested structs below are the only two states this axis
takes; see `BaseBatchSweep.hpp` for how they combine with `YbusPolicy` and
`BatchInitKind` into the four algorithms built on top of this template.

Deliberately free of any back-reference to the owning `BaseBatchSweep`: every method
takes the grid / solver-space state it needs as an explicit parameter, so this struct
can be embedded as a plain member with no CRTP-style cycle.
 **/
struct LS2G_API SbusPolicy
{
    /** the injection is fixed for the whole batch (ContingencyAnalysis): no state, no
        work -- the `BaseBatchSweep` per-step loop simply reuses the fixed member
        `ac_cache_.inj`/`dc_cache_.inj` for every row. **/
    struct NOOP {
        static constexpr bool supports_vary = false;
        // no-op: lets BaseBatchSweep::clear() (a plain, always-compiled member) call
        // sbus_policy_.clear() unconditionally regardless of instantiation.
        void clear() {}
    };

    /** each step has its own injection, built from the four per-step
        (generator / static generator / load) matrices handed to the `modify_*`
        setters (TimeSeries, InjectionSweep, ScenarioSweep). Holds its own per-step
        Sbus matrix, moved verbatim (de-templated off `BatchInitKind`, which never
        actually affected this code -- only the cosmetic `algo_name` text below) from
        the pre-refactor `BaseInjectionSweep<INIT>`.

        LS2G_API here (not just on the enclosing SbusPolicy): `assemble`/
        `constant_sbus_pu` are defined out-of-line in SbusPolicy.cpp (part of
        lightsim2grid_core) -- MSVC's dllexport does not propagate from an enclosing
        class to a nested one, so without this the symbols would be missing from
        lightsim2grid_core.dll's export table for any *other* DLL that ever calls
        into this struct directly (see the identical, and triggered, issue on
        YbusPolicy::Contingency). No current call site outside lightsim2grid_core
        happens to reach this one yet, but the codebase's own convention is to mark
        every core-lib class with out-of-line members this way regardless. **/
    struct LS2G_API Vary {
        static constexpr bool supports_vary = true;

        // RowMajor, matching BaseBatchSolverSynch::CplxMat -- NOT the column-major
        // ls2g::CplxMat from Utils.hpp. Nested here (rather than reusing the ambient
        // alias) so the two never collide and this struct stays embeddable anywhere.
        using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using CplxMat = Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        // one row per step, one column per solver-space bus -- built on demand by
        // materialize() for whoever asks for the whole matrix (get_sbuses()). The
        // row loop never reads it: it builds each row into its own buffer as it
        // goes (fill_row), so a year of hourly rows on a 9241-bus grid does not
        // hold 1.3 GB of injections that are read once. Mutable: a cache.
        mutable CplxMat sbuses;

        // per-step raw injection inputs, as set by the modify_gen_p/modify_sgen_p/
        // modify_load_p/modify_load_q setters (one row per step, one column per grid
        // element of the relevant type). Empty (0 rows) means "never set -- default
        // to the grid's own target value, broadcast across every row" (see assemble()
        // below). Row-count / cross-axis locking is NOT this struct's job: it is
        // shared with the Ybus-mask setters on ScenarioSweep, so `BaseBatchSweep`
        // owns that lock and only hands already-validated data here.
        RealMat gen_p, sgen_p, load_p, load_q;

        // per-step generator voltage-magnitude targets (target_vm_pu), one row per
        // step, one column per generator, as set by modify_gen_v. Empty (0 rows)
        // means "never set": the grid's own target_vm_pu_ was already correctly
        // seeded once, before the very first step, by
        // BaseBatchSolverSynch::_finish_preprocessing's unconditional set_vm() call
        // -- nothing to redo. Unlike gen_p/sgen_p/load_p/load_q above, this axis
        // feeds NOTHING into `sbuses`/assemble()/constant_sbus_pu(): it is not part
        // of the injection at all, it only re-seeds |V| at PV buses before each
        // step's solve (see BaseBatchSweep::_apply_step_gen_v). Kept here anyway
        // (rather than as a separate member of BaseBatchSweep) so it shares this
        // struct's existing row-count-lock-free / "0 rows = unset" contract and
        // stays reachable from wherever gen_p etc. are.
        RealMat gen_v;

        // per-step generator contingency mask, one row per step, one column per
        // generator: true means "this generator is disconnected for this step"
        // (ScenarioSweep only, see BaseBatchSweep::set_contingency_gens). Empty
        // (0 rows) means "never set", like every axis above.
        //
        // Two effects, and they live in two different places. The injection side
        // is here: assemble() takes the generator's active power -- and, for a NON
        // voltage-regulating one, its reactive setpoint -- back out of `sbuses`.
        // The labelling side is not: when the last locally voltage-regulating
        // generator of a bus is masked off, that bus turns PQ for the step, which
        // BaseBatchSweep drives through the algorithm's set_pv_pinned_buses.
        BoolMat gen_off;

        void clear() {
            gen_p = RealMat(); sgen_p = RealMat(); load_p = RealMat(); load_q = RealMat();
            gen_v = RealMat();
            gen_off = BoolMat();
            sbuses = CplxMat();
            gen_bus_.clear(); sgen_bus_.clear(); load_bus_.clear();
            gen_target_p_ = RealVect(); sgen_target_p_ = RealVect();
            load_target_p_ = RealVect(); load_target_q_ = RealVect();
            gen_target_q_ = RealVect(); gen_vreg_.clear();
            constant_pu_ = CplxVect();
            nb_buses_solver_ = 0;
            nb_steps_ = 0;
        }

        // Once per compute(): the routing every row needs -- the solver bus of each
        // element that injects (-1 for an inactive one), the grid's own targets for
        // the axes never set, the constant (storage/SVC/HVDC/...) share via
        // constant_sbus_pu(). The validation the per-element loops used to do
        // (an active element on a disconnected bus) is done here, once.
        void prepare(const LSGrid & grid_model,
                     bool ac_solver_used,
                     int nb_buses_solver,
                     const SolverBusIdVect & id_me_to_solver,
                     const Eigen::Ref<const CplxVect> & complete_sbus_pu,
                     real_type sn_mva,
                     Eigen::Index nb_steps,
                     const char * algo_name);

        // Row `i`'s injection (solver numbering, per unit) into `row`, resized as
        // needed: the same accumulation, in the same order and with the same
        // operations, as the whole-matrix build used to do -- generators, static
        // generators, loads (P then Q), the generators a row disconnects (P then
        // Q), the division by sn_mva, the constant share -- so a row here is the
        // row of that matrix, bit for bit.
        void fill_row(Eigen::Index i, CplxVect & row) const;

        // the whole matrix, for whoever wants it (get_sbuses()): every row through
        // fill_row, built at most once per prepare().
        const CplxMat & materialize() const;

        /**
         * Per-unit injection that the four per-step matrices (gen_p/sgen_p/load_p/
         * load_q) do NOT account for, and which is therefore constant across the
         * steps: storage units, REACTIVE_POWER-mode SVCs, the reactive setpoint of a
         * generator whose voltage regulation is off, a static generator's reactive
         * setpoint, the hvdc injections and (dc only) the phase-shifter term.
         *
         * Computed as `complete - accounted_for`: `complete_sbus_pu` is the complete
         * per-unit injection the gridmodel built (the caller's own ac_cache_.inj (ac) /
         * dc_cache_.inj.cast<cplx_type>() (dc), read-only here), minus the same
         * reconstruction the per-step matrices perform, evaluated at the gridmodel's
         * own target values. Deriving it this way -- instead of listing the elements
         * -- is what keeps it from falling out of date when a new element type starts
         * contributing to LSGrid::fillSbus_me.
         *
         * `algo_name` is only used to name the caller in an error message (this
         * struct no longer knows which alias -- TimeSeries / InjectionSweep /
         * ScenarioSweep -- it is embedded in).
         */
        CplxVect constant_sbus_pu(const LSGrid & grid_model,
                                  const Eigen::Ref<const CplxVect> & complete_sbus_pu,
                                  int nb_buses_solver,
                                  const SolverBusIdVect & id_me_to_solver,
                                  const char * algo_name) const;

        // ---- what prepare() leaves for fill_row() ---------------------------------
        std::vector<int> gen_bus_, sgen_bus_, load_bus_;   // solver bus per element, -1 if inactive
        RealVect gen_target_p_, sgen_target_p_, load_target_p_, load_target_q_;  // the grid's own
        RealVect gen_target_q_;          // a non-regulating generator's reactive setpoint
        std::vector<char> gen_vreg_;     // whether each generator regulates voltage
        CplxVect constant_pu_;
        real_type sn_mva_ = 1.;
        int nb_buses_solver_ = 0;
        Eigen::Index nb_steps_ = 0;

        // Adds (add=true) or subtracts (add=false) `temporal_data` (one row per step,
        // one column per element of `structure_data`) into `Sbuses` (one row per
        // step, one column per solver-space bus), routed through
        // `id_me_to_ac_solver`. Static: touches no Vary state, purely a function of
        // its parameters. Moved verbatim from the pre-refactor
        // BaseInjectionSweep::fill_SBus_real/fill_SBus_imag.
        template<class T>
        static void fill_SBus_real(Eigen::Ref<CplxMat> Sbuses,
                                   const T & structure_data,
                                   const Eigen::Ref<const RealMat> & temporal_data,
                                   const SolverBusIdVect & id_me_to_ac_solver,
                                   bool add,  // if true call += else calls -=
                                   const char * algo_name)
        {
            size_t nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            SolverBusId bus_id_solver;
            GlobalBusId bus_id_me;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                if(bus_id_me.cast_int() == BaseConstants::_deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name << "::fill_SBus_real: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                bus_id_solver = id_me_to_ac_solver[bus_id_me.cast_int()];
                if(bus_id_solver.cast_int() == BaseConstants::_deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name << "::fill_SBus_real: the element with id ";
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
        static void fill_SBus_imag(Eigen::Ref<CplxMat> Sbuses,
                                   const T & structure_data,
                                   const Eigen::Ref<const RealMat> & temporal_data,
                                   const SolverBusIdVect & id_me_to_ac_solver,
                                   bool add,  // if true call += else calls -=
                                   const char * algo_name)
        {
            size_t nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            SolverBusId  bus_id_solver;
            GlobalBusId bus_id_me;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                if(bus_id_me.cast_int() == BaseConstants::_deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name << "::fill_SBus_imag: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                bus_id_solver = id_me_to_ac_solver[bus_id_me.cast_int()];
                if(bus_id_solver.cast_int() == BaseConstants::_deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name << "::fill_SBus_imag: the element with id ";
                    exc_ << el_id;
                    exc_ << " is connected to a disconnected bus while being connected";
                    throw std::runtime_error(exc_.str());
                }
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver.cast_int()) += BaseConstants::my_i * tmp;
                else Sbuses.col(bus_id_solver.cast_int()) -= BaseConstants::my_i * tmp;
            }
        }
    };
};

} // namespace ls2g

#endif  // SBUSPOLICY_H
