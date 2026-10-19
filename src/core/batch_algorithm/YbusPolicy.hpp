// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef YBUSPOLICY_H
#define YBUSPOLICY_H

#include "LSGrid.hpp"

#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace ls2g {

/** What one contingency of a ContingencyAnalysis trips: branches (lines then
    transformers, the gridmodel numbering every batch class uses), shunts, loads and
    generators, by their own ids. Ordered field by field, branches first, so a set of
    branch-only contingencies keeps the order a set of branch sets had. **/
struct LS2G_API ContingencyElements {
    std::set<int> branches;
    std::set<int> shunts;
    std::set<int> loads;
    std::set<int> gens;

    bool operator<(const ContingencyElements & other) const {
        return std::tie(branches, shunts, loads, gens) < std::tie(other.branches, other.shunts, other.loads, other.gens);
    }
    bool operator==(const ContingencyElements & other) const {
        return std::tie(branches, shunts, loads, gens) == std::tie(other.branches, other.shunts, other.loads, other.gens);
    }
    bool empty() const { return branches.empty() && shunts.empty() && loads.empty() && gens.empty(); }
};

/**
Policy controlling how (if at all) a `BaseBatchSweep` instantiation's admittance matrix
varies from step to step. The two nested structs below are the only two states this
axis takes; see `BaseBatchSweep.hpp` for how they combine with `SbusPolicy` and
`BatchInitKind` into the four algorithms built on top of this template.

Deliberately free of any back-reference to the owning `BaseBatchSweep`: every method
takes the grid / solver-space state it needs as an explicit parameter, so this struct
can be embedded as a plain member with no CRTP-style cycle.
 **/
struct LS2G_API YbusPolicy
{
    /** the topology is fixed for the whole batch (TimeSeries, InjectionSweep): no
        state, no work -- the `BaseBatchSweep` per-step loop's before_step/after_step
        hooks are empty for this policy and the compiler elides them entirely. **/
    struct NOOP {
        static constexpr bool supports_contingency = false;
    };

    /** each step disconnects a set of lines/trafos (branch ids, numbered lines-then-
        trafos), emulated by editing the admittance matrix coefficients rather than
        rebuilding it from scratch (used by ContingencyAnalysis and ScenarioSweep).
        Holds its own Coeff-list bookkeeping, moved verbatim from the pre-refactor
        `ContingencyAnalysis`.

        LS2G_API here (not just on the enclosing YbusPolicy) is required, not
        decorative: this nested struct's methods are defined out-of-line in
        YbusPolicy.cpp (part of lightsim2grid_core) and called directly from
        BaseBatchSweep's inline (header-defined) is_grid_connected_after_contingency()
        / pick_reference_slack(), which get compiled into lightsim2grid_cpp.pyd -- a
        *different* DLL on Windows. MSVC's dllexport does NOT propagate from an
        enclosing class to a nested one, so without this the symbols are simply
        missing from lightsim2grid_core.dll's export table (LNK2001 "unresolved
        external symbol" when linking lightsim2grid_cpp.pyd) even though every
        GCC/Clang build (default symbol visibility is far more permissive) links
        fine regardless. **/
    struct LS2G_API Contingency {
        static constexpr bool supports_contingency = true;

        // registered contingencies (ContingencyAnalysis): what each one trips, see
        // ContingencyElements. A std::set, NOT unordered: callers (my_defaults_vect(),
        // clean_flows()) rely on the ordering, which is the row order. Mutated directly
        // by BaseBatchSweep's add_n1/add_nk/add_contingency/remove_*/... (those stay on
        // BaseBatchSweep: they need the element counts for bounds-checking, which this
        // policy deliberately does not know about).
        std::set<ContingencyElements> li_defaults;
        // li_defaults in row order (ContingencyAnalysis), refreshed with the
        // coefficients by init_li_coeffs: what row i trips, without walking the set
        std::vector<ContingencyElements> li_rows;
        // for each entry of li_defaults (ContingencyAnalysis) -- or each row of
        // line_mask/trafo_mask (ScenarioSweep) -- the Ybus coefficients to remove/
        // re-add to emulate it, in the same order. Moved verbatim (as a container)
        // from ContingencyAnalysis::_li_coeffs; ScenarioSweep reuses the same field
        // via init_li_coeffs_from_masks below instead of init_li_coeffs, li_defaults
        // staying empty/unused on that path. A shunt the entry trips is one diagonal
        // coefficient (AC: it is an admittance; nothing in DC, see li_sbus_dp).
        std::vector<std::vector<Coeff> > li_coeffs;
        // aligned with li_coeffs: what to ADD to the entry's injection (solver bus id,
        // per unit), what the matrix edit does not carry:
        //  - DC, a phase shifter the entry disconnects: it is not in the matrix, its
        //    shift is a pair of injections at its two buses
        //    (TrafoContainer::hack_Sbus_for_dc_phase_shifter), so taking its
        //    coefficients out of Bbus leaves them behind -- the opposite pair, here;
        //  - DC, a shunt the entry trips: its active part, which DC puts in Pbus;
        //  - a load or a generator the entry trips (ContingencyAnalysis, whose
        //    injection is fixed; ScenarioSweep takes them out of its own per-row
        //    injection, see SbusPolicy::Vary): the opposite of what it stamps, its
        //    active power and, in AC, its reactive power (a generator's only when it
        //    does not regulate voltage, as GeneratorContainer stamps it).
        // Empty for an entry that needs none of that.
        std::vector<std::vector<std::pair<int, cplx_type> > > li_sbus_dp;

        // dense row-aligned contingency masks (ScenarioSweep only): row i, column j
        // True means "deactivate this branch for step i". Empty (0 rows) means "never
        // set" -- treated as all-False. Columns: line_mask has one per powerline,
        // trafo_mask one per trafo; branch ids are numbered lines-then-trafos, same
        // convention as everywhere else in this file.
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> line_mask;
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trafo_mask;
        // same, one column per shunt (ScenarioSweep set_contingency_shunts): a shunt is
        // an admittance, so it lives with the branches, as a diagonal coefficient
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> shunt_mask;

        // per row, the branches (gridmodel numbering, lines then trafos, sorted) that
        // the row's topological action disconnects on top of the masks (ScenarioSweep
        // set_topo_actions; see BaseBatchSweep::_maybe_resolve_topology, which fills it
        // once per compute() and refuses a branch named by both a mask and the action).
        // Empty (0 rows) means "no topological action". Read wherever the masks are:
        // branch_ids_for_row merges the two.
        std::vector<std::vector<int> > topo_branches_off;

        // per row, the branches the row's action places itself: an end moved to another
        // busbar, or a branch disconnected in the base grid reconnected (ScenarioSweep
        // set_topo_actions). Solver numbering. `base_on`: whether the branch is
        // connected in the base grid (its base contribution is then taken out first).
        // The row's contribution goes in with the branch's RAW block (yac_11..22 /
        // ydc_*: both ends closed), at the row's buses -- entries the base pattern
        // may not hold, reserved up front by BaseBatchSweep::_maybe_reserve_union_pattern.
        struct BranchPlacement {
            int branch_id;
            int bus1_solver;
            int bus2_solver;
            bool base_on;
        };
        std::vector<std::vector<BranchPlacement> > topo_branches_moved;

        // builds li_coeffs from li_defaults (see the pre-refactor
        // ContingencyAnalysis::init_li_coeffs). `grid_model`/`n_line` supply what used
        // to come from the owning class's _grid_model/n_line_ members.
        void init_li_coeffs(const LSGrid & grid_model,
                            bool ac_solver_used,
                            const SolverBusIdVect & id_me_to_solver,
                            size_t n_line);

        // builds li_coeffs from line_mask/trafo_mask instead of li_defaults: one entry
        // per row of `nb_steps`, in row order (ScenarioSweep's per-scenario mask is
        // row-aligned with its injections, unlike ContingencyAnalysis's set of
        // distinct scenarios). A row with both masks empty/all-False produces an
        // empty coefficient list for that step (no branch touched).
        void init_li_coeffs_from_masks(const LSGrid & grid_model,
                                       bool ac_solver_used,
                                       const SolverBusIdVect & id_me_to_solver,
                                       size_t n_line,
                                       Eigen::Index nb_steps);

        // per-row branch ids (gridmodel numbering, lines then trafos) implied by
        // line_mask/trafo_mask row `row` -- the ScenarioSweep counterpart of
        // ContingencyAnalysis's my_defaults_vect()-derived skip list (used by
        // BaseBatchSweep::_row_skip_branch_ids to exclude a row's own disconnected
        // branches from that row's current-limit checks), read directly off the row
        // instead of a precomputed cache, merged with topo_branches_off[row]. Empty
        // if neither mask was ever set and no action disconnects a branch on that row.
        std::vector<int> branch_ids_for_row(Eigen::Index row, size_t n_line) const;
        // per-row shunt ids implied by shunt_mask row `row` (sorted, empty if never set)
        std::vector<int> shunt_ids_for_row(Eigen::Index row) const;
        // whether row `row` trips that shunt / load / generator: its contingency
        // (ContingencyAnalysis, li_rows) or its shunt mask (ScenarioSweep; its loads and
        // generators are SbusPolicy::Vary's)
        bool shunt_off_in_row(size_t row, int shunt_id) const {
            if(row < li_rows.size()) return li_rows[row].shunts.count(shunt_id) > 0;
            return shunt_mask.rows() > 0 && static_cast<Eigen::Index>(row) < shunt_mask.rows() &&
                   shunt_id < shunt_mask.cols() && shunt_mask(static_cast<Eigen::Index>(row), shunt_id);
        }
        bool load_off_in_row(size_t row, int load_id) const {
            return row < li_rows.size() && li_rows[row].loads.count(load_id) > 0;
        }
        bool gen_off_in_row(size_t row, int gen_id) const {
            return row < li_rows.size() && li_rows[row].gens.count(gen_id) > 0;
        }
        // whether any contingency trips a generator (ContingencyAnalysis)
        bool has_gen_off() const {
            for(const ContingencyElements & el : li_defaults) if(!el.gens.empty()) return true;
            return false;
        }

        // remove / re-add one contingency's coefficients from/to Ybus (AC) or the
        // algo's own internal Ybus (DC). Whether the grid stays connected after the
        // removal is NOT decided here: BaseBatchSweep settles it once per
        // contingency, before the row loop, from one depth-first search of the base
        // graph (see BusGraph and BaseBatchSweep::_prepare_connectivity) -- a
        // search per row, on a matrix that is the base one but for four entries, was
        // what these two calls used to cost. Neither reads nor writes any
        // Contingency state, so both are static.
        static void remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const std::vector<Coeff> & coeffs,
                                     bool ac_solver_used,
                                     AlgorithmSelector & algo);
        static void readd_to_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                  const std::vector<Coeff> & coeffs,
                                  bool ac_solver_used,
                                  AlgorithmSelector & algo);

        private:
            // the row's own contribution of a branch placed by its action, as
            // coefficients to REMOVE (negated: remove_from_Ybus's -= adds them)
            static void _append_placement_coeffs(const BranchPlacement & placement,
                                                 const LSGrid & grid_model,
                                                 bool ac_solver_used,
                                                 size_t n_line,
                                                 std::vector<Coeff> & out);
            // what tripping these shunts / loads / generators adds to an entry: a
            // diagonal coefficient per shunt (AC) into `coeffs`, and into `sbus_dp`
            // the injection terms li_sbus_dp describes. An element off in the base grid
            // is not in the matrix nor in the injection: nothing to take out.
            static void _append_shunt_trips(const std::vector<int> & shunt_ids,
                                            const LSGrid & grid_model,
                                            bool ac_solver_used,
                                            const SolverBusIdVect & id_me_to_solver,
                                            std::vector<Coeff> & coeffs,
                                            std::vector<std::pair<int, cplx_type> > & sbus_dp);
            static void _append_load_trips(const std::vector<int> & load_ids,
                                           const LSGrid & grid_model,
                                           bool ac_solver_used,
                                           const SolverBusIdVect & id_me_to_solver,
                                           std::vector<std::pair<int, cplx_type> > & sbus_dp);
            static void _append_gen_trips(const std::vector<int> & gen_ids,
                                          const LSGrid & grid_model,
                                          bool ac_solver_used,
                                          const SolverBusIdVect & id_me_to_solver,
                                          std::vector<std::pair<int, cplx_type> > & sbus_dp);
            // shared by init_li_coeffs and init_li_coeffs_from_masks: the Ybus
            // coefficients (up to 4 per branch) that emulate disconnecting exactly
            // `branch_ids` (gridmodel numbering, lines then trafos), and -- in DC -- the
            // phase-shift injections that go with them (see li_sbus_dp), into
            // `dc_shift_dp` (overwritten).
            static std::vector<Coeff> _coeffs_for_branch_ids(const std::vector<int> & branch_ids,
                                                              const LSGrid & grid_model,
                                                              bool ac_solver_used,
                                                              const SolverBusIdVect & id_me_to_solver,
                                                              size_t n_line,
                                                              std::vector<std::pair<int, cplx_type> > & dc_shift_dp);
    };
};

} // namespace ls2g

#endif  // YBUSPOLICY_H
