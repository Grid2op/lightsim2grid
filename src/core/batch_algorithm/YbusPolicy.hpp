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
#include <vector>

namespace ls2g {

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

        // registered contingencies, gridmodel branch ids (lines then trafos). A
        // std::set<std::set<int>>, NOT unordered: callers (my_defaults_vect(),
        // clean_flows()) rely on the ordering. Moved verbatim from
        // ContingencyAnalysis::_li_defaults. Mutated directly by BaseBatchSweep's
        // add_n1/add_nk/remove_n1/... (those stay on BaseBatchSweep: they need
        // n_total_ from BaseBatchSolverSynch for bounds-checking, which this policy
        // deliberately does not know about).
        std::set<std::set<int> > li_defaults;
        // for each entry of li_defaults (ContingencyAnalysis) -- or each row of
        // line_mask/trafo_mask (ScenarioSweep) -- the Ybus coefficients to remove/
        // re-add to emulate it, in the same order. Moved verbatim (as a container)
        // from ContingencyAnalysis::_li_coeffs; ScenarioSweep reuses the same field
        // via init_li_coeffs_from_masks below instead of init_li_coeffs, li_defaults
        // staying empty/unused on that path.
        std::vector<std::vector<Coeff> > li_coeffs;

        // dense row-aligned contingency masks (ScenarioSweep only): row i, column j
        // True means "deactivate this branch for step i". Empty (0 rows) means "never
        // set" -- treated as all-False. Columns: line_mask has one per powerline,
        // trafo_mask one per trafo; branch ids are numbered lines-then-trafos, same
        // convention as everywhere else in this file.
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> line_mask;
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trafo_mask;

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
                                       size_t n_trafo,
                                       Eigen::Index nb_steps);

        // remove / re-add one contingency's coefficients from/to Ybus (AC) or the
        // algo's own internal Ybus (DC). Returns (remove_from_Ybus only) whether the
        // grid stays connected after removal -- see check_invertible below. Neither
        // reads nor writes any Contingency state, so both are static (moved verbatim
        // from the pre-refactor ContingencyAnalysis::remove_from_Ybus/readd_to_Ybus).
        static bool remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const std::vector<Coeff> & coeffs,
                                     bool ac_solver_used,
                                     AlgorithmSelector & algo);
        static void readd_to_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                  const std::vector<Coeff> & coeffs,
                                  bool ac_solver_used,
                                  AlgorithmSelector & algo);

        // BFS connectivity check on Ybus's sparsity pattern (see the pre-refactor
        // ContingencyAnalysis::check_invertible). Static: uses no Contingency state.
        static bool check_invertible(const Eigen::Ref<const Eigen::SparseMatrix<cplx_type> > & Ybus);

        private:
            // shared by init_li_coeffs and init_li_coeffs_from_masks: the Ybus
            // coefficients (up to 4 per branch) that emulate disconnecting exactly
            // `branch_ids` (gridmodel numbering, lines then trafos).
            static std::vector<Coeff> _coeffs_for_branch_ids(const std::vector<int> & branch_ids,
                                                              const LSGrid & grid_model,
                                                              bool ac_solver_used,
                                                              const SolverBusIdVect & id_me_to_solver,
                                                              size_t n_line);
    };
};

} // namespace ls2g

#endif  // YBUSPOLICY_H
