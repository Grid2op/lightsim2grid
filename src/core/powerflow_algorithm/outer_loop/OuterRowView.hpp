// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_ROW_VIEW_H
#define OUTER_ROW_VIEW_H

#include <array>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "Utils.hpp"
#include "TaggedIdVec.hpp"

namespace ls2g {

class LSGrid;
class BaseOuterLoop;
struct SolverBusLayout;

/**
 * A solve's outer loops as they ended, for other solves to resume from: a batch's "n" solve,
 * which every row seeded from it starts from (see BaseAlgo::get_outer_resume_point and
 * OuterRowView::resume_from). Immutable once built, and shared read-only by every worker of
 * the batch, each of which takes clones of its loops.
 */
struct OuterResumePoint
{
    /// a clone of each loop, in the algorithm's order, with the state the solve left it in
    /// (the buses it froze, the SVCs it switched on, the hvdc regimes, ...)
    std::vector<std::shared_ptr<const BaseOuterLoop> > loops;
    /// the positions the loops settled on, by element id, INT_MIN where none handled it
    std::vector<int> phase_tap;
    std::vector<int> ratio_tap;
    std::vector<int> sections;
    /// the active targets a loop moved, MW, NaN elsewhere (OuterInjections' convention)
    std::vector<real_type> gen_target_p;
    std::vector<real_type> storage_target_p;
    /// the solve's voltages, solver numbering: what a resuming loop reads of the solve before
    CplxVect V;
    RealVect Va;
    RealVect Vm;
};

/**
 * What one row of a batch changes compared with the grid, as the outer loops read it.
 *
 * A batch (TimeSeries, InjectionSweep, ContingencyAnalysis, ScenarioSweep) never modifies its
 * LSGrid: a row is the grid plus a few overrides -- a generator or a branch disconnected, a
 * generator placed on another bus, other targets. A loop sees the row through this, so it
 * acts (and detects) on the row's own state rather than on the base grid's.
 *
 * Plain data: one container per element kind, dense and indexed by element id, EMPTY when
 * the row changes nothing of that kind (the usual case, and then nothing is allocated). The
 * owner (a batch worker) refills it before each row. `layout` and `masked` point to the
 * owner's own data, which outlives the solve.
 *
 * Read it through OuterGridView, which falls back on the grid wherever the row says nothing.
 */
struct OuterRowView
{
    /// the solver labelling of the solve: id_me_to_solver, the slack, the pv / pq split and
    /// the voltage-control plan (null: the grid's AC one)
    const SolverBusLayout * layout = nullptr;
    /// the solver buses this row leaves out (a contingency strands them), sorted (null: none)
    const std::vector<int> * masked = nullptr;
    /// the loops start from where the previous row left them (a TimeSeries row: its equations
    /// are the previous row's), not from the grid (see OuterLoopDriverParams::resume_outer_state)
    bool resume = false;
    /// with `resume`, the solve they start from instead of the algorithm's own last one: the
    /// batch's "n" solve, for a row seeded from it (init_from_n_powerflow). Null: the last one.
    const OuterResumePoint * resume_from = nullptr;
    /// the active power the row's contingency takes out of the balance (MW, generator
    /// convention: lost generation minus lost load), at the units' row or grid targets: the
    /// DistributedSlack loop shares it on its participants before the first solve, as
    /// OpenLoadFlow's security analysis does (a batch's `redistribute_slack`). NaN: nothing to
    /// share.
    real_type lost_p_mw = std::numeric_limits<real_type>::quiet_NaN();
    /// detection only: the pi blocks (y11, y12, y21, y22, pu) of the transformers the solve's
    /// outer loops left at another tap, by increasing transformer id (TrafoContainer::
    /// moved_tap_blocks): the state that solve ended in (null: none moved, or not a detection)
    const std::vector<std::pair<int, std::array<cplx_type, 4> > > * moved_trafo_blocks = nullptr;

    /// One element kind's overrides. Each vector is empty, or sized by the element count.
    struct Units
    {
        static constexpr int GRID_BUS = -2;  ///< bus_solver: where the grid has it
        /// 1: this row disconnects it
        std::vector<char> off;
        /// GRID_BUS, or the solver bus this row places it on (a topological action, which may
        /// reconnect an element the grid has off)
        std::vector<int> bus_solver;
        /// NaN where the row keeps the grid's: active target (MW, the container's convention),
        /// voltage target (pu), reactive target (MVar, the container's convention)
        std::vector<real_type> target_p;
        std::vector<real_type> target_vm;
        std::vector<real_type> target_q;

        bool empty() const {
            return off.empty() && bus_solver.empty() && target_p.empty() && target_vm.empty() && target_q.empty();
        }
        void clear() { off.clear(); bus_solver.clear(); target_p.clear(); target_vm.clear(); target_q.clear(); }
    };
    /// One branch kind's overrides. Each vector is empty, or sized by the branch count.
    struct Branches
    {
        /// 1: this row opens the branch (a contingency)
        std::vector<char> off;
        /// per side, Units::GRID_BUS, -1 (this row disconnects that side) or the solver bus
        /// this row moves it to (a topological action)
        std::vector<int> bus1_solver;
        std::vector<int> bus2_solver;

        bool empty() const { return off.empty() && bus1_solver.empty() && bus2_solver.empty(); }
        void clear() { off.clear(); bus1_solver.clear(); bus2_solver.clear(); }
    };
    Units gens;
    Units storages;
    Units loads;
    Branches lines;
    Branches trafos;

    void clear_overrides() {
        gens.clear();
        storages.clear();
        loads.clear();
        lines.clear();
        trafos.clear();
    }
    /// lay the active targets a solve's outer loops moved (OuterInjections' convention: NaN
    /// where a loop kept the unit's) over this row's: the state that solve ended in, which a
    /// detection after it reads (a unit's reactive limits follow its target P)
    void take_moved_targets(const std::vector<real_type> & gen_target_p,
                            const std::vector<real_type> & storage_target_p);
};

/**
 * The grid as one row of a batch sees it: the grid's data, with the row's overrides on top
 * (none for a single solve, where `row` is null). Everything an outer loop or a trigger reads
 * that a batch row may change goes through here; what no row changes (limits, curves, tap
 * tables, nominal voltages) is read on the grid directly.
 *
 * Cheap to copy: a reference and a pointer.
 */
class LS2G_API OuterGridView
{
    public:
        OuterGridView(const LSGrid & grid, const OuterRowView * row) : grid_(grid), row_(row) {}

        const LSGrid & grid() const { return grid_; }
        const OuterRowView * row() const { return row_; }

        /// the labelling of the solve (the row's, the grid's AC one otherwise)
        const SolverBusLayout & layout() const;
        const SolverBusIdVect & id_me_to_solver() const;
        const GlobalBusIdVect & id_solver_to_me() const;
        /// the row's stranded solver buses, sorted; null: none
        const std::vector<int> * masked() const { return row_ != nullptr ? row_->masked : nullptr; }
        bool is_masked(int bus_solver) const;

        // ----- generators ----------------------------------------------------------------
        /// connected in this row
        bool gen_on(int gen_id) const;
        /// its solver bus in this row, -1 when it is off or its bus is not solved
        int gen_bus_solver(int gen_id) const;
        /// the grid bus it stands on in this row (a placed generator's), -1 when off
        int gen_bus_grid(int gen_id) const;
        /// connected, regulating its own bus through the PV path, in this row
        bool gen_is_local_voltage_controller(int gen_id) const;
        real_type gen_target_p(int gen_id) const;   ///< MW
        /// whether the row sets the unit's active target (else it is the grid's)
        bool row_sets_gen_target_p(int gen_id) const;
        real_type gen_target_vm(int gen_id) const;  ///< pu

        // ----- storage units (load convention for P, as the container) ------------------
        bool storage_on(int storage_id) const;
        int storage_bus_solver(int storage_id) const;
        bool storage_is_local_voltage_controller(int storage_id) const;
        real_type storage_target_p(int storage_id) const;
        bool row_sets_storage_target_p(int storage_id) const;
        real_type storage_target_vm(int storage_id) const;

        // ----- loads (load convention) --------------------------------------------------
        bool load_on(int load_id) const;
        int load_bus_solver(int load_id) const;
        real_type load_target_q(int load_id) const;

        // ----- branches -----------------------------------------------------------------
        /// whether the branch is closed at both ends in this row, and then the solver buses
        /// of its two ends (both solved, possibly equal)
        bool line_closed(int line_id, int & bus1_solver, int & bus2_solver) const;
        bool trafo_closed(int trafo_id, int & bus1_solver, int & bus2_solver) const;
        /// its pi block (y11, y12, y21, y22, pu): at the tap the solve's loops left it at (see
        /// OuterRowView::moved_trafo_blocks), the one the grid stamps otherwise
        std::array<cplx_type, 4> trafo_block(int trafo_id) const;

        /// a generator's target P where a loop or the row moved it, NaN where neither did: what
        /// LSGrid::gen_limits_at_outer_target reads. `loops` (per generator, may be empty) wins
        /// over the row.
        real_type gen_target_p_override(const std::vector<real_type> & loops, int gen_id) const;

    private:
        static int _placed(const OuterRowView::Units * u, int id) {
            if (u == nullptr || id < 0 || static_cast<std::size_t>(id) >= u->bus_solver.size()) return OuterRowView::Units::GRID_BUS;
            return u->bus_solver[static_cast<std::size_t>(id)];
        }
        static bool _value(const std::vector<real_type> & v, int id, real_type & res);

        const LSGrid & grid_;
        const OuterRowView * row_;
};

}  // namespace ls2g

#endif  // OUTER_ROW_VIEW_H
