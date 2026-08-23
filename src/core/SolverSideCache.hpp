// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_SOLVER_SIDE_CACHE_H
#define LS2G_SOLVER_SIDE_CACHE_H

#include <type_traits>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Sparse"

#include "TaggedIdVec.hpp"
#include "Utils.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

// ---------------------------------------------------------------------------
// WHY A SOLVER-SIDE CACHE IS ONE OBJECT
//
// The bus labelling, the matrix, the injections, the slack and the pv-pq split
// are not independent things. They are one picture of the grid taken at one
// instant, and all of it is expressed in ONE bus labelling -- the compact solver
// numbering `id_me_to_solver` defines. Mix two instants, or two labellings, and
// nothing complains: the sizes still match (the bus count rarely changes), the
// powerflow still converges, and the answer is quietly wrong. No assertion can
// catch it after the fact, because every number is individually plausible.
//
// They used to be nine (then twelve, then eighteen) separate members of LSGrid,
// built by a function that took six of them as parameters and reached for the
// other three itself. That is how a batch algorithm building into its own
// vectors ended up writing its pv-pq split into the grid's cache, next to the
// grid's own matrix. One object makes that class of mistake unspeakable rather
// than merely unlikely: you cannot hand a function half of a struct.
//
// The split below is by TYPE, not by lifetime: SolverBusLayout is the part whose
// type does not mention the family's scalar, SolverSideCache<MatScalar> adds the
// two that do. Both halves are per family and are built, reused and retired as
// one unit.
// ---------------------------------------------------------------------------

/**
 * The half of a solver-side cache whose TYPE does not depend on the family: the
 * bus labelling, the slack, the pv-pq split, and the bookkeeping that says
 * whether any of it is still valid.
 *
 * Split out from SolverSideCache so that code which is generic over the family
 * -- the batch algorithms, which run AC or DC but never both in one sweep -- can
 * name and pass it without knowing which of Ybus/Bbus sits next to it. The VALUES
 * are still per family: an AC and a DC cache each own a full one of these, they
 * are never shared. That is the whole point; sharing them across families is what
 * used to make a "nothing changed" claim unsound.
 *
 * New per-family state that is not a matrix or an injection vector (remote /
 * shared voltage-control layouts, HVDC droop data, ...) belongs here.
 */
struct SolverBusLayout
{
    // ---- the bus labelling everything is expressed in ---------------------
    /// grid ("me") bus id -> solver bus id; sized by the grid's total bus count,
    /// `_deactivated_bus_id` wherever a bus carries nothing
    SolverBusIdVect id_me_to_solver;
    /// solver bus id -> grid bus id; sized by the number of buses actually solved
    GlobalBusIdVect id_solver_to_me;

    // ---- the slack ---------------------------------------------------------
    GlobalBusIdVect slack_bus_id_me;      ///< slack buses, grid numbering
    SolverBusIdVect slack_bus_id_solver;  ///< the same, solver numbering
    RealVect slack_weights;               ///< distributed-slack share per solver bus

    // ---- the pv / pq split -------------------------------------------------
    SolverBusIdVect bus_pv;  ///< solver ids, NOT grid ids
    SolverBusIdVect bus_pq;  ///< solver ids, NOT grid ids

    // ---- what makes the above reusable, or not -----------------------------
    /**
     * The grid's bus connectivity at the moment this cache was built. A powerflow
     * compares the grid's current connectivity against this to decide whether the
     * labelling above still describes it; each family keeps its own, because a
     * family that has not solved in a while may be several grid modifications
     * behind the other one. Empty means "nothing built" -- which is also how a
     * cache is retired without throwing its contents away (see the publication
     * step at the end of LSGrid::_pre_process_solver_impl).
     */
    std::vector<bool> last_bus_status;
    /**
     * Durable "never reuse this" switch (LSGrid::allow_ac_cache_reuse etc.). The
     * answer is identical either way, so this is a debugging switch and an escape
     * hatch for code that mutates the grid behind LSGrid's back.
     */
    bool allow_reuse = true;
    /**
     * Set when this family's previous powerflow diverged and its algorithm was
     * reset: the DATA is still a correct picture of the grid (divergence is a
     * numerical failure, not a data one) but the algorithm must rebuild its own
     * internals from it rather than skip on a "nothing changed" it cannot honour.
     * Cleared once the algorithm has actually run.
     */
    bool algo_needs_rebuild = false;

    /// number of buses in the solved system (0 when nothing was ever built)
    [[nodiscard]] Eigen::Index nb_bus_solver() const noexcept {
        return static_cast<Eigen::Index>(id_solver_to_me.size());
    }

    /// Back to "this family has never solved", for everything in this half.
    void clear_layout(){
        id_me_to_solver = SolverBusIdVect();
        id_solver_to_me = GlobalBusIdVect();
        slack_bus_id_me = GlobalBusIdVect();
        slack_bus_id_solver = SolverBusIdVect();
        slack_weights = RealVect();
        bus_pv = SolverBusIdVect();
        bus_pq = SolverBusIdVect();
        last_bus_status.clear();
    }

    /// the family-agnostic half of is_usable(); see there for why this exists
    [[nodiscard]] bool layout_is_usable(std::size_t nb_bus_grid) const noexcept {
        const Eigen::Index nb_solver = nb_bus_solver();
        if(!allow_reuse) return false;              // told never to reuse
        if(nb_solver == 0) return false;            // never built
        if(id_me_to_solver.size() != nb_bus_grid) return false;  // built for another grid size
        if(slack_weights.size() != nb_solver) return false;
        // every bus is pv, pq, or slack: the split can never outnumber the system
        if(bus_pv.size() + bus_pq.size() > static_cast<std::size_t>(nb_solver)) return false;
        // a snapshot that does not cover the grid means "retired", see last_bus_status
        if(last_bus_status.size() != nb_bus_grid) return false;
        return true;
    }
};

/**
 * One solver family's complete cache: the layout above, plus the two things whose
 * type depends on the family.
 *
 * The template parameter picks the family: `cplx_type` is AC (complex admittance
 * Ybus, complex injection Sbus), `real_type` is DC (real susceptance Bbus, real
 * active injection Pbus). Everything else is inherited, which is why the two
 * families share one definition instead of two hand-kept-in-sync copies.
 */
template<class MatScalar>
struct SolverSideCache final : SolverBusLayout
{
    static_assert(std::is_same<MatScalar, cplx_type>::value ||
                  std::is_same<MatScalar, real_type>::value,
                  "SolverSideCache: MatScalar must be cplx_type (AC) or real_type (DC)");

    /// complex Sbus for the AC family, real Pbus for the DC one
    using InjVect = typename std::conditional<std::is_same<MatScalar, cplx_type>::value,
                                              CplxVect, RealVect>::type;
    using Matrix = Eigen::SparseMatrix<MatScalar>;

    /// true for the AC family. A compile-time constant, usable wherever the old
    /// runtime `is_ac` flag was threaded through by hand.
    static constexpr bool is_ac = std::is_same<MatScalar, cplx_type>::value;

    /// Ybus (AC, complex) / Bbus (DC, real), square, nb_bus_solver()
    Matrix mat;
    /// Sbus (AC, complex) / Pbus (DC, real), per unit, nb_bus_solver()
    InjVect inj;

    /// Back to "this family has never solved". Everything that describes a built
    /// system goes, the two policy flags stay: `allow_reuse` is the caller's
    /// standing choice, not derived state, and `algo_needs_rebuild` is about the
    /// algorithm, which a caller resets separately when it means to.
    void clear(){
        clear_layout();
        mat = Matrix();
        inj = InjVect();
    }

    /**
     * Is there really a system in here, of the shape the caller is about to claim
     * is up to date, for a grid with `nb_bus_grid` buses?
     *
     * This is the one guard between a "nothing changed" claim and an out-of-bounds
     * read. A powerflow's change flags only record what happened SINCE that family
     * last solved; they cannot say whether it ever solved at all. Nothing stops a
     * caller from asserting "nothing changed" -- `LSGrid::unset_changes()` does
     * exactly that, for both families at once -- about a family whose data is still
     * default-constructed. The "nothing to rebuild" path was then taken with an
     * empty `id_me_to_solver` / `mat` / `inj`, and everything downstream indexes
     * them with bus ids in the hundreds. Release wheels are -O3 -DNDEBUG, so that
     * is a segfault, not an error.
     *
     * So the flags are never trusted alone: check that the data they describe is
     * really there, and rebuild when it is not. A wrong "nothing changed" then
     * costs time, never memory safety. A handful of integer comparisons, off any
     * inner loop -- every size here is O(1) to read.
     */
    [[nodiscard]] bool is_usable(std::size_t nb_bus_grid) const noexcept {
        if(!layout_is_usable(nb_bus_grid)) return false;
        const Eigen::Index nb_solver = nb_bus_solver();
        if(mat.rows() != nb_solver) return false;
        if(mat.cols() != mat.rows()) return false;
        if(inj.size() != nb_solver) return false;
        return true;
    }
};

/// the AC family: complex Ybus / Sbus
using AcSolverCache = SolverSideCache<cplx_type>;
/// the DC family: real Bbus / Pbus
using DcSolverCache = SolverSideCache<real_type>;

} // namespace ls2g

#endif // LS2G_SOLVER_SIDE_CACHE_H
