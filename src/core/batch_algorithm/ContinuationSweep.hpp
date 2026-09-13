// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CONTINUATIONSWEEP_H
#define CONTINUATIONSWEEP_H

#include "BaseBatchSolverSynch.hpp"
#include "SbusPolicy.hpp"

#include <string>

namespace ls2g {

/**
Continuation powerflow (CPF): traces the solution curve from a BASE injection state
to a TARGET one, following the classical predictor / corrector scheme, and stops at
the "nose" (the steady-state loading limit) or at a requested lambda.

--- Why this is a batch algorithm -------------------------------------------------

A continuation is a sequence of powerflows on a FIXED Ybus where only the injection
moves -- exactly the premise the batch layer is built on. Deriving from
`BaseBatchSolverSynch` (rather than reimplementing a loop around `LSGrid::ac_pf`) is
what buys the whole point of doing this in C++: `_finish_preprocessing` runs the base
case with `tell_all_changed()`, which builds the Jacobian sparsity and pays the ONE
symbolic factorization, and leaves the control in `tell_none_changed()`; every
corrector after that refactorizes numerically only. For KLU on a non-trivial grid the
symbolic analysis dominates, so a curve of 100 points costs one analyze, not 100.

It is a SIBLING of `BaseBatchSweep`, not a fifth instantiation of it, because what
distinguishes a continuation is not which of (Ybus, Sbus) varies -- in that taxonomy
it is exactly TimeSeries' cell -- but WHO DECIDES THE NEXT ROW. Three consequences
the four-algorithm template cannot express:
  - the number of points is not known in advance (the step adapts, and the curve
    stops at the nose);
  - a diverged solve is not a result but a RETRY signal (halve the step and try the
    same point again), where `BaseBatchSweep::_store_row_status` would record it as a
    failed row and a chained instantiation would abandon the range;
  - the rows are strictly sequential, so there is nothing to thread.

--- The algorithm -----------------------------------------------------------------

The parametrised system is

    F(x, lam) = Scomp(x) - Sbus_base - lam . dir = 0,   dir = Sbus_target - Sbus_base

so lam = 0 is the base case and lam = 1 the target case (MATPOWER's `runcpf`
convention, and the reason the input is a target state rather than a scalar factor).

PREDICTOR. Differentiating, J . (dx/dlam) = dir projected onto the mismatch rows (see
NRSystem::cpf_rhs_into for the sign, which is derived from the residual convention
rather than assumed). Writing z = dx/dlam, the tangent is normalised:

    norm  = sqrt(1 + z.z),   t_lam = 1 / norm,   t_x = z / norm

and the step is taken along that unit tangent, so `step` is an arc length and stays
meaningful as the curve turns. The solve reuses the standing factorization -- one
triangular solve, no refactorization (see BaseAlgo::cpf_tangent).

Note that with this (natural) parameterisation t_lam is 1/norm and so is strictly
POSITIVE by construction: it cannot change sign at the nose, it tends to zero there,
as ||z|| -> infinity when J becomes singular. `nose_tol` is therefore a threshold on
t_lam, not a sign test.

CORRECTOR. An ordinary Newton-Raphson at the predicted, FIXED lambda. This is what
makes this a phase-1 continuation: at the nose the corrector's system is singular, so
the curve stops there and the lower branch is not traced. Rounding the nose needs the
lambda unknown and a parameterisation equation inside the Jacobian (MATPOWER's
`cpf.parameterization` 2 / 3), i.e. an NRSystem extension; that is deliberately left
for later.

--- Defaults ----------------------------------------------------------------------

The option names and default values follow MATPOWER's `cpf.*` so that a user who
knows `runcpf` is not faced with new concepts: step 0.05, step_min 1e-4, step_max
0.2, adapt_step off, adapt_step_damping 0.7, adapt_step_tol 1e-3, nose_tol 1e-5,
and "stop at the nose" unless a target lambda is given. The one place this cannot
follow MATPOWER is the parameterisation, whose default there is pseudo arc length
(see above).
**/
class LS2G_API ContinuationSweep final : public BaseBatchSolverSynch
{
    public:
        explicit ContinuationSweep(const LSGrid & init_grid_model):
            BaseBatchSolverSynch(init_grid_model),
            _step(static_cast<real_type>(0.05)),
            _step_min(static_cast<real_type>(1e-4)),
            _step_max(static_cast<real_type>(0.2)),
            _adapt_step(false),
            _adapt_step_damping(static_cast<real_type>(0.7)),
            _adapt_step_tol(static_cast<real_type>(1e-3)),
            _nose_tol(static_cast<real_type>(1e-5)),
            _stop_at_lam(-1.),
            _max_steps(1000),
            _exact_tangent(false),
            _status(0),
            _nb_points(0),
            _nb_retries(0)
            {}

        // Strictly sequential: point k is predicted from point k-1's tangent. Reported
        // here so set_nb_thread(>1) refuses with the base class' own message.
        bool supports_multithread() const override {return false;}

        // ---- the target injection state ----------------------------------------
        // One entry per grid element of the relevant type, in MW / MVAr, i.e. the
        // same units and ordering as LSGrid::get_gen_target_p() and friends. An axis
        // never set (or set empty) means "same as the base case", so it contributes
        // nothing to the direction. Setting none of them at all is an error: the
        // direction would be zero and the continuation would walk lambda along a
        // curve that does not move.
        void set_target_gen_p(const Eigen::Ref<const RealVect> & values)  {_target_gen_p = values;}
        void set_target_sgen_p(const Eigen::Ref<const RealVect> & values) {_target_sgen_p = values;}
        void set_target_load_p(const Eigen::Ref<const RealVect> & values) {_target_load_p = values;}
        void set_target_load_q(const Eigen::Ref<const RealVect> & values) {_target_load_q = values;}
        void clear_target(){
            _target_gen_p = RealVect(); _target_sgen_p = RealVect();
            _target_load_p = RealVect(); _target_load_q = RealVect();
        }

        // ---- options (MATPOWER cpf.* names) ------------------------------------
        real_type get_step() const {return _step;}
        void set_step(real_type x) {_step = _checked_positive(x, "step");}
        real_type get_step_min() const {return _step_min;}
        void set_step_min(real_type x) {_step_min = _checked_positive(x, "step_min");}
        real_type get_step_max() const {return _step_max;}
        void set_step_max(real_type x) {_step_max = _checked_positive(x, "step_max");}
        bool get_adapt_step() const {return _adapt_step;}
        void set_adapt_step(bool on) {_adapt_step = on;}
        real_type get_adapt_step_damping() const {return _adapt_step_damping;}
        void set_adapt_step_damping(real_type x) {_adapt_step_damping = _checked_positive(x, "adapt_step_damping");}
        real_type get_adapt_step_tol() const {return _adapt_step_tol;}
        void set_adapt_step_tol(real_type x) {_adapt_step_tol = _checked_positive(x, "adapt_step_tol");}
        real_type get_nose_tol() const {return _nose_tol;}
        void set_nose_tol(real_type x) {_nose_tol = _checked_positive(x, "nose_tol");}

        // Stop when lambda reaches this value (MATPOWER's numeric `cpf.stop_at`).
        // A non-positive value means "trace until the nose" ('NOSE').
        real_type get_stop_at_lam() const {return _stop_at_lam;}
        void set_stop_at_lam(real_type x) {_stop_at_lam = x;}

        // Hard cap on the number of traced points, so a pathological curve cannot
        // spin forever. Also the number of rows reserved up front.
        int get_max_steps() const {return _max_steps;}
        void set_max_steps(int x);

        // Rebuild and refactorize J at each converged point before taking its
        // tangent. Off by default: compute_pf leaves a factorization from its last
        // refactorizing iteration, which is one iterate behind the converged point --
        // approximate, but the corrector absorbs it. Worth turning on with a lazy
        // RefactorPolicy (Chord), where "one iterate behind" can be several.
        bool get_exact_tangent() const {return _exact_tangent;}
        void set_exact_tangent(bool on) {_exact_tangent = on;}

        // ---- run ----------------------------------------------------------------
        void compute(const Eigen::Ref<const CplxVect> & Vinit, int max_iter, real_type tol);

        // ---- results ------------------------------------------------------------
        // 1 if the curve was traced to its requested end (the nose, or stop_at_lam),
        // 0 otherwise (base case diverged, or the run stopped early -- see get_msg).
        int get_status() const {return _status;}
        const std::string & get_msg() const {return _msg;}

        // Number of traced points, base case included. get_voltages() has exactly
        // this many rows (it is truncated at the end of compute()).
        Eigen::Index nb_points() const {return _nb_points;}

        // lambda at each traced point; lam(0) == 0 (the base case).
        Eigen::Ref<const RealVect> get_lam() const {return _lam;}
        // the tangent's lambda component at each traced point, in (0, 1]. It tends to
        // zero at the nose -- that is the collapse indicator, see the class comment.
        // The last point has no tangent computed after it and reads 0.
        Eigen::Ref<const RealVect> get_tangent_lam() const {return _tangent_lam;}
        // largest lambda reached (0 if only the base case was solved).
        real_type get_lam_max() const {return _nb_points > 0 ? _lam(_nb_points - 1) : static_cast<real_type>(0.);}
        // how many times a corrector failed and the step had to be reduced.
        int nb_retries() const {return _nb_retries;}

        // The direction actually used, in solver ordering and per unit (empty before
        // the first compute()). Exposed for tests and for anyone who wants to check
        // what a steering vector resolved to.
        Eigen::Ref<const CplxVect> get_direction_solver() const {return _direction;}

        // Branch flows at every traced point, same contract as the other batch
        // algorithms' (one row per point, one column per powerline then trafo).
        Eigen::Ref<RealMat> compute_flows() {
            compute_flows_from_Vs();
            return _amps_flows;
        }
        Eigen::Ref<RealMat> compute_power_flows() {
            compute_flows_from_Vs(false);
            return _active_power_flows;
        }

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}

        // L3 (see BaseBatchSolverSynch's block comment on the three cache levels):
        // the curve this class traces, on top of the voltages the base holds.
        void clear_batch_outputs() override {
            _direction = CplxVect();
            _lam = RealVect();
            _tangent_lam = RealVect();
            _nb_points = 0;
            _nb_retries = 0;
            _status = 0;
            _msg.clear();
            BaseBatchSolverSynch::clear_batch_outputs();
        }

        void clear() override {
            BaseBatchSolverSynch::clear();   // -> L1 -> L2 -> L3
            clear_target();                  // the target is a registration, not a level
        }

    protected:
        // Sbus_target - Sbus_base, in solver ordering and per unit, built from the
        // four target axes through SbusPolicy::Vary's own element -> solver-bus
        // scatter helpers (the same ones TimeSeries uses, so a load / generator maps
        // to its bus, aggregates with the others on it, and gets divided by sn_mva
        // exactly as in an ordinary solve). Everything the four axes do NOT cover is
        // identical in the base and the target state and therefore cancels in the
        // difference -- which is why this needs no equivalent of
        // SbusPolicy::Vary::constant_sbus_pu.
        void _build_direction(int nb_buses_solver, const SolverBusIdVect & id_me_to_solver);

        // one target axis as a 1-row matrix of DELTAS (target - base), or a 1 x n
        // zero matrix when the axis was never set. `base` is the grid's own target
        // vector for that element type.
        static SbusPolicy::Vary::RealMat _delta_row(const Eigen::Ref<const RealVect> & target,
                                                    const Eigen::Ref<const RealVect> & base,
                                                    const char * axis_name);

        static real_type _checked_positive(real_type x, const char * name);

        // Re-establish the algorithm's state at the last converged point after a
        // corrector diverged: the failed solve left (Va, Vm) and the standing
        // factorization at some diverged iterate, and both feed the next tangent.
        // Re-solving from V_last converges immediately -- and precisely BECAUSE it
        // converges immediately it never enters the NR loop and so never
        // refactorizes, which is why the refactorization is asked for explicitly.
        bool _restore_at(const Eigen::Ref<const CplxVect> & V_last, real_type lam_last,
                         int max_iter, real_type tol);

        // MATPOWER's predictor-error step control (cpf.adapt_step):
        //   err        = ||[theta(pv,pq); |V|(pq)]_corrected - [..]_predicted||_inf
        //   step_scale = min(2, 1 + damping . (tol/err - 1))
        // clipped to [step_min, step_max]. MATPOWER's error also carries a lambda
        // term; here it is structurally zero (see the definition).
        real_type _adapted_step(real_type step,
                                const Eigen::Ref<const CplxVect> & V_corr,
                                const Eigen::Ref<const CplxVect> & V_pred) const;

    private:
        // target state (empty axis == unchanged from base)
        RealVect _target_gen_p, _target_sgen_p, _target_load_p, _target_load_q;

        // options
        real_type _step;
        real_type _step_min;
        real_type _step_max;
        bool      _adapt_step;
        real_type _adapt_step_damping;
        real_type _adapt_step_tol;
        real_type _nose_tol;
        real_type _stop_at_lam;
        int       _max_steps;
        bool      _exact_tangent;

        // results
        int         _status;
        std::string _msg;
        Eigen::Index _nb_points;
        int         _nb_retries;
        RealVect    _lam;
        RealVect    _tangent_lam;
        CplxVect    _direction;
};

} // namespace ls2g

#endif // CONTINUATIONSWEEP_H
