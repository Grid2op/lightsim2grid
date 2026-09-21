// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ContinuationSweep.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace ls2g {

real_type ContinuationSweep::_checked_positive(real_type x, const char * name)
{
    if(!std::isfinite(x) || (x <= 0.)){
        std::ostringstream exc_;
        exc_ << "ContinuationSweep: " << name << " must be a finite, strictly positive "
             << "number, got " << x << ".";
        throw std::runtime_error(exc_.str());
    }
    return x;
}

void ContinuationSweep::set_max_steps(int x)
{
    if(x < 1){
        std::ostringstream exc_;
        exc_ << "ContinuationSweep: max_steps must be >= 1, got " << x << ".";
        throw std::runtime_error(exc_.str());
    }
    _max_steps = x;
}

SbusPolicy::Vary::RealMat ContinuationSweep::_delta_row(const Eigen::Ref<const RealVect> & target,
                                                        const Eigen::Ref<const RealVect> & base,
                                                        const char * axis_name)
{
    SbusPolicy::Vary::RealMat res = SbusPolicy::Vary::RealMat::Zero(1, base.size());
    if(target.size() == 0) return res;  // axis never set: no contribution
    if(target.size() != base.size()){
        std::ostringstream exc_;
        exc_ << "ContinuationSweep::set_target_" << axis_name << ": expected one value per "
             << "element (" << base.size() << "), got " << target.size() << ".";
        throw std::runtime_error(exc_.str());
    }
    res.row(0) = (target - base).transpose();
    return res;
}

void ContinuationSweep::_build_direction(int nb_buses_solver,
                                         const SolverBusIdVect & id_me_to_solver)
{
    const auto & generators = _grid_model.get_generators_as_data();
    const auto & s_generators = _grid_model.get_static_generators_as_data();
    const auto & loads = _grid_model.get_loads_as_data();
    const char * algo_name = "ContinuationSweep";

    const SbusPolicy::Vary::RealMat d_gen_p =
        _delta_row(_target_gen_p, _grid_model.get_gen_target_p(), "gen_p");
    const SbusPolicy::Vary::RealMat d_sgen_p =
        _delta_row(_target_sgen_p, _grid_model.get_sgen_target_p(), "sgen_p");
    const SbusPolicy::Vary::RealMat d_load_p =
        _delta_row(_target_load_p, _grid_model.get_load_target_p(), "load_p");
    const SbusPolicy::Vary::RealMat d_load_q =
        _delta_row(_target_load_q, loads.get_target_q(), "load_q");

    // Same scatter, same signs (generation adds, load subtracts) and same sn_mva
    // division as an ordinary Sbus build -- these are SbusPolicy::Vary's own helpers,
    // used here on the DELTAS rather than on two absolute states, which is the same
    // thing by linearity and saves assembling the base and target injections only to
    // subtract them.
    SbusPolicy::Vary::CplxMat dir = SbusPolicy::Vary::CplxMat::Zero(1, nb_buses_solver);
    bool add_ = true;
    SbusPolicy::Vary::fill_SBus_real(dir, generators, d_gen_p, id_me_to_solver, add_, algo_name);
    SbusPolicy::Vary::fill_SBus_real(dir, s_generators, d_sgen_p, id_me_to_solver, add_, algo_name);
    add_ = false;
    SbusPolicy::Vary::fill_SBus_real(dir, loads, d_load_p, id_me_to_solver, add_, algo_name);
    SbusPolicy::Vary::fill_SBus_imag(dir, loads, d_load_q, id_me_to_solver, add_, algo_name);

    const real_type sn_mva = _grid_model.get_sn_mva();
    if(std::abs(sn_mva - 1.0) > BaseConstants::_tol_equal_float){
        dir.array() /= static_cast<cplx_type>(sn_mva);
    }
    _direction = dir.row(0).transpose();
}

bool ContinuationSweep::_restore_at(const Eigen::Ref<const CplxVect> & V_last, real_type lam_last,
                                    int max_iter, real_type tol)
{
    CplxVect V = V_last;
    const CplxVect Sbus = ac_cache_.inj + lam_last * _direction;
    const bool conv = compute_one_powerflow(
        ac_cache_.mat, V, Sbus,
        active_layout().slack_bus_id_solver.as_eigen(), active_layout().slack_weights,
        active_layout().bus_pv.as_eigen(), active_layout().bus_pq.as_eigen(),
        max_iter, tol);
    if(!conv) return false;
    // V_last is already a solution, so the solve above converged on its very first
    // convergence check and never entered the NR loop -- which means it never
    // refactorized either, and the standing factorization is still the diverged one
    // this restore exists to discard. Rebuild it explicitly.
    return _algo.cpf_refactorize_at_current();
}

real_type ContinuationSweep::_adapted_step(real_type step,
                                           const Eigen::Ref<const CplxVect> & V_corr,
                                           const Eigen::Ref<const CplxVect> & V_pred) const
{
    // MATPOWER's cpf_error: the infinity norm of (corrected - predicted) over the
    // free state variables -- angles at pv and pq buses, magnitudes at pq buses --
    // and lambda. The lambda term is omitted here because it is structurally zero:
    // this parameterisation's corrector holds lambda FIXED at its predicted value
    // (that is what makes it "natural"), so lam_corrected == lam_predicted exactly.
    // It reappears with a parameterisation that solves for lambda too.
    real_type err = 0.;
    const auto pv = active_layout().bus_pv.as_eigen();
    const auto pq = active_layout().bus_pq.as_eigen();
    for(Eigen::Index k = 0; k < pv.size(); ++k){
        const int b = pv(k);
        err = std::max(err, std::abs(std::arg(V_corr(b)) - std::arg(V_pred(b))));
    }
    for(Eigen::Index k = 0; k < pq.size(); ++k){
        const int b = pq(k);
        err = std::max(err, std::abs(std::arg(V_corr(b)) - std::arg(V_pred(b))));
        err = std::max(err, std::abs(std::abs(V_corr(b)) - std::abs(V_pred(b))));
    }
    if(!(err > 0.)) return _step_max;  // a perfect prediction: take the largest step

    const real_type scale = std::min(static_cast<real_type>(2.),
                                     static_cast<real_type>(1.) +
                                     _adapt_step_damping * (_adapt_step_tol / err - static_cast<real_type>(1.)));
    real_type res = step * scale;
    if(res > _step_max) res = _step_max;
    if(res < _step_min) res = _step_min;
    return res;
}

void ContinuationSweep::compute(const Eigen::Ref<const CplxVect> & Vinit,
                                int max_iter, real_type tol)
{
    auto timer = CustTimer();
    auto timer_preproc = CustTimer();

    const size_t nb_total_bus = _reset_data_and_check_vinit(Vinit);
    _status = 0;
    _msg.clear();
    _nb_points = 0;
    _nb_retries = 0;
    _timer_thread_init = 0.;

    if(!_algo.ac_solver_used()){
        throw std::runtime_error("ContinuationSweep::compute: a continuation powerflow is an AC "
                                 "method; the current algorithm is a DC one. Pick an AC "
                                 "(Newton-Raphson) algorithm with change_algorithm().");
    }
    if(!_algo.supports_cpf()){
        std::ostringstream exc_;
        exc_ << "ContinuationSweep::compute: the algorithm '" << _algo.get_name() << "' is not "
             << "Newton-Raphson based, so it has no Jacobian to take a tangent from. Pick one of "
             << "the NR algorithms (eg NR_KLU / NR_SparseLU).";
        throw std::runtime_error(exc_.str());
    }
    if(_step_min > _step_max){
        std::ostringstream exc_;
        exc_ << "ContinuationSweep::compute: step_min (" << _step_min << ") is larger than "
             << "step_max (" << _step_max << ").";
        throw std::runtime_error(exc_.str());
    }

    const bool ac_solver_used = true;
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);
    _build_direction(nb_buses_solver_, active_layout().id_me_to_solver);

    // A zero direction is not a degenerate curve, it is a meaningless run: J . z = 0
    // gives z = 0, the tangent normalises to (0, ..., 1), and lambda would march to
    // its target along a curve on which nothing moves -- reporting success for a
    // continuation that continued nothing. Refuse it instead.
    if(_direction.size() == 0 || _direction.cwiseAbs().maxCoeff() <= 0.){
        throw std::runtime_error("ContinuationSweep::compute: the direction is zero -- the target "
                                 "state is identical to the base state (or no target was set, or "
                                 "every steered element is disconnected). Set at least one of "
                                 "set_target_gen_p / set_target_sgen_p / set_target_load_p / "
                                 "set_target_load_q to something that differs from the grid's own "
                                 "values.");
    }

    // Reserve the whole curve up front: the base case plus at most _max_steps points.
    // _voltages is truncated to the points actually traced at the end of this
    // function, so get_voltages() / compute_flows_from_Vs() see exactly nb_points()
    // rows and nothing has to know about the reservation.
    const Eigen::Index nb_rows = static_cast<Eigen::Index>(_max_steps) + 1;
    _lam = RealVect::Zero(nb_rows);
    _tangent_lam = RealVect::Zero(nb_rows);

    // The base case ("n" powerflow): this is what builds the Jacobian sparsity and
    // pays the ONE symbolic factorization for the whole curve. It leaves the control
    // in "nothing changed", so every corrector below refactorizes numerically only.
    const bool base_conv = _finish_preprocessing(
        static_cast<size_t>(nb_rows), nb_total_bus, Vinit_solver, static_cast<size_t>(max_iter),
        tol, timer_preproc);

    if(!base_conv){
        std::ostringstream msg_;
        msg_ << "the base case (lambda = 0) did not converge (error: " << _algo.get_error()
             << "); nothing can be continued from it.";
        _msg = msg_.str();
        _voltages = CplxMat::Zero(0, static_cast<Eigen::Index>(nb_total_bus));
        _lam = RealVect(); _tangent_lam = RealVect();
        _timer_total = timer.duration();
        return;
    }

    // ---- point 0: the base case ------------------------------------------------
    CplxVect V_last = _algo.get_V();
    real_type lam = 0.;
    _voltages.row(0)(active_layout().id_solver_to_me.as_eigen()) = V_last.array();
    _lam(0) = 0.;
    _nb_points = 1;

    real_type step = _step;
    // Every buffer the loop needs, allocated ONCE here and assigned into afterwards:
    // declared inside, each iteration would construct and free a fresh nb_bus vector,
    // and a curve is hundreds of iterations long. They are four distinct vectors
    // because each is read after the next is written: V_pred feeds the step adaptation
    // AFTER the corrector has overwritten V_corr, and V_last must survive a corrector
    // failure so the retry can be re-seeded from it (see _restore_at).
    RealVect z;
    CplxVect V_pred;
    CplxVect V_corr;
    CplxVect Sbus;

    for(int it = 0; it < _max_steps; ++it){
        // ---- predictor ---------------------------------------------------------
        if(_exact_tangent && !_algo.cpf_refactorize_at_current()){
            _msg = "could not refactorize the Jacobian at a converged point (exact_tangent).";
            break;
        }
        if(!_algo.cpf_tangent(_direction, z)){
            std::ostringstream msg_;
            msg_ << "the tangent system could not be solved at lambda = " << lam
                 << " (error: " << _algo.get_error() << "); the Jacobian is most likely singular, "
                 << "which happens AT the nose point.";
            _msg = msg_.str();
            break;
        }

        // unit tangent: t = [z; 1] / ||[z; 1]||, so `step` is an arc length and t_lam
        // is 1/norm -- strictly positive, tending to 0 as the curve turns (see the
        // class comment: this parameterisation cannot show a sign change).
        const real_type norm = std::sqrt(static_cast<real_type>(1.) + z.squaredNorm());
        const real_type t_lam = static_cast<real_type>(1.) / norm;
        _tangent_lam(_nb_points - 1) = t_lam;

        if(t_lam < _nose_tol){
            std::ostringstream msg_;
            msg_ << "nose point reached at lambda = " << lam << " (the tangent's lambda component "
                 << "fell to " << t_lam << ", below nose_tol = " << _nose_tol << ").";
            _msg = msg_.str();
            _status = 1;
            break;
        }

        real_type this_step = step;
        real_type lam_pred = lam + this_step * t_lam;
        bool last_point = false;
        if(_stop_at_lam > 0. && lam_pred >= _stop_at_lam){
            // do not overshoot the requested lambda: shrink this step so the
            // predictor lands exactly on it.
            this_step = (_stop_at_lam - lam) / t_lam;
            lam_pred = _stop_at_lam;
            last_point = true;
        }
        _algo.cpf_predict(z, this_step / norm, V_pred);

        // A direction the slack absorbs entirely. Unlike a zero direction this is a
        // perfectly legal thing to ask for -- scaling the machine at a single slack bus
        // moves a real input, it is simply an input the powerflow does not read: that
        // bus' P equation exists and takes the direction, but the slack-absorbed unknown
        // cancels it one for one, so the machine's actual output is unchanged and no bus
        // moves. So this is NOT refused. What is refused is pretending it went somewhere:
        // continuing would march lambda to its target over hundreds of identical points
        // and report a loading margin that is really an artefact of max_steps. Stop here
        // instead, with the base case as the only point and a message saying why.
        //
        // Checked on the first predictor only: which unknowns a direction reaches is a
        // property of it and of the topology, and neither changes along the curve.
        if (it == 0 &&
            (V_pred - V_last).cwiseAbs().maxCoeff() <= BaseConstants::_tol_equal_float) {
            _msg = "the direction moves no bus voltage: everything it asks for is absorbed by "
                   "the slack. Scaling the machine at a single slack bus does this -- it "
                   "changes a real input, but not one the powerflow reads, since that bus' "
                   "generation is whatever balances the grid rather than what was asked of "
                   "it. There is no loading margin along such a direction, so nothing is "
                   "traced past the base case; steer something whose power the powerflow "
                   "actually has to route.";
            break;
        }

        // ---- corrector ---------------------------------------------------------
        V_corr = V_pred;   // compute_one_powerflow overwrites it with the solution
        Sbus = ac_cache_.inj + lam_pred * _direction;
        const bool conv = compute_one_powerflow(
            ac_cache_.mat, V_corr, Sbus,
            active_layout().slack_bus_id_solver.as_eigen(), active_layout().slack_weights,
            active_layout().bus_pv.as_eigen(), active_layout().bus_pq.as_eigen(),
            max_iter, tol);

        if(!conv){
            // Not a failed row: a step that was too long. Halve it and try the same
            // point again -- unless we are already at step_min, in which case there
            // is no solution just ahead and this IS the end of the curve (in
            // practice, the nose).
            ++_nb_retries;
            if(this_step <= _step_min){
                std::ostringstream msg_;
                msg_ << "the corrector stopped converging at lambda = " << lam
                     << " with the smallest allowed step (" << _step_min << "): this is the end "
                     << "of the traceable curve, ie the nose point (to go past it, a "
                     << "parameterisation that keeps the corrector non-singular is needed).";
                _msg = msg_.str();
                _status = 1;
                break;
            }
            step = std::max(_step_min, this_step * static_cast<real_type>(0.5));
            // the diverged solve left (Va, Vm) and the standing factorization at some
            // diverged iterate; both feed the next tangent, so put them back.
            if(!_restore_at(V_last, lam, max_iter, tol)){
                _msg = "could not re-establish the last converged point after a corrector failure.";
                break;
            }
            continue;
        }

        // ---- accept the point --------------------------------------------------
        if(_adapt_step) step = _adapted_step(step, V_corr, V_pred);

        lam = lam_pred;
        V_last = V_corr;
        _voltages.row(_nb_points)(active_layout().id_solver_to_me.as_eigen()) = V_last.array();
        _lam(_nb_points) = lam;
        ++_nb_points;

        if(last_point){
            std::ostringstream msg_;
            msg_ << "reached the requested lambda = " << _stop_at_lam << ".";
            _msg = msg_.str();
            _status = 1;
            break;
        }
    }

    if(_msg.empty()){
        std::ostringstream msg_;
        msg_ << "stopped after the maximum number of steps (" << _max_steps << ") at lambda = "
             << lam << "; raise max_steps, or step, to go further.";
        _msg = msg_.str();
    }

    // truncate every result to the points actually traced
    _voltages.conservativeResize(_nb_points, _voltages.cols());
    _lam.conservativeResize(_nb_points);
    _tangent_lam.conservativeResize(_nb_points);

    _timer_total = timer.duration();
}

} // namespace ls2g
