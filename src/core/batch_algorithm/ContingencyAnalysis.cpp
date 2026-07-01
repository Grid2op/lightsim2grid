// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ContingencyAnalysis.hpp"

#include <queue>
#include <algorithm>
#include <math.h>       /* isfinite */
#include <thread>
#include <vector>
#include <memory>

namespace ls2g {

bool ContingencyAnalysis::check_invertible(const Eigen::SparseMatrix<cplx_type> & Ybus) const{
    std::vector<bool> visited(Ybus.cols(), false); 
    std::vector<bool> already_added(Ybus.cols(), false);
    std::queue<Eigen::Index> neighborhood;
    size_t col_id = 0;  // start by node 0, why not
    while (true)
    {
        visited[col_id] = true;
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(Ybus, col_id); it; ++it)
        {
            // add in the queue all my neighbor (if the coefficient is big enough)
            if(!visited[it.row()] && !already_added[it.row()] && abs(it.value()) > 1e-8){
                neighborhood.push(it.row());
                already_added[it.row()] = true;
            }
        }
        if(neighborhood.empty()) break;
        col_id = neighborhood.front();
        neighborhood.pop();
    }
    
    bool ok = true;
    for(auto el: visited){
        if(!el)
        {
            // this node has not been visited, there is an error
            ok=false;
            break;
        }
    }
    return ok;
}

template<typename T>
std::vector<int> ContingencyAnalysis::disconnected_buses(const Eigen::SparseMatrix<T> & mat) const{
    const int n = static_cast<int>(mat.cols());
    std::vector<int> comp_of_bus(n, -1);  // connected-component label of each bus
    int nb_comp = 0;
    std::queue<int> neighborhood;
    for(int start = 0; start < n; ++start){
        if(comp_of_bus[start] != -1) continue;  // already labelled
        // BFS labelling the component that contains "start"
        comp_of_bus[start] = nb_comp;
        neighborhood.push(start);
        while(!neighborhood.empty()){
            const int col_id = neighborhood.front();
            neighborhood.pop();
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, col_id); it; ++it){
                const int row = static_cast<int>(it.row());
                if(comp_of_bus[row] == -1 && std::abs(it.value()) > 1e-8){
                    comp_of_bus[row] = nb_comp;
                    neighborhood.push(row);
                }
            }
        }
        ++nb_comp;
    }

    if(nb_comp <= 1) return std::vector<int>();  // grid is connected: nothing masked

    // find the largest component (by number of buses)
    std::vector<int> nb_bus_per_comp(nb_comp, 0);
    for(int bus_id = 0; bus_id < n; ++bus_id) nb_bus_per_comp[comp_of_bus[bus_id]] += 1;
    const int main_comp = static_cast<int>(std::distance(
        nb_bus_per_comp.begin(),
        std::max_element(nb_bus_per_comp.begin(), nb_bus_per_comp.end())));

    // every bus not in the main component is masked
    std::vector<int> masked;
    masked.reserve(n - nb_bus_per_comp[main_comp]);
    for(int bus_id = 0; bus_id < n; ++bus_id){
        if(comp_of_bus[bus_id] != main_comp) masked.push_back(bus_id);
    }
    return masked;
}

// explicit instantiations: AC works on the complex Ybus_, DC on the real Bbus_
template std::vector<int> ContingencyAnalysis::disconnected_buses<cplx_type>(const Eigen::SparseMatrix<cplx_type> &) const;
template std::vector<int> ContingencyAnalysis::disconnected_buses<real_type>(const Eigen::SparseMatrix<real_type> &) const;

void ContingencyAnalysis::select_ref_slack_and_masks(){
    const size_t nb_cont = _li_coeffs.size();
    _li_masked.assign(nb_cont, std::vector<int>());
    _skip_mask.assign(nb_cont, 0);

    // 1) per-contingency masked bus set (largest component is kept). Work on a copy of
    // the admittance matrix so the member used by the n-powerflow is left untouched.
    // AC uses the complex Ybus_; DC uses the real Bbus_ (Ybus_ is empty in DC mode).
    if(_algo.ac_solver_used()){
        Eigen::SparseMatrix<cplx_type> Ybus = Ybus_;
        size_t cont_id = 0;
        for(const auto & coeffs_modif: _li_coeffs){
            for(const auto & c: coeffs_modif) Ybus.coeffRef(c.row_id, c.col_id) -= c.value;
            _li_masked[cont_id] = disconnected_buses(Ybus);
            for(const auto & c: coeffs_modif) Ybus.coeffRef(c.row_id, c.col_id) += c.value;
            ++cont_id;
        }
    } else {
        Eigen::SparseMatrix<real_type> Bbus = Bbus_;
        size_t cont_id = 0;
        for(const auto & coeffs_modif: _li_coeffs){
            for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) -= std::real(c.value);
            _li_masked[cont_id] = disconnected_buses(Bbus);
            for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) += std::real(c.value);
            ++cont_id;
        }
    }

    // 2) candidate reference slacks (treated as solver bus ids, consistent with the
    // rest of the batch path). Prefer the ones that actually carry slack power; if
    // none, fall back to any slack (a 0-weight slack can still anchor the angle).
    const Eigen::Index nb_slack = slack_ids_me_.size();
    std::vector<int> candidates;
    for(Eigen::Index i = 0; i < nb_slack; ++i){
        const int b = slack_ids_me_[static_cast<int>(i)].cast_int();
        if(b >= 0 && b < slack_weights_.size() && slack_weights_(b) > 0.) candidates.push_back(b);
    }
    if(candidates.empty()){
        for(Eigen::Index i = 0; i < nb_slack; ++i) candidates.push_back(slack_ids_me_[static_cast<int>(i)].cast_int());
    }
    if(candidates.empty()) return;  // no slack at all: nothing we can do

    // 3) choose the candidate stranded by the fewest contingencies (ties: larger
    // slack weight, then smaller bus id)
    auto is_masked = [](const std::vector<int> & masked, int bus){
        return std::find(masked.begin(), masked.end(), bus) != masked.end();
    };
    int best_bus = candidates[0];
    int best_strand = -1;
    real_type best_weight = -1.;
    for(int bus : candidates){
        int strand = 0;
        for(const auto & masked : _li_masked) if(is_masked(masked, bus)) ++strand;
        const real_type weight = (bus < slack_weights_.size()) ? slack_weights_(bus) : 0.;
        const bool better = (best_strand < 0) ||
                            (strand < best_strand) ||
                            (strand == best_strand && weight > best_weight) ||
                            (strand == best_strand && weight == best_weight && bus < best_bus);
        if(better){ best_strand = strand; best_weight = weight; best_bus = bus; }
    }

    // 4) reorder slack_ids_me_ so the chosen reference is index 0 (the NR uses
    // slack_ids[0] as the angle reference). Done before the (n-)powerflow so the
    // symbolic factorization is built once with this reference.
    {
        std::vector<int> reordered;
        reordered.reserve(static_cast<size_t>(nb_slack));
        reordered.push_back(best_bus);
        for(Eigen::Index i = 0; i < nb_slack; ++i){
            const int b = slack_ids_me_[static_cast<int>(i)].cast_int();
            if(b != best_bus) reordered.push_back(b);
        }
        slack_ids_me_ = GlobalBusIdVect(reordered);
    }

    // 5) contingencies that still strand the chosen reference are skipped
    for(size_t cont_id = 0; cont_id < nb_cont; ++cont_id){
        if(is_masked(_li_masked[cont_id], best_bus)) _skip_mask[cont_id] = 1;
    }
}

RealVect ContingencyAnalysis::masked_slack_weights(const std::vector<int> & masked) const{
    RealVect w = slack_weights_;
    if(masked.empty()) return w;
    const real_type orig_sum = w.sum();
    for(int b : masked) if(b >= 0 && b < w.size()) w(b) = 0.;
    const real_type new_sum = w.sum();
    // rescale so the slack power is shared only among the live slacks (sum preserved).
    // if every slack-carrying bus is masked the vector stays as-is (degenerate: the
    // solve will simply fail to converge and the contingency is reported as skipped).
    if(new_sum > 1e-12 && orig_sum > 1e-12) w *= (orig_sum / new_sum);
    return w;
}

void ContingencyAnalysis::init_li_coeffs(
    bool ac_solver_used,
    const SolverBusIdVect &id_me_to_solver)
{
    _li_coeffs.clear();
    _li_coeffs.reserve(_li_defaults.size());
    const auto & powerlines = _grid_model.get_powerlines_as_data();
    const auto & trafos = _grid_model.get_trafos_as_data();
    // const auto & id_me_to_solver = ac_solver_used ? _grid_model.id_me_to_ac_solver(): _grid_model.id_me_to_dc_solver();
    int bus_1_id, bus_2_id;
    cplx_type y_ff, y_ft, y_tf, y_tt;
    bool status;
    for(const auto & this_cont_id: _li_defaults){
        std::vector<Coeff> this_cont_coeffs;
        this_cont_coeffs.reserve(this_cont_id.size() * 4);  // usually there are 4 coeffs per powerlines / trafos
        for(auto br_id : this_cont_id){
            int el_id;
            const TwoSidesContainer_rxh_A<OneSideContainer_ForBranch> *p_branch;
            if(br_id < n_line_)
            {
                // this is a powerline
                el_id = br_id;
                p_branch = & powerlines;
            }else{
                // this is a trafo
                el_id = br_id - n_line_;
                p_branch = & trafos;
            }
            
            GlobalBusId glob_bus_1 = p_branch->get_bus_side_1(el_id);
            GlobalBusId glob_bus_2 = p_branch->get_bus_side_2(el_id);
            bus_1_id = glob_bus_1.cast_int() == GenericContainer::_deactivated_bus_id ? 
                GenericContainer::_deactivated_bus_id : 
                id_me_to_solver[glob_bus_1.cast_int()].cast_int();
            bus_2_id = glob_bus_2.cast_int() == GenericContainer::_deactivated_bus_id ? 
                GenericContainer::_deactivated_bus_id : 
                id_me_to_solver[glob_bus_2.cast_int()].cast_int();
            status = p_branch->get_status_global()[el_id];
            if(ac_solver_used){
                y_ff = p_branch->yac_eff_11()[el_id];
                y_ft = p_branch->yac_eff_12()[el_id];
                y_tf = p_branch->yac_eff_21()[el_id];
                y_tt = p_branch->yac_eff_22()[el_id];
            }else{
                y_ff = p_branch->ydc_11()[el_id];
                y_ft = p_branch->ydc_12()[el_id];
                y_tf = p_branch->ydc_21()[el_id];
                y_tt = p_branch->ydc_22()[el_id];
            }

            if(status)
            {
                // element is connected, update coeffs based on status of each powerlines
                if((bus_1_id != GenericContainer::_deactivated_bus_id)) this_cont_coeffs.push_back({bus_1_id, bus_1_id, y_ff});
                if((bus_2_id != GenericContainer::_deactivated_bus_id)) this_cont_coeffs.push_back({bus_2_id, bus_2_id, y_tt});
                if((bus_1_id != GenericContainer::_deactivated_bus_id) && (bus_2_id != GenericContainer::_deactivated_bus_id)){
                    this_cont_coeffs.push_back({bus_1_id, bus_2_id, y_ft});
                    this_cont_coeffs.push_back({bus_2_id, bus_1_id, y_tf});
                }
            }
        }
        _li_coeffs.push_back(this_cont_coeffs);
    }
}

bool ContingencyAnalysis::remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                                           const std::vector<Coeff> & coeffs,
                                           bool ac_solver_used,
                                           AlgorithmSelector & algo)
{
    if(ac_solver_used)
    {
        for(const auto & coeff_to_remove: coeffs){
            Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) -= coeff_to_remove.value;
        }
        return check_invertible(Ybus);
    } else{
        // DC solver stores the ybus internally, I update it
        // instead of building it over and over
        for(const Coeff& coeff : coeffs){
            algo.update_internal_Ybus(coeff, false);  // false => remove the coeff (using -= )
        }
        // in DC mode the solver takes the responsibility
        // so Ybus is always "connected".
        return true;
    }
}

IntVect ContingencyAnalysis::is_grid_connected_after_contingency(){
    const bool ac_solver_used = _algo.ac_solver_used();
    // Build the solver inputs (the AC complex Ybus_ or the DC real Bbus_, plus
    // id_me_to_solver_) and the per-contingency coefficients if they are not available
    // yet (i.e. compute() was not called).
    // NB: we use the (correctly indexed) member Ybus_ / Bbus_, NOT
    // _grid_model.get_Ybus_solver(): the latter is the grid model's own Ybus, never built
    // by this class (it works on Ybus_ / Bbus_) and would be an empty 0x0 matrix here ->
    // out-of-bounds coeffRef.
    const bool inputs_ready = (ac_solver_used ? Ybus_.cols() != 0 : Bbus_.cols() != 0);
    if(!inputs_ready || _li_coeffs.size() != _li_defaults.size()){
        const size_t nb_total_bus = _grid_model.total_bus();
        CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                            {_grid_model.get_init_vm_pu(), 0.});
        prepare_solver_input_base(Vinit, ac_solver_used);
        init_li_coeffs(ac_solver_used, id_me_to_solver_);
    }
    IntVect res = IntVect::Constant(_li_coeffs.size(), 0);
    if(ac_solver_used){
        Eigen::SparseMatrix<cplx_type> Ybus = Ybus_;  // correctly-indexed copy
        int cont_id = 0;
        for(const auto & coeffs_modif: _li_coeffs){
            if(remove_from_Ybus(Ybus, coeffs_modif, true, _algo)) res(cont_id) = 1;
            else res(cont_id) = 0;
            readd_to_Ybus(Ybus, coeffs_modif, true, _algo);
            ++cont_id;
        }
    } else {
        // DC: BFS the (real) Bbus_ directly. Unlike the AC `remove_from_Ybus`, this does
        // not touch the solver's internal dc matrix, so it is side-effect free.
        Eigen::SparseMatrix<real_type> Bbus = Bbus_;  // correctly-indexed copy
        int cont_id = 0;
        for(const auto & coeffs_modif: _li_coeffs){
            for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) -= std::real(c.value);
            res(cont_id) = disconnected_buses(Bbus).empty() ? 1 : 0;
            for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) += std::real(c.value);
            ++cont_id;
        }
    }
    return res;
}

int ContingencyAnalysis::pick_reference_slack(){
    const bool ac_solver_used = _algo.ac_solver_used();
    // Build the solver inputs + per-contingency coefficients on demand (same
    // pattern as is_grid_connected_after_contingency / compute).
    const bool inputs_ready = (ac_solver_used ? Ybus_.cols() != 0 : Bbus_.cols() != 0);
    if(!inputs_ready || _li_coeffs.size() != _li_defaults.size()){
        const size_t nb_total_bus = _grid_model.total_bus();
        CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                            {_grid_model.get_init_vm_pu(), 0.});
        prepare_solver_input_base(Vinit, ac_solver_used);
        init_li_coeffs(ac_solver_used, id_me_to_solver_);
    }
    // select_ref_slack_and_masks() scores every candidate slack and reorders
    // slack_ids_me_ so the stranded-by-fewest reference is index 0.
    select_ref_slack_and_masks();
    if(slack_ids_me_.size() == 0) return -1;
    return slack_ids_me_[static_cast<int>(0)].cast_int();  // gridmodel bus id
}

void ContingencyAnalysis::readd_to_Ybus(
    Eigen::SparseMatrix<cplx_type> & Ybus,
    const std::vector<Coeff> & coeffs,
    bool ac_solver_used,
    AlgorithmSelector & algo)
{
    if(ac_solver_used){
        for(const Coeff & coeff_to_remove: coeffs){
            Ybus.coeffRef(coeff_to_remove.row_id, coeff_to_remove.col_id) += coeff_to_remove.value;
        }
    } else {
        // DC solver stores the ybus internally, I update it
        // instead of building it over and over
        for(const Coeff& coeff : coeffs){
            algo.update_internal_Ybus(coeff, true);  // true => add back the coeff (using += )
        }
    }
}

void ContingencyAnalysis::compute(const CplxVect & Vinit, int max_iter, real_type tol)
{
    auto timer = CustTimer();
    auto timer_preproc = CustTimer();

    // perform some initial checks and reset timers
    size_t nb_total_bus = _reset_data_and_check_vinit(Vinit);
    _timer_modif_Ybus = 0.;
    _timer_thread_init = 0.;

    // read from the grid the usefull information
    const auto & sn_mva = _grid_model.get_sn_mva();
    const bool ac_solver_used = _algo.ac_solver_used();
    size_t nb_steps = _li_defaults.size();

    // prepare the gridmodel (compute Ybus, Sbus etc.)
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);

    ///////////////////////////////////////
    // initialize what will change (here Ybus)
    // initialize properly the coefficients that I will need to remove
    init_li_coeffs(ac_solver_used, id_me_to_solver_);
    ////////////////////////////////////

    // "handle disconnected grid" mode: pre-compute the masked bus set of each
    // contingency and choose the reference slack (BEFORE the n-powerflow, so the
    // symbolic factorization is built once with that reference and reused). Supported
    // by the Newton-Raphson family (AC) and the native DC solver; both expose bus
    // masking. A non-NR AC algorithm (eg Gauss-Seidel / Fast-Decoupled) is rejected.
    const bool mask_mode = _handle_disconnected_grid;
    if(mask_mode){
        if(!_algo.supports_bus_masking()){
            throw std::runtime_error("ContingencyAnalysis: the `handle_disconnected_grid` mode "
                                     "requires a Newton-Raphson algorithm (AC) or the DC solver "
                                     "(the active algorithm does not support bus masking). Use "
                                     "`change_algorithm` to select an NR solver (e.g. NR_KLU / "
                                     "NR_SLU) or the DC solver.");
        }
        select_ref_slack_and_masks();
    }

    bool n_powerflow_has_conv = _finish_preprocessing(
        nb_steps,
        nb_total_bus,
        Vinit_solver,  // is modified if _init_from_n_powerflow is true
        max_iter,
        tol,
        timer_preproc
    );

    if(!n_powerflow_has_conv) return;

    auto timer_thread = CustTimer();
    // now perform the security analysis, possibly split across several threads
    const size_t nb_cont = _li_coeffs.size();
    const int nb_thread = std::min(static_cast<int>(nb_cont), std::max(1, _nb_thread));

    if(nb_thread <= 1){
        // single-threaded path: reuse the (already warmed-up) member solver,
        // member Ybus_ and the member accumulators -> identical to the legacy code.
        std::exception_ptr err;
        run_contingency_range(
            0, nb_cont,
            _algo, _algo_controler, Ybus_, Vinit_solver,
            ac_solver_used, mask_mode, max_iter, tol, sn_mva,
            _timer_modif_Ybus, _nb_solved, _timer_solver, err, false);
        if(err) std::rethrow_exception(err);
        _timer_total = timer.duration();
        return;
    }

    // multi-threaded path: every thread owns its solver / control / Ybus copy.
    std::vector<std::unique_ptr<AlgorithmSelector> > algos(nb_thread);
    std::vector<AlgoControl> controls(nb_thread);
    std::vector<Eigen::SparseMatrix<cplx_type> > ybus_copies(nb_thread);
    std::vector<double> th_timer_modif(nb_thread, 0.);
    std::vector<int> th_nb_solved(nb_thread, 0);
    std::vector<double> th_timer_solver(nb_thread, 0.);
    std::vector<std::exception_ptr> th_err(nb_thread);

    // each worker builds its own solver / control / Ybus copy IN-THREAD so this
    // (non-trivial) setup -- in particular the full sparse Ybus_ copy -- runs in
    // parallel instead of serially on the main thread. Only reads shared state
    // (_grid_model, _algo config, get_algo_type()); each thread writes solely its
    // own slot of the pre-sized vectors, so no locking is needed.
    auto init_thread = [&](int t){
        algos[t] = std::make_unique<AlgorithmSelector>();
        AlgorithmSelector & algo = *algos[t];
        algo.set_lsgrid(&_grid_model);
        algo.change_algorithm(get_algo_type());
        algo.set_config(_algo.get_config());  // match the member solver's configuration
        controls[t] = _algo_controler;
        ybus_copies[t] = Ybus_;  // thread-local working copy of the admittance
    };

    // contiguous split of [0, nb_cont) across the threads
    auto algo_for = [&](int t) -> AlgorithmSelector & {
        return *algos[t];
    };
    auto ybus_for = [&](int t) -> Eigen::SparseMatrix<cplx_type> & {
        return ybus_copies[t];
    };
    const size_t base = nb_cont / static_cast<size_t>(nb_thread);
    const size_t rem = nb_cont % static_cast<size_t>(nb_thread);
    auto range_of = [&](int t, size_t & b, size_t & e){
        // contiguous split: the first `rem` threads get one extra contingency
        b = static_cast<size_t>(t) * base + std::min(static_cast<size_t>(t), rem);
        e = b + base + (static_cast<size_t>(t) < rem ? 1 : 0);
    };

    // spawn threads 1..nb_thread-1, then run thread 0's share inline
    std::vector<std::thread> threads;
    threads.reserve(nb_thread - 1);
    for(int t = 1; t < nb_thread; ++t){
        size_t b, e;
        range_of(t, b, e);
        threads.emplace_back([=, &init_thread, &algo_for, &ybus_for, &controls, &th_timer_modif,
                              &th_nb_solved, &th_timer_solver, &th_err, &Vinit_solver](){
            init_thread(t);  // build this worker's solver / control / Ybus copy in parallel
            run_contingency_range(
                b, e,
                algo_for(t), controls[t], ybus_for(t), Vinit_solver,
                ac_solver_used, mask_mode, max_iter, tol, sn_mva,
                th_timer_modif[t], th_nb_solved[t], th_timer_solver[t], th_err[t], true);
        });
    }
    _timer_thread_init = timer_thread.duration();

    {
        size_t b, e;
        range_of(0, b, e);
        init_thread(0);  // thread 0 runs inline: build its solver / control / Ybus copy here
        run_contingency_range(
            b, e,
            algo_for(0), controls[0], ybus_for(0), Vinit_solver,
            ac_solver_used, mask_mode, max_iter, tol, sn_mva,
            th_timer_modif[0], th_nb_solved[0], th_timer_solver[0], th_err[0], true);
    }

    for(auto & th : threads) th.join();

    // surface the first error (if any) and merge the per-thread accumulators
    for(int t = 0; t < nb_thread; ++t){
        if(th_err[t]) std::rethrow_exception(th_err[t]);
    }
    for(int t = 0; t < nb_thread; ++t){
        _timer_modif_Ybus += th_timer_modif[t];
        _nb_solved += th_nb_solved[t];
        _timer_solver += th_timer_solver[t];
    }

    _timer_total = timer.duration();
}

void ContingencyAnalysis::run_contingency_range(
    size_t cont_begin,
    size_t cont_end,
    AlgorithmSelector & algo,
    AlgoControl & control,
    Eigen::SparseMatrix<cplx_type> & Ybus,
    const CplxVect & Vinit_solver,
    bool ac_solver_used,
    bool mask_mode,
    int max_iter,
    real_type tol,
    real_type sn_mva,
    double & timer_modif_Ybus_acc,
    int & nb_solved,
    double & timer_solver,
    std::exception_ptr & err,
    bool needs_solver_init)
{
    try {
        CplxVect V;
        if(needs_solver_init){
            // Fresh per-thread solvers must analyze / factorize on their first actual solve.
            control.tell_all_changed();
        }

        for(size_t cont_id = cont_begin; cont_id < cont_end; ++cont_id){
            const std::vector<Coeff> & coeffs_modif = _li_coeffs[cont_id];
            auto timer_modif_Ybus = CustTimer();
            bool do_store = false;
            bool conv = false;

            if(mask_mode){
                // disconnected-grid handling: solve the largest component while masking
                // the rest, unless the contingency strands the chosen reference slack.
                if(!_skip_mask[cont_id]){
                    const std::vector<int> & masked = _li_masked[cont_id];
                    remove_from_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);  // invertibility handled by masking
                    timer_modif_Ybus_acc += timer_modif_Ybus.duration();

                    algo.set_masked_buses(masked);
                    const RealVect sw = masked.empty() ? slack_weights_ : masked_slack_weights(masked);
                    V = Vinit_solver; // Vinit is reused for each contingencies
                    conv = compute_one_powerflow(
                        algo, control, nb_solved, timer_solver,
                        Ybus, V, Sbus_,
                        slack_ids_me_.as_eigen(), sw,
                        bus_pv_.as_eigen(), bus_pq_.as_eigen(),
                        max_iter, tol / sn_mva);
                    if(needs_solver_init){
                        control.tell_none_changed();
                        needs_solver_init = false;
                    }
                    algo.set_masked_buses(std::vector<int>());  // reset for the next contingency

                    timer_modif_Ybus = CustTimer();
                    readd_to_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                    timer_modif_Ybus_acc += timer_modif_Ybus.duration();

                    if(conv){
                        // masked buses are not simulated: report 0 voltage for them
                        for(int b : masked) if(b >= 0 && b < V.size()) V(b) = cplx_type(0., 0.);
                        do_store = true;
                    }
                }
                // skipped contingency: conv stays false, voltages stay 0
            } else {
                // legacy behaviour: skip the contingency if it disconnects the grid
                bool invertible = remove_from_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                timer_modif_Ybus_acc += timer_modif_Ybus.duration();

                if(invertible)
                {
                    V = Vinit_solver; // Vinit is reused for each contingencies
                    conv = compute_one_powerflow(
                        algo, control, nb_solved, timer_solver,
                        Ybus,
                        V,
                        Sbus_,
                        slack_ids_me_.as_eigen(),
                        slack_weights_,
                        bus_pv_.as_eigen(),
                        bus_pq_.as_eigen(),
                        max_iter,
                        tol / sn_mva);
                    if(needs_solver_init){
                        control.tell_none_changed();
                        needs_solver_init = false;
                    }
                }

                timer_modif_Ybus = CustTimer();
                readd_to_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                timer_modif_Ybus_acc += timer_modif_Ybus.duration();
                do_store = conv && invertible;
            }

            if (do_store) _voltages.row(cont_id)(id_solver_to_me_.as_eigen()) = V.array();
        }
    } catch(...) {
        err = std::current_exception();
    }
}

// by default the flows are not 0 when the powerline is connected in the original topology
// this function sorts this out
void ContingencyAnalysis::clean_flows(bool is_amps)
{
    auto timer = CustTimer();
    size_t cont_id = 0;
    for(const auto & l_id_this_cont: _li_defaults){
        for(auto l_id : l_id_this_cont){
            real_type & el = is_amps ? _amps_flows(cont_id, l_id): _active_power_flows(cont_id, l_id);
            if(isfinite(el)) el = 0.;
        }
        ++cont_id;
    }
    if (is_amps) _timer_compute_A += timer.duration();
    else _timer_compute_P += timer.duration();
}

} // namespace ls2g
