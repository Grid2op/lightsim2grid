// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GridModel.hpp"
#include "ChooseSolver.hpp"  // to avoid circular references

#include <queue>


GridModel::GridModel(const GridModel & other) noexcept
{
    init_vm_pu_ = other.init_vm_pu_;
    sn_mva_ = other.sn_mva_;
    compute_results_ = other.compute_results_;

    // copy the powersystem representation
    // 1. bus
    last_bus_status_saved_ = other.last_bus_status_saved_;
    substations_ = other.substations_;
    max_nb_bus_per_sub_ = substations_.nmax_busbar_per_sub();
    n_sub_ = substations_.nb_sub();

    set_ls_to_orig_internal(other._ls_to_orig);  // sets also orig_to_ls

    // 2. powerline
    powerlines_ = other.powerlines_;

    // 3. shunt
    shunts_ = other.shunts_;

    // 4. transformers
    // have the r, x, h and ratio
    // ratio is computed from the tap, so maybe store tap num and tap_step_pct
    trafos_ = other.trafos_;

    // 5. generators
    total_q_min_per_bus_ = RealVect();
    total_q_max_per_bus_ = RealVect();
    generators_ = other.generators_;

    // 6. loads
    loads_ = other.loads_;

    // 7. static generators
    sgens_ = other.sgens_;

    // 8. storage units
    storages_ = other.storages_;

    // dc lines
    dc_lines_ = other.dc_lines_;

    // assign the right solver
    reset(true, true, true);
    _solver.change_solver(other.get_solver_type());
    _dc_solver.change_solver(other.get_dc_solver_type());
}

//pickle
GridModel::StateRes GridModel::get_state() const 
{
    std::vector<int> ls_to_orig(_ls_to_orig.begin(), _ls_to_orig.end());
    int version_major = VERSION_MAJOR;
    int version_medium = VERSION_MEDIUM;
    int version_minor = VERSION_MINOR;
    auto res_substation = substations_.get_state();
    auto res_line = powerlines_.get_state();
    auto res_shunt = shunts_.get_state();
    auto res_trafo = trafos_.get_state();
    auto res_gen = generators_.get_state();
    auto res_load = loads_.get_state();
    auto res_sgen = sgens_.get_state();
    auto res_storage = storages_.get_state();
    auto res_dc_line = dc_lines_.get_state();

    GridModel::StateRes res(version_major,
                            version_medium,
                            version_minor,
                            ls_to_orig,
                            init_vm_pu_,
                            sn_mva_,
                            last_bus_status_saved_,
                            res_substation,
                            res_line,
                            res_shunt,
                            res_trafo,
                            res_gen,
                            res_load,
                            res_sgen,
                            res_storage,
                            res_dc_line,
                            get_solver_type(),
                            get_dc_solver_type()
                            );
    return res;
};

void GridModel::set_state(GridModel::StateRes & my_state)
{
    // after loading back, the instance need to be reset anyway
    // TODO see if it's worth the trouble NOT to do it
    solver_control_.tell_all_changed();
    compute_results_ = true;

    // extract data from the state
    int version_major = std::get<0>(my_state);
    int version_medium = std::get<1>(my_state);
    int version_minor = std::get<2>(my_state);
    if((version_major != VERSION_MAJOR )| (version_medium != VERSION_MEDIUM) | (version_minor != VERSION_MINOR))
    {
        std::ostringstream exc_;
        exc_ << "GridModel::set_state: Wrong version. You tried to load a lightsim2grid model saved with version ";
        exc_ << version_major << "." << version_medium << "." << version_minor;
        exc_ << " while currently using the package on version ";
        exc_ << VERSION_MAJOR << "." << VERSION_MEDIUM << "." << VERSION_MINOR;
        exc_ << "It is not possible. Please reinstall it.";
        throw std::runtime_error(exc_.str());
    }
    const std::vector<int> & ls_to_pp = std::get<3>(my_state);
    init_vm_pu_ = std::get<4>(my_state);
    sn_mva_ = std::get<5>(my_state);
    // const std::vector<real_type> & bus_vn_kv = std::get<6>(my_state);
    const std::vector<bool> & last_bus_status_saved = std::get<6>(my_state);
    SubstationContainer::StateRes & state_substations = std::get<7>(my_state);
    // powerlines
    LineContainer::StateRes & state_lines = std::get<8>(my_state);
    // shunts
    ShuntContainer::StateRes & state_shunts = std::get<9>(my_state);
    // trafos
    TrafoContainer::StateRes & state_trafos = std::get<10>(my_state);
    // generators
    // total_q_min_per_bus_;
    // total_q_max_per_bus_;
    // total_gen_per_bus_;
    GeneratorContainer::StateRes & state_gens = std::get<11>(my_state);
    // loads
    LoadContainer::StateRes & state_loads = std::get<12>(my_state);
    // static gen
    SGenContainer::StateRes & state_sgens= std::get<13>(my_state);
    // storage units
    LoadContainer::StateRes & state_storages = std::get<14>(my_state);
    // dc lines
    DCLineContainer::StateRes & state_dc_lines = std::get<15>(my_state);

    // assign it to this instance
    set_ls_to_orig(IntVect::Map(ls_to_pp.data(), ls_to_pp.size()));  // set also _orig_to_ls

    // substations
    last_bus_status_saved_ = last_bus_status_saved;
    substations_.set_state(state_substations);
    max_nb_bus_per_sub_ = substations_.nmax_busbar_per_sub();
    n_sub_ = substations_.nb_sub();

    // elements
    // 1. powerlines
    powerlines_.set_state(state_lines);
    // 2. shunts
    shunts_.set_state(state_shunts);
    // 3. trafos
    trafos_.set_state(state_trafos);
    // 4. gen
    total_q_min_per_bus_ = RealVect();
    total_q_max_per_bus_ = RealVect();
    generators_.set_state(state_gens);
    // 5. loads
    loads_.set_state(state_loads);
    // 6. static generators
    sgens_.set_state(state_sgens);
    // 7. storage units
    storages_.set_state(state_storages);
    // dc lines
    dc_lines_.set_state(state_dc_lines);

    // handle the solver
    reset(true, true, true);
    _solver.change_solver(std::get<16>(my_state));
    _dc_solver.change_solver(std::get<17>(my_state));
};

void GridModel::set_ls_to_orig(const IntVect & ls_to_orig){
    if(ls_to_orig.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }

    if(ls_to_orig.size() != substations_.nb_bus()) 
        throw std::runtime_error("Impossible to set the converter ls_to_orig: the provided vector has not the same size as the number of bus on the grid.");
    set_ls_to_orig_internal(ls_to_orig);
}

void GridModel::set_orig_to_ls(const IntVect & orig_to_ls){
    if(orig_to_ls.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }
    _orig_to_ls = orig_to_ls;
    Eigen::Index nb_bus_ls = 0;
    for(const auto el : orig_to_ls){
        if (el != -1) nb_bus_ls += 1;
    }
    if(nb_bus_ls != substations_.nb_bus()) 
        throw std::runtime_error("Impossible to set the converter orig_to_ls: the number of 'non -1' component in the provided vector does not match the number of buses on the grid.");
    _ls_to_orig = IntVect::Constant(nb_bus_ls, -1);
    Eigen::Index ls2or_ind = 0;
    for(auto or2ls_ind = 0; or2ls_ind < nb_bus_ls; ++or2ls_ind){
        const auto my_ind = _orig_to_ls[or2ls_ind];
        if(my_ind >= 0){
            _ls_to_orig[ls2or_ind] = my_ind;
            ls2or_ind++;
        }
    }
}

void GridModel::set_ls_to_orig_internal(const IntVect & ls_to_orig){
    if(ls_to_orig.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }
    
    _ls_to_orig = ls_to_orig;
    const auto size = ls_to_orig.lpNorm<Eigen::Infinity>();
    _orig_to_ls = IntVect::Constant(size + 1, -1);
    int i = 0;
    for(auto el : _ls_to_orig){
        if(el != -1) _orig_to_ls[el] = i;
        ++i;
    }
}

//init
void GridModel::init_bus(unsigned int n_sub,
                         unsigned int n_busbar_per_sub,
                         const RealVect & bus_vn_kv,
                         int nb_line,
                         int nb_trafo){
    /**
    initialize the bus_vn_kv_ member
    and
    initialize the Ybus_ matrix at the proper shape
    **/
    n_sub_ = n_sub;
    max_nb_bus_per_sub_ = n_busbar_per_sub;
    substations_.init_bus(n_sub_, max_nb_bus_per_sub_,  bus_vn_kv);
    _orig_to_ls = IntVect();
    _ls_to_orig = IntVect();
}

void GridModel::reset(bool reset_solver, bool reset_ac, bool reset_dc)
{
    if(reset_ac){
        id_me_to_ac_solver_ = std::vector<SolverBusId>();
        id_ac_solver_to_me_ = std::vector<GlobalBusId>();
        slack_bus_id_ac_solver_ = SolverBusIdVect();
        Ybus_ac_ = Eigen::SparseMatrix<cplx_type>();
    }

    if(reset_dc){
        id_me_to_dc_solver_ = std::vector<SolverBusId>();
        id_dc_solver_to_me_ = std::vector<GlobalBusId>();
        slack_bus_id_dc_solver_ = SolverBusIdVect();
        Ybus_dc_ = Eigen::SparseMatrix<cplx_type>();
    }

    timer_last_ac_pf_= 0.;
    timer_last_dc_pf_ = 0.;

    acSbus_ = CplxVect();
    dcSbus_ = CplxVect();
    bus_pv_ = SolverBusIdVect();
    bus_pq_ = SolverBusIdVect();

    solver_control_.tell_all_changed();
    tell_solver_need_reset(); // also handles last_bus_status_saved_
    
    slack_bus_id_ac_me_ = GlobalBusIdVect();  // slack bus id, gridmodel number
    slack_bus_id_ac_solver_ = SolverBusIdVect();  // slack bus id, solver number
    slack_bus_id_dc_me_ = GlobalBusIdVect();
    slack_bus_id_dc_solver_ = SolverBusIdVect();
    slack_weights_ = RealVect();

    // reset the solvers
    if (reset_solver){
        _solver.reset();
        _solver.set_gridmodel(this);
        _solver.tell_solver_control(solver_control_);

        _dc_solver.reset();
        _dc_solver.set_gridmodel(this);
        _dc_solver.tell_solver_control(solver_control_);
    }
}

CplxVect GridModel::ac_pf(const CplxVect & Vinit,
                          int max_iter,
                          real_type tol)
{
    auto timer = CustTimer();
    const int nb_bus = static_cast<int>(substations_.nb_bus());
    if(Vinit.size() != nb_bus){
        std::ostringstream exc_;
        exc_ << "GridModel::ac_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    bool conv = false;
    CplxVect res = CplxVect();

    // reset_results();  // clear the results  No need to do it, results are neceassirly set or reset in post process

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    bool is_ac = true;
    CplxVect V = pre_process_solver(Vinit, 
                                    acSbus_,
                                    Ybus_ac_,
                                    id_me_to_ac_solver_,
                                    id_ac_solver_to_me_,
                                    slack_bus_id_ac_me_,
                                    slack_bus_id_ac_solver_,
                                    is_ac,
                                    solver_control_);

    // start the solver
    conv = _solver.compute_pf(
        Ybus_ac_,
        V,
        acSbus_,
        _to_intvect(slack_bus_id_ac_solver_),
        slack_weights_,
        _to_intvect(bus_pv_),
        _to_intvect(bus_pq_),
        max_iter,
        tol / sn_mva_);

    // store results (in ac mode) 
    process_results(conv, res, Vinit, true, id_me_to_ac_solver_);

    timer_last_ac_pf_ = timer.duration();
    // return the vector of complex voltage at each bus
    return res;
};

void GridModel::check_solution_q_values_onegen(CplxVect & res,
                                               const GenInfo& gen,
                                               bool check_q_limits) const{
    if(check_q_limits)
    {
        // i need to check the reactive can be absorbed / produced by the generator
        real_type new_q = BaseConstants::my_zero_;
        real_type react_this_bus = std::imag(res.coeff(gen.bus_id));
        if((react_this_bus >= gen.min_q_mvar) && (react_this_bus <= gen.max_q_mvar))
        {
            // this generator is able to handle all reactive
            new_q = BaseConstants::my_zero_;
        }else if(react_this_bus < gen.min_q_mvar){
            // generator cannot absorb enough reactive power
            new_q = react_this_bus - gen.min_q_mvar; //ex. need -50, qmin is -30, remains: (-50) - (-30) = -20 MVAr
        }else{
            // generator cannot produce enough reactive power
            new_q = react_this_bus - gen.max_q_mvar;  // ex. need 50, qmax is 30, remains: 50 - 30 = 20 MVAr
        }
        res.coeffRef(gen.bus_id) = {std::real(res.coeff(gen.bus_id)), new_q};
    }else{
        // the q value for the bus at which the generator is connected will be 0
        res.coeffRef(gen.bus_id) = {std::real(res.coeff(gen.bus_id)), BaseConstants::my_zero_};
    }
}

void GridModel::check_solution_q_values(CplxVect & res, bool check_q_limits) const{
    // test for iterator though generators
    for(const auto & gen: generators_)
    {
        if(!gen.connected)
        {
            // the generator is disconnected, I do nothing
            continue;
        }
        check_solution_q_values_onegen(res, gen, check_q_limits);

        // if(gen.id == gen_slackbus_)
        if(gen.is_slack)
        {
            // slack bus, by definition, can handle all active value
            // This is probably not the case with distributed slack !
            res.coeffRef(gen.bus_id) = {BaseConstants::my_zero_, std::imag(res.coeff(gen.bus_id))};
        }
    }

    // then do the same for dc powerlines
    for(const auto & dcline: dc_lines_)
    {
        if(!dcline.connected_global)
        {
            // the generator is disconnected, I do nothing
            continue;
        }
        check_solution_q_values_onegen(res, dcline.gen_side_1, check_q_limits);
        check_solution_q_values_onegen(res, dcline.gen_side_2, check_q_limits);
    }
}

CplxVect GridModel::check_solution(const CplxVect & V_proposed, bool check_q_limits)
{
    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    const int nb_bus = static_cast<int>(V_proposed.size());
    bool is_ac = true;
    SolverControl reset_solver;
    reset_solver.tell_none_changed();  // TODO reset solver
    CplxVect V = pre_process_solver(V_proposed, 
                                    acSbus_,
                                    Ybus_ac_, 
                                    id_me_to_ac_solver_,
                                    id_ac_solver_to_me_,
                                    slack_bus_id_ac_me_,
                                    slack_bus_id_ac_solver_,
                                    is_ac, reset_solver);

    // compute the mismatch
    CplxVect tmp = Ybus_ac_ * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - acSbus_.array();  // TODO ac or dc here

    // store results
    CplxVect res = _get_results_back_to_orig_nodes(mis,
                                                   id_me_to_ac_solver_,
                                                   static_cast<int>(V_proposed.size())
                                                   );
    if(abs(sn_mva_- 1.) > BaseConstants::_tol_equal_float) res *= sn_mva_;

    // now check reactive values for buses where there are generators and active values of slack bus
    check_solution_q_values(res, check_q_limits);

    // set to 0 the error on the disconnected bus (it is not initialized at 0.0 in _get_results_back_to_orig_nodes)
    for(int bus_id = 0; bus_id < nb_bus; ++bus_id)
    {
        if(substations_.is_bus_connected(GlobalBusId(bus_id))) continue;
        res.coeffRef(bus_id) = BaseConstants::my_zero_;
    }
    return res;
};

CplxVect GridModel::pre_process_solver(
    const CplxVect & Vinit, 
    CplxVect & Sbus,
    Eigen::SparseMatrix<cplx_type> & Ybus,
    std::vector<SolverBusId> & id_me_to_solver,
    std::vector<GlobalBusId> & id_solver_to_me,
    GlobalBusIdVect & slack_bus_id_me,
    SolverBusIdVect & slack_bus_id_solver,
    bool is_ac,
    const SolverControl & solver_control)
{
    // TODO get rid of the "is_ac" argument: this info is available in the _solver already
    if(is_ac){
        if(solver_control.need_reset_solver()){   
            _solver.reset();
        }
    } else {
        if(solver_control.need_reset_solver()){
            _dc_solver.reset();
        }
    }

    bool redo_all = 
            solver_control.need_reset_solver() || 
            solver_control.has_dimension_changed();

    if (redo_all ||
        solver_control.has_slack_participate_changed()){
            slack_bus_id_me = generators_.get_slack_bus_id();
            // this is the slack bus ids with the gridmodel ordering, not the solver ordering.
            // conversion to solver ordering is done in init_slack_bus
        }
    if (redo_all || solver_control.has_one_el_changed_bus()){
        init_bus_status();
    }
    
    // init_bus_status can set the flag "has_dimension_change"
    // so I need to redo this here
    redo_all = 
            solver_control.need_reset_solver() || 
            solver_control.has_dimension_changed();
    bool converter_changed = false;
    if (redo_all ||
        solver_control.ybus_change_sparsity_pattern()){
            init_converter_bus_id(id_me_to_solver, id_solver_to_me);
            const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
            init_Ybus(Ybus, nb_bus_solver);
            converter_changed = true;
        }
    if (redo_all ||
        converter_changed || 
        solver_control.need_recompute_ybus()){
            fillYbus(Ybus, is_ac, id_me_to_solver);
        }
    if (redo_all || converter_changed || solver_control.need_recompute_sbus()) {
            // init Sbus
            Sbus = CplxVect::Constant(id_solver_to_me.size(), 0.);
        }
    if (redo_all || converter_changed ||
        solver_control.has_slack_participate_changed() || 
        solver_control.has_pv_changed() || 
        solver_control.has_pq_changed()) {
            init_slack_bus(Sbus, id_me_to_solver, id_solver_to_me, slack_bus_id_me, slack_bus_id_solver);
            fillpv_pq(
                id_me_to_solver,
                id_solver_to_me,
                slack_bus_id_solver,
                solver_control);
        }
    
    if (is_ac && (redo_all ||
                  solver_control.need_recompute_sbus() ||  // TODO do we need it ?
                  solver_control.has_slack_participate_changed() || 
                  solver_control.has_pv_changed() || 
                  solver_control.has_pq_changed())  // TODO do we need it ?
        ){
        int nb_bus_total = static_cast<int>(substations_.nb_bus());
        total_q_min_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
        total_q_max_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
        total_gen_per_bus_ = Eigen::VectorXi::Constant(nb_bus_total, 0);
        generators_.init_q_vector(nb_bus_total, total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
        dc_lines_.init_q_vector(nb_bus_total, total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
    }

    if (redo_all || converter_changed ||
        solver_control.has_slack_participate_changed() || 
        solver_control.has_pv_changed() || 
        solver_control.has_pq_changed() ||
        solver_control.need_recompute_sbus()) {
            fillSbus_me(Sbus, is_ac, id_me_to_solver);
        }
    
    const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
    CplxVect V = CplxVect::Constant(nb_bus_solver, init_vm_pu_);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        GlobalBusId bus_me_id = id_solver_to_me[bus_solver_id]; 
        if(bus_me_id.cast_int() == BaseConstants::_deactivated_bus_id){
            //TODO DEBUG MODE : only in debug mode
            std::ostringstream exc_;
            exc_ << "GridModel::pre_process_solver: the bus with solver id ";
            exc_ << bus_solver_id;
            exc_ << " is connected, but mapped (in id_solver_to_me) to a disconnected bus (global / gridmodel id)";
            throw std::runtime_error(exc_.str());
        }
        cplx_type tmp = Vinit(bus_me_id.cast_int());
        V(bus_solver_id) = tmp;
    }
    generators_.set_vm(V, id_me_to_solver);
    dc_lines_.set_vm(V, id_me_to_solver);

    if(solver_control_.need_reset_solver() || 
       solver_control_.has_dimension_changed() ||
       solver_control_.has_slack_participate_changed() || 
       solver_control_.has_pv_changed() || 
       solver_control_.has_slack_weight_changed()){
        slack_weights_ = generators_.get_slack_weights_solver(Ybus.rows(), id_me_to_solver); 
    }

    if(is_ac) _solver.tell_solver_control(solver_control_);
    else _dc_solver.tell_solver_control(solver_control_);
    return V;
}

CplxVect GridModel::_get_results_back_to_orig_nodes(const CplxVect & res_tmp, 
                                                    std::vector<SolverBusId> & id_me_to_solver,
                                                    int size)
{
    CplxVect res = CplxVect::Constant(size, {init_vm_pu_, BaseConstants::my_zero_});
    const int nb_bus = static_cast<int>(substations_.nb_bus());
    for (int bus_id_me=0; bus_id_me < nb_bus; ++bus_id_me){
        if(!substations_.is_bus_connected(GlobalBusId(bus_id_me))) continue;  // nothing is done if the bus is connected
        SolverBusId bus_id_solver = id_me_to_solver[bus_id_me];
        if(bus_id_solver.cast_int() == BaseConstants::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "GridModel::_get_results_back_to_orig_nodes: the bus with id ";
            exc_ << bus_id_me;
            exc_ << " is connected to a disconnected bus (solver side)";
            throw std::runtime_error(exc_.str());
        }
        res(bus_id_me) = res_tmp(static_cast<int>(bus_id_solver));
    }
    return res;
}

void GridModel::process_results(bool conv,
                                CplxVect & res,
                                const CplxVect & Vinit,
                                bool ac,
                                std::vector<SolverBusId> & id_me_to_solver)
{
    if (conv){
        if(compute_results_){
            // compute the results of the flows, P,Q,V of loads etc.
            compute_results(ac);
        }
        // solver_control_.tell_none_changed();  // todo automatically set for ac / dc the `tell_none_changed()`
        const CplxVect & res_tmp = ac ? _solver.get_V(): _dc_solver.get_V() ;

        // convert back the results to "big" vector
        res = _get_results_back_to_orig_nodes(res_tmp,
                                              id_me_to_solver,
                                              static_cast<int>(Vinit.size()));
    } else {
        //powerflow diverge
        reset_results();
        // TODO solver control ??? something to do here ?
    }
}

void GridModel::init_converter_bus_id(std::vector<SolverBusId>& id_me_to_solver,
                                      std::vector<GlobalBusId>& id_solver_to_me){

    //TODO get disconnected bus !!! (and have some conversion for it)
    //1. init the conversion bus
    const int nb_bus_init = static_cast<int>(substations_.nb_bus());
    id_me_to_solver = std::vector<SolverBusId>(nb_bus_init, SolverBusId(BaseConstants::_deactivated_bus_id));  // by default, if a bus is disconnected, then it has a -1 there
    id_solver_to_me = std::vector<GlobalBusId>();
    id_solver_to_me.reserve(nb_bus_init);
    int bus_id_solver = 0;
    for(int bus_id_me=0; bus_id_me < nb_bus_init; ++bus_id_me){
        if(substations_.is_bus_connected(GlobalBusId(bus_id_me))){
            // bus is connected
            id_solver_to_me.push_back(GlobalBusId(bus_id_me));
            id_me_to_solver[bus_id_me] = SolverBusId(bus_id_solver);
            ++bus_id_solver;
        }
    }
}

void GridModel::init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                          int nb_bus_solver){
    Ybus = Eigen::SparseMatrix<cplx_type>(nb_bus_solver, nb_bus_solver);
    Ybus.reserve(nb_bus_solver + 4*powerlines_.nb() + 4*trafos_.nb() + 2 * shunts_.nb());
}

void GridModel::init_slack_bus(const CplxVect & Sbus,
                               const std::vector<SolverBusId>& id_me_to_solver,
                               const std::vector<GlobalBusId>& id_solver_to_me,
                               const GlobalBusIdVect & slack_bus_id_me,
                               SolverBusIdVect & slack_bus_id_solver)
{
    slack_bus_id_solver = SolverBusIdVect::Constant(slack_bus_id_me.size(), SolverBusId(BaseConstants::_deactivated_bus_id));

    size_t i = 0;
    for(const GlobalBusId & el: slack_bus_id_me) {
        SolverBusId tmp = id_me_to_solver[el.cast_int()];
        if(tmp.cast_int() == BaseConstants::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "GridModel::init_slack_bus: One of the slack bus is disconnected.";
            exc_ << " You can check bus with global id GlobalBusId : ";
            exc_ << el.cast_int();
            exc_ << ": [";
            for(const auto & el2 : slack_bus_id_me) exc_ << el2.cast_int() << ", ";
            exc_ << "].";
            throw std::out_of_range(exc_.str());
        }
        slack_bus_id_solver(i) = tmp;
        ++i;
    }
    
    if(GenericContainer::is_in_vect(BaseConstants::_deactivated_bus_id, slack_bus_id_solver)){
        // TODO improve error message with the gen_id
        // TODO DEBUG MODE: only check that in debug mode
        throw std::runtime_error("GridModel::init_Sbus: One of the slack bus is disconnected !");
    }
}
void GridModel::fillYbus(
    Eigen::SparseMatrix<cplx_type> & res,
    bool ac,
    const std::vector<SolverBusId>& id_me_to_solver){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/

    // init the Ybus matrix
    res.setZero();  // it should not be needed but might not hurt too much either.
    std::vector<Eigen::Triplet<cplx_type> > tripletList;
    tripletList.reserve(substations_.nb_bus() + 4*powerlines_.nb() + 4*trafos_.nb() + shunts_.nb());
    powerlines_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    trafos_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    loads_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    sgens_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    storages_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    generators_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    dc_lines_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    res.setFromTriplets(tripletList.begin(), tripletList.end());  // works because  "The initial contents of *this is destroyed"
    res.makeCompressed();
}

void GridModel::fillSbus_me(CplxVect & Sbus, bool ac, const std::vector<SolverBusId>& id_me_to_solver)
{
    // init the Sbus 
    Sbus.array() = 0.;  // reset to 0.
    powerlines_.fillSbus(Sbus, id_me_to_solver, ac);  // TODO have a function to dispatch that to all type of elements
    trafos_.fillSbus(Sbus, id_me_to_solver, ac);
    shunts_.fillSbus(Sbus, id_me_to_solver, ac);
    loads_.fillSbus(Sbus, id_me_to_solver, ac);
    sgens_.fillSbus(Sbus, id_me_to_solver, ac);
    storages_.fillSbus(Sbus, id_me_to_solver, ac);
    generators_.fillSbus(Sbus, id_me_to_solver, ac);
    dc_lines_.fillSbus(Sbus, id_me_to_solver, ac);
    if (abs(sn_mva_ - 1.0) > BaseConstants::_tol_equal_float) Sbus /= sn_mva_;
    // in dc mode, this is used for the phase shifter, this should not be divided by sn_mva_ !
    trafos_.hack_Sbus_for_dc_phase_shifter(Sbus, ac, id_me_to_solver);
}

void GridModel::fillpv_pq(const std::vector<SolverBusId>& id_me_to_solver,
                          const std::vector<GlobalBusId>& id_solver_to_me,
                          const SolverBusIdVect & slack_bus_id_solver,
                          const SolverControl & solver_control)
{
    // Nothing to do if neither pv, nor pq nor the dimension of the problem has changed

    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    const int nb_bus = static_cast<int>(id_solver_to_me.size());  // number of bus in the solver!
    std::vector<int> bus_pq;
    bus_pq.reserve(nb_bus);
    std::vector<int> bus_pv;
    bus_pv.reserve(nb_bus);
    std::vector<bool> has_bus_been_added(nb_bus, false);

    bus_pv_ = SolverBusIdVect();
    bus_pq_ = SolverBusIdVect();
    powerlines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    trafos_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    loads_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    storages_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    sgens_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    generators_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    dc_lines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);

    for(int bus_id = 0; bus_id< nb_bus; ++bus_id){
        if(GenericContainer::is_in_vect(bus_id, slack_bus_id_solver)) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_ = SolverBusIdVect(bus_pv.size());
    for(int i = 0; i < static_cast<int>(bus_pv.size()); ++i){
        bus_pv_(i) = SolverBusId(bus_pv[i]);
    }
    bus_pq_ = SolverBusIdVect(bus_pq.size());
    for(int i = 0; i< static_cast<int>(bus_pq.size()); ++i){
        bus_pq_(i) = SolverBusId(bus_pq[i]);
    }
}

void GridModel::compute_results(bool ac){
    // retrieve results from powerflow
    const auto & Va = ac ? _solver.get_Va() : _dc_solver.get_Va();
    const auto & Vm = ac ? _solver.get_Vm() : _dc_solver.get_Vm();
    const auto & V = ac ? _solver.get_V() : _dc_solver.get_V();

    const std::vector<SolverBusId> & id_me_to_solver = ac ? id_me_to_ac_solver_ : id_me_to_dc_solver_;
    // for powerlines
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);  // TODO have a function to dispatch that to all type of elements
    // for trafo
    trafos_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for loads
    loads_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for static gen
    sgens_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for storage units
    storages_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for prods
    generators_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for dclines
    dc_lines_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);

    //handle_slack_bus active power
    CplxVect mismatch;  // power mismatch at each bus (SOLVER BUS !!!)
    RealVect reactive_mismatch;  // not used in dc mode (DO NOT ATTEMPT TO USE IT THERE)
    RealVect active_mismatch;
    if(ac){
        // In AC mode i am not forced to run through all the grid
        // auto tmp = (Ybus_ac_ * V).conjugate();
        mismatch = V.array() * (Ybus_ac_ * V).conjugate().array() - acSbus_.array();
        active_mismatch = mismatch.real() * sn_mva_;
    } else{
        active_mismatch = RealVect::Zero(V.size());
        //TODO SLACK: improve distributed slack for DC mode !
        // it is possible to know in advance the contribution of each slack generators (sum(Sbus) MW 
        // to split among the contributing generators) so it's possible to "mess with" Sbus 
        // for such purpose
        const SolverBusId id_slack = slack_bus_id_dc_solver_(0);
        active_mismatch(id_slack.cast_int()) = -dcSbus_.real().sum() * sn_mva_;
    }
    generators_.set_p_slack(active_mismatch, id_me_to_solver);

    if(ac) reactive_mismatch = mismatch.imag() * sn_mva_;
    // mainly to initialize the Q value of the generators in dc (just fill it with 0.)
    generators_.set_q(reactive_mismatch, id_me_to_solver, ac,
                      total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
    dc_lines_.set_q(reactive_mismatch, id_me_to_solver, ac,
                    total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
}

void GridModel::reset_results(){
    powerlines_.reset_results();  // TODO have a function to dispatch that to all type of elements
    shunts_.reset_results();
    trafos_.reset_results();
    loads_.reset_results();
    sgens_.reset_results();
    storages_.reset_results();
    generators_.reset_results();
    dc_lines_.reset_results();
}

CplxVect GridModel::dc_pf(const CplxVect & Vinit,
                          int max_iter,  // not used for DC
                          real_type tol  // not used for DC
                          )
{
    //TODO SLACK: improve distributed slack for DC mode !
    // the idea is to "mess" with the Sbus beforehand to split the "losses"
    // ie fake the action of generators to adjust Sbus such that sum(Sbus) = 0
    // and the slack contribution factors are met.
    auto timer = CustTimer();

    const int nb_bus = static_cast<int>(substations_.nb_bus());
    if(Vinit.size() != nb_bus){
        //TODO DEBUG MODE: 
        std::ostringstream exc_;
        exc_ << "GridModel::dc_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    bool conv = false;
    CplxVect res = CplxVect();

    // reset_results();  // clear the results  No need to do it, results are neceassirly set or reset in post process

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    bool is_ac = false;
    CplxVect V = pre_process_solver(Vinit,
                                    dcSbus_,
                                    Ybus_dc_,
                                    id_me_to_dc_solver_,
                                    id_dc_solver_to_me_,
                                    slack_bus_id_dc_me_,
                                    slack_bus_id_dc_solver_,
                                    is_ac,
                                    solver_control_);
    // start the solver
    conv = _dc_solver.compute_pf(
        Ybus_dc_,
        V,
        dcSbus_,
        _to_intvect(slack_bus_id_dc_solver_),
        slack_weights_,
        _to_intvect(bus_pv_),
        _to_intvect(bus_pq_),
        max_iter,
        tol);
    // store results (fase -> because I am in dc mode)
    process_results(conv, res, Vinit, is_ac, id_me_to_dc_solver_);
    timer_last_dc_pf_ = timer.duration();
    return res;
}

RealMat GridModel::get_ptdf_solver(){
    if(Ybus_dc_.size() == 0){
        throw std::runtime_error("GridModel::get_ptdf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    const RealMat & PTDF_solver = _dc_solver.get_ptdf();
    return PTDF_solver;
}


RealMat GridModel::get_ptdf(){
    if(Ybus_dc_.size() == 0){
        throw std::runtime_error("GridModel::get_ptdf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    const RealMat & PTDF_solver = get_ptdf_solver();
    RealMat PTDF_grid =  RealMat::Zero(powerlines_.nb() + trafos_.nb(), total_bus());  // , std::numeric_limits<real_type>::quiet_NaN()
    int solver_col = 0;
    for(const GlobalBusId & my_col: id_dc_solver_to_me()){
        PTDF_grid.col(my_col.cast_int()) = PTDF_solver.col(solver_col);
        ++solver_col;
    }
    return PTDF_grid;
}

RealMat GridModel::get_lodf(){
    if(Ybus_dc_.size() == 0){
        throw std::runtime_error("GridModel::get_lodf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    const int nb_el = powerlines_.nb() + trafos_.nb();
    GlobalBusIdVect from_bus(nb_el);
    GlobalBusIdVect to_bus(nb_el);
    // retrieve the from_bus / to_bus from the grid
    from_bus << powerlines_.get_bus_id_side_1(), trafos_.get_bus_id_side_1();
    to_bus << powerlines_.get_bus_id_side_2(), trafos_.get_bus_id_side_2();

    // convert it to solver bus id
    IntVect from_bus_solver(nb_el);  // TODO : SolverBusIdVect here
    IntVect to_bus_solver(nb_el);
    for(int el_id = 0; el_id < nb_el; ++el_id){
        // from side
        GlobalBusId f_grid_bus = from_bus[el_id];
        SolverBusId f_solver_bus = id_me_to_dc_solver_[f_grid_bus.cast_int()];
        from_bus_solver[el_id] = f_solver_bus.cast_int();
        // to side
        GlobalBusId t_grid_bus = to_bus[el_id];
        SolverBusId t_solver_bus = id_me_to_dc_solver_[t_grid_bus.cast_int()];
        to_bus_solver[el_id] = t_solver_bus.cast_int();
    }
    return _dc_solver.get_lodf(from_bus_solver, to_bus_solver);
}

Eigen::SparseMatrix<real_type> GridModel::get_Bf_solver(){
    if(Ybus_dc_.size() == 0){
        throw std::runtime_error("GridModel::get_Bf_solver: Cannot get the Bf matrix without having first computed a DC powerflow.");
    }
    Eigen::SparseMatrix<real_type> Bf;
    fillBf_for_PTDF(Bf);
    return Bf;
}

Eigen::SparseMatrix<real_type> GridModel::get_Bf(){
    if(Ybus_dc_.size() == 0){
        throw std::runtime_error("GridModel::get_Bf: Cannot get the Bf matrix without having first computed a DC powerflow.");
    }
    Eigen::SparseMatrix<real_type> Bf_solver = get_Bf_solver();
    return _relabel_matrix(Bf_solver, id_dc_solver_to_me_, false);
}

void GridModel::add_gen_slackbus(int gen_id, real_type weight){
    if(gen_id < 0)
    {
        std::ostringstream exc_;
        exc_ << "GridModel::add_gen_slackbus: Slack bus should be an id of a generator, thus positive. You provided: ";
        exc_ << gen_id;
        throw std::runtime_error(exc_.str());
    }
    if(gen_id >= generators_.nb())
    {
        std::ostringstream exc_;
        exc_ << "GridModel::add_gen_slackbus: There are only " << generators_.nb() << " generators on the grid. ";
        exc_ << "Generator with id " << gen_id << " does not exist and can't be the slack bus";
        throw std::runtime_error(exc_.str());
    }
    if(weight <= 0.){
        std::ostringstream exc_;
        exc_ << "GridModel::add_gen_slackbus: please enter a valid weight for the slack bus (> 0.)";
        throw std::runtime_error(exc_.str());
    }
    generators_.add_slackbus(gen_id, weight, solver_control_);
}

void GridModel::remove_gen_slackbus(int gen_id){
    if(gen_id < 0)
    {
        // TODO DEBUG MODE: only check when in debug mode
        std::ostringstream exc_;
        exc_ << "GridModel::remove_gen_slackbus: Slack bus should be an id of a generator, thus positive. You provided: ";
        exc_ << gen_id;
        throw std::runtime_error(exc_.str());
    }
    if(gen_id >= generators_.nb())
    {
        // TODO DEBUG MODE: only check when in debug mode
        std::ostringstream exc_;
        exc_ << "GridModel::remove_gen_slackbus: There are only " << generators_.nb() << " generators on the grid. ";
        exc_ << "Generator with id " << gen_id << " does not exist and can't be the slack bus";
        throw std::runtime_error(exc_.str());
    }
    generators_.remove_slackbus(gen_id, solver_control_);
}

/** GRID2OP SPECIFIC REPRESENTATION **/
void GridModel::update_gens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_p_gen);
}

void GridModel::update_sgens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_p_sgen);
}

void GridModel::update_gens_v(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_v_gen);
}

void GridModel::update_loads_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_p_load);
}

void GridModel::update_loads_q(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_q_load);
}

void GridModel::update_storages_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_p_storage);
}

void GridModel::update_topo(Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<const Eigen::Array<int,  Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    loads_.update_topo(has_changed, new_values, solver_control_, substations_);
    generators_.update_topo(has_changed, new_values, solver_control_, substations_);
    storages_.update_topo(has_changed, new_values, solver_control_, substations_);
    // shunts are not in "topo" in grid2op

    // NB we suppose that if a powerline (or a trafo) is disconnected, then both its ends are
    // and same for trafo, obviously
    powerlines_.update_topo(has_changed, new_values, solver_control_, substations_);
    trafos_.update_topo(has_changed, new_values, solver_control_, substations_);
}

// for FDPF (implementation of the alg 2 method FDBX (FDXB will follow)  // TODO FDPF
void GridModel::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, 
                           Eigen::SparseMatrix<real_type> & Bpp, 
                           FDPFMethod xb_or_bx) const
{
    // clear the matrices
    const int nb_bus_solver = static_cast<int>(id_ac_solver_to_me_.size());
    Bp = Eigen::SparseMatrix<real_type>(nb_bus_solver, nb_bus_solver);
    Bpp = Eigen::SparseMatrix<real_type>(nb_bus_solver, nb_bus_solver);

    // init the Bp and Bpp matrices for Fast Decoupled Powerflow  (TODO FDPF: optim when it's NOT needed just like for Ybus)
    std::vector<Eigen::Triplet<real_type> > tripletList_Bp;
    std::vector<Eigen::Triplet<real_type> > tripletList_Bpp;
    tripletList_Bp.reserve(substations_.nb_bus() + 4 * powerlines_.nb() + 4 * trafos_.nb() + shunts_.nb());
    tripletList_Bpp.reserve(substations_.nb_bus() + 4 * powerlines_.nb() + 4 * trafos_.nb() + shunts_.nb());
    // run through the grid and get the parameters to fill them
    powerlines_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    trafos_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    loads_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    sgens_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    storages_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    generators_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    dc_lines_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    // now make the matrices effectively
    Bp.setFromTriplets(tripletList_Bp.begin(), tripletList_Bp.end());
    Bp.makeCompressed();
    Bpp.setFromTriplets(tripletList_Bpp.begin(), tripletList_Bpp.end());
    Bpp.makeCompressed();
}


void GridModel::fillBf_for_PTDF(Eigen::SparseMatrix<real_type> & Bf, bool transpose) const
{
    const int nb_bus_solver = static_cast<int>(id_dc_solver_to_me_.size());
    // TODO DEBUG MODE
    if(nb_bus_solver == 0) throw std::runtime_error("GridModel::fillBf_for_PTDF: it appears no DC powerflow has run on your grid.");
    
    if(transpose){
        Bf = Eigen::SparseMatrix<real_type>(nb_bus_solver, powerlines_.nb() + trafos_.nb());
    }else{
        Bf = Eigen::SparseMatrix<real_type>(powerlines_.nb() + trafos_.nb(), nb_bus_solver);
    }
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(substations_.nb_bus() + 2 * powerlines_.nb() + 2 * trafos_.nb());
    
    powerlines_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    trafos_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    loads_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    sgens_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    storages_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    generators_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    dc_lines_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);

    Bf.setFromTriplets(tripletList.begin(), tripletList.end());
    Bf.makeCompressed();
}

// returns only the gen_id with the highest p that is connected to this bus !
// returns bus_id, gen_bus_id
std::tuple<int, int> GridModel::assign_slack_to_most_connected(){
    auto res = std::tuple<int, int>(-1, -1);
    int res_bus_id = -1;
    int res_gen_id = -1;
    int max_line = -1;
    const unsigned int nb_busbars = substations_.nb_bus();
    std::vector<real_type> gen_p_per_bus(nb_busbars, 0.);
    std::vector<int> nb_line_end_per_bus(nb_busbars, 0);

    // computes the total amount of power produce at each nodes
    powerlines_.gen_p_per_bus(gen_p_per_bus);  // TODO have a function to dispatch that to all type of elements
    shunts_.gen_p_per_bus(gen_p_per_bus);
    trafos_.gen_p_per_bus(gen_p_per_bus);
    loads_.gen_p_per_bus(gen_p_per_bus);
    sgens_.gen_p_per_bus(gen_p_per_bus);
    storages_.gen_p_per_bus(gen_p_per_bus);
    generators_.gen_p_per_bus(gen_p_per_bus);
    dc_lines_.gen_p_per_bus(gen_p_per_bus);

    // computes the total number of "neighbors" (extremity of connected powerlines and trafo, not real neighbors)
    powerlines_.nb_line_end(nb_line_end_per_bus);  // TODO have a function to dispatch that to all type of elements
    shunts_.nb_line_end(nb_line_end_per_bus);
    trafos_.nb_line_end(nb_line_end_per_bus);
    loads_.nb_line_end(nb_line_end_per_bus);
    sgens_.nb_line_end(nb_line_end_per_bus);
    storages_.nb_line_end(nb_line_end_per_bus);
    generators_.nb_line_end(nb_line_end_per_bus);
    dc_lines_.nb_line_end(nb_line_end_per_bus);
    
    // now find the most connected buses
    for(unsigned int bus_id = 0; bus_id < nb_busbars; ++bus_id)
    {
        const auto & nb_lines_this = nb_line_end_per_bus[bus_id];
        if((nb_lines_this > max_line) && (gen_p_per_bus[bus_id] > 0.)){
            res_bus_id = bus_id;
            max_line = nb_lines_this;
        }
    }
    // TODO DEBUG MODE
    if(res_bus_id == -1) throw std::runtime_error("GridModel::assign_slack_to_most_connected: impossible to find anything connected to a node.");
    std::get<0>(res) = res_bus_id;

    // and reset the slack bus
    generators_.remove_all_slackbus();
    res_gen_id = generators_.assign_slack_bus(res_bus_id, gen_p_per_bus, solver_control_);
    std::get<1>(res) = res_gen_id;
    slack_bus_id_ac_solver_ = SolverBusIdVect();
    slack_bus_id_dc_solver_ = SolverBusIdVect();
    slack_weights_ = RealVect();
    return res;
}

// TODO DC LINE: one side might be in the connected comp and not the other !
void GridModel::consider_only_main_component(){
    const auto & slack_buses_id = generators_.get_slack_bus_id();

    // TODO DEBUG MODE
    if(slack_buses_id.size() == 0) throw std::runtime_error("GridModel::consider_only_main_component: no slack is defined on your grid. This function cannot be used.");
    
    // build the graph
    const auto nb_busbars = substations_.nb_bus();
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(2 * powerlines_.nb() + 2 * trafos_.nb());
    powerlines_.get_graph(tripletList);  // TODO have a function to dispatch that to all type of elements
    shunts_.get_graph(tripletList);
    trafos_.get_graph(tripletList);
    loads_.get_graph(tripletList);
    sgens_.get_graph(tripletList);
    storages_.get_graph(tripletList);
    generators_.get_graph(tripletList);
    dc_lines_.get_graph(tripletList);
    Eigen::SparseMatrix<real_type> graph = Eigen::SparseMatrix<real_type>(nb_busbars, nb_busbars);
    graph.setFromTriplets(tripletList.begin(), tripletList.end());
    graph.makeCompressed();

    // find the connected buses
    // TODO copy paste from SecurityAnalysis
    std::vector<bool> tmp_visited(nb_busbars, false);
    std::vector<int> conn_comp(nb_busbars, -1);
    std::vector<bool> already_added(nb_busbars, false);

    int connected_comp = 0;
    std::queue<GlobalBusId> neighborhood;
    while(true)
    {
        neighborhood = std::queue<GlobalBusId>();

        // choose bus id (one of the slack) to start
        bool one_added = false;
        for(const auto & el : slack_buses_id){
            if(!tmp_visited[el.cast_int()] && !already_added[el.cast_int()])
            {
                one_added = true;
                neighborhood.push(el);
                break;
            }
        }
        
        if(!one_added) break; // no more slack bus, I stop

        // start the bfs
        while (true)
        {
            const GlobalBusId col_id = neighborhood.front();
            neighborhood.pop();
            tmp_visited[col_id.cast_int()] = true;
            conn_comp[col_id.cast_int()] = connected_comp;
            for (Eigen::SparseMatrix<real_type>::InnerIterator it(graph, col_id.cast_int()); it; ++it)
            {
                // add in the queue all my neighbor
                if(!tmp_visited[it.row()] && !already_added[it.row()]){
                    neighborhood.push(GlobalBusId(it.row()));
                    already_added[it.row()] = true;
                }
            }
            if(neighborhood.empty()) break;  // no more neighbors
        }

        // go to the next connected comp
        ++connected_comp;
    }

    // TODO speed optim: if connected_comp == 1 => don't do the following 2 steps
    
    // find the connected comp with the most buses
    int main_cc_id = -1;
    std::vector<int> nb_bus_per_cc(connected_comp, 0);
    for(const auto el : conn_comp){
        if(el == -1) continue;
        nb_bus_per_cc[el] += 1;
    }
    main_cc_id = std::distance(nb_bus_per_cc.begin(),
                               std::max_element(nb_bus_per_cc.begin(), nb_bus_per_cc.end()));

    // mark as visited the element in this cc
    std::vector<bool> bus_in_main_cc(nb_busbars, false);
    for(unsigned int bus_id = 0; bus_id < nb_busbars; ++bus_id){
        if(conn_comp[bus_id] == main_cc_id) bus_in_main_cc[bus_id] = true;
    }

    // disconnected elements not in main component
    powerlines_.disconnect_if_not_in_main_component(bus_in_main_cc);
    shunts_.disconnect_if_not_in_main_component(bus_in_main_cc);
    trafos_.disconnect_if_not_in_main_component(bus_in_main_cc);
    loads_.disconnect_if_not_in_main_component(bus_in_main_cc);
    sgens_.disconnect_if_not_in_main_component(bus_in_main_cc);
    storages_.disconnect_if_not_in_main_component(bus_in_main_cc);
    generators_.disconnect_if_not_in_main_component(bus_in_main_cc);
    dc_lines_.disconnect_if_not_in_main_component(bus_in_main_cc);
    // and finally deal with the buses
    init_bus_status();
}
