// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GridModel.h"

GridModel::GridModel(const GridModel & other)
{
    reset(true, true, true);

    init_vm_pu_ = other.init_vm_pu_;
    sn_mva_ = other.sn_mva_;

    // copy the powersystem representation
    // 1. bus
    bus_vn_kv_ = other.bus_vn_kv_;
    bus_status_ = other.bus_status_;

    // 2. powerline
    powerlines_ = other.powerlines_;

    // 3. shunt
    shunts_ = other.shunts_;

    // 4. transformers
    // have the r, x, h and ratio
    // ratio is computed from the tap, so maybe store tap num and tap_step_pct
    trafos_ = other.trafos_;

    // 5. generators
    generators_ = other.generators_;

    // 6. loads
    loads_ = other.loads_;

    // 7. static generators
    sgens_ = other.sgens_;

    // 8. storage units
    storages_ = other.storages_;

    // copy the attributes specific grid2op (speed optimization)
    n_sub_ = other.n_sub_;
    load_pos_topo_vect_ = other.load_pos_topo_vect_;
    gen_pos_topo_vect_ = other.gen_pos_topo_vect_;
    line_or_pos_topo_vect_ = other.line_or_pos_topo_vect_;
    line_ex_pos_topo_vect_ = other.line_ex_pos_topo_vect_;
    trafo_hv_pos_topo_vect_ = other.trafo_hv_pos_topo_vect_;
    trafo_lv_pos_topo_vect_ = other.trafo_lv_pos_topo_vect_;
    storage_pos_topo_vect_ = other.storage_pos_topo_vect_;

    load_to_subid_ = other.load_to_subid_;
    gen_to_subid_ = other.gen_to_subid_;
    line_or_to_subid_ = other.line_or_to_subid_;
    line_ex_to_subid_ = other.line_ex_to_subid_;
    trafo_hv_to_subid_ = other.trafo_hv_to_subid_;
    trafo_lv_to_subid_ = other.trafo_lv_to_subid_;
    storage_to_subid_ = other.storage_to_subid_;

    // assign the right solver
    _solver.change_solver(other._solver.get_type());
    _dc_solver.change_solver(other._dc_solver.get_type());
    compute_results_ = other.compute_results_;
}

//pickle
GridModel::StateRes GridModel::get_state() const
{
    std::vector<real_type> bus_vn_kv(bus_vn_kv_.begin(), bus_vn_kv_.end());
    int version_major = VERSION_MAJOR;
    int version_medium = VERSION_MEDIUM;
    int version_minor = VERSION_MINOR;
    auto res_line = powerlines_.get_state();
    auto res_shunt = shunts_.get_state();
    auto res_trafo = trafos_.get_state();
    auto res_gen = generators_.get_state();
    auto res_load = loads_.get_state();
    auto res_sgen = sgens_.get_state();
    auto res_storage = storages_.get_state();

    GridModel::StateRes res(version_major,
                            version_medium,
                            version_minor,
                            init_vm_pu_,
                            sn_mva_,
                            bus_vn_kv,
                            bus_status_,
                            res_line,
                            res_shunt,
                            res_trafo,
                            res_gen,
                            res_load,
                            res_sgen,
                            res_storage
                            );
    return res;
};

void GridModel::set_state(GridModel::StateRes & my_state)
{
    // after loading back, the instance need to be reset anyway
    // TODO see if it's worth the trouble NOT to do it
    reset(true, true, true);
    need_reset_ = true;
    compute_results_ = true;
    topo_changed_ = true;

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
    init_vm_pu_ = std::get<3>(my_state);
    sn_mva_ = std::get<4>(my_state);
    std::vector<real_type> & bus_vn_kv = std::get<5>(my_state);
    std::vector<bool> & bus_status = std::get<6>(my_state);

    // powerlines
    DataLine::StateRes & state_lines = std::get<7>(my_state);
    // shunts
    DataShunt::StateRes & state_shunts = std::get<8>(my_state);
    // trafos
    DataTrafo::StateRes & state_trafos = std::get<9>(my_state);
    // generators
    DataGen::StateRes & state_gens = std::get<10>(my_state);
    // loads
    DataLoad::StateRes & state_loads = std::get<11>(my_state);
    // static gen
    DataSGen::StateRes & state_sgens= std::get<12>(my_state);
    // storage units
    DataLoad::StateRes & state_storages = std::get<13>(my_state);

    // assign it to this instance

    // buses
    // 1. bus_vn_kv_
    bus_vn_kv_ = RealVect::Map(&bus_vn_kv[0], bus_vn_kv.size());
    // 2. bus status
    bus_status_ = bus_status;

    // elements
    // 1. powerlines
    powerlines_.set_state(state_lines);
    // 2. shunts
    shunts_.set_state(state_shunts);
    // 3. trafos
    trafos_.set_state(state_trafos);
    // 4. gen
    generators_.set_state(state_gens);
    // 5. loads
    loads_.set_state(state_loads);
    // 6. static generators
    sgens_.set_state(state_sgens);
    // 7. storage units
    storages_.set_state(state_storages);
};

//init
void GridModel::init_bus(const RealVect & bus_vn_kv, int nb_line, int nb_trafo){
    /**
    initialize the bus_vn_kv_ member
    and
    initialize the Ybus_ matrix at the proper shape
    **/
    const int nb_bus = static_cast<int>(bus_vn_kv.size());
    bus_vn_kv_ = bus_vn_kv;  // base_kv

    bus_status_ = std::vector<bool>(nb_bus, true); // by default everything is connected
}

void GridModel::reset(bool reset_solver, bool reset_ac, bool reset_dc)
{
    if(reset_ac){
        id_me_to_ac_solver_ = std::vector<int>();
        id_ac_solver_to_me_ = std::vector<int>();
        slack_bus_id_ac_solver_ = Eigen::VectorXi();
        Ybus_ac_ = Eigen::SparseMatrix<cplx_type>();
    }

    if(reset_dc){
        id_me_to_dc_solver_ = std::vector<int>();
        id_dc_solver_to_me_ = std::vector<int>();
        slack_bus_id_dc_solver_ = Eigen::VectorXi();
        Ybus_dc_ = Eigen::SparseMatrix<cplx_type>();
    }

    Sbus_ = CplxVect();
    bus_pv_ = Eigen::VectorXi();
    bus_pq_ = Eigen::VectorXi();
    slack_weights_ = RealVect();
    need_reset_ = true;
    topo_changed_ = true;

    // reset the solvers
    if (reset_solver){
        _solver.reset();
        _dc_solver.reset();
    }
    // std::cout << "GridModel::reset called" << std::endl;
}

CplxVect GridModel::ac_pf(const CplxVect & Vinit,
                          int max_iter,
                          real_type tol)
{
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
    if(Vinit.size() != nb_bus){
        std::ostringstream exc_;
        exc_ << "GridModel::ac_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    bool conv = false;
    CplxVect res = CplxVect();

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    bool is_ac = true;
    bool reset_solver = topo_changed_;   // I reset the solver only if the topology change
    CplxVect V = pre_process_solver(Vinit, Ybus_ac_,
                                    id_me_to_ac_solver_, id_ac_solver_to_me_, slack_bus_id_ac_solver_,
                                    is_ac, reset_solver);

    // start the solver
    slack_weights_ = generators_.get_slack_weights(Ybus_ac_.rows(), id_me_to_ac_solver_); 
    conv = _solver.compute_pf(Ybus_ac_, V, Sbus_, slack_bus_id_ac_solver_, slack_weights_, bus_pv_, bus_pq_, max_iter, tol / sn_mva_);

    // store results (in ac mode)
    process_results(conv, res, Vinit, true, id_me_to_ac_solver_);

    // return the vector of complex voltage at each bus
    return res;
};

CplxVect GridModel::check_solution(const CplxVect & V_proposed, bool check_q_limits)
{
    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    const int nb_bus = static_cast<int>(V_proposed.size());
    bool is_ac = true;
    bool reset_solver = false;
    CplxVect V = pre_process_solver(V_proposed, Ybus_ac_, 
                                    id_me_to_ac_solver_, id_ac_solver_to_me_, slack_bus_id_ac_solver_,
                                    is_ac, reset_solver);

    // compute the mismatch
    CplxVect tmp = Ybus_ac_ * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - Sbus_.array();

    // store results
    CplxVect res = _get_results_back_to_orig_nodes(mis,
                                                   id_me_to_ac_solver_,
                                                   static_cast<int>(V_proposed.size())
                                                   );
    if(sn_mva_ != 1.) res *= sn_mva_;

    // now check reactive values for buses where there are generators and active values of slack bus
    // test for iterator though generator
    for(const auto & gen: generators_)
    {
        if(!gen.connected)
        {
            // the generator is disconnected, I do nothing
            continue;
        }
        if(check_q_limits)
        {
            // i need to check the reactive can be absorbed / produced by the generator
            real_type new_q = my_zero_;
            real_type react_this_bus = std::imag(res.coeff(gen.bus_id));
            if((react_this_bus >= gen.min_q_mvar) && (react_this_bus <= gen.max_q_mvar))
            {
                // this generator is able to handle all reactive
                new_q = my_zero_;
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
            res.coeffRef(gen.bus_id) = {std::real(res.coeff(gen.bus_id)), my_zero_};
        }

        // if(gen.id == gen_slackbus_)
        if(gen.is_slack)
        {
            // slack bus, by definition, can handle all active value
            // This is probably not the case with distributed slack !
            res.coeffRef(gen.bus_id) = {my_zero_, std::imag(res.coeff(gen.bus_id))};
        }
    }

    // set to 0 the error on the disconnected bus (it is not initialized at 0.0 in _get_results_back_to_orig_nodes)
    for(int bus_id = 0; bus_id < nb_bus; ++bus_id)
    {
        if(bus_status_[bus_id]) continue;
        res.coeffRef(bus_id) = my_zero_;
    }
    return res;
};

CplxVect GridModel::pre_process_solver(const CplxVect & Vinit, 
                                       Eigen::SparseMatrix<cplx_type> & Ybus,
                                       std::vector<int> & id_me_to_solver,
                                       std::vector<int> & id_solver_to_me,
                                       Eigen::VectorXi & slack_bus_id_solver,
                                       bool is_ac,
                                       bool reset_solver)
{
    // TODO get rid of the "is_ac" argument: this info is available in the _solver already
    // std::cout << "GridModel::pre_process_solver : topo_changed_ " << topo_changed_ << std::endl;
    // std::cout << "GridModel::pre_process_solver : reset_solver " << reset_solver << std::endl;

    bool reset_ac = topo_changed_ && is_ac;
    bool reset_dc = topo_changed_ && !is_ac;
    // if(need_reset_){ // TODO optimization when it's not mandatory to start from scratch
    if(topo_changed_) reset(reset_solver, reset_ac, reset_dc);  // TODO what if pv and pq changed ? :O
    else{
        // topo is not changed, but i can still reset the solver (TODO: no necessarily needed !)
        if (reset_solver)
        {
            if(is_ac) _solver.reset();
            else _dc_solver.reset();
        }
    }
    slack_bus_id_ = generators_.get_slack_bus_id();
    if(topo_changed_){
        // TODO do not reinit Ybus if the topology does not change
        init_Ybus(Ybus, id_me_to_solver, id_solver_to_me);
        fillYbus(Ybus, is_ac, id_me_to_solver);
    }
    init_Sbus(Sbus_, id_me_to_solver, id_solver_to_me, slack_bus_id_solver);
    fillpv_pq(id_me_to_solver, id_solver_to_me, slack_bus_id_solver); // TODO what if pv and pq changed ? :O
    
    generators_.init_q_vector(static_cast<int>(bus_vn_kv_.size()));
    fillSbus_me(Sbus_, is_ac, id_me_to_solver, slack_bus_id_solver);

    const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
    CplxVect V = CplxVect::Constant(nb_bus_solver, init_vm_pu_);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        int bus_me_id = id_solver_to_me[bus_solver_id];  //TODO DEBUG MODE : POSSIBLE SEGFAULT
        cplx_type tmp = Vinit(bus_me_id);
        V(bus_solver_id) = tmp;
    }
    generators_.set_vm(V, id_me_to_solver);
    return V;
}

CplxVect GridModel::_get_results_back_to_orig_nodes(const CplxVect & res_tmp, 
                                                    std::vector<int> & id_me_to_solver,
                                                    int size)
{
    CplxVect res = CplxVect::Constant(size, {init_vm_pu_, my_zero_});
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
    for (int bus_id_me=0; bus_id_me < nb_bus; ++bus_id_me){
        if(!bus_status_[bus_id_me]) continue;  // nothing is done if the bus is connected
        int bus_id_solver = id_me_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "GridModel::_get_results_back_to_orig_nodes: the bus with id ";
            exc_ << bus_id_me;
            exc_ << " is connected to a disconnected bus (solver side)";
            throw std::runtime_error(exc_.str());
        }
        res(bus_id_me) = res_tmp(bus_id_solver);
    }
    return res;
}

void GridModel::process_results(bool conv, CplxVect & res, const CplxVect & Vinit, bool ac,
                                std::vector<int> & id_me_to_solver)
{
    if (conv){
        if(compute_results_){
            // compute the results of the flows, P,Q,V of loads etc.
            compute_results(ac);
        }
        need_reset_ = false;
        const CplxVect & res_tmp = ac ? _solver.get_V(): _dc_solver.get_V() ;

        // convert back the results to "big" vector
        res = _get_results_back_to_orig_nodes(res_tmp,
                                              id_me_to_solver,
                                              static_cast<int>(Vinit.size()));
    } else {
        //powerflow diverge
        reset_results();
        need_reset_ = true;  // in this case, the powerflow diverge, so i need to recompute Ybus next time
    }
}

void GridModel::init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                          std::vector<int>& id_me_to_solver,
                          std::vector<int>& id_solver_to_me){
    //TODO get disconnected bus !!! (and have some conversion for it)
    //1. init the conversion bus
    const int nb_bus_init = static_cast<int>(bus_vn_kv_.size());
    id_me_to_solver = std::vector<int>(nb_bus_init, _deactivated_bus_id);  // by default, if a bus is disconnected, then it has a -1 there
    id_solver_to_me = std::vector<int>();
    id_solver_to_me.reserve(nb_bus_init);
    int bus_id_solver=0;
    for(int bus_id_me=0; bus_id_me < nb_bus_init; ++bus_id_me){
        if(bus_status_[bus_id_me]){
            // bus is connected
            id_solver_to_me.push_back(bus_id_me);
            id_me_to_solver[bus_id_me] = bus_id_solver;
            ++bus_id_solver;
        }
    }
    const int nb_bus = static_cast<int>(id_solver_to_me.size());

    Ybus = Eigen::SparseMatrix<cplx_type>(nb_bus, nb_bus);
    Ybus.reserve(nb_bus + 2*powerlines_.nb() + 2*trafos_.nb());


}

void GridModel::init_Sbus(CplxVect & Sbus,
                          std::vector<int>& id_me_to_solver,
                          std::vector<int>& id_solver_to_me,
                          Eigen::VectorXi & slack_bus_id_solver){
    
    const int nb_bus = static_cast<int>(id_solver_to_me.size());                 
    Sbus = CplxVect::Constant(nb_bus, 0.);
    slack_bus_id_solver = Eigen::VectorXi::Zero(slack_bus_id_.size());

    size_t i = 0;
    for(auto el: slack_bus_id_) {
        slack_bus_id_solver(i) = id_me_to_solver[el];
        ++i;
    }
    
    if(is_in_vect(_deactivated_bus_id, slack_bus_id_solver)){
        // TODO improve error message with the gen_id
        // TODO DEBUG MODE: only check that in debug mode
        throw std::runtime_error("One of the slack bus is disconnected !");
    }
}
void GridModel::fillYbus(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int>& id_me_to_solver){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/

    // init the Ybus matrix
    std::vector<Eigen::Triplet<cplx_type> > tripletList;
    tripletList.reserve(bus_vn_kv_.size() + 4*powerlines_.nb() + 4*trafos_.nb() + shunts_.nb());
    powerlines_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    shunts_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    trafos_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    loads_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    sgens_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    storages_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    generators_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    res.makeCompressed();
}

void GridModel::fillSbus_me(CplxVect & Sbus, bool ac, const std::vector<int>& id_me_to_solver, Eigen::VectorXi & slack_bus_id_solver)
{
    // init the Sbus vector
    powerlines_.fillSbus(Sbus, true, id_me_to_solver);
    trafos_.fillSbus(Sbus, ac, id_me_to_solver);
    shunts_.fillSbus(Sbus, true, id_me_to_solver);
    loads_.fillSbus(Sbus, true, id_me_to_solver);
    sgens_.fillSbus(Sbus, true, id_me_to_solver);
    storages_.fillSbus(Sbus, true, id_me_to_solver);
    generators_.fillSbus(Sbus, true, id_me_to_solver);

    if (sn_mva_ != 1.0) Sbus /= sn_mva_;
    // in dc mode, this is used for the phase shifter, this should not be divided by sn_mva_ !
    trafos_.hack_Sbus_for_dc_phase_shifter(Sbus, ac, id_me_to_solver);
}

void GridModel::fillpv_pq(const std::vector<int>& id_me_to_solver,
                          std::vector<int>& id_solver_to_me,
                          Eigen::VectorXi & slack_bus_id_solver)
{
    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    const int nb_bus = static_cast<int>(id_solver_to_me.size());  // number of bus in the solver!
    std::vector<int> bus_pq;
    std::vector<int> bus_pv;
    std::vector<bool> has_bus_been_added(nb_bus, false);

    // std::cout << "id_me_to_solver.size(): " << id_me_to_solver.size() << std::endl;

    bus_pv_ = Eigen::VectorXi();
    bus_pq_ = Eigen::VectorXi();
    powerlines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    shunts_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    trafos_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    loads_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    storages_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    sgens_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    generators_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);

    for(int bus_id = 0; bus_id< nb_bus; ++bus_id){
        if(is_in_vect(bus_id, slack_bus_id_solver)) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pv.data(), bus_pv.size());
    bus_pq_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pq.data(), bus_pq.size());
}
void GridModel::compute_results(bool ac){
    // retrieve results from powerflow
    const auto & Va = ac ? _solver.get_Va() : _dc_solver.get_Va();
    const auto & Vm = ac ? _solver.get_Vm() : _dc_solver.get_Vm();
    const auto & V = ac ? _solver.get_V() : _dc_solver.get_V();

    const std::vector<int> & id_me_to_solver = ac ? id_me_to_ac_solver_ : id_me_to_dc_solver_;
    // for powerlines
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for trafo
    trafos_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for loads
    loads_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for static gen
    sgens_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for storage units
    storages_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);
    // for prods
    generators_.compute_results(Va, Vm, V, id_me_to_solver, bus_vn_kv_, sn_mva_, ac);

    //handle_slack_bus active power
    CplxVect mismatch;  // power mismatch at each bus (SOLVER BUS !!!)
    RealVect ractive_mismatch;  // not used in dc mode (DO NOT ATTEMPT TO USE IT THERE)
    RealVect active_mismatch;
    if(ac){
        // In AC mode i am not forced to run through all the grid
        auto tmp = (Ybus_ac_ * V).conjugate();
        mismatch = V.array() * tmp.array() - Sbus_.array();
        active_mismatch = mismatch.real() * sn_mva_;
    } else{
        active_mismatch = RealVect::Zero(V.size());
        //TODO SLACK: improve distributed slack for DC mode !
        // it is possible to know in advance the contribution of each slack generators (sum(Sbus) MW 
        // to split among the contributing generators) so it's possible to "mess with" Sbus 
        // for such purpose
        const auto id_slack = slack_bus_id_dc_solver_(0);
        active_mismatch(id_slack) = -Sbus_.real().sum() * sn_mva_;
    }
    generators_.set_p_slack(active_mismatch, id_me_to_solver);

    if(ac) ractive_mismatch = mismatch.imag() * sn_mva_;
    // mainly to initialize the Q value of the generators in dc (just fill it with 0.)
    generators_.set_q(ractive_mismatch, id_me_to_solver, ac);
}

void GridModel::reset_results(){
    powerlines_.reset_results();
    shunts_.reset_results();
    trafos_.reset_results();
    loads_.reset_results();
    sgens_.reset_results();
    storages_.reset_results();
    generators_.reset_results();
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
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
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

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    bool is_ac = false;
    bool reset_solver = topo_changed_;  // I reset the solver only if the topology change
    CplxVect V = pre_process_solver(Vinit, Ybus_dc_,
                                    id_me_to_dc_solver_, id_dc_solver_to_me_, slack_bus_id_dc_solver_,
                                    is_ac, reset_solver);

    // start the solver
    slack_weights_ = generators_.get_slack_weights(Ybus_dc_.rows(), id_me_to_dc_solver_);
    conv = _dc_solver.compute_pf(Ybus_dc_, V, Sbus_, slack_bus_id_dc_solver_, slack_weights_, bus_pv_, bus_pq_, max_iter, tol);

    // store results (fase -> because I am in dc mode)
    process_results(conv, res, Vinit, false, id_me_to_dc_solver_);
    return res;
}

/**
Retrieve the number of connected buses
**/
int GridModel::nb_bus() const
{
    int res = 0;
    for(const auto & el : bus_status_)
    {
        if(el) ++res;
    }
    return res;
}

void GridModel::add_gen_slackbus(int gen_id, real_type weight){
    if(gen_id < 0)
    {
        std::ostringstream exc_;
        exc_ << "GridModel::add_gen_slackbus: Slack bus should be an id of a generator, thus positive. You provided: ";
        exc_ << gen_id;
        throw std::runtime_error(exc_.str());
    }
    if(gen_id > generators_.nb())
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
    generators_.add_slackbus(gen_id, weight);
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
    if(gen_id > generators_.nb())
    {
        // TODO DEBUG MODE: only check when in debug mode
        std::ostringstream exc_;
        exc_ << "GridModel::remove_gen_slackbus: There are only " << generators_.nb() << " generators on the grid. ";
        exc_ << "Generator with id " << gen_id << " does not exist and can't be the slack bus";
        throw std::runtime_error(exc_.str());
    }
    generators_.remove_slackbus(gen_id);
}

/** GRID2OP SPECIFIC REPRESENTATION **/
void GridModel::update_bus_status(int nb_bus_before,
                                  Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> > active_bus)
{
    for(int bus_id = 0; bus_id < active_bus.rows(); ++bus_id)
    {
        if(active_bus(bus_id, 0)){
            reactivate_bus(bus_id);
        }else{
            deactivate_bus(bus_id);
        }
        if(active_bus(bus_id, 1)){
            reactivate_bus(bus_id + nb_bus_before);
        }else{
            deactivate_bus(bus_id + nb_bus_before);
        }
    }
}

void GridModel::update_gens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &GridModel::change_p_gen);
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

void GridModel::update_topo(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<Eigen::Array<int,  Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    const int nb_bus = static_cast<int>(bus_status_.size());
    for(int i = 0; i < nb_bus; ++i) bus_status_[i] = false;

    update_topo_generic(has_changed, new_values,
                        load_pos_topo_vect_, load_to_subid_,
                        &GridModel::reactivate_load,
                        &GridModel::change_bus_load,
                        &GridModel::deactivate_load
                        );
    update_topo_generic(has_changed, new_values,
                        gen_pos_topo_vect_, gen_to_subid_,
                        &GridModel::reactivate_gen,
                        &GridModel::change_bus_gen,
                        &GridModel::deactivate_gen
                        );
    update_topo_generic(has_changed, new_values,
                        storage_pos_topo_vect_, storage_to_subid_,
                        &GridModel::reactivate_storage,
                        &GridModel::change_bus_storage,
                        &GridModel::deactivate_storage
                        );

    // NB we suppose that if a powerline (or a trafo) is disconnected, then both its ends are
    // and same for trafo, obviously
    update_topo_generic(has_changed, new_values,
                        line_or_pos_topo_vect_, line_or_to_subid_,
                        &GridModel::reactivate_powerline,
                        &GridModel::change_bus_powerline_or,
                        &GridModel::deactivate_powerline
                        );
    update_topo_generic(has_changed, new_values,
                        line_ex_pos_topo_vect_, line_ex_to_subid_,
                        &GridModel::reactivate_powerline,
                        &GridModel::change_bus_powerline_ex,
                        &GridModel::deactivate_powerline
                        );
    update_topo_generic(has_changed, new_values,
                        trafo_hv_pos_topo_vect_, trafo_hv_to_subid_,
                        &GridModel::reactivate_trafo,
                        &GridModel::change_bus_trafo_hv,
                        &GridModel::deactivate_trafo
                        );
    update_topo_generic(has_changed, new_values,
                        trafo_lv_pos_topo_vect_, trafo_lv_to_subid_,
                        &GridModel::reactivate_trafo,
                        &GridModel::change_bus_trafo_lv,
                        &GridModel::deactivate_trafo
                        );
}
