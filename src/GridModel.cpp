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
    reset(true);

    init_vm_pu_ = other.init_vm_pu_;
    sn_mva_ = other.sn_mva_;

    // assign the right solver
    _solver.change_solver(other._solver.get_type());
    compute_results_ = other.compute_results_;

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

    // 6. static generators
    sgens_ = other.sgens_;

    // 6. storage units
    storages_ = other.storages_;

    // 7. slack bus
    gen_slackbus_ = other.gen_slackbus_;
    slack_bus_id_ = other.slack_bus_id_;

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
                            res_storage,
                            gen_slackbus_
                            );
    return res;
};

void GridModel::set_state(GridModel::StateRes & my_state)
{
    // after loading back, the instance need to be reset anyway
    // TODO see if it's worth the trouble NOT to do it
    reset(true);
    need_reset_ = true;
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
    int gen_slackbus = std::get<14>(my_state);

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

    // other stuff
    gen_slackbus_ = gen_slackbus;

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

void GridModel::reset(bool reset_solver)
{
    Ybus_ = Eigen::SparseMatrix<cplx_type>();
    Sbus_ = CplxVect();
    id_me_to_solver_ = std::vector<int>();
    id_solver_to_me_ = std::vector<int>();
    slack_bus_id_solver_ = -1;
    bus_pv_ = Eigen::VectorXi();
    bus_pq_ = Eigen::VectorXi();
    need_reset_ = true;

    // reset the solvers
    if (reset_solver) _solver.reset();
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
    CplxVect V = pre_process_solver(Vinit, true, true);

    // start the solver
    conv = _solver.compute_pf(Ybus_, V, Sbus_, bus_pv_, bus_pq_, max_iter, tol / sn_mva_);

    // store results
    process_results(conv, res, Vinit);

    // return the vector of complex voltage at each bus
    return res;
};

CplxVect GridModel::check_solution(const CplxVect & V_proposed, bool check_q_limits)
{
    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    const int nb_bus = static_cast<int>(V_proposed.size());
    CplxVect V = pre_process_solver(V_proposed, true, false);

    // compute the mismatch
    CplxVect tmp = Ybus_ * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - Sbus_.array();

    // store results
    CplxVect res = _get_results_back_to_orig_nodes(mis, static_cast<int>(V_proposed.size()));
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

        if(gen.id == gen_slackbus_)
        {
            // slack bus, by definition, can handle all active value
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


CplxVect GridModel::pre_process_solver(const CplxVect & Vinit, bool is_ac, bool reset_solver)
{
    // TODO get rid of the "is_ac" argument: this info is available in the _solver already

    // if(need_reset_){ // TODO optimization when it's not mandatory to start from scratch
    reset(reset_solver);
    slack_bus_id_ = generators_.get_slack_bus_id(gen_slackbus_);
    init_Ybus(Ybus_, Sbus_, id_me_to_solver_, id_solver_to_me_, slack_bus_id_solver_);
    fillYbus(Ybus_, is_ac, id_me_to_solver_);
    fillpv_pq(id_me_to_solver_);
    generators_.init_q_vector(static_cast<int>(bus_vn_kv_.size()));
    // }
    fillSbus_me(Sbus_, is_ac, id_me_to_solver_, slack_bus_id_solver_);

    const int nb_bus_solver = static_cast<int>(id_solver_to_me_.size());
    CplxVect V = CplxVect::Constant(nb_bus_solver, init_vm_pu_);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        int bus_me_id = id_solver_to_me_[bus_solver_id];  //POSSIBLE SEGFAULT
        cplx_type tmp = Vinit(bus_me_id);
        V(bus_solver_id) = tmp;
        // TODO save this V somewhere
    }
    generators_.set_vm(V, id_me_to_solver_);
    return V;
}

CplxVect GridModel::_get_results_back_to_orig_nodes(const CplxVect & res_tmp, int size)
{
    CplxVect res = CplxVect::Constant(size, {init_vm_pu_, my_zero_});
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
    for (int bus_id_me=0; bus_id_me < nb_bus; ++bus_id_me){
        if(!bus_status_[bus_id_me]) continue;  // nothing is done if the bus is connected
        int bus_id_solver = id_me_to_solver_[bus_id_me];
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

void GridModel::process_results(bool conv, CplxVect & res, const CplxVect & Vinit)
{
    if (conv){
        if(compute_results_){
            // compute the results of the flows, P,Q,V of loads etc.
            compute_results();
        }
        need_reset_ = false;
        const CplxVect & res_tmp = _solver.get_V();
        // convert back the results to "big" vector
        res = _get_results_back_to_orig_nodes(res_tmp, static_cast<int>(Vinit.size()));
    } else {
        //powerflow diverge
        reset_results();
        need_reset_ = true;  // in this case, the powerflow diverge, so i need to recompute Ybus next time
    }
}

void GridModel::init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                          CplxVect & Sbus,
                          std::vector<int>& id_me_to_solver,
                          std::vector<int>& id_solver_to_me,
                          int & slack_bus_id_solver){
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

    Sbus = CplxVect::Constant(nb_bus, 0.);
    slack_bus_id_solver = id_me_to_solver[slack_bus_id_];
    if(slack_bus_id_solver == _deactivated_bus_id){
        //TODO improve error message with the gen_id
        throw std::runtime_error("The slack bus is disconnected.");
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

void GridModel::fillSbus_me(CplxVect & Sbus, bool ac, const std::vector<int>& id_me_to_solver, int slack_bus_id_solver)
{
    // init the Sbus vector
    powerlines_.fillSbus(Sbus, true, id_me_to_solver);
    shunts_.fillSbus(Sbus, true, id_me_to_solver);
    trafos_.fillSbus(Sbus, ac, id_me_to_solver);
    loads_.fillSbus(Sbus, true, id_me_to_solver);
    sgens_.fillSbus(Sbus, true, id_me_to_solver);
    storages_.fillSbus(Sbus, true, id_me_to_solver);
    generators_.fillSbus(Sbus, true, id_me_to_solver);

    // handle slack bus (in ac only)
    if(ac)
    {
        real_type sum_active = Sbus.sum().real();
        Sbus.coeffRef(slack_bus_id_solver) -= sum_active;
    }
    if (sn_mva_ != 1.0) Sbus /= sn_mva_;
}

void GridModel::fillpv_pq(const std::vector<int>& id_me_to_solver)
{
    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    const int nb_bus = static_cast<int>(id_solver_to_me_.size());  // number of bus in the solver!
    std::vector<int> bus_pq;
    std::vector<int> bus_pv;
    std::vector<bool> has_bus_been_added(nb_bus, false);

    bus_pv_ = Eigen::VectorXi();
    bus_pq_ = Eigen::VectorXi();
    powerlines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    shunts_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    trafos_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    loads_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    storages_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    sgens_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    generators_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);

    for(int bus_id = 0; bus_id< nb_bus; ++bus_id){
        if(bus_id == slack_bus_id_solver_) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pv.data(), bus_pv.size());
    bus_pq_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pq.data(), bus_pq.size());
}
void GridModel::compute_results(){
    // TODO "deactivate" the Q value for DC

    // retrieve results from powerflow
    const auto & Va = _solver.get_Va();
    const auto & Vm = _solver.get_Vm();
    const auto & V = _solver.get_V();

    // for powerlines
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for trafo
    trafos_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for loads
    loads_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for static gen
    sgens_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for storage units
    storages_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);
    // for prods
    generators_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_, sn_mva_);

    //handle_slack_bus
    real_type p_slack = powerlines_.get_p_slack(slack_bus_id_);
    p_slack += trafos_.get_p_slack(slack_bus_id_);
    p_slack += loads_.get_p_slack(slack_bus_id_);
    p_slack += sgens_.get_p_slack(slack_bus_id_);
    p_slack += storages_.get_p_slack(slack_bus_id_);
    p_slack += shunts_.get_p_slack(slack_bus_id_);
    generators_.set_p_slack(gen_slackbus_, p_slack);

    // handle gen_q now
    std::vector<real_type> q_by_bus = std::vector<real_type>(bus_vn_kv_.size(), 0.);
    powerlines_.get_q(q_by_bus);
    trafos_.get_q(q_by_bus);
    loads_.get_q(q_by_bus);
    storages_.get_q(q_by_bus);
    sgens_.get_q(q_by_bus);
    shunts_.get_q(q_by_bus);

    generators_.set_q(q_by_bus);
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

CplxVect GridModel::dc_pf_old(const CplxVect & Vinit,
                              int max_iter,  // not used for DC
                              real_type tol  // not used for DC
                              )
{
    // TODO refactor that with ac pf, this is mostly done, but only mostly...
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
    if(static_cast<int>(Vinit.size()) != nb_bus){
        std::ostringstream exc_;
        exc_ << "GridModel::dc_pf_old: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    Eigen::SparseMatrix<cplx_type> dcYbus_tmp;
    CplxVect Sbus_tmp;
    std::vector<int> id_me_to_solver;
    std::vector<int> id_solver_to_me;
    int slack_bus_id_solver;

    //if(need_reset_){
    slack_bus_id_ = generators_.get_slack_bus_id(gen_slackbus_);
    init_Ybus(dcYbus_tmp, Sbus_tmp, id_me_to_solver, id_solver_to_me, slack_bus_id_solver);
    fillYbus(dcYbus_tmp, false, id_me_to_solver);
    // fillpv_pq(id_me_to_solver);
    //}
    fillSbus_me(Sbus_tmp, false, id_me_to_solver, slack_bus_id_solver);


    // extract only connected bus from Vinit
    const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());

    // DC SOLVER STARTS HERE
    // TODO all this should rather be one in a "dc solver" instead of here
    // remove the slack bus

    // remove the slack bus from Ybus
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    dcYbus_tmp.makeCompressed();
    Eigen::SparseMatrix<real_type> dcYbus = Eigen::SparseMatrix<real_type>(nb_bus_solver - 1, nb_bus_solver - 1);
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(dcYbus_tmp.nonZeros());
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(dcYbus_tmp, k); it; ++it)
        {
            int row_res = static_cast<int>(it.row());
            if(row_res == slack_bus_id_solver) continue;
            row_res = row_res > slack_bus_id_solver ? row_res - 1 : row_res;
            int col_res = static_cast<int>(it.col());
            col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
            tripletList.push_back(Eigen::Triplet<real_type> (row_res, col_res, std::real(it.value())));
        }
    }
    dcYbus.setFromTriplets(tripletList.begin(), tripletList.end());
    dcYbus.makeCompressed();

    // initialize the solver
    Eigen::SparseLU<Eigen::SparseMatrix<real_type>, Eigen::COLAMDOrdering<int> >   solver;
    solver.analyzePattern(dcYbus);
    solver.factorize(dcYbus);
    if(solver.info() != Eigen::Success) {
        // matrix is not connected
        return CplxVect();
    }

    // remove the slack bus from Sbus
    RealVect Sbus = RealVect::Constant(nb_bus_solver - 1, 0.);
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        int col_res = k;
        col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
        Sbus(col_res) = std::real(Sbus_tmp(k));
    }

    // solve for theta: Sbus = dcY . theta
    RealVect Va_dc = solver.solve(Sbus);
    if(solver.info() != Eigen::Success) {
        // solving failed, this should not happen in dc ...
        return CplxVect();
    }

    // retrieve back the results in the proper shape
    const int nb_bus_me = static_cast<int>(bus_vn_kv_.size());
    int bus_id_solver;
    RealVect Va = RealVect::Constant(nb_bus_me, 0.);
    // fill Va from dc approx
    for (int bus_id_me=0; bus_id_me < nb_bus_me; ++bus_id_me){
        if(bus_id_me == slack_bus_id_) continue;  // slack bus is handled elsewhere
        if(!bus_status_[bus_id_me]) continue;  // nothing is done if the bus is not connected

        bus_id_solver = id_me_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One bus is both connected and disconnected");
        }
        bus_id_solver = bus_id_solver > slack_bus_id_solver ? bus_id_solver - 1 : bus_id_solver;
        Va(bus_id_me) = Va_dc(bus_id_solver);
    }
    Va.array() +=  std::arg(Vinit(slack_bus_id_));

    // fill Vm either Vinit if pq or Vm if pv (TODO)
    RealVect Vm;
    if(false){
        Vm = Vinit.array().abs();  // fill Vm = Vinit for all
        // put Vm = 0. for disconnected bus
        for (int bus_id_me=0; bus_id_me < nb_bus_me; ++bus_id_me){
            if(bus_status_[bus_id_me]) continue;  // nothing is done if the bus is connected
            Vm(bus_id_me) = 0.;
        }
        // put Vm = Vm of turned on gen
        // generators_.get_vm_for_dc(Vm);
        // assign vm of the slack bus
        Vm(slack_bus_id_) =  std::abs(Vinit(slack_bus_id_));
    }
    else{
        Vm = RealVect::Constant(Vinit.size(), 1.0);
        for (int bus_id_me=0; bus_id_me < nb_bus_me; ++bus_id_me){
            if(bus_status_[bus_id_me]) continue;  // nothing is done if the bus is connected
            Vm(bus_id_me) = 0.;
        }
        generators_.get_vm_for_dc(Vm);
    }
    //END of the SOLVER PART

    //TODO handle Vm = Vm (gen) for connected generators
    return Vm.array() * (Va.array().cos().cast<cplx_type>() + my_i * Va.array().sin().cast<cplx_type>());
}

CplxVect GridModel::dc_pf(const CplxVect & Vinit,
                          int max_iter,  // not used for DC
                          real_type tol  // not used for DC
                          )
{
    const int nb_bus = static_cast<int>(bus_vn_kv_.size());
    if(Vinit.size() != nb_bus){
        std::ostringstream exc_;
        exc_ << "GridModel::dc_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    SolverType solver_type = _solver.get_type();
    _solver.change_solver(SolverType::DC);

    bool conv = false;
    CplxVect res = CplxVect();

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    CplxVect V = pre_process_solver(Vinit, false, true);

    // start the solver
    conv = _solver.compute_pf(Ybus_, V, Sbus_, bus_pv_, bus_pq_, max_iter, tol);

    // store results
    process_results(conv, res, Vinit);

    // put back the solver to its original state
    // TODO add a better handling of this!
    _solver.change_solver(solver_type);

    // return the vector of complex voltage at each bus
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

void GridModel::add_gen_slackbus(int gen_id){
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
    gen_slackbus_ = gen_id;
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
