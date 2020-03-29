#include "GridModel.h"

// const int GridModel::_deactivated_bus_id = -1;

void GridModel::init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo){
    /**
    initialize the bus_vn_kv_ member
    and
    initialize the Ybus_ matrix at the proper shape
    **/
    int nb_bus = bus_vn_kv.size();
    bus_vn_kv_ = bus_vn_kv;  // base_kv

    bus_status_ = std::vector<bool>(nb_bus, true); // by default everything is connected
}

void GridModel::init_powerlines(const Eigen::VectorXd & branch_r,
                                const Eigen::VectorXd & branch_x,
                                const Eigen::VectorXcd & branch_h,
                                const Eigen::VectorXi & branch_from_id,
                                const Eigen::VectorXi & branch_to_id
                                )
{
    powerlines_.init(branch_r, branch_x, branch_h, branch_from_id, branch_to_id);
}

void GridModel::init_shunt(const Eigen::VectorXd & shunt_p_mw,
                           const Eigen::VectorXd & shunt_q_mvar,
                           const Eigen::VectorXi & shunt_bus_id)
{
    shunts_.init(shunt_p_mw, shunt_q_mvar, shunt_bus_id);
}

void GridModel::init_trafo(const Eigen::VectorXd & trafo_r,
                           const Eigen::VectorXd & trafo_x,
                           const Eigen::VectorXcd & trafo_b,
                           const Eigen::VectorXd & trafo_tap_step_pct,
            //                        const Eigen::VectorXd & trafo_tap_step_degree,
                           const Eigen::VectorXd & trafo_tap_pos,
                           const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                           const Eigen::VectorXi & trafo_hv_id,
                           const Eigen::VectorXi & trafo_lv_id
                           )
{
    /**
    INPUT DATA ARE ALREADY PAIR UNIT !!
    DOES NOT WORK WITH POWERLINES
    **/
    /**
    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!
    Eigen::VectorXd ratio = 1.0 + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array() * (2*trafo_tap_hv.array().cast<double>() - 1.0);

    transformers_r_ = trafo_r;
    transformers_x_ = trafo_x;
    transformers_h_ = trafo_b;
    transformers_ratio_ = ratio;
    transformers_bus_hv_id_ = trafo_hv_id;
    transformers_bus_lv_id_ = trafo_lv_id;
    transformers_status_ = std::vector<bool>(trafo_r.size(), true);
    **/
    trafos_.init(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_tap_hv, trafo_hv_id, trafo_lv_id);
}


void GridModel::init_generators(const Eigen::VectorXd & generators_p,
                     const Eigen::VectorXd & generators_v,
                     const Eigen::VectorXi & generators_bus_id)
{
    generators_p_ = generators_p;
    generators_v_ = generators_v;
    generators_bus_id_ = generators_bus_id;
    generators_status_ = std::vector<bool>(generators_p.size(), true);
}

void GridModel::init_loads(const Eigen::VectorXd & loads_p,
                const Eigen::VectorXd & loads_q,
                const Eigen::VectorXi & loads_bus_id)
{
    loads_p_ = loads_p;
    loads_q_ = loads_q;
    loads_bus_id_ = loads_bus_id;
    loads_status_ = std::vector<bool>(loads_p.size(), true);
}

bool GridModel::compute_newton(const Eigen::VectorXcd & Vinit,
                               int max_iter,
                               double tol)
{
    // TODO optimization when it's not mandatory to start from scratch
    bool res = false;
    if(need_reset_){
        fillYbus();
        _solver.reset();
    }

    // extract only connected bus from Vinit
    int nb_bus_solver = id_solver_to_me_.size();
    Eigen::VectorXcd V = Eigen::VectorXcd::Constant(id_solver_to_me_.size(), 1.04);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        int bus_me_id = id_solver_to_me_[bus_solver_id];  //POSSIBLE SEGFAULT
        cdouble tmp = Vinit(bus_me_id);
        V(bus_solver_id) = tmp;
        // TODO save this V somewhere
    }

    res = _solver.do_newton(Ybus_, V, Sbus_, bus_pv_, bus_pq_, max_iter, tol);
    if (res){
        compute_results();
        need_reset_ = false;
    } else {
        //powerflow diverge
        reset_results();
        need_reset_ = true;  // in this case, the powerflow diverge, so i need to recompute Ybus next time
    }
    return res;
};

void GridModel::init_Ybus(){
    // int nb_line = powerlines_r_.size();
    // int nb_trafo = transformers_r_.size();

    //TODO get disconnected bus !!! (and have some conversion for it)
    //1. init the conversion bus
    int nb_bus_init = bus_vn_kv_.size();
    id_me_to_solver_ = std::vector<int>(nb_bus_init, _deactivated_bus_id);  // by default, if a bus is disconnected, then it has a -1 there
    id_solver_to_me_ = std::vector<int>();
    id_solver_to_me_.reserve(bus_vn_kv_.size());
    int bus_id_solver=0;
    for(int bus_id_me=0; bus_id_me < nb_bus_init; ++bus_id_me){
        if(bus_status_[bus_id_me]){
            // bus is connected
            id_solver_to_me_.push_back(bus_id_me);
            id_me_to_solver_[bus_id_me] = bus_id_solver;
            ++bus_id_solver;
        }
    }
    int nb_bus = id_me_to_solver_.size();

    Ybus_ = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);
    Ybus_.reserve(nb_bus + 2*powerlines_.nb() + 2*trafos_.nb());  //TODO optimize with number of connected powerlines or trafo

    // init diagonal coefficients
    for(int bus_id=0; bus_id < nb_bus; ++bus_id){
        // assign diagonal coefficient
        Ybus_.insert(bus_id, bus_id) = 0.;
    }

    Sbus_ = Eigen::VectorXcd::Constant(nb_bus, 0.);
    slack_bus_id_solver_ = id_me_to_solver_[slack_bus_id_];
    if(slack_bus_id_solver_ == _deactivated_bus_id){
        //TODO improve error message with the gen_id
        throw std::runtime_error("The slack bus is disconnected.");
    }
}

void GridModel::fillYbus(){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/
    // TODO split that between fillYbus, computeSbus and initPvPq
    init_Ybus();

    // init the Ybus matrix
    // fillYbusBranch(Ybus_, true);
    powerlines_.fillYbus(Ybus_, true, id_me_to_solver_);
    // fillYbusShunt(Ybus_, true);
    shunts_.fillYbus(Ybus_, true, id_me_to_solver_);
    // fillYbusTrafo(Ybus_, true);
    trafos_.fillYbus(Ybus_, true, id_me_to_solver_);

    // init the Sbus vector
    cdouble my_i = 1.0i;
    int bus_id_me, bus_id_solver;
    int nb_bus = id_solver_to_me_.size();
    int nb_load = loads_p_.size();
    int nb_gen = generators_p_.size();
    double sum_active = 0.;
    double tmp;

    std::vector<bool> has_gen_conn(nb_bus, false);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!generators_status_[gen_id]) continue;

        bus_id_me = generators_bus_id_(gen_id);
        bus_id_solver = id_me_to_solver_[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        tmp = generators_p_(gen_id);
        Sbus_.coeffRef(bus_id_solver) += tmp; // + my_i * generators_p_(gen_id);
        sum_active += tmp;
        has_gen_conn[bus_id_solver] = true;
    }
    for(int load_id = 0; load_id < nb_load; ++load_id){
        //  i don't do anything if the load is disconnected
        if(!loads_status_[load_id]) continue;

        bus_id_me = loads_bus_id_(load_id);
        bus_id_solver = id_me_to_solver_[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One load is connected to a disconnected bus.");
        }
        tmp = loads_p_(load_id);
        Sbus_.coeffRef(bus_id_solver) -= tmp + my_i * loads_q_(load_id);
        sum_active -= tmp;
    }
    Sbus_.coeffRef(slack_bus_id_solver_) -= sum_active;
    //TODO put the shunt here (but test before it can be done)

    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    std::vector<int> bus_pq;
    std::vector<int> bus_pv;
    std::vector<bool> has_bus_been_added(nb_bus, false);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!generators_status_[gen_id]) continue;

        bus_id_me = generators_bus_id_(gen_id);
        bus_id_solver = id_me_to_solver_[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        if(bus_id_solver == slack_bus_id_solver_) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver);
        has_bus_been_added[bus_id_solver] = true;  // don't add it a second time
    }
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
     //TODO check it has converged!

    // retrieve results from powerflow
    const auto & Va = _solver.get_Va();
    const auto & Vm = _solver.get_Vm();
    const auto & V = _solver.get_V();

    // for powerlines
    /**
    int nb_line = powerlines_r_.size();
    Eigen::VectorXd ratio = Eigen::VectorXd::Constant(nb_line, 1.0);
    res_powerlines(Va, Vm, V, powerlines_status_,
                   nb_line,
                   powerlines_r_,
                   powerlines_x_,
                   powerlines_h_,
                   ratio,
                   powerlines_bus_or_id_,
                   powerlines_bus_ex_id_,
                   res_powerline_por_,  // in MW
                   res_powerline_qor_,  // in MVar
                   res_powerline_vor_,  // in kV
                   res_powerline_aor_,  // in kA
                   res_powerline_pex_,  // in MW
                   res_powerline_qex_,  // in MVar
                   res_powerline_vex_,  // in kV
                   res_powerline_aex_  // in kA
                   );
    **/
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);

    // for trafo
    /**
    int nb_trafo = transformers_r_.size();
    res_powerlines(Va, Vm, V, transformers_status_,
                   nb_trafo,
                   transformers_r_,
                   transformers_x_,
                   transformers_h_,
                   transformers_ratio_,
                   transformers_bus_hv_id_,
                   transformers_bus_lv_id_,
                   res_trafo_por_,  // in MW
                   res_trafo_qor_,  // in MVar
                   res_trafo_vor_,  // in kV
                   res_trafo_aor_,  // in kA
                   res_trafo_pex_,  // in MW
                   res_trafo_qex_,  // in MVar
                   res_trafo_vex_,  // in kV
                   res_trafo_aex_  // in kA
                   );
                   **/

    trafos_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);


    //TODO model disconnected stuff here!
    // for loads
    int nb_load = loads_p_.size();
    v_kv_from_vpu(Va, Vm, loads_status_, nb_load, loads_bus_id_, id_me_to_solver_, bus_vn_kv_, res_load_v_);
    res_load_p_ = loads_p_;
    res_load_q_ = loads_q_;

    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);


    //TODO model disconnected stuff here!
    // for prods
    int nb_gen = generators_p_.size();
    res_gen_p_ = generators_p_;
    res_gen_v_ = Eigen::VectorXd::Constant(nb_gen, -1.0);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators_status_[gen_id]) continue;
        int bus_id_me = generators_bus_id_(gen_id);
        int bus_solver_id = id_me_to_solver_[bus_id_me];
        if(bus_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::compute_results: A generator is connected to a disconnected bus.");
        }
        double vm_me = Vm(bus_solver_id);
        res_gen_v_(gen_id) = vm_me * bus_vn_kv_(bus_id_me);
    }
    //TODO for res_gen_q_ !!!
}

void GridModel::reset_results(){
    res_load_p_ = Eigen::VectorXd(); // in MW
    res_load_q_ = Eigen::VectorXd(); // in MVar
    res_load_v_ = Eigen::VectorXd(); // in kV

    res_gen_p_ = Eigen::VectorXd();  // in MW
    res_gen_q_ = Eigen::VectorXd();  // in MVar
    res_gen_v_ = Eigen::VectorXd();  // in kV

    powerlines_.reset_results();
    shunts_.reset_results();
    trafos_.reset_results();
}

Eigen::VectorXcd GridModel::dc_pf(const Eigen::VectorXd & p, const Eigen::VectorXcd Va0){
    //TODO fix that with deactivated bus! taking into account refacto !!!

    // initialize the dc Ybus matrix
    Eigen::SparseMatrix<double> dcYbus;
    init_dcY(dcYbus);
    dcYbus.makeCompressed();

    // get the correct matrix : remove the slack bus
    int nb_bus = bus_vn_kv_.size();
    Eigen::SparseMatrix<double> mat_to_inv = Eigen::SparseMatrix<double>(nb_bus-1, nb_bus-1);
    mat_to_inv.reserve(dcYbus.nonZeros());

    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    std::vector<Eigen::Triplet<double> > tripletList;
    for (int k=0; k < nb_bus; ++k){
        if(k == slack_bus_id_) continue;  // I don't add anything to the slack bus
        for (Eigen::SparseMatrix<double>::InnerIterator it(dcYbus, k); it; ++it)
        {
            int row_res = it.row();
            if(row_res == slack_bus_id_) continue;
            row_res = row_res > slack_bus_id_ ? row_res-1 : row_res;
            int col_res = it.col();
            col_res = col_res > slack_bus_id_ ? col_res-1 : col_res;
            tripletList.push_back(Eigen::Triplet<double> (row_res, col_res, it.value()));
        }
    }
    mat_to_inv.setFromTriplets(tripletList.begin(), tripletList.end());

    // solve for theta: P = dcY . theta
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    solver.analyzePattern(mat_to_inv);
    solver.factorize(mat_to_inv);

    // remove the slack bus from the p
    Eigen::VectorXd p_tmp = Eigen::VectorXd(nb_bus-1);

    // TODO vectorize like this, but be carefull to side effect
    int index_tmp = 0;
    for (int k=0; k < nb_bus-1; ++k, ++index_tmp){
        if(k == slack_bus_id_) ++index_tmp;
        p_tmp(k) = p(index_tmp); // - p0_tmp(k);
    }
    // solve the system
    Eigen::VectorXd theta_tmp = solver.solve(p_tmp);
    if(solver.info()!=Eigen::Success) {
        // solving failed
        return Eigen::VectorXd();
    }

    // re insert everything to place
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(nb_bus, 0.);
    index_tmp = 0;
    for (int k=0; k < nb_bus-1; ++k, ++index_tmp){
        if(k == slack_bus_id_) ++index_tmp;
        theta(index_tmp) = theta_tmp(k);
    }
    theta.array() +=  std::arg(Va0(slack_bus_id_));
    Eigen::VectorXd Vm = Va0.array().abs();

    int nb_gen = generators_p_.size();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        int bus_id = generators_bus_id_(gen_id);
        Vm(bus_id) = generators_v_(gen_id);
    }
    return Vm.array() * (theta.array().cos().cast<cdouble>() + 1.0i * theta.array().sin().cast<cdouble>());
}

void GridModel::init_dcY(Eigen::SparseMatrix<double> & dcYbus){
    //TODO handle dc with missing bus taking into account refacto!
    int nb_bus = bus_vn_kv_.size();

    // init this matrix
    Eigen::SparseMatrix<cdouble> tmp = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);

    // fill it properly
    //fillYbusBranch(tmp, false);
    // fillYbusShunt(tmp, false);

    powerlines_.fillYbus(Ybus_, false, id_me_to_solver_);
    shunts_.fillYbus(Ybus_, false, id_me_to_solver_);
    trafos_.fillYbus(Ybus_, false, id_me_to_solver_);
    //fillYbusTrafo(tmp, false);

    // take only real part
    dcYbus  = tmp.real();
}

// deactivate a bus. Be careful, if a bus is deactivated, but an element is
//still connected to it, it will throw an exception
void GridModel::deactivate_bus(int bus_id)
{
    _deactivate(bus_id, bus_status_, need_reset_);
}
void GridModel::reactivate_bus(int bus_id)
{
    _reactivate(bus_id, bus_status_, need_reset_);
}

// for powerline
void GridModel::deactivate_powerline(int powerline_id)
{
    powerlines_.deactivate(powerline_id, need_reset_);
}
void GridModel::reactivate_powerline(int powerline_id)
{
    powerlines_.reactivate(powerline_id, need_reset_);
}
void GridModel::change_bus_powerline_or(int powerline_id, int new_bus_id)
{
    powerlines_.change_bus_or(powerline_id, new_bus_id, need_reset_, bus_vn_kv_.size());
}
void GridModel::change_bus_powerline_ex(int powerline_id, int new_bus_id)
{
    powerlines_.change_bus_ex(powerline_id, new_bus_id, need_reset_, bus_vn_kv_.size());
}

// for trafos
void GridModel::deactivate_trafo(int trafo_id)
{
     trafos_.deactivate(trafo_id, need_reset_); //_deactivate(trafo_id, transformers_status_, need_reset_);
}
void GridModel::reactivate_trafo(int trafo_id)
{
    trafos_.reactivate(trafo_id, need_reset_); //_reactivate(trafo_id, transformers_status_, need_reset_);
}
void GridModel::change_bus_trafo_hv(int trafo_id, int new_bus_id)
{
    trafos_.change_bus_hv(trafo_id, new_bus_id, need_reset_, bus_vn_kv_.size()); //_change_bus(trafo_id, new_bus_id, transformers_bus_hv_id_, need_reset_, bus_vn_kv_.size());
}
void GridModel::change_bus_trafo_lv(int trafo_id, int new_bus_id)
{
    trafos_.change_bus_lv(trafo_id, new_bus_id, need_reset_, bus_vn_kv_.size()); //_change_bus(trafo_id, new_bus_id, transformers_bus_lv_id_, need_reset_, bus_vn_kv_.size());
}
// for loads
void GridModel::deactivate_load(int load_id)
{
    _deactivate(load_id, loads_status_, need_reset_);
}
void GridModel::reactivate_load(int load_id)
{
    _reactivate(load_id, loads_status_, need_reset_);
}
void GridModel::change_bus_load(int load_id, int new_bus_id)
{
    _change_bus(load_id, new_bus_id, loads_bus_id_, need_reset_, bus_vn_kv_.size());
}
// for generators
void GridModel::deactivate_gen(int gen_id)
{
    _deactivate(gen_id, generators_status_, need_reset_);
}
void GridModel::reactivate_gen(int gen_id)
{
    _reactivate(gen_id, generators_status_, need_reset_);
}
void GridModel::change_bus_gen(int gen_id, int new_bus_id)
{
    _change_bus(gen_id, new_bus_id, generators_bus_id_, need_reset_, bus_vn_kv_.size());
}
//for shunts
void GridModel::deactivate_shunt(int shunt_id)
{
    shunts_.deactivate(shunt_id, need_reset_); //_deactivate(shunt_id, shunts_status_, need_reset_);
}
void GridModel::reactivate_shunt(int shunt_id)
{
    shunts_.reactivate(shunt_id, need_reset_); //_reactivate(shunt_id, shunts_status_, need_reset_);
}
void GridModel::change_bus_shunt(int shunt_id, int new_bus_id)
{
    shunts_.change_bus(shunt_id, new_bus_id, need_reset_, bus_vn_kv_.size()); //_change_bus(shunt_id, new_bus_id, shunts_bus_id_, need_reset_, bus_vn_kv_.size());
}