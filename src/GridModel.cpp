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

Eigen::VectorXcd GridModel::ac_pf(const Eigen::VectorXcd & Vinit,
                                  int max_iter,
                                  double tol)
{
    // TODO optimization when it's not mandatory to start from scratch
    bool conv = false;
    Eigen::VectorXcd res = Eigen::VectorXcd();
    if(need_reset_){
        init_Ybus(Ybus_, Sbus_, id_me_to_solver_, id_solver_to_me_, slack_bus_id_solver_);
        fillYbus(Ybus_, true, id_me_to_solver_);
        fillpv_pq(id_me_to_solver_);
        _solver.reset();
    }

    fillSbus(Sbus_, true, id_me_to_solver_, slack_bus_id_solver_);

    // extract only connected bus from Vinit
    int nb_bus_solver = id_solver_to_me_.size();
    Eigen::VectorXcd V = Eigen::VectorXcd::Constant(id_solver_to_me_.size(), 1.04);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        int bus_me_id = id_solver_to_me_[bus_solver_id];  //POSSIBLE SEGFAULT
        cdouble tmp = Vinit(bus_me_id);
        V(bus_solver_id) = tmp;
        // TODO save this V somewhere
    }

    conv = _solver.do_newton(Ybus_, V, Sbus_, bus_pv_, bus_pq_, max_iter, tol);
    if (conv){
        compute_results();
        need_reset_ = false;
        res = _solver.get_V();
    } else {
        //powerflow diverge
        reset_results();
        need_reset_ = true;  // in this case, the powerflow diverge, so i need to recompute Ybus next time
    }
    return res;
};

void GridModel::init_Ybus(Eigen::SparseMatrix<cdouble> & Ybus, Eigen::VectorXcd & Sbus,
                          std::vector<int>& id_me_to_solver, std::vector<int>& id_solver_to_me,
                          int & slack_bus_id_solver){
    //TODO get disconnected bus !!! (and have some conversion for it)
    //1. init the conversion bus
    int nb_bus_init = bus_vn_kv_.size();
    id_me_to_solver = std::vector<int>(nb_bus_init, _deactivated_bus_id);  // by default, if a bus is disconnected, then it has a -1 there
    id_solver_to_me = std::vector<int>();
    id_solver_to_me.reserve(bus_vn_kv_.size());
    int bus_id_solver=0;
    for(int bus_id_me=0; bus_id_me < nb_bus_init; ++bus_id_me){
        if(bus_status_[bus_id_me]){
            // bus is connected
            id_solver_to_me.push_back(bus_id_me);
            id_me_to_solver[bus_id_me] = bus_id_solver;
            ++bus_id_solver;
        }
    }
    int nb_bus = id_me_to_solver.size();

    Ybus = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);
    Ybus.reserve(nb_bus + 2*powerlines_.nb() + 2*trafos_.nb());

    // init diagonal coefficients
    for(int bus_id=0; bus_id < nb_bus; ++bus_id){
        // assign diagonal coefficient
        Ybus.insert(bus_id, bus_id) = 0.;
    }

    Sbus = Eigen::VectorXcd::Constant(nb_bus, 0.);
    slack_bus_id_solver = id_me_to_solver[slack_bus_id_];
    if(slack_bus_id_solver == _deactivated_bus_id){
        //TODO improve error message with the gen_id
        throw std::runtime_error("The slack bus is disconnected.");
    }
}

void GridModel::fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int>& id_me_to_solver){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/
    // init the Ybus matrix
    powerlines_.fillYbus(res, ac, id_me_to_solver);
    shunts_.fillYbus(res, ac, id_me_to_solver);
    trafos_.fillYbus(res, ac, id_me_to_solver);
    loads_.fillYbus(res, ac, id_me_to_solver);
    generators_.fillYbus(res, ac, id_me_to_solver);
}

void GridModel::fillSbus(Eigen::VectorXcd & res, bool ac, const std::vector<int>& id_me_to_solver, int slack_bus_id_solver)
{
    // init the Sbus vector
    powerlines_.fillSbus(res, ac, id_me_to_solver);
    shunts_.fillSbus(res, ac, id_me_to_solver);
    trafos_.fillSbus(res, ac, id_me_to_solver);
    loads_.fillSbus(res, ac, id_me_to_solver);
    generators_.fillSbus(res, ac, id_me_to_solver);

    // handle slack bus
    double sum_active = res.sum().real();
    res.coeffRef(slack_bus_id_solver) -= sum_active;
}

void GridModel::fillpv_pq(const std::vector<int>& id_me_to_solver)
{
    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    int nb_bus = id_solver_to_me_.size();  // number of bus in the solver!
    std::vector<int> bus_pq;
    std::vector<int> bus_pv;
    std::vector<bool> has_bus_been_added(nb_bus, false);

    powerlines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    shunts_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    trafos_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
    loads_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver_, id_me_to_solver);
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
     //TODO check it has converged!

    // retrieve results from powerflow
    const auto & Va = _solver.get_Va();
    const auto & Vm = _solver.get_Vm();
    const auto & V = _solver.get_V();

    // for powerlines
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);

    // for trafo
    trafos_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);


    //TODO model disconnected stuff here!
    // for loads
    loads_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);

    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);


    //TODO model disconnected stuff here!
    // for prods
    generators_.compute_results(Va, Vm, V, id_me_to_solver_, bus_vn_kv_);
    //TODO for res_gen_q_ !!!
}

void GridModel::reset_results(){
    powerlines_.reset_results();
    shunts_.reset_results();
    trafos_.reset_results();
    loads_.reset_results();
    generators_.reset_results();
}

Eigen::VectorXcd GridModel::dc_pf(const Eigen::VectorXcd & Vinit,
                                  int max_iter,  // not used for DC
                                  double tol  // not used for DC
                                  )
{
    // TODO refactor that with ac pf, this is mostly done, but only mostly...
    int nb_bus = bus_vn_kv_.size();
    Eigen::SparseMatrix<cdouble> dcYbus_tmp;
    Eigen::VectorXcd Sbus_tmp;
    std::vector<int> id_me_to_solver;
    std::vector<int> id_solver_to_me;
    int slack_bus_id_solver;

    if(need_reset_){
        init_Ybus(dcYbus_tmp, Sbus_tmp, id_me_to_solver, id_solver_to_me, slack_bus_id_solver);
        fillYbus(dcYbus_tmp, false, id_me_to_solver);
        fillpv_pq(id_me_to_solver);
    }
    fillSbus(Sbus_tmp, false, id_me_to_solver, slack_bus_id_solver);


    // extract only connected bus from Vinit
    int nb_bus_solver = id_solver_to_me_.size();
    Eigen::VectorXcd V = Eigen::VectorXcd::Constant(id_solver_to_me_.size(), 1.04);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        int bus_me_id = id_solver_to_me_[bus_solver_id];  //POSSIBLE SEGFAULT
        cdouble tmp = Vinit(bus_me_id);
        V(bus_solver_id) = tmp;
    }


    // DC SOLVER STARTS HERE
    // remove the slack bus
    // TODO this should rather be one in a "dc solver" instead of here

    // remove the slack bus from Ybus
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    Eigen::SparseMatrix<double> dcYbus = Eigen::SparseMatrix<double>(nb_bus - 1, nb_bus - 1);
    std::vector<Eigen::Triplet<double> > tripletList;
    for (int k=0; k < nb_bus; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        for (Eigen::SparseMatrix<cdouble>::InnerIterator it(dcYbus_tmp, k); it; ++it)
        {
            int row_res = it.row();
            if(row_res == slack_bus_id_solver) continue;
            row_res = row_res > slack_bus_id_solver ? row_res - 1 : row_res;
            int col_res = it.col();
            col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
            tripletList.push_back(Eigen::Triplet<double> (row_res, col_res, std::real(it.value())));
        }
    }
    dcYbus.setFromTriplets(tripletList.begin(), tripletList.end());
    dcYbus.makeCompressed();

    // remove the slack bus from Sbus
    Eigen::VectorXd Sbus = Eigen::VectorXd(nb_bus - 1);
    for (int k=0; k < nb_bus; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        int col_res = k;
        col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
        Sbus.coeffRef(col_res) = std::real(Sbus_tmp(k));
    }

    // initialize the solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    solver.analyzePattern(dcYbus);
    solver.factorize(dcYbus);


    // solve for theta: Sbus = dcY . theta
    Eigen::VectorXd Va_dc = solver.solve(Sbus);
    if(solver.info() != Eigen::Success) {
        // solving failed, this should not happen in dc ...
        return Eigen::VectorXcd();
    }
    // until there it's ok

    // retrieve back the results in the proper shape
    int nb_bus_me = bus_vn_kv_.size();
    int bus_id_solver;
    Eigen::VectorXd Va = Eigen::VectorXd::Constant(nb_bus_me, 0.);
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
    Eigen::VectorXd Vm = Vinit.array().abs();  // fill Vm = Vinit for all
    // put Vm = 0. for disconnected bus
    for (int bus_id_me=0; bus_id_me < nb_bus_me; ++bus_id_me){
        if(bus_status_[bus_id_me]) continue;  // nothing is done if the bus is connected
        Vm(bus_id_me) = 0.;
    }
    //END of the SOLVER PART

    //TODO handle Vm = Vm (gen) for connected generators
    return Vm.array() * (Va.array().cos().cast<cdouble>() + 1.0i * Va.array().sin().cast<cdouble>());
    //return Eigen::VectorXcd();
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
    loads_.deactivate(load_id, need_reset_); // _deactivate(load_id, loads_status_, need_reset_);
}
void GridModel::reactivate_load(int load_id)
{
    loads_.reactivate(load_id, need_reset_); //_reactivate(load_id, loads_status_, need_reset_);
}
void GridModel::change_bus_load(int load_id, int new_bus_id)
{
    loads_.change_bus(load_id, new_bus_id, need_reset_, bus_vn_kv_.size()); //_change_bus(load_id, new_bus_id, loads_bus_id_, need_reset_, bus_vn_kv_.size());
}
// for generators
void GridModel::deactivate_gen(int gen_id)
{
    generators_.deactivate(gen_id, need_reset_); // _deactivate(gen_id, generators_status_, need_reset_);
}
void GridModel::reactivate_gen(int gen_id)
{
    generators_.reactivate(gen_id, need_reset_);  // _reactivate(gen_id, generators_status_, need_reset_);
}
void GridModel::change_bus_gen(int gen_id, int new_bus_id)
{
    generators_.change_bus(gen_id, new_bus_id, need_reset_, bus_vn_kv_.size()); //_change_bus(gen_id, new_bus_id, generators_bus_id_, need_reset_, bus_vn_kv_.size());
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