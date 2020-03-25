#include "DataModel.h"

void DataModel::init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo){
    /**
    initialize the bus_vn_kv_ member
    and
    initialize the Ybus_ matrix at the proper shape
    **/
    int nb_bus = bus_vn_kv.size();
    bus_vn_kv_ = bus_vn_kv;  // base_kv
    Ybus_ = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);
    Ybus_.reserve(nb_bus + 2*nb_line + 2*nb_trafo);

    // per unit conversion
    // bus_pu_ = Eigen::VectorXd::Constant(nb_bus, 1.0 / sn_mva_);
    // bus_pu_.array() *= bus_vn_kv_.array() * bus_vn_kv_.array(); // np.square(base_kv) / net.sn_mva

    // init diagonal coefficients
    for(int bus_id=0; bus_id < nb_bus; ++bus_id){
        // assign diagonal coefficient
        Ybus_.insert(bus_id, bus_id) = 0.;
    }

    Sbus_ = Eigen::VectorXcd::Constant(nb_bus, 0.);
}

void DataModel::init_powerlines(const Eigen::VectorXd & branch_r,
                                const Eigen::VectorXd & branch_x,
                                const Eigen::VectorXcd & branch_h,
                                const Eigen::VectorXi & branch_from_id,
                                const Eigen::VectorXi & branch_to_id
                                )
{
    /**
    This method initialize the Ybus matrix from the branch matrix.
    It has to be called once when the solver is initialized. Afterwards, a call to
    "updateYbus" should be made for performance optimiaztion instead. //TODO
    **/

    // TODO check what can be checked: branch_* have same size, no id in branch_to_id that are
    // TODO not in [0, .., buv_vn_kv.size()] etc.

    //TODO consistency with trafo: have a converter methods to convert this value into pu, and store the pu
    // in this method

    powerlines_bus_or_id_ = branch_from_id;
    powerlines_bus_ex_id_ = branch_to_id;
    powerlines_h_ = branch_h;
    powerlines_r_ = branch_r;
    powerlines_x_ = branch_x;
}

void DataModel::init_shunt(const Eigen::VectorXd & shunt_p_mw,
                const Eigen::VectorXd & shunt_q_mvar,
                const Eigen::VectorXi & shunt_bus_id)
{
    /**
    supposes the matrix Ybus_ has already been initialized
    **/

    shunts_p_mw_ = shunt_p_mw;
    shunts_q_mvar_ = shunt_q_mvar;
    shunts_bus_id_ = shunt_bus_id;
}

void DataModel::init_trafo(const Eigen::VectorXd & trafo_r,
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
    //TODO "parrallel" in the pandapower dataframe, like for lines, are not handled. Handle it python side!
    Eigen::VectorXd ratio = 1.0 + 0.01 * trafo_tap_step_pct.array() * trafo_tap_pos.array() * (2*trafo_tap_hv.array().cast<double>() - 1.0);

    transformers_r_ = trafo_r;
    transformers_x_ = trafo_x;
    transformers_h_ = trafo_b;
    transformers_ratio_ = ratio;
    transformers_bus_hv_id_ = trafo_hv_id;
    transformers_bus_lv_id_ = trafo_lv_id;
}


void DataModel::init_generators(const Eigen::VectorXd & generators_p,
                     const Eigen::VectorXd & generators_v,
                     const Eigen::VectorXi & generators_bus_id)
{
    generators_p_ = generators_p;
    generators_v_ = generators_v;
    generators_bus_id_ = generators_bus_id;
}

void DataModel::init_loads(const Eigen::VectorXd & loads_p,
                const Eigen::VectorXd & loads_q,
                const Eigen::VectorXi & loads_bus_id)
{
    loads_p_ = loads_p;
    loads_q_ = loads_q;
    loads_bus_id_ = loads_bus_id;
}

void DataModel::fillYbusBranch(Eigen::SparseMatrix<cdouble> & res, bool ac)
{
    // fill the matrix
    //TODO template here instead of "if"
    int nb_line = powerlines_r_.size();
    cdouble my_i = 1.0i;

    //diagonal coefficients
    for(int line_id =0; line_id < nb_line; ++line_id){
        // get the from / to bus id
        int from_id = powerlines_bus_or_id_(line_id);
        int to_id = powerlines_bus_ex_id_(line_id);

        // convert subsceptance to half subsceptance, applied on each ends
        cdouble h = 0.;
        if(ac){
            h = powerlines_h_(line_id); // yes it's the correct one
            h = my_i * 0.5 * h;
        }

        // compute the admittance y
        cdouble y = 0.;
        cdouble z = powerlines_x_(line_id);
        if(ac){
            z *= my_i;
            z += powerlines_r_(line_id);
        }
        if (z !=0. ) y = 1.0 / z;

        // fill non diagonal coefficient
        res.coeffRef(from_id, to_id) -= y; // * base_for_pu_from;
        res.coeffRef(to_id, from_id) -= y; // * base_for_pu_to;

        // fill diagonal coefficient
        cdouble tmp = y + h;
        res.coeffRef(from_id, from_id) += tmp;
        res.coeffRef(to_id, to_id) += tmp;
    }
}

void DataModel::fillYbusShunt(Eigen::SparseMatrix<cdouble> & res, bool ac){
    int nb_shunt = shunts_q_mvar_.size();
    cdouble tmp;
    int bus_id;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // assign diagonal coefficient
        tmp = shunts_p_mw_(shunt_id) + 1.0i * shunts_q_mvar_(shunt_id);
        bus_id = shunts_bus_id_(shunt_id);
        res.coeffRef(bus_id, bus_id) -= tmp;
    }
}

void DataModel::fillYbusTrafo(Eigen::SparseMatrix<cdouble> & res, bool ac){
    //TODO merge that with fillYbusBranch!
    //TODO template here instead of "if"
    int nb_trafo = transformers_bus_hv_id_.size();
    cdouble my_i = 1.0i;
    for(int trafo_id =0; trafo_id < nb_trafo; ++trafo_id){
        // compute from / to
        int hv_id = transformers_bus_hv_id_(trafo_id);
        int lv_id = transformers_bus_lv_id_(trafo_id);

        // get the transformers ratio
        double r = transformers_ratio_(trafo_id);

        // subsecptance
        cdouble h = 0.;
        if(ac){
            h = transformers_h_(trafo_id);
            h = my_i * 0.5 * h;
        }

        // admittance
        cdouble y = 0.;
        cdouble z = transformers_x_(trafo_id);
        if(ac){
            z *= my_i;
            z += transformers_r_(trafo_id);
        }
        if(z != 0.) y = 1.0 / z;

        // fill non diagonal coefficient
        cdouble tmp = y / r;
        res.coeffRef(hv_id, lv_id) -= tmp ;
        res.coeffRef(lv_id, hv_id) -= tmp;

        // fill diagonal coefficient
        if(!ac){
            r = 1.0; // in dc, r = 1.0 here (same voltage both side)
        }
        res.coeffRef(hv_id, hv_id) += (tmp + h)  / r ;
        res.coeffRef(lv_id, lv_id) += (tmp + h) * r;
    }
}

Eigen::VectorXcd DataModel::dc_pf(const Eigen::VectorXd & p, const Eigen::VectorXcd Va0){
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
//    Eigen::VectorXd p0_tmp = mat_to_inv.col(slack_bus_id_) * Va0(slack_bus_id_); //TODO check that in pandapower
    Eigen::VectorXd p_tmp = Eigen::VectorXd(nb_bus-1);

    // TODO vectorize like this, but be carefull to side effect
//    theta_tmp.segment(0, slack_bus_id_-1) = p.segment(0,slack_bus_id_-1);
//    theta_tmp.segment(slack_bus_id_-1, nb_bus-1 - (slack_bus_id_-1)) = p.segment(0,slack_bus_id_-1);
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

bool DataModel::compute_newton(Eigen::VectorXcd & V,
                               int max_iter,
                               double tol)
{
    init_Ybus();
    _solver.reset();
    auto Sbus = get_Sbus();
    auto pv = get_pv();
    auto pq = get_pq();
    auto Ybus = get_Ybus();
    bool res = _solver.do_newton(Ybus, V, Sbus, pv, pq, max_iter, tol);
    if (res) compute_results();
    else reset_results();
    return res;
};

void DataModel::init_Ybus(){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/

    // init the Ybus matrix
    fillYbusBranch(Ybus_, true);
    fillYbusShunt(Ybus_, true);
    fillYbusTrafo(Ybus_, true);

    // init the Sbus vector
    cdouble my_i = 1.0i;
    int bus_id;
    int nb_bus = bus_vn_kv_.size();
    int nb_load = loads_p_.size();
    int nb_gen = generators_p_.size();
    double sum_active = 0.;
    double tmp;

    std::vector<bool> has_gen_conn(nb_bus, false);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        bus_id = generators_bus_id_(gen_id);
        tmp = generators_p_(gen_id);
        Sbus_.coeffRef(bus_id) += tmp; // + my_i * generators_p_(gen_id);
        sum_active += tmp;
        has_gen_conn[bus_id] = true;
    }
    for(int load_id = 0; load_id < nb_load; ++load_id){
        bus_id = loads_bus_id_(load_id);
        tmp = loads_p_(load_id);
        Sbus_.coeffRef(bus_id) -= tmp + my_i * loads_q_(load_id);
        sum_active -= tmp;
    }

    Sbus_.coeffRef(slack_bus_id_) -= sum_active;
    //TODO put the shunt here (but test before it can be done)

    // TODO i probably need to keep track of the order here !!!!
    std::vector<int> bus_pq;
    std::vector<int> bus_pv;
    std::vector<bool> has_bus_been_added(nb_bus, false);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        bus_id = generators_bus_id_(gen_id);
        if(bus_id == slack_bus_id_) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id]) continue; // i already added this bus
        bus_pv.push_back(bus_id);
        has_bus_been_added[bus_id] = true;
    }
    for(int bus_id = 0; bus_id< nb_bus; ++bus_id){
        if(bus_id == slack_bus_id_) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pv.data(), bus_pv.size());
    bus_pq_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(bus_pq.data(), bus_pq.size());
}

void DataModel::init_dcY(Eigen::SparseMatrix<double> & dcYbus){
    int nb_bus = bus_vn_kv_.size();

    // init this matrix
    Eigen::SparseMatrix<cdouble> tmp = Eigen::SparseMatrix<cdouble>(nb_bus, nb_bus);

    // fill it properly
    fillYbusBranch(tmp, false);
    fillYbusShunt(tmp, false);
    fillYbusTrafo(tmp, false);

    // take only real part
    dcYbus  = tmp.real();
}

void DataModel::compute_results(){
     //TODO check it has converged!

    // retrieve results from powerflow
    const auto & Va = _solver.get_Va();
    const auto & Vm = _solver.get_Vm();
    const auto & V = _solver.get_V();

    // for powerlines
    int nb_line = powerlines_r_.size();
    res_powerlines(Va, Vm, V, nb_line,
                   powerlines_r_,
                   powerlines_x_,
                   powerlines_h_,
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

    // for trafo
    int nb_trafo = transformers_r_.size();
    res_powerlines(Va, Vm, V, nb_trafo,
                   transformers_r_,
                   transformers_x_,
                   transformers_h_,
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

    // for loads
    int nb_load = loads_p_.size();
    res_loads(Va, Vm, nb_load, loads_bus_id_, res_load_v_);
    res_load_p_ = loads_p_;
    res_load_q_ = loads_q_;

    // for shunts
    int nb_shunt = shunts_p_mw_.size();
    res_loads(Va, Vm, nb_shunt, shunts_bus_id_, res_shunt_v_);
    res_shunt_p_ = shunts_p_mw_;
    res_shunt_q_ = shunts_q_mvar_;

    // for prods
    res_gen_p_ = generators_p_;
    res_gen_v_ = generators_v_;
    //TODO for res_gen_q_ !!!
}

void DataModel::_get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v){
    const double _1_sqrt_3 = 1.0 / std::sqrt(3.);
    Eigen::VectorXd p2q2 = p.array() * p.array() + q.array() * q.array();
    p2q2 = p2q2.array().cwiseSqrt();
    a = p2q2.array() * _1_sqrt_3 / v.array();
}

void DataModel::res_powerlines(const Eigen::Ref<Eigen::VectorXd> & Va,
                               const Eigen::Ref<Eigen::VectorXd> & Vm,
                               const Eigen::Ref<Eigen::VectorXcd> & V,
                               int nb_element,
                               const Eigen::VectorXd & el_r,
                               const Eigen::VectorXd & el_x,
                               const Eigen::VectorXcd & el_h,
                               const Eigen::VectorXi & bus_or_id_,
                               const Eigen::VectorXi & bus_ex_id_,
                               Eigen::VectorXd & por,  // in MW
                               Eigen::VectorXd & qor,  // in MVar
                               Eigen::VectorXd & vor,  // in kV
                               Eigen::VectorXd & aor,  // in kA
                               Eigen::VectorXd & pex,  // in MW
                               Eigen::VectorXd & qex,  // in MVar
                               Eigen::VectorXd & vex,  // in kV
                               Eigen::VectorXd & aex  // in kA
                              ){
    por = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    qor = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    vor = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    aor = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    pex = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    qex = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    vex = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    aex = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    cdouble my_i = 1.0i;
    for(int line_id = 0; line_id < nb_element; ++line_id){
        //physical properties
        double r = el_r(line_id);
        double x = el_x(line_id);
        cdouble h = my_i * 0.5 * el_h(line_id);
        cdouble y = 1.0 / (r + my_i * x);

        // double g = std::real(tmp);
        // double b = std::imag(tmp);
        // std::cout << " for powerline " << line_id << std::end;
        // std::cout << "\t r " << r << " x " << x << " g " << g << " b " << b << std::endl;

        // connectivity
        int bus_or_id = bus_or_id_(line_id);
        int bus_ex_id = bus_ex_id_(line_id);

        // results of the powerflow
        // double theta_or = Va(bus_or_id);
        // double theta_ex = Va(bus_ex_id);
        // double v_or = Vm(bus_or_id);
        // double v_ex = Vm(bus_ex_id);
        cdouble Eor = V(bus_or_id);
        cdouble Eex = V(bus_ex_id);

        //std::cout << "\t theta_or " << theta_or << " theta_ex " << theta_ex << " v_or " << v_or << " v_ex " << v_ex << std::endl;

        // tmp element (be smart in computations)
        // double vkj = v_or * v_ex * sn_mva_;
        // double cos_thetakj = std::cos(theta_or - theta_ex);
        // double sin_thetakj = std::sin(theta_or - theta_ex);

        // powerline equations
        cdouble I_orex = y * (Eor - Eex) + h * Eor;
        cdouble I_exor = y * (Eex - Eor) + h * Eex;

        I_orex = std::conj(I_orex);
        I_exor = std::conj(I_exor);
        cdouble s_orex = Eor * I_orex;
        cdouble s_exor = Eex * I_exor;
        // TODO these are probably not the right formula, need to check
        por(line_id) = std::real(s_orex);
        qor(line_id) = std::imag(s_orex);
        pex(line_id) = std::real(s_exor);
        qex(line_id) = std::imag(s_exor);

        // retrieve voltages magnitude in kv instead of pu
        double v_or = Vm(bus_or_id);
        double v_ex = Vm(bus_ex_id);
        double bus_vn_kv_or = bus_vn_kv_(bus_or_id);
        double bus_vn_kv_ex = bus_vn_kv_(bus_ex_id);
        vor(line_id) = v_or * bus_vn_kv_or;
        vex(line_id) = v_ex * bus_vn_kv_ex;
    }
    _get_amps(aor, por, qor, vor);
    _get_amps(aex, pex, qex, vex);
}

void DataModel::res_loads(const Eigen::Ref<Eigen::VectorXd> & Va,
                          const Eigen::Ref<Eigen::VectorXd> & Vm,
                          int nb_element,
                          const Eigen::VectorXi & bus_id,
                          Eigen::VectorXd & v){
    v = Eigen::VectorXd::Constant(nb_element, 0.0);
    for(int el_id = 0; el_id < nb_element; ++el_id){
        int el_bus_id = bus_id(el_id);
        double bus_vn_kv = bus_vn_kv_(el_bus_id);
        v(el_id) = Vm(el_id) * bus_vn_kv;
    }
}

void DataModel::reset_results(){
    res_load_p_ = Eigen::VectorXd(); // in MW
    res_load_q_ = Eigen::VectorXd(); // in MVar
    res_load_v_ = Eigen::VectorXd(); // in kV

    res_gen_p_ = Eigen::VectorXd();  // in MW
    res_gen_q_ = Eigen::VectorXd();  // in MVar
    res_gen_v_ = Eigen::VectorXd();  // in kV

    res_powerline_por_ = Eigen::VectorXd();  // in MW
    res_powerline_qor_ = Eigen::VectorXd();  // in MVar
    res_powerline_vor_ = Eigen::VectorXd();  // in kV
    res_powerline_aor_ = Eigen::VectorXd();  // in kA
    res_powerline_pex_ = Eigen::VectorXd();  // in MW
    res_powerline_qex_ = Eigen::VectorXd();  // in MVar
    res_powerline_vex_ = Eigen::VectorXd();  // in kV
    res_powerline_aex_ = Eigen::VectorXd();  // in kA

    res_trafo_por_ = Eigen::VectorXd();  // in MW
    res_trafo_qor_ = Eigen::VectorXd();  // in MVar
    res_trafo_vor_ = Eigen::VectorXd();  // in kV
    res_trafo_aor_ = Eigen::VectorXd();  // in kA
    res_trafo_pex_ = Eigen::VectorXd();  // in MW
    res_trafo_qex_ = Eigen::VectorXd();  // in MVar
    res_trafo_vex_ = Eigen::VectorXd();  // in kV
    res_trafo_aex_ = Eigen::VectorXd();  // in kA

    res_shunt_p_ = Eigen::VectorXd();  // in MW
    res_shunt_q_ = Eigen::VectorXd();  // in MVar
    res_shunt_v_ = Eigen::VectorXd();  // in kV
}
