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
    bus_pu_ = Eigen::VectorXd::Constant(nb_bus, 1.0 / sn_mva_);
    bus_pu_.array() *= bus_vn_kv_.array() * bus_vn_kv_.array(); // np.square(base_kv) / net.sn_mva

    // init diagonal coefficients
    for(int bus_id=0; bus_id < nb_bus; ++bus_id){
        // assign diagonal coefficient
        Ybus_.insert(bus_id, bus_id) = 0.;
    }

    Sbus_ = Eigen::VectorXcd::Constant(nb_bus, 0.);
}

void DataModel::init_powerlines(const Eigen::VectorXd & branch_r,
                     const Eigen::VectorXd & branch_x,
                     const Eigen::VectorXd & branch_c,
//                             const Eigen::VectorXd & branch_g,
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

    int nb_line = branch_r.size();

    powerlines_bus_or_id_ = branch_from_id;
    powerlines_bus_ex_id_ = branch_to_id;

    Eigen::VectorXd bus_pu_from = Eigen::VectorXd(nb_line); // (from_id);
    for(int i = 0; i < nb_line; ++i) {bus_pu_from(i) = bus_pu_(branch_from_id(i));}

    powerlines_r_ = branch_r.array() / bus_pu_from.array();
    powerlines_x_ = branch_x.array() / bus_pu_from.array();

    powerlines_h_ = Eigen::VectorXcd::Constant(nb_line, 2.0 * f_hz_ * M_PI * 1e-9);
    powerlines_h_.array() *= branch_c.array();
    powerlines_h_.array() *=  bus_pu_from.cast<cdouble>().array();
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

std::tuple<Eigen::VectorXd,
           Eigen::VectorXd,
           Eigen::VectorXcd>
           DataModel::get_trafo_param(const Eigen::VectorXd & trafo_vn_hv,
                           const Eigen::VectorXd & trafo_vn_lv,
                           const Eigen::VectorXd & trafo_vk_percent,
                           const Eigen::VectorXd & trafo_vkr_percent,
                           const Eigen::VectorXd & trafo_sn_trafo_mva,
                           const Eigen::VectorXd & trafo_pfe_kw,
                           const Eigen::VectorXd & trafo_i0_pct,
                           const Eigen::VectorXi & trafo_lv_id)
{
    //TODO only for "trafo model = t"
    //TODO supposes that the step start at 0 for "no ratio"

    //TODO consistency: move this class outside of here
    int nb_trafo = trafo_vn_lv.size();

    Eigen::VectorXd vn_trafo_lv = trafo_vn_lv;
    Eigen::VectorXd vn_lv = Eigen::VectorXd(nb_trafo);  //bus_vn_kv_[trafo_lv_id]; // TODO check if it compiles
    for(int i = 0; i<nb_trafo; ++i) {vn_lv(i) = bus_vn_kv_(trafo_lv_id(i));}  //TODO optimize that

    // compute r and x
    Eigen::VectorXd tmp = vn_trafo_lv.array() / vn_lv.array();
    tmp = tmp.array() * tmp.array();
    Eigen::VectorXd tap_lv = tmp * sn_mva_;
    Eigen::VectorXd _1_sn_trafo_mva = 1.0 / trafo_sn_trafo_mva.array();
    Eigen::VectorXd z_sc = 0.01 * trafo_vk_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    Eigen::VectorXd r_sc = 0.01 * trafo_vkr_percent.array() * _1_sn_trafo_mva.array() * tap_lv.array();
    Eigen::VectorXd tmp2 = z_sc.array()*z_sc.array() - r_sc.array() * r_sc.array();
    Eigen::VectorXd x_sc = z_sc.cwiseSign().array() * tmp2.cwiseSqrt().array();

    // compute h, the subsceptance
    Eigen::VectorXd baseR = Eigen::VectorXd(nb_trafo);
    for(int i = 0; i<nb_trafo; ++i) {baseR(i) = bus_pu_(trafo_lv_id(i));}  //TODO optimize that
    Eigen::VectorXd pfe =  trafo_pfe_kw.array() * 1e-3;

    // Calculate subsceptance ###
    Eigen::VectorXd vnl_squared = trafo_vn_lv.array() * trafo_vn_lv.array();
    Eigen::VectorXd b_real = pfe.array() / vnl_squared.array() * baseR.array();
    tmp2 = (trafo_i0_pct.array() * 0.01 * trafo_sn_trafo_mva.array());
    Eigen::VectorXd b_img =  tmp2.array() * tmp2.array() - pfe.array() * pfe.array();

    for(int i = 0; i<nb_trafo; ++i) {if (b_img(i) < 0.)  b_img(i) = 0.;}
    b_img = b_img.cwiseSqrt();
    b_img.array() *= baseR.array() / vnl_squared.array();
    Eigen::VectorXcd y = - 1.0i * b_real.array().cast<cdouble>() - b_img.array().cast<cdouble>() * trafo_i0_pct.cwiseSign().array();
    Eigen::VectorXcd b_sc = y.array() / tmp.array();

    //transform trafo from t model to pi model, of course...
    // (remove that if trafo model is not t, but directly pi)
    cdouble my_i = 1.0i;
    for(int i = 0; i<nb_trafo; ++i){
        if(b_sc(i) == 0.) continue;
        cdouble za_star = 0.5 * (r_sc(i) + my_i * x_sc(i));
        cdouble zc_star = - my_i / b_sc(i);
        cdouble zSum_triangle = za_star * za_star + 2.0 * za_star * zc_star;
        cdouble zab_triangle = zSum_triangle / zc_star;
        cdouble zbc_triangle = zSum_triangle / za_star;

        r_sc(i) = zab_triangle.real();
        x_sc(i) = zab_triangle.imag();
        b_sc(i) = -2.0 * my_i / zbc_triangle;
    }

    std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd> res =
        std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXcd>(std::move(r_sc), std::move(x_sc), std::move(b_sc));
    return res;
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

void DataModel::compute_newton(){
    init_Ybus();
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
    double sum_active = 0.;
    double tmp;
    for(int gen_id = 0; gen_id < generators_p_.size(); ++gen_id){
        bus_id = generators_bus_id_(gen_id);
        tmp = generators_p_(gen_id);
        Sbus_.coeffRef(bus_id) += tmp; // + my_i * generators_p_(gen_id);
        sum_active += tmp;
    }
    for(int load_id = 0; load_id < loads_p_.size(); ++load_id){
        bus_id = loads_bus_id_(load_id);
        tmp = loads_p_(load_id);
        Sbus_.coeffRef(bus_id) -= tmp + my_i * loads_q_(load_id);
        sum_active -= tmp;
    }

    Sbus_.coeffRef(slack_bus_id_) -= sum_active;
    //TODO put the shunt here (but test before it can be done)
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

