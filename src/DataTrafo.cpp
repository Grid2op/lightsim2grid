// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of PyKLU2Grid, PyKLU2Grid a implements a c++ backend targeting the Grid2Op platform.

#include "DataTrafo.h"

void DataTrafo::init(const Eigen::VectorXd & trafo_r,
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

    r_ = trafo_r;
    x_ = trafo_x;
    h_ = trafo_b;
    ratio_ = ratio;
    bus_hv_id_ = trafo_hv_id;
    bus_lv_id_ = trafo_lv_id;
    status_ = std::vector<bool>(trafo_r.size(), true);
}

void DataTrafo::fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver)
{
    //TODO merge that with fillYbusBranch!
    //TODO template here instead of "if"
    int nb_trafo = nb();
    cdouble my_i = 1.0i;
    for(int trafo_id =0; trafo_id < nb_trafo; ++trafo_id){
        // i don't do anything if the trafo is disconnected
        if(!status_[trafo_id]) continue;

        // compute from / to
        int bus_hv_id_me = bus_hv_id_(trafo_id);
        int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
        if(bus_hv_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::fillYbusTrafo: A trafo is connected (hv) to a disconnected bus.");
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::fillYbusTrafo: A trafo is connected (lv) to a disconnected bus.");
        }

        // get the transformers ratio
        double r = ratio_(trafo_id);

        // subsecptance
        cdouble h = 0.;
        if(ac){
            h = h_(trafo_id);
            h = my_i * 0.5 * h;
        }

        // admittance
        cdouble y = 0.;
        cdouble z = x_(trafo_id);
        if(ac){
            z *= my_i;
            z += r_(trafo_id);
        }
        if(z != 0.) y = 1.0 / z;

        // fill non diagonal coefficient
        cdouble tmp = y / r;
        res.coeffRef(bus_hv_solver_id, bus_lv_solver_id) -= tmp ;
        res.coeffRef(bus_lv_solver_id, bus_hv_solver_id) -= tmp;

        // fill diagonal coefficient
        if(!ac){
            r = 1.0; // in dc, r = 1.0 here (same voltage both side)
        }
        tmp += h;
        res.coeffRef(bus_hv_solver_id, bus_hv_solver_id) += tmp / r ;
        res.coeffRef(bus_lv_solver_id, bus_lv_solver_id) += tmp * r;
    }
}

void DataTrafo::fillYbus(std::vector<Eigen::Triplet<cdouble> > & res, bool ac, const std::vector<int> & id_grid_to_solver)
{
    //TODO merge that with fillYbusBranch!
    //TODO template here instead of "if"
    int nb_trafo = nb();
    cdouble my_i = 1.0i;
    for(int trafo_id =0; trafo_id < nb_trafo; ++trafo_id){
        // i don't do anything if the trafo is disconnected
        if(!status_[trafo_id]) continue;

        // compute from / to
        int bus_hv_id_me = bus_hv_id_(trafo_id);
        int bus_hv_solver_id = id_grid_to_solver[bus_hv_id_me];
        if(bus_hv_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::fillYbusTrafo: A trafo is connected (hv) to a disconnected bus.");
        }
        int bus_lv_id_me = bus_lv_id_(trafo_id);
        int bus_lv_solver_id = id_grid_to_solver[bus_lv_id_me];
        if(bus_lv_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::fillYbusTrafo: A trafo is connected (lv) to a disconnected bus.");
        }

        // get the transformers ratio
        double r = ratio_(trafo_id);

        // subsecptance
        cdouble h = 0.;
        if(ac){
            h = h_(trafo_id);
            h = my_i * 0.5 * h;
        }

        // admittance
        cdouble y = 0.;
        cdouble z = x_(trafo_id);
        if(ac){
            z *= my_i;
            z += r_(trafo_id);
        }
        if(z != 0.) y = 1.0 / z;

        // fill non diagonal coefficient
        cdouble tmp = y / r;
        res.push_back(Eigen::Triplet<cdouble> (bus_hv_solver_id, bus_lv_solver_id, -tmp));
        res.push_back(Eigen::Triplet<cdouble> (bus_lv_solver_id, bus_hv_solver_id, -tmp));

        // fill diagonal coefficient
        if(!ac){
            r = 1.0; // in dc, r = 1.0 here (same voltage both side)
        }
        tmp += h;
        res.push_back(Eigen::Triplet<cdouble>(bus_hv_solver_id, bus_hv_solver_id, tmp / r));
        res.push_back(Eigen::Triplet<cdouble>(bus_lv_solver_id, bus_lv_solver_id, tmp * r));
    }
}

void DataTrafo::compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv
                              )
{
    // it needs to be initialized at 0.
    int nb_element = nb();
    res_p_hv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    res_q_hv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    res_v_hv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    res_a_hv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    res_p_lv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MW
    res_q_lv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in MVar
    res_v_lv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kV
    res_a_lv_ = Eigen::VectorXd::Constant(nb_element, 0.0);  // in kA
    cdouble my_i = 1.0i;
    for(int line_id = 0; line_id < nb_element; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;

        //physical properties
        double r = r_(line_id);
        double x = x_(line_id);
        double ratio_me = ratio_(line_id);
        cdouble h = my_i * 0.5 * h_(line_id);
        cdouble y = 1.0 / (r + my_i * x);
        y /= ratio_me;

        // connectivity
        int bus_or_id_me = bus_hv_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataTrafo::compute_results: A trafo is connected (hv) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_lv_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataTrafo::compute_results: A trafo is connected (lv) to a disconnected bus.");
        }

        // results of the powerflow
        cdouble Eor = V(bus_or_solver_id);
        cdouble Eex = V(bus_ex_solver_id);

        // powerline equations
        cdouble I_orex = (y + h) / ratio_me * Eor - y * Eex;
        cdouble I_exor = (y + h) * ratio_me * Eex - y * Eor;

        I_orex = std::conj(I_orex);
        I_exor = std::conj(I_exor);
        cdouble s_orex = Eor * I_orex;
        cdouble s_exor = Eex * I_exor;

        res_p_hv_(line_id) = std::real(s_orex);
        res_q_hv_(line_id) = std::imag(s_orex);
        res_p_lv_(line_id) = std::real(s_exor);
        res_q_lv_(line_id) = std::imag(s_exor);

        // retrieve voltages magnitude in kv instead of pu
        double v_or = Vm(bus_or_solver_id);
        double v_ex = Vm(bus_ex_solver_id);
        double bus_vn_kv_or = bus_vn_kv(bus_or_id_me);
        double bus_vn_kv_ex = bus_vn_kv(bus_ex_id_me);
        res_v_hv_(line_id) = v_or * bus_vn_kv_or;
        res_v_lv_(line_id) = v_ex * bus_vn_kv_ex;
    }
    _get_amps(res_a_hv_, res_p_hv_, res_q_hv_, res_v_hv_);
    _get_amps(res_a_lv_, res_p_lv_, res_q_lv_, res_v_lv_);
}

void DataTrafo::reset_results(){
    res_p_hv_ = Eigen::VectorXd();  // in MW
    res_q_hv_ = Eigen::VectorXd();  // in MVar
    res_v_hv_ = Eigen::VectorXd();  // in kV
    res_a_hv_ = Eigen::VectorXd();  // in kA
    res_p_lv_ = Eigen::VectorXd();  // in MW
    res_q_lv_ = Eigen::VectorXd();  // in MVar
    res_v_lv_ = Eigen::VectorXd();  // in kV
    res_a_lv_ = Eigen::VectorXd();  // in kA
}

double DataTrafo::get_p_slack(int slack_bus_id)
{
    int nb_element = nb();
    double res = 0.;
    for(int line_id = 0; line_id < nb_element; ++line_id)
    {
        if(!status_[line_id]) continue;
        if(bus_hv_id_(line_id) == slack_bus_id) res += res_p_hv_(line_id);
        if(bus_lv_id_(line_id) == slack_bus_id) res += res_p_lv_(line_id);
    }
    return res;
}

void DataTrafo::get_q(std::vector<double>& q_by_bus)
{
    int nb_element = nb();
    for(int el_id = 0; el_id < nb_element; ++el_id)
    {
        if(!status_[el_id]) continue;
        int bus_id_hv = bus_hv_id_[el_id];
        int bus_id_lv = bus_lv_id_[el_id];
        q_by_bus[bus_id_hv] += res_q_hv_(el_id);
        q_by_bus[bus_id_lv] += res_q_lv_(el_id);
    }
}
