// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataLine.h"

void DataLine::init(const RealVect & branch_r,
                    const RealVect & branch_x,
                    const CplxVect & branch_h,
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

    bus_or_id_ = branch_from_id;
    bus_ex_id_ = branch_to_id;
    powerlines_h_ = branch_h;
    powerlines_r_ = branch_r;
    powerlines_x_ = branch_x;
    status_ = std::vector<bool>(branch_r.size(), true); // by default everything is connected
}

DataLine::StateRes DataLine::get_state() const
{
     std::vector<real_type> branch_r(powerlines_r_.begin(), powerlines_r_.end());
     std::vector<real_type> branch_x(powerlines_x_.begin(), powerlines_x_.end());
     std::vector<cplx_type > branch_h(powerlines_h_.begin(), powerlines_h_.end());
     std::vector<int > branch_from_id(bus_or_id_.begin(), bus_or_id_.end());
     std::vector<int > branch_to_id(bus_ex_id_.begin(), bus_ex_id_.end());
     std::vector<bool> status = status_;
     DataLine::StateRes res(branch_r, branch_x, branch_h, branch_from_id, branch_to_id, status);
     return res;
}
void DataLine::set_state(DataLine::StateRes & my_state)
{
    reset_results();

    std::vector<real_type> & branch_r = std::get<0>(my_state);
    std::vector<real_type> & branch_x = std::get<1>(my_state);
    std::vector<cplx_type > & branch_h = std::get<2>(my_state);
    std::vector<int> & branch_from_id = std::get<3>(my_state);
    std::vector<int> & branch_to_id = std::get<4>(my_state);
    std::vector<bool> & status = std::get<5>(my_state);
    // TODO check sizes

    // now assign the values
    powerlines_r_ = RealVect::Map(&branch_r[0], branch_r.size());
    powerlines_x_ = RealVect::Map(&branch_x[0], branch_x.size());
    powerlines_h_ = CplxVect::Map(&branch_h[0], branch_h.size());

    // input data
    bus_or_id_ = Eigen::VectorXi::Map(&branch_from_id[0], branch_from_id.size());
    bus_ex_id_ = Eigen::VectorXi::Map(&branch_to_id[0], branch_to_id.size());
    status_ = status;
}

void DataLine::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                        bool ac,
                        const std::vector<int> & id_grid_to_solver,
                        real_type sn_mva)
{
    // fill the matrix
    //TODO template here instead of "if" for ac / dc
    int nb_line = powerlines_r_.size();

    //diagonal coefficients
    for(int line_id =0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        // compute from / to
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }

        // convert subsceptance to half subsceptance, applied on each ends
        cplx_type h = 0.;
        if(ac){
            h = powerlines_h_(line_id); // yes it's the correct one
            h = my_i * my_half_ * h;
        }

        // compute the admittance y
        cplx_type y = 0.;
        cplx_type z = powerlines_x_(line_id);
        if(ac){
            z *= my_i;
            z += powerlines_r_(line_id);
        }
        if (z != my_zero_) y = my_one_ / z;

        // fill non diagonal coefficient
        res.push_back(Eigen::Triplet<cplx_type> (bus_or_solver_id, bus_ex_solver_id, -y));
        res.push_back(Eigen::Triplet<cplx_type> (bus_ex_solver_id, bus_or_solver_id, -y));

        // fill diagonal coefficient
        cplx_type tmp = y;
        if(ac) tmp += h;
        // else tmp += std::imag(h); // todo ???? still don't know !

        res.push_back(Eigen::Triplet<cplx_type> (bus_or_solver_id, bus_or_solver_id, tmp));
        res.push_back(Eigen::Triplet<cplx_type> (bus_ex_solver_id, bus_ex_solver_id, tmp));
    }
}
void DataLine::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver)
{

    //TODO this is no more used!!!! see the other fillYbus
    // fill the matrix
    //TODO template here instead of "if" for ac / dc
    int nb_line = powerlines_r_.size();

    //diagonal coefficients
    for(int line_id =0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        // compute from / to
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataLine::fillYbusBranch: A line is connected (or) to a disconnected bus.");
        }

       //TODO this is no more used!!!! see the other fillYbus
        // convert subsceptance to half subsceptance, applied on each ends
        cplx_type h = 0.;
        if(ac){
            h = powerlines_h_(line_id); // yes it's the correct one
            h = my_i * my_half_ * h;
        }

        // compute the admittance y
        cplx_type y = 0.;
        cplx_type z = powerlines_x_(line_id);
        if(ac){
            z *= my_i;
            z += powerlines_r_(line_id);
        }
        if (z != my_zero_ ) y = my_one_ / z;

        //TODO this is no more used!!!! see the other fillYbus
        // fill non diagonal coefficient
        res.coeffRef(bus_or_solver_id, bus_ex_solver_id) -= y; // * base_for_pu_from;
        res.coeffRef(bus_ex_solver_id, bus_or_solver_id) -= y; // * base_for_pu_to;

        // fill diagonal coefficient
        cplx_type tmp = y;
        if(ac){
            tmp += h;
        }

        //TODO this is no more used!!!! see the other fillYbus
        res.coeffRef(bus_or_solver_id, bus_or_solver_id) += tmp;
        res.coeffRef(bus_ex_solver_id, bus_ex_solver_id) += tmp;

        //TODO this is no more used!!!! see the other fillYbus
    }
}

void DataLine::reset_results()
{
    res_powerline_por_ = RealVect();  // in MW
    res_powerline_qor_ = RealVect();  // in MVar
    res_powerline_vor_ = RealVect();  // in kV
    res_powerline_aor_ = RealVect();  // in kA
    res_powerline_pex_ = RealVect();  // in MW
    res_powerline_qex_ = RealVect();  // in MVar
    res_powerline_vex_ = RealVect();  // in kV
    res_powerline_aex_ = RealVect();  // in kA
}


void DataLine::compute_results(const Eigen::Ref<RealVect> & Va,
                               const Eigen::Ref<RealVect> & Vm,
                               const Eigen::Ref<CplxVect> & V,
                               const std::vector<int> & id_grid_to_solver,
                               const RealVect & bus_vn_kv,
                               real_type sn_mva)
{
    // it needs to be initialized at 0.
    int nb_element = nb();
    res_powerline_por_ = RealVect::Constant(nb_element, my_zero_);  // in MW
    res_powerline_qor_ = RealVect::Constant(nb_element, my_zero_);  // in MVar
    res_powerline_vor_ = RealVect::Constant(nb_element, my_zero_);  // in kV
    res_powerline_aor_ = RealVect::Constant(nb_element, my_zero_);  // in kA
    res_powerline_pex_ = RealVect::Constant(nb_element, my_zero_);  // in MW
    res_powerline_qex_ = RealVect::Constant(nb_element, my_zero_);  // in MVar
    res_powerline_vex_ = RealVect::Constant(nb_element, my_zero_);  // in kV
    res_powerline_aex_ = RealVect::Constant(nb_element, my_zero_);  // in kA
    res_powerline_thetaor_ = RealVect::Constant(nb_element, my_zero_);  // in kV
    res_powerline_thetaex_ = RealVect::Constant(nb_element, my_zero_);  // in kV
    for(int line_id = 0; line_id < nb_element; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;

        //physical properties
        real_type r = powerlines_r_(line_id);
        real_type x = powerlines_x_(line_id);
        cplx_type h = my_i * my_half_ * powerlines_h_(line_id);
        cplx_type y = my_one_ / (r + my_i * x);

        // connectivity
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::res_powerlines: A powerline or a trafo is connected (or) to a disconnected bus.");
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataModel::res_powerlines: A powerline or a trafo is connected (ex) to a disconnected bus.");
        }

        // results of the powerflow
        cplx_type Eor = V(bus_or_solver_id);
        cplx_type Eex = V(bus_ex_solver_id);

        // powerline equations
        cplx_type I_orex = (y + h) * Eor - y * Eex;
        cplx_type I_exor = (y + h) * Eex - y * Eor;

        I_orex = std::conj(I_orex);
        I_exor = std::conj(I_exor);
        cplx_type s_orex = Eor * I_orex;
        cplx_type s_exor = Eex * I_exor;

        res_powerline_por_(line_id) = std::real(s_orex) * sn_mva;
        res_powerline_qor_(line_id) = std::imag(s_orex) * sn_mva;
        res_powerline_pex_(line_id) = std::real(s_exor) * sn_mva;
        res_powerline_qex_(line_id) = std::imag(s_exor) * sn_mva;

        // retrieve voltages magnitude in kv instead of pu
        real_type v_or = Vm(bus_or_solver_id);
        real_type v_ex = Vm(bus_ex_solver_id);
        real_type bus_vn_kv_or = bus_vn_kv(bus_or_id_me);
        real_type bus_vn_kv_ex = bus_vn_kv(bus_ex_id_me);
        res_powerline_vor_(line_id) = v_or * bus_vn_kv_or;
        res_powerline_vex_(line_id) = v_ex * bus_vn_kv_ex;

        res_powerline_thetaor_(line_id) = Va(bus_or_solver_id) * 180. / my_pi;
        res_powerline_thetaex_(line_id) = Va(bus_ex_solver_id) * 180. / my_pi;
    }
    _get_amps(res_powerline_aor_, res_powerline_por_, res_powerline_qor_, res_powerline_vor_);
    _get_amps(res_powerline_aex_, res_powerline_pex_, res_powerline_qex_, res_powerline_vex_);
}

real_type DataLine::get_p_slack(int slack_bus_id)
{
    int nb_element = nb();
    real_type res = 0.;
    for(int line_id = 0; line_id < nb_element; ++line_id)
    {
        if(!status_[line_id]) continue;
        if(bus_or_id_(line_id) == slack_bus_id) res += res_powerline_por_(line_id);
        if(bus_ex_id_(line_id) == slack_bus_id) res += res_powerline_pex_(line_id);
    }
    return res;
}

void DataLine::get_q(std::vector<real_type>& q_by_bus)
{
    int nb_element = nb();
    for(int el_id = 0; el_id < nb_element; ++el_id)
    {
        if(!status_[el_id]) continue;
        int bus_id_ex = bus_ex_id_[el_id];
        int bus_id_or = bus_or_id_[el_id];
        q_by_bus[bus_id_or] += res_powerline_qor_(el_id);
        q_by_bus[bus_id_ex] += res_powerline_qex_(el_id);
    }
}
