// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATATRAFO_H
#define DATATRAFO_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

class DataTrafo : public DataGeneric
{
    public:
    typedef std::tuple<
               std::vector<double>, // branch_r
               std::vector<double>, // branch_x
               std::vector<std::complex<double> >, // branch_h
               std::vector<int>, // branch_from_id
               std::vector<int>, // branch_to_id
               std::vector<bool> , // status_
               std::vector<double> // ratio_
           >  StateRes;

    DataTrafo() {};

    void init(const Eigen::VectorXd & trafo_r,
                           const Eigen::VectorXd & trafo_x,
                           const Eigen::VectorXcd & trafo_b,
                           const Eigen::VectorXd & trafo_tap_step_pct,
            //                        const Eigen::VectorXd & trafo_tap_step_degree,
                           const Eigen::VectorXd & trafo_tap_pos,
                           const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                           const Eigen::VectorXi & trafo_hv_id,
                           const Eigen::VectorXi & trafo_lv_id
              );
    DataTrafo::StateRes get_state() const;
    void set_state(DataTrafo::StateRes & my_state );

    int nb() { return r_.size(); }

    void deactivate(int trafo_id, bool & need_reset) {_deactivate(trafo_id, status_, need_reset);}
    void reactivate(int trafo_id, bool & need_reset) {_reactivate(trafo_id, status_, need_reset);}
    void change_bus_hv(int trafo_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(trafo_id, new_bus_id, bus_hv_id_, need_reset, nb_bus);}
    void change_bus_lv(int trafo_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(trafo_id, new_bus_id, bus_lv_id_, need_reset, nb_bus);}
    int get_bus_hv(int trafo_id) {return _get_bus(trafo_id, status_, bus_hv_id_);}
    int get_bus_lv(int trafo_id) {return _get_bus(trafo_id, status_, bus_lv_id_);}

    virtual void fillYbus_spmat(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver);
    virtual void fillYbus(std::vector<Eigen::Triplet<cdouble> > & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv);
    void reset_results();
    virtual double get_p_slack(int slack_bus_id);
    virtual void get_q(std::vector<double>& q_by_bus);

    tuple4d get_res_hv() const {return tuple4d(res_p_hv_, res_q_hv_, res_v_hv_, res_a_hv_);}
    tuple4d get_res_lv() const {return tuple4d(res_p_lv_, res_q_lv_, res_v_lv_, res_a_lv_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties
        Eigen::VectorXd r_;
        Eigen::VectorXd x_;
        Eigen::VectorXcd h_;

        // input data
        Eigen::VectorXi bus_hv_id_;
        Eigen::VectorXi bus_lv_id_;
        std::vector<bool> status_;
        Eigen::VectorXd ratio_;

        //output data
        Eigen::VectorXd res_p_hv_;  // in MW
        Eigen::VectorXd res_q_hv_;  // in MVar
        Eigen::VectorXd res_v_hv_;  // in kV
        Eigen::VectorXd res_a_hv_;  // in kA
        Eigen::VectorXd res_p_lv_;  // in MW
        Eigen::VectorXd res_q_lv_;  // in MVar
        Eigen::VectorXd res_v_lv_;  // in kV
        Eigen::VectorXd res_a_lv_;  // in kA
};

#endif  //DATATRAFO_H
