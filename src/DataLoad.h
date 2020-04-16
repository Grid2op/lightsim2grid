// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATALOAD_H
#define DATALOAD_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

class DataLoad : public DataGeneric
{
    public:
    DataLoad() {};

    void init(const Eigen::VectorXd & loads_p,
              const Eigen::VectorXd & loads_q,
              const Eigen::VectorXi & loads_bus_id
              );

    int nb() { return p_mw_.size(); }

    void deactivate(int load_id, bool & need_reset) {_deactivate(load_id, status_, need_reset);}
    void reactivate(int load_id, bool & need_reset) {_reactivate(load_id, status_, need_reset);}
    void change_bus(int load_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(load_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    int get_bus(int load_id) {return _get_bus(load_id, status_, bus_id_);}
    void change_p(int load_id, double new_p, bool & need_reset);
    void change_q(int load_id, double new_q, bool & need_reset);

    virtual void fillSbus(Eigen::VectorXcd & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                         const Eigen::Ref<Eigen::VectorXd> & Vm,
                         const Eigen::Ref<Eigen::VectorXcd> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const Eigen::VectorXd & bus_vn_kv);
    void reset_results();
    virtual double get_p_slack(int slack_bus_id);
    virtual void get_q(std::vector<double>& q_by_bus);

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties

        // input data
        Eigen::VectorXd p_mw_;
        Eigen::VectorXd q_mvar_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        Eigen::VectorXd res_p_;  // in MW
        Eigen::VectorXd res_q_;  // in MVar
        Eigen::VectorXd res_v_;  // in kV
};

#endif  //DATALOAD_H
