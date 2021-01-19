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

/**
This class is a container for all loads on the grid.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/load.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/load.html#electric-model

NOTE: this class is also used for the storage units! So storage units are modeled as load
which entails that negative storage: the unit is discharging, power is injected in the grid,
positive storage: the unit is charging, power is taken from the grid.
**/
class DataLoad : public DataGeneric
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    public:
    typedef std::tuple<
       std::vector<real_type>, // p_mw
       std::vector<real_type>, // q_mvar
       std::vector<int>, // bus_id
       std::vector<bool> // status
       >  StateRes;

    DataLoad() {};

    // pickle (python)
    DataLoad::StateRes get_state() const;
    void set_state(DataLoad::StateRes & my_state );


    void init(const RealVect & loads_p,
              const RealVect & loads_q,
              const Eigen::VectorXi & loads_bus_id
              );

    int nb() { return p_mw_.size(); }

    void deactivate(int load_id, bool & need_reset) {_deactivate(load_id, status_, need_reset);}
    void reactivate(int load_id, bool & need_reset) {_reactivate(load_id, status_, need_reset);}
    void change_bus(int load_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(load_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    int get_bus(int load_id) {return _get_bus(load_id, status_, bus_id_);}
    void change_p(int load_id, real_type new_p, bool & need_reset);
    void change_q(int load_id, real_type new_q, bool & need_reset);

    virtual void fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<RealVect> & Va,
                         const Eigen::Ref<RealVect> & Vm,
                         const Eigen::Ref<CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv);
    void reset_results();
    virtual real_type get_p_slack(int slack_bus_id);
    virtual void get_q(std::vector<real_type>& q_by_bus);

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties

        // input data
        RealVect p_mw_;
        RealVect q_mvar_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        RealVect res_p_;  // in MW
        RealVect res_q_;  // in MVar
        RealVect res_v_;  // in kV
};

#endif  //DATALOAD_H
