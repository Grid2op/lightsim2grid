// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATASGEN_H
#define DATASGEN_H

#include "Utils.h"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "DataGeneric.h"

/**
This class is a container for all static generator (PQ generators) on the grid.
They are given in the generator convention: positive sign for P,Q => the power is produced.

The convention used for the static is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html#electric-model
**/
class DataSGen: public DataGeneric
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    public:
    typedef std::tuple<
       std::vector<real_type>, // p_mw
       std::vector<real_type>, // q_mvar
       std::vector<real_type>, // p_min
       std::vector<real_type>, //  p_max
       std::vector<real_type>, //  q_min
       std::vector<real_type>, //  q_max
       std::vector<int>, // bus_id
       std::vector<bool> // status
       >  StateRes;

    DataSGen() {};

    // pickle (python)
    DataSGen::StateRes get_state() const;
    void set_state(DataSGen::StateRes & my_state );


    void init(const RealVect & sgen_p,
              const RealVect & sgen_q,
              const RealVect & sgen_pmin,
              const RealVect & sgen_pmax,
              const RealVect & sgen_qmin,
              const RealVect & sgen_qmax,
              const Eigen::VectorXi & sgen_bus_id
              );

    int nb() const { return static_cast<int>(p_mw_.size()); }

    void deactivate(int sgen_id, bool & need_reset) {_deactivate(sgen_id, status_, need_reset);}
    void reactivate(int sgen_id, bool & need_reset) {_reactivate(sgen_id, status_, need_reset);}
    void change_bus(int sgen_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(sgen_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    int get_bus(int sgen_id) {return _get_bus(sgen_id, status_, bus_id_);}
    void change_p(int sgen_id, real_type new_p, bool & need_reset);
    void change_q(int sgen_id, real_type new_q, bool & need_reset);

    virtual void fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<const RealVect> & Va,
                         const Eigen::Ref<const RealVect> & Vm,
                         const Eigen::Ref<const CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva,
                         bool ac);
    void reset_results();
    // TODO SLACK real_type p_slack OR  std::set<real_type> p_slack ???
    virtual real_type get_p_slack(const std::vector<int>& slack_bus_id) const;
    virtual void get_q(std::vector<real_type>& q_by_bus);

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
    const std::vector<bool>& get_status() const {return status_;}
    const Eigen::VectorXi & get_bus_id() const {return bus_id_;}

    protected:
        // physical properties
        RealVect p_min_mw_;
        RealVect p_max_mw_;
        RealVect q_min_mvar_;
        RealVect q_max_mvar_;

        // input data
        RealVect p_mw_;
        RealVect q_mvar_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        RealVect res_p_;  // in MW
        RealVect res_q_;  // in MVar
        RealVect res_v_;  // in kV
        RealVect res_theta_;  // in degree
};

#endif  //DATASGEN_H
