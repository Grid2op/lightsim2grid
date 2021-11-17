// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATASHUNT_H
#define DATASHUNT_H

#include "Utils.h"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"


#include "DataGeneric.h"

/**
This class is a container for all shunts on the grid.

The convention used for the shunt is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html#electric-model
**/
class DataShunt : public DataGeneric
{
    public:
    typedef std::tuple<
           std::vector<real_type>, // p_mw
           std::vector<real_type>, // q_mvar
           std::vector<int>, // bus_id
           std::vector<bool> // status
           >  StateRes;

    DataShunt() {};

    void init(const RealVect & shunt_p_mw,
                     const RealVect & shunt_q_mvar,
                     const Eigen::VectorXi & shunt_bus_id
              );

    // pickle (python)
    DataShunt::StateRes get_state() const;
    void set_state(DataShunt::StateRes & my_state );


    int nb() const { return static_cast<int>(p_mw_.size()); }

    void deactivate(int shunt_id, bool & need_reset) {_deactivate(shunt_id, status_, need_reset);}
    void reactivate(int shunt_id, bool & need_reset) {_reactivate(shunt_id, status_, need_reset);}
    void change_bus(int shunt_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(shunt_id, new_bus_id, bus_id_, need_reset, nb_bus);}
    void change_p(int shunt_id, real_type new_p, bool & need_reset);
    void change_q(int shunt_id, real_type new_q, bool & need_reset);
    int get_bus(int shunt_id) {return _get_bus(shunt_id, status_, bus_id_);}

    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                          bool ac,
                          const std::vector<int> & id_grid_to_solver,
                          real_type sn_mva);
    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);
    virtual void fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);  // in DC i need that

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
        RealVect res_theta_;  // in kV
};

#endif  //DATASHUNT_H
