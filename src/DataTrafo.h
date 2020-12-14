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

/**
This class is a container for all transformers on the grid.
Transformers are modeled "in pi" here. If your trafo are given in a "t" model (like in pandapower
for example) use the DataConverter class.

The convention used for the transformer is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html#electric-model
**/
class DataTrafo : public DataGeneric
{
    public:
    typedef std::tuple<
               std::vector<real_type>, // branch_r
               std::vector<real_type>, // branch_x
               std::vector<cplx_type >, // branch_h
               std::vector<int>, // branch_from_id
               std::vector<int>, // branch_to_id
               std::vector<bool> , // status_
               std::vector<real_type>, // ratio_
               std::vector<bool> , // is_tap_hv_side
               std::vector<real_type> // shift_
           >  StateRes;

    DataTrafo() {};

    void init(const RealVect & trafo_r,
                           const RealVect & trafo_x,
                           const CplxVect & trafo_b,
                           const RealVect & trafo_tap_step_pct,
            //                        const RealVect & trafo_tap_step_degree,
                           const RealVect & trafo_tap_pos,
                           const RealVect & trafo_shift_degree,
                           const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
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

    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);
    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<RealVect> & Va,
                         const Eigen::Ref<RealVect> & Vm,
                         const Eigen::Ref<CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv);
    void reset_results();
    virtual real_type get_p_slack(int slack_bus_id);
    virtual void get_q(std::vector<real_type>& q_by_bus);

    tuple4d get_res_hv() const {return tuple4d(res_p_hv_, res_q_hv_, res_v_hv_, res_a_hv_);}
    tuple4d get_res_lv() const {return tuple4d(res_p_lv_, res_q_lv_, res_v_lv_, res_a_lv_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties
        RealVect r_;
        RealVect x_;
        CplxVect h_;
        std::vector<bool> is_tap_hv_side_;  // whether the tap is hav side or not

        // input data
        Eigen::VectorXi bus_hv_id_;
        Eigen::VectorXi bus_lv_id_;
        std::vector<bool> status_;
        RealVect ratio_;  // transformer ratio
        RealVect shift_;  // phase shifter (in radian !)

        //output data
        RealVect res_p_hv_;  // in MW
        RealVect res_q_hv_;  // in MVar
        RealVect res_v_hv_;  // in kV
        RealVect res_a_hv_;  // in kA
        RealVect res_p_lv_;  // in MW
        RealVect res_q_lv_;  // in MVar
        RealVect res_v_lv_;  // in kV
        RealVect res_a_lv_;  // in kA
};

#endif  //DATATRAFO_H
