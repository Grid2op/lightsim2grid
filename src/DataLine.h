// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATALINE_H
#define DATALINE_H

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "DataGeneric.h"

/**
This class is a container for all the powerlines on the grid.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/line.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/line.html#electric-model

**/
class DataLine : public DataGeneric
{
    public:
    typedef std::tuple<
               std::vector<real_type>, // branch_r
               std::vector<real_type>, // branch_x
               std::vector<cplx_type>, // branch_h
               std::vector<int>, // branch_from_id
               std::vector<int>, // branch_to_id
               std::vector<bool> // status_
               >  StateRes;

    DataLine() {};

    void init(const RealVect & branch_r,
              const RealVect & branch_x,
              const CplxVect & branch_h,
              const Eigen::VectorXi & branch_from_id,
              const Eigen::VectorXi & branch_to_id
              );

    // pickle
    DataLine::StateRes get_state() const;
    void set_state(DataLine::StateRes & my_state );
    template<class T>
    void check_size(const T& my_state)
    {
        //currently unused
        unsigned int size_th = 6;
        if (my_state.size() != size_th)
        {
            std::cout << "LightSim::DataLine state size " << my_state.size() << " instead of "<< size_th << std::endl;
            // TODO more explicit error message
            throw std::runtime_error("Invalid state when loading LightSim::DataLine");
        }
    }

    int nb() { return powerlines_r_.size(); }

    void deactivate(int powerline_id, bool & need_reset) {_deactivate(powerline_id, status_, need_reset);}
    void reactivate(int powerline_id, bool & need_reset) {_reactivate(powerline_id, status_, need_reset);}
    void change_bus_or(int powerline_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(powerline_id, new_bus_id, bus_or_id_, need_reset, nb_bus);}
    void change_bus_ex(int powerline_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(powerline_id, new_bus_id, bus_ex_id_, need_reset, nb_bus);}
    int get_bus_or(int powerline_id) {return _get_bus(powerline_id, status_, bus_or_id_);}
    int get_bus_ex(int powerline_id) {return _get_bus(powerline_id, status_, bus_ex_id_);}
    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                          bool ac,
                          const std::vector<int> & id_grid_to_solver,
                          real_type sn_mva
                          );
    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<RealVect> & Va,
                         const Eigen::Ref<RealVect> & Vm,
                         const Eigen::Ref<CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva);
    void reset_results();
    virtual real_type get_p_slack(int slack_bus_id);
    virtual void get_q(std::vector<real_type>& q_by_bus);

    tuple4d get_lineor_res() const {return tuple4d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_);}
    tuple4d get_lineex_res() const {return tuple4d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_);}
    const std::vector<bool>& get_status() const {return status_;}

    protected:
        // physical properties
        RealVect powerlines_r_;
        RealVect powerlines_x_;
        CplxVect powerlines_h_;

        // input data
        Eigen::VectorXi bus_or_id_;
        Eigen::VectorXi bus_ex_id_;
        std::vector<bool> status_;

        //output data
        RealVect res_powerline_por_;  // in MW
        RealVect res_powerline_qor_;  // in MVar
        RealVect res_powerline_vor_;  // in kV
        RealVect res_powerline_aor_;  // in kA
        RealVect res_powerline_pex_;  // in MW
        RealVect res_powerline_qex_;  // in MVar
        RealVect res_powerline_vex_;  // in kV
        RealVect res_powerline_aex_;  // in kA
};

#endif  //DATALINE_H
