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
        class TrafoInfo
        {
            public:
                // members
                int id;  // id of the generator
                bool connected;
                int bus_hv_id;
                int bus_lv_id;
                real_type r;
                real_type x;
                cplx_type h;
                bool is_tap_hv_side;
                real_type ratio;
                real_type shift;

                bool has_res;
                real_type res_p_hv_mw;
                real_type res_q_hv_mvar;
                real_type res_v_hv_kv;
                real_type res_a_hv_a;
                real_type res_p_lv_mw;
                real_type res_q_lv_mvar;
                real_type res_v_lv_kv;
                real_type res_a_lv_a;

                TrafoInfo(const DataTrafo & r_data_trafo, int my_id):
                id(-1),
                connected(false),
                bus_hv_id(-1),
                bus_lv_id(-1),
                r(-1.0),
                x(-1.0),
                h({0., 0.}),
                is_tap_hv_side(true),
                ratio(-1.0),
                shift(-1.0),
                has_res(false),
                res_p_hv_mw(0.),
                res_q_hv_mvar(0.),
                res_v_hv_kv(0.),
                res_a_hv_a(0.),
                res_p_lv_mw(0.),
                res_q_lv_mvar(0.),
                res_v_lv_kv(0.),
                res_a_lv_a(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_trafo.nb()))
                    {
                        id = my_id;
                        connected = r_data_trafo.status_[my_id];
                        bus_hv_id = r_data_trafo.bus_hv_id_.coeff(my_id);
                        bus_lv_id = r_data_trafo.bus_lv_id_.coeff(my_id);
                        r = r_data_trafo.r_.coeff(my_id);
                        x = r_data_trafo.x_.coeff(my_id);
                        h = r_data_trafo.h_.coeff(my_id);
                        is_tap_hv_side = r_data_trafo.is_tap_hv_side_[my_id];
                        ratio = r_data_trafo.ratio_.coeff(my_id);
                        shift = r_data_trafo.shift_.coeff(my_id);

                        has_res = r_data_trafo.res_p_hv_.size() > 0;
                        if(has_res)
                        {
                            res_p_hv_mw = r_data_trafo.res_p_hv_.coeff(my_id);
                            res_q_hv_mvar = r_data_trafo.res_q_hv_.coeff(my_id);
                            res_v_hv_kv = r_data_trafo.res_v_hv_.coeff(my_id);
                            res_a_hv_a = r_data_trafo.res_a_hv_.coeff(my_id);
                            res_p_lv_mw = r_data_trafo.res_p_lv_.coeff(my_id);
                            res_q_lv_mvar = r_data_trafo.res_q_lv_.coeff(my_id);
                            res_v_lv_kv = r_data_trafo.res_v_lv_.coeff(my_id);
                            res_a_lv_a = r_data_trafo.res_a_lv_.coeff(my_id);
                        }
                    }
                }
        };
        typedef TrafoInfo DataInfo;

    private:
        typedef DataConstIterator<DataTrafo> DataTrafoConstIterator;

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
    //pickle
    DataTrafo::StateRes get_state() const;
    void set_state(DataTrafo::StateRes & my_state );

    int nb() const { return r_.size(); }

    // make it iterable
    typedef DataTrafoConstIterator const_iterator_type;
    const_iterator_type begin() const {return DataTrafoConstIterator(this, 0); }
    const_iterator_type end() const {return DataTrafoConstIterator(this, nb()); }
    TrafoInfo operator[](int id) const
    {
        if(id < 0)
        {
            throw std::range_error("You cannot ask for a negative generator");
        }
        if(id >= nb())
        {
            throw std::range_error("Generator out of bound. Not enough generator on the grid.");
        }
        return TrafoInfo(*this, id);
    }

    // method used within lightsim
    void deactivate(int trafo_id, bool & need_reset) {_deactivate(trafo_id, status_, need_reset);}
    void reactivate(int trafo_id, bool & need_reset) {_reactivate(trafo_id, status_, need_reset);}
    void change_bus_hv(int trafo_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(trafo_id, new_bus_id, bus_hv_id_, need_reset, nb_bus);}
    void change_bus_lv(int trafo_id, int new_bus_id, bool & need_reset, int nb_bus) {_change_bus(trafo_id, new_bus_id, bus_lv_id_, need_reset, nb_bus);}
    int get_bus_hv(int trafo_id) {return _get_bus(trafo_id, status_, bus_hv_id_);}
    int get_bus_lv(int trafo_id) {return _get_bus(trafo_id, status_, bus_lv_id_);}

    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);
    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                          bool ac,
                          const std::vector<int> & id_grid_to_solver,
                          real_type sn_mva);
    virtual void fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver);  // needed for dc mode

    void compute_results(const Eigen::Ref<RealVect> & Va,
                         const Eigen::Ref<RealVect> & Vm,
                         const Eigen::Ref<CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva);
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
