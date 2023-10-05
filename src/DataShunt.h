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
    // iterators part
    public:
        class ShuntInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                bool connected;
                int bus_id;

                real_type target_p_mw;
                real_type target_q_mvar;
                bool has_res;
                real_type res_p_mw;
                real_type res_q_mvar;
                real_type res_v_kv;
                real_type res_theta_deg;

                ShuntInfo(const DataShunt & r_data_shunt, int my_id):
                id(-1),
                connected(false),
                bus_id(-1),
                target_p_mw(0.),
                target_q_mvar(0.),
                has_res(false),
                res_p_mw(0.),
                res_q_mvar(0.),
                res_v_kv(0.),
                res_theta_deg(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_shunt.nb()))
                    {
                        id = my_id;
                        connected = r_data_shunt.status_[my_id];
                        bus_id = r_data_shunt.bus_id_[my_id];

                        target_p_mw = r_data_shunt.p_mw_.coeff(my_id);
                        target_q_mvar = r_data_shunt.q_mvar_.coeff(my_id);

                        has_res = r_data_shunt.res_p_.size() > 0;
                        if(has_res)
                        {
                            res_p_mw = r_data_shunt.res_p_.coeff(my_id);
                            res_q_mvar = r_data_shunt.res_q_.coeff(my_id);
                            res_v_kv = r_data_shunt.res_v_.coeff(my_id);
                            res_theta_deg = r_data_shunt.res_theta_.coeff(my_id);
                        }
                    }
                }
        };
        typedef ShuntInfo DataInfo;

    private:
        typedef DataConstIterator<DataShunt> DataShuntConstIterator;

    public:
        typedef DataShuntConstIterator const_iterator_type;
        const_iterator_type begin() const {return DataShuntConstIterator(this, 0); }
        const_iterator_type end() const {return DataShuntConstIterator(this, nb()); }
        ShuntInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative load id.");
            }
            if(id >= nb())
            {
                throw std::range_error("Generator out of bound. Not enough loads on the grid.");
            }
            return ShuntInfo(*this, id);
        }

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
                          real_type sn_mva) const;
    virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                            std::vector<Eigen::Triplet<real_type> > & Bpp,
                            const std::vector<int> & id_grid_to_solver,
                            real_type sn_mva,
                            FDPFMethod xb_or_bx) const;
    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);
    virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const;  // in DC i need that

    void compute_results(const Eigen::Ref<const RealVect> & Va,
                         const Eigen::Ref<const RealVect> & Vm,
                         const Eigen::Ref<const CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva,
                         bool ac);
    void reset_results();

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
