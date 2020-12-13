// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DATAGENERIC_H
#define DATAGENERIC_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "BaseConstants.h"


// TODO make this class iterable ! with operator begin, end and an iterator
/**
Base class for every object that can be manipulated
**/
class DataGeneric : public BaseConstants
{
    public:

        virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res, bool ac, const std::vector<int> & id_grid_to_solver) {};
        virtual void fillYbus(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver) {};
        virtual void fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){};
        virtual void fillpv(std::vector<int>& bus_pv,
                            std::vector<bool> & has_bus_been_added,
                            int slack_bus_id_solver,
                            const std::vector<int> & id_grid_to_solver) {};
        virtual real_type get_p_slack(int slack_bus_id) {return my_zero_;}
        virtual void set_p_slack(int gen_slackbus, real_type p_slack) {};
        virtual void get_q(std::vector<real_type>& q_by_bus) {};

    protected:
        static const int _deactivated_bus_id;

        /**
        activation / deactivation of elements
        **/
        void _reactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _deactivate(int el_id, std::vector<bool> & status, bool & need_reset);
        void _change_bus(int el_id, int new_bus_me_id, Eigen::VectorXi & el_bus_ids, bool & need_reset, int nb_bus);
        int _get_bus(int el_id, const std::vector<bool> & status_, const Eigen::VectorXi & bus_id_);

        /**
        compute the amps from the p, the q and the v (v should NOT be pair unit)
        **/
        void _get_amps(RealVect & a, const RealVect & p, const RealVect & q, const RealVect & v);

        /**
        **/
        void v_kv_from_vpu(const Eigen::Ref<RealVect> & Va,
                           const Eigen::Ref<RealVect> & Vm,
                           const std::vector<bool> & status,
                           int nb_element,
                           const Eigen::VectorXi & bus_me_id,
                           const std::vector<int> & id_grid_to_solver,
                           const RealVect & bus_vn_kv,
                           RealVect & v);
        /**
        check the size of the elements
        **/
        template<class T>
        void check_size(const T & container, int size, const std::string & container_name)
        {
            if(container.size() != size) throw std::runtime_error(container_name + " do not have the proper size");
        }
};

#endif // DATAGENERIC_H

