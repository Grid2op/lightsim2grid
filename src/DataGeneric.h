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

#include "Utils.h"

/**
Base class for every object that can be manipulated
**/
class DataGeneric
{
    public:

        virtual void fillYbus(std::vector<Eigen::Triplet<cdouble> > & res, bool ac, const std::vector<int> & id_grid_to_solver) {};
        virtual void fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver) {};
        virtual void fillSbus(Eigen::VectorXcd & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){};
        virtual void fillpv(std::vector<int>& bus_pv,
                            std::vector<bool> & has_bus_been_added,
                            int slack_bus_id_solver,
                            const std::vector<int> & id_grid_to_solver) {};
        virtual double get_p_slack(int slack_bus_id) {return 0.;}
        virtual void set_p_slack(int gen_slackbus, double p_slack) {};
        virtual void get_q(const std::vector<double>& q_by_bus) {};

    protected:
        static const int _deactivated_bus_id;
        static const cdouble my_i;

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
        void _get_amps(Eigen::VectorXd & a, const Eigen::VectorXd & p, const Eigen::VectorXd & q, const Eigen::VectorXd & v);

        /**
        **/
        void v_kv_from_vpu(const Eigen::Ref<Eigen::VectorXd> & Va,
                           const Eigen::Ref<Eigen::VectorXd> & Vm,
                           const std::vector<bool> & status,
                           int nb_element,
                           const Eigen::VectorXi & bus_me_id,
                           const std::vector<int> & id_grid_to_solver,
                           const Eigen::VectorXd & bus_vn_kv,
                           Eigen::VectorXd & v);

};

#endif // DATAGENERIC_H

