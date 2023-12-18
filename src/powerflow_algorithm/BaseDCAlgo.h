// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_DC_ALGO_H
#define BASE_DC_ALGO_H

#include "BaseAlgo.h"

template<class LinearSolver>
class BaseDCAlgo: public BaseAlgo
{
    public:
        BaseDCAlgo():
            BaseAlgo(false),
            _linear_solver(),
            need_factorize_(true),
            sizeYbus_with_slack_(0),
            sizeYbus_without_slack_(0){};

        ~BaseDCAlgo(){}

        virtual void reset();

        // TODO SLACK : this should be handled in Sbus by the gridmodel maybe ?
        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                        CplxVect & V,
                        const CplxVect & Sbus,
                        const Eigen::VectorXi & slack_ids,
                        const RealVect & slack_weights,  // currently unused
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        );

        virtual RealMat get_ptdf(const Eigen::SparseMatrix<cplx_type> & dcYbus);  // TODO PTDF
        virtual Eigen::SparseMatrix<real_type> get_lodf();  // TODO PTDF
        virtual Eigen::SparseMatrix<real_type> get_bsdf();  // TODO PTDF

    private:
        // no copy allowed
        BaseDCAlgo( const BaseDCAlgo & ) =delete ;
        BaseDCAlgo & operator=( const BaseDCAlgo & ) =delete;

    protected:
        void fill_mat_bus_id(int nb_bus_solver);
        void fill_dcYbus_noslack(int nb_bus_solver, const Eigen::SparseMatrix<cplx_type> & ref_mat);

        // remove_slack_buses: res_mat is initialized and make_compressed in this function
        template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
        void remove_slack_buses(int nb_bus_solver, const Eigen::SparseMatrix<ref_mat_type> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat);

    protected:
        LinearSolver  _linear_solver;
        bool need_factorize_;

        // save this not to recompute them when not needed
        int sizeYbus_with_slack_;
        int sizeYbus_without_slack_;
        RealVect dcSbus_noslack_;
        Eigen::SparseMatrix<real_type> dcYbus_noslack_;
        Eigen::VectorXi my_pv_;
        Eigen::VectorXi slack_buses_ids_solver_;
        // -1 if bus is slack , else the id of the row / column used in the linear solver representing this bus
        Eigen::VectorXi mat_bus_id_;   // formerly `ybus_to_me`

};

#include "BaseDCAlgo.tpp"

#endif // BASE_DC_ALGO_H
