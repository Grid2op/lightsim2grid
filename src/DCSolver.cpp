// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DCSolver.h"

bool DCSolver::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                          CplxVect & V,
                          const CplxVect & Sbus,
                          const Eigen::VectorXi & pv,
                          const Eigen::VectorXi & pq,
                          int max_iter,
                          real_type tol
                          )
{
    // max_iter is ignored
    // tol is ignored
    // V is used the following way: at pq buses it's completely ignored. For pv bus only the magnitude is used,
    //   and for the slack bus both the magnitude and the angle are used.

    auto timer = CustTimer();
    const int nb_bus_solver = static_cast<int>(Ybus.rows());

    Eigen::SparseMatrix<real_type> dcYbus = Eigen::SparseMatrix<real_type>(nb_bus_solver - 1, nb_bus_solver - 1);

    Eigen::SparseMatrix<cplx_type> dcYbus_tmp = Ybus;
    dcYbus_tmp.makeCompressed();
    const CplxVect & Sbus_tmp = Sbus;

    // find the slack bus
    int slack_bus_id_solver = extract_slack_bus_id(pv, pq, nb_bus_solver);

    // remove the slack bus from Ybus
    // and extract only real part
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(dcYbus_tmp.nonZeros());
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(dcYbus_tmp, k); it; ++it)
        {
            int row_res = static_cast<int>(it.row());
            if(row_res == slack_bus_id_solver) continue;
            row_res = row_res > slack_bus_id_solver ? row_res - 1 : row_res;
            int col_res = static_cast<int>(it.col());
            col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
            tripletList.push_back(Eigen::Triplet<real_type> (row_res, col_res, std::real(it.value())));
        }
    }
    dcYbus.setFromTriplets(tripletList.begin(), tripletList.end());
    dcYbus.makeCompressed();

    // initialize the solver
    Eigen::SparseLU<Eigen::SparseMatrix<real_type>, Eigen::COLAMDOrdering<int> >  dc_solver;
    dc_solver.analyzePattern(dcYbus);
    dc_solver.factorize(dcYbus);
    if(dc_solver.info() != Eigen::Success) {
        // matrix is not connected
        timer_total_nr_ += timer.duration();
        err_ = 1;
        return false;
    }

    // remove the slack bus from Sbus
    RealVect dcSbus = RealVect::Constant(nb_bus_solver - 1, my_zero_);
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        int col_res = k;
        col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
        dcSbus(col_res) = std::real(Sbus_tmp(k));
    }

    // solve for theta: Sbus = dcY . theta
    RealVect Va_dc_without_slack = dc_solver.solve(dcSbus);
    if(dc_solver.info() != Eigen::Success) {
        // solving failed, this should not happen in dc ...
        // matrix is not connected
        timer_total_nr_ += timer.duration();
        err_ = 3;
        return false;
    }

    // retrieve back the results in the proper shape (add back the slack bus)
    // TODO have a better way for this, for example using `.segment(0,npv)`
    // see the BaseSolver.cpp: _evaluate_Fx
    RealVect Va_dc = RealVect::Constant(nb_bus_solver, my_zero_);
    // fill Va from dc approx
    for (int bus_id_with_slack=0; bus_id_with_slack < nb_bus_solver; ++bus_id_with_slack){
        if(bus_id_with_slack == slack_bus_id_solver) continue;  // slack bus is handled elsewhere
        int bus_id_without_slack = bus_id_with_slack > slack_bus_id_solver ? bus_id_with_slack - 1 : bus_id_with_slack;
        Va_dc(bus_id_with_slack) = Va_dc_without_slack(bus_id_without_slack);
    }
    Va_dc.array() += std::arg(V(slack_bus_id_solver));

    // save the results
    Va_ = Va_dc;
    // add the Voltage setpoints of the generator
    Vm_ = V.array().abs(); // RealVect::Constant(V.size(), my_one_);
    //    Vm_(pv) = V(pv).array().abs();
    Vm_(slack_bus_id_solver) = std::abs(V(slack_bus_id_solver));

    // now compute the resulting complex voltage
    V_ = (Va_.array().cos().cast<cplx_type>() + my_i * Va_.array().sin().cast<cplx_type>());
    V_.array() *= Vm_.array();
    nr_iter_ = 1;
    V = V_;
    timer_total_nr_ += timer.duration();
    return true;
}


