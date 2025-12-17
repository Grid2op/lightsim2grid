// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LineContainer.hpp"

#include <sstream>

LineInfo::LineInfo(const LineContainer & r_data, int my_id):
TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::TwoSidesContainer_rxh_AInfo(r_data, my_id) {}

void LineContainer::init(
    const RealVect & branch_r,
    const RealVect & branch_x,
    const CplxVect & branch_h,
    const Eigen::VectorXi & branch_from_id,
    const Eigen::VectorXi & branch_to_id
)
{
    /**
    This method initialize the Ybus matrix from the branch matrix.
    It has to be called once when the solver is initialized. Afterwards, a call to
    "updateYbus" should be made for performance optimiaztion instead. //TODO
    **/

   init(
    branch_r, 
    branch_x,
    0.5 * branch_h,
    0.5 * branch_h,
    branch_from_id,
    branch_to_id);
}

void LineContainer::init(const RealVect & branch_r,
                         const RealVect & branch_x,
                         const CplxVect & branch_h_or,
                         const CplxVect & branch_h_ex,
                         const Eigen::VectorXi & branch_from_id,
                         const Eigen::VectorXi & branch_to_id
                         )
{
    /**
    This method initialize the Ybus matrix from the branch matrix.
    It has to be called once when the solver is initialized. Afterwards, a call to
    "updateYbus" should be made for performance optimiaztion instead. //TODO
    **/

    // TODO check what can be checked: branch_* have same size, no id in branch_to_id that are
    // TODO not in [0, .., buv_vn_kv.size()] etc.

    //TODO consistency with trafo: have a converter methods to convert this value into pu, and store the pu
    // in this method

    int size = static_cast<int>(branch_r.size());
    GenericContainer::check_size(branch_r, size, "branch_r");
    GenericContainer::check_size(branch_x, size, "branch_x");
    GenericContainer::check_size(branch_h_or, size, "branch_h_or");
    GenericContainer::check_size(branch_h_ex, size, "branch_h_ex");
    GenericContainer::check_size(branch_from_id, size, "branch_from_id");
    GenericContainer::check_size(branch_to_id, size, "branch_to_id");

    // bus_or_id_ = branch_from_id;
    // bus_ex_id_ = branch_to_id;
    // status_ = std::vector<bool>(branch_r.size(), true); // by default everything is connected
    h_side_1_ = branch_h_or;
    h_side_2_ = branch_h_ex;
    r_ = branch_r;
    x_ = branch_x;
    init_tsc(branch_from_id, branch_to_id, "trafo");
    _update_model_coeffs();
    reset_results();
}

void LineContainer::_update_model_coeffs()
{
    const auto my_size = r_.size();

    yac_11_ = CplxVect::Zero(my_size);
    yac_12_ = CplxVect::Zero(my_size);
    yac_21_ = CplxVect::Zero(my_size);
    yac_22_ = CplxVect::Zero(my_size);

    ydc_11_ = CplxVect::Zero(my_size);
    ydc_12_ = CplxVect::Zero(my_size);
    ydc_21_ = CplxVect::Zero(my_size);
    ydc_22_ = CplxVect::Zero(my_size);
    for(int i = 0; i < my_size; ++i)
    {
        // for AC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.2
        const cplx_type ys = 1. / cplx_type(r_(i), x_(i));
        const cplx_type h_or = h_side_1_(i);
        const cplx_type h_ex = h_side_2_(i);
        yac_11_(i) = (ys + h_or);
        yac_22_(i) = (ys + h_ex);
        yac_12_(i) = -ys;
        yac_21_(i) = -ys;

        // for DC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
        // except here I only care about the real part, so I remove the "1/j"
        cplx_type tmp = 1. / cplx_type(x_(i), 0.);
        ydc_11_(i) = tmp;
        ydc_22_(i) = tmp;
        ydc_21_(i) = -tmp;
        ydc_12_(i) = -tmp;
    }
}

void LineContainer::fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                               std::vector<Eigen::Triplet<real_type> > & Bpp,
                               const std::vector<int> & id_grid_to_solver,
                               real_type sn_mva,
                               FDPFMethod xb_or_bx) const
{

    // For Bp
    // temp_branch[:, BR_B] = zeros(nl)           ## zero out line charging shunts
    // temp_branch[:, TAP] = ones(nl)             ## cancel out taps
    // if alg == 2:                               ## if XB method
    //    temp_branch[:, BR_R] = zeros(nl)       ## zero out line resistance

    // For Bpp
    // temp_branch[:, SHIFT] = zeros(nl)          ## zero out phase shifters
    // if alg == 3:                               ## if BX method
    //     temp_branch[:, BR_R] = zeros(nl)    ## zero out line resistance
    const Eigen::Index nb_line = static_cast<int>(r_.size());
    real_type yft_bp, ytf_bp, yff_bp, ytt_bp;
    real_type yft_bpp, ytf_bpp, yff_bpp, ytt_bpp;
    //diagonal coefficients
    for(Eigen::Index line_id=0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_global_[line_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = get_bus_side_1(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBp_Bpp: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (or side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = get_bus_side_2(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBp_Bpp: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (ex side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }

        // get the coefficients
        cplx_type ys_bp, ys_bpp;
        if(xb_or_bx==FDPFMethod::XB){
            ys_bp = 1. / (0. + my_i * x_(line_id));
            ys_bpp = 1. / (r_(line_id) + my_i * x_(line_id));
        }else if (xb_or_bx==FDPFMethod::BX){
            ys_bp = 1. / (r_(line_id) + my_i * x_(line_id));
            ys_bpp = 1. / (0. + my_i * x_(line_id));
        }else{
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBp_Bpp: unknown method for the FDPF powerflow for line id ";
            exc_ << line_id;
            throw std::runtime_error(exc_.str());            
        }
        const real_type ys_bp_r = std::imag(ys_bp); 
        yff_bp = ys_bp_r;
        ytt_bp = ys_bp_r;
        yft_bp = -ys_bp_r;
        ytf_bp = -ys_bp_r;
        const real_type ys_bpp_r = std::imag(ys_bpp); 
        yff_bpp = ys_bpp_r + std::imag(h_side_1_(line_id));
        ytt_bpp = ys_bpp_r + std::imag(h_side_2_(line_id));
        yft_bpp = -ys_bpp_r;
        ytf_bpp = -ys_bpp_r;

        // and now add them
        Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_ex_solver_id, -yft_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_or_solver_id, -ytf_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_or_solver_id, -yff_bp));
        Bp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_ex_solver_id, -ytt_bp));

        Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_ex_solver_id, -yft_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_or_solver_id, -ytf_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, bus_or_solver_id, -yff_bpp));
        Bpp.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, bus_ex_solver_id, -ytt_bpp));

    }
}
