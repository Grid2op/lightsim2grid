// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LineContainer.h"

#include <sstream>

void LineContainer::init(const RealVect & branch_r,
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

    // TODO check what can be checked: branch_* have same size, no id in branch_to_id that are
    // TODO not in [0, .., buv_vn_kv.size()] etc.

    //TODO consistency with trafo: have a converter methods to convert this value into pu, and store the pu
    // in this method

    int size = static_cast<int>(branch_r.size());
    GenericContainer::check_size(branch_r, size, "branch_r");
    GenericContainer::check_size(branch_x, size, "branch_x");
    GenericContainer::check_size(branch_h, size, "branch_h");
    GenericContainer::check_size(branch_from_id, size, "branch_from_id");
    GenericContainer::check_size(branch_to_id, size, "branch_to_id");

    bus_or_id_ = branch_from_id;
    bus_ex_id_ = branch_to_id;
    powerlines_h_or_ = 0.5 * branch_h;
    powerlines_h_ex_ = 0.5 * branch_h;
    powerlines_r_ = branch_r;
    powerlines_x_ = branch_x;
    status_ = std::vector<bool>(branch_r.size(), true); // by default everything is connected
    _update_model_coeffs();
    reset_results();
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

    bus_or_id_ = branch_from_id;
    bus_ex_id_ = branch_to_id;
    powerlines_h_or_ = branch_h_or;
    powerlines_h_ex_ = branch_h_ex;
    powerlines_r_ = branch_r;
    powerlines_x_ = branch_x;
    status_ = std::vector<bool>(branch_r.size(), true); // by default everything is connected
    _update_model_coeffs();
    reset_results();
}

LineContainer::StateRes LineContainer::get_state() const
{
     std::vector<real_type> branch_r(powerlines_r_.begin(), powerlines_r_.end());
     std::vector<real_type> branch_x(powerlines_x_.begin(), powerlines_x_.end());
     std::vector<cplx_type > branch_hor(powerlines_h_or_.begin(), powerlines_h_or_.end());
     std::vector<cplx_type > branch_hex(powerlines_h_ex_.begin(), powerlines_h_ex_.end());
     std::vector<int > branch_from_id(bus_or_id_.begin(), bus_or_id_.end());
     std::vector<int > branch_to_id(bus_ex_id_.begin(), bus_ex_id_.end());
     std::vector<bool> status = status_;
     LineContainer::StateRes res(names_, branch_r, branch_x, branch_hor, branch_hex, branch_from_id, branch_to_id, status);
     return res;
}
void LineContainer::set_state(LineContainer::StateRes & my_state)
{
    names_ = std::get<0>(my_state);
    std::vector<real_type> & branch_r = std::get<1>(my_state);
    std::vector<real_type> & branch_x = std::get<2>(my_state);
    std::vector<cplx_type > & branch_h_or = std::get<3>(my_state);
    std::vector<cplx_type > & branch_h_ex = std::get<4>(my_state);
    std::vector<int> & branch_from_id = std::get<5>(my_state);
    std::vector<int> & branch_to_id = std::get<6>(my_state);
    std::vector<bool> & status = std::get<7>(my_state);
    // TODO check sizes

    // now assign the values
    powerlines_r_ = RealVect::Map(&branch_r[0], branch_r.size());
    powerlines_x_ = RealVect::Map(&branch_x[0], branch_x.size());
    powerlines_h_or_ = CplxVect::Map(&branch_h_or[0], branch_h_or.size());
    powerlines_h_ex_ = CplxVect::Map(&branch_h_ex[0], branch_h_ex.size());

    // input data
    bus_or_id_ = Eigen::VectorXi::Map(&branch_from_id[0], branch_from_id.size());
    bus_ex_id_ = Eigen::VectorXi::Map(&branch_to_id[0], branch_to_id.size());
    status_ = status;

    _update_model_coeffs();
    reset_results();
}

void LineContainer::_update_model_coeffs()
{
    const auto my_size = powerlines_r_.size();

    yac_ff_ = CplxVect::Zero(my_size);
    yac_ft_ = CplxVect::Zero(my_size);
    yac_tf_ = CplxVect::Zero(my_size);
    yac_tt_ = CplxVect::Zero(my_size);

    ydc_ff_ = CplxVect::Zero(my_size);
    ydc_ft_ = CplxVect::Zero(my_size);
    ydc_tf_ = CplxVect::Zero(my_size);
    ydc_tt_ = CplxVect::Zero(my_size);
    for(int i = 0; i < my_size; ++i)
    {
        // for AC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.2
        const cplx_type ys = 1. / (powerlines_r_(i) + my_i * powerlines_x_(i));
        const cplx_type h_or = my_i * powerlines_h_or_(i);
        const cplx_type h_ex = my_i * powerlines_h_ex_(i);
        yac_ff_(i) = (ys + h_or);
        yac_tt_(i) = (ys + h_ex);
        yac_tf_(i) = -ys;
        yac_ft_(i) = -ys;

        // for DC
        // see https://matpower.org/docs/MATPOWER-manual.pdf eq. 3.21
        // except here I only care about the real part, so I remove the "1/j"
        cplx_type tmp = 1. / (powerlines_x_(i));
        ydc_ff_(i) = tmp;
        ydc_tt_(i) = tmp;
        ydc_tf_(i) = -tmp;
        ydc_ft_(i) = -tmp;
    }
}

void LineContainer::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver)
{
    throw std::runtime_error("You should not use that!");
}

void LineContainer::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                             bool ac,
                             const std::vector<int> & id_grid_to_solver,
                             real_type sn_mva) const
{
    // fill the matrix
    //TODO template here instead of "if" for ac / dc
    const Eigen::Index nb_line = static_cast<int>(powerlines_r_.size());
    cplx_type yft, ytf, yff, ytt;

    //diagonal coefficients
    for(Eigen::Index line_id =0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        // compute from / to
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "Line::fillYbusBranch: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (or side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "Line::fillYbusBranch: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (ex side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        
        if(ac){
            // ac mode
            yft = yac_ft_(line_id);
            ytf = yac_tf_(line_id);
            yff = yac_ff_(line_id);
            ytt = yac_tt_(line_id);
        }else{
            // dc mode
            yft = ydc_ft_(line_id);
            ytf = ydc_tf_(line_id);
            yff = ydc_ff_(line_id);
            ytt = ydc_tt_(line_id);
        }
        res.push_back(Eigen::Triplet<cplx_type> (bus_or_solver_id, bus_ex_solver_id, yft));
        res.push_back(Eigen::Triplet<cplx_type> (bus_ex_solver_id, bus_or_solver_id, ytf));
        res.push_back(Eigen::Triplet<cplx_type> (bus_or_solver_id, bus_or_solver_id, yff));
        res.push_back(Eigen::Triplet<cplx_type> (bus_ex_solver_id, bus_ex_solver_id, ytt));
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
    const Eigen::Index nb_line = static_cast<int>(powerlines_r_.size());
    real_type yft_bp, ytf_bp, yff_bp, ytt_bp;
    real_type yft_bpp, ytf_bpp, yff_bpp, ytt_bpp;
    //diagonal coefficients
    for(Eigen::Index line_id=0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBp_Bpp: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (or side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
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
            ys_bp = 1. / (0. + my_i * powerlines_x_(line_id));
            ys_bpp = 1. / (powerlines_r_(line_id) + my_i * powerlines_x_(line_id));
        }else if (xb_or_bx==FDPFMethod::BX){
            ys_bp = 1. / (powerlines_r_(line_id) + my_i * powerlines_x_(line_id));
            ys_bpp = 1. / (0. + my_i * powerlines_x_(line_id));
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
        yff_bpp = ys_bpp_r + std::imag(my_i * powerlines_h_or_(line_id));
        ytt_bpp = ys_bpp_r + std::imag(my_i * powerlines_h_ex_(line_id));
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


void LineContainer::fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & Bf,
                                    const std::vector<int> & id_grid_to_solver,
                                    real_type sn_mva,
                                    int nb_powerline,
                                    bool transpose) const
{
    const Eigen::Index nb_line = powerlines_r_.size();

    for(Eigen::Index line_id=0; line_id < nb_line; ++line_id){
        // i only add this if the powerline is connected
        if(!status_[line_id]) continue;

        // get the from / to bus id
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBf_for_PTDF: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (or side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::fillBf_for_PTDF: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (ex side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        real_type x = powerlines_x_(line_id);
        
        // TODO
        // Bf (nb_branch, nb_bus) : en dc un truc du genre 1 / x / tap for (1..nb_branch, from_bus)
        // and -1. / x / tap for (1..nb_branch, to_bus) 
        if(transpose){
            Bf.push_back(Eigen::Triplet<real_type> (bus_or_solver_id, line_id, 1. / x));
            Bf.push_back(Eigen::Triplet<real_type> (bus_ex_solver_id, line_id, -1. / x));
        }else{
            Bf.push_back(Eigen::Triplet<real_type> (line_id, bus_or_solver_id, 1. / x));
            Bf.push_back(Eigen::Triplet<real_type> (line_id, bus_ex_solver_id, -1. / x));
        }
    }

}


void LineContainer::reset_results()
{
    res_powerline_por_ = RealVect(nb());  // in MW
    res_powerline_qor_ = RealVect(nb());  // in MVar
    res_powerline_vor_ = RealVect(nb());  // in kV
    res_powerline_aor_ = RealVect(nb());  // in kA
    res_powerline_pex_ = RealVect(nb());  // in MW
    res_powerline_qex_ = RealVect(nb());  // in MVar
    res_powerline_vex_ = RealVect(nb());  // in kV
    res_powerline_aex_ = RealVect(nb());  // in kA
    res_powerline_thetaor_ = RealVect(nb());
    res_powerline_thetaex_ = RealVect(nb());
}

void LineContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                    const Eigen::Ref<const RealVect> & Vm,
                                    const Eigen::Ref<const CplxVect> & V,
                                    const std::vector<int> & id_grid_to_solver,
                                    const RealVect & bus_vn_kv,
                                    real_type sn_mva,
                                    bool ac)
{
    // it needs to be initialized at 0.
    Eigen::Index nb_element = nb();
    for(Eigen::Index line_id = 0; line_id < nb_element; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) {
            res_powerline_por_(line_id) = 0.0;  // in MW
            res_powerline_qor_(line_id) = 0.0;  // in MVar
            res_powerline_vor_(line_id) = 0.0;  // in kV
            res_powerline_aor_(line_id) = 0.0;  // in kA
            res_powerline_pex_(line_id) = 0.0;  // in MW
            res_powerline_qex_(line_id) = 0.0;  // in MVar
            res_powerline_vex_(line_id) = 0.0;  // in kV
            res_powerline_aex_(line_id) = 0.0;  // in kA
            res_powerline_thetaor_(line_id) = 0.0;  // in kV
            res_powerline_thetaex_(line_id) = 0.0;  // in kV
            continue;
        }

        // connectivity
        int bus_or_id_me = bus_or_id_(line_id);
        int bus_or_solver_id = id_grid_to_solver[bus_or_id_me];
        if(bus_or_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::compute_results: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (or side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        int bus_ex_id_me = bus_ex_id_(line_id);
        int bus_ex_solver_id = id_grid_to_solver[bus_ex_id_me];
        if(bus_ex_solver_id == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LineContainer::compute_results: the line with id ";
            exc_ << line_id;
            exc_ << " is connected (ex side) to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }

        // retrieve voltages magnitude in kv instead of pu
        real_type v_or = Vm(bus_or_solver_id);
        real_type v_ex = Vm(bus_ex_solver_id);
        real_type bus_vn_kv_or = bus_vn_kv(bus_or_id_me);
        real_type bus_vn_kv_ex = bus_vn_kv(bus_ex_id_me);
        
        // for the voltage
        res_powerline_vor_(line_id) = v_or * bus_vn_kv_or;
        res_powerline_vex_(line_id) = v_ex * bus_vn_kv_ex;

        // retrieve the voltage angle in degree (instead of radian)
        res_powerline_thetaor_(line_id) = Va(bus_or_solver_id) * my_180_pi_;
        res_powerline_thetaex_(line_id) = Va(bus_ex_solver_id) * my_180_pi_;

        // results of the powerflow
        cplx_type Eor = V(bus_or_solver_id);
        cplx_type Eex = V(bus_ex_solver_id);

        if(ac){
            // result of the ac powerflow
            cplx_type I_orex =  yac_ff_(line_id) * Eor + yac_ft_(line_id) * Eex;
            cplx_type I_exor =  yac_tt_(line_id) * Eex + yac_tf_(line_id) * Eor;

            I_orex = std::conj(I_orex);
            I_exor = std::conj(I_exor);
            cplx_type s_orex = Eor * I_orex;
            cplx_type s_exor = Eex * I_exor;

            res_powerline_por_(line_id) = std::real(s_orex) * sn_mva;
            res_powerline_qor_(line_id) = std::imag(s_orex) * sn_mva;
            res_powerline_pex_(line_id) = std::real(s_exor) * sn_mva;
            res_powerline_qex_(line_id) = std::imag(s_exor) * sn_mva;

        }else{
            // result of the dc powerflow
            res_powerline_por_(line_id) = (std::real(ydc_ff_(line_id)) * Va(bus_or_solver_id) + std::real(ydc_ft_(line_id)) * Va(bus_ex_solver_id)) * sn_mva;
            res_powerline_pex_(line_id) = (std::real(ydc_tt_(line_id)) * Va(bus_ex_solver_id) + std::real(ydc_tf_(line_id)) * Va(bus_or_solver_id)) * sn_mva;   

            res_powerline_qor_(line_id) = 0.;
            res_powerline_qex_(line_id) = 0.;
            // for the voltage (by hypothesis vm = 1)
            // res_powerline_vor_(line_id) = bus_vn_kv_or;
            // res_powerline_vex_(line_id) = bus_vn_kv_ex;        
        }
    }
    _get_amps(res_powerline_aor_, res_powerline_por_, res_powerline_qor_, res_powerline_vor_);
    _get_amps(res_powerline_aex_, res_powerline_pex_, res_powerline_qex_, res_powerline_vex_);
}

void LineContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const{

    const auto my_size = powerlines_r_.size();
    for(Eigen::Index line_id = 0; line_id < my_size; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;
        
        const auto bus_or_id_me = bus_or_id_(line_id);        
        if(bus_or_id_me == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "LineContainer::reconnect_connected_buses: Line with id ";
            exc_ << line_id;
            exc_ << " is connected (origin) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_powerline(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[bus_or_id_me] = true;

        const auto bus_ex_id_me = bus_ex_id_(line_id);        
        if(bus_ex_id_me == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "LineContainer::reconnect_connected_buses: Line with id ";
            exc_ << line_id;
            exc_ << " is connected (ext) to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_powerline(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[bus_ex_id_me] = true;
    }
}

void LineContainer::nb_line_end(std::vector<int> & res) const{
    const auto my_size = powerlines_r_.size();
    for(Eigen::Index line_id = 0; line_id < my_size; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;
        const auto bus_or = bus_or_id_(line_id);
        const auto bus_ex = bus_ex_id_(line_id);
        res[bus_or] += 1;
        res[bus_ex] += 1;
    }
}

void LineContainer::get_graph(std::vector<Eigen::Triplet<real_type> > & res) const
{
    const auto my_size = powerlines_r_.size();
    for(Eigen::Index line_id = 0; line_id < my_size; ++line_id){
        // don't do anything if the element is disconnected
        if(!status_[line_id]) continue;
        const auto bus_or = bus_or_id_(line_id);
        const auto bus_ex = bus_ex_id_(line_id);
        res.push_back(Eigen::Triplet<real_type>(bus_or, bus_ex, 1.));
        res.push_back(Eigen::Triplet<real_type>(bus_ex, bus_or, 1.));
    }
}

void LineContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component)
{
    const Eigen::Index nb_line = nb();
    SolverControl unused_solver_control;
    for(Eigen::Index i = 0; i < nb_line; ++i){
        if(!status_[i]) continue;
        auto bus_or = bus_or_id_(i);
        auto bus_ex = bus_ex_id_(i);
        if(!busbar_in_main_component[bus_or] || !busbar_in_main_component[bus_ex]){
            deactivate(i, unused_solver_control);
        }
    }
}
