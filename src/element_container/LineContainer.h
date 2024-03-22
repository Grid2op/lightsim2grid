// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LINE_CONTAINER_H
#define LINE_CONTAINER_H

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "GenericContainer.h"

/**
This class is a container for all the powerlines on the grid.

**/
class LineContainer : public GenericContainer
{
    public:
        class LineInfo
        {
            public:
                // members
                int id;  // id of the line
                std::string name;
                bool connected;
                int bus_or_id;
                int bus_ex_id;
                real_type r_pu;
                real_type x_pu;
                cplx_type h_pu;
                cplx_type h_or_pu;
                cplx_type h_ex_pu;

                bool has_res;
                real_type res_p_or_mw;
                real_type res_q_or_mvar;
                real_type res_v_or_kv;
                real_type res_a_or_ka;
                real_type res_theta_or_deg;
                real_type res_p_ex_mw;
                real_type res_q_ex_mvar;
                real_type res_v_ex_kv;
                real_type res_a_ex_ka;
                real_type res_theta_ex_deg;

                LineInfo(const LineContainer & r_data_line, int my_id):
                id(my_id),
                name(""),
                connected(false),
                bus_or_id(-1),
                bus_ex_id(-1),
                r_pu(-1.0),
                x_pu(-1.0),
                h_pu(0., 0.),
                h_or_pu(0., 0.),
                h_ex_pu(0., 0.),
                has_res(false),
                res_p_or_mw(0.),
                res_q_or_mvar(0.),
                res_v_or_kv(0.),
                res_a_or_ka(0.),
                res_theta_or_deg(0.),
                res_p_ex_mw(0.),
                res_q_ex_mvar(0.),
                res_v_ex_kv(0.),
                res_a_ex_ka(0.),
                res_theta_ex_deg(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_line.nb()))
                    {
                        id = my_id;
                        if(r_data_line.names_.size()){
                            name = r_data_line.names_[my_id];
                        }
                        connected = r_data_line.status_[my_id];
                        bus_or_id = r_data_line.bus_or_id_.coeff(my_id);
                        bus_ex_id = r_data_line.bus_ex_id_.coeff(my_id);
                        r_pu = r_data_line.powerlines_r_.coeff(my_id);
                        x_pu = r_data_line.powerlines_x_.coeff(my_id);
                        h_or_pu = r_data_line.powerlines_h_or_.coeff(my_id);
                        h_ex_pu = r_data_line.powerlines_h_ex_.coeff(my_id);
                        h_pu = (h_or_pu + h_ex_pu);

                        has_res = r_data_line.res_powerline_por_.size() > 0;
                        if(has_res)
                        {
                            res_p_or_mw = r_data_line.res_powerline_por_.coeff(my_id);
                            res_q_or_mvar = r_data_line.res_powerline_qor_.coeff(my_id);
                            res_v_or_kv = r_data_line.res_powerline_vor_.coeff(my_id);
                            res_a_or_ka = r_data_line.res_powerline_aor_.coeff(my_id);
                            res_p_ex_mw = r_data_line.res_powerline_pex_.coeff(my_id);
                            res_q_ex_mvar = r_data_line.res_powerline_qex_.coeff(my_id);
                            res_v_ex_kv = r_data_line.res_powerline_vex_.coeff(my_id);
                            res_a_ex_ka = r_data_line.res_powerline_aex_.coeff(my_id);
                            res_theta_or_deg = r_data_line.res_powerline_thetaor_.coeff(my_id);
                            res_theta_ex_deg = r_data_line.res_powerline_thetaex_.coeff(my_id);
                        }
                    }
                }
        };
        typedef LineInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<LineContainer> LineConstIterator;

    public:
    typedef std::tuple<
               std::vector<std::string>,
               std::vector<real_type>, // branch_r
               std::vector<real_type>, // branch_x
               std::vector<cplx_type>, // branch_h
               std::vector<cplx_type>, // branch_h
               std::vector<int>, // branch_from_id
               std::vector<int>, // branch_to_id
               std::vector<bool> // status_
               >  StateRes;

    LineContainer() {};

    void init(const RealVect & branch_r,
              const RealVect & branch_x,
              const CplxVect & branch_h,
              const Eigen::VectorXi & branch_from_id,
              const Eigen::VectorXi & branch_to_id
              );

    void init(const RealVect & branch_r,
              const RealVect & branch_x,
              const CplxVect & branch_h_or,
              const CplxVect & branch_h_ex,
              const Eigen::VectorXi & branch_from_id,
              const Eigen::VectorXi & branch_to_id
              );

    // pickle
    LineContainer::StateRes get_state() const;
    void set_state(LineContainer::StateRes & my_state );
    template<class T>
    void check_size(const T& my_state)
    {
        //currently unused
        unsigned int size_th = 6;
        if (my_state.size() != size_th)
        {
            std::cout << "LightSim::LineContainer state size " << my_state.size() << " instead of "<< size_th << std::endl;
            // TODO more explicit error message
            throw std::runtime_error("Invalid state when loading LightSim::LineContainer");
        }
    }

    int nb() const { return static_cast<int>(powerlines_r_.size()); }

    // make it iterable
    typedef LineConstIterator const_iterator_type;
    const_iterator_type begin() const {return LineConstIterator(this, 0); }
    const_iterator_type end() const {return LineConstIterator(this, nb()); }
    LineInfo operator[](int id) const
    {
        if(id < 0)
        {
            throw std::range_error("You cannot ask for line with negative id");
        }
        if(id >= nb())
        {
            throw std::range_error("Generator out of bound. Not enough powerlines on the grid.");
        }
        return LineInfo(*this, id);
    }
    virtual void reconnect_connected_buses(std::vector<bool> & bus_status) const;
    virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component);
    virtual void nb_line_end(std::vector<int> & res) const;
    virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const;
    virtual void update_bus_status(std::vector<bool> & bus_status) const {
        const int nb_ = nb();
        for(int el_id = 0; el_id < nb_; ++el_id)
        {
            if(!status_[el_id]) continue;
            bus_status[bus_or_id_[el_id]] = true;
            bus_status[bus_ex_id_[el_id]] = true;
        }
    }    

    void deactivate(int powerline_id, SolverControl & solver_control) {
        // std::cout << "line: deactivate called\n";
        if(status_[powerline_id]){
            solver_control.tell_recompute_ybus();
            // but sparsity pattern do not change here (possibly one more coeff at 0.)
            solver_control.tell_ybus_some_coeffs_zero();
            solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        }
        _deactivate(powerline_id, status_);
    }
    void reactivate(int powerline_id, SolverControl & solver_control) {
        if(!status_[powerline_id]){
            solver_control.tell_recompute_ybus();
            solver_control.tell_ybus_change_sparsity_pattern();  // sparsity pattern might change: a non zero coeff can pop up
            solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        }
        _reactivate(powerline_id, status_);
    }
    void change_bus_or(int powerline_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        // std::cout << "line: change_bus_or called\n";
        _change_bus(powerline_id, new_bus_id, bus_or_id_, solver_control, nb_bus);
        }
    void change_bus_ex(int powerline_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        // std::cout << "line: change_bus_or called\n";
        _change_bus(powerline_id, new_bus_id, bus_ex_id_, solver_control, nb_bus);
        }
    int get_bus_or(int powerline_id) {return _get_bus(powerline_id, status_, bus_or_id_);}
    int get_bus_ex(int powerline_id) {return _get_bus(powerline_id, status_, bus_ex_id_);}
    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                          bool ac,
                          const std::vector<int> & id_grid_to_solver,
                          real_type sn_mva
                          ) const;
    virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                            std::vector<Eigen::Triplet<real_type> > & Bpp,
                            const std::vector<int> & id_grid_to_solver,
                            real_type sn_mva,
                            FDPFMethod xb_or_bx) const;
    virtual void fillBf_for_PTDF(std::vector<Eigen::Triplet<real_type> > & Bf,
                                 const std::vector<int> & id_grid_to_solver,
                                 real_type sn_mva,
                                 int nb_powerline,
                                 bool transpose) const;
                                 
    virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);

    void compute_results(const Eigen::Ref<const RealVect> & Va,
                         const Eigen::Ref<const RealVect> & Vm,
                         const Eigen::Ref<const CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva,
                         bool ac);
    void reset_results();

    tuple4d get_lineor_res() const {return tuple4d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_);}
    tuple4d get_lineex_res() const {return tuple4d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_);}
    tuple5d get_res_or_full() const {return tuple5d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_, res_powerline_thetaor_);}
    tuple5d get_res_ex_full() const {return tuple5d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_, res_powerline_thetaex_);}
    Eigen::Ref<const RealVect> get_theta_or() const {return res_powerline_thetaor_;}
    Eigen::Ref<const RealVect> get_theta_ex() const {return res_powerline_thetaex_;}
    const std::vector<bool>& get_status() const {return status_;}
    Eigen::Ref<const Eigen::VectorXi> get_bus_from() const {return bus_or_id_;}
    Eigen::Ref<const Eigen::VectorXi> get_bus_to() const {return bus_ex_id_;}

    // model paramters
    Eigen::Ref<const CplxVect> yac_ff() const {return yac_ff_;}
    Eigen::Ref<const CplxVect> yac_ft() const {return yac_ft_;}
    Eigen::Ref<const CplxVect> yac_tf() const {return yac_tf_;}
    Eigen::Ref<const CplxVect> yac_tt() const {return yac_tt_;}

    Eigen::Ref<const CplxVect> ydc_ff() const {return ydc_ff_;}
    Eigen::Ref<const CplxVect> ydc_ft() const {return ydc_ft_;}
    Eigen::Ref<const CplxVect> ydc_tf() const {return ydc_tf_;}
    Eigen::Ref<const CplxVect> ydc_tt() const {return ydc_tt_;}
    // for consistency with trafo, when used for example in BaseMultiplePowerflow...
    Eigen::Ref<const RealVect> dc_x_tau_shift() const {return RealVect();}

    protected:
        void _update_model_coeffs();

    protected:
        // physical properties
        RealVect powerlines_r_;
        RealVect powerlines_x_;
        CplxVect powerlines_h_or_;
        CplxVect powerlines_h_ex_;

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
        RealVect res_powerline_thetaor_; // in degree
        RealVect res_powerline_thetaex_; // in degree

        // model coefficients
        CplxVect yac_ff_;
        CplxVect yac_ft_;
        CplxVect yac_tf_;
        CplxVect yac_tt_;

        CplxVect ydc_ff_;
        CplxVect ydc_ft_;
        CplxVect ydc_tf_;
        CplxVect ydc_tt_;
};

#endif  //LINE_CONTAINER_H
