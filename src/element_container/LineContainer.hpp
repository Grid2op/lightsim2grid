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

#include "Utils.hpp"
#include "BaseSubstation.hpp"
#include "OneSideContainer.hpp"
#include "TwoSidesContainer_rxh_A.hpp"

/**
This class is a container for all the powerlines on the grid.

**/
class LineContainer : public TwoSidesContainer_rxh_A<OneSideContainer>
{
    public:
        class LineInfo : public TwoSidesContainer_rxh_A<OneSideContainer>::TwoSidesContainer_rxh_AInfo
        {
            public:
                LineInfo(const LineContainer & r_data, int my_id):
                TwoSidesContainer_rxh_A<OneSideContainer>::TwoSidesContainer_rxh_AInfo(r_data, my_id) {}
        };
        typedef LineInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<LineContainer> LineConstIterator;
    
    public:
        // make it iterable
        typedef LineConstIterator const_iterator_type;
        LineConstIterator begin() const {return LineConstIterator(this, 0); }
        LineConstIterator end() const {return LineConstIterator(this, nb()); }
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

    public:
        typedef std::tuple<
                   TwoSidesContainer_rxh_A<OneSideContainer>::StateRes
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
        StateRes get_state() const
        {
            StateRes res(get_tsc_rxha_state());
            return res;
        }
        void set_state(LineContainer::StateRes & my_state )
        {
            set_tsc_rxha_state(std::get<0>(my_state));
            _update_model_coeffs();
            reset_results();
        }
        // template<class T>
        // void check_size(const T& my_state)
        // {
        //     //currently unused
        //     unsigned int size_th = 6;
        //     if (my_state.size() != size_th)
        //     {
        //         std::cout << "LightSim::LineContainer state size " << my_state.size() << " instead of "<< size_th << std::endl;
        //         // TODO more explicit error message
        //         throw std::runtime_error("Invalid state when loading LightSim::LineContainer");
        //     }
        // }
    
        // int nb() const { return static_cast<int>(powerlines_r_.size()); }
    
        // virtual void reconnect_connected_buses(Substation & substation) const;
        // virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component);
        // virtual void nb_line_end(std::vector<int> & res) const;
        // virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const;
        // virtual void update_bus_status(Substation & substation) const {
        //     const int nb_ = nb();
        //     for(int el_id = 0; el_id < nb_; ++el_id)
        //     {
        //         if(!status_[el_id]) continue;
        //         substation.reconnect_bus(bus_or_id_[el_id]);
        //         substation.reconnect_bus(bus_ex_id_[el_id]);
        //     }
        // }    
    
        // void update_topo(
        //     Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
        //     Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
        //     SolverControl & solver_control,
        //     Substation & substations
        // ) {}
        // void deactivate(int powerline_id, SolverControl & solver_control) {
        //     // std::cout << "line: deactivate called\n";
        //     if(status_[powerline_id]){
        //         solver_control.tell_recompute_ybus();
        //         // but sparsity pattern do not change here (possibly one more coeff at 0.)
        //         solver_control.tell_ybus_some_coeffs_zero();
        //         solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        //     }
        //     _generic_deactivate(powerline_id, status_);
        // }
        // void reactivate(int powerline_id, SolverControl & solver_control) {
        //     if(!status_[powerline_id]){
        //         solver_control.tell_recompute_ybus();
        //         solver_control.tell_ybus_change_sparsity_pattern();  // sparsity pattern might change: a non zero coeff can pop up
        //         solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
        //     }
        //     _generic_reactivate(powerline_id, status_);
        // }
        // void change_bus_or(int powerline_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        //     // std::cout << "line: change_bus_or called\n";
        //     _generic_change_bus(powerline_id, new_bus_id, bus_or_id_, solver_control, nb_bus);
        //     }
        // void change_bus_ex(int powerline_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        //     // std::cout << "line: change_bus_or called\n";
        //     _generic_change_bus(powerline_id, new_bus_id, bus_ex_id_, solver_control, nb_bus);
        //     }
        // int get_bus_or(int powerline_id) {return _get_bus(powerline_id, status_, bus_or_id_);}
        // int get_bus_ex(int powerline_id) {return _get_bus(powerline_id, status_, bus_ex_id_);}
        // virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
        //                       bool ac,
        //                       const std::vector<int> & id_grid_to_solver,
        //                       real_type sn_mva
        //                       ) const;
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
                                     
        // virtual void fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver);
        
        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const std::vector<int> & id_grid_to_solver,
                             const RealVect & bus_vn_kv,
                             real_type sn_mva,
                             bool ac){
            compute_results_tsc_rxha(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        }

        void reset_results() {reset_results_tsc_rxha();}
        
        // tuple4d get_lineor_res() const {return tuple4d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_);}
        // tuple4d get_lineex_res() const {return tuple4d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_);}
        // tuple5d get_res_or_full() const {return tuple5d(res_powerline_por_, res_powerline_qor_, res_powerline_vor_, res_powerline_aor_, res_powerline_thetaor_);}
        // tuple5d get_res_ex_full() const {return tuple5d(res_powerline_pex_, res_powerline_qex_, res_powerline_vex_, res_powerline_aex_, res_powerline_thetaex_);}
        // Eigen::Ref<const RealVect> get_theta_or() const {return res_powerline_thetaor_;}
        // Eigen::Ref<const RealVect> get_theta_ex() const {return res_powerline_thetaex_;}
        // const std::vector<bool>& get_status() const {return status_;}
        // Eigen::Ref<const Eigen::VectorXi> get_bus_from() const {return bus_or_id_;}
        // Eigen::Ref<const Eigen::VectorXi> get_bus_to() const {return bus_ex_id_;}
        
        // model paramters
        // Eigen::Ref<const CplxVect> yac_ff() const {return yac_ff_;}
        // Eigen::Ref<const CplxVect> yac_ft() const {return yac_ft_;}
        // Eigen::Ref<const CplxVect> yac_tf() const {return yac_tf_;}
        // Eigen::Ref<const CplxVect> yac_tt() const {return yac_tt_;}
        
        // Eigen::Ref<const CplxVect> ydc_ff() const {return ydc_ff_;}
        // Eigen::Ref<const CplxVect> ydc_ft() const {return ydc_ft_;}
        // Eigen::Ref<const CplxVect> ydc_tf() const {return ydc_tf_;}
        // Eigen::Ref<const CplxVect> ydc_tt() const {return ydc_tt_;}
        // for consistency with trafo, when used for example in BaseMultiplePowerflow...
        Eigen::Ref<const RealVect> dc_x_tau_shift() const {return RealVect();}
        
        // void set_or_pos_topo_vect(Eigen::Ref<const IntVect> pos_topo_vect)
        // {
        //     or_pos_topo_vect_.array() = pos_topo_vect;
        // }
        // void set_ex_pos_topo_vect(Eigen::Ref<const IntVect> pos_topo_vect)
        // {
        //     ex_pos_topo_vect_.array() = pos_topo_vect;
        // }
        
        // void set_or_subid(Eigen::Ref<const IntVect> subid)
        // {
        //     or_to_subid_.array() = subid;
        // }
        // void set_ex_subid(Eigen::Ref<const IntVect> subid)
        // {
        //     ex_to_subid_.array() = subid;
        // }
    
    protected:
        void _update_model_coeffs();

    protected:
        // physical properties
        // RealVect powerlines_r_;
        // RealVect powerlines_x_;
        // CplxVect powerlines_h_or_;
        // CplxVect powerlines_h_ex_;

        // specific grid2op
        // IntVect or_pos_topo_vect_;
        // IntVect ex_pos_topo_vect_;
        // IntVect or_to_subid_;
        // IntVect ex_to_subid_;

        // input data
        // Eigen::VectorXi bus_or_id_;
        // Eigen::VectorXi bus_ex_id_;
        // std::vector<bool> status_;

        //output data
        // RealVect res_powerline_por_;  // in MW
        // RealVect res_powerline_qor_;  // in MVar
        // RealVect res_powerline_vor_;  // in kV
        // RealVect res_powerline_aor_;  // in kA
        // RealVect res_powerline_pex_;  // in MW
        // RealVect res_powerline_qex_;  // in MVar
        // RealVect res_powerline_vex_;  // in kV
        // RealVect res_powerline_aex_;  // in kA
        // RealVect res_powerline_thetaor_; // in degree
        // RealVect res_powerline_thetaex_; // in degree

        // model coefficients
        // CplxVect yac_ff_;
        // CplxVect yac_ft_;
        // CplxVect yac_tf_;
        // CplxVect yac_tt_;

        // CplxVect ydc_ff_;
        // CplxVect ydc_ft_;
        // CplxVect ydc_tf_;
        // CplxVect ydc_tt_;
};

#endif  //LINE_CONTAINER_H
