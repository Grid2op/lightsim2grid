// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef COMPUTERS_H
#define COMPUTERS_H

#include "GridModel.h"

class Computers
{
    public:
        typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RealMat;
        typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CplxMat;
        
        Computers(const GridModel & init_grid_model):
                _grid_model(init_grid_model),
                _Sbuses(),
                _amps_flows(),
                _voltages(),
                _nb_solved(0),
                _status(1), // 1: success, 0: failure
                _compute_flows(true),
                _solver(),
                _timer_total(0.) ,
                _timer_solver(0.) ,
                _timer_pre_proc(0.) 
                {
                    // make sure that my "grid_model" is ready to be used (for ac and dc)
                    Eigen::Index nb_bus = init_grid_model.total_bus();
                    CplxVect V = CplxVect::Constant(nb_bus, 1.);
                    const auto & Vtmp = init_grid_model.get_V();
                    for(int i = 0; i < Vtmp.size(); ++i) V[i] = Vtmp[i];
                    _grid_model.dc_pf(V, 10, 1e-5);
                    _grid_model.ac_pf(V, 10, 1e-5);
                };

        Computers(const Computers&) = delete;

        // control on whether I compute the flows or not
        void deactivate_flow_computations() {_compute_flows = false;}
        void activate_flow_computations() {_compute_flows = true;}

        // solver "control"
        void change_solver(const SolverType & type){
            _solver.change_solver(type);
        }
        std::vector<SolverType> available_solvers() {return _solver.available_solvers(); }
        SolverType get_solver_type() {return _solver.get_type(); }

        // timers
        double total_time() const {return _timer_total;}
        double solver_time() const {return _timer_solver;}
        double preprocessing_time() const {return _timer_pre_proc;}
        int nb_solved() const {return _nb_solved;}

        // status
        int get_status() const {return _status;}

        /**
        This function computes the results of running as many powerflow when varying the 
        injection (Sbus). 

        Each line of `Sbuses` will be a time step, and each column with 
        **/
        int compute_Sbuses(Eigen::Ref<const RealMat> gen_p,
                           Eigen::Ref<const RealMat> sgen_p,
                           Eigen::Ref<const RealMat> load_p,
                           Eigen::Ref<const RealMat> load_q,
                           const CplxVect & Vinit,
                           int max_iter,
                           real_type tol);

        Eigen::Ref<const RealMat > get_flows() const {return _amps_flows;}
        Eigen::Ref<const CplxMat > get_voltages() const {return _voltages;}
        Eigen::Ref<const CplxMat > get_sbuses() const {return _Sbuses;}

    protected:
        template<class T>
        void fill_SBus_real(CplxMat & Sbuses,
                            const T & structure_data,
                            const RealMat & temporal_data,
                            const std::vector<int> & id_me_to_ac_solver,
                            bool add  // if true call += else calls -=
                            ) const 
        {
            auto nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            int  bus_id_solver, bus_id_me;
            for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                bus_id_solver = id_me_to_ac_solver[bus_id_me];
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver) += tmp;
                else Sbuses.col(bus_id_solver) -= tmp;
            }
        }

        template<class T>
        void fill_SBus_imag(CplxMat & Sbuses,
                            const T & structure_data,
                            const RealMat & temporal_data,
                            const std::vector<int> & id_me_to_ac_solver,
                            bool add  // if true call += else calls -=
                            ) const 
        {
            auto nb_el = structure_data.nb();
            const auto & el_status = structure_data.get_status();
            const auto & el_bus_id = structure_data.get_bus_id();
            int  bus_id_solver, bus_id_me;
            cplx_type tmp;
            for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                bus_id_solver = id_me_to_ac_solver[bus_id_me];
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver) += BaseConstants::my_i * tmp;
                else Sbuses.col(bus_id_solver) -= BaseConstants::my_i * tmp;;
            }
        }

        bool compute_one_powerflow(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                   CplxVect & V,
                                   const CplxVect & Sbus,
                                   const Eigen::VectorXi & bus_pv,
                                   const Eigen::VectorXi & bus_pq,
                                   const std::vector<int> & id_ac_solver_to_me,
                                   int max_iter,
                                   double tol
                                   );
    private:
        // inputs
        GridModel _grid_model;
        CplxMat _Sbuses;
        // outputs
        RealMat _amps_flows;
        CplxMat _voltages;
        int _nb_solved;
        int _status;

        // parameters
        bool _compute_flows;
        ChooseSolver _solver;

        //timers
        double _timer_total;
        double _timer_solver;
        double _timer_pre_proc;
};
#endif  //COMPUTERS_H
