// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEMULTIPLEPOWERFLOW_H
#define BASEMULTIPLEPOWERFLOW_H

#include "GridModel.h"

/**
This is a utility class, used for Computers and SecurityAnalysis that abstract some computations when
the same solver is re used multiple times.
 **/
class BaseMultiplePowerflow
{
    public:
        typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RealMat;
        typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CplxMat;
        
        BaseMultiplePowerflow(const GridModel & init_grid_model):
            _grid_model(init_grid_model),
            n_line_(init_grid_model.nb_powerline()),
            n_trafos_(init_grid_model.nb_trafo()),
            n_total_(n_line_ + n_trafos_),
            _solver(),
            _amps_flows(),
            _voltages(),
            _nb_solved(0),
            _timer_compute_A(0.),
            _timer_solver(0.)
            {
                // make sure that my "grid_model" is ready to be used (for ac and dc)
                Eigen::Index nb_bus = init_grid_model.total_bus();
                CplxVect V = CplxVect::Constant(nb_bus, 1.04);
                // const auto & Vtmp = init_grid_model.get_V();
                // for(int i = 0; i < Vtmp.size(); ++i) V[i] = Vtmp[i];
                _grid_model.dc_pf(V, 10, 1e-5);
                _grid_model.ac_pf(V, 10, 1e-5);
                
                // assign the right solver type
                _solver.change_solver(_grid_model.get_solver_type());
            }

        BaseMultiplePowerflow(const BaseMultiplePowerflow&) = delete;
    
        // solver "control"
        void change_solver(const SolverType & type){
            _solver.change_solver(type);
        }
        std::vector<SolverType> available_solvers() {return _solver.available_solvers(); }
        SolverType get_solver_type() {return _solver.get_type(); }

        // utlities informations
        double amps_computation_time() const {return _timer_compute_A;}
        double solver_time() const {return _timer_solver;}
        int nb_solved() const {return _nb_solved;}

        // results
        Eigen::Ref<const RealMat > get_flows() const {return _amps_flows;}
        Eigen::Ref<const CplxMat > get_voltages() const {return _voltages;}
        
    protected:
        template<class T>
        void compute_amps_flows(const T & structure_data,
                                real_type sn_mva,
                                Eigen::Index lag_id) 
        {
            const auto & bus_vn_kv = _grid_model.get_bus_vn_kv();
            const auto & el_status = structure_data.get_status();
            const auto & bus_from = structure_data.get_bus_from();
            const auto & bus_to = structure_data.get_bus_to();
            const auto & v_yac_ff = structure_data.yac_ff();
            const auto & v_yac_ft = structure_data.yac_ft();
            Eigen::Index nb_el = structure_data.nb();
            real_type sqrt_3 = sqrt(3.);

            for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                // retrieve which buses are used
                int bus_from_me = bus_from(el_id);
                int bus_to_me = bus_to(el_id);

                // retrieve voltages
                const auto Efrom = _voltages.col(bus_from_me);  // vector (one voltages per step)
                const auto Eto = _voltages.col(bus_to_me);

                const real_type bus_vn_kv_f = bus_vn_kv(bus_from_me);
                const RealVect v_f_kv = Efrom.array().abs() * bus_vn_kv_f;

                // retrieve physical parameters
                const cplx_type yac_ff = v_yac_ff(el_id);  // scalar
                const cplx_type yac_ft = v_yac_ft(el_id);

                // trafo equations (to get the power at the "from" side)
                CplxVect I_ft =  yac_ff * Efrom + yac_ft * Eto;
                I_ft = I_ft.array().conjugate();
                const CplxVect S_ft = Efrom.array() * I_ft.array();

                // now compute the current flow
                RealVect res = S_ft.array().abs() * sn_mva;
                res.array() /= sqrt_3 * v_f_kv.array();
                _amps_flows.col(el_id + lag_id) = res;
            }
        }

        bool compute_one_powerflow(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                   CplxVect & V,
                                   const CplxVect & Sbus,
                                   const Eigen::VectorXi & bus_pv,
                                   const Eigen::VectorXi & bus_pq,
                                   int max_iter,
                                   double tol
                                   );

        void compute_flows_from_Vs();

        CplxVect extract_Vsolver_from_Vinit(const CplxVect& Vinit,
                                            Eigen::Index nb_buses_solver,
                                            Eigen::Index nb_total_bus,
                                            const std::vector<int> & id_me_to_ac_solver){
            // extract V solver from the given V
            CplxVect Vinit_solver = CplxVect::Constant(nb_buses_solver, {_grid_model.get_init_vm_pu(), 0.});
            Eigen::Index tmp;
            for(Eigen::Index bus_id_grid = 0; bus_id_grid < nb_total_bus; ++bus_id_grid){
                tmp = id_me_to_ac_solver[bus_id_grid];
                if(tmp == GridModel::_deactivated_bus_id) continue;
                Vinit_solver[tmp] = Vinit[bus_id_grid];
            }
            return Vinit_solver;
        }
    protected:
        // inputs
        GridModel _grid_model;

        // properties of the grid
        const Eigen::Index n_line_;
        const Eigen::Index n_trafos_;
        const Eigen::Index n_total_;

        // solver
        ChooseSolver _solver;

        // outputs
        RealMat _amps_flows;
        CplxMat _voltages;
        
        // timers
        int _nb_solved;
        double _timer_compute_A;
        double _timer_solver;

};
#endif // BASEMULTIPLEPOWERFLOW_H
