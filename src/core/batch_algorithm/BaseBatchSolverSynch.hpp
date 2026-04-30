// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEMULTIPLEPOWERFLOW_H
#define BASEMULTIPLEPOWERFLOW_H

#include "LSGrid.hpp"

namespace ls2g {

/**
This is a utility class, used for TimeSeries and SecurityAnalysis that abstract some computations when
the same solver is re used multiple times.

It allows to perform "batch" powerflow one a time in a synchronous manner.

The "solver" of the gridmodel is never really used to perform powerflows.

**/
class LS2G_API BaseBatchSolverSynch : protected BaseConstants
{
    public:
        using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using CplxMat =  Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        
        explicit BaseBatchSolverSynch(const LSGrid & init_grid_model) noexcept:
            _grid_model(init_grid_model),
            n_line_(init_grid_model.nb_powerline()),
            n_trafos_(init_grid_model.nb_trafo()),
            n_total_(n_line_ + n_trafos_),
            _algo(),
            _voltages(),
            _amps_flows(),
            _active_power_flows()
            {
            }
        virtual ~BaseBatchSolverSynch() noexcept = default;  // to avoid warning about overload virtual
        BaseBatchSolverSynch(const BaseBatchSolverSynch&) = delete;
        BaseBatchSolverSynch(BaseBatchSolverSynch&&) = delete;
        BaseBatchSolverSynch & operator=(BaseBatchSolverSynch&&) = delete;
        BaseBatchSolverSynch & operator=(const BaseBatchSolverSynch&) = delete;
    
        bool get_init_from_n_powerflow() const noexcept {return _init_from_n_powerflow;}
        void set_init_from_n_powerflow(bool do_it) noexcept {_init_from_n_powerflow = do_it;}

        // solver "control"
        virtual void change_algorithm(const AlgorithmType & type){
            _algo.change_algorithm(type);
            this->clear();
        }
        
        // Returns the enum-typed solvers available in this build.
        // Does not include plugin (Custom) solvers; use AlgorithmRegistry::available_default_algorithms()
        // for the full list.
        std::vector<AlgorithmType> available_default_algorithms() const {return _algo.available_default_algorithms(); }
        
        // Returns all solver names currently registered (built-in + plugins).
        std::vector<std::string> available_algorithm_names() const {
            return AlgorithmRegistry::instance().available_algorithm_names();
        }

        AlgorithmType get_algo_type() const {return _algo.get_type(); }

        // // TODO
        // void change_gridmodel(const GridModel & new_grid_model){
        //     CplxVect V = CplxVect::Constant(new_grid_model.total_bus(), 1.04);
        //     _grid_model = new_grid_model.copy(); // not implemented !!!
        // }

        // utlities informations
        double amps_computation_time() const {return _timer_compute_A;}
        double solver_time() const {return _timer_solver;}
        int nb_solved() const {return _nb_solved;}
        virtual void clear() {
            _algo.reset();
            _amps_flows = RealMat();
            _active_power_flows = RealMat();
            _voltages = CplxMat();
            _nb_solved = 0;
            _timer_compute_A = 0.;
            _timer_compute_P = 0.;
            _timer_solver = 0.;
        }

        // results 
        // this should not be const, see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html#pass-by-reference
        // tl;dr: const can make copies ! OR NOT I AM LOST
        const RealMat & get_flows() const {return _amps_flows;}
        const RealMat & get_power_flows() const {return _active_power_flows;}
        const CplxMat & get_voltages() const {return _voltages;}
        
    protected:
        template<class T>
        void compute_amps_flows(const T & structure_data,
                                real_type sn_mva,
                                size_t lag_id,
                                bool is_trafo) 
        {
            const auto & bus_vn_kv = _grid_model.get_bus_vn_kv();
            const auto & el_status = structure_data.get_status_global();
            const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
            const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
            bool is_ac = _algo.ac_solver_used();

            const auto & vect_y_ff = is_ac ? structure_data.yac_eff_11() : structure_data.ydc_11();
            const auto & vect_y_ft = is_ac ? structure_data.yac_eff_12() : structure_data.ydc_12();
            Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift(); // not used in AC nor if it's powerline anyway

            size_t nb_el = structure_data.nb();
            real_type sqrt_3 = sqrt(3.);

            RealVect res;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                // retrieve which buses are used
                GlobalBusId bus_from_me = bus_from(el_id);
                GlobalBusId bus_to_me = bus_to(el_id);

                // retrieve voltages
                const auto Efrom = _voltages.col(bus_from_me.cast_int());  // vector (one voltages per step)
                const auto Eto = _voltages.col(bus_to_me.cast_int());

                const real_type bus_vn_kv_f = bus_vn_kv(bus_from_me.cast_int());
                const RealVect v_f_kv = Efrom.array().abs() * bus_vn_kv_f;

                // retrieve physical parameters
                const cplx_type y_ff = vect_y_ff(el_id);  // scalar
                const cplx_type y_ft = vect_y_ft(el_id);

                if(is_ac){
                    // trafo equations (to get the power at the "from" side)
                    CplxVect I_ft =  y_ff * Efrom + y_ft * Eto;
                    I_ft = I_ft.array().conjugate();
                    const CplxVect S_ft = Efrom.array() * I_ft.array();

                    // now compute the current flow
                    res = S_ft.array().abs() * sn_mva;
                }else{
                    res = (std::real(y_ff) * Efrom.array().arg() + std::real(y_ft) * Eto.array().arg()) * sn_mva;
                    if(is_trafo) res.array() -= dc_x_tau_shift(el_id);
                    res.array() = res.array().abs();
                }
                res.array() /= sqrt_3 * v_f_kv.array();
                _amps_flows.col(el_id + lag_id) = res;
            }
        }
        template<class T>
        void compute_active_power_flows(const T & structure_data,
                                        real_type sn_mva,
                                        size_t lag_id,
                                        bool is_trafo) 
        {
            const auto & bus_vn_kv = _grid_model.get_bus_vn_kv();
            const auto & el_status = structure_data.get_status_global();
            const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
            const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
            const bool is_ac = _algo.ac_solver_used();

            Eigen::Ref<const CplxVect> vect_y_ff = is_ac ? structure_data.yac_eff_11() : structure_data.ydc_11();
            Eigen::Ref<const CplxVect> vect_y_ft = is_ac ? structure_data.yac_eff_12() : structure_data.ydc_12();
            Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift(); // not used in AC nor if it's powerline anyway

            Eigen::Index nb_el = structure_data.nb();

            RealVect res;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                // retrieve which buses are used
                GlobalBusId bus_from_me = bus_from(el_id);
                GlobalBusId bus_to_me = bus_to(el_id);

                // retrieve voltages
                const auto & Efrom = _voltages.col(bus_from_me.cast_int());  // vector (one voltages per step)
                const auto & Eto = _voltages.col(bus_to_me.cast_int());

                const real_type bus_vn_kv_f = bus_vn_kv(bus_from_me.cast_int());
                const RealVect v_f_kv = Efrom.array().abs() * bus_vn_kv_f;

                // retrieve physical parameters
                const cplx_type y_ff = vect_y_ff(el_id);  // scalar
                const cplx_type y_ft = vect_y_ft(el_id);

                // trafo equations (to get the power at the "from" side)
                if(is_ac){
                    CplxVect I_ft = y_ff * Efrom + y_ft * Eto;
                    I_ft = I_ft.array().conjugate();
                    const CplxVect S_ft = Efrom.array() * I_ft.array();

                    // now compute the active flow
                    res = S_ft.array().real() * sn_mva;
                }else{
                    res = (std::real(y_ff) * Efrom.array().arg() + std::real(y_ft) * Eto.array().arg()) * sn_mva;
                    if(is_trafo) res.array() -= dc_x_tau_shift(el_id);
                }
                _active_power_flows.col(el_id + lag_id) = res;
            }
        }

        bool compute_one_powerflow(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                   CplxVect & V,
                                   const CplxVect & Sbus,
                                   Eigen::Ref<const IntVect> slack_ids,
                                   const RealVect & slack_weights,
                                   Eigen::Ref<const IntVect> bus_pv,
                                   Eigen::Ref<const IntVect> bus_pq,
                                   int max_iter,
                                   double tol
                                   );

        void compute_flows_from_Vs(bool amps=true);

        CplxVect extract_Vsolver_from_Vinit(const CplxVect& Vinit,
                                            size_t nb_buses_solver,
                                            size_t nb_total_bus,
                                            const SolverBusIdVect & id_me_to_ac_solver){
            // extract V solver from the given V
            CplxVect Vinit_solver = CplxVect::Constant(nb_buses_solver, {_grid_model.get_init_vm_pu(), 0.});
            int tmp;
            for(size_t bus_id_grid = 0; bus_id_grid < nb_total_bus; ++bus_id_grid){
                tmp = static_cast<int>(id_me_to_ac_solver[bus_id_grid]);
                if(tmp == BaseConstants::_deactivated_bus_id) continue;
                Vinit_solver[tmp] = Vinit[bus_id_grid];
            }
            return Vinit_solver;
        }
    protected:

        CplxVect prepare_solver_input_base(const CplxVect & Vinit, bool ac_solver_used){
            // clear previous data
            Sbus_ = CplxVect();
            Ybus_ = Eigen::SparseMatrix<cplx_type>();
            id_solver_to_me_.clear();
            id_me_to_solver_.clear();
            slack_ids_solver_ = SolverBusIdVect();
            slack_ids_me_ = GlobalBusIdVect();
            nb_buses_solver_ = -1;

            // fill the data correctly
            _algo_controler.tell_all_changed();
            CplxVect res = _grid_model.pre_process_solver(
                Vinit, 
                Sbus_,
                Ybus_,
                id_me_to_solver_,
                id_solver_to_me_,
                slack_ids_me_,
                slack_ids_solver_,
                ac_solver_used,
                _algo_controler);

            // extract relevant information
            nb_buses_solver_ = static_cast<int>(Ybus_.cols());
            const SolverBusIdVect & gm_bus_pv = _grid_model.get_pv_solver();
            const SolverBusIdVect & gm_bus_pq = _grid_model.get_pq_solver();
            const RealVect & gm_bus_sw = _grid_model.get_slack_weights_solver();

            // TODO copies are made here, which is not ideal
            bus_pv_ = gm_bus_pv;  // was _to_intvect()
            bus_pq_ = gm_bus_pq;  // was _to_intvect()
            slack_weights_ = gm_bus_sw;
            return res;
        }

        size_t _reset_data_and_check_vinit(const CplxVect & Vinit){
            const size_t nb_total_bus = _grid_model.total_bus();
            if(Vinit.size() != nb_total_bus){
                std::ostringstream exc_;
                exc_ << "TimeSeries::compute_Sbuses: Size of the Vinit should be the same as the total number of buses. Currently:  ";
                exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_total_bus << " buses.";
                exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
                throw std::runtime_error(exc_.str());
            }

            // reset timers
            _nb_solved = 0;
            _timer_pre_proc = 0.;
            _timer_total = 0.;
            _timer_solver = 0.;
            return nb_total_bus;
        }

        bool _finish_preprocessing(
            size_t nb_steps,
            size_t nb_total_bus,
            CplxVect & Vinit_solver,  // is modified if _init_from_n_powerflow is true !
            size_t max_iter,
            real_type tol,
            CustTimer  & timer_preproc  // non const because double duration() is not const
        ){

                // init the results matrices
                _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus); 
                _amps_flows = RealMat::Zero(0, n_total_);
                _active_power_flows = RealMat::Zero(0, n_total_);

                // reset the solver
                _algo.reset();

                // perform the initial powerflow / "powerflow in n"
                // (needed to init the underlying solver with the correct sparsity pattern in particular)
                _algo_controler.tell_all_changed();
                _algo.tell_solver_control(_algo_controler);
                _grid_model.get_generators().set_vm(Vinit_solver, id_me_to_solver_);
                CplxVect Vinit_solver2 = Vinit_solver;
                bool conv = _algo.compute_pf(
                    Ybus_,
                    Vinit_solver2,
                    Sbus_,
                    slack_ids_me_.as_eigen(),
                    slack_weights_,
                    bus_pv_.as_eigen(),
                    bus_pq_.as_eigen(),
                    max_iter,
                    tol);

                // check if we init the n-1 cases with results from the n cases
                // or not
                if(_init_from_n_powerflow) Vinit_solver = _algo.get_V();

                // everything init from n-case above
                _algo_controler.tell_none_changed();
                
                // end of pre processing
                _timer_pre_proc = timer_preproc.duration();
                return conv;
        }

    protected:
        bool _init_from_n_powerflow = false;
        //timers
        double _timer_total = 0.;
        double _timer_pre_proc = 0.;

        // inputs
        LSGrid _grid_model;

        // properties of the grid
        const size_t n_line_;
        const size_t n_trafos_;
        const size_t n_total_;

        // algo
        AlgorithmSelector _algo;

        // outputs
        CplxMat _voltages;
        RealMat _amps_flows;
        RealMat _active_power_flows;
        
        // timers
        int _nb_solved = 0;
        double _timer_compute_A = 0.;
        double _timer_compute_P = 0.;
        double _timer_solver = 0.;

        // solver control
        AlgoControl _algo_controler;

        // internal data
        CplxVect Sbus_;
        Eigen::SparseMatrix<cplx_type> Ybus_;
        GlobalBusIdVect id_solver_to_me_;
        SolverBusIdVect id_me_to_solver_;
        SolverBusIdVect slack_ids_solver_;
        GlobalBusIdVect slack_ids_me_;
        int nb_buses_solver_;

        // TODO everything is copied here
        // not optimal...
        SolverBusIdVect bus_pv_;
        SolverBusIdVect bus_pq_;
        RealVect slack_weights_;

};


} // namespace ls2g

#endif // BASEMULTIPLEPOWERFLOW_H