// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef COMPUTERS_H
#define COMPUTERS_H

#include "BaseBatchSolverSynch.h"

/**
Allws the computation of time series, that is, the same grid topology is used along with time
series of injections (productions and loads) to compute powerflows/
 **/
class TimeSeries: public BaseBatchSolverSynch
{
    public:
        TimeSeries(const GridModel & init_grid_model):
            BaseBatchSolverSynch(init_grid_model),
            _Sbuses(),
            _status(1), // 1: success, 0: failure
            _compute_flows(true),
            _timer_total(0.) ,
            _timer_pre_proc(0.)
            {}

        TimeSeries(const TimeSeries&) = delete;

        // control on whether I compute the flows or not
        void deactivate_flow_computations() {_compute_flows = false;}
        void activate_flow_computations() {_compute_flows = true;}

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}

        // status
        int get_status() const {return _status;}

        /**
        This function computes the results of running as many powerflow when varying the 
        injection (Sbus). 

        Each line of `Sbuses` will be a time step, and each column with 
        **/
        int compute_Vs(Eigen::Ref<const RealMat> gen_p,
                       Eigen::Ref<const RealMat> sgen_p,
                       Eigen::Ref<const RealMat> load_p,
                       Eigen::Ref<const RealMat> load_q,
                       const CplxVect & Vinit,
                       const int max_iter,
                       const real_type tol);
                       
        virtual void clear(){
            BaseBatchSolverSynch::clear();
            _Sbuses = CplxMat();
            _status = 1;
            _compute_flows = true;
            _timer_total = 0.;
            _timer_pre_proc = 0.;
        }

        Eigen::Ref<const CplxMat > get_sbuses() const {return _Sbuses;}
        Eigen::Ref<const RealMat > compute_flows() {
            compute_flows_from_Vs();
            return _amps_flows;
        }
        Eigen::Ref<const RealMat > compute_power_flows() {
            compute_flows_from_Vs(false);
            return _active_power_flows;
        }

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
            for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;
                bus_id_me = el_bus_id(el_id);
                bus_id_solver = id_me_to_ac_solver[bus_id_me];
                const auto & tmp = temporal_data.col(el_id).cast<cplx_type>();
                if(add) Sbuses.col(bus_id_solver) += BaseConstants::my_i * tmp;
                else Sbuses.col(bus_id_solver) -= BaseConstants::my_i * tmp;
            }
        }

    private:
        // inputs
        CplxMat _Sbuses;

        // outputs
        int _status;

        // parameters
        bool _compute_flows;

        //timers
        double _timer_total;
        double _timer_pre_proc;
};
#endif  //COMPUTERS_H
