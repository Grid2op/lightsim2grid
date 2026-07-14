// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TRAFO_CONTAINER_H
#define TRAFO_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SubstationContainer.hpp"
#include "OneSideContainer_forBranch.hpp"
#include "TwoSidesContainer_rxh_A.hpp"

#include <algorithm>

namespace ls2g {

class TrafoContainer;
class LS2G_API TrafoInfo : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::TwoSidesContainer_rxh_AInfo
{
    public:
        // members
        real_type ratio;
        real_type shift_rad;
        bool is_tap_side1;

        inline TrafoInfo(const TrafoContainer & r_data_trafo, int my_id) noexcept;
};

/**
This class is a container for all transformers on the grid.
Transformers are modeled "in pi" here. If your trafo are given in a "t" model (like in pandapower
for example) use the DataConverter class.

The convention used for the transformer is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/trafo.html#electric-model
**/
class LS2G_API TrafoContainer final : public TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>, public IteratorAdder<TrafoContainer, TrafoInfo>
{
    //////////////////////////////
    // access data from base class
    public:
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::get_buses_side_1;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::get_buses_side_2;

    protected:
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::side_1_;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::side_2_;
        using TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::status_global_;
    //////////////////////////////

    friend class TrafoInfo;

    public:
        using DataInfo = TrafoInfo;

    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
                   TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>::StateRes,
                   std::vector<real_type>, // ratio_
                   std::vector<bool> , // is_tap_hv_side
                   std::vector<real_type>, // shift_
                   bool,  // ignore_tap_side_for_shift_
                   bool,  // shift_dependent_rx_
                   std::vector<real_type>,  // base_r_
                   std::vector<real_type>,  // base_x_
                   std::vector<std::vector<real_type> >,  // rx_corr_alpha_
                   std::vector<std::vector<real_type> >   // rx_corr_pct_
               >;

        TrafoContainer() noexcept = default;
        virtual ~TrafoContainer() noexcept = default;

        void init(const Eigen::Ref<const RealVect> & trafo_r,
                  const Eigen::Ref<const RealVect> & trafo_x,
                  const Eigen::Ref<const CplxVect> & trafo_b,
                  const Eigen::Ref<const RealVect> & trafo_tap_step_pct,
                  const Eigen::Ref<const RealVect> & trafo_tap_pos,
                  const Eigen::Ref<const RealVect> & trafo_shift_degree,
                  const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                  const Eigen::Ref<const Eigen::VectorXi> & trafo_hv_id,
                  const Eigen::Ref<const Eigen::VectorXi> & trafo_lv_id,
                  bool ignore_tap_side_for_shift
                  );

        void init(const Eigen::Ref<const RealVect> & trafo_r,
                  const Eigen::Ref<const RealVect> & trafo_x,
                  const Eigen::Ref<const CplxVect> & trafo_b,
                  const Eigen::Ref<const RealVect> & trafo_ratio,
                  const Eigen::Ref<const RealVect> & trafo_shift_degree,
                  const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                  const Eigen::Ref<const Eigen::VectorXi> & trafo_hv_id,
                  const Eigen::Ref<const Eigen::VectorXi> & trafo_lv_id,
                  bool ignore_tap_side_for_shift
                  );

        //pickle
        StateRes get_state() const;
        void set_state(StateRes & my_state );

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static TrafoContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "TrafoContainer"; }  // written into / checked against the binary file header

        bool ignore_tap_side_for_shift() const { return ignore_tap_side_for_shift_; }

        /**
         * Declare that the series impedance (r, x) of (some) transformers depends on
         * the phase-shift angle `alpha` (= `shift_`), and supply that dependency as a
         * per-transformer table of sample points `alpha (rad) -> r/x correction (%)`
         * (the per-step r/x deltas of a pypowsybl phase-tap-changer; r% == x%). The
         * effective impedance is `base * (1 + corr(shift_) / 100)`, recomputed (by
         * interpolation on `shift_`) every time the coefficients are rebuilt -- in
         * particular whenever `change_shift` / `change_ratio` is called -- so there is
         * NO notion of a discrete "tap" in lightsim2grid. Pass an empty inner vector
         * for a transformer that has no such dependency. `enable` is the master flag
         * (kept false for pandapower, which has no such data).
         */
        void set_shift_dependent_rx(bool enable,
                                    const std::vector<std::vector<real_type> > & alpha_rad,
                                    const std::vector<std::vector<real_type> > & rx_corr_pct,
                                    DualAlgoControl & solver_control);

        void hack_Sbus_for_dc_phase_shifter(
            CplxVect & Sbus,
            bool ac,
            const SolverBusIdVect & id_grid_to_solver);  // needed for dc mode

        void compute_results(const Eigen::Ref<const RealVect> & Va,
                             const Eigen::Ref<const RealVect> & Vm,
                             const Eigen::Ref<const CplxVect> & V,
                             const SolverBusIdVect & id_grid_to_solver,
                             const Eigen::Ref<const RealVect> & bus_vn_kv,
                             real_type sn_mva,
                             bool ac)
        {
            // compute base values
            compute_results_tsc_rxha_no_amps(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
            // adjust for phase shifters
            if(!ac){
                Eigen::Ref<RealVect> res_p_side_1 = get_res_p_side_1();
                Eigen::Ref<RealVect> res_p_side_2 = get_res_p_side_2();
                const std::vector<bool> & status1 = side_1_.get_status();
                const std::vector<bool> & status2 = side_2_.get_status();

                const int nb_element = nb();
                for(int el_id = 0; el_id < nb_element; ++el_id){
                    if(status_global_[el_id] && status1[el_id] && status2[el_id]){
                        res_p_side_1(el_id) += dc_x_tau_shift_(el_id) * sn_mva;
                        res_p_side_2(el_id) -= dc_x_tau_shift_(el_id) * sn_mva;
                    }
                }
            }
            // compute amps flow
            compute_amps_after_all_set();
        }
        
        void reset_results(){
            reset_results_tsc_rxha();
        }

        Eigen::Ref<const RealVect> dc_x_tau_shift() const {return dc_x_tau_shift_;}

        void change_ratio(
            int el_id,
            real_type new_ratio,
            DualAlgoControl & solver_control){
                // el_id indexes ratio_ with an unchecked Eigen operator() below (OOB write).
                _check_in_range(el_id, ratio_, "change_ratio");
                if(std::abs(ratio_(el_id) - new_ratio) >_tol_equal_float){
                    ratio_(el_id) = new_ratio;
                    // TODO speed: only some part needs to be recomputed
                    _update_internal_coeffs(el_id); 
                    solver_control.ac_algo_controler().tell_recompute_ybus(); solver_control.dc_algo_controler().tell_recompute_ybus();
                }
        }
        
        /**
         * The shift is in radian (not degree !)
         * 
         * It is the shift on the "side 1" (regardless of the value of "is_tap_hv_side").
         * If the tap is on the other side, the user has the reponsibility to
         * take the opposite (ie -0.1 instead of +0.1)
         */
        void change_shift(
            int el_id,
            real_type new_shift_rad,
            DualAlgoControl & solver_control){
                // el_id indexes shift_ with an unchecked Eigen operator() below (OOB write).
                _check_in_range(el_id, shift_, "change_shift");
                if(std::abs(shift_(el_id) - new_shift_rad) >_tol_equal_float){
                    shift_(el_id) = new_shift_rad;
                    // TODO speed: only some part needs to be recomputed
                    _update_internal_coeffs(el_id); 
                    solver_control.ac_algo_controler().tell_recompute_ybus(); solver_control.dc_algo_controler().tell_recompute_ybus();
                    solver_control.dc_algo_controler().tell_recompute_sbus();  // only in DC however
                }
        }
        
    protected:
        // void _update_model_coeffs();
        virtual void _update_model_coeffs_one_el(int el_id) override;
        virtual void _update_other_model_coeffs() override {
            dc_x_tau_shift_ = RealVect::Zero(nb());
        }

        virtual bool _deactivate(int el_id, DualAlgoControl & solver_control) override {
            bool has_been_changed = TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>:: _deactivate(el_id, solver_control);
            if(has_been_changed){
                // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
            return has_been_changed;
        }

        virtual bool _reactivate(int el_id, DualAlgoControl & solver_control) override {
            bool has_been_changed = TwoSidesContainer_rxh_A<OneSideContainer_ForBranch>:: _reactivate(el_id, solver_control);
            if(has_been_changed){
                // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
            return has_been_changed;
        }

        virtual void _change_bus_side_1(
            int /*el_id*/,
            GridModelBusId /*new_gridmodel_bus_id*/,
            DualAlgoControl & solver_control,
            const SubstationContainer & /*substation*/,
            bool has_effectively_changed
        ) override {
            if(has_effectively_changed){
                // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
        }

        virtual void _change_bus_side_2(
            int /*el_id*/,
            GridModelBusId /*new_gridmodel_bus_id*/,
            DualAlgoControl & solver_control,
            const SubstationContainer & /*substation*/,
            bool has_effectively_changed
        ) override {
            if(has_effectively_changed){
                // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
        }

        virtual void _update_topo(
            DualAlgoControl & solver_control,
            SubstationContainer & /*substations*/,
            const std::vector<bool> & side1_changed,
            const std::vector<bool> & side2_changed
        ) override
        {
            bool onechanged_1 = std::any_of(side1_changed.begin(), side1_changed.end(), [](bool v) { return v; });
            bool onechanged_2 = std::any_of(side2_changed.begin(), side2_changed.end(), [](bool v) { return v; });
            if(onechanged_1 || onechanged_2){
                // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
        }
    private:
        /**
         * whether to ignore the tap position for phase shifter (alpha).
         * 
         * This is the default behaviour in pandapower, where the phase shifter
         * is always assigned to side 1.
         *
         * Default-initialized: a container whose init_trafo() is never called
         * (grid without trafos) is still copied and serialized, and an
         * indeterminate bool there is undefined behavior (found by valgrind
         * over the C++ unit tests).
         */
        bool ignore_tap_side_for_shift_ = false;
        
        // physical properties
        std::vector<bool> is_tap_side1_;  // whether the tap is hav side or not

        // input data
        RealVect ratio_;  // transformer ratio (no unit) (depends on is_tap_side1_)
        RealVect shift_;  // phase shifter (in radian !) (might depends on is_tap_side1, if ignore_tap_side_for_shift_ is true, then it is the shift side1)

        // alpha-dependent series-impedance correction (phase-shifting transformers).
        // lightsim2grid has no "tap" concept: the per-step r/x correction is stored
        // as a function of the phase-shift `alpha` and looked up by the current
        // `shift_`. Disabled (flag false, empty tables) for pandapower.
        // Default-initialized for the same reason as ignore_tap_side_for_shift_.
        bool shift_dependent_rx_ = false;
        RealVect base_r_;  // neutral (uncorrected) r, per trafo
        RealVect base_x_;  // neutral (uncorrected) x, per trafo
        std::vector<std::vector<real_type> > rx_corr_alpha_;  // per trafo: alpha (rad), ascending
        std::vector<std::vector<real_type> > rx_corr_pct_;    // per trafo: r/x correction (%) at each alpha (r% == x%)

        //output data

        // model coefficients
        RealVect dc_x_tau_shift_;

        // r/x correction (%) at the current `shift_(el_id)`, linearly interpolated on
        // the stored `alpha -> correction` samples (clamped outside the range). 0 if
        // the transformer carries no such dependency.
        real_type _shift_rx_corr_pct(int el_id) const {
            const std::vector<real_type> & xs = rx_corr_alpha_[el_id];
            const std::vector<real_type> & ys = rx_corr_pct_[el_id];
            const std::size_t n = xs.size();
            if(n == 0) return my_zero_;
            const real_type a = shift_(el_id);
            if(a <= xs.front()) return ys.front();
            if(a >= xs.back()) return ys.back();
            std::size_t hi = 1;
            while(hi < n && xs[hi] < a) ++hi;
            const real_type t = (a - xs[hi - 1]) / (xs[hi] - xs[hi - 1]);
            return ys[hi - 1] + t * (ys[hi] - ys[hi - 1]);
        }

    protected:

        virtual real_type fillBf_for_PTDF_coeff(int tr_id) const override {
            real_type res = x_(tr_id);
            real_type tau = is_tap_side1_[tr_id] ? ratio_(tr_id) : 1. / ratio_(tr_id);
            return res * tau;
        }

        virtual int fillBf_for_PTDF_id(int tr_id, int nb_powerline) const override {
            return tr_id + nb_powerline;
        }

        virtual FDPFCoeffs get_fdpf_coeffs(int tr_id, FDPFMethod xb_or_bx) const override;
};

inline TrafoInfo::TrafoInfo(const TrafoContainer & r_data_trafo, int my_id) noexcept:
TwoSidesContainer_rxh_AInfo(r_data_trafo, my_id),
ratio(-1.0),
shift_rad(-1.0),
is_tap_side1(true)
{
    if(my_id < 0) return;
    if(static_cast<size_t>(my_id) >= r_data_trafo.nb()) return;
    is_tap_side1 = r_data_trafo.is_tap_side1_[my_id];
    ratio = r_data_trafo.ratio_.coeff(my_id);
    shift_rad = r_data_trafo.shift_.coeff(my_id);
}


} // namespace ls2g

#endif  //TRAFO_CONTAINER_H
