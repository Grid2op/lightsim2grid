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
#include "BranchContainer.hpp"
#include "TapChangers.hpp"

#include <algorithm>
#include <array>
#include <utility>

namespace ls2g {

class TrafoContainer;
class LS2G_API TrafoInfo : public BranchContainer::BranchInfo
{
    public:
        // members
        real_type ratio;
        real_type shift_rad;
        bool is_tap_side1;
        // its ratio tap changer (has_ratio_tap_changer false: none, the rest is meaningless)
        bool has_ratio_tap_changer;
        int ratio_tap_position;
        int ratio_low_tap;
        int ratio_high_tap;
        RegulationMode ratio_regulation_mode;
        bool ratio_regulating;
        real_type ratio_target;  // a voltage in pu of the regulated bus' nominal voltage
        real_type ratio_deadband;
        int ratio_regulated;  // the regulated grid bus id
        // its phase tap changer, the same
        bool has_phase_tap_changer;
        int phase_tap_position;
        int phase_low_tap;
        int phase_high_tap;
        RegulationMode phase_regulation_mode;
        bool phase_regulating;
        real_type phase_target;  // MW (ACTIVE_POWER) or A (CURRENT_LIMITER)
        real_type phase_deadband;
        int phase_regulated;  // the side of the transformer it regulates (1 or 2)
        // the position of its phase tap changer in the last results (the input one unless an
        // outer loop moved it, see TrafoContainer::set_results_phase_tap_override)
        int res_phase_tap_position;

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
class LS2G_API TrafoContainer final : public BranchContainer, public IteratorAdder<TrafoContainer, TrafoInfo>
{
    friend class TrafoInfo;

    public:
        using DataInfo = TrafoInfo;

    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
                   BranchContainer::StateRes,
                   std::vector<real_type>, // ratio_
                   std::vector<bool> , // is_tap_hv_side
                   std::vector<real_type>, // shift_
                   bool,  // ignore_tap_side_for_shift_
                   bool,  // shift_dependent_rx_
                   std::vector<real_type>,  // base_r_
                   std::vector<real_type>,  // base_x_
                   std::vector<std::vector<real_type> >,  // rx_corr_alpha_
                   std::vector<std::vector<real_type> >,  // rx_corr_pct_
                   TapChangers::StateRes,  // ratio_taps_
                   TapChangers::StateRes,  // phase_taps_
                   std::vector<real_type>,  // base_ratio_
                   std::vector<cplx_type>,  // base_h1_
                   std::vector<cplx_type>   // base_h2_
               >;
        enum StateResIdx {
            BRANCH_STATE = 0,
            RATIO,
            IS_TAP_SIDE1,
            SHIFT,
            IGNORE_TAP_SIDE_FOR_SHIFT,
            SHIFT_DEPENDENT_RX,
            BASE_R,
            BASE_X,
            RX_CORR_ALPHA,
            RX_CORR_PCT,
            RATIO_TAPS,
            PHASE_TAPS,
            BASE_RATIO,
            BASE_H1,
            BASE_H2,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "TrafoContainer::StateRes and StateResIdx do not match");

        TrafoContainer() noexcept = default;
        ~TrafoContainer() noexcept override = default;

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

        /**
         * The ratio (`phase` false) or phase (`phase` true) tap changer of transformer `el`:
         * its steps, one per position from `low_tap` (`rho`, `alpha_deg` -- in degree, empty
         * for a ratio changer --, and the corrections of r, x, g, b in %), and its current
         * position. The transformer's pi model is then taken at its tap positions, as
         * OpenLoadFlow takes it (Transformers.getTapCharacteristics): its r, x, g, b are the
         * neutral ones times (1 + step % / 100) of each changer, its ratio the neutral one
         * times the rho of each, its shift the alpha of its phase changer.
         *
         * The neutral values are the r, x, h given at init (or update_physical_parameters)
         * and the ratio the transformer has when this is called divided by the rho of its
         * changers at their positions: the ratio is taken to be at the taps already, this
         * changer's included when it is the first of its kind (give init the ratio at the
         * current taps, as pypowsybl's `rho` is), the replaced one's otherwise.
         */
        void set_tap_changer(bool phase, int el, int low_tap, int position,
                             const std::vector<real_type> & rho,
                             const std::vector<real_type> & alpha_deg,
                             const std::vector<real_type> & r_pct,
                             const std::vector<real_type> & x_pct,
                             const std::vector<real_type> & g_pct,
                             const std::vector<real_type> & b_pct,
                             DualAlgoControl & solver_control);
        /// what that tap changer regulates (see RegulationMode for `target`'s unit), data only
        void set_tap_regulation(bool phase, int el, RegulationMode mode, bool regulating,
                                real_type target, real_type deadband, int regulated);
        /// move that tap changer to `position`: the pi model follows, ratio and shift included
        void change_tap_position(bool phase, int el, int position, DualAlgoControl & solver_control);
        /**
         * The position of that tap changer whose ratio (`phase` false) or shift (`phase` true,
         * rad) is closest to `value`, the other changer staying where it is: OpenLoadFlow's
         * ClosestTapPositionFinder (PiModelArray), the current position kept unless another
         * is strictly closer, the lowest such position otherwise.
         */
        int closest_tap_position(bool phase, int el, real_type value) const;
        const TapChangers & get_tap_changers(bool phase) const { return phase ? phase_taps_ : ratio_taps_; }

        /**
         * The pi model coefficients (y11, y12, y21, y22, pu, the raw block of a transformer
         * connected at both ends) of transformer `el` with its phase tap changer at
         * `phase_position` and the shift `shift_rad` (that position's alpha, or a continuous
         * one: OpenLoadFlow's PiModelArray overrides alpha only, r, x, g, b stay the tap's);
         * its ratio changer where it is. At the current position and shift, exactly the block
         * the container stamps (the shift-dependent impedance of set_shift_dependent_rx
         * aside, which is not OpenLoadFlow's).
         */
        std::array<cplx_type, 4> pi_block_at(int el, int phase_position, real_type shift_rad) const;
        /// the sign the shift enters the block with: d y12 / d shift = j s y12, d y21 / d shift = -j s y21
        real_type shift_sign(int el) const { return (!is_tap_side1_[el] && !ignore_tap_side_for_shift_) ? -1. : 1.; }

        /**
         * Positions of the phase tap changers to compute the next results at (one per
         * transformer, TapChangers' range; anything else keeps the input position): the
         * flows are published at those taps, the inputs are not modified. An outer loop's
         * (see LSGrid::compute_results); empty clears it.
         */
        void set_results_phase_tap_override(const std::vector<int> & positions) { results_phase_tap_override_ = positions; }
        /// the positions of the last results, see TrafoInfo::res_phase_tap_position
        const std::vector<int> & get_res_phase_tap_position() const { return res_phase_tap_position_; }

        void hack_Sbus_for_dc_phase_shifter(
            Eigen::Ref<CplxVect> Sbus,
            bool ac,
            const SolverBusIdVect & id_grid_to_solver);  // needed for dc mode

    protected:
        bool _in_topo_vect() const override { return true; }

        void _compute_results(const Eigen::Ref<const RealVect> & Va,
                              const Eigen::Ref<const RealVect> & Vm,
                              const Eigen::Ref<const CplxVect> & V,
                              const SolverBusIdVect & id_grid_to_solver,
                              const Eigen::Ref<const RealVect> & bus_vn_kv,
                              real_type sn_mva,
                              bool ac) override
        {
            // an outer loop's taps: the flows at those, the inputs untouched
            const std::vector<int> moved = _apply_results_phase_tap_override();
            // compute base values
            _compute_branch_results_no_amps(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
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
            _compute_amps();
            _restore_results_phase_tap_override(moved);
        }

        // see set_results_phase_tap_override: moves the overridden taps (returns the
        // transformers moved, with their input positions), and records the positions
        std::vector<int> _apply_results_phase_tap_override();
        void _restore_results_phase_tap_override(const std::vector<int> & moved);

    public:
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
                    solver_control.tell_recompute_ybus();
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
                    solver_control.tell_recompute_ybus();
                    solver_control.dc_algo_controler().tell_recompute_sbus();  // only in DC however
                }
        }
        
    protected:
        void _update_model_coeffs_one_el(int el_id) override;

        // Any connectivity change moves the phase shifter's DC Sbus term
        // (hack_Sbus_for_dc_phase_shifter only stamps a transformer connected at
        // both ends), on top of the Kron update the branch does.
        // TODO speed: only when dc_x_tau_shift_ is not 0, but be carefull, dc_x_tau_shift_ can be changed later
        void _on_connectivity_changed(int el_id, DualAlgoControl & solver_control) override {
            BranchContainer::_on_connectivity_changed(el_id, solver_control);
            solver_control.dc_algo_controler().tell_recompute_sbus();
        }

        // new r / x / h are the neutral (uncorrected) values. The DC Sbus term of a phase
        // shifter depends on x, so it has to be recomputed too.
        void _on_physical_parameters_updated(DualAlgoControl & solver_control) override {
            base_r_ = r_;
            base_x_ = x_;
            base_h1_ = h_side_1_;
            base_h2_ = h_side_2_;
            solver_control.dc_algo_controler().tell_recompute_sbus();
        }

        // ratio_ and shift_ at the tap positions (the r, x and h follow in
        // _update_model_coeffs_one_el)
        void _apply_tap_position(int el_id) {
            ratio_(el_id) = base_ratio_(el_id) * ratio_taps_.rho(el_id) * phase_taps_.rho(el_id);
            if(phase_taps_.has(el_id)) shift_(el_id) = phase_taps_.alpha(el_id);
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
        CplxVect base_h1_;  // neutral h on side 1, per trafo
        CplxVect base_h2_;  // neutral h on side 2, per trafo
        RealVect base_ratio_;  // ratio with no tap changer step (ratio_ at the taps, see _apply_tap_position)

        // the tap changers, see set_tap_changer (none: neutral)
        TapChangers ratio_taps_;
        TapChangers phase_taps_;
        // see set_results_phase_tap_override, get_res_phase_tap_position (not serialized: results)
        std::vector<int> results_phase_tap_override_;
        std::vector<int> res_phase_tap_position_;
        std::vector<std::pair<int, std::array<real_type, 2> > > results_saved_;  // ratio, shift before the override

        // the AC pi block (y11, y12, y21, y22) of the given parameters
        static std::array<cplx_type, 4> _pi_coeffs(real_type r, real_type x, cplx_type h1, cplx_type h2,
                                                   real_type ratio, real_type shift, bool tap_side1,
                                                   bool ignore_tap_side_for_shift);
        std::vector<std::vector<real_type> > rx_corr_alpha_;  // per trafo: alpha (rad), ascending
        std::vector<std::vector<real_type> > rx_corr_pct_;    // per trafo: r/x correction (%) at each alpha (r% == x%)

        //output data

        // model coefficients
        RealVect dc_x_tau_shift_;

        // r/x correction (%) at the current `shift_(el_id)`, linearly interpolated on
        // the stored `alpha -> correction` samples (clamped outside the range). 0 if
        // the transformer carries no such dependency.
        real_type _shift_rx_corr_pct(int el_id) const { return _rx_corr_pct_at(el_id, shift_(el_id)); }
        real_type _rx_corr_pct_at(int el_id, real_type a) const {
            const std::vector<real_type> & xs = rx_corr_alpha_[el_id];
            const std::vector<real_type> & ys = rx_corr_pct_[el_id];
            const std::size_t n = xs.size();
            if(n == 0) return my_zero_;
            if(a <= xs.front()) return ys.front();
            if(a >= xs.back()) return ys.back();
            std::size_t hi = 1;
            while(hi < n && xs[hi] < a) ++hi;
            const real_type t = (a - xs[hi - 1]) / (xs[hi] - xs[hi - 1]);
            return ys[hi - 1] + t * (ys[hi] - ys[hi - 1]);
        }

    protected:

        real_type _ptdf_x(int tr_id) const override {
            real_type res = x_(tr_id);
            real_type tau = is_tap_side1_[tr_id] ? ratio_(tr_id) : 1. / ratio_(tr_id);
            return res * tau;
        }

        int _ptdf_row(int tr_id, int nb_powerline) const override {
            return tr_id + nb_powerline;
        }

        FDPFCoeffs _fdpf_coeffs(int tr_id, FDPFMethod xb_or_bx) const override;
};

inline TrafoInfo::TrafoInfo(const TrafoContainer & r_data_trafo, int my_id) noexcept:
BranchInfo(r_data_trafo, my_id),
ratio(-1.0),
shift_rad(-1.0),
is_tap_side1(true),
has_ratio_tap_changer(false),
ratio_tap_position(0),
ratio_low_tap(0),
ratio_high_tap(0),
ratio_regulation_mode(RegulationMode::FIXED),
ratio_regulating(false),
ratio_target(0.),
ratio_deadband(0.),
ratio_regulated(-1),
has_phase_tap_changer(false),
phase_tap_position(0),
phase_low_tap(0),
phase_high_tap(0),
phase_regulation_mode(RegulationMode::FIXED),
phase_regulating(false),
phase_target(0.),
phase_deadband(0.),
phase_regulated(-1),
res_phase_tap_position(0)
{
    if(my_id < 0) return;
    if(my_id >= r_data_trafo.nb()) return;
    is_tap_side1 = r_data_trafo.is_tap_side1_[my_id];
    ratio = r_data_trafo.ratio_.coeff(my_id);
    shift_rad = r_data_trafo.shift_.coeff(my_id);
    const TapChangers & rtc = r_data_trafo.ratio_taps_;
    if(rtc.nb() > my_id && rtc.has(my_id)){
        has_ratio_tap_changer = true;
        ratio_tap_position = rtc.position(my_id);
        ratio_low_tap = rtc.low_tap(my_id);
        ratio_high_tap = rtc.high_tap(my_id);
    }
    if(rtc.nb() > my_id){
        ratio_regulation_mode = rtc.mode(my_id);
        ratio_regulating = rtc.regulating(my_id);
        ratio_target = rtc.target(my_id);
        ratio_deadband = rtc.deadband(my_id);
        ratio_regulated = rtc.regulated(my_id);
    }
    const TapChangers & ptc = r_data_trafo.phase_taps_;
    if(ptc.nb() > my_id && ptc.has(my_id)){
        has_phase_tap_changer = true;
        phase_tap_position = ptc.position(my_id);
        phase_low_tap = ptc.low_tap(my_id);
        phase_high_tap = ptc.high_tap(my_id);
    }
    if(ptc.nb() > my_id){
        phase_regulation_mode = ptc.mode(my_id);
        phase_regulating = ptc.regulating(my_id);
        phase_target = ptc.target(my_id);
        phase_deadband = ptc.deadband(my_id);
        phase_regulated = ptc.regulated(my_id);
    }
    const std::vector<int> & res_pos = r_data_trafo.get_res_phase_tap_position();
    res_phase_tap_position = static_cast<std::size_t>(my_id) < res_pos.size() ? res_pos[static_cast<std::size_t>(my_id)]
                                                                              : phase_tap_position;
}


} // namespace ls2g

#endif  //TRAFO_CONTAINER_H
