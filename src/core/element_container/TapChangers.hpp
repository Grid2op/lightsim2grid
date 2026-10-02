// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TAP_CHANGERS_H
#define TAP_CHANGERS_H

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Utils.hpp"

namespace ls2g {

/**
 * What a discrete control (a transformer's tap changer, a shunt's sections) regulates. The
 * unit of its target and deadband follows: a voltage in pu of the regulated bus' nominal
 * voltage, a reactive or active power in MVar / MW, a current in A.
 */
enum class LS2G_API RegulationMode : int {
    FIXED = 0,            ///< nothing: the position only moves when someone moves it
    VOLTAGE = 1,          ///< the voltage of a bus (`regulated` the grid bus id)
    REACTIVE_POWER = 2,   ///< the reactive power through one side of the branch (`regulated` the side)
    CURRENT_LIMITER = 3,  ///< the current through one side of the branch, kept below the target
    ACTIVE_POWER = 4      ///< the active power through one side of the branch (`regulated` the side)
};

/**
 * The tap changers of ONE kind (ratio or phase) of every transformer of a container: per
 * transformer, an optional step table and its current position, and what the changer
 * regulates. A transformer without a table has none of that kind and contributes a neutral
 * step (rho 1, alpha 0, no impedance correction).
 *
 * A step is what IIDM (and OpenLoadFlow's Transformers.getTapCharacteristics) gives: a
 * ratio factor rho, a phase shift alpha (rad, 0 for a ratio changer) and corrections, in %,
 * of the transformer's r, x, g and b. The position runs from `low_tap` to
 * `low_tap + n_steps - 1`, like IIDM's.
 */
class LS2G_API TapChangers
{
    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
                   std::vector<int>,  // low_tap
                   std::vector<int>,  // position
                   std::vector<std::vector<real_type> >,  // rho
                   std::vector<std::vector<real_type> >,  // alpha (rad)
                   std::vector<std::vector<real_type> >,  // r (%)
                   std::vector<std::vector<real_type> >,  // x (%)
                   std::vector<std::vector<real_type> >,  // g (%)
                   std::vector<std::vector<real_type> >,  // b (%)
                   std::vector<int>,  // regulation mode (RegulationMode)
                   std::vector<bool>,  // regulating
                   std::vector<real_type>,  // target
                   std::vector<real_type>,  // deadband
                   std::vector<int>   // regulated: grid bus id (VOLTAGE) or side (otherwise), -1 none
               >;
        enum StateResIdx {
            LOW_TAP = 0,
            POSITION,
            RHO,
            ALPHA,
            R_PCT,
            X_PCT,
            G_PCT,
            B_PCT,
            MODE,
            REGULATING,
            TARGET,
            DEADBAND,
            REGULATED,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "TapChangers::StateRes and StateResIdx do not match");

        /// `nb` transformers, none with a tap changer
        void resize(int nb)
        {
            const std::size_t n = static_cast<std::size_t>(nb);
            low_tap_.assign(n, 0);
            position_.assign(n, 0);
            rho_.assign(n, std::vector<real_type>());
            alpha_.assign(n, std::vector<real_type>());
            r_pct_.assign(n, std::vector<real_type>());
            x_pct_.assign(n, std::vector<real_type>());
            g_pct_.assign(n, std::vector<real_type>());
            b_pct_.assign(n, std::vector<real_type>());
            mode_.assign(n, static_cast<int>(RegulationMode::FIXED));
            regulating_.assign(n, false);
            target_.assign(n, 0.);
            deadband_.assign(n, 0.);
            regulated_.assign(n, -1);
        }
        int nb() const { return static_cast<int>(low_tap_.size()); }
        bool has(int el) const { return !rho_[static_cast<std::size_t>(el)].empty(); }
        bool any() const
        {
            for (const auto & steps : rho_) if (!steps.empty()) return true;
            return false;
        }

        /**
         * The step table of transformer `el` (all vectors of the same, non-zero length; `alpha`
         * may be empty for a ratio changer, meaning 0 everywhere) and its current position.
         */
        void set(int el, int low_tap, int position,
                 const std::vector<real_type> & rho,
                 const std::vector<real_type> & alpha_rad,
                 const std::vector<real_type> & r_pct,
                 const std::vector<real_type> & x_pct,
                 const std::vector<real_type> & g_pct,
                 const std::vector<real_type> & b_pct,
                 const std::string & where)
        {
            _check_el(el, where);
            const std::size_t n = rho.size();
            if (n == 0 || r_pct.size() != n || x_pct.size() != n || g_pct.size() != n || b_pct.size() != n ||
                (!alpha_rad.empty() && alpha_rad.size() != n)) {
                std::ostringstream exc_;
                exc_ << where << ": the step table of transformer " << el << " has columns of different "
                     << "lengths (or none at all): every column has one value per position.";
                throw std::runtime_error(exc_.str());
            }
            const std::size_t k = static_cast<std::size_t>(el);
            low_tap_[k] = low_tap;
            rho_[k] = rho;
            alpha_[k] = alpha_rad.empty() ? std::vector<real_type>(n, 0.) : alpha_rad;
            r_pct_[k] = r_pct;
            x_pct_[k] = x_pct;
            g_pct_[k] = g_pct;
            b_pct_[k] = b_pct;
            set_position(el, position, where);
        }

        /// what the changer of transformer `el` regulates, see RegulationMode
        void set_regulation(int el, RegulationMode mode, bool regulating, real_type target,
                            real_type deadband, int regulated, const std::string & where)
        {
            _check_el(el, where);
            const std::size_t k = static_cast<std::size_t>(el);
            mode_[k] = static_cast<int>(mode);
            regulating_[k] = regulating;
            target_[k] = target;
            deadband_[k] = deadband;
            regulated_[k] = regulated;
        }

        int low_tap(int el) const { return low_tap_[static_cast<std::size_t>(el)]; }
        int high_tap(int el) const { return low_tap(el) + static_cast<int>(rho_[static_cast<std::size_t>(el)].size()) - 1; }
        int position(int el) const { return position_[static_cast<std::size_t>(el)]; }
        void set_position(int el, int position, const std::string & where)
        {
            _check_el(el, where);
            if (!has(el) || position < low_tap(el) || position > high_tap(el)) {
                std::ostringstream exc_;
                exc_ << where << ": position " << position << " is not one of transformer " << el << "'s";
                if (has(el)) exc_ << " (from " << low_tap(el) << " to " << high_tap(el) << ")";
                else exc_ << " (it has no tap changer of that kind)";
                throw std::runtime_error(exc_.str());
            }
            position_[static_cast<std::size_t>(el)] = position;
        }

        // the step at `position` (the current one by default); a transformer without a
        // changer of this kind gives the neutral step
        real_type rho_at(int el, int position) const { return _at(rho_, el, position, 1.); }
        real_type alpha_at(int el, int position) const { return _at(alpha_, el, position, 0.); }
        real_type r_factor_at(int el, int position) const { return 1. + _at(r_pct_, el, position, 0.) / 100.; }
        real_type x_factor_at(int el, int position) const { return 1. + _at(x_pct_, el, position, 0.) / 100.; }
        real_type g_factor_at(int el, int position) const { return 1. + _at(g_pct_, el, position, 0.) / 100.; }
        real_type b_factor_at(int el, int position) const { return 1. + _at(b_pct_, el, position, 0.) / 100.; }
        real_type rho(int el) const { return rho_at(el, position(el)); }
        real_type alpha(int el) const { return alpha_at(el, position(el)); }
        real_type r_factor(int el) const { return r_factor_at(el, position(el)); }
        real_type x_factor(int el) const { return x_factor_at(el, position(el)); }
        real_type g_factor(int el) const { return g_factor_at(el, position(el)); }
        real_type b_factor(int el) const { return b_factor_at(el, position(el)); }

        RegulationMode mode(int el) const { return static_cast<RegulationMode>(mode_[static_cast<std::size_t>(el)]); }
        bool regulating(int el) const { return regulating_[static_cast<std::size_t>(el)]; }
        real_type target(int el) const { return target_[static_cast<std::size_t>(el)]; }
        real_type deadband(int el) const { return deadband_[static_cast<std::size_t>(el)]; }
        int regulated(int el) const { return regulated_[static_cast<std::size_t>(el)]; }

        StateRes get_state() const
        {
            return StateRes(low_tap_, position_, rho_, alpha_, r_pct_, x_pct_, g_pct_, b_pct_,
                            mode_, regulating_, target_, deadband_, regulated_);
        }

        /// `nb` the number of transformers the state must describe
        void set_state(const StateRes & state, int nb, const std::string & where)
        {
            const std::size_t n = static_cast<std::size_t>(nb);
            const auto & low = std::get<LOW_TAP>(state);
            const auto & pos = std::get<POSITION>(state);
            const auto & rho = std::get<RHO>(state);
            const auto & alpha = std::get<ALPHA>(state);
            const auto & r = std::get<R_PCT>(state);
            const auto & x = std::get<X_PCT>(state);
            const auto & g = std::get<G_PCT>(state);
            const auto & b = std::get<B_PCT>(state);
            if (low.size() != n || pos.size() != n || rho.size() != n || alpha.size() != n || r.size() != n ||
                x.size() != n || g.size() != n || b.size() != n || std::get<MODE>(state).size() != n ||
                std::get<REGULATING>(state).size() != n || std::get<TARGET>(state).size() != n ||
                std::get<DEADBAND>(state).size() != n || std::get<REGULATED>(state).size() != n) {
                std::ostringstream exc_;
                exc_ << where << ": the tap changer state does not describe " << nb << " transformer(s).";
                throw std::runtime_error(exc_.str());
            }
            // the tables are read by position: each must be complete and the position in it
            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t steps = rho[k].size();
                const bool consistent = alpha[k].size() == steps && r[k].size() == steps && x[k].size() == steps &&
                                        g[k].size() == steps && b[k].size() == steps;
                const bool in_range = steps == 0 ||
                    (pos[k] >= low[k] && pos[k] - low[k] < static_cast<int>(steps));
                if (!consistent || !in_range) {
                    std::ostringstream exc_;
                    exc_ << where << ": the tap changer of transformer " << k << " is inconsistent (columns of "
                         << "different lengths, or a position outside its table).";
                    throw std::runtime_error(exc_.str());
                }
            }
            low_tap_ = low;
            position_ = pos;
            rho_ = rho;
            alpha_ = alpha;
            r_pct_ = r;
            x_pct_ = x;
            g_pct_ = g;
            b_pct_ = b;
            mode_ = std::get<MODE>(state);
            regulating_ = std::get<REGULATING>(state);
            target_ = std::get<TARGET>(state);
            deadband_ = std::get<DEADBAND>(state);
            regulated_ = std::get<REGULATED>(state);
        }

    private:
        void _check_el(int el, const std::string & where) const
        {
            if (el < 0 || el >= nb()) {
                std::ostringstream exc_;
                exc_ << where << ": no transformer with id " << el << " (there are " << nb() << ").";
                throw std::out_of_range(exc_.str());
            }
        }
        real_type _at(const std::vector<std::vector<real_type> > & column, int el, int position, real_type neutral) const
        {
            const std::vector<real_type> & steps = column[static_cast<std::size_t>(el)];
            if (steps.empty()) return neutral;
            return steps[static_cast<std::size_t>(position - low_tap(el))];
        }

        std::vector<int> low_tap_;
        std::vector<int> position_;
        std::vector<std::vector<real_type> > rho_;
        std::vector<std::vector<real_type> > alpha_;
        std::vector<std::vector<real_type> > r_pct_;
        std::vector<std::vector<real_type> > x_pct_;
        std::vector<std::vector<real_type> > g_pct_;
        std::vector<std::vector<real_type> > b_pct_;
        std::vector<int> mode_;
        std::vector<bool> regulating_;
        std::vector<real_type> target_;
        std::vector<real_type> deadband_;
        std::vector<int> regulated_;
};

}  // namespace ls2g

#endif  // TAP_CHANGERS_H
