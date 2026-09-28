// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LIGHT_ENV_H
#define LIGHT_ENV_H

#include "LSGrid.hpp"
#include "CustTimer.hpp"

#include "topo_action.hpp"
#include "inj_action.hpp"
#include "protections.hpp"

#include <type_traits>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace ls2g {

class LightEnv;

/**
 * The observation of a `LightEnv`: a read-only view on the environment's current state, it
 * holds nothing but a pointer to its environment and copies nothing. Every getter returns an
 * `Eigen::Ref` on memory the environment (or its grid) owns, so the values follow the
 * environment as it steps: this is not a snapshot, copy what you want to keep.
 *
 * A view is valid until the next `reset()` of its environment, which rebuilds the grid and the
 * protections, or until a step that ends the episode by a divergence (the grid then drops its
 * results). Once a step returned `done`, the values are not meaningful.
 *
 * Lines use grid2op numbering (powerlines then transformers), the "or" side of a transformer
 * being its hv side. Powers are in MW / MVAr, currents in kA (the unit of the thermal limits
 * of `Protections`), `topo_vect` holds local busbar ids (-1 for a disconnected element).
 */
class LightEnvObservation
{
    public:
        explicit LightEnvObservation(const LightEnv & env): env_(&env) {}

        Eigen::Ref<const RealVect> get_rho() const;

        Eigen::Ref<const RealVect> get_p_or() const;
        Eigen::Ref<const RealVect> get_q_or() const;
        Eigen::Ref<const RealVect> get_a_or() const;
        Eigen::Ref<const RealVect> get_p_ex() const;
        Eigen::Ref<const RealVect> get_q_ex() const;
        Eigen::Ref<const RealVect> get_a_ex() const;

        Eigen::Ref<const RealVect> get_load_p() const;
        Eigen::Ref<const RealVect> get_gen_p() const;

        Eigen::Ref<const IntVect> get_topo_vect() const;
        Eigen::Ref<const IntVect> get_time_before_cooldown_line() const;
        Eigen::Ref<const IntVect> get_time_before_cooldown_sub() const;

        int get_current_step() const;

    private:
        const LightEnv * env_;
};

/**
 * A `std::unique_ptr` whose copy is a deep copy (`LSGrid` has a copy constructor but no copy
 * assignment, so the live grid of a `LightEnv` is held through a pointer).
 */
template<class T>
class DeepCopyPtr
{
    public:
        explicit DeepCopyPtr(T * ptr): ptr_(ptr) {}
        DeepCopyPtr(const DeepCopyPtr & other): ptr_(other.ptr_ ? new T(*other.ptr_) : nullptr) {}
        DeepCopyPtr & operator=(const DeepCopyPtr & other){
            if(this != &other) ptr_.reset(other.ptr_ ? new T(*other.ptr_) : nullptr);
            return *this;
        }
        DeepCopyPtr(DeepCopyPtr &&) noexcept = default;
        DeepCopyPtr & operator=(DeepCopyPtr &&) noexcept = default;

        T & operator*() const {return *ptr_;}
        T * operator->() const {return ptr_.get();}
        explicit operator bool() const noexcept {return static_cast<bool>(ptr_);}
        void reset(T * ptr) {ptr_.reset(ptr);}

    private:
        std::unique_ptr<T> ptr_;
};

/**
 * Everything a `LightEnv` holds but its observation. Kept apart so that it is copied member
 * by member, while `LightEnv` re-points its observation to itself: adding a member here
 * cannot be forgotten by the copy.
 *
 * What never changes during an episode (the initial grid, the time series, the registered
 * actions) is shared between copies, read-only: `assign_time_series` / `init_actions` on
 * one env replace its own pointer and leave the others alone. The rest (the live grid, the
 * protections, the cooldowns, the step...) is copied, so a copy is an independent env at the
 * same point of the same episode.
 */
class LightEnvState
{
    protected:
        typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RealMat;

        // one row per step, one column per element, NaN meaning "unchanged"
        struct TimeSeries
        {
            RealMat load_p;
            RealMat load_q;
            RealMat gen_p;
            RealMat gen_v;
            RealMat storage_p;
            RealMat shunt_p;
            RealMat shunt_q;
            RealMat sgen_p;
            RealMat sgen_q;
        };

        explicit LightEnvState(const LSGrid & gridmodel):
            has_been_checked_(false),
            max_step_(0),
            actions_(std::make_shared<const std::vector<TopoAction> >()),
            init_grid_(std::make_shared<const LSGrid>(gridmodel)),
            grid_(new LSGrid(gridmodel)),
            max_iter_(10),
            tol_(1e-8),
            current_step_(0),
            nb_timestep_cooldown_sub_(0),
            nb_timestep_cooldown_line_(0),
            nb_timestep_reconnection_(10),
            timer_step_(0.),
            timer_reset_(0.),
            timer_obs_(0.),
            timer_update_gridmodel_(0.)
            {}

        // consistency
        bool has_been_checked_;

        // time series (shared between copies), nullptr until assign_time_series
        std::shared_ptr<const TimeSeries> time_series_;
        int max_step_;  // size of the input data

        // protections and thermal limits
        Protections protections_;

        // the actions the agent can take (shared between copies)
        std::shared_ptr<const std::vector<TopoAction> > actions_;

        // powergrid state: the initial one (restored at reset, shared between copies) and the
        // live one
        std::shared_ptr<const LSGrid> init_grid_;
        DeepCopyPtr<LSGrid> grid_;
        CplxVect V_;
        int max_iter_;
        real_type tol_;

        // time related information
        int current_step_;
        int nb_timestep_cooldown_sub_;
        int nb_timestep_cooldown_line_;
        int nb_timestep_reconnection_;
        Eigen::VectorXi time_step_sub_cooldown_;
        Eigen::VectorXi time_step_line_cooldown_;

        // the buffers the observation views (grid2op order, see LightEnv::extract_observation)
        RealVect p_or_;
        RealVect q_or_;
        RealVect a_or_;
        RealVect p_ex_;
        RealVect q_ex_;
        RealVect a_ex_;
        IntVect topo_vect_;

        // timers
        double timer_step_;
        double timer_reset_;
        double timer_obs_;
        double timer_update_gridmodel_;
};

// what makes the moves of LightEnv noexcept (LSGrid itself is not nothrow-movable, it is
// held through DeepCopyPtr): a member added to LightEnvState must keep these true, on every
// standard library -- MSVC's std::unordered_map / std::map / std::list / std::deque, for
// instance, do not have a noexcept move constructor (their sentinel node is allocated)
static_assert(std::is_nothrow_move_constructible<LightEnvState>::value,
              "LightEnvState should be nothrow move constructible");
static_assert(std::is_nothrow_move_assignable<LightEnvState>::value,
              "LightEnvState should be nothrow move assignable");

/**
 * A (very) limited grid2op environment in pure c++.
 *
 * The grid given at construction is the initial state of every episode: `reset()` restores
 * its topology. Injections are replayed from the time series given to `assign_time_series`
 * (one row per step, NaN meaning "unchanged"). The agent acts through `step(act_id)`, where
 * `act_id` indexes the actions given to `init_actions` (a `TopoAction` each). Without any
 * action initialised, `act_id = 0` is the only valid id and does nothing.
 *
 * Cooldowns follow grid2op: after an action affects a substation (resp. a line) that
 * substation (resp. line) cannot be acted on for `nb_timestep_cooldown_sub` (resp.
 * `nb_timestep_cooldown_line`) steps, and a line disconnected by the protections cannot be
 * reconnected for `nb_timestep_reconnection` steps. An action that touches an element still
 * in cooldown is illegal: it is replaced by "do nothing" and `info["is_illegal"]` is "true".
 *
 * A `LightEnv` can be copied: the copy is an independent env at the same point of the same
 * episode, see `LightEnvState` for what is shared and what is copied. Its observation views the
 * copy, not the original.
 *
 * It can be moved too (noexcept, so a `std::vector<LightEnv>` moves rather than copies when it
 * grows): nothing is copied, the observation views the env moved to. A moved-from env can only
 * be destroyed or assigned to; anything else (reset, step, the grid, load_p / gen_p of its
 * observation) throws `std::logic_error`.
 */
class LightEnv : protected LightEnvState
{
    public:
        typedef std::unordered_map<std::string, std::string> InfoReturnedType;
        typedef std::tuple<const LightEnvObservation &, double, bool, bool, InfoReturnedType > StepReturnedType;
        typedef std::tuple<const LightEnvObservation &, InfoReturnedType > ResetReturnedType;

        explicit LightEnv(const LSGrid & gridmodel):
            LightEnvState(gridmodel),
            observation_(*this)
            {
                // every buffer an observation views is allocated once, here, and only ever
                // written in place: a numpy view on it stays valid as the env steps
                const int nb_line = static_cast<int>(grid_->nb_powerline() + grid_->nb_trafo());
                p_or_ = RealVect::Zero(nb_line);
                q_or_ = RealVect::Zero(nb_line);
                a_or_ = RealVect::Zero(nb_line);
                p_ex_ = RealVect::Zero(nb_line);
                q_ex_ = RealVect::Zero(nb_line);
                a_ex_ = RealVect::Zero(nb_line);
                time_step_sub_cooldown_ = Eigen::VectorXi::Zero(grid_->get_n_sub());
                time_step_line_cooldown_ = Eigen::VectorXi::Zero(nb_line);
                topo_vect_ = IntVect::Constant(aux_dim_topo(*grid_), BaseConstants::_deactivated_bus_id);
            }

        // the state is copied (or moved) member by member, the observation keeps viewing
        // the env it belongs to
        LightEnv(const LightEnv & other): LightEnvState(other), observation_(*this) {}
        LightEnv(LightEnv && other) noexcept: LightEnvState(std::move(other)), observation_(*this) {
            other.has_been_checked_ = false;
        }
        LightEnv & operator=(const LightEnv & other){
            LightEnvState::operator=(other);
            return *this;
        }
        LightEnv & operator=(LightEnv && other) noexcept {
            if(this != &other){
                LightEnvState::operator=(std::move(other));
                other.has_been_checked_ = false;
            }
            return *this;
        }

        void assign_time_series(const Eigen::Ref<const RealMat> & load_p,
                                const Eigen::Ref<const RealMat> & load_q,
                                const Eigen::Ref<const RealMat> & gen_p,
                                const Eigen::Ref<const RealMat> & gen_v,
                                const Eigen::Ref<const RealMat> & storage_p,
                                const Eigen::Ref<const RealMat> & shunt_p,
                                const Eigen::Ref<const RealMat> & shunt_q,
                                const Eigen::Ref<const RealMat> & sgen_p,
                                const Eigen::Ref<const RealMat> & sgen_q){
            has_been_checked_ = false;
            time_series_ = std::make_shared<const TimeSeries>(
                TimeSeries{load_p, load_q, gen_p, gen_v, storage_p, shunt_p, shunt_q, sgen_p, sgen_q});
        }

        void assign_protections(const Protections & protections){
            has_been_checked_ = false;
            protections_ = protections;
        }

        /**
         * Register the actions the agent can take: `step(i)` will play `actions[i]`.
         *
         * Every action is checked against the initial grid (see `TopoAction::check_validity`);
         * if one is invalid, `std::invalid_argument` is thrown naming it and nothing is
         * registered (the previous actions, if any, are kept).
         */
        void init_actions(const std::vector<TopoAction> & actions){
            std::vector<TopoAction> checked;
            checked.reserve(actions.size());
            for(size_t i = 0; i < actions.size(); ++i){
                TopoAction act = actions[i];
                try{
                    act.check_validity(*init_grid_);
                }catch(const std::exception & exc_){
                    std::ostringstream msg;
                    msg << "LightEnv::init_actions: action " << i << " is invalid: " << exc_.what();
                    throw std::invalid_argument(msg.str());
                }
                checked.push_back(act);
            }
            actions_ = std::make_shared<const std::vector<TopoAction> >(std::move(checked));
        }
        int nb_actions() const {return static_cast<int>(actions_->size());}
        const std::vector<TopoAction> & get_actions() const {return *actions_;}

        const Protections & get_protections() const {return protections_;}
        const LSGrid & get_grid() const {aux_check_not_moved_from("get_grid"); return *grid_;}
        double get_step_time() const {return timer_step_;}
        double get_reset_time() const {return timer_reset_;}
        double get_obs_time() const {return timer_obs_;}
        double get_update_gridmodel_time() const {return timer_update_gridmodel_;}

        int get_max_iter() const {return max_iter_;}
        void set_max_iter(int max_iter) {max_iter_ = max_iter;}

        real_type get_tol() const {return tol_;}
        void set_tol(real_type tol) {tol_ = tol;}

        // cooldown parameters (grid2op NB_TIMESTEP_COOLDOWN_SUB / _LINE / NB_TIMESTEP_RECONNECTION)
        int get_nb_timestep_cooldown_sub() const {return nb_timestep_cooldown_sub_;}
        void set_nb_timestep_cooldown_sub(int val) {aux_check_positive(val, "nb_timestep_cooldown_sub"); nb_timestep_cooldown_sub_ = val;}
        int get_nb_timestep_cooldown_line() const {return nb_timestep_cooldown_line_;}
        void set_nb_timestep_cooldown_line(int val) {aux_check_positive(val, "nb_timestep_cooldown_line"); nb_timestep_cooldown_line_ = val;}
        int get_nb_timestep_reconnection() const {return nb_timestep_reconnection_;}
        void set_nb_timestep_reconnection(int val) {aux_check_positive(val, "nb_timestep_reconnection"); nb_timestep_reconnection_ = val;}

        // current cooldowns (grid2op obs.time_before_cooldown_sub / _line)
        Eigen::Ref<const Eigen::VectorXi> get_time_before_cooldown_sub() const {return time_step_sub_cooldown_;}
        Eigen::Ref<const Eigen::VectorXi> get_time_before_cooldown_line() const {return time_step_line_cooldown_;}

        const LightEnvObservation & get_obs() const {return observation_;}
        int get_max_step() const {return max_step_;}
        int get_current_step() const {return current_step_;}

        ResetReturnedType reset(){
            aux_check_not_moved_from("reset");
            auto timer_reset = CustTimer();
            reset_timers();

            // back to the initial topology
            grid_.reset(new LSGrid(*init_grid_));

            if(!has_been_checked_){
                perform_internal_checks();
                has_been_checked_ = true;
            }else{
                protections_.reset_episode();
            }

            current_step_ = 0;

            // reset the cooldowns
            reset_cooldowns();

            // set the correct injection to the grid
            apply_injections(current_step_);

            // run the powerflow
            V_ = protections_.run_powerflow(*grid_, max_iter_, tol_);
            // TODO if divergence !!!!
            if(V_.size() != 0) protections_.update_rho(*grid_);

            // extract the observation
            extract_observation();

            timer_reset_ += timer_reset.duration();
            return ResetReturnedType(observation_, InfoReturnedType());
        }

        StepReturnedType step(int act_id){
            aux_check_not_moved_from("step");
            auto timer_step = CustTimer();
            if(!has_been_checked_){
                throw std::runtime_error("Environment cannot be used, you most likely need to call env.reset() ");
            }
            // fails (leaving the env untouched) on an unknown action id
            const TopoAction * action = aux_get_action(act_id);

            current_step_ += 1;
            InfoReturnedType info;
            info["is_illegal"] = "false";
            if (current_step_ >= max_step_){
                info["success"] = "true";
                info["failure"] = "false";
                info["survival_time"] = std::to_string(survival_ratio());
                has_been_checked_ = false;
                timer_step_ += timer_step.duration();
                return StepReturnedType(observation_, 1., true, false, info);
            }

            // apply the topology, if legal (cooldowns)
            bool action_applied = false;
            std::vector<bool> subs_impacted;
            std::vector<bool> lines_impacted;
            if(action != nullptr && !action->is_do_nothing()){
                action->compute_impact(*grid_, subs_impacted, lines_impacted);
                const bool is_illegal = aux_is_illegal(subs_impacted, lines_impacted);
                if(is_illegal){
                    info["is_illegal"] = "true";
                }else{
                    action->apply_to_gridmodel(*grid_);
                    action_applied = true;
                }
            }

            // apply the injection
            auto timer_update_gridmodel = CustTimer();
            apply_injections(current_step_);
            timer_update_gridmodel_ += timer_update_gridmodel.duration();
            // TODO apply redispatching, storage, curtailment and other !

            V_ = protections_.next_grid_state(*grid_, max_iter_, tol_);
            if(V_.size() == 0){
                // divergence
                info["success"] = "false";
                info["failure"] = "true";
                info["survival_time"] = std::to_string(survival_ratio());
                has_been_checked_ = false;
                timer_step_ += timer_step.duration();
                return StepReturnedType(observation_, 0., true, true, info);
            }

            // extract the observation
            extract_observation();

            // update time dependant information
            update_cooldowns(action_applied, subs_impacted, lines_impacted);

            timer_step_ += timer_step.duration();
            return StepReturnedType(observation_, survival_ratio(), false, false, info);
        }

    protected:
        // the live grid is only ever null in a moved-from env
        void aux_check_not_moved_from(const char * where) const {
            if(!grid_){
                std::ostringstream exc_;
                exc_ << "LightEnv::" << where << ": this env has been moved from, it can only be destroyed or "
                     << "assigned to.";
                throw std::logic_error(exc_.str());
            }
        }

        const TopoAction * aux_get_action(int act_id) const {
            const std::vector<TopoAction> & actions = *actions_;
            if(actions.empty()){
                if(act_id != 0){
                    std::ostringstream exc_;
                    exc_ << "LightEnv::step: no action has been initialised (see init_actions), "
                         << "only act_id = 0 (do nothing) is valid, you provided " << act_id << ".";
                    throw std::out_of_range(exc_.str());
                }
                return nullptr;
            }
            if(act_id < 0 || act_id >= static_cast<int>(actions.size())){
                std::ostringstream exc_;
                exc_ << "LightEnv::step: unknown action id " << act_id << ", "
                     << actions.size() << " actions have been initialised (valid ids: 0 to "
                     << actions.size() - 1 << ").";
                throw std::out_of_range(exc_.str());
            }
            return &actions[act_id];
        }

        bool aux_is_illegal(const std::vector<bool> & subs_impacted,
                            const std::vector<bool> & lines_impacted) const {
            for(size_t sub_id = 0; sub_id < subs_impacted.size(); ++sub_id){
                if(subs_impacted[sub_id] && time_step_sub_cooldown_(sub_id) > 0) return true;
            }
            for(size_t line_id = 0; line_id < lines_impacted.size(); ++line_id){
                if(lines_impacted[line_id] && time_step_line_cooldown_(line_id) > 0) return true;
            }
            return false;
        }

        /**
         * Same order as grid2op (`BaseEnv._aux_register_env_converged`): one step passed for
         * the lines, the lines the protections disconnected get the reconnection cooldown,
         * the lines the action touched get the line cooldown (unless already longer), then
         * one step passed for the substations and the ones the action touched get the
         * substation cooldown.
         */
        void update_cooldowns(bool action_applied,
                              const std::vector<bool> & subs_impacted,
                              const std::vector<bool> & lines_impacted){
            for(int i = 0; i < time_step_line_cooldown_.size(); ++i){
                if(time_step_line_cooldown_(i) > 0) time_step_line_cooldown_(i) -= 1;
            }
            for(int line_id : protections_.lines_disconnected_this_step()){
                time_step_line_cooldown_(line_id) = nb_timestep_reconnection_;
            }
            if(action_applied && nb_timestep_cooldown_line_ > 0){
                for(size_t i = 0; i < lines_impacted.size(); ++i){
                    if(lines_impacted[i] && time_step_line_cooldown_(i) < nb_timestep_cooldown_line_){
                        time_step_line_cooldown_(i) = nb_timestep_cooldown_line_;
                    }
                }
            }
            for(int i = 0; i < time_step_sub_cooldown_.size(); ++i){
                if(time_step_sub_cooldown_(i) > 0) time_step_sub_cooldown_(i) -= 1;
            }
            if(action_applied && nb_timestep_cooldown_sub_ > 0){
                for(size_t i = 0; i < subs_impacted.size(); ++i){
                    if(subs_impacted[i]) time_step_sub_cooldown_(i) = nb_timestep_cooldown_sub_;
                }
            }
        }

        void apply_injections(int step_id){
            const TimeSeries & ts = *time_series_;
            const InjAction inj_action(ts.load_p.row(step_id),
                                       ts.load_q.row(step_id),
                                       ts.gen_p.row(step_id),
                                       ts.gen_v.row(step_id),
                                       ts.storage_p.row(step_id),
                                       ts.shunt_p.row(step_id),
                                       ts.shunt_q.row(step_id),
                                       ts.sgen_p.row(step_id),
                                       ts.sgen_q.row(step_id)
                                       );
            inj_action.apply_to_gridmodel(*grid_);
        }

        // fraction of the episode survived, in [0, 1]: the reward of a step and the
        // "survival_time" of the info at the end of an episode
        double survival_ratio() const {
            if(max_step_ <= 0) return 1.;
            return static_cast<double>(current_step_) / static_cast<double>(max_step_);
        }

        // what an observation reads that is not stored in grid2op order by the grid: the flows
        // (lines then trafos) and the topology vector, written in place in buffers allocated
        // by the constructor. rho, load_p, gen_p and the cooldowns are read where they live.
        void extract_observation(){
            auto timer_obs = CustTimer();
            const LSGrid & grid = *grid_;
            const int nb_line = static_cast<int>(grid.nb_powerline());
            const int nb_trafo = static_cast<int>(grid.nb_trafo());
            const tuple4d line_or = grid.get_line_res1();
            const tuple4d line_ex = grid.get_line_res2();
            const tuple4d trafo_or = grid.get_trafo_res1();
            const tuple4d trafo_ex = grid.get_trafo_res2();
            p_or_.head(nb_line) = std::get<0>(line_or);
            q_or_.head(nb_line) = std::get<1>(line_or);
            a_or_.head(nb_line) = std::get<3>(line_or);
            p_ex_.head(nb_line) = std::get<0>(line_ex);
            q_ex_.head(nb_line) = std::get<1>(line_ex);
            a_ex_.head(nb_line) = std::get<3>(line_ex);
            p_or_.tail(nb_trafo) = std::get<0>(trafo_or);
            q_or_.tail(nb_trafo) = std::get<1>(trafo_or);
            a_or_.tail(nb_trafo) = std::get<3>(trafo_or);
            p_ex_.tail(nb_trafo) = std::get<0>(trafo_ex);
            q_ex_.tail(nb_trafo) = std::get<1>(trafo_ex);
            a_ex_.tail(nb_trafo) = std::get<3>(trafo_ex);

            if(topo_vect_.size() > 0){
                const SubstationContainer & subs = grid.get_substations();
                aux_fill_topo_vect(subs, grid.get_loads().get_pos_topo_vect(),
                                   [&grid](int el_id){return grid.get_loads().get_bus(el_id);});
                aux_fill_topo_vect(subs, grid.get_generators().get_pos_topo_vect(),
                                   [&grid](int el_id){return grid.get_generators().get_bus(el_id);});
                aux_fill_topo_vect(subs, grid.get_storages().get_pos_topo_vect(),
                                   [&grid](int el_id){return grid.get_storages().get_bus(el_id);});
                aux_fill_topo_vect(subs, grid.get_lines().get_pos_topo_vect_side_1(),
                                   [&grid](int el_id){return grid.get_lines().get_bus_side_1(el_id);});
                aux_fill_topo_vect(subs, grid.get_lines().get_pos_topo_vect_side_2(),
                                   [&grid](int el_id){return grid.get_lines().get_bus_side_2(el_id);});
                aux_fill_topo_vect(subs, grid.get_trafos().get_pos_topo_vect_side_1(),
                                   [&grid](int el_id){return grid.get_trafos().get_bus_side_1(el_id);});
                aux_fill_topo_vect(subs, grid.get_trafos().get_pos_topo_vect_side_2(),
                                   [&grid](int el_id){return grid.get_trafos().get_bus_side_2(el_id);});
            }
            timer_obs_ += timer_obs.duration();
        }

        template<class BusGetter>
        void aux_fill_topo_vect(const SubstationContainer & subs,
                                const Eigen::Ref<const IntVect> & pos_topo_vect,
                                BusGetter get_bus){
            for(int el_id = 0; el_id < pos_topo_vect.size(); ++el_id){
                topo_vect_(pos_topo_vect(el_id)) = subs.gridmodel_to_local(get_bus(el_id)).cast_int();
            }
        }

        // size of the grid2op topology vector, 0 if the grid does not carry the positions
        // (they are set by LightSimBackend, not by the grid converters)
        static int aux_dim_topo(const LSGrid & grid){
            const IntVect * positions[] = {
                &grid.get_loads().get_pos_topo_vect(),
                &grid.get_generators().get_pos_topo_vect(),
                &grid.get_storages().get_pos_topo_vect(),
                &grid.get_lines().get_pos_topo_vect_side_1(),
                &grid.get_lines().get_pos_topo_vect_side_2(),
                &grid.get_trafos().get_pos_topo_vect_side_1(),
                &grid.get_trafos().get_pos_topo_vect_side_2()};
            const int nb_els[] = {
                grid.get_loads().nb(),
                grid.get_generators().nb(),
                grid.get_storages().nb(),
                static_cast<int>(grid.nb_powerline()),
                static_cast<int>(grid.nb_powerline()),
                static_cast<int>(grid.nb_trafo()),
                static_cast<int>(grid.nb_trafo())};
            int dim_topo = 0;
            for(int i = 0; i < 7; ++i){
                if(positions[i]->size() != nb_els[i]) return 0;
                dim_topo += nb_els[i];
            }
            return dim_topo;
        }

        void reset_cooldowns(){
            time_step_sub_cooldown_.setZero();
            time_step_line_cooldown_.setZero();
        }

        void reset_timers(){
            timer_step_ = 0.;
            timer_reset_ = 0.;
            timer_obs_ = 0.;
            timer_update_gridmodel_ = 0.;
        }

        static void aux_check_positive(int val, const std::string & name){
            if(val < 0){
                std::ostringstream exc_;
                exc_ << "LightEnv: " << name << " should be >= 0, you provided " << val << ".";
                throw std::invalid_argument(exc_.str());
            }
        }

        void perform_internal_checks(){
            // correct number of rows (steps)
            if(!time_series_){
                throw std::runtime_error("LightEnv: no time series, call assign_time_series before reset.");
            }
            const TimeSeries & ts = *time_series_;
            const int n_ts = static_cast<int>(ts.load_p.rows());
            aux_check_row(ts.load_p, n_ts, "perform_internal_checks (load_p)");
            aux_check_row(ts.load_q, n_ts, "perform_internal_checks (load_q)");
            aux_check_row(ts.gen_p, n_ts, "perform_internal_checks (gen_p)");
            aux_check_row(ts.gen_v, n_ts, "perform_internal_checks (gen_v)");
            aux_check_row(ts.storage_p, n_ts, "perform_internal_checks (storage_p)");
            aux_check_row(ts.shunt_p, n_ts, "perform_internal_checks (shunt_p)");
            aux_check_row(ts.shunt_q, n_ts, "perform_internal_checks (shunt_q)");
            aux_check_row(ts.sgen_p, n_ts, "perform_internal_checks (sgen_p)");
            aux_check_row(ts.sgen_q, n_ts, "perform_internal_checks (sgen_q)");
            max_step_ = n_ts;

            // correct number of columns (number of elelments)
            aux_check_col(ts.load_p, grid_->get_loads().nb(), "perform_internal_checks (load_p)");
            aux_check_col(ts.load_q, grid_->get_loads().nb(), "perform_internal_checks (load_q)");
            aux_check_col(ts.gen_p, grid_->get_generators().nb(), "perform_internal_checks (gen_p)");
            aux_check_col(ts.gen_v, grid_->get_generators().nb(), "perform_internal_checks (gen_v)");
            aux_check_col(ts.storage_p, grid_->get_storages().nb(), "perform_internal_checks (storage_p)");
            aux_check_col(ts.shunt_p, grid_->get_shunts().nb(), "perform_internal_checks (shunt_p)");
            aux_check_col(ts.shunt_q, grid_->get_shunts().nb(), "perform_internal_checks (shunt_q)");
            aux_check_col(ts.sgen_p, grid_->get_static_generators().nb(), "perform_internal_checks (sgen_p)");
            aux_check_col(ts.sgen_q, grid_->get_static_generators().nb(), "perform_internal_checks (sgen_q)");

            // protections
            protections_.check_validity(*grid_);

            // cooldowns
            reset_cooldowns();
        }

        template<class EigenType>
        void aux_check_row(const EigenType & new_vect,
                           int size_th,
                           const std::string & error_detailed) const{
            if(new_vect.rows() != size_th)
            {
                std::ostringstream exc_;
                exc_ << "LightEnv::" << error_detailed << ".";
                exc_ << "Theorical size is: ";
                exc_ << size_th;
                exc_ << " but your provided data with ";
                exc_ << new_vect.rows() << " rows.";
                throw std::runtime_error(exc_.str());
            }
        }
        template<class EigenType>
        void aux_check_col(const EigenType & new_vect,
                           int size_th,
                           const std::string & error_detailed) const{
            if(new_vect.cols() != size_th)
            {
                std::ostringstream exc_;
                exc_ << "LightEnv::" << error_detailed << ".";
                exc_ << "Theorical size is: ";
                exc_ << size_th;
                exc_ << " but your provided data with ";
                exc_ << new_vect.cols() << " columns.";
                throw std::runtime_error(exc_.str());
            }
        }

    protected:
        LightEnvObservation observation_;

        friend class LightEnvObservation;
};

inline Eigen::Ref<const RealVect> LightEnvObservation::get_rho() const {return env_->protections_.get_rho();}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_p_or() const {return env_->p_or_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_q_or() const {return env_->q_or_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_a_or() const {return env_->a_or_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_p_ex() const {return env_->p_ex_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_q_ex() const {return env_->q_ex_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_a_ex() const {return env_->a_ex_;}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_load_p() const {
    env_->aux_check_not_moved_from("get_obs().get_load_p");
    return std::get<0>(env_->grid_->get_loads_res());
}
inline Eigen::Ref<const RealVect> LightEnvObservation::get_gen_p() const {
    env_->aux_check_not_moved_from("get_obs().get_gen_p");
    return std::get<0>(env_->grid_->get_gen_res());
}
inline Eigen::Ref<const IntVect> LightEnvObservation::get_topo_vect() const {
    if(env_->topo_vect_.size() == 0){
        throw std::runtime_error("LightEnvObservation::get_topo_vect: the grid of this environment has no position "
                                 "in the grid2op topology vector (set_*_pos_topo_vect, done by LightSimBackend).");
    }
    return env_->topo_vect_;
}
inline Eigen::Ref<const IntVect> LightEnvObservation::get_time_before_cooldown_line() const {return env_->time_step_line_cooldown_;}
inline Eigen::Ref<const IntVect> LightEnvObservation::get_time_before_cooldown_sub() const {return env_->time_step_sub_cooldown_;}
inline int LightEnvObservation::get_current_step() const {return env_->current_step_;}

} // namespace ls2g

#endif // LIGHT_ENV_H
