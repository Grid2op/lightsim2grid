// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LSGRID_H
#define LSGRID_H

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <tuple>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI

#include "Utils.hpp"

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

// import data classes
#include "SubstationContainer.hpp"
#include "element_container/GenericContainer.hpp"
#include "element_container/LineContainer.hpp"
#include "element_container/ShuntContainer.hpp"
#include "element_container/TrafoContainer.hpp"
#include "element_container/LoadContainer.hpp"
#include "element_container/StorageContainer.hpp"
#include "element_container/GeneratorContainer.hpp"
#include "element_container/SGenContainer.hpp"
#include "element_container/SvcContainer.hpp"
#include "element_container/HvdcLineContainer.hpp"
#include "HvdcDroopData.hpp"
#include "VoltageControlData.hpp"

// import newton raphson solvers using different linear algebra solvers
#include "AlgorithmSelector.hpp"

namespace ls2g {

//TODO implement a BFS check to make sure the Ymatrix is "connected" [one single component]
class LS2G_API LSGrid final
{
    public:
        using IntVectRowMaj = Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor>;

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
                std::string, // version major
                std::string, // version medium
                std::string, // version minor
                std::vector<int>, // ls_to_orig
                real_type,  // init_vm_pu
                real_type, //sn_mva
                std::vector<bool>,  // bus_status
                SubstationContainer::StateRes,
                // powerlines
                LineContainer::StateRes ,
                // shunts
                ShuntContainer::StateRes,
                // trafos
                TrafoContainer::StateRes,
                // gens
                GeneratorContainer::StateRes,
                // loads
                LoadContainer::StateRes,
                // static generators
                SGenContainer::StateRes,
                // storage units
                StorageContainer::StateRes,
                //hvdc lines (was the "dc lines" before lightsim2grid 0.12)
                HvdcLineContainer::StateRes,
                // static var compensators (appended; old pickles are version-gated)
                SvcContainer::StateRes,
                // algo *registry names* (not AlgorithmType): the enum collapses
                // every external/plugin solver onto AlgorithmType::Custom, which
                // loses the information needed to restore it (and made loading
                // such a grid fail unconditionally). The name round-trips, and is
                // also immune to a renumbering of the enum.
                std::string, // ac_algo name
                std::string, // dc_algo name
                // algo config (scaling/refactor/line-search params, appended;
                // old pickles are version-gated): int_params, real_params
                std::tuple<std::vector<int>, std::vector<double> >,  // ac_algo_config
                std::tuple<std::vector<int>, std::vector<double> >,  // dc_algo_config
                // relevant kwargs the grid was built with (eg by init_from_pypowsybl), as
                // a string->string map flattened into parallel key/value vectors (appended;
                // old pickles are version-gated). See get_init_kwargs()/set_init_kwargs().
                std::vector<std::string>,  // init_kwargs keys
                std::vector<std::string>,  // init_kwargs values
                // fused-bus representative lookup (appended; old pickles/binary
                // formats are version-gated). See get_bus_fusion_rep().
                std::vector<int>  // bus_fusion_rep
                >;

        // named indices into the StateRes tuple above (get_state()/set_state()
        // use these instead of raw std::get<N> literals so the two stay in
        // sync when fields are reordered or appended)
        static const std::size_t VERSION_MAJOR_ID = 0;
        static const std::size_t VERSION_MEDIUM_ID = 1;
        static const std::size_t VERSION_MINOR_ID = 2;
        static const std::size_t LS_TO_ORIG_ID = 3;
        static const std::size_t INIT_VM_PU_ID = 4;
        static const std::size_t SN_MVA_ID = 5;
        static const std::size_t BUS_STATUS_ID = 6;
        static const std::size_t SUBSTATION_ID = 7;
        static const std::size_t LINE_ID = 8;
        static const std::size_t SHUNT_ID = 9;
        static const std::size_t TRAFO_ID = 10;
        static const std::size_t GEN_ID = 11;
        static const std::size_t LOAD_ID = 12;
        static const std::size_t SGEN_ID = 13;
        static const std::size_t STORAGE_ID = 14;
        static const std::size_t HVDC_ID = 15;
        static const std::size_t SVC_ID = 16;
        static const std::size_t AC_ALGO_NAME_ID = 17;
        static const std::size_t DC_ALGO_NAME_ID = 18;
        static const std::size_t AC_ALGO_CONFIG_ID = 19;
        static const std::size_t DC_ALGO_CONFIG_ID = 20;
        static const std::size_t INIT_KWARGS_KEYS_ID = 21;
        static const std::size_t INIT_KWARGS_VALUES_ID = 22;
        static const std::size_t BUS_FUSION_REP_ID = 23;

        LSGrid():
          timer_last_ac_pf_(0.),
          timer_last_dc_pf_(0.),
          algo_controler_(),
          compute_results_(true),
          init_vm_pu_(1.04),
          sn_mva_(1.0),
          max_nb_bus_per_sub_(2){
            _algo.change_algorithm(AlgorithmType::NR_SparseLU);
            _dc_algo.change_algorithm(AlgorithmType::DC_SparseLU);
            _algo.set_lsgrid(this);
            _dc_algo.set_lsgrid(this);
            algo_controler_.ac_algo_controler().tell_all_changed();
            algo_controler_.dc_algo_controler().tell_all_changed();
        }
        LSGrid(const LSGrid & other);
        // LSGrid(LSGrid && other) noexcept = default;  // TODO
        LSGrid copy() const {
            LSGrid res(*this);
            return res;
        }
        ~LSGrid() noexcept = default;

        void set_ls_to_orig(const Eigen::Ref<const IntVect> & ls_to_orig);  // set both _ls_to_orig and _orig_to_ls
        void set_orig_to_ls(const Eigen::Ref<const IntVect> & orig_to_ls);  // set both _orig_to_ls and _ls_to_orig
        [[nodiscard]] const IntVect & get_ls_to_orig(void) const {return _ls_to_orig;}
        [[nodiscard]] const IntVect & get_orig_to_ls(void) const {return _orig_to_ls;}
        double timer_last_ac_pf() const {return timer_last_ac_pf_;}
        double timer_last_dc_pf() const {return timer_last_dc_pf_;}

        /**
         * @brief Return the total number of buses (both connected and disconnected)
         * 
         * For the number of connected buses see `nb_bus()`
         * 
         * @return Eigen::Index 
         */
        size_t total_bus() const {return substations_.nb_bus();}

        const SolverBusIdVect & id_me_to_ac_solver() const {return id_me_to_ac_solver_;}
        const GlobalBusIdVect & id_ac_solver_to_me() const {return id_ac_solver_to_me_;}
        const SolverBusIdVect & id_me_to_dc_solver() const {return id_me_to_dc_solver_;}
        const GlobalBusIdVect & id_dc_solver_to_me() const {return id_dc_solver_to_me_;}

        std::vector<int> id_me_to_ac_solver_numpy() const {
            return id_me_to_ac_solver_.to_int_vector();
        }
        std::vector<int> id_ac_solver_to_me_numpy() const {
            return id_ac_solver_to_me_.to_int_vector();
        }
        std::vector<int> id_me_to_dc_solver_numpy() const {
            return id_me_to_dc_solver_.to_int_vector();
        }
        std::vector<int> id_dc_solver_to_me_numpy() const {
            return id_dc_solver_to_me_.to_int_vector();
        }

        // retrieve the underlying data (raw class)
        [[nodiscard]] const GeneratorContainer & get_generators_as_data() const {return generators_;}
        // turned off generators are not pv
        void turnedoff_no_pv(){
            algo_controler_.ac_algo_controler().has_pv_changed(); 
            algo_controler_.dc_algo_controler().has_pv_changed();
            generators_.turnedoff_no_pv(algo_controler_);
        }
        // turned off generators are pv
        void turnedoff_pv(){
            algo_controler_.ac_algo_controler().has_pv_changed();
            algo_controler_.dc_algo_controler().has_pv_changed();
            generators_.turnedoff_pv(algo_controler_);
        }
        [[nodiscard]] bool get_turnedoff_gen_pv() const {return generators_.get_turnedoff_gen_pv();}
        void update_slack_weights(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & could_be_slack){
            generators_.update_slack_weights(could_be_slack, algo_controler_);
        }
        void update_slack_weights_by_id(const Eigen::Ref<const IntVect> & slack_ids){
            generators_.update_slack_weights_by_id(slack_ids, algo_controler_);
        }

        // Force a given (gridmodel) bus to be the angle reference among the slack
        // buses, WITHOUT changing the slack set or weights: the slack list is
        // reordered at solve time so this bus becomes slack_ids[0] (the NR
        // reference). Pass -1 to clear (default: natural generator order). Used to
        // align the base ac_pf with ContingencyAnalysis::pick_reference_slack() so
        // the GPU companion inherits a reference stranded by the fewest
        // contingencies. Triggers a slack re-evaluation on the next solve.
        void set_reference_slack_bus(int bus_id){
            _forced_ref_slack_bus_id = bus_id;
            algo_controler_.ac_algo_controler().tell_slack_participate_changed();
            algo_controler_.dc_algo_controler().tell_slack_participate_changed();
        }
        [[nodiscard]] int get_reference_slack_bus() const {return _forced_ref_slack_bus_id;}

        [[nodiscard]] const SGenContainer & get_static_generators_as_data() const {return sgens_;}
        [[nodiscard]] const LoadContainer & get_loads_as_data() const {return loads_;}
        [[nodiscard]] const LineContainer & get_powerlines_as_data() const {return powerlines_;}
        [[nodiscard]] const TrafoContainer & get_trafos_as_data() const {return trafos_;}
        [[nodiscard]] const HvdcLineContainer & get_dclines_as_data() const {return hvdc_lines_;}
        [[nodiscard]] Eigen::Ref<const RealVect> get_bus_vn_kv() const {return substations_.get_bus_vn_kv();}

        // per-bus min/max operating voltage (kV), optional: empty if never set
        void set_bus_voltage_limits(const Eigen::Ref<const RealVect> & bus_vmin_kv, const Eigen::Ref<const RealVect> & bus_vmax_kv){
            substations_.init_bus_voltage_limits(bus_vmin_kv, bus_vmax_kv);
        }
        [[nodiscard]] Eigen::Ref<const RealVect> get_bus_vmin_kv() const {return substations_.get_bus_vmin_kv();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_bus_vmax_kv() const {return substations_.get_bus_vmax_kv();}

        std::tuple<int, int> assign_slack_to_most_connected();
        void consider_only_main_component();
        /**
         * Not relevant for dc lines, which always have the default to 
         * synch both sides !
         */
        void set_ignore_status_global(bool ignore_status_global){
            powerlines_.set_ignore_status_global(ignore_status_global);
            trafos_.set_ignore_status_global(ignore_status_global);
        }
        [[nodiscard]] bool get_ignore_status_global() const{
            return powerlines_.get_ignore_status_global();
        }
        /**
         * Not relevant for dc lines, which always have the default to 
         * synch both sides !
         */
        void set_synch_status_both_side(bool synch_status_both_side){
            powerlines_.set_synch_status_both_side(synch_status_both_side);
            trafos_.set_synch_status_both_side(synch_status_both_side);
        }
        [[nodiscard]] bool get_synch_status_both_side() const{
            return powerlines_.get_synch_status_both_side();
        }

        // solver "control"
        void change_algorithm(const AlgorithmType & type){

            if(_algo.is_fdpf(type)) init_fdpf_coeffs();

            if(_algo.is_dc(type)){
                _dc_algo.change_algorithm(type);
                algo_controler_.dc_algo_controler().tell_all_changed();
            } else {
                _algo.change_algorithm(type);
                algo_controler_.ac_algo_controler().tell_all_changed();
            }
        }
        // String-based overload: looks up the solver by registry name.
        // For known built-in names the enum-based routing applies;
        // for plugin names the solver always goes to the AC solver slot.
        void change_algorithm(const std::string& name) {
            // Peek at IS_AC to decide which slot to update.
            std::unique_ptr<BaseAlgo> tmp = AlgorithmRegistry::instance().make(name);
            if (tmp->IS_AC) {
                _algo.change_algorithm(name);
                algo_controler_.ac_algo_controler().tell_all_changed();
            }
            else {
                _dc_algo.change_algorithm(name);
                algo_controler_.dc_algo_controler().tell_all_changed();
            }
        }
        std::vector<AlgorithmType> available_default_algorithms() {return _algo.available_default_algorithms(); }
        // Returns all solver names currently registered (built-in + plugins).
        std::vector<std::string> available_algorithm_names() const {
            return AlgorithmRegistry::instance().available_algorithm_names();
        }
        [[nodiscard]] AlgorithmType get_algo_type() const {return _algo.get_type(); }
        [[nodiscard]] AlgorithmType get_dc_algo_type() const {return _dc_algo.get_type(); }
        [[nodiscard]] const AlgorithmSelector & get_algo() const {return _algo;}
        [[nodiscard]] const AlgorithmSelector & get_dc_algo() const {return _dc_algo;}

        // do i compute the results (in terms of P,Q,V or loads, generators and flows on lines
        void deactivate_result_computation(){compute_results_=false;}
        void reactivate_result_computation(){compute_results_=true;}

        // All methods to init this data model, all need to be pair unit when applicable
        void init_bus(unsigned int n_sub, unsigned int n_busbar_per_sub, const Eigen::Ref<const RealVect> & bus_vn_kv, int nb_line, int nb_trafo);
        // Both scalars below must be finite and strictly positive: they are physical
        // per-unit quantities, and a degenerate value does not fail loudly, it produces
        // a *confidently wrong* answer. `sn_mva_` is the base power of the entire
        // per-unit system (Sbus is divided by it, every MW/MVar result multiplied back
        // by it, and even the solver tolerance is scaled by it in ac_pf), and
        // `init_vm_pu_` is the flat-start voltage magnitude. See check_positive_finite.
        void set_init_vm_pu(real_type init_vm_pu) {
            check_positive_finite(init_vm_pu, "init_vm_pu");
            init_vm_pu_ = init_vm_pu;
        }
        [[nodiscard]] real_type get_init_vm_pu() const {return init_vm_pu_;}
        void set_sn_mva(real_type sn_mva) {
            check_positive_finite(sn_mva, "sn_mva");
            if(sn_mva == sn_mva_) return;
            sn_mva_ = sn_mva;
            // `sn_mva_` is the base power of the whole per-unit system: it is passed
            // to fillYbus / fillBdc AND divides fillSbus_me, so changing it
            // invalidates every cached matrix and injection vector of both families.
            // Grids are normally built with their final sn_mva and never touch it
            // again, so the blunt (and obviously correct) invalidation is free in
            // practice.
            prevent_cache_reuse();
        }
        [[nodiscard]] real_type get_sn_mva() const {return sn_mva_;}

        // Throws std::runtime_error unless `value` is finite and > 0. Written as
        // `!(value > 0.)` on purpose: `value <= 0.` is FALSE for NaN, so the naive
        // form lets a NaN straight through.
        static void check_positive_finite(real_type value, const char * name);

        void init_powerlines(const Eigen::Ref<const RealVect> & branch_r,
                             const Eigen::Ref<const RealVect> & branch_x,
                             const Eigen::Ref<const CplxVect> & branch_h,
                             const Eigen::Ref<const Eigen::VectorXi> & branch_from_id,
                             const Eigen::Ref<const Eigen::VectorXi> & branch_to_id
                             ){
            powerlines_.init(branch_r, branch_x, branch_h, branch_from_id, branch_to_id);
        }
        void init_powerlines_full(const Eigen::Ref<const RealVect> & branch_r,
                                  const Eigen::Ref<const RealVect> & branch_x,
                                  const Eigen::Ref<const CplxVect> & branch_h1,
                                  const Eigen::Ref<const CplxVect> & branch_h2,
                                  const Eigen::Ref<const Eigen::VectorXi> & branch_from_id,
                                  const Eigen::Ref<const Eigen::VectorXi> & branch_to_id
                             ){
            powerlines_.init(branch_r, branch_x, branch_h1,
                             branch_h2, branch_from_id, 
                             branch_to_id);
        }

        void init_shunt(const Eigen::Ref<const RealVect> & shunt_p_mw,
                        const Eigen::Ref<const RealVect> & shunt_q_mvar,
                        const Eigen::Ref<const Eigen::VectorXi> & shunt_bus_id){
            shunts_.init(shunt_p_mw, shunt_q_mvar, shunt_bus_id);
        }
        void init_trafo_pandapower(const Eigen::Ref<const RealVect> & trafo_r,
                                   const Eigen::Ref<const RealVect> & trafo_x,
                                   const Eigen::Ref<const CplxVect> & trafo_b,
                                   const Eigen::Ref<const RealVect> & trafo_tap_step_pct,
                                   const Eigen::Ref<const RealVect> & trafo_tap_pos,
                                   const Eigen::Ref<const RealVect> & trafo_shift_degree,
                                   const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                                   const Eigen::Ref<const Eigen::VectorXi> & bus1_id,
                                   const Eigen::Ref<const Eigen::VectorXi> & bus2_id,
                                   bool ignore_tap_side_for_shift
                                   ){
            trafos_.init(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_shift_degree,
                         trafo_tap_hv, bus1_id, bus2_id, ignore_tap_side_for_shift);
        }
        void init_trafo(const Eigen::Ref<const RealVect> & trafo_r,
                        const Eigen::Ref<const RealVect> & trafo_x,
                        const Eigen::Ref<const CplxVect> & trafo_b,
                        const Eigen::Ref<const RealVect> & trafo_ratio,
                        const Eigen::Ref<const RealVect> & trafo_shift_degree,
                        const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::Ref<const Eigen::VectorXi> & bus1_id,
                        const Eigen::Ref<const Eigen::VectorXi> & bus2_id,
                        bool ignore_tap_side_for_shift
                           ){
            trafos_.init(trafo_r, trafo_x, trafo_b, trafo_ratio, trafo_shift_degree,
                         trafo_tap_hv, bus1_id, bus2_id, ignore_tap_side_for_shift);
        }

        void init_generators(const Eigen::Ref<const RealVect> & generators_p,
                             const Eigen::Ref<const RealVect> & generators_v,
                             const Eigen::Ref<const RealVect> & generators_min_q,
                             const Eigen::Ref<const RealVect> & generators_max_q,
                             const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id){
            generators_.init(generators_p, generators_v, generators_min_q, generators_max_q, generators_bus_id);
        }
        void init_generators_full(const Eigen::Ref<const RealVect> & generators_p,
                                  const Eigen::Ref<const RealVect> & generators_v,
                                  const Eigen::Ref<const RealVect> & generators_q,
                                  const std::vector<bool> & voltage_regulator_on,
                                  const Eigen::Ref<const RealVect> & generators_min_q,
                                  const Eigen::Ref<const RealVect> & generators_max_q,
                                  const Eigen::Ref<const Eigen::VectorXi> & generators_bus_id){
            generators_.init_full(generators_p, generators_v, generators_q, voltage_regulator_on,
                                  generators_min_q, generators_max_q, generators_bus_id);
        }
        void init_loads(const Eigen::Ref<const RealVect> & loads_p,
                        const Eigen::Ref<const RealVect> & loads_q,
                        const Eigen::Ref<const Eigen::VectorXi> & loads_bus_id){
            loads_.init(loads_p, loads_q, loads_bus_id);
        }
        void init_sgens(const Eigen::Ref<const RealVect> & sgen_p,
                        const Eigen::Ref<const RealVect> & sgen_q,
                        const Eigen::Ref<const RealVect> & sgen_pmin,
                        const Eigen::Ref<const RealVect> & sgen_pmax,
                        const Eigen::Ref<const RealVect> & sgen_qmin,
                        const Eigen::Ref<const RealVect> & sgen_qmax,
                        const Eigen::Ref<const Eigen::VectorXi> & sgen_bus_id){
            sgens_.init(sgen_p, sgen_q, sgen_pmin, sgen_pmax, sgen_qmin, sgen_qmax, sgen_bus_id);
        }
        void init_svcs(const std::vector<int> & regulation_mode,
                       const Eigen::Ref<const RealVect> & target_vm_pu,
                       const Eigen::Ref<const RealVect> & q_setpoint_mvar,
                       const Eigen::Ref<const RealVect> & slope_pu,
                       const Eigen::Ref<const RealVect> & b_min,
                       const Eigen::Ref<const RealVect> & b_max,
                       const Eigen::Ref<const Eigen::VectorXi> & regulated_bus_id,
                       const Eigen::Ref<const Eigen::VectorXi> & svc_bus_id){
            svcs_.init(regulation_mode, target_vm_pu, q_setpoint_mvar, slope_pu,
                       b_min, b_max, regulated_bus_id, svc_bus_id);
        }
        void init_storages(const Eigen::Ref<const RealVect> & storages_p,
                           const Eigen::Ref<const RealVect> & storages_q,
                           const Eigen::Ref<const Eigen::VectorXi> & storages_bus_id){
            storages_.init(storages_p, storages_q, storages_bus_id);
        }
        void init_dclines(const Eigen::Ref<const Eigen::VectorXi> & branch_from_id,
                          const Eigen::Ref<const Eigen::VectorXi> & branch_to_id,
                          const Eigen::Ref<const RealVect> & p_mw,
                          const Eigen::Ref<const RealVect> & loss_percent,
                          const Eigen::Ref<const RealVect> & loss_mw,
                          const Eigen::Ref<const RealVect> & vm1_pu,
                          const Eigen::Ref<const RealVect> & vm2_pu,
                          const Eigen::Ref<const RealVect> & min_q1,
                          const Eigen::Ref<const RealVect> & max_q1,
                          const Eigen::Ref<const RealVect> & min_q2,
                          const Eigen::Ref<const RealVect> & max_q2){
            hvdc_lines_.init_legacy(branch_from_id, branch_to_id, p_mw,
                                    loss_percent, loss_mw, vm1_pu, vm2_pu,
                                    min_q1, max_q1, min_q2, max_q2);
        }
        /**
         * Full IIDM-style hvdc initialization (two converter stations - VSC or
         * LCC - per line, optional angle-droop). See HvdcLineContainer::init.
         */
        void init_hvdc_lines(const Eigen::Ref<const Eigen::VectorXi> & bus1_id,
                             const Eigen::Ref<const Eigen::VectorXi> & bus2_id,
                             const std::vector<int> & type1,
                             const std::vector<int> & type2,
                             const Eigen::Ref<const RealVect> & loss_factor1,
                             const Eigen::Ref<const RealVect> & loss_factor2,
                             const std::vector<bool> & vreg1_on,
                             const std::vector<bool> & vreg2_on,
                             const Eigen::Ref<const RealVect> & vm1_pu,
                             const Eigen::Ref<const RealVect> & vm2_pu,
                             const Eigen::Ref<const RealVect> & q1_setpoint_mvar,
                             const Eigen::Ref<const RealVect> & q2_setpoint_mvar,
                             const Eigen::Ref<const RealVect> & min_q1,
                             const Eigen::Ref<const RealVect> & max_q1,
                             const Eigen::Ref<const RealVect> & min_q2,
                             const Eigen::Ref<const RealVect> & max_q2,
                             const Eigen::Ref<const RealVect> & power_factor1,
                             const Eigen::Ref<const RealVect> & power_factor2,
                             const std::vector<int> & converters_mode,
                             const Eigen::Ref<const RealVect> & p_setpoint_mw,
                             const Eigen::Ref<const RealVect> & r_ohm,
                             const Eigen::Ref<const RealVect> & nominal_v_kv,
                             const std::vector<bool> & droop_enabled,
                             const Eigen::Ref<const RealVect> & droop_p0_mw,
                             const Eigen::Ref<const RealVect> & droop_mw_per_deg,
                             const Eigen::Ref<const RealVect> & pmax_1to2_mw,
                             const Eigen::Ref<const RealVect> & pmax_2to1_mw){
            const RealVect no_legacy_loss = RealVect::Zero(bus1_id.size());
            hvdc_lines_.init(bus1_id, bus2_id, type1, type2,
                             loss_factor1, loss_factor2, vreg1_on, vreg2_on,
                             vm1_pu, vm2_pu, q1_setpoint_mvar, q2_setpoint_mvar,
                             min_q1, max_q1, min_q2, max_q2,
                             power_factor1, power_factor2,
                             converters_mode, p_setpoint_mw, r_ohm, nominal_v_kv,
                             no_legacy_loss, no_legacy_loss,
                             droop_enabled, droop_p0_mw, droop_mw_per_deg,
                             pmax_1to2_mw, pmax_2to1_mw);
        }

        void init_bus_status(){
            substations_.disconnect_all_buses();

            powerlines_.reconnect_connected_buses(substations_);
            shunts_.reconnect_connected_buses(substations_);
            trafos_.reconnect_connected_buses(substations_);
            generators_.reconnect_connected_buses(substations_);
            loads_.reconnect_connected_buses(substations_);
            sgens_.reconnect_connected_buses(substations_);
            storages_.reconnect_connected_buses(substations_);
            hvdc_lines_.reconnect_connected_buses(substations_);
            const std::vector<bool> new_status = substations_.get_bus_status();

            // Each family is compared against the connectivity ITS OWN cache was
            // built with: an AC powerflow that finds the AC snapshot up to date
            // must not conclude anything about DC, which may be several grid
            // modifications behind. A snapshot is refreshed only by
            // _mark_cache_valid(), ie when that family has actually rebuilt.
            _flag_dimension_change(new_status, last_bus_status_saved_, algo_controler_.ac_algo_controler());
            _flag_dimension_change(new_status, last_bus_status_dc_, algo_controler_.dc_algo_controler());
        }
        void set_substation_names(const std::vector<std::string> & sub_names){
            substations_.init_sub_names(sub_names);
        }
        [[nodiscard]] const std::vector<std::string> & get_substation_names()const {
            return substations_.get_sub_names();
        }

        void add_gen_slackbus(int gen_id, real_type weight);
        void remove_gen_slackbus(int gen_id);

        //pickle
        LSGrid::StateRes get_state() const ;
        // `restore_algorithm == true` (the default) also re-selects the AC / DC
        // solver the grid was saved with, and re-applies its configuration. It
        // throws if that solver is not registered here (a plugin that has not
        // been loaded, or a built-in needing an optional backend this build
        // lacks). Pass false to load the grid *data* only and keep this object's
        // current (default) solvers -- see load_binary_without_algorithm().
        void set_state(LSGrid::StateRes & my_state, bool restore_algorithm = true) ;

        // Whole-grid consistency check. Verifies that every index the grid carries
        // (element bus ids, substation ids, position in the topology vector,
        // generator slack and remote-regulated bus references) is in range for this
        // grid. Throws std::out_of_range on an out-of-range index and
        // std::runtime_error on a structural error.
        //
        // Called automatically at the end of set_state() (so loading a pickle or a
        // binary file cannot leave an out-of-range index that would later cause an
        // out-of-bounds read/write during a powerflow), and exposed to Python so the
        // grid loaders (pandapower / pypowsybl / matpower / powermodels) can call it
        // right after building a grid. Runs in O(number of elements), off the solver
        // hot path.
        void check_grid() const;

    private:
        // Re-select, by registry name, the algorithm a grid was saved with.
        // Throws a message naming the missing solver (and how to get it) when it
        // is not registered here -- typically a plugin that has not been loaded,
        // or a built-in needing an optional backend this build lacks.
        static void _restore_algorithm(AlgorithmSelector & algo_selector,
                                       const std::string & name,
                                       const char * ac_or_dc);
    public:

        // fast binary serialization (additive alternative to pickle, see
        // BinaryArchive.hpp -- readable by any lightsim2grid version sharing
        // the same BINARY_FORMAT_VERSION)
        void save_binary(const std::string & path, bool atomic = true) const;
        static LSGrid load_binary(const std::string & path);
        static const char * binary_type_tag() { return "LSGrid"; }  // written into / checked against the binary file header
        // binary-load hook (detected by load_binary_generic): rewrites the
        // version strings embedded in a StateRes read from a binary file to
        // this build's values, so the pickle-shared exact-version check in
        // set_state() does not reject a cross-version load that the binary
        // format number allows
        static void fixup_binary_state(LSGrid::StateRes & state);

        // algo config (scaling/refactor policy params) — also part of StateRes,
        // see LSGrid::get_state()/set_state()
        [[nodiscard]] AlgoConfig get_ac_algo_config() const { return _algo.get_config(); }
        void set_ac_algo_config(const AlgoConfig& cfg) { _algo.set_config(cfg); }
        [[nodiscard]] AlgoConfig get_dc_algo_config() const { return _dc_algo.get_config(); }
        void set_dc_algo_config(const AlgoConfig& cfg) { _dc_algo.set_config(cfg); }
        template<class T>
        void check_size(const T& my_state)
        {
            // currently un used
            unsigned int size_th = 6;
            if (my_state.size() != size_th)
            {
                // TODO more explicit error message
                throw std::runtime_error("Invalid state when loading LightSim::LSGrid");
            }
        }

        //powerflows
        // control the need to refactorize the topology
        // ---- cache reuse -------------------------------------------------------
        // A powerflow reuses what the previous one of the SAME family built --
        // the bus labelling, Ybus / Sbus, the pv-pq split, the slack weights --
        // and only re-stamps what the grid reports as modified since. This is on
        // by default and needs nothing from the caller: every powerflow that
        // converges marks its own family "in sync" on the way out (see
        // process_results / _mark_cache_valid). Every grid modification made
        // through LSGrid's own API sets the flags that undo that, per family.
        //
        // Turning a family off makes it rebuild everything, every solve. The
        // answer is identical either way, so this is a debugging switch (is a
        // wrong result a caching bug?) or a safety net for code that mutates the
        // containers behind LSGrid's back -- not something a normal user needs.
        void allow_ac_cache_reuse(bool allowed){allow_ac_cache_reuse_ = allowed;}
        void allow_dc_cache_reuse(bool allowed){allow_dc_cache_reuse_ = allowed;}
        void allow_cache_reuse(bool allowed){allow_ac_cache_reuse_ = allowed; allow_dc_cache_reuse_ = allowed;}
        [[nodiscard]] bool get_allow_ac_cache_reuse() const {return allow_ac_cache_reuse_;}
        [[nodiscard]] bool get_allow_dc_cache_reuse() const {return allow_dc_cache_reuse_;}
        // true only when BOTH families are allowed to reuse their cache
        [[nodiscard]] bool get_allow_cache_reuse() const {return allow_ac_cache_reuse_ && allow_dc_cache_reuse_;}

        /**
         * Throw away everything a family cached and start its next powerflow from
         * scratch (matrices, bus labelling, factorization, and the connectivity
         * snapshot used to detect topology changes).
         *
         * Needed only when the grid was modified through a path that bypasses
         * LSGrid's own `change_*` / `deactivate_*` / `reactivate_*` methods --
         * those already set the narrower, per-family flags themselves -- or when
         * in doubt after a change of unclear scope. Always correct, just more
         * expensive than letting the narrower flags do their job.
         *
         * Formerly `tell_solver_need_reset()`, still available under that name.
         */
        void prevent_ac_cache_reuse(){
            last_bus_status_saved_ = std::vector<bool>(substations_.nb_bus(), false);
            algo_controler_.ac_algo_controler().tell_solver_need_reset();
        }
        void prevent_dc_cache_reuse(){
            last_bus_status_dc_ = std::vector<bool>(substations_.nb_bus(), false);
            algo_controler_.dc_algo_controler().tell_solver_need_reset();
        }
        void prevent_cache_reuse(){prevent_ac_cache_reuse(); prevent_dc_cache_reuse();}
        // backward-compatible name of `prevent_cache_reuse()`
        void tell_solver_need_reset(){prevent_cache_reuse();}

        /**
         * Inform the grid that a powerflow has been run and that every "past changes"
         * can be "forgotten" : if a solver is re run, then some things are cached to
         * avoid un necessary computation.
         *
         * Historical, manual counterpart of the automatic marking described above.
         * When both families cache automatically (the default) there is nothing
         * left for it to do and it returns immediately; it only has an effect on a
         * family whose automatic reuse was turned off -- and even then that
         * family's `allow_*_cache_reuse(false)` still wins, so the claim is
         * recorded but the next powerflow of that family rebuilds anyway.
         *
         * Kept for backward compatibility: new code should simply not call it.
         */
        void unset_changes(){
            if(allow_ac_cache_reuse_ && allow_dc_cache_reuse_) return;  // already done, after every powerflow
            _mark_cache_valid(true);
            _mark_cache_valid(false);
        }  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.

        void tell_recompute_ybus(){algo_controler_.ac_algo_controler().tell_recompute_ybus(); algo_controler_.dc_algo_controler().tell_recompute_ybus();}
        void tell_recompute_sbus(){algo_controler_.ac_algo_controler().tell_recompute_sbus(); algo_controler_.dc_algo_controler().tell_recompute_sbus();}
        void tell_ybus_change_sparsity_pattern(){algo_controler_.ac_algo_controler().tell_ybus_change_sparsity_pattern(); algo_controler_.dc_algo_controler().tell_ybus_change_sparsity_pattern();}
        [[nodiscard]] const AlgoControl & get_ac_algo_controler() const {return algo_controler_.ac_algo_controler();}
        [[nodiscard]] const AlgoControl & get_dc_algo_controler() const {return algo_controler_.dc_algo_controler();}

        // dc powerflow
        CplxVect dc_pf(const Eigen::Ref<const CplxVect> & Vinit,
                       int max_iter,  // not used for DC
                       real_type tol  // not used for DC
                       );

        /**
         * @brief Retrieve the PTDF matrice, with bus labeled in the "gridmodel" format.
         * Deactivated buses are represented by a column of 0.
         * 
         * It has the size (nb_line + nb_trafo, total_bus)
         * 
         * @return RealMat 
        */
        RealMat get_ptdf();

        /**
         * @brief Retrieve the LODF matrice.
         * 
         * 
         * It has the size (nb_line + nb_trafo, nb_line + nb_trafo)
         * 
         * @return RealMat 
        */
        RealMat get_lodf();

        /**
         * @brief Retrieve the PTDF matrice, with bus labeled in the "solver" format.
         * 
         * Deactivated buses are represented by a column of 0.
         * 
         * It has the size (nb_line + nb_trafo, nb_bus)
         * NB nb_bus is the number of activated buses !
         * 
         * @return RealMat 
        */
        RealMat get_ptdf_solver();
        
        Eigen::SparseMatrix<real_type> get_Bf_solver();
        Eigen::SparseMatrix<real_type> get_Bf();

        // ac powerflow
        CplxVect ac_pf(const Eigen::Ref<const CplxVect> & Vinit,
                       int max_iter,
                       real_type tol);

        // check the kirchoff law
        CplxVect check_solution(const Eigen::Ref<const CplxVect> & V, bool check_q_limits);

        // deactivate a bus. Be careful, if a bus is deactivated, but an element is
        // still connected to it, it will throw an exception
        void deactivate_bus(GlobalBusId global_bus_id) {
            if(substations_.is_bus_connected(global_bus_id)){
                // bus was connected, dim of matrix change
                algo_controler_.ac_algo_controler().need_reset_solver();
                algo_controler_.ac_algo_controler().need_recompute_sbus();
                algo_controler_.ac_algo_controler().need_recompute_ybus();
                algo_controler_.ac_algo_controler().ybus_change_sparsity_pattern();
                algo_controler_.dc_algo_controler().need_reset_solver();
                algo_controler_.dc_algo_controler().need_recompute_sbus();
                algo_controler_.dc_algo_controler().need_recompute_ybus();
                algo_controler_.dc_algo_controler().ybus_change_sparsity_pattern();
                GenericContainer::_generic_deactivate(global_bus_id, substations_);
            }
        }
        void deactivate_bus_python(int global_bus_id) {
            deactivate_bus(GlobalBusId(global_bus_id));
        }

        // if a bus is connected, but isolated, it will make the powerflow diverge
        void reactivate_bus(GlobalBusId global_bus_id) {
            if(!substations_.is_bus_connected(global_bus_id)){
                // bus was not connected, dim of matrix change
                algo_controler_.ac_algo_controler().need_reset_solver();
                algo_controler_.ac_algo_controler().need_recompute_sbus();
                algo_controler_.ac_algo_controler().need_recompute_ybus();
                algo_controler_.ac_algo_controler().ybus_change_sparsity_pattern();
                algo_controler_.dc_algo_controler().need_reset_solver();
                algo_controler_.dc_algo_controler().need_recompute_sbus();
                algo_controler_.dc_algo_controler().need_recompute_ybus();
                algo_controler_.dc_algo_controler().ybus_change_sparsity_pattern();
                GenericContainer::_generic_reactivate(global_bus_id, substations_); 
            }
        }
        void reactivate_bus_python(int global_bus_id) {
            reactivate_bus(GlobalBusId(global_bus_id));
        }

        /**
         * @brief Return the total number of connected buses !
         * 
         * For the total number of buses, see `total_bus()`
         * 
         * @return int 
         */
        int nb_connected_bus() const {return substations_.nb_connected_bus();}
        size_t nb_powerline() const {return powerlines_.nb();}
        size_t nb_trafo() const {return trafos_.nb();}

        // read only data accessor
        [[nodiscard]] const SubstationContainer & get_substations() const {return substations_;}
        [[nodiscard]] const LineContainer & get_lines() const {return powerlines_;}
        [[nodiscard]] const HvdcLineContainer & get_dclines() const {return hvdc_lines_;}
        [[nodiscard]] const TrafoContainer & get_trafos() const {return trafos_;}
        [[nodiscard]] const GeneratorContainer & get_generators() const {return generators_;}
        [[nodiscard]] const LoadContainer & get_loads() const {return loads_;}
        [[nodiscard]] const StorageContainer & get_storages() const {return storages_;}
        [[nodiscard]] const SGenContainer & get_static_generators() const {return sgens_;}
        [[nodiscard]] const SvcContainer & get_svcs() const {return svcs_;}
        [[nodiscard]] const ShuntContainer & get_shunts() const {return shunts_;}
        [[nodiscard]] const std::vector<bool> & get_bus_status() const {return substations_.get_bus_status();}
        
        void set_line_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, powerlines_.nb(), "set_line_names");
            powerlines_.set_names(names);
        }
        [[nodiscard]] const std::vector<std::string> & get_line_names() const {return powerlines_.get_names();}
        void set_dcline_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, hvdc_lines_.nb(), "set_dcline_names");
            hvdc_lines_.set_names(names);
        }
        void set_trafo_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, trafos_.nb(), "set_trafo_names");
            trafos_.set_names(names);
        }
        [[nodiscard]] const std::vector<std::string> & get_trafo_names() const {return trafos_.get_names();}
        // per-side current limit, in kA, optional: empty if never set
        void set_line_current_limit_side1(const Eigen::Ref<const RealVect> & limit_a1_ka){
            GenericContainer::check_size(limit_a1_ka, powerlines_.nb(), "set_line_current_limit_side1");
            powerlines_.set_limit_a1_ka(limit_a1_ka);
        }
        void set_line_current_limit_side2(const Eigen::Ref<const RealVect> & limit_a2_ka){
            GenericContainer::check_size(limit_a2_ka, powerlines_.nb(), "set_line_current_limit_side2");
            powerlines_.set_limit_a2_ka(limit_a2_ka);
        }
        void set_trafo_current_limit_side1(const Eigen::Ref<const RealVect> & limit_a1_ka){
            GenericContainer::check_size(limit_a1_ka, trafos_.nb(), "set_trafo_current_limit_side1");
            trafos_.set_limit_a1_ka(limit_a1_ka);
        }
        void set_trafo_current_limit_side2(const Eigen::Ref<const RealVect> & limit_a2_ka){
            GenericContainer::check_size(limit_a2_ka, trafos_.nb(), "set_trafo_current_limit_side2");
            trafos_.set_limit_a2_ka(limit_a2_ka);
        }
        void set_gen_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, generators_.nb(), "set_gen_names");
            generators_.set_names(names);
        }
        void set_load_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, loads_.nb(), "set_load_names");
            loads_.set_names(names);
        }
        void set_storage_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, storages_.nb(), "set_storage_names");
            storages_.set_names(names);
        }
        void set_sgen_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, sgens_.nb(), "set_sgen_names");
            sgens_.set_names(names);
        }
        void set_shunt_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, shunts_.nb(), "set_shunt_names");
            shunts_.set_names(names);
        }

        //deactivate a powerline (disconnect it)
        void deactivate_powerline(int powerline_id) {
            powerlines_.deactivate(powerline_id, algo_controler_);
        }
        void reactivate_powerline(int powerline_id) {
            powerlines_.reactivate(powerline_id, algo_controler_);
        }
        // per-side (de)activation: model a powerline connected on one terminal only
        // ("half-open"). With set_synch_status_both_side(false), the other side stays
        // as is; with the default (true) the other side follows (whole-line behaviour).
        void deactivate_powerline_side1(int powerline_id) { powerlines_.deactivate_side_1(powerline_id, algo_controler_); }
        void deactivate_powerline_side2(int powerline_id) { powerlines_.deactivate_side_2(powerline_id, algo_controler_); }
        void reactivate_powerline_side1(int powerline_id) { powerlines_.reactivate_side_1(powerline_id, algo_controler_); }
        void reactivate_powerline_side2(int powerline_id) { powerlines_.reactivate_side_2(powerline_id, algo_controler_); }

        /**
         * Change the bus on the "side 1" of the powerline powerline_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id". **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus1_powerline(int powerline_id, GridModelBusId new_gridmodel_bus_id) {
            powerlines_.change_bus_side_1(
                powerline_id,
                new_gridmodel_bus_id,
                algo_controler_,    
                substations_);
        }
        void change_bus1_powerline_python(int powerline_id, int new_gridmodel_bus_id) {
            change_bus1_powerline(powerline_id, GridModelBusId(new_gridmodel_bus_id));
        }

        /**
         * Change the bus on the "side 2" of the powerline powerline_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus2_powerline(int powerline_id, GridModelBusId new_gridmodel_bus_id) {
            powerlines_.change_bus_side_2(
                powerline_id,
                new_gridmodel_bus_id,
                algo_controler_,
                substations_);
        }
        void change_bus2_powerline_python(int powerline_id, int new_gridmodel_bus_id) {
            change_bus2_powerline(powerline_id, GridModelBusId(new_gridmodel_bus_id));
        }
        int get_bus1_powerline(int powerline_id) const {return powerlines_.get_bus_side_1(powerline_id).cast_int();}
        int get_bus2_powerline(int powerline_id) const {return powerlines_.get_bus_side_2(powerline_id).cast_int();}

        //deactivate trafo
        void deactivate_trafo(int trafo_id) {trafos_.deactivate(trafo_id, algo_controler_); }
        void reactivate_trafo(int trafo_id) {trafos_.reactivate(trafo_id, algo_controler_); }
        // per-side (de)activation of a transformer terminal ("half-open"), see the
        // powerline equivalents above.
        void deactivate_trafo_side1(int trafo_id) { trafos_.deactivate_side_1(trafo_id, algo_controler_); }
        void deactivate_trafo_side2(int trafo_id) { trafos_.deactivate_side_2(trafo_id, algo_controler_); }
        void reactivate_trafo_side1(int trafo_id) { trafos_.reactivate_side_1(trafo_id, algo_controler_); }
        void reactivate_trafo_side2(int trafo_id) { trafos_.reactivate_side_2(trafo_id, algo_controler_); }

        /**
         * Change the bus on the "side 1" of the trafo trafo_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus1_trafo(int trafo_id, GridModelBusId new_gridmodel_bus_id) {
            trafos_.change_bus_side_1(
                trafo_id,
                new_gridmodel_bus_id,
                algo_controler_,
                substations_); 
        }
        void change_bus1_trafo_python(int trafo_id, int new_gridmodel_bus_id) {
            change_bus1_trafo(trafo_id, GridModelBusId(new_gridmodel_bus_id));
        }

        /**
         * Change the bus on the "side 2" of the trafo trafo_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus2_trafo(int trafo_id, GridModelBusId new_gridmodel_bus_id) {
            trafos_.change_bus_side_2(trafo_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus2_trafo_python(int trafo_id, int new_gridmodel_bus_id) {
            change_bus2_trafo(trafo_id, GridModelBusId(new_gridmodel_bus_id));
        }
        int get_bus1_trafo(int trafo_id) const {
            return trafos_.get_bus_side_1(trafo_id).cast_int();
        }
        int get_bus2_trafo(int trafo_id) const {
            return trafos_.get_bus_side_2(trafo_id).cast_int();
        }
        void change_ratio_trafo(int trafo_id, real_type new_ratio){
            trafos_.change_ratio(trafo_id, new_ratio, algo_controler_);
        }

        /**
         * The shift is in radian (not degree !)
         */
        void change_shift_trafo(int trafo_id, real_type new_shift_rad){
            trafos_.change_shift(trafo_id, new_shift_rad, algo_controler_);
        }
        void change_shift_trafo_deg(int trafo_id, real_type new_shift_deg){
            real_type new_shift_rad = new_shift_deg / BaseConstants::my_180_pi_;
            change_shift_trafo(trafo_id, new_shift_rad);
        }
        /**
         * Declare the alpha (phase-shift) dependence of the transformers' series
         * impedance: per transformer, sample points `alpha (rad) -> r/x correction (%)`.
         * The effective r / x is then refreshed from these samples (interpolated on the
         * current shift) whenever the coefficients are rebuilt, so changing the shift /
         * ratio of a phase-shifting transformer keeps its impedance correct -- without
         * any "tap" concept. `enable` is kept false for pandapower (no such data).
         */
        void set_trafo_shift_dependent_rx(bool enable,
                                          const std::vector<std::vector<real_type> > & alpha_rad,
                                          const std::vector<std::vector<real_type> > & rx_corr_pct){
            trafos_.set_shift_dependent_rx(enable, alpha_rad, rx_corr_pct, algo_controler_);
        }

        //load
        void deactivate_load(int load_id) {loads_.deactivate(load_id, algo_controler_); }
        void reactivate_load(int load_id) {loads_.reactivate(load_id, algo_controler_); }

        /**
         * Change the bus on the load load_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus_load(int load_id, GridModelBusId new_gridmodel_bus_id) {
            loads_.change_bus(load_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus_load_python(int load_id, int new_gridmodel_bus_id) {
            change_bus_load(load_id, GridModelBusId(new_gridmodel_bus_id));
        }
        void change_p_load(int load_id, real_type new_p) {loads_.change_p_nothrow(load_id, new_p, algo_controler_); }
        void change_q_load(int load_id, real_type new_q) {loads_.change_q_nothrow(load_id, new_q, algo_controler_); }
        [[nodiscard]] int get_bus_load(int load_id) const {return loads_.get_bus(load_id).cast_int();}

        //generator
        void deactivate_gen(int gen_id) {generators_.deactivate(gen_id, algo_controler_); }
        void reactivate_gen(int gen_id) {generators_.reactivate(gen_id, algo_controler_); }

        /**
         * Change the bus on the generator generator_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus_gen(int gen_id, GridModelBusId new_gridmodel_bus_id) {
            generators_.change_bus(gen_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus_gen_python(int gen_id, int new_gridmodel_bus_id) {
            change_bus_gen(gen_id, GridModelBusId(new_gridmodel_bus_id));
        }
        void change_p_gen(int gen_id, real_type new_p) {generators_.change_p_nothrow(gen_id, new_p, algo_controler_); }
        void change_q_gen(int gen_id, real_type new_q) {generators_.change_q_nothrow(gen_id, new_q, algo_controler_); }
        void change_v_gen(int gen_id, real_type new_v_pu) {generators_.change_v_nothrow(gen_id, new_v_pu, algo_controler_); }
        [[nodiscard]] int get_bus_gen(int gen_id) const {return generators_.get_bus(gen_id).cast_int();}

        //static var compensator (SVC)
        void deactivate_svc(int svc_id) {svcs_.deactivate(svc_id, algo_controler_); }
        void reactivate_svc(int svc_id) {svcs_.reactivate(svc_id, algo_controler_); }
        void change_bus_svc(int svc_id, GridModelBusId new_gridmodel_bus_id) {
            svcs_.change_bus(svc_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus_svc_python(int svc_id, int new_gridmodel_bus_id) {
            change_bus_svc(svc_id, GridModelBusId(new_gridmodel_bus_id));
        }
        [[nodiscard]] int get_bus_svc(int svc_id) const {return svcs_.get_bus(svc_id).cast_int();}
        void set_svc_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, svcs_.nb(), "set_svc_names");
            svcs_.set_names(names);
        }

        //shunt
        void deactivate_shunt(int shunt_id) {shunts_.deactivate(shunt_id, algo_controler_); }
        void reactivate_shunt(int shunt_id) {shunts_.reactivate(shunt_id, algo_controler_); }
        /**
         * Change the bus on the shunt shunt_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus_shunt(int shunt_id, GridModelBusId new_gridmodel_bus_id) {
            shunts_.change_bus(shunt_id, new_gridmodel_bus_id, algo_controler_, substations_);  
        }
        void change_bus_shunt_python(int shunt_id, int new_gridmodel_bus_id) {
            change_bus_shunt(shunt_id, GridModelBusId(new_gridmodel_bus_id));  
        }
        void change_p_shunt(int shunt_id, real_type new_p) {shunts_.change_p_nothrow(shunt_id, new_p, algo_controler_); }
        void change_q_shunt(int shunt_id, real_type new_q) {shunts_.change_q_nothrow(shunt_id, new_q, algo_controler_); }
        [[nodiscard]] int get_bus_shunt(int shunt_id) const {return shunts_.get_bus(shunt_id).cast_int();}

        //static gen
        void deactivate_sgen(int sgen_id) {sgens_.deactivate(sgen_id, algo_controler_); }
        void reactivate_sgen(int sgen_id) {sgens_.reactivate(sgen_id, algo_controler_); }
        /**
         * Change the bus on the static generator sgen_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus_sgen(int sgen_id, GridModelBusId new_gridmodel_bus_id) {
            sgens_.change_bus(sgen_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus_sgen_python(int sgen_id, int new_gridmodel_bus_id) {
            change_bus_sgen(sgen_id, GridModelBusId(new_gridmodel_bus_id));
        }
        void change_p_sgen(int sgen_id, real_type new_p) {sgens_.change_p_nothrow(sgen_id, new_p, algo_controler_); }
        void change_q_sgen(int sgen_id, real_type new_q) {sgens_.change_q_nothrow(sgen_id, new_q, algo_controler_); }
        [[nodiscard]] int get_bus_sgen(int sgen_id) const {return sgens_.get_bus(sgen_id).cast_int();}

        //storage units
        void deactivate_storage(int storage_id) {storages_.deactivate(storage_id, algo_controler_); }
        void reactivate_storage(int storage_id) {storages_.reactivate(storage_id, algo_controler_); }
        /**
         * Change the bus on the storage storage_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus_storage(int storage_id, GridModelBusId new_gridmodel_bus_id) {
            storages_.change_bus(storage_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus_storage_python(int sgen_id, int new_gridmodel_bus_id) {
            change_bus_storage(sgen_id, GridModelBusId(new_gridmodel_bus_id));
        }
        void change_p_storage(int storage_id, real_type new_p) {
               storages_.change_p_nothrow(storage_id, new_p, algo_controler_);
            }
        void change_q_storage(int storage_id, real_type new_q) {storages_.change_q_nothrow(storage_id, new_q, algo_controler_); }
        [[nodiscard]] int get_bus_storage(int storage_id) const {return storages_.get_bus(storage_id).cast_int();}

        //deactivate a powerline (disconnect it)
        void deactivate_dcline(int dcline_id) {hvdc_lines_.deactivate(dcline_id, algo_controler_); }
        void reactivate_dcline(int dcline_id) {hvdc_lines_.reactivate(dcline_id, algo_controler_); }
        // Disconnect only one converter station of an HVDC line ("half-open"): the
        // other station keeps injecting its scheduled P / regulating Q-V, matching
        // OpenLoadFlow which treats a VSC station with its DC partner switched off as
        // a still-active local reactive/voltage-support device (not a dead branch).
        // Unlike powerlines/trafos this is unconditional (HvdcLineContainer's
        // synch_status_both_side_ defaults to false), no `keep_half_open_lines` needed.
        void deactivate_dcline_side1(int dcline_id) {
            hvdc_lines_.deactivate_side_1(dcline_id, algo_controler_);
            hvdc_lines_.disable_droop(dcline_id);  // remote angle is gone, see disable_droop's doc
        }
        void deactivate_dcline_side2(int dcline_id) {
            hvdc_lines_.deactivate_side_2(dcline_id, algo_controler_);
            hvdc_lines_.disable_droop(dcline_id);
        }
        void change_p_dcline(int dcline_id, real_type new_p) {hvdc_lines_.change_p(dcline_id, new_p, algo_controler_); }
        void change_v1_dcline(int dcline_id, real_type new_v_pu) {hvdc_lines_.change_v_side_1(dcline_id, new_v_pu, algo_controler_); }
        void change_v2_dcline(int dcline_id, real_type new_v_pu) {hvdc_lines_.change_v_side_2(dcline_id, new_v_pu, algo_controler_); }
        /**
         * Set the droop regime of an angle-droop (AC emulation) hvdc line.
         * This is an INPUT of the solver: 0 = linear, +1 = saturated 1->2,
         * -1 = saturated 2->1. The saturation logic (against pmax_1to2 /
         * pmax_2to1) is meant to be run between two solves, in Python.
         */
        void set_status_droop_hvdc(int dcline_id, int status) {hvdc_lines_.set_status_droop(dcline_id, status, algo_controler_); }
        [[nodiscard]] int get_status_droop_hvdc(int dcline_id) const {return hvdc_lines_.get_status_droop(dcline_id);}
        [[nodiscard]] Eigen::Ref<const IntVect> get_status_droop_hvdc_vect() const {return hvdc_lines_.get_status_droop_vect();}
        /**
         * Per-solve data of the connected angle-droop (AC emulation) hvdc
         * lines, in solver bus labelling and per-unit. Consumed by the Hvdc
         * extension of the Newton-Raphson system (ac = true) and by the DC
         * algorithm (ac = false). Only valid once `pre_process_solver` ran
         * (ie from within a solver's compute_pf).
         */
        void fill_hvdc_droop_solver_data(HvdcDroopSolverData & data, bool ac) const;

        // Python-facing wrapper around fill_hvdc_droop_solver_data(ac=true):
        // (bus1, bus2, status, p0, k, lf1, lf2, r, pmax12, pmax21), one entry
        // per CONNECTED droop-enabled hvdc line, solver bus numbering, pu.
        // Ground truth for external solvers re-deriving the theta-dependent
        // droop flow contribution to F independently.
        std::tuple<Eigen::VectorXi, Eigen::VectorXi, Eigen::VectorXi,
                   RealVect, RealVect, RealVect, RealVect, RealVect, RealVect, RealVect>
        get_hvdc_droop_data_solver() const {
            HvdcDroopSolverData d;
            fill_hvdc_droop_solver_data(d, true);
            return {d.bus1, d.bus2, d.status, d.p0, d.k, d.lf1, d.lf2, d.r, d.pmax12, d.pmax21};
        }
        /**
         * Per-solve data of the ACTIVE voltage-mode controllers (remote-regulating
         * generators and, later, voltage-mode SVCs), grouped by regulated solver
         * bus, in solver bus labelling and per-unit. Consumed by the VoltageControl
         * extension of the Newton-Raphson system (ac = true only; empty in DC).
         * Throws a clear error on the singular-for-us configurations (see Phase 0
         * probe #3). Only valid once `pre_process_solver` ran.
         */
        void fill_voltage_control_solver_data(VoltageControlSolverData & data, bool ac) const;
        /**
         * Solver-bus ids of the slack buses that need a free Vm unknown and a Q
         * equation (added by the MultiSlack NR extension), i.e. every slack bus
         * whose magnitude is NOT pinned by a local voltage-regulating generator.
         * This covers distributed-slack participants whose generator is PQ
         * (voltage_regulator_on == false), slack buses hosting a remote-voltage
         * controller, and SVC-regulated slack buses. A slack bus that DOES host a
         * local PV generator stays Vm-fixed (PV-like) with no Q equation. AC
         * labelling. Only valid once `pre_process_solver` ran (it needs
         * `id_me_to_ac_solver_` / `slack_bus_id_ac_solver_`).
         */
        std::set<int> get_free_vm_slack_solver_buses() const;
        /**
         * GRID-bus ids whose voltage magnitude is set by a VoltageControl group
         * rather than by the classical PV treatment. A bus lands here as soon as
         * an ACTIVE remote-regulating generator or an active voltage-mode SVC aims
         * at it, because the group's bordered voltage row needs that bus to keep a
         * Vm unknown (and hence a Q equation).
         *
         * Buses regulated ONLY by things standing on them -- local generators, or
         * voltage-regulating hvdc converter stations -- are deliberately NOT
         * included: several machines sharing a bus and all regulating it locally is
         * the ordinary PV case, handled exactly as before by the per-bus reactive
         * redistribution (`GeneratorContainer::set_q`). It is the arrival of a
         * remote controller (or an SVC, which is always a group controller) that
         * switches the bus over to the bordered formulation -- and then every local
         * regulator on that bus, generator or converter station, joins the group as
         * a co-controller instead of pinning it.
         *
         * Works in grid ids and reads only input data, so unlike the solver-side
         * accessors it is valid before / independently of `pre_process_solver`
         * (`fillpv_pq` needs it while it is still building that very labelling).
         */
        std::set<int> get_group_controlled_buses() const;
        /**
         * Set the grid bus whose voltage generator `gen_id` regulates (remote
         * voltage control). `bus_id` == the generator's own bus restores ordinary
         * local PV control.
         */
        void set_gen_regulated_bus(int gen_id, int bus_id) {
            // bus_id is stored verbatim and later dereferenced during solve, so reject
            // an out-of-range bus here (gen_id is validated inside set_regulated_bus).
            if(bus_id < 0 || bus_id >= static_cast<int>(total_bus())){
                std::ostringstream exc_;
                exc_ << "LSGrid::set_gen_regulated_bus: regulated bus id should be >= 0 and < ";
                exc_ << total_bus() << " (number of buses on the grid), you provided " << bus_id << ".";
                throw std::out_of_range(exc_.str());
            }
            generators_.set_regulated_bus(gen_id, bus_id, algo_controler_);
        }
        /**
         * Change the bus on the dc line "side 1" dcline_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus1_dcline(int dcline_id, GridModelBusId new_gridmodel_bus_id) {
            hvdc_lines_.change_bus_side_1(dcline_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus1_dcline_python(int dcline_id, int new_gridmodel_bus_id) {
            change_bus1_dcline(dcline_id, GridModelBusId(new_gridmodel_bus_id)); 
        }
        /**
         * Change the bus on the dc line "side 2" dcline_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */
        void change_bus2_dcline(int dcline_id, GridModelBusId new_gridmodel_bus_id) {
            hvdc_lines_.change_bus_side_2(dcline_id, new_gridmodel_bus_id, algo_controler_, substations_);
        }
        void change_bus2_dcline_python(int dcline_id, int new_gridmodel_bus_id) {
            change_bus2_dcline(dcline_id, GridModelBusId(new_gridmodel_bus_id)); 
        }
        [[nodiscard]] int get_bus1_dcline(int dcline_id) const {return hvdc_lines_.get_bus_side_1(dcline_id).cast_int();}
        [[nodiscard]] int get_bus2_dcline(int dcline_id) const {return hvdc_lines_.get_bus_side_2(dcline_id).cast_int();}

        // All results access
        [[nodiscard]] tuple3d get_loads_res() const {return loads_.get_res();}
        [[nodiscard]] const std::vector<bool>& get_loads_status() const { return loads_.get_status();}
        [[nodiscard]] tuple3d get_shunts_res() const {return shunts_.get_res();}
        [[nodiscard]] const std::vector<bool>& get_shunts_status() const { return shunts_.get_status();}
        [[nodiscard]] tuple3d get_gen_res() const {return generators_.get_res();}
        [[nodiscard]] const std::vector<bool>& get_gen_status() const { return generators_.get_status();}
        [[nodiscard]] tuple4d get_line_res1() const {return powerlines_.get_res_side_1();}
        [[nodiscard]] tuple4d get_line_res2() const {return powerlines_.get_res_side_2();}
        [[nodiscard]] const std::vector<bool>& get_lines_status() const { return powerlines_.get_status_global();}
        // per-side status (relevant for half-open lines, see `set_synch_status_both_side` /
        // `keep_half_open_lines`): `get_lines_status()` is the *global* status (both sides
        // disconnected), these report each side independently.
        [[nodiscard]] const std::vector<bool>& get_lines_status_side1() const { return powerlines_.get_status_side_1();}
        [[nodiscard]] const std::vector<bool>& get_lines_status_side2() const { return powerlines_.get_status_side_2();}
        [[nodiscard]] tuple4d get_trafo_res1() const {return trafos_.get_res_side_1();}
        [[nodiscard]] tuple4d get_trafo_res2() const {return trafos_.get_res_side_2();}
        [[nodiscard]] const std::vector<bool>& get_trafo_status() const { return trafos_.get_status_global();}
        [[nodiscard]] const std::vector<bool>& get_trafo_status_side1() const { return trafos_.get_status_side_1();}
        [[nodiscard]] const std::vector<bool>& get_trafo_status_side2() const { return trafos_.get_status_side_2();}
        [[nodiscard]] tuple3d get_storages_res() const {return storages_.get_res();}
        [[nodiscard]] const std::vector<bool>& get_storages_status() const { return storages_.get_status();}
        [[nodiscard]] tuple3d get_sgens_res() const {return sgens_.get_res();}
        [[nodiscard]] const std::vector<bool>& get_sgens_status() const { return sgens_.get_status();}
        [[nodiscard]] tuple3d get_dcline_res1() const {return hvdc_lines_.get_res_side_1();}
        [[nodiscard]] tuple3d get_dcline_res2() const {return hvdc_lines_.get_res_side_2();}
        [[nodiscard]] const std::vector<bool>& get_dclines_status() const { return hvdc_lines_.get_status_global();}

        [[nodiscard]] Eigen::Ref<const RealVect> get_gen_theta() const  {return generators_.get_theta();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_load_theta() const  {return loads_.get_theta();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_shunt_theta() const  {return shunts_.get_theta();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_storage_theta() const  {return storages_.get_theta();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_line_theta1() const {return powerlines_.get_theta_side_1();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_line_theta2() const {return powerlines_.get_theta_side_2();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_trafo_theta1() const {return trafos_.get_theta_side_1();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_trafo_theta2() const {return trafos_.get_theta_side_2();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_dcline_theta1() const {return hvdc_lines_.get_theta_side_1();}
        [[nodiscard]] Eigen::Ref<const RealVect> get_dcline_theta2() const {return hvdc_lines_.get_theta_side_2();}

        [[nodiscard]] const GlobalBusIdVect & get_all_shunt_buses() const {return shunts_.get_buses();}
        [[nodiscard]] Eigen::Ref<const IntVect> get_all_shunt_buses_numpy() const {return shunts_.get_bus_id_numpy();}

        [[nodiscard]] Eigen::Ref<const RealVect>  get_shunt_target_p() const {return shunts_.get_target_p();}
        [[nodiscard]] Eigen::Ref<const RealVect>  get_load_target_p() const {return loads_.get_target_p();}
        [[nodiscard]] Eigen::Ref<const RealVect>  get_gen_target_p() const {return generators_.get_target_p();}
        [[nodiscard]] Eigen::Ref<const RealVect>  get_sgen_target_p() const {return sgens_.get_target_p();}
        [[nodiscard]] Eigen::Ref<const RealVect>  get_storage_target_p() const {return storages_.get_target_p();}

        // complete results (with theta)
        [[nodiscard]] tuple4d get_loads_res_full() const {return loads_.get_res_full();}
        [[nodiscard]] tuple4d get_shunts_res_full() const {return shunts_.get_res_full();}
        [[nodiscard]] tuple4d get_gen_res_full() const {return generators_.get_res_full();}
        [[nodiscard]] tuple5d get_line_res1_full() const {return powerlines_.get_res_full_side_1();}
        [[nodiscard]] tuple5d get_line_res2_full() const {return powerlines_.get_res_full_side_2();}
        [[nodiscard]] tuple5d get_trafo_res1_full() const {return trafos_.get_res_full_side_1();}
        [[nodiscard]] tuple5d get_trafo_res2_full() const {return trafos_.get_res_full_side_2();}
        [[nodiscard]] tuple4d get_storages_res_full() const {return storages_.get_res_full();}
        [[nodiscard]] tuple4d get_sgens_res_full() const {return sgens_.get_res_full();}
        [[nodiscard]] tuple4d get_dcline_res1_full() const {return hvdc_lines_.get_res_full_side_1();}
        [[nodiscard]] tuple4d get_dcline_res2_full() const {return hvdc_lines_.get_res_full_side_2();}

        /**
         * @brief Get the Ybus solver object (AC)
         * 
         * This function allows to retrieve the Ybus passed to the AC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_Ybus` before that)
         * 
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        Eigen::SparseMatrix<cplx_type> get_Ybus_solver() const {
            return Ybus_ac_;  // This is copied to python
        }

        /**
         * @brief Get the Ybus solver object (DC)
         * 
         * This function allows to retrieve the Ybus passed to the DC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_dcYbus` before that)
         * 
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        Eigen::SparseMatrix<real_type> get_dcYbus_solver() const {
            return Bbus_dc_;  // This is copied to python (DC admittance matrix is real)
        }

        /**
         * @brief Get the Sbus solver object (AC)
         * 
         * This function allows to retrieve the Sbus passed to the AC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_Sbus` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] Eigen::Ref<const CplxVect> get_Sbus_solver() const{
            return acSbus_;
        }

        /**
         * @brief Get the Sbus solver object (DC)
         * 
         * This function allows to retrieve the Sbus passed to the DC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_dcSbus` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] Eigen::Ref<const RealVect> get_dcSbus_solver() const{
            return dcPbus_;  // DC power injection is real (active power P)
        }

        /**
         * @brief Get the Ybus gridmodel object (AC powerflow)
         * 
         * This function allows to retrieve the Ybus as seen by the gridmodel.
         * 
         * It may contain empty rows / columns for disconnected buses.
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of active buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Ybus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] Eigen::SparseMatrix<cplx_type> get_Ybus() const {
            return _relabel_matrix(Ybus_ac_, id_ac_solver_to_me_);
        }

        /**
         * @brief Get the Ybus gridmodel object (DC powerflow)
         * 
         * This function allows to retrieve the Ybus (DC) as seen by the gridmodel.
         * 
         * It may contain empty rows / columns for disconnected buses.
         * 
         * It has the "grimodel bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of active buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_dcYbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] Eigen::SparseMatrix<real_type> get_dcYbus() const {
            return _relabel_matrix(Bbus_dc_, id_dc_solver_to_me_);
        }

        /**
         * @brief Get the Sbus solver object (AC)
         * 
         * This function allows to retrieve the Sbus as represented by the "gridmodel"
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "solver" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of total buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Sbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] CplxVect get_Sbus() const {
            return _relabel_vector(acSbus_, id_ac_solver_to_me_);
        }

        /**
         * @brief Get the Sbus solver object (DC)
         * 
         * This function allows to retrieve the Sbus as represented by the "gridmodel"
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "solver" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of total buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Sbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        [[nodiscard]] RealVect get_dcSbus() const {
            return _relabel_vector(dcPbus_, id_dc_solver_to_me_);
        }

        /**
         * @brief Get the vector (list) of pv buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        [[nodiscard]] const SolverBusIdVect & get_pv_solver() const{
            return _has_ac_cache() ? bus_pv_ac_ : bus_pv_dc_;
        }
        /**
         * @brief Get the vector (list) of pv buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        [[nodiscard]] Eigen::Ref<const IntVect> get_pv_solver_numpy() const{
            return get_pv_solver().as_eigen();  // was _to_intvect()
        }
        // Same, but for one explicitly-named solver family. The AC and the DC
        // solver each build their own pv / pq split and their own slack weights
        // (they do not share a bus labelling, and either may be one powerflow
        // behind the other), so ask for the one you mean: the family-less
        // accessors above answer for AC when an AC powerflow has run, and only
        // fall back to DC otherwise.
        [[nodiscard]] const SolverBusIdVect & get_ac_pv_solver() const {return bus_pv_ac_;}
        [[nodiscard]] const SolverBusIdVect & get_dc_pv_solver() const {return bus_pv_dc_;}
        [[nodiscard]] Eigen::Ref<const IntVect> get_ac_pv_solver_numpy() const {return bus_pv_ac_.as_eigen();}
        [[nodiscard]] Eigen::Ref<const IntVect> get_dc_pv_solver_numpy() const {return bus_pv_dc_.as_eigen();}

        /**
         * @brief Get the vector (list) of pv buses, gridmodel labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return const Eigen::VectorXi
         */
        [[nodiscard]] GlobalBusIdVect get_pv() const{
            if(_has_ac_cache()) return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(bus_pv_ac_, id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(bus_pv_dc_, id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_pv: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }
        [[nodiscard]] IntVect get_pv_numpy() const{
            return get_pv().as_eigen();  // was _to_intvect()
        }

        /**
         * @brief Get the vector (list) of pq buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        [[nodiscard]] const SolverBusIdVect & get_pq_solver() const{
            return _has_ac_cache() ? bus_pq_ac_ : bus_pq_dc_;
        }
        /**
         * @brief Get the vector (list) of pq buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        [[nodiscard]] const Eigen::Ref<const IntVect> get_pq_solver_numpy() const{
            return get_pq_solver().as_eigen();  // was _to_intvect()
        }
        // per-family variants, see get_ac_pv_solver / get_dc_pv_solver above
        [[nodiscard]] const SolverBusIdVect & get_ac_pq_solver() const {return bus_pq_ac_;}
        [[nodiscard]] const SolverBusIdVect & get_dc_pq_solver() const {return bus_pq_dc_;}
        [[nodiscard]] Eigen::Ref<const IntVect> get_ac_pq_solver_numpy() const {return bus_pq_ac_.as_eigen();}
        [[nodiscard]] Eigen::Ref<const IntVect> get_dc_pq_solver_numpy() const {return bus_pq_dc_.as_eigen();}

        /**
         * @brief Get the vector (list) of pq buses, grimodel labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return const Eigen::VectorXi
         */
        [[nodiscard]] GlobalBusIdVect get_pq() const{
            if(_has_ac_cache()) return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(bus_pq_ac_, id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(bus_pq_dc_, id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_pq: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }
        [[nodiscard]] IntVect get_pq_numpy() const{
            return get_pq().as_eigen();  // was _to_intvect()
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (AC), solver labelling
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        [[nodiscard]] const SolverBusIdVect & get_slack_ids_solver() const{
            return slack_bus_id_ac_solver_;
        }
        /**
         * @brief Get the vector (list) of pq buses, grimodel labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return const Eigen::VectorXi
         */
        [[nodiscard]] Eigen::Ref<const IntVect> get_slack_ids_solver_numpy() const{
            return slack_bus_id_ac_solver_.as_eigen();  // was _to_intvect()
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (AC), gridmodel labelling
         * 
         * @return const Eigen::VectorXi 
         */
        [[nodiscard]] GlobalBusIdVect get_slack_ids() const {
            return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(slack_bus_id_ac_solver_, id_ac_solver_to_me_);
        }
        [[nodiscard]] IntVect get_slack_ids_numpy() const {
            return get_slack_ids().as_eigen();  // was _to_intvect()
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (DC), solver labelling
         * 
         * @return const SolverBusIdVect &
         */
        [[nodiscard]] const SolverBusIdVect & get_slack_ids_dc_solver() const{
            return slack_bus_id_dc_solver_;
        }
        /**
         * @brief Get the ids of the buses that participate to the slack (DC), solver labelling
         * 
         * @return Eigen::Ref<const IntVect> 
         */
        [[nodiscard]] Eigen::Ref<const IntVect> get_slack_ids_dc_solver_numpy() const{
            return slack_bus_id_dc_solver_.as_eigen();  // was _to_intvect()
        }

        [[nodiscard]] GlobalBusIdVect get_slack_ids_dc() const{
            return _relabel_vector2<SolverBusIdVect, GlobalBusIdVect>(
                slack_bus_id_dc_solver_,
                id_dc_solver_to_me_);
        }
        /**
         * @brief Get the ids of the buses that participate to the slack (DC), gridmodel labelling
         * 
         * @return const Eigen::VectorXi 
         */
        [[nodiscard]] IntVect get_slack_ids_dc_numpy() const{
            return get_slack_ids_dc().as_eigen();  // was _to_intvect()
        }

        /**
         * @brief Get the slack weights for each buses (solver labelling)
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        [[nodiscard]] Eigen::Ref<const RealVect> get_slack_weights_solver() const{
            return _has_ac_cache() ? slack_weights_ac_ : slack_weights_dc_;
        }
        // per-family variants, see get_ac_pv_solver / get_dc_pv_solver above
        [[nodiscard]] Eigen::Ref<const RealVect> get_ac_slack_weights_solver() const {return slack_weights_ac_;}
        [[nodiscard]] Eigen::Ref<const RealVect> get_dc_slack_weights_solver() const {return slack_weights_dc_;}

        /**
         * @brief Get the slack weights for each buses (gridmodel labelling)
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        [[nodiscard]] RealVect get_slack_weights() const{
            if(_has_ac_cache()) return _relabel_vector(slack_weights_ac_, id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(slack_weights_dc_, id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_slack_weights: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (complex) voltage angles for each buses (solver labelling)
         * 
         * @return Eigen::Ref<const CplxVect> 
         */
        [[nodiscard]] Eigen::Ref<const CplxVect> get_V_solver() const{
            return _algo.get_V();
        }

        /**
         * @brief Get the (complex) voltage angles for each buses (grimodel labelling)
         * 
         * @return CplxVect
         */
        [[nodiscard]] CplxVect get_V() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_V_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_V_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_V: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (real) voltage angle for each buses of the grid (solver labelling)
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        [[nodiscard]] Eigen::Ref<const RealVect> get_Va_solver() const{
            return _algo.get_Va();
        }

        /**
         * @brief Get the (real) voltage angle for each buses of the grid (grimodel labelling)
         * 
         * @return const RealVect
         */
        [[nodiscard]] RealVect get_Va() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_Va_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_Va_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_Va: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (real) voltage magnitude for each buses of the grid (solver labelling)
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        [[nodiscard]] Eigen::Ref<const RealVect> get_Vm_solver() const{
            return _algo.get_Vm();
        }

        /**
         * @brief Get the (real) voltage magnitude for each buses of the grid (grimodel labelling)
         * 
         * @return const RealVect
         */
        [[nodiscard]] RealVect get_Vm() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_Vm_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_Vm_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("LSGrid::get_Vm: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        [[nodiscard]] Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J_solver() const{
            return _algo.get_J();
        }

        /**
         * @brief Returns the last computed jacobian matrix (solver labelling)
         * 
         * @return Eigen::SparseMatrix<real_type> 
         */
        [[nodiscard]] Eigen::SparseMatrix<real_type> get_J_python_solver() const{
            return _algo.get_J_python();  // This is copied to python
        }

        // Ledger maps of the augmented NR Jacobian, in solver bus numbering
        // (size n_bus, -1 when the bus owns no such row/column). They describe
        // the layout of get_J_solver() and let external batched solvers (e.g.
        // gpusim2grid) rebuild the dS scatter / residual layout without
        // re-deriving it from Ybus. Only valid for Newton-Raphson algorithms
        // after a solve; throw otherwise (via the active algo).
        [[nodiscard]] IntVect get_theta_to_J_col_solver() const { return _algo.get_theta_to_J_col_python(); }
        [[nodiscard]] IntVect get_vm_to_J_col_solver()    const { return _algo.get_vm_to_J_col_python(); }
        [[nodiscard]] IntVect get_q_to_J_col_solver()     const { return _algo.get_q_to_J_col_python(); }
        [[nodiscard]] IntVect get_p_to_J_row_solver()     const { return _algo.get_p_to_J_row_python(); }
        [[nodiscard]] IntVect get_q_to_J_row_solver()     const { return _algo.get_q_to_J_row_python(); }

        // Compact (bus, row/col) registration pair lists -- the row/col
        // counterpart of the *_to_J_col / *_to_J_row maps above, but preserving
        // EVERY registration in order (a bus may appear more than once, or be
        // absent from the bus-keyed map's current value if a later
        // registration shadowed it there -- see NRLedger's "Multiplicity
        // rules"). NRSystem::_residual() itself iterates these, not the
        // bus-keyed maps; external batched solvers must do the same to
        // reproduce every contribution.
        [[nodiscard]] IntVect get_p_buses_solver()     const { return _algo.get_p_buses_python(); }
        [[nodiscard]] IntVect get_p_rows_solver()      const { return _algo.get_p_rows_python(); }
        [[nodiscard]] IntVect get_q_buses_solver()     const { return _algo.get_q_buses_python(); }
        [[nodiscard]] IntVect get_q_rows_solver()      const { return _algo.get_q_rows_python(); }
        [[nodiscard]] IntVect get_theta_buses_solver() const { return _algo.get_theta_buses_python(); }
        [[nodiscard]] IntVect get_theta_cols_solver()  const { return _algo.get_theta_cols_python(); }
        [[nodiscard]] IntVect get_vm_buses_solver()    const { return _algo.get_vm_buses_python(); }
        [[nodiscard]] IntVect get_vm_cols_solver()     const { return _algo.get_vm_cols_python(); }
        // MultiSlack slack_absorbed J column (-1 when distributed slack inactive).
        [[nodiscard]] int     get_slack_col_solver()      const { return _algo.get_slack_col(); }
        // MultiSlack converged slack_absorbed VALUE (pu; 0 when inactive). This
        // is the ground-truth state after convergence -- NOT the 0 initial
        // guess an external solver's own linearized derivation starts from.
        [[nodiscard]] real_type get_slack_absorbed_solver() const { return _algo.get_slack_absorbed(); }
        // VoltageControl (gen / SVC / hvdc station) converged reactive injection per
        // controller (pu), + its (kind: 0=GEN,1=SVC; element id) identity, in
        // controller registration order (empty when the extension is inactive).
        // Ground truth for external solvers deriving their own controller_q.
        [[nodiscard]] RealVect get_controller_q_solver()       const { return _algo.get_controller_q(); }
        [[nodiscard]] IntVect  get_controller_kind_solver()    const { return _algo.get_controller_kind(); }
        [[nodiscard]] IntVect  get_controller_elem_id_solver() const { return _algo.get_controller_elem_id(); }
        // J column of each controller's own Q unknown, controller registration
        // order -- NOT the bus-keyed get_q_to_J_col_solver() (that map only keeps
        // the LAST controller registered at a given bus, see NRLedger::
        // add_q_unknown's own doc). External solvers rebuilding this bordered
        // block (e.g. gpusim2grid) MUST use this whenever two controllers can
        // share a bus, exactly like p_buses()/p_rows() etc. must be used instead
        // of the bus-keyed maps for the base P/Q block.
        [[nodiscard]] IntVect  get_controller_q_col_solver()   const { return _algo.get_controller_q_col(); }

        [[nodiscard]] real_type get_computation_time() const{ return _algo.get_computation_time();}
        [[nodiscard]] real_type get_dc_computation_time() const{ return _dc_algo.get_computation_time();}

        // part dedicated to grid2op backend, optimized for grid2op data representation (for speed)
        // this is not recommended to use it outside of its intended usage within grid2op !
        void update_gens_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                           const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        void update_sgens_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                           const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        void update_gens_v(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                           const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        void update_loads_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                            const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        void update_loads_q(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                            const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        /**
         * Update the topology based on the topology vector id.
         *
         * The new_values are given in LocalBusId (-1, 1, 2 etc.) and not in
         * SolverBusId nor GridModelBusId
         */
        void update_topo(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                         const Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values);
        void update_storages_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                               const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values);

        void set_load_pos_topo_vect(const Eigen::Ref<const IntVect> & load_pos_topo_vect)
        {
            loads_.set_pos_topo_vect(load_pos_topo_vect);
        }
        void set_gen_pos_topo_vect(const Eigen::Ref<const IntVect> & gen_pos_topo_vect)
        {
            generators_.set_pos_topo_vect(gen_pos_topo_vect);
        }
        void set_storage_pos_topo_vect(const Eigen::Ref<const IntVect> & sto_pos_topo_vect)
        {
            storages_.set_pos_topo_vect(sto_pos_topo_vect);
        }
        void set_line_pos1_topo_vect(const Eigen::Ref<const IntVect> & line_or_pos_topo_vect)
        {
            powerlines_.set_pos_topo_vect_side_1(line_or_pos_topo_vect);
        }
        void set_line_pos2_topo_vect(const Eigen::Ref<const IntVect> & line_ex_pos_topo_vect)
        {
            powerlines_.set_pos_topo_vect_side_2(line_ex_pos_topo_vect);
        }
        void set_trafo_pos1_topo_vect(const Eigen::Ref<const IntVect> & trafo_hv_pos_topo_vect)
        {
            trafos_.set_pos_topo_vect_side_1(trafo_hv_pos_topo_vect);
        }
        void set_trafo_pos2_topo_vect(const Eigen::Ref<const IntVect> & trafo_lv_pos_topo_vect)
        {
            trafos_.set_pos_topo_vect_side_2(trafo_lv_pos_topo_vect);
        }

        void set_load_to_subid(const Eigen::Ref<const IntVect> & load_to_subid)
        {
            loads_.set_subid(load_to_subid);
        }
        void set_gen_to_subid(const Eigen::Ref<const IntVect> & gen_to_subid)
        {
            generators_.set_subid(gen_to_subid);
        }
        void set_storage_to_subid(const Eigen::Ref<const IntVect> & storage_to_subid)
        {
            storages_.set_subid(storage_to_subid);
        }
        void set_shunt_to_subid(const Eigen::Ref<const IntVect> & shunt_to_subid)
        {
            shunts_.set_subid(shunt_to_subid);
        }
        void set_line_to_sub1_id(const Eigen::Ref<const IntVect> & line_or_to_subid)
        {
            powerlines_.set_subid_side_1(line_or_to_subid);
        }
        void set_line_to_sub2_id(const Eigen::Ref<const IntVect> & line_ex_to_subid)
        {
            powerlines_.set_subid_side_2(line_ex_to_subid);
        }
        void set_trafo_to_sub1_id(const Eigen::Ref<const IntVect> & trafo_hv_to_subid)
        {
            trafos_.set_subid_side_1(trafo_hv_to_subid);
        }
        void set_trafo_to_sub2_id(const Eigen::Ref<const IntVect> & trafo_lv_to_subid)
        {
            trafos_.set_subid_side_2(trafo_lv_to_subid);
        }
        void set_n_sub(int n_sub)
        {
            n_sub_ = n_sub;
        }
        int get_n_sub() const
        {
            return n_sub_;
        }
        void set_max_nb_bus_per_sub(int max_nb_bus_per_sub)
        {
            if(substations_.nb_bus() !=  static_cast<unsigned int>(n_sub_ * max_nb_bus_per_sub)){
                std::ostringstream exc_;
                exc_ << "LSGrid::set_max_nb_bus_per_sub: ";
                exc_ << "your model counts ";
                exc_ << substations_.nb_bus()  << " buses according to `substations_.nb_bus()` but ";
                exc_ << n_sub_ * max_nb_bus_per_sub_ << " according to n_sub_ * max_nb_bus_per_sub_.";
                exc_ << "Both should match: either reinit it with another call to `init_bus` or set properly the number of ";
                exc_ << "substations / buses per substations with `set_n_sub` / `set_max_nb_bus_per_sub`";
                throw std::runtime_error(exc_.str());
            }
            max_nb_bus_per_sub_ = max_nb_bus_per_sub;
        }
        [[nodiscard]] int get_max_nb_bus_per_sub() const { return max_nb_bus_per_sub_;}

        /**
         * Relevant kwargs the grid was built with (eg by `init_from_pypowsybl`), as a
         * string->string map (eg {"sort_index": "True", "buses_for_sub": "False"}).
         * Empty for a grid not built that way (eg pandapower/powermodels), or a
         * default-constructed one. Set by the Python-side converter, never read by any
         * C++ logic itself -- purely a way for downstream Python consumers (eg a
         * pypowsybl-shaped result view) to recover conversion-time settings that are
         * otherwise plain Python arguments lost after `init()` returns.
         */
        [[nodiscard]] const std::map<std::string, std::string> & get_init_kwargs() const { return init_kwargs_;}
        void set_init_kwargs(const std::map<std::string, std::string> & init_kwargs) { init_kwargs_ = init_kwargs;}

        /**
         * Fused-bus lookup, size `total_bus()` (empty if unset / not built with bus
         * fusion). For each ls bus id, gives the ls bus id of the "representative" bus
         * it was merged into by `fuse_zero_impedance_branches` (identity for a bus not
         * involved in any fusion). Set by the Python-side converter (`init_from_pypowsybl`,
         * see `_from_pypowsybl.py`); never read by any C++ powerflow logic -- purely so a
         * downstream Python consumer (eg `LightsimResultNetwork`) can recover the voltage
         * of a bus whose own lightsim bus ended up with no elements (and so disconnected)
         * after fusion, by redirecting to its representative, which carries the real
         * solved voltage.
         */
        [[nodiscard]] const IntVect & get_bus_fusion_rep() const { return _bus_fusion_rep;}
        void set_bus_fusion_rep(const Eigen::Ref<const IntVect> & bus_fusion_rep){
            if(bus_fusion_rep.size() != 0 && static_cast<size_t>(bus_fusion_rep.size()) != total_bus()){
                std::ostringstream exc_;
                exc_ << "LSGrid::set_bus_fusion_rep: the provided vector has size ";
                exc_ << bus_fusion_rep.size() << " but this grid counts " << total_bus() << " buses.";
                throw std::runtime_error(exc_.str());
            }
            _bus_fusion_rep = bus_fusion_rep;
        }

        void fillSbus_other(Eigen::Ref<CplxVect> res, bool ac, const SolverBusIdVect& id_me_to_solver){
            fillSbus_me(res, ac, id_me_to_solver);
        }

        /**
        computes Ybus_ and Sbus_. It has different flags to have more control on the purpose for this "computation"
        is_ac indicates if you want to perform and AC powerflow or a DC powerflow and reset_solver indicates
        if you will perform a powerflow after it or not. (usually put ``true`` here).

        This is use internally by ac_pf or dc_pf but also when doing batched solvers (*eg* TimeSeries or Contingency analysis)
        **/
        // `init_pv_vm_targets`: when true (the default, used by the real ac_pf / dc_pf
        // solves), voltage-controlled elements with no droop/slope (generators, HVDC
        // converters, zero-slope SVCs -- see SvcContainer::set_vm) have the proposed
        // voltage MAGNITUDE at their regulated bus snapped to their own target, a
        // reasonable NR-initialization heuristic. `check_solution` passes `false`: it
        // is testing a caller-supplied (eg OLF's) voltage as-is, and silently
        // overwriting it there defeats the whole point of the check -- even a tiny,
        // physically-correct gap between that voltage and the local target (the
        // regulator doing its job) can look like a large spurious power mismatch at a
        // strongly-meshed bus.
        // Sbus stays a plain reference (not Eigen::Ref): forwarded into
        // _pre_process_solver_impl's inj, which forwards into prepare_injection,
        // which reassigns/resizes it.
        // `slack_weights` / `bus_pv` / `bus_pq` complete the group: they used to be
        // taken from this grid's own members whatever the caller passed for the rest,
        // which meant a foreign builder silently wrote its pv-pq split and its slack
        // weights into the grid's cache while leaving the grid's Ybus / Sbus behind.
        // The whole group now travels together -- see the ownership check at the top
        // of _pre_process_solver_impl, and BaseBatchSolverSynch::prepare_solver_input_base
        // for the caller that owns its own.
        CplxVect pre_process_solver(const Eigen::Ref<const CplxVect> & Vinit,
                                    CplxVect & Sbus,
                                    Eigen::SparseMatrix<cplx_type> & Ybus,
                                    SolverBusIdVect & id_me_to_solver,
                                    GlobalBusIdVect & id_solver_to_me,
                                    GlobalBusIdVect & slack_bus_id_me,
                                    SolverBusIdVect & slack_bus_id_solver,
                                    RealVect & slack_weights,
                                    SolverBusIdVect & bus_pv,
                                    SolverBusIdVect & bus_pq,
                                    bool is_ac,
                                    const AlgoControl & solver_control,
                                    bool init_pv_vm_targets = true);

        // DC-specific pre processing: builds the real Bbus (admittance) matrix and the
        // real Pbus (active power) vector, reusing the shared bus-mapping helpers. Mirrors
        // pre_process_solver but keeps the whole DC path real (no complex Ybus / Sbus).
        // Pbus stays a plain reference (not Eigen::Ref): same resize-forwarding
        // reason as pre_process_solver's Sbus above.
        CplxVect pre_process_dc_solver(const Eigen::Ref<const CplxVect> & Vinit,
                                       RealVect & Pbus,
                                       Eigen::SparseMatrix<real_type> & Bbus,
                                       SolverBusIdVect & id_me_to_solver,
                                       GlobalBusIdVect & id_solver_to_me,
                                       GlobalBusIdVect & slack_bus_id_me,
                                       SolverBusIdVect & slack_bus_id_solver,
                                       RealVect & slack_weights,
                                       SolverBusIdVect & bus_pv,
                                       SolverBusIdVect & bus_pq,
                                       const AlgoControl & solver_control);

        //for FDPF
        void fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, 
                        Eigen::SparseMatrix<real_type> & Bpp, 
                        FDPFMethod xb_or_bx) const;
        void init_fdpf_coeffs(){
            powerlines_.init_fdpf_coeffs();
            trafos_.init_fdpf_coeffs();
        }
        void fillBf_for_PTDF(Eigen::SparseMatrix<real_type> & Bf, bool transpose=false) const;

        Eigen::SparseMatrix<real_type> debug_get_Bp_python(FDPFMethod xb_or_bx){
            Eigen::SparseMatrix<real_type> Bp;
            Eigen::SparseMatrix<real_type> Bpp;
            init_fdpf_coeffs();
            fillBp_Bpp(Bp, Bpp, xb_or_bx);
            return Bp;
        }
        Eigen::SparseMatrix<real_type> debug_get_Bpp_python(FDPFMethod xb_or_bx){
            Eigen::SparseMatrix<real_type> Bp;
            Eigen::SparseMatrix<real_type> Bpp;
            init_fdpf_coeffs();
            fillBp_Bpp(Bp, Bpp, xb_or_bx);
            return Bpp;
        }

    protected:
        // Set both _ls_to_orig and _orig_to_ls.
        //
        // NOT noexcept (it used to be declared so, wrongly): it assigns one Eigen
        // vector and allocates another, either of which throws std::bad_alloc on
        // failure -- and a throw out of a noexcept function is not an error the
        // caller can handle, it is an immediate std::terminate(). Nothing needs the
        // guarantee (the only callers are the copy constructor and set_ls_to_orig,
        // neither of them noexcept), so the honest declaration is the plain one.
        //
        // It indexes _orig_to_ls with the *values* of `ls_to_orig`: only ever call
        // it on a vector check_ls_to_orig_values() has already accepted, as
        // set_ls_to_orig does.
        void set_ls_to_orig_internal(const Eigen::Ref<const IntVect> & ls_to_orig);
        // throws std::out_of_range unless every entry is -1 or a sane, non-negative
        // original-grid bus id -- see the definition in LSGrid.cpp for why.
        static void check_ls_to_orig_values(const Eigen::Ref<const IntVect> & ls_to_orig);

        // init the Ybus matrix (its size, it is filled up elsewhere) and also the 
        // converter from "my bus id" to the "solver bus id" (id_me_to_solver and id_solver_to_me)
        void init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                       int nb_bus_solver);
        void init_Bbus(Eigen::SparseMatrix<real_type> & Bbus,
                       int nb_bus_solver);  // DC: real admittance matrix
        void init_converter_bus_id(SolverBusIdVect& id_me_to_solver,
                                   GlobalBusIdVect& id_solver_to_me);

        // converts the slack_bus_id from gridmodel ordering into solver ordering
        void init_slack_bus(const SolverBusIdVect & id_me_to_solver,
                            const GlobalBusIdVect& id_solver_to_me,
                            const GlobalBusIdVect & slack_bus_id_me,
                            SolverBusIdVect & slack_bus_id_solver
                        );

        /**
         * @brief Build the result matrix (eg Ybus) (labelled using the gridmodel) 
         * from the input same matrix (eg Ybus) but labelled with the solver convention
         * 
         * @param Ybus : solver labelling
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @param relabel_row : whether to relabel also the row id
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        // Ybus stays a plain reference: T is deduced from the caller's argument, and
        // Eigen::Ref<const Eigen::SparseMatrix<T>> is a non-deduced context (callers
        // pass a concrete Eigen::SparseMatrix<T>, so deduction would fail).
        template<typename T>
        Eigen::SparseMatrix<T> _relabel_matrix(const Eigen::SparseMatrix<T> & Ybus,
                                               const GlobalBusIdVect & id_solver_to_me,
                                               bool relabel_row=true) const {
            // TODO optim : if relabel_row is false, then we can just copy
            // paste the columns easily in the target matrix, which should be
            // way faster than this function.
            using index_type = typename Eigen::SparseMatrix<T>::StorageIndex ;
            const int nb_conn_bus = nb_connected_bus();
            if(id_solver_to_me.size() == 0) throw std::runtime_error("LSGrid::_relabel_matrix: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Ybus.cols() != nb_conn_bus) throw std::runtime_error("LSGrid::_relabel_matrix: impossible to retrieve the `gridmodel`: the input matrix has not the right number of columns, (.., nb connected bus) expected");
            if(relabel_row & (Ybus.rows() != nb_conn_bus)) throw std::runtime_error("LSGrid::_relabel_matrix: impossible to retrieve the `gridmodel`: the input matrix has not the right number of columnd (nb connected bus, ...) expected");
            Eigen::SparseMatrix<T> res(relabel_row ? total_bus() : Ybus.rows(), total_bus());
            res.reserve(Ybus.nonZeros());
            std::vector<Eigen::Triplet<T> > tripletList;
            tripletList.reserve(Ybus.nonZeros());
            const auto n_col = Ybus.cols();
            for (Eigen::Index col_=0; col_ < n_col; ++col_){
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(Ybus, col_); it; ++it)
                {
                    if(relabel_row) tripletList.push_back({static_cast<index_type>(id_solver_to_me[static_cast<size_t>(it.row())]),
                                                           static_cast<index_type>(id_solver_to_me[static_cast<size_t>(it.col())]),
                                                           it.value()});
                    else tripletList.push_back({static_cast<index_type>(it.row()), 
                                                static_cast<index_type>(id_solver_to_me[static_cast<size_t>(it.col())]),
                                                it.value()});
                }
            }
            res.setFromTriplets(tripletList.begin(), tripletList.end());
            res.makeCompressed();
            return res;
        }

        /**
         * @brief Build the Sbus (or any other vector labelled using the gridmodel convention)
         * from the same vector (input) that uses the solver convention.
         *
         * Overloaded on purpose (not "copy paste, find a better way"): callers pass either
         * a genuine Eigen::Matrix<T,...> lvalue (eg `acSbus_`) or an Eigen::Ref<const
         * Matrix<T,...>> prvalue returned by a `get_..._solver()` accessor (eg
         * `get_V_solver()`). Both need their own overload because T only appears deduced
         * from the argument's own type: a Ref argument cannot deduce T through this
         * plain-Matrix parameter (Ref isn't-a Matrix), and conversely a genuine Matrix
         * argument cannot deduce T through the other overload's Eigen::Ref<const
         * Matrix<T,...>> parameter (a non-deduced context) -- confirmed by deleting this
         * overload and observing every `get_..._solver()`-fed call site fail to compile.
         *
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1> > & Sbus,
                                                            const GlobalBusIdVect & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("LSGrid::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Sbus.size() != nb_connected_bus()) throw std::runtime_error("LSGrid::_relabel_vector: impossible to retrieve the `gridmodel` input solver has not the right size, expected (nb connected bus, ).");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(total_bus());
            for(auto solver_id = 0; solver_id < Sbus.size(); ++solver_id){
                res[id_solver_to_me[solver_id].cast_int()] = Sbus[solver_id];
            }
            return res;
        }

        /**
         * @brief Build the Sbus (or any other vector labelled using the gridmodel convention)
         * from the same vector (input) that uses the solver convention.
         *
         * Sibling overload of the one above, for callers passing a genuine Eigen::Matrix<T,...>
         * lvalue instead of an Eigen::Ref -- see the comment above for why both are needed.
         * Kept as a plain reference (not Eigen::Ref): T is deduced from this argument, and
         * Eigen::Ref<const Eigen::Matrix<T,...>> is a non-deduced context here.
         *
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Sbus,
                                                            const GlobalBusIdVect & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("LSGrid::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Sbus.size() != nb_connected_bus()) throw std::runtime_error("LSGrid::_relabel_vector: impossible to retrieve the `gridmodel` input solver has not the right size, expected (nb connected bus, ).");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(total_bus());
            for(int solver_id = 0; solver_id < Sbus.size(); ++solver_id){
                
                res[id_solver_to_me[solver_id].cast_int()] = Sbus[solver_id];
            }
            return res;
        }

        /**
         * @brief Build the pv; pq or slack ids 
         * (or any other vector labelled using the gridmodel convention) 
         * from the same vector (input) that uses the solver convention.
         */
        template<class InputBusIdVect, class OutputBusIdVect>
        OutputBusIdVect _relabel_vector2(
            const InputBusIdVect & solver_indexed_bus,
            const GlobalBusIdVect & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("LSGrid::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            OutputBusIdVect res = OutputBusIdVect(solver_indexed_bus.size(), static_cast<typename OutputBusIdVect::value_type>(-1));  // TODO or 0...
            size_t pos_id = 0;
            for(const auto & el_id : solver_indexed_bus){
                res[pos_id] = id_solver_to_me[el_id.cast_int()];
                ++ pos_id;
            }
            return res;
        }
    
    protected:
        // Shared implementation of pre_process_solver (AC: complex Ybus / Sbus) and
        // pre_process_dc_solver (DC: real Bbus / Pbus). The matrix scalar type selects the
        // family (cplx_type => AC, real_type => DC); the type-specific steps (matrix init / fill,
        // injection assembly) are tag-dispatched to the overloads just below.
        // inj stays a plain reference (not Eigen::Ref): it is forwarded straight into
        // prepare_injection (CplxVect&/RealVect& overload), which reassigns/resizes it
        // -- not a template-deduction issue here (MatScalar/InjVect are always given
        // explicitly at the two call sites below), just resize-forwarding.
        template<class MatScalar, class InjVect>
        CplxVect _pre_process_solver_impl(const Eigen::Ref<const CplxVect> & Vinit,
                                          InjVect & inj,
                                          Eigen::SparseMatrix<MatScalar> & mat,
                                          SolverBusIdVect & id_me_to_solver,
                                          GlobalBusIdVect & id_solver_to_me,
                                          GlobalBusIdVect & slack_bus_id_me,
                                          SolverBusIdVect & slack_bus_id_solver,
                                          RealVect & slack_weights,
                                          SolverBusIdVect & bus_pv,
                                          SolverBusIdVect & bus_pq,
                                          const AlgoControl & solver_control,
                                          bool init_pv_vm_targets);

        // How many of the nine solver-side containers `_pre_process_solver_impl` was
        // handed are this LSGrid's own cache for that family. Nine means the caller is
        // the grid itself (ac_pf / dc_pf / check_solution, which pass their members);
        // zero means a foreign builder (the batch algorithms, which own theirs). Any
        // number in between is a programming error -- see the call site.
        // Everything is compared as `const void *` so the two families' differently
        // typed members can be handled by one runtime branch (C++14: no `if constexpr`).
        // Nine pointer comparisons summed branchlessly, once per powerflow.
        [[nodiscard]] int _nb_own_cache_args(bool is_ac,
                                             const void * inj,
                                             const void * mat,
                                             const void * id_me_to_solver,
                                             const void * id_solver_to_me,
                                             const void * slack_bus_id_me,
                                             const void * slack_bus_id_solver,
                                             const void * slack_weights,
                                             const void * bus_pv,
                                             const void * bus_pq) const noexcept {
            if(is_ac){
                return (inj                 == (const void *) &acSbus_)
                     + (mat                 == (const void *) &Ybus_ac_)
                     + (id_me_to_solver     == (const void *) &id_me_to_ac_solver_)
                     + (id_solver_to_me     == (const void *) &id_ac_solver_to_me_)
                     + (slack_bus_id_me     == (const void *) &slack_bus_id_ac_me_)
                     + (slack_bus_id_solver == (const void *) &slack_bus_id_ac_solver_)
                     + (slack_weights       == (const void *) &slack_weights_ac_)
                     + (bus_pv              == (const void *) &bus_pv_ac_)
                     + (bus_pq              == (const void *) &bus_pq_ac_);
            }
            return (inj                 == (const void *) &dcPbus_)
                 + (mat                 == (const void *) &Bbus_dc_)
                 + (id_me_to_solver     == (const void *) &id_me_to_dc_solver_)
                 + (id_solver_to_me     == (const void *) &id_dc_solver_to_me_)
                 + (slack_bus_id_me     == (const void *) &slack_bus_id_dc_me_)
                 + (slack_bus_id_solver == (const void *) &slack_bus_id_dc_solver_)
                 + (slack_weights       == (const void *) &slack_weights_dc_)
                 + (bus_pv              == (const void *) &bus_pv_dc_)
                 + (bus_pq              == (const void *) &bus_pq_dc_);
        }
        // matrix (re)initialization, overloaded per family (no `if constexpr`, C++14)
        void init_solver_matrix(Eigen::SparseMatrix<cplx_type> & mat, int nb_bus_solver){ init_Ybus(mat, nb_bus_solver); }
        void init_solver_matrix(Eigen::SparseMatrix<real_type> & mat, int nb_bus_solver){ init_Bbus(mat, nb_bus_solver); }
        void fill_solver_matrix(Eigen::SparseMatrix<cplx_type> & mat, const SolverBusIdVect & id_me_to_solver){ fillYbus(mat, true, id_me_to_solver); }
        void fill_solver_matrix(Eigen::SparseMatrix<real_type> & mat, const SolverBusIdVect & id_me_to_solver){ fillBdc(mat, id_me_to_solver); }
        // injection assembly, overloaded per family (AC fills complex Sbus + reactive vectors; DC fills real Pbus)
        void prepare_injection(CplxVect & Sbus, bool redo_all, bool converter_changed,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const AlgoControl & solver_control);
        void prepare_injection(RealVect & Pbus, bool redo_all, bool converter_changed,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const AlgoControl & solver_control);

        void fillYbus(Eigen::SparseMatrix<cplx_type> & res, bool ac, const SolverBusIdVect& id_me_to_solver);
        void fillBdc(Eigen::SparseMatrix<real_type> & res, const SolverBusIdVect& id_me_to_solver);  // DC: real admittance matrix
        void fillSbus_me(Eigen::Ref<CplxVect> res, bool ac, const SolverBusIdVect& id_me_to_solver);
        // writes the pv / pq split into the caller-supplied vectors: the AC and the
        // DC family each own theirs (bus_pv_ac_ / bus_pv_dc_, ...), and they are
        // expressed in that family's own solver labelling.
        void fillpv_pq(const SolverBusIdVect& id_me_to_solver,
                       const GlobalBusIdVect& id_solver_to_me,
                       const SolverBusIdVect & slack_bus_id_solver,
                       SolverBusIdVect & bus_pv_out,
                       SolverBusIdVect & bus_pq_out);

        // results
        /**process the results from the solver to this instance
        **/
        // res stays a plain reference (not Eigen::Ref): it is reassigned below
        // (res = _get_results_back_to_orig_nodes(...)) to total_bus() size, which can
        // differ from the caller's initial (often empty) res.
        void process_results(bool conv, CplxVect & res, const Eigen::Ref<const CplxVect> & Vinit, bool ac,
                             SolverBusIdVect & id_me_to_solver);

        // Sanity-check the voltages a solver returned before they are consumed.
        // A wrong-sized V/Va/Vm is a contract violation -> throws
        // std::runtime_error. Non-finite values -> the solve is marked as
        // non-converged (ErrorType::InifiniteValue) and this returns false, so no
        // NaN/Inf propagates and no out-of-bounds access happens downstream.
        // Returns true when the outputs are usable.
        // NB only called for external (plugin) solvers, ie those whose type is
        // AlgorithmType::Custom -- built-in solvers pay nothing for this.
        bool _check_solver_output(bool ac);

        /**
        Compute the results vector from the Va, Vm post powerflow
        **/
        void compute_results(bool ac);
        /**
        reset the results in case of divergence of the powerflow.
        **/
        void reset_results();

        /**
        reset the solver, and all its results
        **/
        void reset(bool reset_solver, bool reset_ac, bool reset_dc);

        /**
        optimization for grid2op
        **/
        template<class T>
        void update_continuous_values(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                                      const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
                                      T fun)
        {
            // new_values is indexed by has_changed's length below; a shorter new_values
            // would over-read the buffer (Eigen operator[] is unchecked in release).
            if(new_values.rows() != has_changed.rows()){
                std::ostringstream exc_;
                exc_ << "LSGrid::update_continuous_values: 'has_changed' (size " << has_changed.rows();
                exc_ << ") and 'new_values' (size " << new_values.rows() << ") must have the same size.";
                throw std::runtime_error(exc_.str());
            }
            for(int el_id = 0; el_id < has_changed.rows(); ++el_id)
            {
                if(has_changed(el_id))
                {
                    (this->*fun)(el_id, static_cast<real_type>(new_values[el_id]));
                }
            }
        }

        CplxVect _get_results_back_to_orig_nodes(const Eigen::Ref<const CplxVect> & res_tmp,
                                                 SolverBusIdVect & id_me_to_solver,
                                                 int size);

        void check_solution_q_values( Eigen::Ref<CplxVect> res, bool check_q_limits) const;
        void check_solution_q_values_onegen(Eigen::Ref<CplxVect> res,
                                            int bus_id,
                                            real_type min_q_mvar,
                                            real_type max_q_mvar,
                                            bool check_q_limits) const;

    private:
        // ---- cache-reuse bookkeeping (see the public API above) ---------------
        // Has an AC powerflow ever built the AC solver-side data on this grid?
        // Used by the family-less get_pv_solver() / get_pq_solver() /
        // get_slack_weights_solver() to pick which family to answer for.
        [[nodiscard]] bool _has_ac_cache() const noexcept {return id_ac_solver_to_me_.size() > 0;}

        // Record that `family_control`'s cached data now matches the grid: its
        // flags go back to "nothing changed" and its connectivity snapshot is
        // refreshed. Called after every powerflow of that family (converged or
        // not: pre_process built the data against the current grid either way),
        // unless that family's reuse was turned off.
        void _mark_cache_valid(bool ac){
            if(ac){
                algo_controler_.ac_algo_controler().tell_none_changed();
                last_bus_status_saved_ = substations_.get_bus_status();
            } else {
                algo_controler_.dc_algo_controler().tell_none_changed();
                last_bus_status_dc_ = substations_.get_bus_status();
            }
        }

        // One family's half of init_bus_status(): flag a dimension change iff the
        // connectivity differs from the one that family's cache was built with.
        static void _flag_dimension_change(const std::vector<bool> & new_status,
                                           const std::vector<bool> & last_status,
                                           AlgoControl & family_control){
            if(new_status.size() != last_status.size()){
                // the snapshot was never taken (fresh grid, or a family that has
                // not solved since a prevent_*_cache_reuse), or the grid changed
                // size: either way nothing may be reused
                family_control.tell_dimension_changed();
                return;
            }
            for(std::size_t global_bus = 0; global_bus < new_status.size(); ++global_bus){
                if(last_status[global_bus] != new_status[global_bus]){
                    family_control.tell_dimension_changed();
                    return;
                }
            }
        }

        // memory for the import
        // TODO switches: move to BaseSubstation
        /**
         * _ls_to_orig: has the size of the number of possible buses in lightsim2grid 
         * (*ie* `n_sub_ * max_nb_bus_per_sub_` ) and gives the id of the corresponding
         * bus in the original grid (pandapower or pypowsybl).
         * 
         * If a "-1" is present, then this bus does not exist in the original grid, 
         * it is only present in the lightsim2grid gridmodel.
         */
        IntVect _ls_to_orig; 
        /**
         * Opposite to _ls_to_orig. The vector _orig_to_ls has the size of the number
         * of buses in the original grid (pandapower or pypowsybl) and tells 
         * to which bus of lightsim2grid it corresponds. It should be a >= integer
         * between 0 and `n_sub_ * max_nb_bus_per_sub_`
         */
        IntVect _orig_to_ls;
        IntVect _bus_fusion_rep;  // see get_bus_fusion_rep()

        // member of the grid
        double timer_last_ac_pf_;
        double timer_last_dc_pf_;

        DualAlgoControl algo_controler_;  // independent change tracking for the AC and DC solver families
        bool compute_results_;
        real_type init_vm_pu_;  // default vm initialization, mainly for dc powerflow
        real_type sn_mva_;

        // powersystem representation
        // 1. bus
        int n_sub_;
        int max_nb_bus_per_sub_;
        SubstationContainer substations_;
        std::vector<bool> last_bus_status_saved_;
        std::map<std::string, std::string> init_kwargs_;  // see get_init_kwargs()

        // always have the length of the number of buses,
        // id_me_to_model_[id_me] gives -1 if the bus "id_me" is deactivated, or "id_model" if it is activated.
        SolverBusIdVect id_me_to_ac_solver_;
        // convert the bus id from the model to the bus id of me.
        // it has a variable size, that depends on the number of connected bus. if "id_model" is an id of a bus
        // sent to the solver, then id_model_to_me_[id_model] is the bus id of this model of the grid.
        GlobalBusIdVect id_ac_solver_to_me_;

        SolverBusIdVect id_me_to_dc_solver_;
        GlobalBusIdVect id_dc_solver_to_me_;

        // 2. powerline
        LineContainer powerlines_;

        // 3. shunt
        ShuntContainer shunts_;

        // 4. transformers
        // have the r, x, h and ratio
        // ratio is computed from the tap, so maybe store tap num and tap_step_pct
        TrafoContainer trafos_;

        // 5. generators
        RealVect total_q_min_per_bus_;  // TODO switches: move to BaseSubstation
        RealVect total_q_max_per_bus_;  // TODO switches: move to BaseSubstation
        Eigen::VectorXi total_gen_per_bus_;
        GeneratorContainer generators_;

        // 6. loads
        LoadContainer loads_;

        // 7. static generators (P,Q generators)
        SGenContainer sgens_;

        // 8. storage units
        StorageContainer storages_;

        // 9. hvdc (converter stations + lines, exposed through the "dcline" API)
        HvdcLineContainer hvdc_lines_;

        // 9b. static var compensators (SVC)
        SvcContainer svcs_;

        // 10. slack bus
        // std::vector<int> slack_bus_id_;
        GlobalBusIdVect slack_bus_id_ac_me_;  // slack bus id, gridmodel number
        SolverBusIdVect slack_bus_id_ac_solver_;  // slack bus id, solver number
        GlobalBusIdVect slack_bus_id_dc_me_;
        SolverBusIdVect slack_bus_id_dc_solver_;
        // AC family's slack weights (solver labelling). The DC twin is appended at
        // the end of the class -- see the ABI note there.
        RealVect slack_weights_ac_;

        // as matrix, for the solver
        Eigen::SparseMatrix<cplx_type> Ybus_ac_;
        Eigen::SparseMatrix<real_type> Bbus_dc_;  // DC admittance matrix is real (susceptance B)
        CplxVect acSbus_;
        RealVect dcPbus_;  // DC power injection is real (active power P)
        SolverBusIdVect bus_pv_ac_;  // id are the solver internal id and NOT the initial id
        SolverBusIdVect bus_pq_ac_;  // id are the solver internal id and NOT the initial id

        // TODO have version of the stuff above for the public api, indexed with "me" and not "solver"

        // to solve the newton raphson
        AlgorithmSelector _algo;
        AlgorithmSelector _dc_algo;

        // forced angle-reference slack bus (gridmodel id, -1 = none). Declared
        // LAST so existing member offsets are unchanged (ABI-stable for the
        // gpusim2grid cross-module LSGrid cast). See set_reference_slack_bus.
        int _forced_ref_slack_bus_id = -1;

        // ---- DC twins of the solver-side data ---------------------------------
        // The AC and the DC solver each cache their own bus labelling, matrix and
        // injection vector; these three used to be the exception -- one shared
        // copy for both families, silently overwritten by whichever solved last.
        // That is what made `unset_changes()` unsafe across families (a family
        // could reuse its own maps while reading the other's pv / pq split) and
        // it is what makes per-family cache reuse possible now that it is gone.
        // Appended here, after _forced_ref_slack_bus_id, for the same ABI reason
        // as above: the AC copies keep the offsets they always had.
        RealVect slack_weights_dc_;
        SolverBusIdVect bus_pv_dc_;  // id are the solver internal id and NOT the initial id
        SolverBusIdVect bus_pq_dc_;  // id are the solver internal id and NOT the initial id
        // bus connectivity each family's cache was last built with (the AC twin is
        // `last_bus_status_saved_`, kept under its historical name because it is
        // the one written to / read from StateRes -- see get_state / set_state,
        // which restores it into BOTH families).
        std::vector<bool> last_bus_status_dc_;

        // May a powerflow reuse what the previous one of the same family built
        // (Ybus / Sbus, the bus labelling, the pv-pq split, the slack weights)?
        // ON by default: after every powerflow that converges, that family's grid
        // is marked "nothing changed since" automatically, so the next one only
        // re-stamps what the grid actually reports as modified. Turn a family off
        // to force it to rebuild everything on every solve -- the answer is the
        // same either way, so this is a debugging / paranoia switch, not a
        // correctness one. See allow_cache_reuse() / prevent_cache_reuse().
        bool allow_ac_cache_reuse_ = true;
        bool allow_dc_cache_reuse_ = true;

        // Set when a family's powerflow diverged. The DATA it built (Ybus / Sbus,
        // the labelling, the pv-pq split) is still a correct picture of the grid
        // -- divergence is a numerical failure, not a data one -- so it is kept
        // and stays reusable. What is not reusable is the ALGORITHM's own state:
        // a half-converged iterate, a factorization of a matrix it gave up on.
        // That is reset on the spot, and this flag makes the next powerflow hand
        // the algorithm an "everything changed" control so it rebuilds its
        // internals (Jacobian sparsity, Bp/Bpp, the linear solver) from the
        // cached matrices instead of assuming they are still there. The built-in
        // algorithms would manage without it -- their reset() raises their own
        // `need_factorize_`, which gates the same rebuilds -- but an external
        // (plugin) solver is promised nothing beyond the AlgoControl it is given.
        bool ac_algo_needs_rebuild_ = false;
        bool dc_algo_needs_rebuild_ = false;
};


} // namespace ls2g

#endif  //LSGRID_H