// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SHUNT_CONTAINER_H
#define SHUNT_CONTAINER_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"

namespace ls2g {

class ShuntContainer;
class LS2G_API ShuntInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        // its sections (has_sections false: none, the rest is meaningless)
        bool has_sections;
        int section_count;
        int max_section_count;
        // what they regulate, data only (see ShuntContainer::set_section_regulation)
        bool regulating;
        real_type target_vm_pu;
        real_type target_deadband_pu;
        int regulated_bus;

        inline ShuntInfo(const ShuntContainer & r_data_shunt, int my_id) noexcept;
};

/**
This class is a container for all shunts on the grid.

The convention used for the shunt is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html#electric-model
**/
class LS2G_API ShuntContainer final: public OneSideContainer_PQ, public IteratorAdder<ShuntContainer, ShuntInfo>
{
    friend class ShuntInfo;
    public:
        using DataInfo = ShuntInfo;

    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
                   OneSideContainer_PQ::StateRes,
                   std::vector<std::vector<real_type> >,  // section_p_mw_
                   std::vector<std::vector<real_type> >,  // section_q_mvar_
                   std::vector<int>,  // section_count_
                   std::vector<bool>,  // regulating_
                   std::vector<real_type>,  // target_vm_pu_
                   std::vector<real_type>,  // target_deadband_pu_
                   std::vector<int>   // regulated_bus_
               >;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
            SECTION_P,
            SECTION_Q,
            SECTION_COUNT,
            REGULATING,
            TARGET_VM,
            TARGET_DEADBAND,
            REGULATED_BUS,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "ShuntContainer::StateRes and StateResIdx do not match");
        
        ShuntContainer() noexcept = default;
        ~ShuntContainer() noexcept override = default;
        
        
        void init(const Eigen::Ref<const RealVect> & shunt_p_mw,
                  const Eigen::Ref<const RealVect> & shunt_q_mvar,
                  const Eigen::Ref<const Eigen::VectorXi> & shunt_bus_id
                  )
        {
            init_osc_pq(shunt_p_mw,
                        shunt_q_mvar,
                        shunt_bus_id,
                        "shunts");
            const std::size_t n = static_cast<std::size_t>(nb());
            section_p_mw_.assign(n, std::vector<real_type>());
            section_q_mvar_.assign(n, std::vector<real_type>());
            section_count_.assign(n, 0);
            regulating_.assign(n, false);
            target_vm_pu_.assign(n, 0.);
            target_deadband_pu_.assign(n, 0.);
            regulated_bus_.assign(n, -1);
            reset_results();
        }

        /**
         * The sections of shunt `el`: for k = 1 .. max_section_count, the active and reactive
         * power (MW, MVar at 1 pu, the same convention as init) of the shunt with k sections
         * on -- cumulative, IIDM's non-linear model (a linear one is k times its
         * per-section value) -- and the count on now, from 0 (nothing) to max_section_count.
         * The shunt's p and q become those of that count.
         */
        void set_sections(int el, int section_count,
                          const std::vector<real_type> & p_mw,
                          const std::vector<real_type> & q_mvar,
                          DualAlgoControl & solver_control);
        /// the voltage the sections regulate (data only): `target_vm_pu` and the deadband in
        /// pu of `regulated_bus`' nominal voltage, `regulated_bus` a grid bus id
        void set_section_regulation(int el, bool regulating, real_type target_vm_pu,
                                    real_type target_deadband_pu, int regulated_bus);
        /// switch `section_count` sections of shunt `el` on: its p and q follow
        void change_section_count(int el, int section_count, DualAlgoControl & solver_control);
        bool has_sections(int el) const { return !section_q_mvar_[static_cast<std::size_t>(el)].empty(); }
        int get_section_count(int el) const { return section_count_[static_cast<std::size_t>(el)]; }
        int get_max_section_count(int el) const { return static_cast<int>(section_q_mvar_[static_cast<std::size_t>(el)].size()); }
    
        // pickle (python)
        ShuntContainer::StateRes get_state() const;
        void set_state(ShuntContainer::StateRes & my_state );

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static ShuntContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "ShuntContainer"; }  // written into / checked against the binary file header
        
        void _fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                       bool ac,
                       const SolverBusIdVect & id_grid_to_solver,
                       real_type sn_mva) const override;
        void _fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                         std::vector<Eigen::Triplet<real_type> > & Bpp,
                         const SolverBusIdVect & id_grid_to_solver,
                         real_type sn_mva,
                         FDPFMethod xb_or_bx) const override;
    protected:
        void _fillSbus(Eigen::Ref<CplxVect> Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;  // in DC i need that
        
    protected:
        // a shunt is in Ybus (AC) AND in Sbus (DC only, its active part: _fillSbus
        // stamps nothing in AC), so the Sbus flag is raised on the DC family alone
        void _on_change_p(int shunt_id, real_type new_p, DualAlgoControl & solver_control) override
        {
            if(abs(target_p_mw_(shunt_id) - new_p) > _tol_equal_float){
                solver_control.tell_recompute_ybus();
                solver_control.dc_algo_controler().tell_recompute_sbus();
            }
        }
        void _on_change_q(int shunt_id, real_type new_q, DualAlgoControl & solver_control) override
        {
            if(abs(target_q_mvar_(shunt_id) - new_q) > _tol_equal_float){
                solver_control.tell_recompute_ybus();
            }
        }
        void _on_change_bus(int /*el_id*/, GridModelBusId /*new_bus_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_ybus();
            solver_control.tell_one_el_changed_bus();
            solver_control.dc_algo_controler().tell_recompute_sbus();
        }
        void _on_deactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_ybus();
            solver_control.tell_one_el_changed_bus();
            solver_control.dc_algo_controler().tell_recompute_sbus();
        }
        void _on_reactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_ybus();
            solver_control.tell_one_el_changed_bus();
            solver_control.dc_algo_controler().tell_recompute_sbus();
        }

        void _compute_res_pq(
            const Eigen::Ref<const RealVect> & Va,
            const Eigen::Ref<const RealVect> & Vm,
            const Eigen::Ref<const CplxVect> & V,
            const SolverBusIdVect & id_grid_to_solver,
            const Eigen::Ref<const RealVect> & bus_vn_kv,
            real_type sn_mva,
            bool ac) override;

    protected:
        // the p and q (MW, MVar) of shunt `el` with `count` sections on
        real_type _section_p(int el, int count) const {
            return count == 0 ? 0. : section_p_mw_[static_cast<std::size_t>(el)][static_cast<std::size_t>(count - 1)];
        }
        real_type _section_q(int el, int count) const {
            return count == 0 ? 0. : section_q_mvar_[static_cast<std::size_t>(el)][static_cast<std::size_t>(count - 1)];
        }
        void _check_section_count(int el, int section_count, const char * where) const;

        // the sections, see set_sections (empty: none)
        std::vector<std::vector<real_type> > section_p_mw_;
        std::vector<std::vector<real_type> > section_q_mvar_;
        std::vector<int> section_count_;
        // what they regulate, see set_section_regulation
        std::vector<bool> regulating_;
        std::vector<real_type> target_vm_pu_;
        std::vector<real_type> target_deadband_pu_;
        std::vector<int> regulated_bus_;
};

inline ShuntInfo::ShuntInfo(const ShuntContainer & r_data_shunt, int my_id) noexcept:
OneSidePQInfo(r_data_shunt, my_id),
has_sections(false),
section_count(0),
max_section_count(0),
regulating(false),
target_vm_pu(0.),
target_deadband_pu(0.),
regulated_bus(-1)
{
    if(my_id < 0) return;
    if(my_id >= static_cast<int>(r_data_shunt.section_count_.size())) return;
    const std::size_t k = static_cast<std::size_t>(my_id);
    has_sections = r_data_shunt.has_sections(my_id);
    section_count = r_data_shunt.section_count_[k];
    max_section_count = r_data_shunt.get_max_section_count(my_id);
    regulating = r_data_shunt.regulating_[k];
    target_vm_pu = r_data_shunt.target_vm_pu_[k];
    target_deadband_pu = r_data_shunt.target_deadband_pu_[k];
    regulated_bus = r_data_shunt.regulated_bus_[k];
}


} // namespace ls2g

#endif  //SHUNT_CONTAINER_H
