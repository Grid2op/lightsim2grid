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
        // no members
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
        using StateRes = std::tuple<OneSideContainer_PQ::StateRes >;
        enum StateResIdx {
            OSC_PQ_STATE = 0,
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
            reset_results();
        }
    
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
        // physical properties

        // input data

        //output data

};

inline ShuntInfo::ShuntInfo(const ShuntContainer & r_data_shunt, int my_id) noexcept:
OneSidePQInfo(r_data_shunt, my_id){}


} // namespace ls2g

#endif  //SHUNT_CONTAINER_H
