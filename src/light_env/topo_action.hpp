// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TOPO_ACTION_H
#define TOPO_ACTION_H

#include "GridModel.hpp"

#include "light_env_utils.hpp"


class TopoAction
{
    protected:
        // typedef std::vector<std::pair<int, int> > ElBusMapping;
        typedef std::map<int, int> ElBusMapping;

    public:
        // TODO
        TopoAction() {}

        void add_element(ElementType el_type, int el_id, int bus_global_id){
            ElBusMapping & buses_tab = choose_bus_table(el_type);
            buses_tab.emplace(std::pair<int, int>(el_id, bus_global_id));
        }

        // TODO
        void apply_to_gridmodel(GridModel & grid) const {
            update_grid_values(grid,
                               loads_bus,
                               &GridModel::change_bus_load,
                               &GridModel::deactivate_load,
                               &GridModel::reactivate_load);
            update_grid_values(grid,
                               gens_bus,
                               &GridModel::change_bus_gen,
                               &GridModel::deactivate_gen,
                               &GridModel::reactivate_gen);
            update_grid_values(grid,
                               lines_or_bus,
                               &GridModel::change_bus1_powerline,
                               &GridModel::deactivate_powerline,
                               &GridModel::reactivate_powerline);
            update_grid_values(grid,
                               lines_ex_bus,
                               &GridModel::change_bus2_powerline,
                               &GridModel::deactivate_powerline,
                               &GridModel::reactivate_powerline);
            update_grid_values(grid,
                               trafo_hv_bus,
                               &GridModel::change_bus1_trafo,
                               &GridModel::deactivate_trafo,
                               &GridModel::reactivate_trafo);
            update_grid_values(grid,
                               trafo_lv_bus,
                               &GridModel::change_bus2_trafo,
                               &GridModel::deactivate_trafo,
                               &GridModel::reactivate_trafo);
            update_grid_values(grid,
                               storages_bus,
                               &GridModel::change_bus_storage,
                               &GridModel::deactivate_storage,
                               &GridModel::reactivate_storage);

            // TODO element not supported in grid2op:
            // - hvdc
            // - shunts
            // - static gen
        }

        // TODO check compliance (correct sub, correct bus)
        void check_validity(const GridModel & grid) const {

        }
    
    protected:
        template<class FunctorTypeChange, class FunctorTypeDisco, class FunctorTypeReco>
        void update_grid_values(GridModel & grid,
                                const ElBusMapping & new_values,
                                FunctorTypeChange fun_change,
                                FunctorTypeDisco fun_disc,
                                FunctorTypeReco fun_reco) const
        {
            for(const auto & el : new_values)
            {
                const int el_id = el.first;
                const int bus_global_id = el.second;
                if(bus_global_id == 0) continue;  // element not changed
                if(bus_global_id == -1)
                {
                    // element disconnected
                    (grid.*fun_disc)(el_id);
                }else{
                    // element changed bus (reconnected or simply changed)
                    (grid.*fun_reco)(el_id);
                    (grid.*fun_change)(el_id, GlobalBusId(bus_global_id));
                }
            }
        }

        ElBusMapping & choose_bus_table(ElementType el_type){
            switch (el_type)
            {
            case ElementType::load:
                return loads_bus;
            case ElementType::gen:
                return gens_bus;
            case ElementType::line_or:
                return lines_or_bus;
            case ElementType::line_ex:
                return lines_ex_bus;
            case ElementType::storage:
                return storages_bus;
            case ElementType::trafo_hv:
                return trafo_hv_bus;
            case ElementType::trafo_lv:
                return trafo_lv_bus;
            default:
                std::ostringstream exc_;
                exc_ << "TopoAction modification of element type ";
                exc_ << el_type << " is not supported at the moment.";
                throw std::runtime_error(exc_.str());
            }
        }

    protected:
        // bus id should be global (lightsim2grid / gridmodel) bus and
        // not local bus id
        ElBusMapping loads_bus;
        ElBusMapping gens_bus;
        ElBusMapping storages_bus;
        ElBusMapping static_gens_bus;
        ElBusMapping shunts_bus;
        ElBusMapping lines_or_bus;
        ElBusMapping lines_ex_bus;
        ElBusMapping trafo_hv_bus;
        ElBusMapping trafo_lv_bus;
        ElBusMapping dc_lines_or_bus;
        ElBusMapping dc_lines_ex_bus;
        
};

#endif // TOPO_ACTION_H
