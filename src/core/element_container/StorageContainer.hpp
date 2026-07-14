// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef STORAGE_CONTAINER_H
#define STORAGE_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"

namespace ls2g {


class StorageContainer;
class LS2G_API StorageInfo : public OneSideContainer_PQ::OneSidePQInfo
{
    public:
        inline StorageInfo(const StorageContainer & r_data_storage, int my_id) noexcept;
};


/**
This class is a container for all storage units (batteries) on the grid.

Storage units are modeled as PQ injections using the **load convention** (same as
in pandapower and grid2op): positive `target_p` means the unit is charging (power is
taken from the grid), negative `target_p` means the unit is discharging (power is
injected in the grid).

NOTE: PowSyBl / IIDM batteries use the opposite (generator) convention for their
`target_p` / `target_q` setpoints; the pypowsybl converter
(`init_from_pypowsybl`) negates them before feeding this container.

This is a dedicated container (rather than reusing :class:`LoadContainer`) so that
storage-specific behaviour can diverge in the future without affecting loads.
**/
class LS2G_API StorageContainer final: public OneSideContainer_PQ, public IteratorAdder<StorageContainer, StorageInfo>
{
    friend class StorageInfo;

    public:
        using DataInfo = StorageInfo;

    // regular implementation
    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
           OneSideContainer_PQ::StateRes  // state of the base class
           > ;

        StorageContainer() noexcept = default;
        virtual ~StorageContainer() noexcept = default;

        // pickle (python)
        StorageContainer::StateRes get_state() const;
        void set_state(StorageContainer::StateRes & my_state);

        // fast binary serialization (additive alternative to pickle, see BinaryArchive.hpp)
        void save_binary(const std::string & path, bool atomic = true) const;
        static StorageContainer load_binary(const std::string & path);
        static const char * binary_type_tag() { return "StorageContainer"; }  // written into / checked against the binary file header

        void init(const Eigen::Ref<const RealVect> & storage_p_mw,
                  const Eigen::Ref<const RealVect> & storage_q_mvar,
                  const Eigen::Ref<const Eigen::VectorXi> & storage_bus_id
                  )
        {
            init_osc_pq(storage_p_mw,
                        storage_q_mvar,
                        storage_bus_id,
                        "storages");
            reset_results();
        }

        void fillSbus(CplxVect & Sbus, const SolverBusIdVect & id_grid_to_solver, bool ac) const override;

    protected:
        void _compute_results(const Eigen::Ref<const RealVect> & /*Va*/,
                                    const Eigen::Ref<const RealVect> & /*Vm*/,
                                    const Eigen::Ref<const CplxVect> & /*V*/,
                                    const SolverBusIdVect & /*id_grid_to_solver*/,
                                    const Eigen::Ref<const RealVect> & /*bus_vn_kv*/,
                                    real_type /*sn_mva*/,
                                    bool ac) override
                                    {

                                            set_osc_pq_res_p();
                                            set_osc_pq_res_q(ac);
                                    }
};

inline StorageInfo::StorageInfo(const StorageContainer & r_data_storage, int my_id) noexcept:
        OneSidePQInfo(r_data_storage, my_id) {}


} // namespace ls2g

#endif  //STORAGE_CONTAINER_H
