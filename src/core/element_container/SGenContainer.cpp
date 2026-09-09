// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SGenContainer.hpp"
#include "BinaryArchive.hpp"

#include <iostream>

namespace ls2g {

void SGenContainer::init(const Eigen::Ref<const RealVect> & sgen_p,
                         const Eigen::Ref<const RealVect> & sgen_q,
                         const Eigen::Ref<const RealVect> & sgen_pmin,
                         const Eigen::Ref<const RealVect> & sgen_pmax,
                         const Eigen::Ref<const RealVect> & sgen_qmin,
                         const Eigen::Ref<const RealVect> & sgen_qmax,
                         const Eigen::Ref<const Eigen::VectorXi> & sgen_bus_id)
{
    init_osc_pq(sgen_p, sgen_q, sgen_bus_id, "static_generators");
    
    int size = nb();
    check_size(sgen_pmin, size, "sgen_pmin");
    check_size(sgen_pmax, size, "sgen_pmax");
    check_size(sgen_qmin, size, "sgen_qmin");
    check_size(sgen_qmax, size, "sgen_qmax");

    p_min_mw_ = sgen_pmin;
    p_max_mw_ = sgen_pmax;
    q_min_mvar_ = sgen_qmin;
    q_max_mvar_ = sgen_qmax;
    reset_results();
}

SGenContainer::StateRes SGenContainer::get_state() const
{
     std::vector<real_type> p_min(p_min_mw_.begin(), p_min_mw_.end());
     std::vector<real_type> p_max(p_max_mw_.begin(), p_max_mw_.end());
     std::vector<real_type> q_min(q_min_mvar_.begin(), q_min_mvar_.end());
     std::vector<real_type> q_max(q_max_mvar_.begin(), q_max_mvar_.end());
     SGenContainer::StateRes res(get_osc_pq_state(), p_min, p_max, q_min, q_max);
     return res;
}

void SGenContainer::set_state(SGenContainer::StateRes & my_state )
{    
    set_osc_pq_state(std::get<0>(my_state));

    std::vector<real_type> & p_min = std::get<1>(my_state);
    std::vector<real_type> & p_max = std::get<2>(my_state);
    std::vector<real_type> & q_min = std::get<3>(my_state);
    std::vector<real_type> & q_max = std::get<4>(my_state);
    const auto size = nb();

    GenericContainer::check_size(p_min, size, "p_min");
    GenericContainer::check_size(p_max, size, "p_max");
    GenericContainer::check_size(q_min, size, "q_min");
    GenericContainer::check_size(q_max, size, "q_max");

    p_min_mw_ = RealVect::Map(p_min.data(), size);
    p_max_mw_ = RealVect::Map(p_max.data(), size);
    q_min_mvar_ = RealVect::Map(q_min.data(), size);
    q_max_mvar_ = RealVect::Map(q_max.data(), size);
    reset_results();
}

void SGenContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

SGenContainer SGenContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<SGenContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
