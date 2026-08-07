// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LSGrid.hpp"
#include "AlgorithmSelector.hpp"  // to avoid circular references
#include "BinaryArchive.hpp"

#include <cmath>      // std::isfinite (check_positive_finite)
#include <queue>
#include <sstream>
#include <stdexcept>

namespace ls2g {


LSGrid::LSGrid(const LSGrid & other)
{
    init_vm_pu_ = other.init_vm_pu_;
    sn_mva_ = other.sn_mva_;
    compute_results_ = other.compute_results_;
    init_kwargs_ = other.init_kwargs_;
    _bus_fusion_rep = other._bus_fusion_rep;

    // copy the powersystem representation
    // 1. bus
    last_bus_status_saved_ = other.last_bus_status_saved_;
    substations_ = other.substations_;
    max_nb_bus_per_sub_ = substations_.nmax_busbar_per_sub();
    n_sub_ = substations_.nb_sub();

    set_ls_to_orig_internal(other._ls_to_orig);  // sets also orig_to_ls

    // 2. powerline
    powerlines_ = other.powerlines_;

    // 3. shunt
    shunts_ = other.shunts_;

    // 4. transformers
    // have the r, x, h and ratio
    // ratio is computed from the tap, so maybe store tap num and tap_step_pct
    trafos_ = other.trafos_;

    // 5. generators
    total_q_min_per_bus_ = RealVect();
    total_q_max_per_bus_ = RealVect();
    generators_ = other.generators_;

    // 6. loads
    loads_ = other.loads_;

    // 7. static generators
    sgens_ = other.sgens_;

    // 8. storage units
    storages_ = other.storages_;

    // hvdc lines
    hvdc_lines_ = other.hvdc_lines_;

    // static var compensators
    svcs_ = other.svcs_;

    // assign the right solver. By *name*, not by AlgorithmType: an external
    // (plugin) solver has type AlgorithmType::Custom, which the type-based
    // overload rejects -- copying a grid using a plugin used to throw.
    reset(true, true, true);
    _algo.change_algorithm(other._algo.get_name());
    _algo.set_config(other.get_algo().get_config());
    _dc_algo.change_algorithm(other._dc_algo.get_name());
    _dc_algo.set_config(other.get_dc_algo().get_config());
}

//pickle
LSGrid::StateRes LSGrid::get_state() const 
{
    std::vector<int> ls_to_orig(_ls_to_orig.begin(), _ls_to_orig.end());
    std::string version_major = VERSION_MAJOR;
    std::string version_medium = VERSION_MEDIUM;
    std::string version_minor = VERSION_MINOR;
    auto res_substation = substations_.get_state();
    auto res_line = powerlines_.get_state();
    auto res_shunt = shunts_.get_state();
    auto res_trafo = trafos_.get_state();
    auto res_gen = generators_.get_state();
    auto res_load = loads_.get_state();
    auto res_sgen = sgens_.get_state();
    auto res_storage = storages_.get_state();
    auto res_hvdc_line = hvdc_lines_.get_state();
    auto res_svc = svcs_.get_state();

    AlgoConfig ac_algo_cfg = get_ac_algo_config();
    AlgoConfig dc_algo_cfg = get_dc_algo_config();
    auto res_ac_algo_cfg = std::make_tuple(ac_algo_cfg.int_params, ac_algo_cfg.real_params);
    auto res_dc_algo_cfg = std::make_tuple(dc_algo_cfg.int_params, dc_algo_cfg.real_params);

    std::vector<std::string> init_kwargs_keys;
    std::vector<std::string> init_kwargs_values;
    init_kwargs_keys.reserve(init_kwargs_.size());
    init_kwargs_values.reserve(init_kwargs_.size());
    for (const auto & kv : init_kwargs_) {
        init_kwargs_keys.push_back(kv.first);
        init_kwargs_values.push_back(kv.second);
    }
    std::vector<int> bus_fusion_rep(_bus_fusion_rep.begin(), _bus_fusion_rep.end());

    LSGrid::StateRes res(version_major,
                            version_medium,
                            version_minor,
                            ls_to_orig,
                            init_vm_pu_,
                            sn_mva_,
                            last_bus_status_saved_,
                            res_substation,
                            res_line,
                            res_shunt,
                            res_trafo,
                            res_gen,
                            res_load,
                            res_sgen,
                            res_storage,
                            res_hvdc_line,
                            res_svc,
                            _algo.get_name(),
                            _dc_algo.get_name(),
                            res_ac_algo_cfg,
                            res_dc_algo_cfg,
                            init_kwargs_keys,
                            init_kwargs_values,
                            bus_fusion_rep
                            );
    return res;
};

void LSGrid::set_state(LSGrid::StateRes & my_state, bool restore_algorithm)
{
    // after loading back, the instance need to be reset anyway
    // TODO see if it's worth the trouble NOT to do it
    algo_controler_.ac_algo_controler().tell_all_changed();
    algo_controler_.dc_algo_controler().tell_all_changed();
    compute_results_ = true;

    // extract data from the state
    std::string version_major = std::get<VERSION_MAJOR_ID>(my_state);
    std::string version_medium = std::get<VERSION_MEDIUM_ID>(my_state);
    std::string version_minor = std::get<VERSION_MINOR_ID>(my_state);
    if((version_major != VERSION_MAJOR )| (version_medium != VERSION_MEDIUM) | (version_minor != VERSION_MINOR))
    {
        std::ostringstream exc_;
        exc_ << "LSGrid::set_state: Wrong version. You tried to load a lightsim2grid model saved with version ";
        exc_ << version_major << "." << version_medium << "." << version_minor;
        exc_ << " while currently using the package on version ";
        exc_ << VERSION_MAJOR << "." << VERSION_MEDIUM << "." << VERSION_MINOR;
        exc_ << "It is not possible. Please reinstall it.";
        throw std::runtime_error(exc_.str());
    }
    const std::vector<int> & ls_to_pp = std::get<LS_TO_ORIG_ID>(my_state);
    init_vm_pu_ = std::get<INIT_VM_PU_ID>(my_state);
    sn_mva_ = std::get<SN_MVA_ID>(my_state);
    const std::vector<bool> & last_bus_status_saved = std::get<BUS_STATUS_ID>(my_state);
    SubstationContainer::StateRes & state_substations = std::get<SUBSTATION_ID>(my_state);
    // powerlines
    LineContainer::StateRes & state_lines = std::get<LINE_ID>(my_state);
    // shunts
    ShuntContainer::StateRes & state_shunts = std::get<SHUNT_ID>(my_state);
    // trafos
    TrafoContainer::StateRes & state_trafos = std::get<TRAFO_ID>(my_state);
    // generators
    // total_q_min_per_bus_;
    // total_q_max_per_bus_;
    // total_gen_per_bus_;
    GeneratorContainer::StateRes & state_gens = std::get<GEN_ID>(my_state);
    // loads
    LoadContainer::StateRes & state_loads = std::get<LOAD_ID>(my_state);
    // static gen
    SGenContainer::StateRes & state_sgens= std::get<SGEN_ID>(my_state);
    // storage units
    StorageContainer::StateRes & state_storages = std::get<STORAGE_ID>(my_state);
    // hvdc lines
    HvdcLineContainer::StateRes & state_hvdc_lines = std::get<HVDC_ID>(my_state);
    // static var compensators
    SvcContainer::StateRes & state_svcs = std::get<SVC_ID>(my_state);
    // algo configs (scaling/refactor/line-search params)
    auto & state_ac_algo_cfg = std::get<AC_ALGO_CONFIG_ID>(my_state);
    auto & state_dc_algo_cfg = std::get<DC_ALGO_CONFIG_ID>(my_state);
    // relevant kwargs the grid was built with (eg by init_from_pypowsybl)
    const std::vector<std::string> & init_kwargs_keys = std::get<INIT_KWARGS_KEYS_ID>(my_state);
    const std::vector<std::string> & init_kwargs_values = std::get<INIT_KWARGS_VALUES_ID>(my_state);
    // fused-bus representative lookup (see get_bus_fusion_rep())
    const std::vector<int> & bus_fusion_rep = std::get<BUS_FUSION_REP_ID>(my_state);

    // substations
    last_bus_status_saved_ = last_bus_status_saved;
    substations_.set_state(state_substations);
    max_nb_bus_per_sub_ = substations_.nmax_busbar_per_sub();
    n_sub_ = substations_.nb_sub();

    // assign it to this instance -- must run *after* substations_ is restored above,
    // since set_ls_to_orig() validates the vector size against substations_.nb_bus()
    // (on a freshly default-constructed LSGrid, as pickle/binary restore does, that
    // would otherwise still be 0 and any non-empty ls_to_orig would wrongly be
    // rejected as a size mismatch).
    set_ls_to_orig(IntVect::Map(ls_to_pp.data(), ls_to_pp.size()));  // set also _orig_to_ls

    // elements
    // 1. powerlines
    powerlines_.set_state(state_lines);
    // 2. shunts
    shunts_.set_state(state_shunts);
    // 3. trafos
    trafos_.set_state(state_trafos);
    // 4. gen
    total_q_min_per_bus_ = RealVect();
    total_q_max_per_bus_ = RealVect();
    generators_.set_state(state_gens);
    // 5. loads
    loads_.set_state(state_loads);
    // 6. static generators
    sgens_.set_state(state_sgens);
    // 7. storage units
    storages_.set_state(state_storages);
    // hvdc lines
    hvdc_lines_.set_state(state_hvdc_lines);
    svcs_.set_state(state_svcs);

    // handle the solver
    reset(true, true, true);
    if (restore_algorithm) {
        _restore_algorithm(_algo, std::get<AC_ALGO_NAME_ID>(my_state), "AC");
        _restore_algorithm(_dc_algo, std::get<DC_ALGO_NAME_ID>(my_state), "DC");

        // algo configs -- must be restored *after* change_algorithm() above,
        // since set_config() operates on the currently-selected concrete solver
        // (same order as the copy constructor). They describe the tuning of the
        // solver we just re-selected, so they are skipped together with it.
        AlgoConfig ac_algo_cfg;
        ac_algo_cfg.int_params = std::get<0>(state_ac_algo_cfg);
        ac_algo_cfg.real_params = std::get<1>(state_ac_algo_cfg);
        set_ac_algo_config(ac_algo_cfg);

        AlgoConfig dc_algo_cfg;
        dc_algo_cfg.int_params = std::get<0>(state_dc_algo_cfg);
        dc_algo_cfg.real_params = std::get<1>(state_dc_algo_cfg);
        set_dc_algo_config(dc_algo_cfg);
    }

    // relevant kwargs the grid was built with (eg by init_from_pypowsybl).
    // The map is serialized as two parallel vectors, whose lengths are stored
    // independently (both in a pickle and in the binary format): a file declaring
    // more keys than values makes the loop below read past the end of
    // init_kwargs_values -- an out-of-bounds read over std::string objects, ie
    // constructing a std::string from whatever the heap holds there. Check first.
    if (init_kwargs_keys.size() != init_kwargs_values.size()) {
        std::ostringstream exc_;
        exc_ << "LSGrid::set_state: the serialized `_init_kwargs` is inconsistent: "
             << init_kwargs_keys.size() << " keys but " << init_kwargs_values.size()
             << " values. They are two parallel vectors and must have the same length.";
        throw std::runtime_error(exc_.str());
    }
    init_kwargs_.clear();
    for (std::size_t i = 0; i < init_kwargs_keys.size(); ++i) {
        init_kwargs_[init_kwargs_keys[i]] = init_kwargs_values[i];
    }

    // fused-bus representative lookup -- must run after substations_ is restored
    // above, same reasoning as set_ls_to_orig() (validates against total_bus()).
    set_bus_fusion_rep(IntVect::Map(bus_fusion_rep.data(), bus_fusion_rep.size()));

    // Now that every container has been restored, validate the whole grid: a
    // pickle or binary file is only length-checked while being read, so an
    // out-of-range bus / substation / topo-vector index (or a NaN electrical
    // value) would otherwise slip through and cause an out-of-bounds access on
    // the next powerflow. check_grid() turns that into a clean exception here.
    check_grid();
};

void LSGrid::_restore_algorithm(AlgorithmSelector & algo_selector,
                                const std::string & name,
                                const char * ac_or_dc)
{
    // A name that cannot even be a solver name never was one: the state is
    // corrupted (or was not produced by lightsim2grid). Say so, rather than
    // sending the reader looking for a plugin that never existed.
    if (!is_valid_solver_name(name)) {
        std::ostringstream exc_;
        exc_ << "LSGrid::set_state: the " << ac_or_dc << " solver name read from this state, '"
             << printable(name) << "', is not a valid solver name (" << solver_name_rule()
             << "). This state is corrupted.";
        throw std::runtime_error(exc_.str());
    }
    if (AlgorithmRegistry::instance().is_registered(name)) {
        algo_selector.change_algorithm(name);
        return;
    }
    // The name was resolvable when the grid was saved but is not now. Either the
    // solver comes from a plugin that has not been loaded in this process, or it
    // needs an optional linear-algebra backend (KLU / NICSLU / CKTSO) this build
    // was not compiled with. Say so, and say what to do about it -- the grid data
    // itself is perfectly fine, only the solver choice cannot be honoured.
    // `name` comes straight from the file: escape it (a corrupted one can hold
    // arbitrary bytes, which would make pybind11 raise UnicodeDecodeError while
    // converting what() instead of the RuntimeError we mean to report).
    std::ostringstream exc_;
    exc_ << "LSGrid::set_state: this grid was saved using the " << ac_or_dc
         << " solver '" << printable(name) << "', which is not available here. ";
    if (name.rfind("NR_", 0) == 0 || name.rfind("NRSing_", 0) == 0 ||
        name.rfind("DC_", 0) == 0 || name.rfind("FDPF_", 0) == 0) {
        exc_ << "It looks like a built-in solver relying on an optional linear-algebra "
             << "backend (KLU / NICSLU / CKTSO) that this build of lightsim2grid does not "
             << "include: reinstall lightsim2grid with that support, or re-save the grid "
             << "after selecting a solver available everywhere (eg 'NR_SparseLU'). ";
    } else {
        exc_ << "It is most likely provided by a solver plugin: load it first with "
             << "lightsim2grid.load_algorithm_plugin(<path to the plugin>), then load this "
             << "grid again. ";
    }
    exc_ << "Solvers currently available: ";
    const std::vector<std::string> available = AlgorithmRegistry::instance().available_algorithm_names();
    for (std::size_t i = 0; i < available.size(); ++i) {
        if (i) exc_ << ", ";
        exc_ << "'" << available[i] << "'";
    }
    exc_ << ".";
    throw std::runtime_error(exc_.str());
}

void LSGrid::check_grid() const
{
    // The two grid-wide scalars, before anything else. `set_state` assigns them
    // straight from the serialized state (bypassing the setters), and the loaders take
    // them verbatim from their source file -- `init_from_powermodels` /
    // `init_from_matpower` do `set_sn_mva(float(network["baseMVA"]))`, and
    // `init_from_pandapower` `set_sn_mva(pp_net.sn_mva)`.
    //
    // sn_mva_ is the base power of the whole per-unit system: Sbus is divided by it,
    // every MW / MVar result multiplied back by it, and ac_pf even scales the solver
    // tolerance with it. A degenerate value does not make the powerflow fail, it makes
    // it *quietly wrong*: with sn_mva_ = NaN the DC powerflow reports CONVERGENCE and
    // hands back NaN flows (the built-in solvers' own finiteness guards see a perfectly
    // finite Va -- Bbus does not involve sn_mva -- and the size/finiteness check on the
    // solver output is deliberately reserved for external solvers), and with a negative
    // one it reports a plausible-looking but sign-inverted per-unit system.
    check_positive_finite(sn_mva_, "sn_mva");
    check_positive_finite(init_vm_pu_, "init_vm_pu");

    // The substation container FIRST: it defines nb_bus / nb_sub, the bounds every
    // per-element check below is expressed against, and it carries the vector
    // (bus_status_) those very ids are used to index. Validating elements against a
    // self-inconsistent substation container would prove nothing.
    substations_.check_valid();

    const int nb_bus = static_cast<int>(substations_.nb_bus());
    const int nb_sub = substations_.nb_sub();

    // Per-container range + finiteness checks. Each container appends the
    // pos_topo_vect entries it carries (an optional field) to the accumulator it is
    // given, for the global permutation check below.
    //
    // Only the containers update_topo() actually walks may contribute to
    // `all_pos_topo_vect`: that vector's own length is the `dim_topo` the permutation
    // is proved against, and update_topo() sizes its caller arrays as
    // `nb loads + nb gens + nb storages + 2 * nb lines + 2 * nb trafos`. Letting a
    // shunt / sgen / svc / hvdc position into the same pot inflates that bound, so a
    // *validated* load position could still be >= the length of the arrays
    // update_topo() indexes with it. Those containers have no position in the grid2op
    // topology vector to begin with (there is no setter for one), so a state carrying
    // one is inconsistent: collect them apart and reject.
    std::vector<int> all_pos_topo_vect;   // elements update_topo() drives
    std::vector<int> pos_topo_vect_not_in_topo;  // must stay empty
    powerlines_.check_valid(nb_bus, nb_sub, substations_, all_pos_topo_vect);
    trafos_.check_valid(nb_bus, nb_sub, substations_, all_pos_topo_vect);
    generators_.check_valid(nb_bus, nb_sub, substations_, all_pos_topo_vect);
    loads_.check_valid(nb_bus, nb_sub, substations_, all_pos_topo_vect);
    storages_.check_valid(nb_bus, nb_sub, substations_, all_pos_topo_vect);
    shunts_.check_valid(nb_bus, nb_sub, substations_, pos_topo_vect_not_in_topo);
    sgens_.check_valid(nb_bus, nb_sub, substations_, pos_topo_vect_not_in_topo);
    hvdc_lines_.check_valid(nb_bus, nb_sub, substations_, pos_topo_vect_not_in_topo);
    svcs_.check_valid(nb_bus, nb_sub, substations_, pos_topo_vect_not_in_topo);
    if(!pos_topo_vect_not_in_topo.empty())
    {
        throw std::runtime_error(
            "LSGrid::check_grid: a shunt, static generator, svc or hvdc line declares a position "
            "in the grid2op topology vector. Only loads, generators, storage units, powerlines and "
            "transformers are part of that vector (it is what sizes the arrays passed to "
            "update_topo), so this state is inconsistent.");
    }

    // pos_topo_vect is grid2op-specific and optional: it is either set on every
    // topology-participating element or on none. When set, the collected values
    // must be a permutation of [0, dim_topo) -- that is exactly what makes it
    // safe to index the (dim_topo-sized) arrays passed to update_topo(). K here
    // is dim_topo; K distinct values all in [0, K) is precisely a permutation.
    if(!all_pos_topo_vect.empty())
    {
        const int dim_topo = static_cast<int>(all_pos_topo_vect.size());
        std::vector<char> seen(dim_topo, 0);
        for(int pos : all_pos_topo_vect)
        {
            if((pos < 0) || (pos >= dim_topo))
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: a position in the topology vector (" << pos
                     << ") is out of range [0, " << dim_topo << "). The pos_topo_vect of all "
                     << "elements must form a permutation of [0, dim_topo).";
                throw std::out_of_range(exc_.str());
            }
            if(seen[pos])
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: the position " << pos << " in the topology vector is "
                     << "assigned to more than one element (pos_topo_vect values must be unique).";
                throw std::runtime_error(exc_.str());
            }
            seen[pos] = 1;
        }
    }

    // Bus-id mapping vectors carried alongside the grid. They are never read by the
    // C++ powerflow itself (only by downstream python result views), but they are
    // part of the serialized state, they are sized against the grid, and
    // set_ls_to_orig_internal() indexes _orig_to_ls with the *values* of
    // _ls_to_orig -- so an inconsistent pair is exactly the kind of thing this
    // function exists to reject rather than discover later.
    if(_ls_to_orig.size() != 0)
    {
        if(static_cast<size_t>(_ls_to_orig.size()) != substations_.nb_bus())
        {
            std::ostringstream exc_;
            exc_ << "LSGrid::check_grid: _ls_to_orig has " << _ls_to_orig.size()
                 << " entries while the grid has " << substations_.nb_bus()
                 << " buses (it is indexed by lightsim bus id).";
            throw std::runtime_error(exc_.str());
        }
        check_ls_to_orig_values(_ls_to_orig);
    }
    if(_orig_to_ls.size() != 0)
    {
        const int nb_bus_ls = static_cast<int>(substations_.nb_bus());
        for(Eigen::Index i = 0; i < _orig_to_ls.size(); ++i)
        {
            const int ls_id = _orig_to_ls(i);
            if(ls_id == GenericContainer::_deactivated_bus_id) continue;
            if((ls_id < 0) || (ls_id >= nb_bus_ls))
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: _orig_to_ls[" << i << "] = " << ls_id
                     << " is out of range: it must be -1 or a lightsim bus id in [0, "
                     << nb_bus_ls << ").";
                throw std::out_of_range(exc_.str());
            }
        }
    }
    if(_bus_fusion_rep.size() != 0)
    {
        const int nb_bus_ls = static_cast<int>(substations_.nb_bus());
        if(static_cast<size_t>(_bus_fusion_rep.size()) != substations_.nb_bus())
        {
            std::ostringstream exc_;
            exc_ << "LSGrid::check_grid: _bus_fusion_rep has " << _bus_fusion_rep.size()
                 << " entries while the grid has " << substations_.nb_bus() << " buses.";
            throw std::runtime_error(exc_.str());
        }
        for(Eigen::Index i = 0; i < _bus_fusion_rep.size(); ++i)
        {
            const int rep = _bus_fusion_rep(i);
            if((rep < 0) || (rep >= nb_bus_ls))
            {
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: _bus_fusion_rep[" << i << "] = " << rep
                     << " is out of range [0, " << nb_bus_ls << "): every bus must be merged into "
                     << "an existing bus (identity for a bus involved in no fusion).";
                throw std::out_of_range(exc_.str());
            }
        }
    }
}

void LSGrid::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

LSGrid LSGrid::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<LSGrid>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

void LSGrid::fixup_binary_state(LSGrid::StateRes & state) {
    // must be defined in this translation unit: the VERSION_* macros here are
    // exactly the ones get_state()/set_state() above are compiled with, so
    // the rewritten fields always pass set_state()'s equality check
    std::get<VERSION_MAJOR_ID>(state) = VERSION_MAJOR;
    std::get<VERSION_MEDIUM_ID>(state) = VERSION_MEDIUM;
    std::get<VERSION_MINOR_ID>(state) = VERSION_MINOR;
}

void LSGrid::set_ls_to_orig(const Eigen::Ref<const IntVect> & ls_to_orig){
    if(ls_to_orig.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }

    if(static_cast<size_t>(ls_to_orig.size()) != substations_.nb_bus())
        throw std::runtime_error("Impossible to set the converter ls_to_orig: the provided vector has not the same size as the number of bus on the grid.");
    check_ls_to_orig_values(ls_to_orig);
    set_ls_to_orig_internal(ls_to_orig);
}

// Every entry of _ls_to_orig is either -1 ("this lightsim bus has no counterpart in
// the original grid") or an original-grid bus id, which set_ls_to_orig_internal uses
// *as an index* into the _orig_to_ls it allocates. That allocation is sized from
// `lpNorm<Infinity>()`, ie the max ABSOLUTE value, so a negative entry other than -1
// (say -5) sizes the vector from |−5| and then writes at index −5: an out-of-bounds
// heap write. A huge positive entry is the other end of the same problem (`size + 1`
// overflows int for INT_MAX, and even short of that it asks for a multi-gigabyte
// allocation for a handful of buses). Both are rejected here, before any allocation.
// This runs on the python property setter AND on set_state(), ie on every pickle and
// every binary file.
void LSGrid::check_positive_finite(real_type value, const char * name)
{
    // NB `!(value > 0.)`, not `value <= 0.`: every comparison with NaN is false, so the
    // naive form silently accepts a NaN -- which is exactly the value that does the most
    // damage here (see check_grid()).
    if(std::isfinite(value) && (value > 0.)) return;
    std::ostringstream exc_;
    exc_ << "LSGrid: '" << name << "' is " << value
         << "; it must be a finite, strictly positive number.";
    throw std::runtime_error(exc_.str());
}

void LSGrid::check_ls_to_orig_values(const Eigen::Ref<const IntVect> & ls_to_orig)
{
    // an original grid can never have more buses than the biggest vector we could
    // possibly want to allocate for it; keep the bound generous but finite.
    constexpr int max_orig_bus_id = 100000000;  // 100M buses, ~400 MB for _orig_to_ls
    for(Eigen::Index i = 0; i < ls_to_orig.size(); ++i){
        const int el = ls_to_orig(i);
        if(el == GenericContainer::_deactivated_bus_id) continue;  // -1: no counterpart, legal
        if(el < 0){
            std::ostringstream exc_;
            exc_ << "LSGrid::set_ls_to_orig: ls_to_orig[" << i << "] = " << el
                 << " is negative. Entries must be either -1 (this bus has no counterpart in the "
                 << "original grid) or a valid original-grid bus id >= 0 (it is used as an index "
                 << "into the reverse mapping).";
            throw std::out_of_range(exc_.str());
        }
        if(el > max_orig_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::set_ls_to_orig: ls_to_orig[" << i << "] = " << el
                 << " exceeds the maximum supported original-grid bus id (" << max_orig_bus_id
                 << "). The reverse mapping is sized from the largest id, so this would ask for "
                 << "an unreasonable allocation.";
            throw std::out_of_range(exc_.str());
        }
    }
}

void LSGrid::set_orig_to_ls(const Eigen::Ref<const IntVect> & orig_to_ls){
    if(orig_to_ls.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }
    size_t nb_bus_ls = 0;
    for(const auto el : orig_to_ls){
        if (el != -1) nb_bus_ls += 1;
    }
    if(nb_bus_ls != substations_.nb_bus())
        throw std::runtime_error("Impossible to set the converter orig_to_ls: the number of 'non -1' component in the provided vector does not match the number of buses on the grid.");
    // orig_to_ls[orig_bus_id] is a LIGHTSIM bus id, used below as an index into
    // _ls_to_orig (which has one slot per lightsim bus): validate the range before
    // writing anything, an unchecked entry here is an out-of-bounds heap write.
    for(Eigen::Index i = 0; i < orig_to_ls.size(); ++i){
        const int ls_id = orig_to_ls(i);
        if(ls_id == GenericContainer::_deactivated_bus_id) continue;  // -1: this original bus has no lightsim bus
        if((ls_id < 0) || (static_cast<size_t>(ls_id) >= nb_bus_ls)){
            std::ostringstream exc_;
            exc_ << "LSGrid::set_orig_to_ls: orig_to_ls[" << i << "] = " << ls_id
                 << " is out of range: entries must be either -1 (this original bus has no "
                 << "counterpart in lightsim2grid) or a lightsim bus id in [0, " << nb_bus_ls << ").";
            throw std::out_of_range(exc_.str());
        }
    }
    // _ls_to_orig is the INVERSE of _orig_to_ls: _ls_to_orig[ls_id] == orig_id iff
    // _orig_to_ls[orig_id] == ls_id. Walk the whole input (not just its first
    // nb_bus_ls entries -- the non -1 ones are not necessarily at the front) and
    // index the result by the lightsim id, so the two arrays really are inverses.
    _orig_to_ls = orig_to_ls;
    _ls_to_orig = IntVect::Constant(nb_bus_ls, -1);
    for(Eigen::Index orig_id = 0; orig_id < _orig_to_ls.size(); ++orig_id){
        const int ls_id = _orig_to_ls(orig_id);
        if(ls_id >= 0) _ls_to_orig[ls_id] = static_cast<int>(orig_id);
    }
}

// NB deliberately NOT noexcept -- it allocates, see the declaration in LSGrid.hpp.
void LSGrid::set_ls_to_orig_internal(const Eigen::Ref<const IntVect> & ls_to_orig){
    if(ls_to_orig.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }
    
    _ls_to_orig = ls_to_orig;
    const auto size = ls_to_orig.lpNorm<Eigen::Infinity>();
    _orig_to_ls = IntVect::Constant(size + 1, -1);
    int i = 0;
    for(auto el : _ls_to_orig){
        if(el != -1) _orig_to_ls[el] = i;
        ++i;
    }
}

//init
void LSGrid::init_bus(unsigned int n_sub,
                         unsigned int n_busbar_per_sub,
                         const Eigen::Ref<const RealVect> & bus_vn_kv,
                         int /*nb_line*/,
                         int /*nb_trafo*/){
    /**
    initialize the bus_vn_kv_ member
    and
    initialize the Ybus_ matrix at the proper shape
    **/
    n_sub_ = n_sub;
    max_nb_bus_per_sub_ = n_busbar_per_sub;
    substations_.init_bus(n_sub_, max_nb_bus_per_sub_,  bus_vn_kv);
    _orig_to_ls = IntVect();
    _ls_to_orig = IntVect();
}

void LSGrid::reset(bool reset_solver, bool reset_ac, bool reset_dc)
{
    if(reset_ac){
        id_me_to_ac_solver_ = SolverBusIdVect();
        id_ac_solver_to_me_ = GlobalBusIdVect();
        slack_bus_id_ac_solver_ = SolverBusIdVect();
        Ybus_ac_ = Eigen::SparseMatrix<cplx_type>();
    }

    if(reset_dc){
        id_me_to_dc_solver_ = SolverBusIdVect();
        id_dc_solver_to_me_ = GlobalBusIdVect();
        slack_bus_id_dc_solver_ = SolverBusIdVect();
        Bbus_dc_ = Eigen::SparseMatrix<real_type>();
    }

    timer_last_ac_pf_= 0.;
    timer_last_dc_pf_ = 0.;

    acSbus_ = CplxVect();
    dcPbus_ = RealVect();
    bus_pv_ = SolverBusIdVect();
    bus_pq_ = SolverBusIdVect();

    algo_controler_.ac_algo_controler().tell_all_changed();
    algo_controler_.dc_algo_controler().tell_all_changed();
    tell_solver_need_reset(); // also handles last_bus_status_saved_
    
    slack_bus_id_ac_me_ = GlobalBusIdVect();  // slack bus id, gridmodel number
    slack_bus_id_ac_solver_ = SolverBusIdVect();  // slack bus id, solver number
    slack_bus_id_dc_me_ = GlobalBusIdVect();
    slack_bus_id_dc_solver_ = SolverBusIdVect();
    slack_weights_ = RealVect();

    // reset the solvers
    if (reset_solver){
        _algo.reset();
        _algo.set_lsgrid(this);
        _algo.tell_solver_control(algo_controler_.ac_algo_controler());

        _dc_algo.reset();
        _dc_algo.set_lsgrid(this);
        _dc_algo.tell_solver_control(algo_controler_.dc_algo_controler());
    }
}

CplxVect LSGrid::ac_pf(const Eigen::Ref<const CplxVect> & Vinit,
                          int max_iter,
                          real_type tol)
{
    auto timer = CustTimer();
    const int nb_bus = static_cast<int>(substations_.nb_bus());
    if(Vinit.size() != nb_bus){
        std::ostringstream exc_;
        exc_ << "LSGrid::ac_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    BaseAlgo::check_iter_tol("LSGrid::ac_pf", max_iter, tol);
    if(hvdc_lines_.has_droop_active() && !_algo.supports_hvdc_droop()){
        std::ostringstream exc_;
        exc_ << "LSGrid::ac_pf: the grid counts hvdc lines with the angle-droop (AC emulation) enabled, ";
        exc_ << "which is only supported by the Newton-Raphson algorithms (not by the Gauss-Seidel / ";
        exc_ << "fast-decoupled ones). Please `change_algorithm` or disable the droop.";
        throw std::runtime_error(exc_.str());
    }
    bool conv = false;
    CplxVect res = CplxVect();

    // reset_results();  // clear the results  No need to do it, results are neceassirly set or reset in post process

    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    bool is_ac = true;
    CplxVect V = pre_process_solver(Vinit, 
                                    acSbus_,
                                    Ybus_ac_,
                                    id_me_to_ac_solver_,
                                    id_ac_solver_to_me_,
                                    slack_bus_id_ac_me_,
                                    slack_bus_id_ac_solver_,
                                    is_ac,
                                    algo_controler_.ac_algo_controler());

    // start the solver
    conv = _algo.compute_pf(
        Ybus_ac_,
        V,
        acSbus_,
        slack_bus_id_ac_solver_.as_eigen(),  // was _to_intvect()
        slack_weights_,
        bus_pv_.as_eigen(),  // was _to_intvect()
        bus_pq_.as_eigen(),  // was _to_intvect()
        max_iter,
        tol / sn_mva_);

    // store results (in ac mode) 
    process_results(conv, res, Vinit, true, id_me_to_ac_solver_);

    timer_last_ac_pf_ = timer.duration();
    // return the vector of complex voltage at each bus
    return res;
};

void LSGrid::fill_hvdc_droop_solver_data(HvdcDroopSolverData & data, bool ac) const
{
    data.clear();
    const int nb_hvdc = static_cast<int>(hvdc_lines_.nb());
    if(nb_hvdc == 0) return;
    const SolverBusIdVect & id_me_to_solver = ac ? id_me_to_ac_solver_ : id_me_to_dc_solver_;
    // upper bound of the solver bus ids emitted below: they are handed straight
    // to the solver, which uses them to index its n_bus-sized tables with no
    // bounds check of its own (release wheels are -O3 -DNDEBUG)
    const int nb_bus_solver = static_cast<int>((ac ? id_ac_solver_to_me_ : id_dc_solver_to_me_).size());
    const std::vector<bool> & droop_enabled = hvdc_lines_.get_droop_enabled();
    const std::vector<bool> & status_global = hvdc_lines_.get_status_global();
    std::vector<int> indices;
    indices.reserve(nb_hvdc);
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        if(!droop_enabled[hvdc_id] || !status_global[hvdc_id]) continue;
        indices.push_back(hvdc_id);
    }
    const int nb_droop = static_cast<int>(indices.size());
    if(nb_droop == 0) return;
    data.bus1 = Eigen::VectorXi(nb_droop);
    data.bus2 = Eigen::VectorXi(nb_droop);
    data.status = Eigen::VectorXi(nb_droop);
    data.p0 = RealVect(nb_droop);
    data.k = RealVect(nb_droop);
    data.lf1 = RealVect(nb_droop);
    data.lf2 = RealVect(nb_droop);
    data.r = RealVect(nb_droop);
    data.pmax12 = RealVect(nb_droop);
    data.pmax21 = RealVect(nb_droop);
    data.connected1.assign(nb_droop, true);
    data.connected2.assign(nb_droop, true);
    for(int pos = 0; pos < nb_droop; ++pos){
        const int hvdc_id = indices[pos];
        // angle-droop ("AC emulation") needs both remote angles: this must never
        // happen (both call sites that can half-open an hvdc line --
        // LSGrid::deactivate_dcline_side1/2 and
        // HvdcLineContainer::disconnect_if_not_in_main_component -- call
        // disable_droop at the same time), but enforce it explicitly rather than
        // silently indexing id_me_to_solver with the open side's bus id (-1) below.
        if(!hvdc_lines_.get_connected_side_1(hvdc_id) || !hvdc_lines_.get_connected_side_2(hvdc_id)){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_hvdc_droop_solver_data: hvdc line with id ";
            exc_ << hvdc_id;
            exc_ << " has angle-droop enabled while half-open (one side "
                    "disconnected) -- this should never happen, disable_droop "
                    "must be called whenever a side is opened.";
            throw std::runtime_error(exc_.str());
        }
        const GlobalBusId bus_1 = hvdc_lines_.get_bus_side_1(hvdc_id);
        const GlobalBusId bus_2 = hvdc_lines_.get_bus_side_2(hvdc_id);
        const SolverBusId bus_1_solver = id_me_to_solver[bus_1.cast_int()];
        const SolverBusId bus_2_solver = id_me_to_solver[bus_2.cast_int()];
        if((bus_1_solver.cast_int() == GenericContainer::_deactivated_bus_id) ||
           (bus_2_solver.cast_int() == GenericContainer::_deactivated_bus_id)){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_hvdc_droop_solver_data: hvdc line with id ";
            exc_ << hvdc_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        // The solver indexes Va / the per-bus mismatch and its own bus-keyed
        // row/column tables with these ids, unchecked. `id_me_to_solver` should
        // never produce one out of range, so this is a postcondition of the bus
        // mapping rather than an input check -- but it is the last place that
        // can still turn a broken mapping into an exception instead of an
        // out-of-bounds access inside the Newton loop.
        if((bus_1_solver.cast_int() < 0) || (bus_1_solver.cast_int() >= nb_bus_solver) ||
           (bus_2_solver.cast_int() < 0) || (bus_2_solver.cast_int() >= nb_bus_solver)){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_hvdc_droop_solver_data: hvdc line with id " << hvdc_id
                 << " resolves to solver bus ids (" << bus_1_solver.cast_int() << ", "
                 << bus_2_solver.cast_int() << ") outside [0, " << nb_bus_solver
                 << "). The grid -> solver bus mapping is inconsistent.";
            throw std::runtime_error(exc_.str());
        }
        data.bus1(pos) = bus_1_solver.cast_int();
        data.bus2(pos) = bus_2_solver.cast_int();
        data.status(pos) = hvdc_lines_.get_status_droop(hvdc_id);
        data.p0(pos) = hvdc_lines_.get_droop_p0_mw(hvdc_id) / sn_mva_;
        data.k(pos) = hvdc_lines_.get_droop_k_mw_per_rad(hvdc_id) / sn_mva_;
        data.lf1(pos) = hvdc_lines_.get_loss_factor_side_1(hvdc_id);
        data.lf2(pos) = hvdc_lines_.get_loss_factor_side_2(hvdc_id);
        const real_type v_nom = hvdc_lines_.get_nominal_v_kv(hvdc_id);
        data.r(pos) = (v_nom > 0.) ? hvdc_lines_.get_r_ohm(hvdc_id) * sn_mva_ / (v_nom * v_nom) : 0.;
        data.pmax12(pos) = hvdc_lines_.get_pmax_1to2_mw(hvdc_id) / sn_mva_;
        data.pmax21(pos) = hvdc_lines_.get_pmax_2to1_mw(hvdc_id) / sn_mva_;
        // status_global[hvdc_id] being true only guarantees at least ONE side
        // is in the main synchronous component (see
        // disconnect_if_not_in_main_component) -- the OTHER side may be
        // individually open. bus1/bus2 solver ids stay valid either way (the
        // AC bus itself, not the converter, is what id_me_to_solver maps),
        // but the theta-dependent flow into an open side is meaningless.
        data.connected1[pos] = hvdc_lines_.get_connected_side_1(hvdc_id);
        data.connected2[pos] = hvdc_lines_.get_connected_side_2(hvdc_id);
    }
}

void fill_hvdc_droop_data_from_grid(const LSGrid * lsgrid_ptr, HvdcDroopSolverData & data, bool ac)
{
    data.clear();
    if(lsgrid_ptr != nullptr) lsgrid_ptr->fill_hvdc_droop_solver_data(data, ac);
}

void LSGrid::fill_voltage_control_solver_data(VoltageControlSolverData & data, bool ac) const
{
    data.clear();
    if(!ac) return;  // DC: no voltage control (no-op, SVC contributes nothing)
    const SolverBusIdVect & id_me_to_solver = id_me_to_ac_solver_;
    const int nb_bus_solver = static_cast<int>(id_ac_solver_to_me_.size());
    if(nb_bus_solver == 0) return;

    // Cheap pre-scan, before anything that allocates. Everything below --- two
    // nb_bus_solver-sized vectors and the `free_vm_slack` set --- is only there
    // to validate controllers, so a grid with none of them (the common case, and
    // every grid in the benchmark suite) has nothing to do here. This function
    // runs on EVERY solve, so without this early exit that allocation was paid
    // per powerflow to produce an empty result.
    {
        bool any_controller = false;
        const int nb_gen_scan = static_cast<int>(generators_.nb());
        for(int gen_id = 0; gen_id < nb_gen_scan; ++gen_id){
            if(generators_.gen_is_voltage_controller(gen_id)){ any_controller = true; break; }
        }
        if(!any_controller){
            const int nb_svc_scan = static_cast<int>(svcs_.nb());
            for(int svc_id = 0; svc_id < nb_svc_scan; ++svc_id){
                if(svcs_.svc_is_voltage_controller(svc_id)){ any_controller = true; break; }
            }
        }
        if(!any_controller) return;
    }

    // PQ membership: a bus owns a Q equation AND a Vm unknown iff it is a PQ bus
    // (PV buses have only theta/P, the slack none). bus_pq_ is set by fillpv_pq.
    std::vector<bool> is_pq(nb_bus_solver, false);
    for(int k = 0; k < static_cast<int>(bus_pq_.size()); ++k){
        const int b = bus_pq_(k).cast_int();
        if(b >= 0 && b < nb_bus_solver) is_pq[b] = true;
    }
    // Slack buses are not PQ in the base block, but a slack bus that is not pinned
    // by a local PV generator is given a Q equation + free Vm by the MultiSlack
    // extension (see LSGrid::get_free_vm_slack_solver_buses), so a controller on
    // such a slack bus is supported even though `is_pq` is false there. A slack
    // bus that IS locally pinned (another generator regulates it directly) gets
    // no such Q equation at all -- checking membership of the whole `slack_bus_
    // id_ac_solver_` list here (as opposed to just this "free" subset) would
    // wrongly accept that case: its Q equation lookup then resolves to -1, the
    // controller's own reactive-injection column ends up with no Jacobian entry
    // anywhere, and the factorization fails with ErrorType::SolverFactor instead
    // of this function's own clear error.
    const std::set<int> free_vm_slack = get_free_vm_slack_solver_buses();
    std::vector<bool> has_free_q(nb_bus_solver, false);
    for(int b : free_vm_slack){
        if(b >= 0 && b < nb_bus_solver) has_free_q[b] = true;
    }

    // 1. collect the active voltage-mode controllers (remote-regulating gens for
    //    now; voltage-mode SVCs are added with the SvcContainer). Per controller:
    //    solver bus, regulated solver bus, v_set (pu), sharing key, kind, elem id.
    struct Raw { int bus; int reg_bus; real_type v_set; real_type slope; real_type weight; int kind; int elem_id; };
    std::vector<Raw> raws;
    const int nb_gen = static_cast<int>(generators_.nb());
    const GlobalBusIdVect & gen_buses = generators_.get_buses();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators_.gen_is_voltage_controller(gen_id)) continue;
        const int ctrl_grid = gen_buses(gen_id).cast_int();
        const int reg_grid  = generators_.get_regulated_bus_id(gen_id);
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        const int reg_solver  = (reg_grid >= 0) ? id_me_to_solver[reg_grid].cast_int()
                                                : GenericContainer::_deactivated_bus_id;
        if(ctrl_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " is a voltage controller but its bus is disconnected.";
            throw std::runtime_error(exc_.str());
        }
        if(reg_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates a disconnected bus.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[ctrl_solver] && !has_free_q[ctrl_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates a remote bus but its OWN bus has no reactive (Q) equation"
                    " (it is a PV bus that is not a slack, or a slack bus already locally"
                    " pinned by another voltage-regulating generator). This is not supported"
                    " in v1.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[reg_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates bus " << reg_grid << " which has no voltage (Vm) unknown"
                    " (it is a slack or an already-PV bus). This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        const real_type w = generators_.get_max_q(gen_id) - generators_.get_min_q(gen_id);
        raws.push_back({ctrl_solver, reg_solver, generators_.get_target_vm_pu(gen_id),
                        static_cast<real_type>(0.), w, VoltageControlSolverData::GEN, gen_id});
    }

    // ... and the active VOLTAGE-mode SVCs (local or remote, with or without slope)
    const int nb_svc = static_cast<int>(svcs_.nb());
    const GlobalBusIdVect & svc_buses = svcs_.get_buses();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svcs_.svc_is_voltage_controller(svc_id)) continue;
        const int ctrl_grid = svc_buses(svc_id).cast_int();
        const int reg_grid  = svcs_.get_regulated_bus_id(svc_id);
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        const int reg_solver  = (reg_grid >= 0) ? id_me_to_solver[reg_grid].cast_int()
                                                : GenericContainer::_deactivated_bus_id;
        if(ctrl_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " is a voltage controller but its bus is disconnected.";
            throw std::runtime_error(exc_.str());
        }
        if(reg_solver == GenericContainer::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " regulates a disconnected bus.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[ctrl_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " is at a bus with no reactive (Q) equation (it is a slack or PV bus)."
                    " This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_pq[reg_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << svc_id
                 << " regulates bus " << reg_grid << " which has no voltage (Vm) unknown"
                    " (it is a slack or an already-PV bus). This is not supported in v1.";
            throw std::runtime_error(exc_.str());
        }
        const real_type w = svcs_.get_b_max(svc_id) - svcs_.get_b_min(svc_id);
        raws.push_back({ctrl_solver, reg_solver, svcs_.get_target_vm_pu(svc_id),
                        svcs_.get_slope_pu(svc_id), w, VoltageControlSolverData::SVC, svc_id});
    }
    if(raws.empty()) return;

    // 2. group by regulated solver bus (merge gens that share a regulated bus),
    //    checking the v_set agree within tolerance.
    std::vector<int> grp_reg;
    std::vector<real_type> grp_vset;
    std::vector<std::vector<int> > grp_members;  // indices into raws
    for(int i = 0; i < static_cast<int>(raws.size()); ++i){
        int g = -1;
        for(int gg = 0; gg < static_cast<int>(grp_reg.size()); ++gg)
            if(grp_reg[gg] == raws[i].reg_bus){ g = gg; break; }
        if(g == -1){
            g = static_cast<int>(grp_reg.size());
            grp_reg.push_back(raws[i].reg_bus);
            grp_vset.push_back(raws[i].v_set);
            grp_members.push_back(std::vector<int>());
        } else if(std::abs(grp_vset[g] - raws[i].v_set) > BaseConstants::_tol_equal_float){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: several controllers regulate the"
                    " same bus with conflicting voltage setpoints (" << grp_vset[g] << " vs "
                 << raws[i].v_set << " pu).";
            throw std::runtime_error(exc_.str());
        }
        grp_members[g].push_back(i);
    }

    // 2b. v1 restriction: an SVC may only be ALONE in its control group. The
    //     cross-weight sharing of an SVC with other controllers (and any sloped
    //     SVC sharing a regulated bus, cf Phase 0 probe #3) is not supported yet.
    for(int g = 0; g < static_cast<int>(grp_members.size()); ++g){
        if(grp_members[g].size() <= 1) continue;
        for(int idx : grp_members[g]){
            if(raws[idx].kind == VoltageControlSolverData::SVC){
                std::ostringstream exc_;
                exc_ << "LSGrid::fill_voltage_control_solver_data: SVC " << raws[idx].elem_id
                     << " shares a regulated bus with other controllers, which is not"
                        " supported in v1 (an SVC must be the only controller of its bus).";
                throw std::runtime_error(exc_.str());
            }
        }
    }

    // 3. emit, controllers grouped contiguously
    const int ng = static_cast<int>(grp_reg.size());
    const int nc = static_cast<int>(raws.size());
    data.bus = Eigen::VectorXi(nc);
    data.kind = Eigen::VectorXi(nc);
    data.elem_id = Eigen::VectorXi(nc);
    data.slope = RealVect(nc);
    data.weight = RealVect(nc);
    data.group = Eigen::VectorXi(nc);
    data.reg_bus = Eigen::VectorXi(ng);
    data.v_set = RealVect(ng);
    data.grp_start = Eigen::VectorXi(ng);
    data.grp_count = Eigen::VectorXi(ng);
    int cursor = 0;
    for(int g = 0; g < ng; ++g){
        data.reg_bus(g) = grp_reg[g];
        data.v_set(g) = grp_vset[g];
        data.grp_start(g) = cursor;
        data.grp_count(g) = static_cast<int>(grp_members[g].size());
        for(int idx : grp_members[g]){
            const Raw & r = raws[idx];
            data.bus(cursor) = r.bus;
            data.kind(cursor) = r.kind;
            data.elem_id(cursor) = r.elem_id;
            data.slope(cursor) = r.slope;
            // floor the sharing key to keep the N>1 sharing rows non-singular
            data.weight(cursor) = (std::abs(r.weight) > BaseConstants::_tol_equal_float) ? r.weight : BaseConstants::_tol_equal_float;
            data.group(cursor) = g;
            ++cursor;
        }
    }

    // Postconditions of this function, checked once here rather than by the
    // solver on every solve. The VoltageControl extension walks each group as
    // [grp_start(g), grp_start(g) + grp_count(g)) and sizes its reactive-sharing
    // vectors from `grp_count(g) - 1` as a std::size_t, then indexes Vm / the
    // per-bus mismatch and its own bus-keyed tables with `bus` / `reg_bus` --
    // all unchecked (release wheels are -O3 -DNDEBUG). An empty group would
    // underflow that subtraction; a bus id out of range would read or write
    // past the ledger. Neither can happen given the construction above, which is
    // exactly why this belongs here: it is cheap, it runs only when the data is
    // actually rebuilt (the extensions cache it otherwise), and it keeps the
    // guarantee at the point the data is produced instead of at the point it is
    // consumed.
    if(cursor != nc){
        std::ostringstream exc_;
        exc_ << "LSGrid::fill_voltage_control_solver_data: the control groups cover " << cursor
             << " controllers but " << nc << " were collected.";
        throw std::runtime_error(exc_.str());
    }
    for(int g = 0; g < ng; ++g){
        if(data.grp_count(g) < 1){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: control group " << g
                 << " has no controller.";
            throw std::runtime_error(exc_.str());
        }
        if((data.reg_bus(g) < 0) || (data.reg_bus(g) >= nb_bus_solver)){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: control group " << g
                 << " regulates solver bus " << data.reg_bus(g) << ", outside [0, "
                 << nb_bus_solver << "). The grid -> solver bus mapping is inconsistent.";
            throw std::runtime_error(exc_.str());
        }
    }
    for(int j = 0; j < nc; ++j){
        if((data.bus(j) < 0) || (data.bus(j) >= nb_bus_solver)){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: controller " << j
                 << " sits on solver bus " << data.bus(j) << ", outside [0, "
                 << nb_bus_solver << "). The grid -> solver bus mapping is inconsistent.";
            throw std::runtime_error(exc_.str());
        }
    }
}

std::set<int> LSGrid::get_free_vm_slack_solver_buses() const
{
    std::set<int> res;
    // solver-bus ids of the slack buses
    std::set<int> slack;
    for(int k = 0; k < static_cast<int>(slack_bus_id_ac_solver_.size()); ++k){
        slack.insert(slack_bus_id_ac_solver_(k).cast_int());
    }
    if(slack.empty()) return res;

    // A slack bus is Vm-fixed (PV-like, no Q equation) only when a LOCAL
    // voltage-regulating generator pins its magnitude. Collect those buses.
    std::set<int> locally_vfixed;
    const SolverBusIdVect & id_me_to_solver = id_me_to_ac_solver_;
    const int nb_gen = static_cast<int>(generators_.nb());
    const GlobalBusIdVect & gen_buses = generators_.get_buses();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!generators_.gen_is_local_voltage_controller(gen_id)) continue;
        const int ctrl_grid = gen_buses(gen_id).cast_int();
        const int ctrl_solver = id_me_to_solver[ctrl_grid].cast_int();
        if(ctrl_solver == GenericContainer::_deactivated_bus_id) continue;
        locally_vfixed.insert(ctrl_solver);
    }

    // Every slack bus whose magnitude is NOT pinned locally needs a free Vm
    // unknown + Q equation: distributed-slack PQ participants (the common case),
    // remote-voltage controllers, and SVC-regulated slack buses all fall here.
    for(int b : slack){
        if(!locally_vfixed.count(b)) res.insert(b);
    }
    return res;
}

void LSGrid::check_solution_q_values_onegen(Eigen::Ref<CplxVect> res,
                                               int bus_id,
                                               real_type min_q_mvar,
                                               real_type max_q_mvar,
                                               bool check_q_limits) const{
    if(check_q_limits)
    {
        // i need to check the reactive can be absorbed / produced by the generator
        real_type new_q = BaseConstants::my_zero_;
        real_type react_this_bus = std::imag(res.coeff(bus_id));
        if((react_this_bus >= min_q_mvar) && (react_this_bus <= max_q_mvar))
        {
            // this generator is able to handle all reactive
            new_q = BaseConstants::my_zero_;
        }else if(react_this_bus < min_q_mvar){
            // generator cannot absorb enough reactive power
            new_q = react_this_bus - min_q_mvar; //ex. need -50, qmin is -30, remains: (-50) - (-30) = -20 MVAr
        }else{
            // generator cannot produce enough reactive power
            new_q = react_this_bus - max_q_mvar;  // ex. need 50, qmax is 30, remains: 50 - 30 = 20 MVAr
        }
        res.coeffRef(bus_id) = {std::real(res.coeff(bus_id)), new_q};
    }else{
        // the q value for the bus at which the generator is connected will be 0
        res.coeffRef(bus_id) = {std::real(res.coeff(bus_id)), BaseConstants::my_zero_};
    }
}

void LSGrid::check_solution_q_values(Eigen::Ref<CplxVect> res, bool check_q_limits) const{
    // test for iterator though generators
    for(const auto & gen: generators_)
    {
        if(!gen.connected)
        {
            // the generator is disconnected, I do nothing
            continue;
        }
        // Only a voltage-regulating generator has a genuinely FREE Q at its own bus in
        // the real NR system (GeneratorContainer::fillSbus never stamps Q there, local
        // or remote control alike -- see fillpv / the VoltageControl extension). A
        // non-regulating (fixed PQ) generator's Q is deterministic and already part of
        // Sbus, so masking it here would hide a real mismatch instead of correctly
        // reporting "this bus's Q is free, don't judge it".
        if(gen.voltage_regulator_on){
            check_solution_q_values_onegen(res, gen.bus_id, gen.min_q_mvar, gen.max_q_mvar, check_q_limits);
        }

        // if(gen.id == gen_slackbus_)
        if(gen.is_slack)
        {
            // slack bus, by definition, can handle all active value
            // This is probably not the case with distributed slack !
            res.coeffRef(gen.bus_id) = {BaseConstants::my_zero_, std::imag(res.coeff(gen.bus_id))};
        }
    }

    // then do the same for the hvdc converter stations
    for(const auto & hvdc: hvdc_lines_)
    {
        if(!hvdc.connected_global)
        {
            // the hvdc line is disconnected, I do nothing
            continue;
        }
        const auto & station_1 = hvdc.station_side_1;
        const auto & station_2 = hvdc.station_side_2;
        // a side may be open while the line is still connected_global (a line whose
        // remote converter is in another synchronous component, see
        // HvdcLineContainer::disconnect_if_not_in_main_component): skip the open side.
        // Same "only mask a genuinely free Q" reasoning as for generators above: a
        // fixed-Q / power-factor (LCC) station's Q is real Sbus data, not a free
        // variable -- see ConverterStationContainer::fillSbus_station.
        if(station_1.connected && station_1.voltage_regulator_on)
            check_solution_q_values_onegen(res, station_1.bus_id, station_1.min_q_mvar, station_1.max_q_mvar, check_q_limits);
        if(station_2.connected && station_2.voltage_regulator_on)
            check_solution_q_values_onegen(res, station_2.bus_id, station_2.min_q_mvar, station_2.max_q_mvar, check_q_limits);
    }

    // ... and the VOLTAGE-mode static var compensators: their reactive injection is
    // solved by the VoltageControl extension (not stamped in Sbus), so the reactive
    // mismatch at their bus must be absorbed. SVC reactive limits are NEVER enforced.
    for(const auto & svc: svcs_)
    {
        if(!svc.connected) continue;
        if(svc.regulation_mode != SvcContainer::RegulationMode::VOLTAGE) continue;  // REACTIVE_POWER is in Sbus
        check_solution_q_values_onegen(res, svc.bus_id, 0., 0., false);  // never enforce SVC limits
    }
}

CplxVect LSGrid::check_solution(const Eigen::Ref<const CplxVect> & V_proposed, bool check_q_limits)
{
    // pre process the data to define a proper jacobian matrix, the proper voltage vector etc.
    const int nb_bus = static_cast<int>(V_proposed.size());
    bool is_ac = true;
    AlgoControl reset_solver;
    // `AlgoControl`'s own default constructor already asks for a full rebuild
    // (`need_reset_solver_=true`) -- only downgrade that to "nothing changed" once the
    // AC solver bus mapping has actually been built at least once (by a prior
    // `ac_pf`/`dc_pf`/`check_solution` call). Calling `check_solution` as the very
    // first operation on a freshly-built model with `tell_none_changed()`
    // unconditionally used to leave `id_me_to_ac_solver_`/`Ybus_ac_` at their
    // default-constructed (empty) size, and `fill_hvdc_droop_solver_data`'s
    // `id_me_to_solver[bus_id]` lookup (bus ids in the hundreds/thousands) then read
    // out of bounds -- a silent, hard-to-reproduce segfault instead of a clean rebuild.
    if(id_me_to_ac_solver_.size() > 0) reset_solver.tell_none_changed();
    CplxVect V = pre_process_solver(V_proposed,
                                    acSbus_,
                                    Ybus_ac_,
                                    id_me_to_ac_solver_,
                                    id_ac_solver_to_me_,
                                    slack_bus_id_ac_me_,
                                    slack_bus_id_ac_solver_,
                                    is_ac, reset_solver,
                                    false);  // do NOT snap regulated buses to their target: we are testing V_proposed as-is

    // compute the mismatch
    CplxVect tmp = Ybus_ac_ * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    CplxVect mis = V.array() * tmp.array() - acSbus_.array();  // TODO ac or dc here

    // the angle-droop (AC emulation) hvdc flows are not part of Sbus: they
    // leave the buses, so they add to the computed power, ie to the mismatch
    HvdcDroopSolverData droop_data;
    fill_hvdc_droop_solver_data(droop_data, true);
    for(int k = 0; k < droop_data.size(); ++k){
        const real_type theta_1 = std::arg(V(droop_data.bus1(k)));
        const real_type theta_2 = std::arg(V(droop_data.bus2(k)));
        const real_type raw = droop_data.p0(k) + droop_data.k(k) * (theta_1 - theta_2);
        real_type p1_flow, p2_flow;
        droop_data.flows_pu(k, raw, p1_flow, p2_flow);
        mis(droop_data.bus1(k)) += p1_flow;
        mis(droop_data.bus2(k)) += p2_flow;
    }

    // store results
    CplxVect res = _get_results_back_to_orig_nodes(mis,
                                                   id_me_to_ac_solver_,
                                                   static_cast<int>(V_proposed.size())
                                                   );
    if(abs(sn_mva_- 1.) > BaseConstants::_tol_equal_float) res *= sn_mva_;

    // now check reactive values for buses where there are generators and active values of slack bus
    check_solution_q_values(res, check_q_limits);

    // set to 0 the error on the disconnected bus (it is not initialized at 0.0 in _get_results_back_to_orig_nodes)
    for(int bus_id = 0; bus_id < nb_bus; ++bus_id)
    {
        if(substations_.is_bus_connected(GlobalBusId(bus_id))) continue;
        res.coeffRef(bus_id) = BaseConstants::my_zero_;
    }
    return res;
};

// AC injection: complex Sbus + the reactive-power vectors (Q limits / gen count per bus)
// Sbus stays a plain reference (not Eigen::Ref): it is reassigned below
// (Sbus = CplxVect::Constant(...)) to a size that can change with topology.
void LSGrid::prepare_injection(CplxVect & Sbus, bool redo_all, bool converter_changed,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const AlgoControl & solver_control)
{
    if (redo_all || converter_changed || solver_control.need_recompute_sbus()) {
        // init Sbus
        Sbus = CplxVect::Constant(id_solver_to_me.size(), 0.);
    }
    if (redo_all ||
        solver_control.need_recompute_sbus() ||  // TODO do we need it ?
        solver_control.has_slack_participate_changed() ||
        solver_control.has_pv_changed() ||
        solver_control.has_pq_changed())  // TODO do we need it ?
    {
        int nb_bus_total = static_cast<int>(substations_.nb_bus());
        total_q_min_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
        total_q_max_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
        total_gen_per_bus_ = Eigen::VectorXi::Constant(nb_bus_total, 0);
        generators_.init_q_vector(nb_bus_total, total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
        hvdc_lines_.init_q_vector(nb_bus_total, total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
    }
    if (redo_all || converter_changed ||
        solver_control.has_slack_participate_changed() ||
        solver_control.has_pv_changed() ||
        solver_control.has_pq_changed() ||
        solver_control.need_recompute_sbus()) {
            fillSbus_me(Sbus, true, id_me_to_solver);
        }
}

// DC injection: real Pbus, assembled by reusing the (complex) Sbus fills and keeping the real part
// Pbus stays a plain reference (not Eigen::Ref): it is reassigned below
// (Pbus = Sbus_tmp.real()) to a size that can change with topology.
void LSGrid::prepare_injection(RealVect & Pbus, bool redo_all, bool converter_changed,
                               const SolverBusIdVect & id_me_to_solver,
                               const GlobalBusIdVect & id_solver_to_me,
                               const AlgoControl & solver_control)
{
    if (redo_all || converter_changed ||
        solver_control.has_slack_participate_changed() ||
        solver_control.has_pv_changed() ||
        solver_control.has_pq_changed() ||
        solver_control.need_recompute_sbus()) {
            CplxVect Sbus_tmp = CplxVect::Constant(id_solver_to_me.size(), 0.);
            fillSbus_me(Sbus_tmp, false, id_me_to_solver);
            Pbus = Sbus_tmp.real();
        }
}

template<class MatScalar, class InjVect>
CplxVect LSGrid::_pre_process_solver_impl(
    const Eigen::Ref<const CplxVect> & Vinit,
    InjVect & inj,
    Eigen::SparseMatrix<MatScalar> & mat,
    SolverBusIdVect & id_me_to_solver,
    GlobalBusIdVect & id_solver_to_me,
    GlobalBusIdVect & slack_bus_id_me,
    SolverBusIdVect & slack_bus_id_solver,
    const AlgoControl & solver_control,
    bool init_pv_vm_targets)
{
    // cplx_type matrix => AC solver family, real_type matrix => DC solver family
    const bool is_ac = std::is_same<MatScalar, cplx_type>::value;
    if(solver_control.need_reset_solver()){
        if(is_ac) _algo.reset();
        else _dc_algo.reset();
    }

    bool redo_all =
            solver_control.need_reset_solver() ||
            solver_control.has_dimension_changed();

    if (redo_all || solver_control.has_slack_participate_changed()){
        slack_bus_id_me = generators_.get_slack_bus_id();
        // this is the slack bus ids with the gridmodel ordering, not the solver ordering.
        // conversion to solver ordering is done in init_slack_bus

        // Optional forced angle reference: move the requested (gridmodel) bus to
        // the front so the NR uses it as slack_ids[0] (the reference) without
        // changing the slack set or weights. See LSGrid::set_reference_slack_bus.
        if (_forced_ref_slack_bus_id >= 0){
            std::vector<int> sids = slack_bus_id_me.to_int_vector();
            for (std::size_t i = 1; i < sids.size(); ++i){
                if (sids[i] == _forced_ref_slack_bus_id){
                    const int ref = sids[i];
                    sids.erase(sids.begin() + static_cast<std::ptrdiff_t>(i));
                    sids.insert(sids.begin(), ref);
                    slack_bus_id_me = GlobalBusIdVect(sids);
                    break;
                }
            }
        }
    }
    if (redo_all || solver_control.has_one_el_changed_bus()){
        init_bus_status();
    }

    // init_bus_status can set the flag "has_dimension_change", so redo this here
    redo_all =
            solver_control.need_reset_solver() ||
            solver_control.has_dimension_changed();
    bool converter_changed = false;
    if (redo_all || solver_control.ybus_change_sparsity_pattern()){
        init_converter_bus_id(id_me_to_solver, id_solver_to_me);
        const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
        init_solver_matrix(mat, nb_bus_solver);
        converter_changed = true;
    }
    if (redo_all || converter_changed || solver_control.need_recompute_ybus()){
        fill_solver_matrix(mat, id_me_to_solver);
    }
    if (redo_all || converter_changed ||
        solver_control.has_slack_participate_changed() ||
        solver_control.has_pv_changed() ||
        solver_control.has_pq_changed()) {
            init_slack_bus(id_me_to_solver, id_solver_to_me, slack_bus_id_me, slack_bus_id_solver);
            fillpv_pq(id_me_to_solver, id_solver_to_me, slack_bus_id_solver, solver_control);
        }

    // type-specific injection assembly (complex Sbus for AC, real Pbus for DC)
    prepare_injection(inj, redo_all, converter_changed, id_me_to_solver, id_solver_to_me, solver_control);

    const int nb_bus_solver = static_cast<int>(id_solver_to_me.size());
    CplxVect V = CplxVect::Constant(nb_bus_solver, init_vm_pu_);
    for(int bus_solver_id = 0; bus_solver_id < nb_bus_solver; ++bus_solver_id){
        GlobalBusId bus_me_id = id_solver_to_me[bus_solver_id];
        if(bus_me_id.cast_int() == BaseConstants::_deactivated_bus_id){
            //TODO DEBUG MODE : only in debug mode
            std::ostringstream exc_;
            exc_ << "LSGrid::pre_process_solver: the bus with solver id ";
            exc_ << bus_solver_id;
            exc_ << " is connected, but mapped (in id_solver_to_me) to a disconnected bus (global / gridmodel id)";
            throw std::runtime_error(exc_.str());
        }
        V(bus_solver_id) = Vinit(bus_me_id.cast_int());
    }
    if(init_pv_vm_targets){
        // NR-initialization heuristic only: snaps regulated buses with no droop/slope
        // to their own target voltage magnitude. Skipped by check_solution, which must
        // evaluate the caller-supplied voltage as given (see the `init_pv_vm_targets`
        // doc on `pre_process_solver`).
        generators_.set_vm(V, id_me_to_solver);
        hvdc_lines_.set_vm(V, id_me_to_solver);
        svcs_.set_vm(V, id_me_to_solver);  // VOLTAGE-mode SVCs (init quality at the regulated bus)
    }

    if(solver_control.need_reset_solver() ||
       solver_control.has_dimension_changed() ||
       solver_control.has_slack_participate_changed() ||
       solver_control.has_pv_changed() ||
       solver_control.has_slack_weight_changed()){
        slack_weights_ = generators_.get_slack_weights_solver(mat.rows(), id_me_to_solver);
    }

    if(is_ac) _algo.tell_solver_control(solver_control);
    else _dc_algo.tell_solver_control(solver_control);

    // keep the member bus mapping in sync with the one we just built. The single-shot
    // ac_pf / dc_pf pass the member itself as `id_me_to_solver` (self-assign, skipped
    // below), but the batch algorithms (TimeSeries / ContingencyAnalysis) pass their own
    // local vector. The hvdc droop extension reads the member through
    // fill_hvdc_droop_solver_data(), so it must reflect the active mapping in every path.
    if(is_ac){
        if(&id_me_to_solver != &id_me_to_ac_solver_) id_me_to_ac_solver_ = id_me_to_solver;
    } else {
        if(&id_me_to_solver != &id_me_to_dc_solver_) id_me_to_dc_solver_ = id_me_to_solver;
    }
    return V;
}

CplxVect LSGrid::pre_process_solver(
    const Eigen::Ref<const CplxVect> & Vinit,
    CplxVect & Sbus,
    Eigen::SparseMatrix<cplx_type> & Ybus,
    SolverBusIdVect & id_me_to_solver,
    GlobalBusIdVect & id_solver_to_me,
    GlobalBusIdVect & slack_bus_id_me,
    SolverBusIdVect & slack_bus_id_solver,
    bool /*is_ac*/,  // kept for API compatibility; DC now goes through pre_process_dc_solver
    const AlgoControl & solver_control,
    bool init_pv_vm_targets)
{
    return _pre_process_solver_impl<cplx_type>(
        Vinit, Sbus, Ybus, id_me_to_solver, id_solver_to_me,
        slack_bus_id_me, slack_bus_id_solver, solver_control, init_pv_vm_targets);
}

CplxVect LSGrid::pre_process_dc_solver(
    const Eigen::Ref<const CplxVect> & Vinit,
    RealVect & Pbus,
    Eigen::SparseMatrix<real_type> & Bbus,
    SolverBusIdVect & id_me_to_solver,
    GlobalBusIdVect & id_solver_to_me,
    GlobalBusIdVect & slack_bus_id_me,
    SolverBusIdVect & slack_bus_id_solver,
    const AlgoControl & solver_control)
{
    return _pre_process_solver_impl<real_type>(
        Vinit, Pbus, Bbus, id_me_to_solver, id_solver_to_me,
        slack_bus_id_me, slack_bus_id_solver, solver_control, true);
}

CplxVect LSGrid::_get_results_back_to_orig_nodes(const Eigen::Ref<const CplxVect> & res_tmp,
                                                    SolverBusIdVect & id_me_to_solver,
                                                    int size)
{
    CplxVect res = CplxVect::Constant(size, {init_vm_pu_, BaseConstants::my_zero_});
    const int nb_bus = static_cast<int>(substations_.nb_bus());
    for (int bus_id_me=0; bus_id_me < nb_bus; ++bus_id_me){
        if(!substations_.is_bus_connected(GlobalBusId(bus_id_me))) continue;  // nothing is done if the bus is connected
        SolverBusId bus_id_solver = id_me_to_solver[bus_id_me];
        if(bus_id_solver.cast_int() == BaseConstants::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::_get_results_back_to_orig_nodes: the bus with id ";
            exc_ << bus_id_me;
            exc_ << " is connected to a disconnected bus (solver side)";
            throw std::runtime_error(exc_.str());
        }
        res(bus_id_me) = res_tmp(static_cast<int>(bus_id_solver));
    }
    return res;
}

void LSGrid::process_results(bool conv,
                                CplxVect & res,
                                const Eigen::Ref<const CplxVect> & Vinit,
                                bool ac,
                                SolverBusIdVect & id_me_to_solver)
{
    if (conv){
        // An external (plugin) solver can claim convergence but return malformed
        // voltages. Validate their size/finiteness before we index them below
        // (compute_results / _get_results_back_to_orig_nodes both index the solver
        // vectors with an unchecked operator()).
        // Only external solvers are checked: the built-in ones are covered by the
        // test suite, so they pay nothing here.
        //
        // "Is it built-in?" is asked of the REGISTRY, which records it at
        // registration time (SolverOrigin), not of AlgorithmType. Gating on
        // `get_type() == AlgorithmType::Custom` was wrong: AlgorithmType is a fixed
        // enum of serialized solver identities, and a built-in only has a member
        // there if one was added for it -- the NRRefactorRetry_* family never got
        // one, so name_to_algo_type() reports those three built-ins as Custom
        // exactly like a plugin, and they were paying for this check on every
        // single solve. The flag is cached by AlgorithmSelector::change_algorithm,
        // so this stays a bool read on the hot path.
        const bool is_external_algo =
            !(ac ? _algo.is_builtin_algo() : _dc_algo.is_builtin_algo());
        if (is_external_algo) conv = _check_solver_output(ac);
    }
    if (conv){
        if(compute_results_){
            // compute the results of the flows, P,Q,V of loads etc.
            compute_results(ac);
        }
        // solver_control_.tell_none_changed();  // todo automatically set for ac / dc the `tell_none_changed()`
        // was `const CplxVect & res_tmp = ...`: get_V() already returns Eigen::Ref<const
        // CplxVect>, so binding that to a concrete CplxVect& forced a full-vector copy
        // here on every single solve.
        const Eigen::Ref<const CplxVect> res_tmp = ac ? _algo.get_V(): _dc_algo.get_V() ;

        // convert back the results to "big" vector
        res = _get_results_back_to_orig_nodes(res_tmp,
                                              id_me_to_solver,
                                              static_cast<int>(Vinit.size()));
    } else {
        //powerflow diverge
        reset_results();
        // TODO solver control ??? something to do here ?
    }
}

bool LSGrid::_check_solver_output(bool ac)
{
    const Eigen::Ref<const CplxVect> V  = ac ? _algo.get_V()  : _dc_algo.get_V();
    const Eigen::Ref<const RealVect> Va = ac ? _algo.get_Va() : _dc_algo.get_Va();
    const Eigen::Ref<const RealVect> Vm = ac ? _algo.get_Vm() : _dc_algo.get_Vm();
    const int nb_bus_solver = ac ? static_cast<int>(id_ac_solver_to_me_.size())
                                 : static_cast<int>(id_dc_solver_to_me_.size());
    const char * algo_name = ac ? "AC" : "DC";

    if((V.size() != nb_bus_solver) || (Va.size() != nb_bus_solver) || (Vm.size() != nb_bus_solver))
    {
        // wrong size = the solver broke its contract. This would cause an
        // out-of-bounds read in compute_results / _get_results_back_to_orig_nodes.
        std::ostringstream exc_;
        exc_ << "LSGrid::process_results: the " << algo_name << " algorithm reported convergence but "
             << "returned voltage vectors of an unexpected size (V: " << V.size() << ", Va: " << Va.size()
             << ", Vm: " << Vm.size() << ", while the solver problem has " << nb_bus_solver
             << " buses). This is a bug in the (possibly plugin) solver.";
        throw std::runtime_error(exc_.str());
    }
    if((!V.allFinite()) || (!Va.allFinite()) || (!Vm.allFinite()))
    {
        // Non-finite voltage: a well-behaved solver reports this itself
        // (ErrorType::InifiniteValue, non-convergence); a misbehaving one may not.
        // Treat it as a non-converged solve so no NaN/Inf propagates to the results.
        (ac ? _algo : _dc_algo).set_error(ErrorType::InifiniteValue);
        return false;
    }
    return true;
}

void LSGrid::init_converter_bus_id(SolverBusIdVect& id_me_to_solver,
                                      GlobalBusIdVect& id_solver_to_me){

    //TODO get disconnected bus !!! (and have some conversion for it)
    //1. init the conversion bus
    const int nb_bus_init = static_cast<int>(substations_.nb_bus());
    id_me_to_solver = SolverBusIdVect(nb_bus_init, SolverBusId(BaseConstants::_deactivated_bus_id));  // by default, if a bus is disconnected, then it has a -1 there
    id_solver_to_me = GlobalBusIdVect();
    id_solver_to_me.reserve(nb_bus_init);
    int bus_id_solver = 0;
    for(int bus_id_me=0; bus_id_me < nb_bus_init; ++bus_id_me){
        if(substations_.is_bus_connected(GlobalBusId(bus_id_me))){
            // bus is connected
            id_solver_to_me.push_back(GlobalBusId(bus_id_me));
            id_me_to_solver[bus_id_me] = SolverBusId(bus_id_solver);
            ++bus_id_solver;
        }
    }
}

void LSGrid::init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                          int nb_bus_solver){
    Ybus = Eigen::SparseMatrix<cplx_type>(nb_bus_solver, nb_bus_solver);
    Ybus.reserve(nb_bus_solver + 4*powerlines_.nb() + 4*trafos_.nb() + 2 * shunts_.nb());
}

void LSGrid::init_Bbus(Eigen::SparseMatrix<real_type> & Bbus,
                          int nb_bus_solver){
    // DC: real admittance matrix, only lines and trafos contribute (no shunt)
    Bbus = Eigen::SparseMatrix<real_type>(nb_bus_solver, nb_bus_solver);
    Bbus.reserve(nb_bus_solver + 4*powerlines_.nb() + 4*trafos_.nb());
}

void LSGrid::init_slack_bus(const SolverBusIdVect& id_me_to_solver,
                               const GlobalBusIdVect& /*id_solver_to_me*/,
                               const GlobalBusIdVect & slack_bus_id_me,
                               SolverBusIdVect & slack_bus_id_solver)
{
    slack_bus_id_solver = SolverBusIdVect(slack_bus_id_me.size(), SolverBusId(BaseConstants::_deactivated_bus_id));

    size_t i = 0;
    for(const GlobalBusId & el: slack_bus_id_me) {
        SolverBusId tmp = id_me_to_solver[el.cast_int()];
        if(tmp.cast_int() == BaseConstants::_deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LSGrid::init_slack_bus: One of the slack bus is disconnected.";
            exc_ << " You can check bus with global id GlobalBusId : ";
            exc_ << el.cast_int();
            exc_ << ": [";
            for(const auto & el2 : slack_bus_id_me) exc_ << el2.cast_int() << ", ";
            exc_ << "].";
            throw std::out_of_range(exc_.str());
        }
        slack_bus_id_solver(i) = tmp;
        ++i;
    }
    
    if(GenericContainer::is_in_vect(BaseConstants::_deactivated_bus_id, slack_bus_id_solver.to_int_vector())){
        // TODO improve error message with the gen_id
        // TODO DEBUG MODE: only check that in debug mode
        throw std::runtime_error("LSGrid::init_Sbus: One of the slack bus is disconnected !");
    }
}
void LSGrid::fillYbus(
    Eigen::SparseMatrix<cplx_type> & res,
    bool ac,
    const SolverBusIdVect& id_me_to_solver){
    /**
    Supposes that the powerlines, shunt and transformers are initialized.
    And it fills the Ybus matrix.
    **/

    // init the Ybus matrix
    res.setZero();  // it should not be needed but might not hurt too much either.
    std::vector<Eigen::Triplet<cplx_type> > tripletList;
    tripletList.reserve(substations_.nb_bus() + 4*powerlines_.nb() + 4*trafos_.nb() + shunts_.nb());
    powerlines_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    trafos_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    loads_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    sgens_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    storages_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    generators_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    hvdc_lines_.fillYbus(tripletList, ac, id_me_to_solver, sn_mva_);
    res.setFromTriplets(tripletList.begin(), tripletList.end());  // works because  "The initial contents of *this is destroyed"
    res.makeCompressed();
}

void LSGrid::fillBdc(
    Eigen::SparseMatrix<real_type> & res,
    const SolverBusIdVect& id_me_to_solver){
    /**
    Real DC admittance (Bbus) matrix. Only lines and transformers contribute to the DC
    admittance matrix (shunts contribute to Pbus, not Bbus). Builds real triplets directly.
    **/
    res.setZero();
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(4*powerlines_.nb() + 4*trafos_.nb());
    powerlines_.fillBdc(tripletList, id_me_to_solver, sn_mva_);
    trafos_.fillBdc(tripletList, id_me_to_solver, sn_mva_);
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    res.makeCompressed();
}

void LSGrid::fillSbus_me(Eigen::Ref<CplxVect> Sbus, bool ac, const SolverBusIdVect& id_me_to_solver)
{
    // init the Sbus 
    Sbus.array() = 0.;  // reset to 0.
    powerlines_.fillSbus(Sbus, id_me_to_solver, ac);  // TODO have a function to dispatch that to all type of elements
    trafos_.fillSbus(Sbus, id_me_to_solver, ac);
    shunts_.fillSbus(Sbus, id_me_to_solver, ac);
    loads_.fillSbus(Sbus, id_me_to_solver, ac);
    sgens_.fillSbus(Sbus, id_me_to_solver, ac);
    storages_.fillSbus(Sbus, id_me_to_solver, ac);
    generators_.fillSbus(Sbus, id_me_to_solver, ac);
    hvdc_lines_.fillSbus(Sbus, id_me_to_solver, ac);
    svcs_.fillSbus(Sbus, id_me_to_solver, ac);  // REACTIVE_POWER-mode SVCs only
    if (abs(sn_mva_ - 1.0) > BaseConstants::_tol_equal_float) Sbus /= sn_mva_;
    // in dc mode, this is used for the phase shifter, this should not be divided by sn_mva_ !
    trafos_.hack_Sbus_for_dc_phase_shifter(Sbus, ac, id_me_to_solver);
}

void LSGrid::fillpv_pq(const SolverBusIdVect& id_me_to_solver,
                          const GlobalBusIdVect& id_solver_to_me,
                          const SolverBusIdVect & slack_bus_id_solver,
                          const AlgoControl & /*solver_control*/)
{
    // Nothing to do if neither pv, nor pq nor the dimension of the problem has changed

    // init pq and pv vector
    // TODO remove the order here..., i could be faster in this piece of code (looping once through the buses)
    const int nb_bus = static_cast<int>(id_solver_to_me.size());  // number of bus in the solver!
    std::vector<int> bus_pq;
    bus_pq.reserve(nb_bus);
    std::vector<int> bus_pv;
    bus_pv.reserve(nb_bus);
    std::vector<bool> has_bus_been_added(nb_bus, false);

    bus_pv_ = SolverBusIdVect();
    bus_pq_ = SolverBusIdVect();
    powerlines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    trafos_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    loads_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    storages_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    sgens_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    generators_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);
    hvdc_lines_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_me_to_solver);

    for(int bus_id = 0; bus_id< nb_bus; ++bus_id){
        if(GenericContainer::is_in_vect(bus_id, slack_bus_id_solver.to_int_vector())) continue;  // slack bus is not PQ either
        if(has_bus_been_added[bus_id]) continue; // a pv bus cannot be PQ
        bus_pq.push_back(bus_id);
        has_bus_been_added[bus_id] = true;  // don't add it a second time
    }
    bus_pv_ = SolverBusIdVect(bus_pv.size(), SolverBusId(0));
    for(int i = 0; i < static_cast<int>(bus_pv.size()); ++i){
        bus_pv_(i) = SolverBusId(bus_pv[i]);
    }
    bus_pq_ = SolverBusIdVect(bus_pq.size(), SolverBusId(0));
    for(int i = 0; i< static_cast<int>(bus_pq.size()); ++i){
        bus_pq_(i) = SolverBusId(bus_pq[i]);
    }
}

void LSGrid::compute_results(bool ac){
    // retrieve results from powerflow
    const auto & Va = ac ? _algo.get_Va() : _dc_algo.get_Va();
    const auto & Vm = ac ? _algo.get_Vm() : _dc_algo.get_Vm();
    const auto & V = ac ? _algo.get_V() : _dc_algo.get_V();

    const SolverBusIdVect & id_me_to_solver = ac ? id_me_to_ac_solver_ : id_me_to_dc_solver_;
    // for powerlines
    powerlines_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);  // TODO have a function to dispatch that to all type of elements
    // for trafo
    trafos_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for loads
    loads_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for static gen
    sgens_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for storage units
    storages_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for shunts
    shunts_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for prods
    generators_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for dclines
    hvdc_lines_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);
    // for static var compensators
    svcs_.compute_results(Va, Vm, V, id_me_to_solver, substations_.get_bus_vn_kv(), sn_mva_, ac);

    //handle_slack_bus active power
    CplxVect mismatch;  // power mismatch at each bus (SOLVER BUS !!!)
    RealVect reactive_mismatch;  // not used in dc mode (DO NOT ATTEMPT TO USE IT THERE)
    RealVect active_mismatch;
    if(ac){
        // In AC mode i am not forced to run through all the grid
        // auto tmp = (Ybus_ac_ * V).conjugate();
        mismatch = V.array() * (Ybus_ac_ * V).conjugate().array() - acSbus_.array();
        active_mismatch = mismatch.real() * sn_mva_;
    } else{
        // distributed slack (DC): the global active power imbalance -sum(Pbus) is shared among the
        // participating slack buses proportionally to their (normalized) slack weights. With a single
        // slack this assigns the whole imbalance to the reference bus (historical behaviour).
        active_mismatch = RealVect::Zero(V.size());
        const real_type imbalance = -dcPbus_.sum() * sn_mva_;
        if(slack_weights_.size() == active_mismatch.size()){
            for(int k=0; k < slack_weights_.size(); ++k){
                if(slack_weights_(k) <= BaseConstants::my_zero_) continue;
                active_mismatch(k) = slack_weights_(k) * imbalance;
            }
        } else {
            // fallback (should not happen): assign the whole imbalance to the reference slack bus
            const SolverBusId id_slack = slack_bus_id_dc_solver_(0);
            active_mismatch(id_slack.cast_int()) = imbalance;
        }
    }
    generators_.set_p_slack(active_mismatch, id_me_to_solver);

    if(ac) reactive_mismatch = mismatch.imag() * sn_mva_;
    // mainly to initialize the Q value of the generators in dc (just fill it with 0.)
    generators_.set_q(reactive_mismatch, id_me_to_solver, ac,
                      total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);
    hvdc_lines_.set_q(reactive_mismatch, id_me_to_solver, ac,
                    total_gen_per_bus_, total_q_min_per_bus_, total_q_max_per_bus_);

    // VoltageControl (remote gen + SVC) write-back: the reactive output of the
    // voltage-mode controllers is solved inside the NR system (not by the per-bus
    // redistribution above, which skips them). Pull it from the AC algorithm and
    // store it (pu -> MVAr). Empty for DC / non-NR algorithms.
    if(ac){
        const RealVect ctrl_q    = _algo.get_controller_q();
        const IntVect  ctrl_kind = _algo.get_controller_kind();
        const IntVect  ctrl_elem = _algo.get_controller_elem_id();
        for(int i = 0; i < static_cast<int>(ctrl_q.size()); ++i){
            const real_type q_mvar = ctrl_q(i) * sn_mva_;
            if(ctrl_kind(i) == VoltageControlSolverData::GEN){
                generators_.set_voltage_control_q(ctrl_elem(i), q_mvar);
            } else {  // VoltageControlSolverData::SVC
                svcs_.set_voltage_control_q(ctrl_elem(i), q_mvar);
            }
        }
    }
}

void LSGrid::reset_results(){
    powerlines_.reset_results();  // TODO have a function to dispatch that to all type of elements
    shunts_.reset_results();
    trafos_.reset_results();
    loads_.reset_results();
    sgens_.reset_results();
    storages_.reset_results();
    generators_.reset_results();
    hvdc_lines_.reset_results();
    svcs_.reset_results();
}

CplxVect LSGrid::dc_pf(const Eigen::Ref<const CplxVect> & Vinit,
                          int max_iter,  // only validated: not used for DC (single linear solve)
                          real_type tol  // only validated: not used for DC (single linear solve)
                          )
{
    //TODO SLACK: improve distributed slack for DC mode !
    // the idea is to "mess" with the Sbus beforehand to split the "losses"
    // ie fake the action of generators to adjust Sbus such that sum(Sbus) = 0
    // and the slack contribution factors are met.
    auto timer = CustTimer();

    const int nb_bus = static_cast<int>(substations_.nb_bus());
    if(Vinit.size() != nb_bus){
        //TODO DEBUG MODE: 
        std::ostringstream exc_;
        exc_ << "LSGrid::dc_pf: Size of the Vinit should be the same as the total number of buses. Currently:  ";
        exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_bus << " buses.";
        exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
        throw std::runtime_error(exc_.str());
    }
    // DC ignores max_iter / tol, but nonsensical values still indicate a bug
    // at the call site: reject them the same way ac_pf does
    BaseAlgo::check_iter_tol("LSGrid::dc_pf", max_iter, tol);
    bool conv = false;
    CplxVect res = CplxVect();

    // reset_results();  // clear the results  No need to do it, results are neceassirly set or reset in post process

    // pre process the data: builds the real DC admittance matrix Bbus_dc_ and real power vector dcPbus_
    bool is_ac = false;
    CplxVect V = pre_process_dc_solver(Vinit,
                                       dcPbus_,
                                       Bbus_dc_,
                                       id_me_to_dc_solver_,
                                       id_dc_solver_to_me_,
                                       slack_bus_id_dc_me_,
                                       slack_bus_id_dc_solver_,
                                       algo_controler_.dc_algo_controler());
    // start the solver (native real DC entry point)
    conv = _dc_algo.compute_pf_dc(
        Bbus_dc_,
        V,
        dcPbus_,
        slack_bus_id_dc_solver_.as_eigen(),  // was _to_intvect()
        slack_weights_,
        bus_pv_.as_eigen(),  // was _to_intvect()
        bus_pq_.as_eigen());  // was _to_intvect()
    // store results (fase -> because I am in dc mode)
    process_results(conv, res, Vinit, is_ac, id_me_to_dc_solver_);
    timer_last_dc_pf_ = timer.duration();
    return res;
}

RealMat LSGrid::get_ptdf_solver(){
    if(Bbus_dc_.size() == 0){
        throw std::runtime_error("LSGrid::get_ptdf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    // return the freshly-computed matrix directly (RVO/move) instead of binding it
    // to a local const-ref first: `const RealMat& x = ...; return x;` defeats RVO
    // and forces an extra full-matrix copy, since a reference can't be moved from.
    return _dc_algo.get_ptdf();
}


RealMat LSGrid::get_ptdf(){
    if(Bbus_dc_.size() == 0){
        throw std::runtime_error("LSGrid::get_ptdf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    const RealMat & PTDF_solver = get_ptdf_solver();
    RealMat PTDF_grid =  RealMat::Zero(powerlines_.nb() + trafos_.nb(), total_bus());  // , std::numeric_limits<real_type>::quiet_NaN()
    int solver_col = 0;
    for(const GlobalBusId & my_col: id_dc_solver_to_me()){
        PTDF_grid.col(my_col.cast_int()) = PTDF_solver.col(solver_col);
        ++solver_col;
    }
    return PTDF_grid;
}

RealMat LSGrid::get_lodf(){
    if(Bbus_dc_.size() == 0){
        throw std::runtime_error("LSGrid::get_lodf: Cannot get the ptdf without having first computed a DC powerflow.");
    }
    const size_t n_line = powerlines_.nb();
    const size_t nb_el = n_line + trafos_.nb();
    // retrieve the from_bus / to_bus from the grid
    GlobalBusIdVect from_bus = GlobalBusIdVect::concat(powerlines_.get_bus_id_side_1(), trafos_.get_bus_id_side_1());
    GlobalBusIdVect to_bus   = GlobalBusIdVect::concat(powerlines_.get_bus_id_side_2(), trafos_.get_bus_id_side_2());
    const auto & status1_line = powerlines_.get_status_side_1();
    const auto & status2_line = powerlines_.get_status_side_2();
    const auto & status1_trafo = trafos_.get_status_side_1();
    const auto & status2_trafo = trafos_.get_status_side_2();

    // convert it to solver bus id
    IntVect from_bus_solver(nb_el);  // TODO : SolverBusIdVect here
    IntVect to_bus_solver(nb_el);
    for(size_t el_id = 0; el_id < nb_el; ++el_id){
        const bool is_dc_connected = el_id < n_line
            ? (status1_line[el_id] && status2_line[el_id])
            : (status1_trafo[el_id - n_line] && status2_trafo[el_id - n_line]);
        if(!is_dc_connected){
            // half-open (see keep_half_open_lines) or fully disconnected: this
            // branch carries no DC flow at all (TwoSidesContainer_rxh_A::fillBdc
            // drops it from Bbus entirely -- "disco on one side == disco on both
            // sides"), and its open/stale bus id must not index
            // id_me_to_dc_solver_ -- propagate the deactivated sentinel instead,
            // so BaseDCAlgo::get_lodf gives it the identity treatment (its
            // "outage" changes nothing, anywhere).
            from_bus_solver[el_id] = BaseConstants::_deactivated_bus_id;
            to_bus_solver[el_id] = BaseConstants::_deactivated_bus_id;
            continue;
        }
        // from side
        GlobalBusId f_grid_bus = from_bus[el_id];
        SolverBusId f_solver_bus = id_me_to_dc_solver_[f_grid_bus.cast_int()];
        from_bus_solver[el_id] = f_solver_bus.cast_int();
        // to side
        GlobalBusId t_grid_bus = to_bus[el_id];
        SolverBusId t_solver_bus = id_me_to_dc_solver_[t_grid_bus.cast_int()];
        to_bus_solver[el_id] = t_solver_bus.cast_int();
    }
    return _dc_algo.get_lodf(from_bus_solver, to_bus_solver);
}

Eigen::SparseMatrix<real_type> LSGrid::get_Bf_solver(){
    if(Bbus_dc_.size() == 0){
        throw std::runtime_error("LSGrid::get_Bf_solver: Cannot get the Bf matrix without having first computed a DC powerflow.");
    }
    Eigen::SparseMatrix<real_type> Bf;
    fillBf_for_PTDF(Bf);
    return Bf;
}

Eigen::SparseMatrix<real_type> LSGrid::get_Bf(){
    if(Bbus_dc_.size() == 0){
        throw std::runtime_error("LSGrid::get_Bf: Cannot get the Bf matrix without having first computed a DC powerflow.");
    }
    Eigen::SparseMatrix<real_type> Bf_solver = get_Bf_solver();
    return _relabel_matrix(Bf_solver, id_dc_solver_to_me_, false);
}

void LSGrid::add_gen_slackbus(int gen_id, real_type weight){
    if(gen_id < 0)
    {
        std::ostringstream exc_;
        exc_ << "LSGrid::add_gen_slackbus: Slack bus should be an id of a generator, thus positive. You provided: ";
        exc_ << gen_id;
        throw std::runtime_error(exc_.str());
    }
    if(gen_id >= generators_.nb())
    {
        std::ostringstream exc_;
        exc_ << "LSGrid::add_gen_slackbus: There are only " << generators_.nb() << " generators on the grid. ";
        exc_ << "Generator with id " << gen_id << " does not exist and can't be the slack bus";
        throw std::runtime_error(exc_.str());
    }
    if(weight <= 0.){
        std::ostringstream exc_;
        exc_ << "LSGrid::add_gen_slackbus: please enter a valid weight for the slack bus (> 0.)";
        throw std::runtime_error(exc_.str());
    }
    generators_.add_slackbus(gen_id, weight, algo_controler_);
}

void LSGrid::remove_gen_slackbus(int gen_id){
    if(gen_id < 0)
    {
        // TODO DEBUG MODE: only check when in debug mode
        std::ostringstream exc_;
        exc_ << "LSGrid::remove_gen_slackbus: Slack bus should be an id of a generator, thus positive. You provided: ";
        exc_ << gen_id;
        throw std::runtime_error(exc_.str());
    }
    if(gen_id >= generators_.nb())
    {
        // TODO DEBUG MODE: only check when in debug mode
        std::ostringstream exc_;
        exc_ << "LSGrid::remove_gen_slackbus: There are only " << generators_.nb() << " generators on the grid. ";
        exc_ << "Generator with id " << gen_id << " does not exist and can't be the slack bus";
        throw std::runtime_error(exc_.str());
    }
    generators_.remove_slackbus(gen_id, algo_controler_);
}

/** GRID2OP SPECIFIC REPRESENTATION **/
void LSGrid::update_gens_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_gen);
}

void LSGrid::update_sgens_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_sgen);
}

void LSGrid::update_gens_v(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_v_gen);
}

void LSGrid::update_loads_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_load);
}

void LSGrid::update_loads_q(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_q_load);
}

void LSGrid::update_storages_p(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                              const Eigen::Ref<const Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_storage);
}

void LSGrid::update_topo(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                            const Eigen::Ref<const Eigen::Array<int,  Eigen::Dynamic, Eigen::RowMajor> > & new_values)
{
    // Both arrays come straight from python and are indexed BY POSITION IN THE
    // TOPOLOGY VECTOR: each container does `has_changed(pos_topo_vect_(el_id))` /
    // `new_values(pos_topo_vect_(el_id))` with an unchecked Eigen operator(). The
    // positions themselves are validated (check_grid() proves they form a
    // permutation of [0, dim_topo)), but nothing checked that the caller's arrays
    // are dim_topo long -- a shorter one reads past its end. dim_topo is exactly
    // the number of topology-participating element sides, so compute it and demand
    // both arrays match it.
    const Eigen::Index dim_topo =
        static_cast<Eigen::Index>(loads_.nb()) +
        static_cast<Eigen::Index>(generators_.nb()) +
        static_cast<Eigen::Index>(storages_.nb()) +
        2 * static_cast<Eigen::Index>(powerlines_.nb()) +
        2 * static_cast<Eigen::Index>(trafos_.nb());
    if((has_changed.rows() != dim_topo) || (new_values.rows() != dim_topo)){
        std::ostringstream exc_;
        exc_ << "LSGrid::update_topo: 'has_changed' (size " << has_changed.rows()
             << ") and 'new_values' (size " << new_values.rows() << ") must both have the size of "
             << "the topology vector (" << dim_topo << " = nb loads + nb gens + nb storages + "
             << "2 * nb lines + 2 * nb trafos). They are indexed by position in the topology "
             << "vector, so a shorter array would be read out of bounds.";
        throw std::runtime_error(exc_.str());
    }
    loads_.update_topo(has_changed, new_values, algo_controler_, substations_);
    generators_.update_topo(has_changed, new_values, algo_controler_, substations_);
    storages_.update_topo(has_changed, new_values, algo_controler_, substations_);
    // shunts are not in "topo" in grid2op

    // NB we suppose that if a powerline (or a trafo) is disconnected, then both its ends are
    // and same for trafo, obviously
    powerlines_.update_topo(has_changed, new_values, algo_controler_, substations_);
    trafos_.update_topo(has_changed, new_values, algo_controler_, substations_);
}

// for FDPF (implementation of the alg 2 method FDBX (FDXB will follow)  // TODO FDPF
void LSGrid::fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, 
                           Eigen::SparseMatrix<real_type> & Bpp, 
                           FDPFMethod xb_or_bx) const
{
    // clear the matrices
    const int nb_bus_solver = static_cast<int>(id_ac_solver_to_me_.size());
    Bp = Eigen::SparseMatrix<real_type>(nb_bus_solver, nb_bus_solver);
    Bpp = Eigen::SparseMatrix<real_type>(nb_bus_solver, nb_bus_solver);

    // init the Bp and Bpp matrices for Fast Decoupled Powerflow  (TODO FDPF: optim when it's NOT needed just like for Ybus)
    std::vector<Eigen::Triplet<real_type> > tripletList_Bp;
    std::vector<Eigen::Triplet<real_type> > tripletList_Bpp;
    tripletList_Bp.reserve(substations_.nb_bus() + 4 * powerlines_.nb() + 4 * trafos_.nb() + shunts_.nb());
    tripletList_Bpp.reserve(substations_.nb_bus() + 4 * powerlines_.nb() + 4 * trafos_.nb() + shunts_.nb());
    // run through the grid and get the parameters to fill them
    powerlines_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    trafos_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    loads_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    sgens_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    storages_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    generators_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    hvdc_lines_.fillBp_Bpp(tripletList_Bp, tripletList_Bpp, id_me_to_ac_solver_, sn_mva_, xb_or_bx);
    // now make the matrices effectively
    Bp.setFromTriplets(tripletList_Bp.begin(), tripletList_Bp.end());
    Bp.makeCompressed();
    Bpp.setFromTriplets(tripletList_Bpp.begin(), tripletList_Bpp.end());
    Bpp.makeCompressed();
}


void LSGrid::fillBf_for_PTDF(Eigen::SparseMatrix<real_type> & Bf, bool transpose) const
{
    const int nb_bus_solver = static_cast<int>(id_dc_solver_to_me_.size());
    // TODO DEBUG MODE
    if(nb_bus_solver == 0) throw std::runtime_error("LSGrid::fillBf_for_PTDF: it appears no DC powerflow has run on your grid.");
    
    if(transpose){
        Bf = Eigen::SparseMatrix<real_type>(nb_bus_solver, powerlines_.nb() + trafos_.nb());
    }else{
        Bf = Eigen::SparseMatrix<real_type>(powerlines_.nb() + trafos_.nb(), nb_bus_solver);
    }
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(substations_.nb_bus() + 2 * powerlines_.nb() + 2 * trafos_.nb());
    
    powerlines_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);  // TODO have a function to dispatch that to all type of elements
    shunts_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    trafos_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    loads_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    sgens_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    storages_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    generators_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);
    hvdc_lines_.fillBf_for_PTDF(tripletList, id_me_to_dc_solver_, sn_mva_, powerlines_.nb(), transpose);

    Bf.setFromTriplets(tripletList.begin(), tripletList.end());
    Bf.makeCompressed();
}

// returns only the gen_id with the highest p that is connected to this bus !
// returns bus_id, gen_bus_id
std::tuple<int, int> LSGrid::assign_slack_to_most_connected(){
    auto res = std::tuple<int, int>(-1, -1);
    int res_bus_id = -1;
    int res_gen_id = -1;
    int max_line = -1;
    const unsigned int nb_busbars = substations_.nb_bus();
    std::vector<real_type> gen_p_per_bus(nb_busbars, 0.);
    std::vector<int> nb_line_end_per_bus(nb_busbars, 0);

    // computes the total amount of power produce at each nodes
    powerlines_.gen_p_per_bus(gen_p_per_bus);  // TODO have a function to dispatch that to all type of elements
    shunts_.gen_p_per_bus(gen_p_per_bus);
    trafos_.gen_p_per_bus(gen_p_per_bus);
    loads_.gen_p_per_bus(gen_p_per_bus);
    sgens_.gen_p_per_bus(gen_p_per_bus);
    storages_.gen_p_per_bus(gen_p_per_bus);
    generators_.gen_p_per_bus(gen_p_per_bus);
    hvdc_lines_.gen_p_per_bus(gen_p_per_bus);

    // computes the total number of "neighbors" (extremity of connected powerlines and trafo, not real neighbors)
    powerlines_.nb_line_end(nb_line_end_per_bus);  // TODO have a function to dispatch that to all type of elements
    shunts_.nb_line_end(nb_line_end_per_bus);
    trafos_.nb_line_end(nb_line_end_per_bus);
    loads_.nb_line_end(nb_line_end_per_bus);
    sgens_.nb_line_end(nb_line_end_per_bus);
    storages_.nb_line_end(nb_line_end_per_bus);
    generators_.nb_line_end(nb_line_end_per_bus);
    hvdc_lines_.nb_line_end(nb_line_end_per_bus);
    
    // now find the most connected buses
    for(unsigned int bus_id = 0; bus_id < nb_busbars; ++bus_id)
    {
        const auto & nb_lines_this = nb_line_end_per_bus[bus_id];
        if((nb_lines_this > max_line) && (gen_p_per_bus[bus_id] > 0.)){
            res_bus_id = bus_id;
            max_line = nb_lines_this;
        }
    }
    // TODO DEBUG MODE
    if(res_bus_id == -1) throw std::runtime_error("LSGrid::assign_slack_to_most_connected: impossible to find anything connected to a node.");
    std::get<0>(res) = res_bus_id;

    // and reset the slack bus
    generators_.remove_all_slackbus();
    res_gen_id = generators_.assign_slack_bus(res_bus_id, gen_p_per_bus, algo_controler_);
    std::get<1>(res) = res_gen_id;
    slack_bus_id_ac_solver_ = SolverBusIdVect();
    slack_bus_id_dc_solver_ = SolverBusIdVect();
    slack_weights_ = RealVect();
    return res;
}

// TODO DC LINE: one side might be in the connected comp and not the other !
void LSGrid::consider_only_main_component(){
    const auto & slack_buses_id = generators_.get_slack_bus_id();

    // TODO DEBUG MODE
    if(slack_buses_id.size() == 0) throw std::runtime_error("LSGrid::consider_only_main_component: no slack is defined on your grid. This function cannot be used.");
    
    // build the graph
    const auto nb_busbars = substations_.nb_bus();
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(2 * powerlines_.nb() + 2 * trafos_.nb());
    powerlines_.get_graph(tripletList);  // TODO have a function to dispatch that to all type of elements
    shunts_.get_graph(tripletList);
    trafos_.get_graph(tripletList);
    loads_.get_graph(tripletList);
    sgens_.get_graph(tripletList);
    storages_.get_graph(tripletList);
    generators_.get_graph(tripletList);
    hvdc_lines_.get_graph(tripletList);
    Eigen::SparseMatrix<real_type> graph = Eigen::SparseMatrix<real_type>(nb_busbars, nb_busbars);
    graph.setFromTriplets(tripletList.begin(), tripletList.end());
    graph.makeCompressed();

    // find the connected buses
    // TODO copy paste from SecurityAnalysis
    std::vector<bool> tmp_visited(nb_busbars, false);
    std::vector<int> conn_comp(nb_busbars, -1);
    std::vector<bool> already_added(nb_busbars, false);

    int connected_comp = 0;
    std::queue<GlobalBusId> neighborhood;
    while(true)
    {
        neighborhood = std::queue<GlobalBusId>();

        // choose bus id (one of the slack) to start
        bool one_added = false;
        for(const auto & el : slack_buses_id){
            if(!tmp_visited[el.cast_int()] && !already_added[el.cast_int()])
            {
                one_added = true;
                neighborhood.push(el);
                already_added[el.cast_int()] = true;
                break;
            }
        }
        
        if(!one_added) break; // no more slack bus, I stop

        // start the bfs
        while (true)
        {
            const GlobalBusId col_id = neighborhood.front();
            neighborhood.pop();
            tmp_visited[col_id.cast_int()] = true;
            conn_comp[col_id.cast_int()] = connected_comp;
            for (Eigen::SparseMatrix<real_type>::InnerIterator it(graph, col_id.cast_int()); it; ++it)
            {
                // add in the queue all my neighbor
                if(!tmp_visited[it.row()] && !already_added[it.row()]){
                    neighborhood.push(GlobalBusId(it.row()));
                    already_added[it.row()] = true;
                }
            }
            if(neighborhood.empty()) break;  // no more neighbors
        }

        // go to the next connected comp
        ++connected_comp;
    }

    // TODO speed optim: if connected_comp == 1 => don't do the following 2 steps
    
    // find the connected comp with the most buses
    int main_cc_id = -1;
    std::vector<int> nb_bus_per_cc(connected_comp, 0);
    for(const auto el : conn_comp){
        if(el == -1) continue;
        nb_bus_per_cc[el] += 1;
    }
    main_cc_id = std::distance(nb_bus_per_cc.begin(),
                               std::max_element(nb_bus_per_cc.begin(), nb_bus_per_cc.end()));

    // mark as visited the element in this cc
    std::vector<bool> bus_in_main_cc(nb_busbars, false);
    for(size_t bus_id = 0; bus_id < nb_busbars; ++bus_id){
        if(conn_comp[bus_id] == main_cc_id) bus_in_main_cc[bus_id] = true;
    }
    // disconnected elements not in main component
    powerlines_.disconnect_if_not_in_main_component(bus_in_main_cc);
    shunts_.disconnect_if_not_in_main_component(bus_in_main_cc);
    trafos_.disconnect_if_not_in_main_component(bus_in_main_cc);
    loads_.disconnect_if_not_in_main_component(bus_in_main_cc);
    sgens_.disconnect_if_not_in_main_component(bus_in_main_cc);
    storages_.disconnect_if_not_in_main_component(bus_in_main_cc);
    generators_.disconnect_if_not_in_main_component(bus_in_main_cc);
    hvdc_lines_.disconnect_if_not_in_main_component(bus_in_main_cc);
    // and finally deal with the buses
    init_bus_status();
}

} // namespace ls2g
