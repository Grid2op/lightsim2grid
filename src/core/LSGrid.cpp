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

#include <queue>

namespace ls2g {


LSGrid::LSGrid(const LSGrid & other)
{
    init_vm_pu_ = other.init_vm_pu_;
    sn_mva_ = other.sn_mva_;
    compute_results_ = other.compute_results_;

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

    // assign the right solver
    reset(true, true, true);
    _algo.change_algorithm(other.get_algo_type());
    _dc_algo.change_algorithm(other.get_dc_algo_type());
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
                            get_algo_type(),
                            get_dc_algo_type(),
                            res_svc
                            );
    return res;
};

void LSGrid::set_state(LSGrid::StateRes & my_state)
{
    // after loading back, the instance need to be reset anyway
    // TODO see if it's worth the trouble NOT to do it
    algo_controler_.ac_algo_controler().tell_all_changed();
    algo_controler_.dc_algo_controler().tell_all_changed();
    compute_results_ = true;

    // extract data from the state
    std::string version_major = std::get<0>(my_state);
    std::string version_medium = std::get<1>(my_state);
    std::string version_minor = std::get<2>(my_state);
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
    const std::vector<int> & ls_to_pp = std::get<3>(my_state);
    init_vm_pu_ = std::get<4>(my_state);
    sn_mva_ = std::get<5>(my_state);
    // const std::vector<real_type> & bus_vn_kv = std::get<6>(my_state);
    const std::vector<bool> & last_bus_status_saved = std::get<6>(my_state);
    SubstationContainer::StateRes & state_substations = std::get<7>(my_state);
    // powerlines
    LineContainer::StateRes & state_lines = std::get<8>(my_state);
    // shunts
    ShuntContainer::StateRes & state_shunts = std::get<9>(my_state);
    // trafos
    TrafoContainer::StateRes & state_trafos = std::get<10>(my_state);
    // generators
    // total_q_min_per_bus_;
    // total_q_max_per_bus_;
    // total_gen_per_bus_;
    GeneratorContainer::StateRes & state_gens = std::get<11>(my_state);
    // loads
    LoadContainer::StateRes & state_loads = std::get<12>(my_state);
    // static gen
    SGenContainer::StateRes & state_sgens= std::get<13>(my_state);
    // storage units
    StorageContainer::StateRes & state_storages = std::get<14>(my_state);
    // hvdc lines
    HvdcLineContainer::StateRes & state_hvdc_lines = std::get<15>(my_state);
    // static var compensators (index 16/17 are the ac/dc algo types)
    SvcContainer::StateRes & state_svcs = std::get<18>(my_state);

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
    _algo.change_algorithm(std::get<16>(my_state));
    _dc_algo.change_algorithm(std::get<17>(my_state));
};

void LSGrid::save_binary(const std::string & path) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

LSGrid LSGrid::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<LSGrid>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

void LSGrid::set_ls_to_orig(const IntVect & ls_to_orig){
    if(ls_to_orig.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }

    if(static_cast<size_t>(ls_to_orig.size()) != substations_.nb_bus())
        throw std::runtime_error("Impossible to set the converter ls_to_orig: the provided vector has not the same size as the number of bus on the grid.");
    set_ls_to_orig_internal(ls_to_orig);
}

void LSGrid::set_orig_to_ls(const IntVect & orig_to_ls){
    if(orig_to_ls.size() == 0){
        _ls_to_orig = IntVect();
        _orig_to_ls = IntVect();
        return;
    }
    _orig_to_ls = orig_to_ls;
    size_t nb_bus_ls = 0;
    for(const auto el : orig_to_ls){
        if (el != -1) nb_bus_ls += 1;
    }
    if(nb_bus_ls != substations_.nb_bus()) 
        throw std::runtime_error("Impossible to set the converter orig_to_ls: the number of 'non -1' component in the provided vector does not match the number of buses on the grid.");
    _ls_to_orig = IntVect::Constant(nb_bus_ls, -1);
    size_t ls2or_ind = 0;
    for(size_t or2ls_ind = 0; or2ls_ind < nb_bus_ls; ++or2ls_ind){
        const auto my_ind = _orig_to_ls[or2ls_ind];
        if(my_ind >= 0){
            _ls_to_orig[ls2or_ind] = my_ind;
            ls2or_ind++;
        }
    }
}

void LSGrid::set_ls_to_orig_internal(const IntVect & ls_to_orig) noexcept{
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
                         const RealVect & bus_vn_kv,
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

CplxVect LSGrid::ac_pf(const CplxVect & Vinit,
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
    if(hvdc_lines_.has_droop_active() && !_algo.supports_hvdc_droop(_algo.get_type())){
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
    for(int pos = 0; pos < nb_droop; ++pos){
        const int hvdc_id = indices[pos];
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
    // such a slack bus is supported even though `is_pq` is false there.
    std::vector<bool> is_slack(nb_bus_solver, false);
    for(int k = 0; k < static_cast<int>(slack_bus_id_ac_solver_.size()); ++k){
        const int b = slack_bus_id_ac_solver_(k).cast_int();
        if(b >= 0 && b < nb_bus_solver) is_slack[b] = true;
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
        if(!is_pq[ctrl_solver] && !is_slack[ctrl_solver]){
            std::ostringstream exc_;
            exc_ << "LSGrid::fill_voltage_control_solver_data: generator " << gen_id
                 << " regulates a remote bus but its OWN bus has no reactive (Q) equation"
                    " (it is a PV bus that is not a slack). This is not supported in v1.";
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

void LSGrid::check_solution_q_values_onegen(CplxVect & res,
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

void LSGrid::check_solution_q_values(CplxVect & res, bool check_q_limits) const{
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

CplxVect LSGrid::check_solution(const CplxVect & V_proposed, bool check_q_limits)
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
    const CplxVect & Vinit,
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
    const CplxVect & Vinit,
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
    const CplxVect & Vinit,
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

CplxVect LSGrid::_get_results_back_to_orig_nodes(const CplxVect & res_tmp,
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
                                const CplxVect & Vinit,
                                bool ac,
                                SolverBusIdVect & id_me_to_solver)
{
    if (conv){
        if(compute_results_){
            // compute the results of the flows, P,Q,V of loads etc.
            compute_results(ac);
        }
        // solver_control_.tell_none_changed();  // todo automatically set for ac / dc the `tell_none_changed()`
        const CplxVect & res_tmp = ac ? _algo.get_V(): _dc_algo.get_V() ;

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

void LSGrid::fillSbus_me(CplxVect & Sbus, bool ac, const SolverBusIdVect& id_me_to_solver)
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

CplxVect LSGrid::dc_pf(const CplxVect & Vinit,
                          int /*max_iter*/,  // not used for DC
                          real_type /*tol*/  // not used for DC
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
    const RealMat & PTDF_solver = _dc_algo.get_ptdf();
    return PTDF_solver;
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
    const size_t nb_el = powerlines_.nb() + trafos_.nb();
    // retrieve the from_bus / to_bus from the grid
    GlobalBusIdVect from_bus = GlobalBusIdVect::concat(powerlines_.get_bus_id_side_1(), trafos_.get_bus_id_side_1());
    GlobalBusIdVect to_bus   = GlobalBusIdVect::concat(powerlines_.get_bus_id_side_2(), trafos_.get_bus_id_side_2());

    // convert it to solver bus id
    IntVect from_bus_solver(nb_el);  // TODO : SolverBusIdVect here
    IntVect to_bus_solver(nb_el);
    for(size_t el_id = 0; el_id < nb_el; ++el_id){
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
void LSGrid::update_gens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_gen);
}

void LSGrid::update_sgens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_sgen);
}

void LSGrid::update_gens_v(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_v_gen);
}

void LSGrid::update_loads_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_load);
}

void LSGrid::update_loads_q(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_q_load);
}

void LSGrid::update_storages_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                              Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
    update_continuous_values(has_changed, new_values, &LSGrid::change_p_storage);
}

void LSGrid::update_topo(Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<const Eigen::Array<int,  Eigen::Dynamic, Eigen::RowMajor> > new_values)
{
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
