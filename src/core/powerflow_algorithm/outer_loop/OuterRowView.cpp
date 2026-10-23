// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "OuterRowView.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "LSGrid.hpp"
#include "BaseOuterLoop.hpp"

namespace ls2g {

// out-of-class definition: the C++14 build odr-uses it (vector::assign takes a const reference)
constexpr int OuterRowView::Units::GRID_BUS;

namespace {

void overlay_targets(std::vector<real_type> & row, const std::vector<real_type> & moved)
{
    for (std::size_t k = 0; k < moved.size(); ++k) {
        if (std::isnan(moved[k])) continue;
        if (row.size() < moved.size()) row.resize(moved.size(), std::numeric_limits<real_type>::quiet_NaN());
        row[k] = moved[k];
    }
}

}  // namespace

void OuterRowView::take_moved_targets(const std::vector<real_type> & gen_target_p,
                                      const std::vector<real_type> & storage_target_p)
{
    overlay_targets(gens.target_p, gen_target_p);
    overlay_targets(storages.target_p, storage_target_p);
}

const SolverBusLayout & OuterGridView::layout() const
{
    return (row_ != nullptr && row_->layout != nullptr) ? *row_->layout : grid_.get_ac_layout();
}

const SolverBusIdVect & OuterGridView::id_me_to_solver() const { return layout().id_me_to_solver; }

const GlobalBusIdVect & OuterGridView::id_solver_to_me() const { return layout().id_solver_to_me; }

bool OuterGridView::is_masked(int bus_solver) const
{
    const std::vector<int> * m = masked();
    return m != nullptr && std::binary_search(m->begin(), m->end(), bus_solver);
}

bool OuterGridView::_value(const std::vector<real_type> & v, int id, real_type & res)
{
    if (id < 0 || static_cast<std::size_t>(id) >= v.size()) return false;
    const real_type x = v[static_cast<std::size_t>(id)];
    if (std::isnan(x)) return false;
    res = x;
    return true;
}

namespace {

// the solver bus of a grid bus, -1 when it is not solved
int to_solver(const SolverBusIdVect & id_me_to_solver, int bus_me)
{
    if (bus_me < 0 || bus_me >= static_cast<int>(id_me_to_solver.size())) return -1;
    const int b = id_me_to_solver[bus_me].cast_int();
    return b == BaseConstants::_deactivated_bus_id ? -1 : b;
}

// whether an element of `container` is connected in a row with these overrides
template<class Container>
bool unit_on(const Container & container, const OuterRowView::Units * u, int id)
{
    if (id < 0 || id >= container.nb()) return false;
    if (u != nullptr) {
        if (static_cast<std::size_t>(id) < u->off.size() && u->off[static_cast<std::size_t>(id)] != 0) return false;
        // a placed element is on, whatever the grid says
        if (static_cast<std::size_t>(id) < u->bus_solver.size() &&
            u->bus_solver[static_cast<std::size_t>(id)] != OuterRowView::Units::GRID_BUS) return true;
    }
    return container.get_status()[id];
}

template<class Container>
int unit_bus_solver(const Container & container, const OuterRowView::Units * u, int id,
                    const SolverBusIdVect & id_me_to_solver)
{
    if (!unit_on(container, u, id)) return -1;
    if (u != nullptr && static_cast<std::size_t>(id) < u->bus_solver.size()) {
        const int placed = u->bus_solver[static_cast<std::size_t>(id)];
        if (placed != OuterRowView::Units::GRID_BUS) return placed;
    }
    return to_solver(id_me_to_solver, container.get_bus_id()(id).cast_int());
}

}  // namespace

// ----- generators ----------------------------------------------------------------------

bool OuterGridView::gen_on(int gen_id) const
{
    return unit_on(grid_.get_generators(), row_ != nullptr ? &row_->gens : nullptr, gen_id);
}

int OuterGridView::gen_bus_solver(int gen_id) const
{
    return unit_bus_solver(grid_.get_generators(), row_ != nullptr ? &row_->gens : nullptr, gen_id, id_me_to_solver());
}

int OuterGridView::gen_bus_grid(int gen_id) const
{
    if (!gen_on(gen_id)) return -1;
    const int placed = _placed(row_ != nullptr ? &row_->gens : nullptr, gen_id);
    if (placed == OuterRowView::Units::GRID_BUS) return grid_.get_generators().get_bus_id()(gen_id).cast_int();
    const GlobalBusIdVect & to_me = id_solver_to_me();
    return (placed >= 0 && placed < static_cast<int>(to_me.size())) ? to_me[placed].cast_int() : -1;
}

bool OuterGridView::gen_is_local_voltage_controller(int gen_id) const
{
    return gen_on(gen_id) && grid_.get_generators().would_be_local_voltage_controller(gen_id);
}

real_type OuterGridView::gen_target_p(int gen_id) const
{
    real_type res;
    if (row_ != nullptr && _value(row_->gens.target_p, gen_id, res)) return res;
    return grid_.get_generators().get_target_p()(gen_id);
}

bool OuterGridView::row_sets_gen_target_p(int gen_id) const
{
    real_type res;
    return row_ != nullptr && _value(row_->gens.target_p, gen_id, res);
}

real_type OuterGridView::gen_target_vm(int gen_id) const
{
    real_type res;
    if (row_ != nullptr && _value(row_->gens.target_vm, gen_id, res)) return res;
    return grid_.get_generators().get_target_vm_pu(gen_id);
}

real_type OuterGridView::gen_target_p_override(const std::vector<real_type> & loops, int gen_id) const
{
    real_type res;
    if (_value(loops, gen_id, res)) return res;
    if (row_ != nullptr && _value(row_->gens.target_p, gen_id, res)) return res;
    return std::numeric_limits<real_type>::quiet_NaN();
}

// ----- storage units -------------------------------------------------------------------

bool OuterGridView::storage_on(int storage_id) const
{
    return unit_on(grid_.get_storages(), row_ != nullptr ? &row_->storages : nullptr, storage_id);
}

int OuterGridView::storage_bus_solver(int storage_id) const
{
    return unit_bus_solver(grid_.get_storages(), row_ != nullptr ? &row_->storages : nullptr, storage_id, id_me_to_solver());
}

bool OuterGridView::storage_is_local_voltage_controller(int storage_id) const
{
    return storage_on(storage_id) && grid_.get_storages().would_be_local_voltage_controller(storage_id);
}

real_type OuterGridView::storage_target_p(int storage_id) const
{
    real_type res;
    if (row_ != nullptr && _value(row_->storages.target_p, storage_id, res)) return res;
    return grid_.get_storages().get_target_p()(storage_id);
}

bool OuterGridView::row_sets_storage_target_p(int storage_id) const
{
    real_type res;
    return row_ != nullptr && _value(row_->storages.target_p, storage_id, res);
}

real_type OuterGridView::storage_target_vm(int storage_id) const
{
    real_type res;
    if (row_ != nullptr && _value(row_->storages.target_vm, storage_id, res)) return res;
    return grid_.get_storages().get_target_vm_pu(storage_id);
}

// ----- loads ---------------------------------------------------------------------------

bool OuterGridView::load_on(int load_id) const
{
    return unit_on(grid_.get_loads(), row_ != nullptr ? &row_->loads : nullptr, load_id);
}

int OuterGridView::load_bus_solver(int load_id) const
{
    return unit_bus_solver(grid_.get_loads(), row_ != nullptr ? &row_->loads : nullptr, load_id, id_me_to_solver());
}

real_type OuterGridView::load_target_q(int load_id) const
{
    real_type res;
    if (row_ != nullptr && _value(row_->loads.target_q, load_id, res)) return res;
    return grid_.get_loads().get_target_q()(load_id);
}

// ----- branches ------------------------------------------------------------------------

namespace {

// one side of a branch in a row: -1 when that side is open or its bus is not solved
template<class Container>
int branch_side(const Container & branches, const std::vector<int> & row_bus, int el, int side,
                const SolverBusIdVect & id_me_to_solver)
{
    if (static_cast<std::size_t>(el) < row_bus.size()) {
        const int moved = row_bus[static_cast<std::size_t>(el)];
        if (moved != OuterRowView::Units::GRID_BUS) return moved < 0 ? -1 : moved;
    }
    const bool on = side == 1 ? branches.get_status_side_1()[el] : branches.get_status_side_2()[el];
    if (!on) return -1;
    const int bus_me = (side == 1 ? branches.get_bus_side_1(el) : branches.get_bus_side_2(el)).cast_int();
    return to_solver(id_me_to_solver, bus_me);
}

template<class Container>
bool branch_closed(const Container & branches, const OuterRowView::Branches * b, int el,
                   const SolverBusIdVect & id_me_to_solver, int & bus1, int & bus2)
{
    bus1 = bus2 = -1;
    if (el < 0 || el >= branches.nb()) return false;
    static const std::vector<int> none;
    if (b != nullptr && static_cast<std::size_t>(el) < b->off.size() && b->off[static_cast<std::size_t>(el)] != 0) return false;
    // a row that moves an end reconnects the branch there, whatever the grid's global status
    const bool moved = b != nullptr &&
        ((static_cast<std::size_t>(el) < b->bus1_solver.size() && b->bus1_solver[static_cast<std::size_t>(el)] != OuterRowView::Units::GRID_BUS) ||
         (static_cast<std::size_t>(el) < b->bus2_solver.size() && b->bus2_solver[static_cast<std::size_t>(el)] != OuterRowView::Units::GRID_BUS));
    if (!moved && !branches.get_status_global()[el]) return false;
    bus1 = branch_side(branches, b != nullptr ? b->bus1_solver : none, el, 1, id_me_to_solver);
    bus2 = branch_side(branches, b != nullptr ? b->bus2_solver : none, el, 2, id_me_to_solver);
    return bus1 >= 0 && bus2 >= 0;
}

}  // namespace

bool OuterGridView::line_closed(int line_id, int & bus1_solver, int & bus2_solver) const
{
    return branch_closed(grid_.get_lines(), row_ != nullptr ? &row_->lines : nullptr, line_id,
                         id_me_to_solver(), bus1_solver, bus2_solver);
}

bool OuterGridView::trafo_closed(int trafo_id, int & bus1_solver, int & bus2_solver) const
{
    return branch_closed(grid_.get_trafos(), row_ != nullptr ? &row_->trafos : nullptr, trafo_id,
                         id_me_to_solver(), bus1_solver, bus2_solver);
}

// ----- OuterInjections ----------------------------------------------------------------

real_type OuterInjections::gen_target_p_in_force(const OuterGridView & view, int gen_id) const
{
    if (gen_id >= 0 && static_cast<std::size_t>(gen_id) < gen_target_p.size() &&
        !std::isnan(gen_target_p[static_cast<std::size_t>(gen_id)])) {
        return gen_target_p[static_cast<std::size_t>(gen_id)];
    }
    return view.gen_target_p(gen_id);
}

real_type OuterInjections::storage_target_p_in_force(const OuterGridView & view, int storage_id) const
{
    if (storage_id >= 0 && static_cast<std::size_t>(storage_id) < storage_target_p.size() &&
        !std::isnan(storage_target_p[static_cast<std::size_t>(storage_id)])) {
        return storage_target_p[static_cast<std::size_t>(storage_id)];
    }
    return view.storage_target_p(storage_id);
}

void OuterInjections::_add_injection(int bus_solver, const cplx_type & ds)
{
    if (ds == cplx_type(0., 0.)) return;
    (*Sbus)(bus_solver) += ds;
    if (Sbus_target != nullptr) (*Sbus_target)(bus_solver) += ds;
}

real_type OuterInjections::move_gen_target_p(const OuterGridView & view, int gen_id, int bus_solver, real_type mw)
{
    const LSGrid & grid = view.grid();
    const real_type sn_mva = grid.get_sn_mva();
    if (gen_target_p.empty()) {
        gen_target_p.assign(static_cast<std::size_t>(grid.get_generators().nb()), std::numeric_limits<real_type>::quiet_NaN());
    }
    const real_type delta = mw - gen_target_p_in_force(view, gen_id);
    cplx_type ds = {delta / sn_mva, 0.};
    const bool pq = !grid.get_generators().get_voltage_regulator_on(gen_id);
    const real_type q_before = pq ? grid.gen_target_q_at_outer_target(
        gen_id, view.gen_target_p_override(gen_target_p, gen_id)) : 0.;
    gen_target_p[static_cast<std::size_t>(gen_id)] = mw;
    if (pq) ds += cplx_type(0., (grid.gen_target_q_at_outer_target(
        gen_id, view.gen_target_p_override(gen_target_p, gen_id)) - q_before) / sn_mva);
    _add_injection(bus_solver, ds);
    return delta;
}

real_type OuterInjections::move_storage_target_p(const OuterGridView & view, int storage_id, int bus_solver, real_type mw)
{
    const LSGrid & grid = view.grid();
    if (storage_target_p.empty()) {
        storage_target_p.assign(static_cast<std::size_t>(grid.get_storages().nb()), std::numeric_limits<real_type>::quiet_NaN());
    }
    // load convention: the injection moves the other way
    const real_type delta = storage_target_p_in_force(view, storage_id) - mw;
    storage_target_p[static_cast<std::size_t>(storage_id)] = mw;
    _add_injection(bus_solver, cplx_type(delta / grid.get_sn_mva(), 0.));
    return delta;
}

}  // namespace ls2g
