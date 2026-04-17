// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TAGGEDIDVEC_HPP
#define TAGGEDIDVEC_HPP

#include <vector>
#include "Eigen/Core"

const int LOCAL_BUS = 0;
const int GLOBAL_BUS = 1;
const int GRIDMODEL_BUS = GLOBAL_BUS;
const int SOLVER_BUS = 2;


// Forward declaration: IntClass<U> is fully defined in Utils.hpp after this
// header is included.  All template bodies below are only instantiated after
// IntClass is complete, so the forward declaration is sufficient.
template<int U> class IntClass;

// Self-contained IntVect typedef (identical to the one in Utils.hpp; identical
// redefinitions are permitted in C++).
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IntVect;


// ============================================================================
// TaggedIdVec<Tag>  —  replaces  Eigen::Matrix<IntClass<Tag>, Dynamic, 1>
//
// Stores plain int internally.  Element access returns IntClass<Tag> by value
// (const) or through a proxy reference (non-const).
//
// Key benefit: .raw() / .int_vect() expose the underlying IntVect through a
// plain int* pointer — no reinterpret_cast, no UB.
// ============================================================================
// template<int Tag>
// class TaggedIdVec {
// public:
//     using value_type = IntClass<Tag>;
//     using Index      = Eigen::Index;

//     // ----- proxy reference (non-const element assignment) ------------------
//     struct reference {
//         int& v_;
//         explicit reference(int& v) noexcept : v_(v) {}
//         reference& operator=(IntClass<Tag> x) noexcept { v_ = x.cast_int(); return *this; }
//         reference& operator=(int         x) noexcept { v_ = x;             return *this; }
//         operator IntClass<Tag>() const noexcept { return IntClass<Tag>(v_); }
//         explicit operator int() const noexcept  { return v_; }
//         int cast_int()          const noexcept  { return v_; }
//     };

//     // ----- forward const iterator (yields IntClass<Tag> by value) ----------
//     struct const_iterator {
//         const int* ptr_;
//         explicit const_iterator(const int* p) noexcept : ptr_(p) {}
//         IntClass<Tag>   operator*()  const noexcept { return IntClass<Tag>(*ptr_); }
//         const_iterator& operator++()       noexcept { ++ptr_; return *this; }
//         const_iterator  operator++(int)    noexcept { auto t = *this; ++ptr_; return t; }
//         bool operator==(const const_iterator& o) const noexcept { return ptr_ == o.ptr_; }
//         bool operator!=(const const_iterator& o) const noexcept { return ptr_ != o.ptr_; }
//     };

//     // ----- constructors ----------------------------------------------------
//     TaggedIdVec() = default;

//     // Factory analogous to Eigen::Matrix<>::Constant(n, val)
//     static TaggedIdVec Constant(Index n, IntClass<Tag> val) {
//         TaggedIdVec r; r.data_ = IntVect::Constant(n, val.cast_int()); return r;
//     }
//     static TaggedIdVec Constant(Index n, int val) {
//         TaggedIdVec r; r.data_ = IntVect::Constant(n, val); return r;
//     }

//     // ----- size / resize ---------------------------------------------------
//     Index size()  const noexcept { return data_.size(); }
//     bool  empty() const noexcept { return data_.size() == 0; }
//     void  resize(Index n)        { data_.resize(n); }

//     // ----- element access --------------------------------------------------
//     IntClass<Tag> operator()(Index i) const noexcept { return IntClass<Tag>(data_(i)); }
//     reference     operator()(Index i)       noexcept { return reference(data_(i)); }
//     IntClass<Tag> operator[](Index i) const noexcept { return IntClass<Tag>(data_(i)); }
//     reference     operator[](Index i)       noexcept { return reference(data_(i)); }

//     // ----- iteration -------------------------------------------------------
//     const_iterator begin() const noexcept { return const_iterator(data_.data()); }
//     const_iterator end()   const noexcept { return const_iterator(data_.data() + data_.size()); }

//     // ----- UB-free accessors -----------------------------------------------
//     // Returns a Map backed by the int storage — int*→int*, no reinterpret_cast
//     Eigen::Map<const IntVect> raw() const noexcept {
//         return Eigen::Map<const IntVect>(data_.data(), data_.size());
//     }
//     // Zero-copy read access to the underlying IntVect
//     const IntVect& int_vect() const noexcept { return data_; }

// private:
//     IntVect data_;
// };


// ============================================================================
// TaggedIdStdVec<Tag>  —  replaces  std::vector<IntClass<Tag>>
//
// Stores plain int internally.  Presents a std::vector-like interface with
// IntClass<Tag> element type via a proxy reference.
//
// Key benefit: .as_eigen() returns Map<const IntVect> for Eigen fancy-index
// operations (voltage scatter, etc.) — int*→int*, no reinterpret_cast.
// ============================================================================
template<int Tag>
class TaggedIdStdVec {
public:
    using value_type = IntClass<Tag>;
    using size_type  = std::size_t;

    // ----- proxy reference -------------------------------------------------
    struct reference {
        int& v_;
        explicit reference(int& v) noexcept : v_(v) {}
        reference& operator=(IntClass<Tag> x) noexcept { v_ = x.cast_int(); return *this; }
        reference& operator=(int         x) noexcept { v_ = x;             return *this; }
        operator IntClass<Tag>() const noexcept { return IntClass<Tag>(v_); }
        explicit operator int() const noexcept  { return v_; }
        int cast_int()          const noexcept  { return v_; }
    };

    // ----- forward const iterator ------------------------------------------
    struct const_iterator {
        const int* ptr_;
        explicit const_iterator(const int* p) noexcept : ptr_(p) {}
        IntClass<Tag>   operator*()  const noexcept { return IntClass<Tag>(*ptr_); }
        const_iterator& operator++()       noexcept { ++ptr_; return *this; }
        const_iterator  operator++(int)    noexcept { auto t = *this; ++ptr_; return t; }
        bool operator==(const const_iterator& o) const noexcept { return ptr_ == o.ptr_; }
        bool operator!=(const const_iterator& o) const noexcept { return ptr_ != o.ptr_; }
    };

    // ----- constructors ----------------------------------------------------
    TaggedIdStdVec() = default;

    // Analogous to std::vector<IntClass<Tag>>(n, val)
    explicit TaggedIdStdVec(size_type n, IntClass<Tag> val)
        : data_(n, val.cast_int()) {}

    // // Construct from a std::vector of tagged int
    // explicit TaggedIdStdVec(const std::vector<IntClass<Tag>>& v) {
    //     data_.reserve(v.size());
    //     for (const auto& el : v) data_.push_back(el.cast_int());
    // }

    // Construct from a std::vector of int
    explicit TaggedIdStdVec(const std::vector<int>& v) {
        data_ = v;
    }

    // Construct from an IntVect (eigen vector of int)
    explicit TaggedIdStdVec(const IntVect & v) {
        data_= std::vector<int>(v.begin(), v.end());
    }

    // ----- std::vector interface -------------------------------------------
    size_type size()  const noexcept { return data_.size(); }
    bool      empty() const noexcept { return data_.empty(); }
    void clear()              { data_.clear(); }
    void reserve(size_type n) { data_.reserve(n); }
    void push_back(IntClass<Tag> val) { data_.push_back(val.cast_int()); }

    // signed index — used pervasively in loops with int variables
    IntClass<Tag> operator[](int i) const noexcept {
        return IntClass<Tag>(data_[static_cast<size_type>(i)]);
    }
    reference operator[](int i) noexcept {
        return reference(data_[static_cast<size_type>(i)]);
    }
    // unsigned index
    IntClass<Tag> operator[](size_type i) const noexcept { return IntClass<Tag>(data_[i]); }
    reference     operator[](size_type i)       noexcept { return reference(data_[i]); }

    // Eigen-style operator() — same semantics as operator[], for compatibility
    // with call sites that use () notation (previously Eigen::Matrix-backed types)
    IntClass<Tag> operator()(int i)       const noexcept { return (*this)[i]; }
    reference     operator()(int i)             noexcept { return (*this)[i]; }
    IntClass<Tag> operator()(size_type i) const noexcept { return (*this)[i]; }
    reference     operator()(size_type i)       noexcept { return (*this)[i]; }

    // ----- iteration -------------------------------------------------------
    const_iterator begin() const noexcept { return const_iterator(data_.data()); }
    const_iterator end()   const noexcept { return const_iterator(data_.data() + data_.size()); }

    // ----- UB-free accessors -----------------------------------------------
    // Returns an Eigen Map for Eigen fancy indexing (voltage scatter, etc.)
    Eigen::Map<const IntVect> as_eigen() const noexcept {
        return Eigen::Map<const IntVect>(data_.data(),
                                        static_cast<Eigen::Index>(data_.size()));
    }
    // Raw pointer for C APIs / pybind11
    const int*       raw_data()      const noexcept { return data_.data(); }
    // Copy to std::vector<int> for numpy() accessors
    const std::vector<int> & to_int_vector() const          { return data_; }

private:
    std::vector<int> data_;
};


// ---- Type aliases ----------------------------------------------------------
// Type aliases
// GlobalBusIdVect: std::vector-backed (TaggedIdStdVec) — used for mappings:
//   id_solver_to_me, id_ac_solver_to_me_, id_dc_solver_to_me_, …
// SolverBusIdVect: std::vector-backed (TaggedIdStdVec) — used for mappings:
//   id_me_to_solver, id_grid_to_solver, id_me_to_ac_solver_, …
// GridModelBusIdVect: Eigen-backed (TaggedIdVec) — alias for TaggedIdVec<1>,
//   used for Eigen-style bus-id lists (bus_id_, slack_bus_id_*_me_, from_bus, …)

using GlobalBusIdVect    = TaggedIdStdVec<GLOBAL_BUS>;
using GridModelBusIdVect = TaggedIdStdVec<GRIDMODEL_BUS>;
using SolverBusIdVect    = TaggedIdStdVec<SOLVER_BUS>;

#endif // TAGGEDIDVEC_HPP
