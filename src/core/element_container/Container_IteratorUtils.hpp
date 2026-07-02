// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CONTAINER_ITERATOR_UTILS_H
#define CONTAINER_ITERATOR_UTILS_H

#include <algorithm>  // for std::find

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "BaseConstants.hpp"

namespace ls2g {

// TODO simplify this iterator classes

// iterator type
template<class DataType>
class GenericContainerConstIterator
{
    protected:
        using DataInfo = typename DataType::DataInfo;

        const DataType * const _p_data_;
        int my_id;

    public:
        DataInfo my_info;

        // functions
        GenericContainerConstIterator(const DataType * const data_, int id):
            _p_data_(data_),
            my_id(id),
            my_info(*data_, id)
            {};

        // Returns a COPY, not a reference: `my_info` is a single member reused/overwritten
        // in place by every `operator++`/`operator--` (see below), so a caller holding onto
        // a `const DataInfo&` across more than one increment would silently see it mutate
        // out from under them. This bit pybind11's `__iter__` binding hard: `py::make_iterator`
        // dereferences+advances the SAME C++ iterator for every `__next__()`, so
        // `list(container)` (which exhausts the iterator before Python is done with any of
        // the yielded objects) ended up handing back N aliases of the SAME final `my_info`
        // state (mostly the past-the-end default-constructed one, id=-1) instead of N
        // independent per-element snapshots. Returning by value forces pybind11 to copy at
        // each step; `model.get_X()[i]` (`operator[]`, unaffected) was always correct.
        DataInfo operator*() const { return my_info; }
        bool operator==(const GenericContainerConstIterator<DataType> & other) const { return (my_id == other.my_id) && (_p_data_ == other._p_data_); }
        bool operator!=(const GenericContainerConstIterator<DataType> & other) const { return !(*this == other); }
        GenericContainerConstIterator<DataType> & operator++()
        {
            ++my_id;
            my_info = DataInfo(*_p_data_, my_id);
            return *this;
        }
        GenericContainerConstIterator<DataType> & operator--()
        {
            --my_id;
            my_info = DataInfo(*_p_data_, my_id);
            return *this;
        }
        int size() const { return _p_data_->nb(); }
};

/**
 * CRTP to avoid to repeat the following piece of code to make all 
 * derived class "iterable"
 * (see https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)
 * 
 * Example of use: 
 * 
 * class SubstationContainer: public IteratorAdder<SubstationContainer, SubstationInfo> 
 * {
 *    ...
 * };
 * 
 * And then you can use: 
 * 
 * substations_.begin(), substations_.end(), substations_[xxx], for(const auto & sub : substations_) etc.
 */
template<class ConcreteContainer, class ConcreteContainerInfo>
class IteratorAdder
{
    private:
        using DataConstIterator = GenericContainerConstIterator<ConcreteContainer>;

    public:
        DataConstIterator begin() const {return DataConstIterator(static_cast<const ConcreteContainer*>(this), 0); }
        DataConstIterator end() const {return DataConstIterator(
            static_cast<const ConcreteContainer*>(this), 
            static_cast<const ConcreteContainer*>(this)->nb()); 
        }
        ConcreteContainerInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative element id.");
            }
            if(id >= static_cast<const ConcreteContainer*>(this)->nb())
            {
                throw std::range_error("Load out of bound. Not enough elements of this type on the grid.");
            }
            return ConcreteContainerInfo(static_cast<const ConcreteContainer&>(*this), id);
        }
};


} // namespace ls2g

#endif // CONTAINER_ITERATOR_UTILS_H
