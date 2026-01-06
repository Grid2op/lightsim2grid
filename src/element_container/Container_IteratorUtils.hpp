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

// TODO simplify this iterator classes

// iterator type
template<class DataType>
class GenericContainerConstIterator
{
    protected:
        typedef typename DataType::DataInfo DataInfo;

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

        const DataInfo& operator*() const { return my_info; }
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
        typedef GenericContainerConstIterator<ConcreteContainer> DataConstIterator;
        // typedef typename ConcreteContainer::DataConstIterator DataConstIterator;
        // typedef typename ConcreteContainer::DataInfo DataInfo;
        // using DataContainerImpl::nb; 

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

#endif // CONTAINER_ITERATOR_UTILS_H

