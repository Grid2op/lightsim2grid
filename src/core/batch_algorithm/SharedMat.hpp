// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LS2G_SHARED_MAT_H
#define LS2G_SHARED_MAT_H

#include <memory>
#include <utility>

namespace ls2g {

/**
 * A result matrix of a batch algorithm, held through a shared_ptr.
 *
 * share() hands out co-ownership of the current buffer: the python bindings put it
 * in the base of the numpy array they return, so that array stays valid for as long
 * as it lives. Replacing the result never writes over, nor frees, a buffer someone
 * else still holds: reset() and overwrite() install a new one instead, and mut()
 * copies before writing. An array therefore keeps the values of the run that
 * produced it.
 *
 * The write accessors, from the cheapest:
 *   - raw():       the current buffer, in place. Only for the code filling a buffer
 *                  that reset() / overwrite() just installed, before anything could
 *                  share it -- the per-row writes of a compute(), from its threads.
 *   - overwrite(): a buffer to replace entirely: the current one if nothing else
 *                  holds it (no allocation), a new empty one otherwise.
 *   - mut():       the current values, to modify: copied first if shared.
 *   - reset():     always a new, empty buffer (what clear() wants: the memory goes
 *                  as soon as no array holds it).
 *
 * Copying a holder copies the values (two algorithms never write into one buffer);
 * a moved-from holder is left empty, not null.
 */
template<class M>
class SharedMat
{
    public:
        SharedMat() : ptr_(std::make_shared<M>()) {}
        SharedMat(const SharedMat & other) : ptr_(std::make_shared<M>(*other.ptr_)) {}
        SharedMat(SharedMat && other) : ptr_(std::move(other.ptr_)) { other.ptr_ = std::make_shared<M>(); }
        SharedMat & operator=(const SharedMat & other){
            if(this != &other) ptr_ = std::make_shared<M>(*other.ptr_);
            return *this;
        }
        SharedMat & operator=(SharedMat && other){
            if(this != &other){
                ptr_ = std::move(other.ptr_);
                other.ptr_ = std::make_shared<M>();
            }
            return *this;
        }

        const M & get() const { return *ptr_; }
        M & raw() { return *ptr_; }
        M & overwrite() {
            if(ptr_.use_count() > 1) ptr_ = std::make_shared<M>();
            return *ptr_;
        }
        M & mut() {
            if(ptr_.use_count() > 1) ptr_ = std::make_shared<M>(*ptr_);
            return *ptr_;
        }
        M & reset() {
            ptr_ = std::make_shared<M>();
            return *ptr_;
        }

        std::shared_ptr<const M> share() const { return ptr_; }

    private:
        std::shared_ptr<M> ptr_;
};

}  // namespace ls2g

#endif  // LS2G_SHARED_MAT_H
