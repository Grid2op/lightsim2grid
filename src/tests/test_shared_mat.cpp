// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// SharedMat: the holder of the batch results that python arrays co-own. Whatever the
// holder does next, a buffer handed out by share() keeps its values.

#include <catch2/catch_test_macros.hpp>

#include "Utils.hpp"
#include "batch_algorithm/SharedMat.hpp"

using namespace ls2g;
using Mat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

TEST_CASE("SharedMat overwrite reuses an unshared buffer and leaves a shared one alone", "[shared_mat]")
{
    SharedMat<Mat> m;
    m.overwrite() = Mat::Constant(2, 3, 1.);
    const real_type * first = m.get().data();
    m.overwrite().setConstant(2.);                // nothing holds it: written in place
    REQUIRE(m.get().data() == first);

    auto held = m.share();
    m.overwrite() = Mat::Constant(2, 3, 3.);      // held: a new buffer
    REQUIRE(m.get().data() != held->data());
    REQUIRE((held->array() == 2.).all());
    REQUIRE((m.get().array() == 3.).all());
}

TEST_CASE("SharedMat mut copies a shared buffer before writing", "[shared_mat]")
{
    SharedMat<Mat> m;
    m.overwrite() = Mat::Constant(2, 2, 1.);
    auto held = m.share();
    m.mut()(0, 0) = 5.;
    REQUIRE((*held)(0, 0) == 1.);
    REQUIRE(m.get()(0, 0) == 5.);
    REQUIRE(m.get()(1, 1) == 1.);                 // the other values came along
}

TEST_CASE("SharedMat reset and destruction keep what was shared", "[shared_mat]")
{
    std::shared_ptr<const Mat> held;
    {
        SharedMat<Mat> m;
        m.overwrite() = Mat::Constant(3, 1, 7.);
        held = m.share();
        m.reset();
        REQUIRE(m.get().size() == 0);
    }
    REQUIRE(held->size() == 3);
    REQUIRE((held->array() == 7.).all());
}

TEST_CASE("SharedMat copies deeply and leaves a moved-from holder empty", "[shared_mat]")
{
    SharedMat<Mat> a;
    a.overwrite() = Mat::Constant(1, 2, 4.);
    SharedMat<Mat> b(a);
    REQUIRE(b.get().data() != a.get().data());
    b.raw()(0, 0) = 0.;
    REQUIRE(a.get()(0, 0) == 4.);

    SharedMat<Mat> c(std::move(a));
    REQUIRE(c.get()(0, 1) == 4.);
    REQUIRE(a.get().size() == 0);                 // usable, not null
}
