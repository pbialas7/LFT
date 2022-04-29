//
// Created by pbialas on 24.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include "Field/Lattice.h"

TEST_CASE("MultiIndex constructors", "[Constructors] [MultiIndex]") {
    const int Lx = 8;
    const int Ly = 8;

    MultiIndex<int8_t, 2> index({Lx, Ly});
    REQUIRE(index.dims[0] == Lx);
    REQUIRE(index.dims[0] == Ly);
    REQUIRE(index.volumes_[0] == 1);
    REQUIRE(index.volumes_[1] == Lx);
    REQUIRE(index.volumes_[2] == Lx * Ly);
    REQUIRE(index.volume() == 64);

}

TEST_CASE("MultiIndex indexing", "[MultiIndex][Indexing]") {
    const int Lx = 3;
    const int Ly = 2;

    MultiIndex<int8_t, 2> index({Lx, Ly});

    REQUIRE(index.volume() == 6);

    // [ 3 4 5]
    // [ 0 1 2]

    for (int i = 0; i < Lx * Ly; ++i) {
        INFO("i = " << i);
        REQUIRE(index.idx() == i);
        ++index;
    }

    {
        MultiIndex<int8_t, 2> index({Lx, Ly}, {2, 1});
        REQUIRE(index.coords()[0] == 2);
        REQUIRE(index.coords()[1] == 1);

        REQUIRE(index.volumes_[0] == 1);
        REQUIRE(index.volumes_[1] == Lx);
        REQUIRE(index.idx() == 5);
    }

    {
        MultiIndex<int8_t, 2> index({Lx, Ly});
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 0);
            REQUIRE(coords[1] == 0);
        }
        ++index;
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 1);
            REQUIRE(coords[1] == 0);
        }
        ++index;
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 2);
            REQUIRE(coords[1] == 0);
        }
        ++index;
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 0);
            REQUIRE(coords[1] == 1);
        }
        ++index;
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 1);
            REQUIRE(coords[1] == 1);
        }
        ++index;
        {
            auto coords = index.coords();
            REQUIRE(coords[0] == 2);
            REQUIRE(coords[1] == 1);
        }

    }
}

TEST_CASE("MultiIndex inc", "[MultiIndex][inc]") {
    const int Lx = 4;
    const int Ly = 4;

    MultiIndex<int8_t, 2> index({Lx, Ly});

    REQUIRE(index.volume() == Lx * Ly);


    for (int i = 0; i < Lx * Ly; ++i, ++index) {
        INFO("i = " << i);
        REQUIRE(index.idx() == i);

    }
}
