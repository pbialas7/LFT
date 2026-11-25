//
// Created by pbialas on 20.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include "Field/Lattice.h"

using namespace lft;

TEST_CASE("Latttice constructors indexing 1D", "[Constructors] [Lattice][1D]{indexing") {
    const uint8_t Lx = 8;

    Lattice<int8_t, 1> lattice({Lx});
    REQUIRE(lattice.n_elements == 8);

    for (int i = 0; i < Lx; i++) {
        INFO("i = " <<i);
        REQUIRE(int(lattice.up(i,0))==(i+1)%Lx);
        REQUIRE(int(lattice.dn(i,0))==(i-1+Lx)%Lx);
    }
}

TEST_CASE("Latttice constructors", "[Constructors] [Lattice]") {
    const int Lx = 8;
    const int Ly = 7;

    Lattice<int8_t> lattice({Lx, Ly});
    REQUIRE(lattice.n_elements == Lx*Ly);

    for (int8_t i = 0; i < Lx; i++) {
        for (int8_t j = 0; j < Ly; j++) {
            INFO("i = " <<i << " j = " << j);
            auto idx = i + j * Lx;

            REQUIRE(lattice.idx({i, j})==idx);
            REQUIRE(lattice.idx({i, j})==idx);
        }
    }
}

TEST_CASE("Lattice indexing", "[Lattice][Indexing]") {
    const int Lx = 2;
    const int Ly = 2;

    Lattice<int8_t> lattice({Lx, Ly});
    REQUIRE(lattice.n_elements == 4);

    // [ 2 3 ]
    // [ 0 1 ]


    REQUIRE(lattice.idx({0, 0}) == 0);
    REQUIRE(lattice.idx({0, 1}) == 2);
    REQUIRE(lattice.idx({1, 0}) == 1);
    REQUIRE(lattice.idx({1, 1}) == 3);
}


TEST_CASE("Latttice up dn", "[Constructors][Lattice]") {
    {
        const int Lx = 2;
        const int Ly = 2;

        Lattice<int8_t> lattice({Lx, Ly});
        REQUIRE(lattice.n_elements == 4);

        // [ 2 3 ]
        // [ 0 1 ]

        REQUIRE(lattice.up(0, 0) == 1);
        REQUIRE(lattice.up(0, 1) == 2);
        REQUIRE(lattice.up(1, 0) == 0);
        REQUIRE(lattice.up(1, 1) == 3);
        REQUIRE(lattice.up(2, 0) == 3);
        REQUIRE(lattice.up(2, 1) == 0);
        REQUIRE(lattice.up(3, 0) == 2);
        REQUIRE(lattice.up(3, 1) == 1);

        REQUIRE(lattice.dn(0, 0) == 1);
        REQUIRE(lattice.dn(0, 1) == 2);
        REQUIRE(lattice.dn(1, 0) == 0);
        REQUIRE(lattice.dn(1, 1) == 3);
        REQUIRE(lattice.dn(2, 0) == 3);
        REQUIRE(lattice.dn(2, 1) == 0);
        REQUIRE(lattice.dn(3, 0) == 2);
        REQUIRE(lattice.dn(3, 1) == 1);
    }
}

TEST_CASE("Latttice up dn 4x4", "[Constructors][Lattice]") {
    const int Lx = 4;
    const int Ly = 4;

    Lattice<int8_t> lattice({Lx, Ly});
    REQUIRE(lattice.n_elements == 16);

    // [12 13 14 15]
    // [ 8  9 10 11]
    // [ 4  5  6  7]
    // [ 0  1  2  3]

    REQUIRE(lattice.up(0, 0) == 1);
    REQUIRE(lattice.up(0, 1) == 4);

    REQUIRE(lattice.up(7, 0) == 4);
    REQUIRE(lattice.up(7, 1) == 11);

    REQUIRE(lattice.dn(12, 0) == 15);
    REQUIRE(lattice.up(12, 1) == 0);
}

TEST_CASE("Lattice even/odd", "[Lattice][Even/Odd]") {
    const int Lx = 4;
    const int Ly = 4;

    Lattice<int8_t> lattice({Lx, Ly});
    REQUIRE(lattice.n_elements == 16);

    // [12 13 14 15]
    // [ 8  9 10 11]
    // [ 4  5  6  7]
    // [ 0  1  2  3]

    REQUIRE(lattice.even(0) == 0);
    REQUIRE(lattice.even(1) == 2);
    REQUIRE(lattice.even(2) == 5);
    REQUIRE(lattice.even(3) == 7);
    REQUIRE(lattice.even(4) == 8);
    REQUIRE(lattice.even(5) == 10);
    REQUIRE(lattice.even(6) == 13);
    REQUIRE(lattice.even(7) == 15);

    REQUIRE(lattice.odd(0) == 1);
    REQUIRE(lattice.odd(1) == 3);
    REQUIRE(lattice.odd(2) == 4);
    REQUIRE(lattice.odd(3) == 6);
    REQUIRE(lattice.odd(4) == 9);
    REQUIRE(lattice.odd(5) == 11);
    REQUIRE(lattice.odd(6) == 12);
    REQUIRE(lattice.odd(7) == 14);
}


TEST_CASE("Lattice checkerboard requires even extents", "[Lattice][Even/Odd]") {
    SECTION("even extents are accepted and neighbours have opposite parity") {
        Lattice<uint32_t> lattice({4, 6});
        REQUIRE(lattice.has_checkerboard());
        REQUIRE_NOTHROW(lattice.require_checkerboard());
        REQUIRE(lattice.n_elements == 24);
        std::vector<bool> is_even(lattice.n_elements, false);
        for (std::size_t i = 0; i < lattice.n_elements / 2; ++i) is_even[lattice.even(i)] = true;
        for (std::size_t i = 0; i < lattice.n_elements / 2; ++i) {
            auto s = lattice.even(i);
            for (int d = 0; d < 2; ++d) {
                REQUIRE_FALSE(is_even[lattice.up(s, d)]);
                REQUIRE_FALSE(is_even[lattice.dn(s, d)]);
            }
        }
    }
    SECTION("any odd extent is rejected") {
        REQUIRE_THROWS_AS(Lattice<uint32_t>({3, 4}).require_checkerboard(), std::invalid_argument);
        REQUIRE_THROWS_AS(Lattice<uint32_t>({4, 3}).require_checkerboard(), std::invalid_argument);
        REQUIRE_THROWS_AS(Lattice<uint32_t>({3, 3}).require_checkerboard(), std::invalid_argument);
        REQUIRE_THROWS_AS((Lattice<uint32_t, 1>({5}).require_checkerboard()), std::invalid_argument);
        REQUIRE_FALSE(Lattice<uint32_t>({3, 4}).has_checkerboard());
    }
}
