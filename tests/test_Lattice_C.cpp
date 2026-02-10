//
// Created by pbialas on 20.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include "Field/Lattice.h"

using namespace lft;

TEST_CASE("Latttice constructors indexing 1D C", "[Constructors] [Lattice][1D]{indexing") {
    const uint8_t Lx = 8;

    Lattice<int8_t, 1> lattice({Lx}, 'C');
    REQUIRE(lattice.n_elements == 8);

    for (int i = 0; i < Lx; i++) {
        INFO("i = " <<i);
        REQUIRE(int(lattice.up(i,0))==(i+1)%Lx);
        REQUIRE(int(lattice.dn(i,0))==(i-1+Lx)%Lx);
    }
}

TEST_CASE("Latttice constructors C", "[Constructors] [Lattice]") {
    const int Lx = 8;
    const int Ly = 7;

    Lattice<int8_t> lattice({Lx, Ly}, 'C');
    REQUIRE(lattice.n_elements == Lx*Ly);

    for (int8_t i = 0; i < Lx; i++) {
        for (int8_t j = 0; j < Ly; j++) {
            INFO("i = " <<i << " j = " << j);
            auto idx = j + i * Ly;

            REQUIRE(lattice.idx({i, j})==idx);
            REQUIRE(lattice.idx({i, j})==idx);
        }
    }
}

TEST_CASE("Lattice indexing C", "[Lattice][Indexing]") {
    const int Lx = 3;
    const int Ly = 2;

    Lattice<int8_t> lattice({Lx, Ly}, 'C');
    REQUIRE(lattice.n_elements == Lx*Ly);


    // [0 1]
    // [2 3]
    // [4 5]


    REQUIRE(lattice.idx({0, 0}) == 0);
    REQUIRE(lattice.idx({0, 1}) == 1);
    REQUIRE(lattice.idx({1, 0}) == 2);
    REQUIRE(lattice.idx({1, 1}) == 3);
    REQUIRE(lattice.idx({2, 0}) == 4);
    REQUIRE(lattice.idx({2, 1}) == 5);
}


TEST_CASE("Latttice up dn C", "[Constructors][Lattice]") {
    {
        const int Lx = 2;
        const int Ly = 3;

        Lattice<int8_t> lattice({Lx, Ly}, 'C');
        REQUIRE(lattice.n_elements == Lx * Ly);

        // [0 1 2 ]
        // [3 4 5 ]

        REQUIRE(lattice.up(0, 0) == 3);
        REQUIRE(lattice.up(0, 1) == 1);
        REQUIRE(lattice.up(1, 0) == 4);
        REQUIRE(lattice.up(1, 1) == 2);
        REQUIRE(lattice.up(2, 0) == 5);
        REQUIRE(lattice.up(2, 1) == 0);
        REQUIRE(lattice.up(3, 0) == 0);
        REQUIRE(lattice.up(3, 1) == 4);
        REQUIRE(lattice.up(4, 0) == 1);
        REQUIRE(lattice.up(4, 1) == 5);
        REQUIRE(lattice.up(5, 0) == 2);
        REQUIRE(lattice.up(5, 1) == 3);

        REQUIRE(lattice.dn(0, 0) == 3);
        REQUIRE(lattice.dn(0, 1) == 2);
        REQUIRE(lattice.dn(1, 0) == 4);
        REQUIRE(lattice.dn(1, 1) == 0);
        REQUIRE(lattice.dn(2, 0) == 5);
        REQUIRE(lattice.dn(2, 1) == 1);
        REQUIRE(lattice.dn(3, 0) == 0);
        REQUIRE(lattice.dn(3, 1) == 5);
        REQUIRE(lattice.dn(4, 0) == 1);
        REQUIRE(lattice.dn(4, 1) == 3);
        REQUIRE(lattice.dn(5, 0) == 2);
        REQUIRE(lattice.dn(5, 1) == 4);
    }
}

TEST_CASE("Latttice up dn 4x4 C", "[Constructors][Lattice]") {
    const int Lx = 3;
    const int Ly = 4;

    Lattice<int8_t> lattice({Lx, Ly}, 'C');
    REQUIRE(lattice.n_elements == Lx*Ly);

    // [ 0  1  2  3 ]
    // [ 4  5  6  7 ]
    // [ 8  9 10 11 ]

    REQUIRE(lattice.up(0, 0) == 4);
    REQUIRE(lattice.up(0, 1) == 1);

    REQUIRE(lattice.up(7, 0) == 11);
    REQUIRE(lattice.up(7, 1) == 4);

    REQUIRE(lattice.up(11, 0) == 3);
    REQUIRE(lattice.up(11, 1) == 8);


    REQUIRE(lattice.dn(0, 0) == 8);
    REQUIRE(lattice.dn(0, 1) == 3);
}

TEST_CASE("Lattice even/odd C", "[Lattice][Even/Odd]") {
    const int Lx = 4;
    const int Ly = 6;

    Lattice<int8_t> lattice({Lx, Ly}, 'C');
    REQUIRE(lattice.n_elements == Lx*Ly);

    // [ 0  1  2  3 4  5 ]
    // [ 6  7  8  9 10 11]
    // [12 13 14 15 16 17]
    // [18 19 20 21 22 23]

    REQUIRE(lattice.even(0) == 0);
    REQUIRE(lattice.even(1) == 2);
    REQUIRE(lattice.even(2) == 4);
    REQUIRE(lattice.even(3) == 7);
    REQUIRE(lattice.even(4) == 9);
    REQUIRE(lattice.even(5) == 11);
    REQUIRE(lattice.even(6) == 12);
    REQUIRE(lattice.even(7) == 14);
    REQUIRE(lattice.even(8) == 16);
    REQUIRE(lattice.even(9) == 19);
    REQUIRE(lattice.even(10) == 21);
    REQUIRE(lattice.even(11) == 23);

    REQUIRE(lattice.odd(0) == 1);
    REQUIRE(lattice.odd(1) == 3);
    REQUIRE(lattice.odd(2) == 5);
    REQUIRE(lattice.odd(3) == 6);
    REQUIRE(lattice.odd(4) == 8);
    REQUIRE(lattice.odd(5) == 10);
    REQUIRE(lattice.odd(6) == 13);
    REQUIRE(lattice.odd(7) == 15);
    REQUIRE(lattice.odd(8) == 17);
    REQUIRE(lattice.odd(9) == 18);
    REQUIRE(lattice.odd(10) == 20);
    REQUIRE(lattice.odd(11) == 22);
}
