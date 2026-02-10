//
// Created by pbialas on 20.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include "Field/Field.h"

TEST_CASE("Field Constructors", "[Field][Constructors]") {
    const int Lx = 8;
    const int Ly = 8;

    Lattice<> lat({Lx, Ly});

    Field<int8_t, decltype(lat)> field1(lat);
    REQUIRE(field1.n_elements == 64);
    Field<int8_t, decltype(lat)> field2(lat, -1);
    REQUIRE(field1.n_elements == 64);
    for (auto i = 0; i < field2.n_elements; i++) {
        REQUIRE(field2[i] == -1);
    }
}

TEST_CASE("Field Corona", "[Field][Corona]") {
    const int Lx = 4;
    const int Ly = 4;

    Lattice<int8_t> lat({Lx, Lx});

    Field<int8_t, decltype(lat)> field(lat);
    REQUIRE(field.n_elements == Lx * Ly);


    // [12 13 14 15]
    // [ 8  9 10 11]
    // [ 4  5  6  7]
    // [ 0  1  2  3]

    for (auto i = 0; i < field.n_elements; i++) {
        field[i] = i;
    }

    REQUIRE(int(field.corona(0)) == 20);
    REQUIRE(int(field.corona(15)) == (11+3+14+12));
    REQUIRE(int(field.corona(5)) == (4+6+1+9));

}
