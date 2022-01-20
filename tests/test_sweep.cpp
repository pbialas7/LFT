//
// Created by pbialas on 20.01.2022.
//

//
// Created by pbialas on 20.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include "Field/Field.h"
#include "MonteCarlo/sweep.h"

template<typename F>
struct MakeOne {
    typename F::lattice_t::size_t operator()(F &field, typename F::lattice_t::size_t site) {
        field[site] = 1;
        return 1;
    }
};

TEST_CASE("sweep", "[sweep]") {
    const int Lx = 8;
    const int Ly = 8;

    Lattice<> lat({Lx, Ly});

    Field<int8_t> field(lat, -1);
    for (auto i = 0; i < field.n_elements; i++) {
        REQUIRE(field[i] == -1);
    }
    MakeOne<Field<int8_t>> update;
    auto acceptance = sweep(field, update);
    REQUIRE(acceptance == Lx * Ly);
    for (auto i = 0; i < field.n_elements; i++) {
        REQUIRE(int(field[i]) == 1);
    }

}