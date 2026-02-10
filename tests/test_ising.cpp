//
// Created by pbialas on 20.01.2022.
//

#include  "catch2/catch_test_macros.hpp"

#include <iostream>

#include "Field/Field.h"
#include "MonteCarlo/sweep.h"
#include "Ising/ising.h"


using namespace lft;

TEST_CASE("ising sweep", "[ising][sweep]") {
    const int Lx = 8;
    const int Ly = 8;

    using lattice_t = Lattice<>;
    lattice_t lat({Lx, Ly});

    std::mt19937_64 rng(4132413);

    {
        ising::IsingField<lattice_t> field(lat, 1);
        ising::HeathBath<lattice_t> update(0.0, rng);
        auto acceptance = sweep(field, update);
        REQUIRE(acceptance == Lx * Ly);

    }

    {
        ising::IsingField<lattice_t> field(lat, 1);
        ising::HeathBath<lattice_t> update(1.0, rng);
        auto acceptance = sweep(field, update);
        REQUIRE(acceptance == Lx * Ly);

    }
    {
        ising::IsingField<lattice_t> field(lat);
        ising::hot_start(field, rng);
        std::cout << mean<double>(field) << " ";
        ising::HeathBath<lattice_t> update(1.0, rng);
        for (auto i = 0; i < 100; ++i) {
            auto acceptance = sweep(field, update);
            REQUIRE(acceptance == Lx * Ly);
        }

    }

}


TEST_CASE("ising sweep bigger lattice", "[ising][sweep]") {
    const int Lx = 64;
    const int Ly = 64;

    using lattice_t = Lattice<>;
    lattice_t lat({Lx, Ly});

    std::mt19937_64 rng(4132413);

    {
        ising::IsingField<lattice_t> field(lat);
        ising::hot_start(field, rng);
        std::cout << mean<double>(field) << " ";
        ising::HeathBath<lattice_t> update(1.0, rng);
        for (auto i = 0; i < 1000; ++i) {
            auto acceptance = sweep(field, update);
            REQUIRE(acceptance == Lx * Ly);
        }

        double M = 0.0, E = 0.0;
        for (auto i = 0; i < 100; ++i) {
            auto acceptance = sweep(field, update);
            M += ising::magnetisation<double>(field);
            E += ising::energy<double>(field);
        }
        M /= 100;
        E /= 100;
        std::cout.precision(12);
        std::cout << M << " " << E << std::endl;
    }
}