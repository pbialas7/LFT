//
// Created by pbialas on 25.11.2025.
//

#include <Field/Lattice.h>

#include "MonteCarlo/sweep.h"
#include <EdwardsAnderson/ea.h>

int main() {
    std::mt19937_64 rng(1423);
    using lattice_t = Lattice<uint32_t>;
    lattice_t lat({16,16});
    ea::SpinField<lattice_t> spin_field(lat,1);

    Lattice<uint32_t,3> j_lat({16,16,2});
    ea::JField<lattice_t> j_field(j_lat,1);

    ea::HeathBath<lattice_t> update(0.1,rng, j_field);

    for (int i =0;i<100000;++i) {
        sweep(spin_field, update);
    }

    return 0;
}