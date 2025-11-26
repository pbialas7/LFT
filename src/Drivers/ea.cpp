//
// Created by pbialas on 25.11.2025.
//

#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include <Field/Lattice.h>

#include <MonteCarlo/sweep.h>
#include <EdwardsAnderson/ea.h>
#include "utils/fs.h"
#include "ising_base_options.h"


int main(int argc, char *argv[]) {
    IsingBaseOptions base_options;

    auto results = base_options.cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return 1;
    }

    spdlog::info("Lx {} Ly {}", base_options.Lx, base_options.Ly);

    std::mt19937_64 rng(base_options.seed);
    using lattice_t = Lattice<uint32_t>;
    lattice_t lat({base_options.Lx, base_options.Ly});
    ea::SpinField<lattice_t> spin_field(lat, 1);

    Lattice<uint32_t, 3> j_lat({lat.dims[0], lat.dims[1], 2});
    ea::JField<lattice_t> j_field(j_lat, 1);

    auto j_path = make_file_path(base_options.data_dir, "j", base_options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    ea::HeathBath<lattice_t> update(base_options.beta, rng, j_field);

    for (int i = 0; i < base_options.n_term; ++i) {
        sweep(spin_field, update);
    }

    for (int i = 0; i < base_options.n_sweeps; ++i) {
        sweep(spin_field, update);
    }

    return 0;
}
