//
// Created by pbialas on 25.11.2025.
//

#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include <Field/Lattice.h>

#include "MonteCarlo/sweep.h"
#include "EdwardsAnderson/ea.h"
#include "utils/fs.h"
#include "ising_base_options.h"


int main(int argc, char *argv[]) {
    IsingBaseOptions base_options;

    int meas_freq = 0;
    bool ising = false;
    bool binary = false;


    base_options.cli |= lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("configuration save frequency");
    base_options.cli |= lyra::opt(ising)["--ising"]("Set J = 1");
    base_options.cli |= lyra::opt(binary)["--binary"]("Sets J =+/-1");

    auto results = base_options.cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return 1;
    }

    spdlog::info("Lx {} Ly {}", base_options.Lx, base_options.Ly);

    std::mt19937_64 rng(base_options.seed);
    std::bernoulli_distribution bern(0.5);
    std::normal_distribution<double> normal(0.0, 1.0);
    using lattice_t = lft::Lattice<uint32_t>;
    lattice_t lat({base_options.Lx, base_options.Ly},'C');
    ea::SpinField<lattice_t> spin_field(lat, 1);

    lft::Lattice<uint32_t, 3> j_lat({2,lat.dims[0], lat.dims[1]},'C');
    auto j_field=lft::make_field(j_lat, 1.0);
    if (!ising) {
        if (binary) {
            for (int i = 0; i < j_field.n_elements; ++i) {
                if (!bern(rng)) {
                    j_field[i] = -1;
                }
            }
        } else {
            for (int i = 0; i < j_field.n_elements; ++i) {
                j_field[i] = normal(rng);
            }
        }
    }

    auto j_path = make_file_path(base_options.data_dir, "j", base_options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    ea::HeathBath<lattice_t> update(base_options.beta, rng, j_field);

    for (int i = 0; i < base_options.n_term; ++i) {
        sweep(spin_field, update);
    }

    auto *em_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "em", base_options.name, "txt"),
        meas_freq > 0, std::fstream::out);
    auto cfg_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "cfg", base_options.name, "bin"),
        base_options.save_freq > 0, std::ios::out | std::ios::binary);


    for (int i = 0; i < base_options.n_sweeps; ++i) {
        sweep(spin_field, update);
        if (meas_freq > 0 && (i % meas_freq) == 0) {
            if (em_stream_ptr) {
                *em_stream_ptr << ea::energy<double>(spin_field, j_field) << " ";
                *em_stream_ptr << ea::magnetisation<double>(spin_field) << "\n";
            }
        }

        if (base_options.save_freq > 0 && (i % base_options.save_freq) == 0) {
            if (cfg_stream_ptr) {
                spin_field.write(*cfg_stream_ptr);
            }
        }
    }
    if (em_stream_ptr)
        delete em_stream_ptr;
    if (cfg_stream_ptr)
        delete cfg_stream_ptr;

    return 0;
}
