//
// Created by pbialas on 25.11.2025.
//

#include <fstream>
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;

#include <omp.h>

#include <Field/Lattice.h>

#include "MonteCarlo/sweep.h"
#include "EdwardsAnderson/ea.h"
#include "utils/fs.h"
#include "ising_base_options.h"


int main(int argc, char *argv[]) {
    auto max_threads = omp_get_max_threads();

    spdlog::info("Max number of threads = {}", max_threads);


    IsingBaseOptions base_options;

    int meas_freq = 0;
    bool ising = false;
    bool binary = false;
    bool two_replicas = false;
    std::string j_file_path;
    int n_replicas = 1;
    int n_threads = 1;

    base_options.cli |= lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measure save frequency");
    base_options.cli |= lyra::opt(ising)["--ising"]("Set J = 1");
    base_options.cli |= lyra::opt(binary)["--binary"]("Sets J =+/-1");
    base_options.cli |= lyra::opt(two_replicas)["-q"]["--two-replicas"]("Simulates two replicas.");
    base_options.cli |= lyra::opt(j_file_path, "J file")["-j"]["--j-file"]("Fiole with link variables");


    auto results = base_options.cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return 1;
    }
    omp_set_num_threads(n_threads);
    if (two_replicas)
        n_replicas = 2;

    spdlog::info("Lx {} Ly {}", base_options.Lx, base_options.Ly);


    using rng_t = std::mt19937_64;

    rng_t rng(base_options.seed);
    std::bernoulli_distribution bern(0.5);
    std::normal_distribution<double> normal(0.0, 1.0);
    using lattice_t = lft::Lattice<uint32_t>;
    lattice_t lat({base_options.Lx, base_options.Ly}, 'C');
    std::array<ea::SpinField<lattice_t> *, 2> replica;

    for (int j = 0; j < n_replicas; ++j) {
        replica[j] = new ea::SpinField<lattice_t>(lat, 1);
    }


    lft::Lattice<uint32_t, 3> j_lat({2, lat.dims[0], lat.dims[1]}, 'C');
    auto j_field = lft::make_field(j_lat, 1.0f);
    if (j_file_path.empty()) {
        if (!ising) {
            if (binary)
                ea::init_bernoulli(j_field, rng);
            else
                ea::init_gaussian(j_field, rng);
        }
    } else {
        std::ifstream ifs(j_file_path, std::ios::in);
        if (!ifs) {
            spdlog::error("Error opening J file : {}", j_file_path);
            exit(1);
        }
        for (std::size_t i = 0; i < j_lat.n_elements; ++i) {
            float val = 7.0;
            if (!(ifs >> val)) {
                spdlog::error("J file  {} too short at {}", j_file_path, i);
                exit(1);
            }
            j_field[i] = val;
        }
    }

    auto j_path = make_file_path(base_options.data_dir, "j", base_options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    ea::HeathBath<float, lattice_t, rng_t> update(base_options.beta, rng, j_field);

    auto start_term = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < base_options.n_term; ++i) {
        for (int j = 0; j < n_replicas; ++j) {
            sweep(*replica[j], update);
        }
    }
    auto end_term = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_term_seconds = end_term - start_term;
    spdlog::info("Thermalisation took {:.3} seconds", elapsed_term_seconds.count());

    auto *em_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "em", base_options.name, "txt"),
        meas_freq > 0, std::fstream::out);
    auto cfg_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "cfg", base_options.name, "bin"),
        base_options.save_freq > 0, std::ios::out | std::ios::binary);


    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < base_options.n_sweeps; ++i) {
        for (int j = 0; j < n_replicas; ++j) {
            sweep(*replica[j], update);
        }
        if (meas_freq > 0 && (i % meas_freq) == 0) {
            if (em_stream_ptr) {
                for (int j = 0; j < n_replicas; ++j) {
                    *em_stream_ptr << ea::energy<double>(*replica[j], j_field) << " ";
                    *em_stream_ptr << ea::magnetisation<double>(*replica[j]) << " ";
                }
                if (two_replicas) {
                    *em_stream_ptr << ea::overlap<double>(*replica[0], *replica[1]) << " ";
                    *em_stream_ptr << ea::link_overlap<double>(*replica[0], *replica[1]) << "\n";
                } else
                    *em_stream_ptr << "\n";
                em_stream_ptr->flush();
            }
        }

        if (base_options.save_freq > 0 && (i % base_options.save_freq) == 0) {
            if (cfg_stream_ptr) {
                for (int j = 0; j < n_replicas; ++j) {
                    replica[j]->write(*cfg_stream_ptr);
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    spdlog::info("Sweeps took {:.3} seconds", elapsed_seconds.count());

    for (int j = 0; j < n_replicas; ++j) {
        delete replica[j];
    }

    if (em_stream_ptr)
        delete em_stream_ptr;
    if (cfg_stream_ptr)
        delete cfg_stream_ptr;

    return 0;
}
