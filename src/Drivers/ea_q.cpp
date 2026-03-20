//
// Created by pbialas on 25.11.2025.
//

#include <fstream>
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;
#include <thread>


#include<omp.h>

#include <Field/Lattice.h>

#include "utils/rand.h"
#include "MonteCarlo/sweep.h"
#include "EdwardsAnderson/ea.h"
#include "EdwardsAnderson/parallel_tempering.h"
#include "utils/fs.h"
#include "utils/hardware.h"
#include "ising_base_options.h"


/**
 * Sets the global log level based on a string name.
 * Handles "trace", "debug", "info", "warn", "err", "critical", "off"
 */
void set_log_level(const std::string &level_name) {
    // from_str returns the level enum; it is case-insensitive by default
    spdlog::level::level_enum level = spdlog::level::from_str(level_name);

    // If the string is invalid (doesn't match any level), from_str returns 'off'
    // but usually, it's better to check if it was a deliberate "off"
    if (level == spdlog::level::off && level_name != "off") {
        spdlog::warn("Unknown log level '{}', defaulting to 'info'", level_name);
        spdlog::set_level(spdlog::level::info);
    } else {
        spdlog::set_level(level);
        spdlog::debug("Log level set to {}", level_name);
    }
}

template<typename Field, typename RNG>
void init_field(Field &field, const std::string &file_path, bool binary, bool ising, RNG &rng) {
    if (file_path.empty()) {
        if (!ising) {
            if (binary)
                lft::ea::init_bernoulli(field, rng);
            else
                lft::ea::init_gaussian(field, rng);
        }
    } else {
        std::ifstream ifs(file_path, std::ios::in);
        if (!ifs) {
            spdlog::error("Error opening file : {}", file_path);
            exit(1);
        }
        for (std::size_t i = 0; i < field.lat.n_elements; ++i) {
            float val = 7.0;
            if (!(ifs >> val)) {
                spdlog::error("File  {} too short at {}", file_path, i);
                exit(1);
            }
            field[i] = val;
        }
    }
}

using lattice_t = lft::Lattice<uint32_t>;

void measure_em(std::fstream *em_stream_ptr,
                lft::ea::Replicas<lattice_t> &replica,
                const lft::ea::JField<float, lattice_t> &j_field) {
    if (em_stream_ptr) {
        for (int j = 0; j < replica.q; ++j) {
            *em_stream_ptr << lft::ea::energy<double>(*replica[j], j_field) << " ";
            *em_stream_ptr << lft::ea::magnetisation<double>(*replica[j]) << " ";
        }
        if (replica.q > 1) {
            *em_stream_ptr << lft::ea::overlap<double>(*replica[0], *replica[1]) << " ";
            *em_stream_ptr << lft::ea::link_overlap<double>(*replica[0], *replica[1]) << "\n";
        } else
            *em_stream_ptr << "\n";
        em_stream_ptr->flush();
    }
}


int main(int argc, char *argv[]) {
    auto max_threads = std::thread::hardware_concurrency();
    IsingBaseOptions base_options;

    int meas_freq = 0;
    bool ising = false;
    bool binary = false;
    bool two_replicas = false;
    std::string j_file_path;
    int n_replicas = 1;
    std::string spdlog_level("info");
    int n_threads = 1;


    base_options.cli |= lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measure save frequency");
    base_options.cli |= lyra::opt(ising)["--ising"]("Set J = 1");
    base_options.cli |= lyra::opt(binary)["--binary"]("Sets J =+/-1");
    base_options.cli |= lyra::opt(two_replicas)["-q"]["--two-replicas"]("Simulates two replicas.");
    base_options.cli |= lyra::opt(j_file_path, "J file")["-j"]["--j-file"]("File with link variables");
    base_options.cli |= lyra::opt(spdlog_level, "spdlog level")["--level"]("Sets the spdlog level");
    base_options.cli |= lyra::opt(n_threads, "N threads")["--n-threads"]("Set number of threads to use");


    auto results = base_options.cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return 1;
    }

    set_log_level(spdlog_level);
    omp_set_num_threads(n_threads);
    spdlog::info("{} threads available, running on {}", max_threads, n_threads);


    spdlog::info("Simulating a {}x{} lattice", base_options.Lx, base_options.Ly);


    // Random number generator for initialising  fields.
    auto rng = std::mt19937(base_options.seed);

    // Random number generator for simulations
    lft::rand::taus_array taus_rng(max_threads);
    taus_rng.gen_seeds(base_options.seed);


    // Creating link variables
    lattice_t lat({base_options.Lx, base_options.Ly}, 'C');
    lft::ea::JLattice<lattice_t> j_lat({2, lat.dims[0], lat.dims[1]}, 'C');
    auto j_field = lft::make_field(j_lat, 1.0f);
    init_field(j_field, j_file_path, binary, ising, rng);

    auto j_path = make_file_path(base_options.data_dir, "j", base_options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    // Creating replicas

    if (two_replicas)
        n_replicas = 2;


    // Creating Parallel tempering updater.
    lft::ea::ParallelTempering<lattice_t> temperer(n_replicas, 1, {0.9f}, j_field);
    for (int j = 0; j < n_replicas; ++j) {
        temperer.replicas[0][j] = new lft::ea::SpinField(lat, 1);
        init_field(*temperer.replicas[0][j], "", true, base_options.cold_start, rng);
    }


    // Thermalisation loop
    auto start_term = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < base_options.n_term; ++i) {
        if (n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);
    }
    auto end_term = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_term_seconds = end_term - start_term;
    spdlog::info("Thermalisation took {:.3} seconds", elapsed_term_seconds.count());


    // Measurements
    auto *em_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "em", base_options.name, "txt"),
        meas_freq > 0, std::fstream::out);
    auto cfg_stream_ptr = otional_fstream_ptr(
        make_file_path(base_options.data_dir, "cfg", base_options.name, "bin"),
        base_options.save_freq > 0, std::ios::out | std::ios::binary);


    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < base_options.n_sweeps; ++i) {
        if (n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);

        if (meas_freq > 0 && (i % meas_freq) == 0) {
            measure_em(em_stream_ptr, temperer.replicas[0], j_field);
        }

        if (base_options.save_freq > 0 && (i % base_options.save_freq) == 0) {
            if (cfg_stream_ptr) {
                for (int j = 0; j < n_replicas; ++j) {
                    temperer.replicas[0][j]->write(*cfg_stream_ptr);
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    spdlog::info("Sweeps took {:.3} seconds", elapsed_seconds.count());


    if (em_stream_ptr)
        delete em_stream_ptr;
    if (cfg_stream_ptr)
        delete cfg_stream_ptr;

    return 0;
}
