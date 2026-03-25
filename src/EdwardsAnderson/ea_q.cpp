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
#include "ea.h"
#include "parallel_tempering.h"
#include "utils/fs.h"
#include "utils/hardware.h"
#include "options.h"


/**
 * Sets the global log level based on a string name.
 * Handles "trace", "debug", "info", "warn", "err", "critical", "off"
 */
void set_log_level(const std::string& level_name) {
    // from_str returns the level enum; it is case-insensitive by default
    spdlog::level::level_enum level = spdlog::level::from_str(level_name);

    // If the string is invalid (doesn't match any level), from_str returns 'off'
    // but usually, it's better to check if it was a deliberate "off"
    if (level == spdlog::level::off && level_name != "off") {
        spdlog::warn("Unknown log level '{}', defaulting to 'info'", level_name);
        spdlog::set_level(spdlog::level::info);
    }
    else {
        spdlog::set_level(level);
        spdlog::debug("Log level set to {}", level_name);
    }
}

template <typename Field, typename RNG>
void init_field(Field& field, const std::string& file_path, bool binary, bool ising, RNG& rng) {
    if (file_path.empty()) {
        if (!ising) {
            if (binary)
                lft::ea::init_bernoulli(field, rng);
            else
                lft::ea::init_gaussian(field, rng);
        }
    }
    else {
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

void measure_em(std::fstream* em_stream_ptr,
                lft::ea::Replicas<lattice_t>& replica,
                const lft::ea::JField<float, lattice_t>& j_field) {
    if (em_stream_ptr) {
        for (int j = 0; j < replica.q; ++j) {
            *em_stream_ptr << lft::ea::energy<double>(*replica[j], j_field) << " ";
            *em_stream_ptr << lft::ea::magnetisation<double>(*replica[j]) << " ";
        }
        if (replica.q > 1) {
            *em_stream_ptr << lft::ea::overlap<double>(*replica[0], *replica[1]) << " ";
            *em_stream_ptr << lft::ea::link_overlap<double>(*replica[0], *replica[1]) << "\n";
        }
        else
            *em_stream_ptr << "\n";
        em_stream_ptr->flush();
    }
}


int main(int argc, char* argv[]) {
    auto max_threads = std::thread::hardware_concurrency();
    lft::ea::Options options(argc, argv);

    int n_replicas = 1;


    auto options_stream = std::fstream(make_file_path(options.data_dir, "opt", options.name, "yaml"), std::ios::out);
    options_stream << options.emit().c_str() << std::endl;
    options_stream.close();

    set_log_level(options.spdlog_level);
    omp_set_num_threads(options.n_threads);
    spdlog::info("{} threads available, running on {}", max_threads, options.n_threads);


    spdlog::info("Simulating a {}x{} lattice", options.Lx, options.Ly);


    // Random number generator for initialising  fields.
    auto rng = std::mt19937(options.seed);

    // Random number generator for simulations
    lft::rand::taus_array taus_rng(max_threads);
    taus_rng.gen_seeds(options.seed);


    // Creating link variables
    lattice_t lat({options.Lx, options.Ly}, 'C');
    lft::ea::JLattice<lattice_t> j_lat({2, lat.dims[0], lat.dims[1]}, 'C');
    auto j_field = lft::make_field(j_lat, 1.0f);
    init_field(j_field, options.j_file_path, options.binary, options.ising, rng);

    auto j_path = make_file_path(options.data_dir, "j", options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    // Creating replicas

    if (options.two_replicas)
        n_replicas = 2;


    // Creating Parallel tempering updater.

    lft::ea::ParallelTempering<lattice_t> temperer(n_replicas, options.n_betas(),
                                                   options.beta, j_field);
    for (int i = 0; i < options.n_betas(); ++i) {
        for (int j = 0; j < n_replicas; ++j) {
            temperer.replicas[i][j] = new lft::ea::SpinField(lat, 1);
            init_field(*temperer.replicas[i][j], "", true, options.cold_start, rng);
        }
    }


    // Thermalization loop
    auto start_term = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < options.n_term; ++i) {
        if (options.n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);

        if ((i > 0) && (options.exchange_freq > 0) && (i % options.exchange_freq) == 0) {
            for (int j = 0; j < n_replicas; ++j)
                temperer.exchange(j, rng);
        }
    }

    auto end_term = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_term_seconds = end_term - start_term;
    spdlog::info("Thermalization took {:.3} seconds", elapsed_term_seconds.count());

    if (options.exchange_freq > 0) {
        spdlog::info("Exchange acceptance rates:");
        for (auto i = 0; i < options.n_betas() - 1; i++) {
            std::cout << std::format("{:.3f}->{:.3f} {:.2f}", options.beta[i], options.beta[i + 1],
                                     (double)temperer.accepted_v[i] / temperer.exchange_v[i]) <<
                std::endl;
        }
    }
    temperer.reset();

    // Measurements
    std::vector<std::fstream*> em_stream_ptrs(options.n_betas(), nullptr);
    for (int i = 0; i < options.n_betas(); ++i) {
        em_stream_ptrs[i] = optional_fstream_ptr(
            make_file_path(options.data_dir, "em",
                           options.name + std::format("_b{:02d}", i), "txt"),
            options.meas_freq > 0, std::fstream::out);
    }
    auto cfg_stream_ptr = optional_fstream_ptr(
        make_file_path(options.data_dir, "cfg", options.name, "bin"),
        options.save_freq > 0, std::ios::out | std::ios::binary);


    //Main loop
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < options.n_sweeps; ++i) {
        if (options.n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);


        if ((i > 0) && (options.exchange_freq > 0) && (i % options.exchange_freq) == 0) {
            for (int j = 0; j < n_replicas; ++j)
                temperer.exchange(j, rng);
        }


        //Measurements
        if (options.meas_freq > 0 && (i % options.meas_freq) == 0) {
            for (int j = 0; j < options.n_betas(); ++j) {
                measure_em(em_stream_ptrs[j], temperer.replicas[j], j_field);
            }
        }

        //Saving configurations
        if (options.save_freq > 0 && (i % options.save_freq) == 0) {
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
    if (options.exchange_freq > 0) {
        spdlog::info("Exchange acceptance rates:");
        for (auto i = 0; i < options.n_betas() - 1; i++) {
            std::cout << std::format("{:.3f}->{:.3f} {:.2f}", options.beta[i], options.beta[i + 1],
                                     (double)temperer.accepted_v[i] / temperer.exchange_v[i]) <<
                std::endl;
        }
    }
    temperer.reset();

    for (int i = 0; i < options.n_betas(); i++)
        if (em_stream_ptrs[i])
            delete em_stream_ptrs[i];

    if (cfg_stream_ptr)
        delete cfg_stream_ptr;

    return 0;
}
