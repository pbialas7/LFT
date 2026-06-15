//
// Created by pbialas on 25.11.2025.
//

#include <fstream>
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;
#include <thread>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include<omp.h>

#include <Field/Lattice.h>

#include "utils/rand.h"
#include "utils/fs.h"
#include "utils/logging.h"

#include "MonteCarlo/sweep.h"
#include "ea.h"
#include "parallel_tempering_mt.h"
#include "options_cli11.h"
#include "sg_file.h"


inline std::string format_hms(double total_seconds) {
    if (total_seconds < 0) total_seconds = 0;
    long long s = static_cast<long long>(total_seconds + 0.5); // round
    long long h = s / 3600;
    s %= 3600;
    long long m = s / 60;
    s %= 60;

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << ":"
            << std::setw(2) << m << ":" << std::setw(2) << s;
    return oss.str();
}

inline void print_progress_bar(
    std::int64_t current_iteration, // 1-based preferred
    std::int64_t final_iterations,
    const std::chrono::high_resolution_clock::time_point &start_time,
    int bar_width = 40) {
    if (final_iterations <= 0) return;

    current_iteration = std::clamp<std::int64_t>(current_iteration, 0, final_iterations);
    const double progress = static_cast<double>(current_iteration) /
                            static_cast<double>(final_iterations);

    const int filled = static_cast<int>(progress * bar_width);

    const auto now = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(now - start_time).count();

    double eta = 0.0;
    if (current_iteration > 0 && current_iteration < final_iterations) {
        const double sec_per_iter = elapsed / static_cast<double>(current_iteration);
        eta = sec_per_iter * static_cast<double>(final_iterations - current_iteration);
    }

    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
        std::cout << (i < filled ? '#' : '-');
    }
    std::cout << "] "
            << std::setw(3) << static_cast<int>(progress * 100.0) << "% "
            << "elapsed " << format_hms(elapsed) << " "
            << "eta " << format_hms(eta)
            << std::flush;

    if (current_iteration == final_iterations) {
        std::cout << "\n";
    }
}

/**
 * Initializes the field with either random values or from a file. If file_path is empty, it initializes randomly.
 * If ising is true, it initializes with 1. If binary is true, it initializes with +/-1.
 * Otherwise, it initializes with Gaussian random numbers.
 */
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
        spdlog::info("Reading link variables from {}", file_path);
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
                lft::ea::Replicas<lattice_t> &replicas,
                const lft::ea::JField<float, lattice_t> &j_field) {
    if (em_stream_ptr) {
        for (int j = 0; j < replicas.n_replicas; ++j) {
            *em_stream_ptr << lft::ea::energy<double>(*replicas[j], j_field) << " ";
            *em_stream_ptr << lft::ea::magnetisation<double>(*replicas[j]) << " ";
        }
        if (replicas.n_replicas > 1) {
            *em_stream_ptr << lft::ea::overlap<double>(*replicas[0], *replicas[1]) << " ";
            *em_stream_ptr << lft::ea::link_overlap<double>(*replicas[0], *replicas[1]) << "\n";
        } else
            *em_stream_ptr << "\n";
        em_stream_ptr->flush();
    }
}


int main(int argc, char *argv[]) {
    auto max_threads = std::thread::hardware_concurrency();

    lft::ea::Options options;
    options.parse(argc, argv);
    spdlog::debug("data_dir after parse: '{}'", options.data_dir);

    if (options.beta.empty()) {
        spdlog::error("No betas specified");
        return 1;
    }

    int n_replicas = options.n_replicas;

    spdlog::info("Writing to data dir = {}", options.data_dir);
    auto options_stream = std::fstream(
        make_file_path(options.data_dir, "opt", options.Lx, options.Ly, options.name, "yaml"),
        std::ios::out
    );
    options_stream << options.emit().c_str() << std::endl;
    options_stream.close();

    auto options_toml_stream = std::fstream(
        make_file_path(options.data_dir, "opt", options.Lx, options.Ly, options.name, "toml"),
        std::ios::out
    );
    options_toml_stream << options.app.config_to_str(CLI::ConfigOutputMode::AllDefaults) << std::endl;
    options_toml_stream.close();

    set_log_level(options.spdlog_level);
    omp_set_num_threads(options.n_threads);
    spdlog::info("{} threads available, using {}", max_threads, options.n_threads);

    spdlog::info("Simulating a {}x{} lattice", options.Lx, options.Ly);

    // Random number generator for initializing fields.
    auto rng = std::mt19937_64(options.seed);
    auto j_rng = std::mt19937_64(options.j_seed);

    // Random number generator for simulations
    lft::rand::taus_array taus_rng(max_threads);
    taus_rng.gen_seeds(options.seed);


    // Creating link variables
    lattice_t lat({options.Lx, options.Ly}, 'C');
    lft::ea::JLattice<lattice_t> j_lat({2, lat.dims[0], lat.dims[1]}, 'C');
    auto j_field = lft::make_field(j_lat, 1.0f);
    init_field(j_field, options.j_file_path, options.binary, options.ising, j_rng);

    // Writing out the J link variables
    auto j_path = make_file_path(options.data_dir, "j", options.Lx, options.Ly, options.name, "txt");
    std::fstream j_file(j_path, std::fstream::out);
    j_file << j_field << "\n";
    j_file.close();

    //Creating chains and replicas

    spdlog::info("Using {} betas {} with replicas", options.n_betas(), n_replicas);

    // Creating Parallel tempering updater.
    // All chains and replicas share same link variables.
    lft::ea::ParallelTemperingMT<lattice_t> temperer(n_replicas, options.n_betas(),
                                                     options.beta, j_field);

    // Created and initializes each replica either randomly with +/-1 (hot start) or 1 (cold start).
    for (int i = 0; i < options.n_betas(); ++i) {
        for (int j = 0; j < n_replicas; ++j) {
            temperer.chains[i][j] = new lft::ea::SpinField(lat, 1);
            init_field(*temperer.chains[i][j], "", true, options.cold_start, rng);
        }
    }


    /*********************
    * Thermalization loop
    *********************/
    spdlog::info("Starting the thermalization loop with {} sweeps and exchange frequency {}", options.n_term,
                 options.exchange_freq);
    auto start_term = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < options.n_term; ++i) {
        if (options.n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);

        if ((i > 0) && (options.exchange_freq > 0) && (i % options.exchange_freq) == 0) {
            if (options.n_threads > 1)
                temperer.exchange_mt(taus_rng);
            else
                temperer.exchange(taus_rng[0]);
        }

        if (options.cluster_freq > 0 && (i > 0) && (i % options.cluster_freq) == 0) {
            if (options.n_threads > 1)
                temperer.houdayer_mt(taus_rng);
            else
                temperer.houdayer(taus_rng[0]);
        }
        print_progress_bar(i + 1, options.n_term, start_term, 40);
    } // End thermalization loop

    auto end_term = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_term_seconds = end_term - start_term;
    spdlog::info("Thermalization took {:.3} seconds", elapsed_term_seconds.count());
    if (options.exchange_freq > 0 && options.n_term > 1) {
        temperer.print_stats(std::cout);
    }
    temperer.reset_stats();

    // Creating output files

    // Energy and magnetization
    std::vector<std::fstream *> em_stream_ptrs(options.n_betas(), nullptr);
    for (int i = 0; i < options.n_betas(); ++i) {
        em_stream_ptrs[i] = optional_fstream_ptr(
            make_file_path(options.data_dir, "em", options.Lx, options.Ly,
                           options.name + std::format("_b{:03d}", i), "txt"),
            options.meas_freq > 0, std::fstream::out);
    }

    // Configuration files
    std::vector<std::fstream *> cfg_stream_ptrs(options.n_betas(), nullptr);
    for (int i = 0; i < options.n_betas(); ++i) {
        cfg_stream_ptrs[i] = optional_fstream_ptr(
            make_file_path(options.data_dir, "cfg", options.Lx, options.Ly, options.name + std::format("_b{:03d}", i),
                           "bin"),
            options.save_freq > 0, std::ios::out | std::ios::binary);
    }


    /***********************
     * Main loop
     ***********************/
    spdlog::info("Starting main loop with {} sweeps and exchange frequency {}", options.n_sweeps,
                 options.exchange_freq);
    spdlog::info("measure frequency {} configurations save frequency {}", options.meas_freq, options.save_freq);
    auto start_main_all = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sweep_time(0.0);
    std::chrono::duration<double> exchange_time(0.0);
    std::chrono::duration<double> cluster_time(0.0);

    for (int i = 0; i < options.n_sweeps; ++i) {
        auto start_main_sweep = std::chrono::high_resolution_clock::now();
        if (options.n_threads > 1)
            temperer.sweep_mt(1, taus_rng);
        else
            temperer.sweep_t1(1, taus_rng[0]);
        auto end_main_sweep = std::chrono::high_resolution_clock::now();
        sweep_time += std::chrono::duration<double>(end_main_sweep - start_main_sweep);

        auto start_main_exchange = std::chrono::high_resolution_clock::now();
        if ((i > 0) && (options.exchange_freq > 0) && (i % options.exchange_freq) == 0) {
            if (options.n_threads > 1)
                temperer.exchange_mt(taus_rng);
            else
                temperer.exchange(taus_rng[0]);
        }
        auto end_main_exchange = std::chrono::high_resolution_clock::now();
        exchange_time += std::chrono::duration<double>(end_main_exchange - start_main_exchange);


        auto start_main_cluster = std::chrono::high_resolution_clock::now();
        if (options.cluster_freq > 0 && (i > 0) && (i % options.cluster_freq) == 0) {
            if (options.n_threads > 1)
                temperer.houdayer_mt(taus_rng);
            else
                temperer.houdayer(taus_rng[0]);
        }
        auto end_main_cluster = std::chrono::high_resolution_clock::now();
        cluster_time += std::chrono::duration<double>(end_main_cluster - start_main_cluster);


        //Measurements
        if (options.meas_freq > 0 && (i % options.meas_freq) == 0) {
            for (int j = 0; j < options.n_betas(); ++j) {
                measure_em(em_stream_ptrs[j], temperer.chains[j], j_field);
            }
        }

        //Saving configurations
        if (options.save_freq > 0 && (i % options.save_freq) == 0) {
            for (int i = 0; i < options.n_betas(); ++i)
                if (cfg_stream_ptrs[i]) {
                    for (int j = 0; j < n_replicas; ++j) {
                        temperer.chains[i][j]->write(*cfg_stream_ptrs[i]);
                    }
                    cfg_stream_ptrs[i]->flush();
                }
        }
        print_progress_bar(i + 1, options.n_sweeps, start_main_all, 40);
    } // End main loop

    auto end_main_all = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_main_all - start_main_all;
    spdlog::info("Main loop  took {:.3} seconds", elapsed_seconds.count());
    spdlog::info("Sweeps loop  took {:.3} seconds", sweep_time.count());
    spdlog::info("Exchange took {:.3} seconds", exchange_time.count());
    spdlog::info("Clusters took {:.3} seconds", cluster_time.count());

    if (options.exchange_freq > 0 && options.n_sweeps > 0) {
        temperer.print_stats(std::cout);
        auto exch_file = make_file_path(options.data_dir, "exch", options.Lx, options.Ly, options.name, "txt");
        spdlog::info("Exchange file {}", exch_file.string());
        std::fstream exch_stream(exch_file, std::fstream::out);
        temperer.print_stats(exch_stream);
        exch_stream.close();
    }
    temperer.reset_stats();

    // Close files (delete the pointers)
    for (int i = 0; i < options.n_betas(); i++)
        if (em_stream_ptrs[i])
            delete em_stream_ptrs[i];

    for (int i = 0; i < options.n_betas(); i++)
        if (cfg_stream_ptrs[i])
            delete cfg_stream_ptrs[i];

    return 0;
}
