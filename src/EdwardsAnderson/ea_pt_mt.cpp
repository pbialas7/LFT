//
// Parallel tempering for the Edwards-Anderson model with replicas and Houdayer cluster moves.
//
// The same source is built twice: `ea_pt_mt` (EA_DIM=2) and `ea_pt_mt_3D` (EA_DIM=3).
//

#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <omp.h>

#include <Field/Lattice.h>

#include "utils/rand.h"
#include "utils/fs.h"
#include "utils/logging.h"
#include "utils/progress_bar.h"

#include "ea.h"
#include "ea_io.h"
#include "parallel_tempering_mt.h"
#include "options_cli11.h"

#ifndef EA_DIM
#error "EA_DIM (2 or 3) must be defined"
#endif

namespace {
    using Duration = std::chrono::duration<double>;

    struct StepTimes {
        Duration sweep{0};
        Duration exchange{0};
        Duration cluster{0};
    };

    template <int DIM>
    std::array<std::uint32_t, DIM> lattice_dims(const lft::ea::Options& options) {
        static_assert(DIM == 2 || DIM == 3);
        if constexpr (DIM == 2)
            return {options.Lx, options.Ly};
        else
            return {options.Lx, options.Ly, options.Lz};
    }

    template <std::size_t N>
    std::string dims_string(const std::array<std::uint32_t, N>& dims) {
        std::string result;
        for (std::size_t d = 0; d < N; ++d)
            result += (d ? "x" : "") + std::to_string(dims[d]);
        return result;
    }

    // Runs op() and adds its duration to *slot, if given.
    template <typename Op>
    void timed(Duration* slot, Op&& op) {
        if (!slot) {
            op();
            return;
        }
        auto start = std::chrono::high_resolution_clock::now();
        op();
        *slot += std::chrono::high_resolution_clock::now() - start;
    }

    /**
     * One Monte Carlo step i: a sweep of all chains, then, every exchange_freq / cluster_freq steps, replica
     * exchange and Houdayer cluster moves. rngs holds one random stream per chain. If times is not null, the time
     * spent in each part is added to it.
     */
    template <typename PT, typename RNGS>
    void mc_step(PT& temperer, RNGS& rngs, const lft::ea::Options& options, int i, StepTimes* times) {
        timed(times ? &times->sweep : nullptr, [&] { temperer.sweep(1, rngs); });

        if ((i > 0) && (options.exchange_freq > 0) && (i % options.exchange_freq) == 0)
            timed(times ? &times->exchange : nullptr, [&] { temperer.exchange(rngs); });

        if (options.cluster_freq > 0 && (i > 0) && (i % options.cluster_freq) == 0)
            timed(times ? &times->cluster : nullptr, [&] { temperer.houdayer(rngs); });
    }

    // Name of the output files of the chain with index i_beta.
    std::string beta_name(const lft::ea::Options& options, int i_beta) {
        return options.name + std::format("_b{:03d}", i_beta);
    }

    // An output file: its path (for error messages) and the stream, which is null if the file is not written.
    struct OutFile {
        fs::path path;
        std::unique_ptr<std::fstream> stream;

        // Throws if a write to the stream has failed. Does nothing for a file that is not written.
        void check() const {
            if (stream)
                check_stream(*stream, path);
        }

        void close() {
            if (stream)
                close_checked(*stream, path);
        }
    };

    // One output file per beta, opened (and checked) only if `create` is true.
    template <std::size_t N>
    std::vector<OutFile> open_beta_files(const lft::ea::Options& options, const std::array<std::uint32_t, N>& dims,
                                         const std::string& prefix, const std::string& ext, bool create,
                                         std::ios::openmode mode) {
        std::vector<OutFile> files;
        for (int i = 0; i < options.n_betas(); ++i) {
            auto path = make_file_path(options.data_dir, prefix, dims, beta_name(options, i), ext);
            auto stream = optional_fstream_ptr(path, create, mode);
            files.push_back({std::move(path), std::move(stream)});
        }
        return files;
    }

    // Creates the data directory if needed; throws if the path exists but is not a directory.
    void ensure_data_dir(const fs::path& data_dir) {
        if (!fs::exists(data_dir)) {
            fs::create_directories(data_dir);
            spdlog::info("Created data dir {}", data_dir.string());
        }
        else if (!fs::is_directory(data_dir)) {
            throw std::runtime_error(std::format("--data-dir {} exists but is not a directory", data_dir.string()));
        }
    }

    // Writes `text` to a new file at `path`, checking the open, the write and the close.
    void write_text_file(const fs::path& path, const std::string& text) {
        auto stream = open_checked(path, std::ios::out);
        stream << text;
        close_checked(stream, path);
    }

    // NaN/inf check on the bit pattern: -ffast-math (-ffinite-math-only) lets the compiler fold std::isfinite to true.
    inline bool is_finite(float x) {
        return (std::bit_cast<std::uint32_t>(x) & 0x7f800000u) != 0x7f800000u;
    }

    /*
     * Checks the options before anything is written. Logs every problem found and returns false if the run must
     * not start. Inconsistent but harmless combinations only produce a warning.
     */
    template <int DIM>
    bool validate(const lft::ea::Options& options) {
        bool ok = true;
        auto fail = [&](const std::string& msg) {
            spdlog::error("{}", msg);
            ok = false;
        };

        const auto dims = lattice_dims<DIM>(options);
        for (auto d : dims) {
            if (d == 0 || d % 2 != 0) {
                fail(std::format("Invalid lattice size {}: every extent must be positive and even "
                                 "(the checkerboard update needs even extents)", dims_string(dims)));
                break;
            }
        }
        if (DIM == 2 && options.Lz != 0)
            fail(std::format("--Lz {} given to the 2D driver; use ea_pt_mt_3D for 3D lattices", options.Lz));

        if (options.beta.empty())
            fail("No betas specified");
        for (std::size_t i = 0; i < options.beta.size(); ++i) {
            auto b = options.beta[i];
            if (!is_finite(b) || b < 0)
                fail(std::format("beta[{}] = {} must be finite and non-negative", i, b));
            if (i > 0 && !(options.beta[i - 1] < b))
                fail(std::format("betas must be strictly increasing, but beta[{}] = {} >= beta[{}] = {}",
                                 i - 1, options.beta[i - 1], i, b));
        }

        if (options.n_replicas < 1)
            fail(std::format("--n-replicas must be at least 1, got {}", options.n_replicas));
        if (options.n_threads < 1)
            fail(std::format("--n-threads must be at least 1, got {}", options.n_threads));

        const std::pair<const char*, int> non_negative[] = {
            {"--n-sweeps", options.n_sweeps}, {"--n-term", options.n_term},
            {"--meas-freq", options.meas_freq}, {"--save-freq", options.save_freq},
            {"--exchange-freq", options.exchange_freq}, {"--cluster-freq", options.cluster_freq}};
        for (const auto& [name, value] : non_negative)
            if (value < 0)
                fail(std::format("{} must not be negative, got {}", name, value));

        if (!options.j_file_path.empty()) {
            if (!std::ifstream(options.j_file_path))
                fail(std::format("Cannot read the coupling file --j-file {}", options.j_file_path));
            if (options.ising || options.binary)
                spdlog::warn("Couplings are read from {}; --ising/--binary are ignored", options.j_file_path);
        }
        else if (options.ising && options.binary) {
            spdlog::warn("Both --ising and --binary given; using --ising (J = 1)");
        }

        return ok;
    }

    template <int DIM>
    int run(lft::ea::Options& options) {
        using lattice_t = lft::Lattice<uint32_t, DIM>;
        auto max_threads = std::thread::hardware_concurrency();

        set_log_level(options.spdlog_level);
        if (!validate<DIM>(options))
            return 1;

        const auto dims = lattice_dims<DIM>(options);
        int n_replicas = options.n_replicas;

        spdlog::info("Writing to data dir = {}", options.data_dir);
        ensure_data_dir(options.data_dir);
        write_text_file(make_file_path(options.data_dir, "opt", dims, options.name, "yaml"),
                        std::string(options.emit().c_str()) + "\n");
        write_text_file(make_file_path(options.data_dir, "opt", dims, options.name, "toml"),
                        options.app.config_to_str(CLI::ConfigOutputMode::AllDefaults) + "\n");

        omp_set_num_threads(options.n_threads);
        spdlog::info("{} threads available, using {}", max_threads, options.n_threads);
        if (options.n_threads > options.n_betas())
            spdlog::info("Only {} chains: at most {} threads are used", options.n_betas(), options.n_betas());

        spdlog::info("Simulating a {} lattice", dims_string(dims));

        // Random number generator for initializing fields.
        auto rng = std::mt19937_64(options.seed);
        auto j_rng = std::mt19937_64(options.j_seed);

        // Random number generators for the simulation: one stream per beta chain, so that the results do not
        // depend on the number of threads.
        lft::rand::taus_array taus_rng(options.n_betas());
        taus_rng.gen_seeds(options.seed);


        // Creating link variables
        lattice_t lat(dims, 'C');
        lat.require_checkerboard();
        typename lft::ea::JLattice<lattice_t>::dim_t j_dims;
        j_dims[0] = DIM;
        for (int d = 0; d < DIM; ++d)
            j_dims[d + 1] = lat.dims[d];
        lft::ea::JLattice<lattice_t> j_lat(j_dims, 'C');
        auto j_field = lft::make_field(j_lat, 1.0f);
        lft::ea::init_j_field(j_field, options.j_file_path, options.binary, options.ising, j_rng);

        // Writing out the J link variables
        auto j_path = make_file_path(options.data_dir, "j", dims, options.name, "txt");
        auto j_file = open_checked(j_path, std::ios::out);
        j_file << j_field << "\n";
        close_checked(j_file, j_path);

        // Output files for the main loop: energy and magnetization, and configurations, one file per beta.
        // They are opened before thermalization so that a bad path fails before any simulation time is spent.
        auto em_files = open_beta_files(options, dims, "em", "txt", options.meas_freq > 0, std::ios::out);
        auto cfg_files = open_beta_files(options, dims, "cfg", "bin", options.save_freq > 0,
                                         std::ios::out | std::ios::binary);

        //Creating chains and replicas

        spdlog::info("Using {} betas with {} replicas", options.n_betas(), n_replicas);

        // Creating Parallel tempering updater.
        // All chains and replicas share same link variables.
        lft::ea::ParallelTemperingMT<lattice_t> temperer(n_replicas, options.n_betas(),
                                                         options.beta, j_field);

        // Created and initializes each replica either randomly with +/-1 (hot start) or 1 (cold start).
        for (int i = 0; i < options.n_betas(); ++i) {
            for (int j = 0; j < n_replicas; ++j) {
                temperer.chains[i][j] = new lft::ea::SpinField(lat, 1);
                lft::ea::init_spin_field(*temperer.chains[i][j], options.cold_start, rng);
            }
        }


        /*********************
        * Thermalization loop
        *********************/
        spdlog::info("Starting the thermalization loop with {} sweeps, exchange frequency {} and cluster frequency {}",
                     options.n_term, options.exchange_freq, options.cluster_freq);
        auto start_term = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < options.n_term; ++i) {
            mc_step(temperer, taus_rng, options, i, nullptr);
            print_progress_bar(i + 1, options.n_term, start_term, 40);
        } // End thermalization loop

        auto end_term = std::chrono::high_resolution_clock::now();
        Duration elapsed_term_seconds = end_term - start_term;
        spdlog::info("Thermalization took {:.3} seconds", elapsed_term_seconds.count());
        if (options.exchange_freq > 0 && options.n_term > 1) {
            temperer.print_stats(std::cout);
        }
        temperer.reset_stats();


        /***********************
         * Main loop
         ***********************/
        spdlog::info("Starting main loop with {} sweeps and exchange frequency {}", options.n_sweeps,
                     options.exchange_freq);
        spdlog::info("measure frequency {} configurations save frequency {}", options.meas_freq, options.save_freq);
        auto start_main_all = std::chrono::high_resolution_clock::now();
        StepTimes times;

        for (int i = 0; i < options.n_sweeps; ++i) {
            mc_step(temperer, taus_rng, options, i, &times);

            //Measurements
            if (options.meas_freq > 0 && (i % options.meas_freq) == 0) {
                for (int j = 0; j < options.n_betas(); ++j) {
                    lft::ea::measure_em(em_files[j].stream.get(), temperer.chains[j], j_field);
                    em_files[j].check();
                }
            }

            //Saving configurations
            if (options.save_freq > 0 && (i % options.save_freq) == 0) {
                for (int i_beta = 0; i_beta < options.n_betas(); ++i_beta)
                    if (auto& cfg = cfg_files[i_beta]; cfg.stream) {
                        for (int j = 0; j < n_replicas; ++j) {
                            temperer.chains[i_beta][j]->write(*cfg.stream);
                        }
                        cfg.stream->flush();
                        cfg.check();
                    }
            }
            print_progress_bar(i + 1, options.n_sweeps, start_main_all, 40);
        } // End main loop

        auto end_main_all = std::chrono::high_resolution_clock::now();
        Duration elapsed_seconds = end_main_all - start_main_all;
        spdlog::info("Main loop  took {:.3} seconds", elapsed_seconds.count());
        spdlog::info("Sweeps loop  took {:.3} seconds", times.sweep.count());
        spdlog::info("Exchange took {:.3} seconds", times.exchange.count());
        spdlog::info("Clusters took {:.3} seconds", times.cluster.count());

        for (auto& f : em_files) f.close();
        for (auto& f : cfg_files) f.close();

        if (options.exchange_freq > 0 && options.n_sweeps > 0) {
            temperer.print_stats(std::cout);
            auto exch_file = make_file_path(options.data_dir, "exch", dims, options.name, "txt");
            spdlog::info("Exchange file {}", exch_file.string());
            auto exch_stream = open_checked(exch_file, std::ios::out);
            temperer.print_stats(exch_stream);
            close_checked(exch_stream, exch_file);
        }
        temperer.reset_stats();

        return 0;
    }
}

int main(int argc, char* argv[]) {
    lft::ea::Options options;
    // Stop after --help (exit code 0) or a parse error (non-zero), before anything is written.
    if (auto exit_code = options.parse(argc, argv))
        return *exit_code;
    spdlog::debug("data_dir after parse: '{}'", options.data_dir);

    try {
        return run<EA_DIM>(options);
    }
    catch (const std::exception& e) {
        spdlog::error("{}", e.what());
        return 1;
    }
}
