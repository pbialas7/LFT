//
// Created by pbialas on 20.01.2022.
//

#include <random>
#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

#include "lyra/lyra.hpp"
#include "spdlog/spdlog.h"

#include "Ising/ising.h"
#include "Ising/correlation_function.h"
#include "MonteCarlo/sweep.h"
#include "Ising/Wolff.h"
#include "utils/fs.h"
#include "utils/extract.h"

int main(int argc, char* argv[]) {
    uint32_t Lx = 0, Ly = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta = 0.0;
    std::string name;
    std::string data_dir{"."};
    int meas_freq = 10;
    int corr_freq = 10;
    int save_freq = 0;
    bool wolff = false;
    int n_clusters = 16;
    int tune_wolff = 0;

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli()
        | lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
        | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
        | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
        | lyra::opt(n_sweeps, "n sweeps ")["-n"]["--n-sweeps"]("number of measurment sweeps")
        | lyra::opt(n_clusters, "n clusters")["-k"]["--n-clusters"](
            "number of clusters in a sweep for Wolff algorithm")
        | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps")
        | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
        | lyra::opt(wolff)["-w"]["--wolff"]("use Wolff algorithm")

        | lyra::opt(seed, "seed")["--seed"]("seed")
        | lyra::opt(name, "name")["--name"]("name")
        | lyra::opt(data_dir, "data_dir")["--data-dir"]("data directory")
        | lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measurment frequency")
        | lyra::opt(corr_freq, "correlation measure frequency")["--corr-freq"]("correlation measure frequency")
        | lyra::opt(tune_wolff, "tune steps")["--tune-wolff"]("number of tuning steps for wolf algorithm")
        | lyra::opt(save_freq, "save frequency")["--save-freq"]("configuration save frequency");


    auto results = cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return 1;
    }
    std::mt19937_64 rng(seed);
    spdlog::default_logger()->set_level(spdlog::level::warn);

    spdlog::debug("Lx {} Ly {} beta {} n-sweeps {} seed {} cold-start {}", Lx, Ly, beta, n_sweeps, seed, cold_start);

    fs::path data_path{data_dir};
    if (!fs::exists(data_path)) {
        std::cerr << "Data directory does `" << data_path << "' not exist, creating." << std::endl;
        fs::create_directories(data_path);
    }

    using lattice_t = lft::Lattice<uint32_t>;
    lattice_t lat({Lx, Ly});
    ising::IsingField<lattice_t> ising(lat, 1);
    // ReSharper disable once CppDFAConstantConditions
    if (!cold_start)
        ising::hot_start(ising, rng);


    // ReSharper disable once CppDFAUnusedValue
    ising::HeathBath<lattice_t> update(beta, rng);
    Wolff wolff_update(ising, beta, rng);

    // ReSharper disable once CppDFAConstantConditions
    double c = 0.0;
    for (std::size_t i = 0; i < n_term; i++) {
        // ReSharper disable once CppDFAUnreachableCode
        if (wolff) {
            c += wolff_update.sweep(n_clusters);
        }
        else {
            c += sweep(ising, update);
        }
    }

    if (n_term > 0)
        std::cout << c / n_term << " after " << n_term << "\n";

    c = 0.0;

    // ReSharper disable once CppDFAUnreachableCode
    if (wolff && tune_wolff > 0) {
        // ReSharper disable once CppDFAUnreachableCode
        for (std::size_t i = 0; i < tune_wolff; i++) {
            c += wolff_update.sweep(n_clusters);
        }

        c /= tune_wolff;
        n_clusters = int(ceil(Lx * Ly / c));

        std::cerr << "tuning " << c << " after " << tune_wolff << " n clusters " << n_clusters << "\n";
    }

    auto em_path = make_file_path(data_dir, "em", name, "txt");
    auto corr_path = make_file_path(data_dir, "cor", name, "bin");
    auto cfg_path = make_file_path(data_dir, "cfg", name, "bin");
    auto edges_path = make_file_path(data_dir, "edges", name, "bin");

    std::fstream energy_mag{em_path, std::fstream::out | std::fstream::trunc};
    std::fstream correlations{
        corr_path, std::fstream::binary | std::fstream::out | std::fstream::trunc
    };
    std::fstream configurations{
        cfg_path, std::fstream::binary | std::fstream::out | std::fstream::trunc
    };

    std::fstream edges_stream{
        edges_path, std::fstream::binary | std::fstream::out | std::fstream::trunc
    };

    std::vector<double> cor_function(Lx, 0.0);
    std::vector<decltype(ising)::field_t> edges(Lx + Ly - 1);
    // ReSharper disable once CppDFAConstantConditions
    c = 0.0;
    for (std::size_t i = 0; i < n_sweeps; i++) {
        // ReSharper disable once CppDFAUnreachableCode

        double c_size = 0.0;
        if (wolff) {
            c_size += wolff_update.sweep(n_clusters);
            c += c_size;
        }
        else
            c += sweep(ising, update);

        if ((meas_freq > 0) && (i + 1) % meas_freq == 0) {
            energy_mag << ising::energy<double>(ising) << " " << ising::magnetisation<double>(ising)
                << " " << c_size << std::endl;
        }
        if ((corr_freq > 0) && (i + 1) % corr_freq == 0) {
            std::fill_n(cor_function.begin(), Lx, 0.0);
            correlation<double>(ising, cor_function);
            correlations.write(reinterpret_cast<char*>(cor_function.data()), cor_function.size() * sizeof(double));
        }


        if ((save_freq > 0) && (i + 1) % save_freq == 0) {
            configurations.write(reinterpret_cast<const char*>(ising.data()),
                                 ising.n_elements * sizeof(ising::IsingField<lattice_t>::field_t));
            extract_edges(ising, edges);
            edges_stream.write(reinterpret_cast<const char*>(edges.data()),
                               edges.size() * sizeof(decltype(ising)::field_t));
        }
    }

    if (n_sweeps > 0)
        std::cout << c / n_sweeps << "\n";

    return 0;
}
