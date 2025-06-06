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

int main(int argc, char *argv[]) {
    int Lx = 0, Ly = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta = 0.0;
    std::string name;
    std::string data_dir{"."};
    int meas_freq = 10;
    int save_freq = 0;
    bool wolff = false;
    int n_clusters = 16;

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
               | lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measurment frquency")
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

    using lattice_t = Lattice<uint32_t>;
    lattice_t lat({Lx, Ly});
    ising::IsingField<lattice_t> ising(lat, 1);
    // ReSharper disable once CppDFAConstantConditions
    if (!cold_start)
        ising::hot_start(ising, rng);

    ising::HeathBath<lattice_t> update(beta, rng);
    Wolff<ising::IsingField<lattice_t>, std::mt19937_64> wolff_update(ising, beta, rng);

    // ReSharper disable once CppDFAConstantConditions
    for (std::size_t i = 0; i < n_term; i++) {
        // ReSharper disable once CppDFAUnreachableCode
        if (wolff) {
            auto c = wolff_update.sweep(n_clusters);
        } else {
            auto c = sweep(ising, update);
        }
    }


    auto em_path = make_file_path(data_dir, "em", name, "txt");
    auto corr_path = make_file_path(data_dir, "cor", name, "bin");
    auto cfg_path = make_file_path(data_dir, "cfg", name, "bin");

    std::fstream energy_mag{em_path, std::fstream::out | std::fstream::trunc};
    std::fstream correlations{
        corr_path, std::fstream::binary | std::fstream::out | std::fstream::trunc
    };
    std::fstream configurations{
        cfg_path, std::fstream::binary | std::fstream::out | std::fstream::trunc
    };
    std::vector<double> cor_function(Lx, 0.0);
    double c = 0.0;
    // ReSharper disable once CppDFAConstantConditions
    for (std::size_t i = 0; i < n_sweeps; i++) {
        // ReSharper disable once CppDFAUnreachableCode

        if (wolff) {
            c += wolff_update.sweep(n_clusters);
        } else
            sweep(ising, update);

        if ((meas_freq > 0) && (i + 1) % meas_freq == 0) {
            energy_mag << ising::energy<double>(ising) << " " << ising::magnetisation<double>(ising) << std::endl;
            std::fill_n(cor_function.begin(), Lx, 0.0);
            correlation<double>(ising, cor_function);
            correlations.write(reinterpret_cast<char *>(cor_function.data()), cor_function.size() * sizeof(double));
        }

        if ((save_freq > 0) && (i + 1) % save_freq == 0) {
            configurations.write(reinterpret_cast<const char *>(ising.data()),
                                 ising.n_elements * sizeof(ising::IsingField<lattice_t>::field_t));
        }
    }

    std::cout << c / n_sweeps << "\n";
    return 0;
}
