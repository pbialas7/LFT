//
// Created by pbialas on 20.01.2022.
//

#include <random>
#include <iostream>
#include <fstream>

#include "lyra/lyra.hpp"
#include "spdlog/spdlog.h"

#include "Ising/ising.h"
#include "Ising/correlation_function.h"
#include "MonteCarlo/sweep.h"

int main(int argc, char *argv[]) {


    int Lx = 0, Ly = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta;
    std::string name = "";
    int meas_freq = 10;

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli()
               | lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
               | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
               | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
               | lyra::opt(n_sweeps, "n sweeps ")["-n"]("number of measurment sweeps")
               | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps")
               | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
               | lyra::opt(seed, "seed")["--seed"]("seed")
               | lyra::opt(name, "name")["--name"]("name")
               | lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measurment frquency");


    auto results = cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return (1);
    }
    std::mt19937_64 rng(seed);
    spdlog::default_logger()->set_level(spdlog::level::warn);

    spdlog::debug("Lx {} Ly {} beta {} n-sweeps {} seed {} cold-start {}", Lx, Ly, beta, n_sweeps, seed, cold_start);


    using lattice_t = Lattice<uint32_t>;
    lattice_t lat({Lx, Ly});
    ising::IsingField<lattice_t> ising(lat, 1);
    if (!cold_start)
        ising::hot_start(ising, rng);

    ising::HeathBath<lattice_t> update(beta, rng);

    for (std::size_t i = 0; i < n_term; i++) {
        sweep(ising, update);
    }

    std::fstream energy_mag{std::string("em_") + name + ".txt", energy_mag.out | energy_mag.trunc};
    std::fstream correlations{std::string("cor_") + name + ".bin",
                              correlations.binary | correlations.out | correlations.trunc};
    std::vector<double> cor_function(Lx/2, 0.0);
    for (std::size_t i = 0; i < n_sweeps; i++) {
        sweep(ising, update);
        if ((i + 1) % meas_freq == 0) {
            energy_mag << ising::energy<double>(ising) << " " << ising::magnetisation<double>(ising) << "\n";
            std::fill_n(cor_function.begin(), Lx / 2, 0.0);
            correlation<double>(ising, cor_function);
            correlations.write(reinterpret_cast<char *>(cor_function.data()), cor_function.size() * sizeof(double));
        }
    }


    return 0;
}