//
// Created by pbialas on 20.01.2022.
//

#include <random>

#include "lyra/lyra.hpp"
#include "spdlog/spdlog.h"

#include "Ising/ising.h"
#include "MonteCarlo/sweep.h"

int main(int argc, char *argv[]) {


    int Lx = 0, Ly = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta;

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli()
               | lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
               | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
               | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
               | lyra::opt(n_sweeps, "n sweeps ")["-n"]("number of measurment sweeps")
               | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps")
               | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
               | lyra::opt(seed, "seed")["--seed"]("seed");


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
    std::cout.precision(12);
    std::cout << ising::M(ising) << "\n";
    return 0;
}