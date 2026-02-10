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


    uint32_t Lx = 0, Ly = 0, Lz = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta;
    std::string name = "";
    int meas_freq = 10;
    int cfg_save_freq = 0;

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli()
               | lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
               | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
               | lyra::opt(Lz, "Lz")["-z"]["--Lz"]("Lz")
               | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
               | lyra::opt(n_sweeps, "n sweeps ")["-n"]["--n-sweeps"]("number of measurment sweeps")
               | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps")
               | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
               | lyra::opt(seed, "seed")["--seed"]("seed")
               | lyra::opt(name, "name")["--name"]("name")
               | lyra::opt(cfg_save_freq, "cfg save freq")["--cfg-save-freq"]("configuration save frequency")
               | lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measurment frquency");


    auto results = cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return (1);
    }
    std::mt19937_64 rng(seed);
    spdlog::default_logger()->set_level(spdlog::level::warn);

    spdlog::debug("Lx {} Ly {} Lz beta {} n-sweeps {} seed {} cold-start {}", Lx, Ly, Lz, beta, n_sweeps, seed,
                  cold_start);


    using lattice_t = lft::Lattice<uint32_t, 3>;
    lattice_t lat({Lx, Ly, Lz});
    ising::IsingField<lattice_t> ising(lat, 1);
    if (!cold_start)
        ising::hot_start(ising, rng);

    ising::HeathBath<lattice_t> update(beta, rng);

    for (std::size_t i = 0; i < n_term; i++) {
        sweep(ising, update);
    }

    std::fstream energy_mag{std::string("em3d_") + name + ".txt", energy_mag.out | energy_mag.trunc};


    std::ostream* out(&std::cout);
    if (cfg_save_freq > 0) {
        out = new std::ofstream(std::string("cfg_") + name + ".txt");
    }
    for (std::size_t i = 0; i < n_sweeps; i++) {
        sweep(ising, update);
        if ((i + 1) % meas_freq == 0) {
            energy_mag << ising::E(ising) << " " << ising::M(ising) << "\n";
        }
        if (cfg_save_freq > 0) {
            if ((i + 1) % cfg_save_freq == 0) {
                for (int i = 0; i < ising.n_elements; i++)
                    *out << int(ising[i]) << " ";
                *out << "\n";
            }
        }
    }

    out->flush();
    if (cfg_save_freq > 0) {
        delete out;
    }
    return 0;
}