//
// Created by pbialas on 20.01.2022.
//

#include <random>
#include <iostream>
#include <fstream>

#include "lyra/lyra.hpp"
#include "spdlog/spdlog.h"

#include "Ising/ising.h"
#include "MonteCarlo/sweep.h"

int main(int argc, char* argv[]) {
    spdlog::default_logger()->set_level(spdlog::level::info);
    uint32_t Lx = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double beta;
    int cfg_save_freq = 0;
    std::string cfg_file_name;

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli()
        | lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
        | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
        | lyra::opt(n_sweeps, "n sweeps ")["-n"]("number of measurment sweeps")
        | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps")
        | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
        | lyra::opt(seed, "seed")["--seed"]("seed")
        | lyra::opt(cfg_save_freq, "cfg save freq")["--cfg-save-freq"]("configuration save frequency")
        | lyra::opt(cfg_file_name, "cfg file name")["-o"]["--cfg-output"]("configurations file name");


    auto results = cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return (1);
    }
    std::mt19937_64 rng(seed);


    spdlog::info("Lx {}  beta {} n-sweeps {} seed {} cold-start {}", Lx, beta, n_sweeps, seed, cold_start);


    using lattice_t = lft::Lattice<uint32_t, 1>;
    lattice_t lat({Lx});
    ising::IsingField<lattice_t> ising(lat, 1);
    if (!cold_start)
        ising::hot_start(ising, rng);

    ising::HeathBath<lattice_t> update(beta, rng);
    spdlog::info("Starting termalisation");
    for (std::size_t i = 0; i < n_term; i++) {
        sweep(ising, update);
    }
    spdlog::info("Finished termalisation");

    spdlog::info("Opening cfg file");
    std::ofstream ofs(cfg_file_name, std::ios::trunc);

    spdlog::info("Starting ising");
    for (std::size_t i = 0; i < n_sweeps; i++) {
        sweep(ising, update);
        if (cfg_save_freq > 0) {
            if ((i + 1) % cfg_save_freq == 0) {
                for (int i = 0; i < ising.n_elements; i++)
                    ofs << int(ising[i]) << " ";
                ofs << "\n";
            }
        }
    }
    spdlog::info("Finished ising");
    ofs.flush();

    spdlog::info("Closed cfg file");
    return 0;
}
