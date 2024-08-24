//
// Created by pbialas on 20.01.2022.
//

#include <random>
#include <fstream>
#include <chrono>
#include <cmath>

#include "lyra/lyra.hpp"


#include "Phi4/phi4.h"
#include "MonteCarlo/sweep.h"

#include "spdlog/spdlog.h"
#include "spdlog/fmt/chrono.h"
#include "spdlog/sinks/stdout_color_sinks.h"

int main(int argc, char *argv[]) {
    auto console = spdlog::stdout_color_mt("console");
    console->set_level(spdlog::level::debug);


    int Lx = 0, Ly = 0;
    std::size_t n_sweeps = 0, n_term = 0;
    bool cold_start = false;
    double kappa = 1.0;
    double M2 = 0.0;
    double lambda = 1.0;
    double eps = 0.5;
    int save_freq = 0;
    std::string out_file_name = "o.bin";
    int n_hits = 4;


    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    auto cli = lyra::cli() |
               lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx") |
               lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly") |
               lyra::opt(kappa, "kappa")["-k"]["--kappa"]("") |
               lyra::opt(M2, "M2")["--M2"]("e") |
               lyra::opt(lambda, "lambda")["-l"]["--lambda"]("") |
               lyra::opt(eps, "eps")["-e"]["--eps"]("") |
               lyra::opt(n_sweeps, "n sweeps ")["-n"]["--n-sweeps"]("number of measurment sweeps") |
               lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of termalisation sweeps") |
               lyra::opt(cold_start)["-c"]["--cold-start"]("cold start") |
               lyra::opt(save_freq, "save freq")["--save-freq"]("configuration save frequency") |
               lyra::opt(seed, "seed")["--seed"]("") |
               lyra::opt(out_file_name, "output")["-o"]["--output"]("congigurations output file") |
               lyra::opt(n_hits, "n hits")["--n-hits"]("Number of hits in multi-hit metropolis (=4)");


    auto results = cli.parse({argc, argv});
    if (!results) {
        std::cerr << "Error in command line: " << results.message() << std::endl;
        return (1);
    }
    std::mt19937_64 rng(seed);
    spdlog::default_logger()->set_level(spdlog::level::warn);

    console->debug("Lx {} Ly {} M2 {} n-sweeps {} seed {} output {} cold-start {}", Lx, Ly, M2, n_sweeps, seed,
                   out_file_name, cold_start);


    using lattice_t = Lattice<uint32_t>;
    lattice_t lat({Lx, Ly});

    using phi_t = phi4::ScalarField<lattice_t>;
    phi_t phi(lat, 1.0);
    if (!cold_start)
        phi4::hot_start(phi, rng);

    phi4::Metropolis<lattice_t> update(kappa, M2, lambda, n_hits, rng);
    update.set_eps(eps);

    auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < n_term; i++) {
        sweep(phi, update);
    }

    std::ofstream out_cfg;
    if (n_sweeps > 0 && save_freq > 0) {
        out_cfg.open(out_file_name.c_str(), std::ios::binary);
    }
    double acceptance = 0.0;
    console->debug("n sweeps {}", n_sweeps);
    for (std::size_t i = 0; i < n_sweeps; i++) {
        auto accepted = sweep(phi, update);
        auto accept = (double) accepted / phi.n_elements;
        acceptance += accept;
        if (save_freq > 0 && ((i + 1) % save_freq) == 0)
            out_cfg.write((const char *) phi.data(), phi.n_elements * sizeof(phi_t::field_t));
    }
    out_cfg.close();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    fmt::print("elapsed time {}\n", elapsed_seconds);

    console->debug("accepted {}", acceptance);
    std::cout << "acceptance = " << acceptance / n_sweeps << "\n";
}