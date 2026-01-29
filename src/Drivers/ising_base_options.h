//
// Created by Piotr Bialas on 26/11/2025.
//

#pragma once

#include <random>

#include "lyra/lyra.hpp"

class IsingBaseOptions {
public:
    IsingBaseOptions() {
        cli |= lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
                | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
                | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
                | lyra::opt(n_sweeps, "n sweeps ")["-n"]["--n-sweeps"]("number of measurement sweeps")
                | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of normalisation sweeps")
                | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
                | lyra::opt(seed, "seed")["--seed"]("seed")
                | lyra::opt(name, "name")["--name"]("name")
                | lyra::opt(data_dir, "data_dir")["--data-dir"]("data directory")
                | lyra::opt(save_freq, "save frequency")["--save-freq"]("configuration save frequency");
    }


    double beta = 0.0;
    int Lx = 0;
    int Ly = 0;
    int save_freq = 0;
    int measure_freq = 0;
    int n_term = 0;
    int n_sweeps = 0;
    bool cold_start = false;
    std::string name{"name"};
    std::string data_dir{"."};

    std::mt19937_64::result_type seed = std::mt19937_64::default_seed;

    lyra::cli cli;
};
