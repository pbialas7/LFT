//
// Created by Piotr Bialas on 21.03.2026.
//

#pragma once

#include <random>

#include <yaml-cpp/yaml.h>

#include "lyra/lyra.hpp"

#define pair(name) YAML::Key<<#name<<YAML::Value<<name

namespace lft::ea {
    class Options {
    public:
        Options() {
            spdlog_level = "info";
            cli |= lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
                | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
                | lyra::opt(beta, "beta")["-b"]["--beta"]("inverse temperature")
                | lyra::opt(n_sweeps, "n sweeps ")["-n"]["--n-sweeps"]("number of measurement sweeps")
                | lyra::opt(n_term, "n term ")["-t"]["--n-term"]("number of normalisation sweeps")
                | lyra::opt(cold_start)["-c"]["--cold-start"]("cold start")
                | lyra::opt(seed, "seed")["--seed"]("seed")
                | lyra::opt(name, "name")["--name"]("name")
                | lyra::opt(data_dir, "data_dir")["--data-dir"]("data directory")
                | lyra::opt(save_freq, "save frequency")["--save-freq"]("configuration save frequency")
                | lyra::opt(meas_freq, "measure frequency")["--meas-freq"]("measure save frequency")
                | lyra::opt(ising)["--ising"]("Set J = 1")
                | lyra::opt(binary)["--binary"]("Sets J =+/-1")
                | lyra::opt(two_replicas)["-q"]["--two-replicas"]("Simulates two replicas.")
                | lyra::opt(j_file_path, "J file")["-j"]["--j-file"]("File with link variables")
                | lyra::opt(spdlog_level, "spdlog level")["--level"]("Sets the spdlog level")
                | lyra::opt(n_threads, "N threads")["--n-threads"]("Set number of threads to use")
                | lyra::opt(exchange_freq, "exchange freq")["--exchange-freq"]("Replicas exchange frequency")
                | lyra::opt(spdlog_level, "spdlog level")["--debug"]("Loging level");
        }

        lyra::cli cli;

        double beta = 0.0;
        uint32_t Lx = 0;
        uint32_t Ly = 0;
        int save_freq = 0;
        int measure_freq = 0;
        int n_term = 0;
        int n_sweeps = 0;
        bool cold_start = false;
        std::string name{"name"};
        std::string data_dir{"."};
        size_t seed = std::mt19937_64::default_seed;
        int meas_freq = 0;
        bool ising = false;
        bool binary = false;
        bool two_replicas = false;
        std::string j_file_path;
        int n_replicas = 1;
        std::string spdlog_level;
        int n_threads = 1;
        int exchange_freq = 0;

        YAML::Emitter& emit() {
            yaml << YAML::BeginMap;
            yaml << pair(beta);
            yaml << pair(Lx);
            yaml << pair(Ly);
            yaml << pair(save_freq);
            yaml << pair(measure_freq);
            yaml << pair(n_term);
            yaml << pair(n_sweeps);
            yaml << pair(cold_start);
            yaml << pair(seed);
            yaml << pair(name);
            yaml << pair(data_dir);
            yaml << pair(spdlog_level);
            yaml << pair(n_threads);
            yaml << pair(exchange_freq);
            yaml << pair(ising);
            yaml << pair(binary);
            yaml << pair(two_replicas);
            yaml << pair(j_file_path);
            return yaml;
        }

    private:
        YAML::Emitter yaml;
    };
}
