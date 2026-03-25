//
// Created by Piotr Bialas on 21.03.2026.
//

#pragma once

#include <random>
#include <vector>
#include <string>
#include <sstream>

#include <spdlog/spdlog.h>

#include <yaml-cpp/yaml.h>

#include "lyra/lyra.hpp"


namespace lft::ea {
    class Options {
        std::vector<float> split_floats(const std::string& s, char delim = ',') {
            std::vector<float> result;
            std::stringstream ss(s);
            std::string token;
            while (std::getline(ss, token, delim))
                result.push_back(std::stof(token));
            return result;
        }

        std::string raw_betas;

    public:
        Options(int argc, char* argv[]) {
            spdlog_level = "info";
            cli |= lyra::opt(Lx, "Lx")["-x"]["--Lx"]("Lx")
                | lyra::opt(Ly, "Ly")["-y"]["--Ly"]("Ly")
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
                | lyra::opt(spdlog_level, "spdlog level")["--debug"]("Loging level")
                | lyra::opt(raw_betas, "beta")["--beta"]("Comma-separated list of betas for parallel tempering");

            auto results = cli.parse({argc, argv});
            if (!results) {
                spdlog::error("Error in command line: {}.", results.message());
                exit(1);
            }

            beta = split_floats(raw_betas);
            if (beta.empty()) {
                spdlog::error("No betas specified.");
                exit(0);
            }
        }

        int n_betas() const { return beta.size(); }

        lyra::cli cli;


        std::vector<float> beta;
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

        YAML::Emitter& emit();

    private:
        YAML::Emitter yaml;
    };
}
