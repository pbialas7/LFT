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
#include <CLI/CLI.hpp>

namespace lft::ea {
    class Options {
        std::vector<float> split_floats(const std::string &s, char delim = ',') {
            auto trim = [](std::string &x) {
                while (!x.empty() && std::isspace(static_cast<unsigned char>(x.front()))) x.erase(x.begin());
                while (!x.empty() && std::isspace(static_cast<unsigned char>(x.back()))) x.pop_back();
            };

            std::vector<float> result;
            std::stringstream ss(s);
            std::string token;

            while (std::getline(ss, token, delim)) {
                trim(token);
                if (token.empty()) continue;
                result.push_back(std::stof(token));
            }

            return result;
        }

    public:
        CLI::App app;
        bool debug = false;

        Options();

        int parse(int argc, char *argv[]) {
            try {
                app.parse(argc, argv); // instead of CLI11_PARSE macro inside a non-main function
            } catch (const CLI::ParseError &e) {
                return app.exit(e);
            }

            if (debug) spdlog_level = "debug";


            return 0;
        }

        int n_betas() const { return beta.size(); }

        std::string raw_betas;
        std::vector<float> beta;
        uint32_t Lx = 0;
        uint32_t Ly = 0;
        uint32_t Lz = 0;
        int save_freq = 0;
        int n_term = 0;
        int n_sweeps = 0;
        bool cold_start = false;
        std::string name{"name"};
        std::string data_dir{"."};
        size_t seed = std::mt19937_64::default_seed;
        size_t j_seed = std::mt19937_64::default_seed;
        int meas_freq = 0;
        bool ising = false;
        bool binary = false;
        std::string j_file_path;
        int n_replicas = 1;
        std::string spdlog_level{"info"};
        int n_threads = 1;
        int exchange_freq = 0;
        int cluster_freq = 0;

        YAML::Emitter &emit();

    private:
        YAML::Emitter yaml;
    };
}
