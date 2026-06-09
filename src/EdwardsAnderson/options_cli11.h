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
        // std::vector<float> split_floats(const std::string &s, char delim = ',') {
        //     std::vector<float> result;
        //     std::stringstream ss(s);
        //     std::string token;
        //     while (std::getline(ss, token, delim))
        //         result.push_back(std::stof(token));
        //     return result;
        // }


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

        Options() {
            //app.add_option("-h,--help", "Print this help documentation");
            app.add_option("--Lx", Lx, "Lx");
            app.add_option("--Ly", Ly, "Ly");
            app.add_option("-n, --n-sweeps", n_sweeps, "number of measurement sweeps");
            app.add_option("-t,--n-term", n_term, "number of thermalization sweeps");
            app.add_option("-c,--cold-start", cold_start, "cold start");
            app.add_option("--seed", seed, "seed");
            app.add_option("--name", name, "name");
            app.add_option("--data-dir", data_dir, "data directory");
            app.add_option("--save-freq", save_freq, "configuration save frequency");
            app.add_option("--meas-freq", meas_freq, "measure save frequency");
            app.add_flag("--ising", ising, "Set J = 1");
            app.add_flag("--binary", binary, "Sets J = +/-1");
            app.add_option("-q,--n-replicas", n_replicas, "Simulates n replicas.");
            app.add_option("-j,--j-file", j_file_path, "File with link variables");
            app.add_option("--level", spdlog_level, "Sets the spdlog level");
            app.add_option("--n-threads", n_threads, "Set number of threads to use");
            app.add_option("--exchange-freq", exchange_freq, "Replicas exchange frequency");
            app.add_option("--beta", raw_betas, "Comma-separated list of betas for parallel tempering");


            // class member
            app.add_flag("--debug", debug, "Enable debug logging");

            app.set_config("--config", "", "Path to optional config file");
        }

        int parse(int argc, char *argv[]) {
            try {
                app.parse(argc, argv); // instead of CLI11_PARSE macro inside a non-main function
            } catch (const CLI::ParseError &e) {
                return app.exit(e);
            }

            if (debug) spdlog_level = "debug";

            try {
                beta = split_floats(raw_betas);
            } catch (const std::exception &e) {
                spdlog::error("Invalid --beta '{}': {}", raw_betas, e.what());
                return 1;
            }

            if (beta.empty()) {
                spdlog::error("No betas specified in '{}'.", raw_betas);
                return 1;
            }

            return 0;
        }

        int n_betas() const { return beta.size(); }

        std::string raw_betas;
        std::vector<float> beta;
        uint32_t Lx = 0;
        uint32_t Ly = 0;
        int save_freq = 0;
        int n_term = 0;
        int n_sweeps = 0;
        bool cold_start = false;
        std::string name{"name"};
        std::string data_dir{"."};
        size_t seed = std::mt19937_64::default_seed;
        int meas_freq = 0;
        bool ising = false;
        bool binary = false;
        std::string j_file_path;
        int n_replicas = 1;
        std::string spdlog_level{"info"};
        int n_threads = 1;
        int exchange_freq = 0;

        YAML::Emitter &emit();

    private:
        YAML::Emitter yaml;
    };
}
