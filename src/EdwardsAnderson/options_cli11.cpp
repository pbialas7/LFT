//
// Created by Piotr Bialas on 21.03.2026.
//

#include "options_cli11.h"

#define pair(name) YAML::Key<<#name<<YAML::Value<<name

lft::ea::Options::Options() {
    //app.add_option("-h,--help", "Print this help documentation");
    app.add_option("--Lx", Lx, "Lx")->capture_default_str();
    app.add_option("--Ly", Ly, "Ly")->capture_default_str();
    app.add_option("--Lz", Lz, "Lz")->capture_default_str();
    app.add_option("-n, --n-sweeps", n_sweeps, "number of measurement sweeps")->capture_default_str();
    app.add_option("-t,--n-term", n_term, "number of thermalization sweeps")->capture_default_str();
    app.add_option("-c,--cold-start", cold_start, "cold start")->capture_default_str();
    app.add_option("--seed", seed, "seed")->capture_default_str();
    app.add_option("--j-seed", j_seed, "j seed")->capture_default_str();
    app.add_option("--name", name, "name")->capture_default_str();
    app.add_option("--data-dir", data_dir, "data directory")->capture_default_str();
    app.add_option("--save-freq", save_freq, "configuration save frequency")->capture_default_str();
    app.add_option("--meas-freq", meas_freq, "measure save frequency")->capture_default_str();
    app.add_flag("--ising", ising, "Set J = 1")->capture_default_str();
    app.add_flag("--binary", binary, "Sets J = +/-1")->capture_default_str();
    app.add_option("-q,--n-replicas", n_replicas, "Simulates n replicas.")->capture_default_str();
    app.add_option("-j,--j-file", j_file_path, "File with link variables")->capture_default_str();
    app.add_option("--level", spdlog_level, "Sets the spdlog level")->capture_default_str();
    app.add_option("--n-threads", n_threads, "Set number of threads to use")->capture_default_str();
    app.add_option("--exchange-freq", exchange_freq, "Replicas exchange frequency")->capture_default_str();
    app.add_option("--cluster-freq", cluster_freq, "Replicas cluster frequency")->capture_default_str();
    app.add_option("--beta", beta, "list of betas for parallel tempering")->capture_default_str();

    // class member
    app.add_flag("--debug", debug, "Enable debug logging")->capture_default_str();

    app.set_config("--config", "", "Path to optional config file");
}

YAML::Emitter& lft::ea::Options::emit() {
    yaml << YAML::BeginMap;
    yaml << YAML::Key << "beta" << YAML::Value << YAML::Flow << beta;
    yaml << pair(Lx);
    yaml << pair(Ly);
    yaml << pair(Lz);
    yaml << pair(save_freq);
    yaml << pair(meas_freq);
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
    yaml << pair(n_replicas);
    yaml << pair(j_file_path);
    yaml << pair(cluster_freq);
    yaml << YAML::EndMap;
    return yaml;
}
