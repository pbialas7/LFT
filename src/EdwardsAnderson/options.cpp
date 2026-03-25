//
// Created by Piotr Bialas on 21.03.2026.
//

#include "options.h"

#define pair(name) YAML::Key<<#name<<YAML::Value<<name

YAML::Emitter& lft::ea::Options::emit() {
    yaml << YAML::BeginMap;
    yaml << YAML::Key << "beta" << YAML::Value << YAML::Flow << beta;
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
    yaml << pair(no_pt);
    yaml << YAML::EndMap;
    return yaml;
}
