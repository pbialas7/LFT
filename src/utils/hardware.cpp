//
// Created by pbialas on 20.03.2026.
//

#include "hardware.h"

#include <fstream>
#include <string>

namespace lft::hdw {
    // Returns frequency in MHz, or -1.0 on failure
    double get_core_mhz(int core_id) {
        std::string path = "/sys/devices/system/cpu/cpu" +
                           std::to_string(core_id) +
                           "/cpufreq/scaling_cur_freq";
        std::ifstream file(path);
        if (!file.is_open()) return -1.0;

        long kHz;
        if (file >> kHz) {
            return kHz / 1000.0; // Convert kHz to MHz
        }
        return -1.0;
    }
}