//
// Created by Piotr Bialas on 30.05.2026.
//

#include "logging.h"

#include <string>

#include <spdlog/spdlog.h>

/**
 * Sets the global log level based on a string name.
 * Handles "trace", "debug", "info", "warn", "err", "critical", "off"
 */
void set_log_level(const std::string& level_name) {
    // from_str returns the level enum; it is case-insensitive by default
    spdlog::level::level_enum level = spdlog::level::from_str(level_name);

    // If the string is invalid (doesn't match any level), from_str returns 'off'
    // but usually, it's better to check if it was a deliberate "off"
    if (level == spdlog::level::off && level_name != "off") {
        spdlog::warn("Unknown log level '{}', defaulting to 'info'", level_name);
        spdlog::set_level(spdlog::level::info);
    }
    else {
        spdlog::set_level(level);
        spdlog::debug("Log level set to {}", level_name);
    }
}
