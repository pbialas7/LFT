// 
// Created by pbialas on 20.06.2026.
//

#pragma once

#include <chrono>
#include <cstdint>
#include <string>

// Formats a duration in seconds as HH:MM:SS.
std::string format_hms(double total_seconds);

void print_progress_bar(
    std::int64_t current_iteration, // 1-based preferred
    std::int64_t final_iterations,
    const std::chrono::high_resolution_clock::time_point& start_time,
    int bar_width = 40);
