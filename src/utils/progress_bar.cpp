//
// Created by Piotr Bialas on 20.06.2026.
//

#include "progress_bar.h"

#include <iomanip>
#include <sstream>
#include <iostream>

std::string format_hms(double total_seconds) {
    if (total_seconds < 0) total_seconds = 0;
    long long s = static_cast<long long>(total_seconds + 0.5); // round
    long long h = s / 3600;
    s %= 3600;
    long long m = s / 60;
    s %= 60;

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << ":"
        << std::setw(2) << m << ":" << std::setw(2) << s;
    return oss.str();
}


void print_progress_bar(
    std::int64_t current_iteration, // 1-based preferred
    std::int64_t final_iterations,
    const std::chrono::high_resolution_clock::time_point& start_time,
    int bar_width) {
    if (final_iterations <= 0) return;

    current_iteration = std::clamp<std::int64_t>(current_iteration, 0, final_iterations);
    const double progress = static_cast<double>(current_iteration) /
        static_cast<double>(final_iterations);

    const int filled = static_cast<int>(progress * bar_width);

    const auto now = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(now - start_time).count();

    double eta = 0.0;
    if (current_iteration > 0 && current_iteration < final_iterations) {
        const double sec_per_iter = elapsed / static_cast<double>(current_iteration);
        eta = sec_per_iter * static_cast<double>(final_iterations - current_iteration);
    }

    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
        std::cout << (i < filled ? '#' : '-');
    }
    std::cout << "] "
        << std::setw(3) << static_cast<int>(progress * 100.0) << "% "
        << "elapsed " << format_hms(elapsed) << " "
        << "eta " << format_hms(eta)
        << std::flush;

    if (current_iteration == final_iterations) {
        std::cout << "\n";
    }
}
