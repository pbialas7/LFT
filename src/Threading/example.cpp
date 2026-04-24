#include <iostream>
#include <chrono>
#include <vector>
#include "thread_pool.hpp"

// A CPU-bound task that sums integers from 1..n
long long sum_to(long long n) {
    long long s = 0;
    for (long long i = 1; i <= n; ++i) s += i;
    return s;
}

int main() {
    // --- Basic usage ---
    ThreadPool pool(4);   // 4 worker threads
    std::cout << "Pool size: " << pool.size() << " threads\n";

    // Submit a lambda
    auto f1 = pool.submit([] {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        return std::string("hello from lambda");
    });

    // Submit a free function with arguments
    auto f2 = pool.submit(sum_to, 1'000'000LL);

    std::cout << "f1 = " << f1.get() << "\n";
    std::cout << "f2 = " << f2.get() << "\n";

    // --- Batch of tasks ---
    constexpr int NUM_TASKS = 20;
    std::vector<std::future<long long>> futures;
    futures.reserve(NUM_TASKS);

    for (int i = 1; i <= NUM_TASKS; ++i)
        futures.push_back(pool.submit(sum_to, i * 100'000LL));

    long long total = 0;
    for (auto& f : futures)
        total += f.get();

    std::cout << "Batch total = " << total << "\n";

    // --- Exception propagation ---
    auto f3 = pool.submit([] -> int {
        throw std::runtime_error("task error");
        return 0;
    });

    try {
        f3.get();
    } catch (const std::exception& e) {
        std::cout << "Caught: " << e.what() << "\n";
    }

    // Pool destructor: drains the queue then joins all threads
    return 0;
}

// Compile with:
//   g++ -std=c++17 -O2 -pthread example.cpp -o example
