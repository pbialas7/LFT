//
// Created by pbialas on 29.04.22.
//
#include <random>
#include <iostream>

int main() {

    std::mt19937_64 rng(770);
    std::normal_distribution<float> norm{0, 1.0};

    size_t n = 100000;
    double mean = 0.0;
    for (size_t i = 0; i < 1000000; i++) {
        mean += norm(rng);
    }

    std::cout << mean / n << "\n";

    return 0;
}