//
// Created by pbialas on 03.06.25.
//

#pragma once
#include <cmath>
#include <random>

#include "Field/Field.h"


template<typename F, typename RNG=std::mt19937_64>
class Wolff {
public:
    using lattice_t = typename F::lattice_t;
    using field_t = typename F::field_t;
    using rng_t = RNG;

    Wolff(F &field, double beta, rng_t &rng): beta_(beta), lat(field.lat), field_(field), cluster_(field.lat),
                                              rng_(rng), s(0, lat.n_elements - 1), u(0.0, 1.0) {
    }

    const lattice_t &lat;

    size_t fill_and_flip(size_t site_idx) {
        auto p = 1.0 - std::exp(-2.0 * beta_);

        std::vector<size_t> stack;
        stack.reserve(1024);

        size_t cluster_size = 1;
        auto dir = field_[site_idx];
        stack.push_back(site_idx);
        field_[site_idx] *= -1; // flip the site
        while (!stack.empty()) {
            auto site = stack.back();
            stack.pop_back();

            for (auto mu = 0; mu < 2 * lattice_t::DIM; ++mu) {
                auto n_ = lat.nn(site, mu);
                if (field_[n_] == dir) {
                    if (u(rng_) < p) {
                        cluster_size++;
                        stack.push_back(n_);
                        field_[n_] *= -1;
                    }
                }
            }
        }

        return cluster_size;
    }

    double sweep(size_t n) {
        size_t flipped = 0;
        for (size_t i = 0; i < n; ++i) {
            auto site = s(rng_);
            flipped += fill_and_flip(site);
        }
        return double(flipped) / n;
    }

private:
    double beta_;
    F &field_;
    lft::Field<std::int8_t, lattice_t> cluster_;
    rng_t &rng_;
    std::uniform_real_distribution<double> u;
    std::uniform_int_distribution<size_t> s;
};
