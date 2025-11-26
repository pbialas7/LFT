//
// Created by pbialas on 25.11.2025.
//

#pragma once


#include <cmath>
#include <random>

#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace ea {
    template<typename L>
    using SpinField = Field<int8_t, L>;
    template<typename L>
    using JLattice = Lattice<typename L::size_t, L::DIM + 1>;
    template<typename L>
    using JField = Field<double, JLattice<L> >;

    template<typename L, typename RNG=std::mt19937_64>
    struct HeathBath {
        using size_t = typename L::size_t;


        using field_class = SpinField<L>;
        using field_t = typename field_class::field_t;
        using rng_t = RNG;


        HeathBath(double beta, rng_t &rng, const JField<L> &j) : beta_(beta), rng_(rng), j_(j), j_lat_(j.lat) {
        }

        size_t operator()(field_class &field, size_t i) {
            double corona = 0.0;
            for (auto d = 0; d < L::DIM; ++d) {
                corona += field[field.lat.up(i, d)] * j_[i + j_lat_.volumes[L::DIM - 1] * d];
            }

            for (auto d = 0; d < L::DIM; ++d) {
                auto dn_site = field.lat.dn(i, d);
                corona += field[dn_site] * j_[dn_site + j_lat_.volumes[L::DIM - 1] * d];
            }
            double p_up = std::exp(beta_ * corona) / (std::exp(beta_ * corona) + std::exp(-beta_ * corona));

            auto r = u_(rng_);
            if (r < p_up)
                field[i] = 1;
            else
                field[i] = -1;
            return 1;
        }

    private:
        double beta_;
        rng_t &rng_;
        std::uniform_real_distribution<double> u_;
        JField<L> j_;
        JLattice<L> j_lat_;
    };
}
