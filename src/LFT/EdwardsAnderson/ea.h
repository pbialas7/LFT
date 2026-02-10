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
    using SpinField = lft::Field<int8_t, L>;
    template<typename L>
    using JLattice = lft::Lattice<typename L::index_t, L::DIM + 1>;
    template<typename L>
    using JField = lft::Field<double, JLattice<L> >;

    template<typename L, typename RNG=std::mt19937_64>
    struct HeathBath {
        using size_t = L::size_t;


        using field_class = SpinField<L>;
        using field_t = field_class::field_t;
        using rng_t = RNG;


        HeathBath(double beta, rng_t &rng, const JField<L> &j) : beta_(beta), rng_(rng), j_(j), j_lat_(j.lat) {
        }

        size_t operator()(field_class &field, size_t i) {
            double corona = 0.0;
            for (auto d = 0; d < L::DIM; ++d) {
                corona += field[field.lat.up(i, d)] * j_[i + j_lat_.n_elements * d];
            }

            for (auto d = 0; d < L::DIM; ++d) {
                auto dn_site = field.lat.dn(i, d);
                corona += field[dn_site] * j_[dn_site + j_lat_.n_elements * d];
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

    template<typename Float, typename F>
    Float magnetisation(const F &f) {
        return mean<Float>(f);
    }

    template<typename Float, typename F, typename JF>
    Float energy(const F &field, const JF &j_) {
        Float e = 0.0;
        for (int i = 0; i < field.n_elements; ++i) {
            double corona = 0.0;
            for (auto d = 0; d < F::DIM; ++d) {
                auto i_up = field.lat.up(i, d);
                auto j_up = j_[i + j_.lat.n_elements * d];
                auto s_up = field[i_up];
                corona += s_up * j_up;
            }
            e += field[i] * corona;
        }

        return -e;
    }
}
