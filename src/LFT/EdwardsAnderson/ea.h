//
// Created by pbialas on 25.11.2025.
//

#pragma once


#include <cmath>
#include <random>

#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace ea {
    template <typename L>
    using SpinField = lft::Field<int8_t, L>;
    template <typename L>
    using JLattice = lft::Lattice<typename L::index_t, L::DIM + 1>;
    template <typename F, typename L>
    using JField = lft::Field<F, JLattice<L>>;


    template <typename F, typename RNG>
    void init_bernoulli(F& j_field, RNG& rng) {
        std::bernoulli_distribution bern(0.5);
        for (int i = 0; i < j_field.n_elements; ++i) {
            if (bern(rng)) {
                j_field[i] = 1;
            }
            else {
                j_field[i] = -1;
            }
        }
    }

    template <typename F, typename RNG>
    void init_gaussian(F& j_field, RNG& rng) {
        std::normal_distribution<float> normal(0.0f, 1.0f);
        for (int i = 0; i < j_field.n_elements; ++i) {
            j_field[i] = normal(rng);
        }
    }


    template <typename F, typename L, typename RNG=std::mt19937_64>
    struct HeathBath {
        using size_t = L::size_t;
        using index_t = L::index_t;


        using field_class = SpinField<L>;
        using field_t = field_class::field_t;
        using rng_t = RNG;


        HeathBath(double beta, rng_t& rng, const JField<F, L>& j) : beta_(beta), rng_(rng), j_(j), j_lat_(j.lat) {
        }

        size_t operator()(field_class& field, size_t i) {
            double corona = 0.0;
            for (auto d = 0; d < L::DIM; ++d) {
                auto idx = i + field.lat.n_elements * d;
                corona += field[field.lat.up(i, d)] * j_[idx];
            }

            for (auto d = 0; d < L::DIM; ++d) {
                auto dn_site = field.lat.dn(i, d);
                auto idx = dn_site + field.lat.n_elements * d;
                corona += field[dn_site] * j_[idx];
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
        rng_t& rng_;
        std::uniform_real_distribution<double> u_;
        JField<F, L> j_;
        JLattice<L> j_lat_;
    };

    template <typename Float, typename F>
    Float magnetisation(const F& f) {
        return mean<Float>(f);
    }

    template <typename Float, typename F, typename JF>
    Float energy(const F& field, const JF& j_) {
        Float e = 0.0;
        for (int i = 0; i < field.n_elements; ++i) {
            double corona = 0.0;
            for (auto d = 0; d < F::DIM; ++d) {
                auto i_up = field.lat.up(i, d);
                auto j_up = j_[i + field.lat.n_elements * d];
                auto s_up = field[i_up];
                corona += s_up * j_up;
            }
            e += field[i] * corona;
        }

        return -e;
    }


    template <typename Float, typename F, typename JF>
    Float energy_dn(const F& field, const JF& j_) {
        Float e = 0.0;
        for (int i = 0; i < field.n_elements; ++i) {
            double corona = 0.0;
            for (auto d = 0; d < F::DIM; ++d) {
                auto i_dn = field.lat.dn(i, d);
                auto j_dn = j_[i_dn + field.lat.n_elements * d];
                auto s_dn = field[i_dn];
                corona += s_dn * j_dn;
            }
            e += field[i] * corona;
        }

        return -e;
    }

    template <typename Float, typename F>
    Float overlap(const F& f1, const F& f2) {
        Float overlap_ = 0.0;
        for (int i = 0; i < f1.n_elements; ++i) {
            overlap_ += f1[i] * f2[i];
        }
        return overlap_ / f1.n_elements;
    }

    template <typename Float, typename F>
    Float link_overlap(const F& field1, const F& field2) {
        Float link_overlap_ = 0.0;
        for (int i = 0; i < field1.n_elements; ++i) {
            for (auto d = 0; d < F::DIM; ++d) {
                auto i_up = field1.lat.up(i, d);
                auto s_up1 = field1[i_up];
                auto s_up2 = field2[i_up];
                link_overlap_ += field1[i] * s_up1 * field2[i] * s_up2;
            }
        }

        return -link_overlap_ / (field1.n_elements * F::DIM);
    }
}
