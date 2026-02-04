//
// Created by pbialas on 20.01.2022.
//

#pragma once

#include <cmath>
#include <random>


#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace ising {
    template<typename L> using IsingField = Field<int8_t, L>;

    template<typename L, typename RNG=std::mt19937_64>
    struct HeathBath {
        using size_t = typename L::size_t;
        using Field = IsingField<L>;
        using field_t = typename Field::field_t;
        using rng_t = RNG;

        HeathBath(double beta, rng_t &rng) : beta_(beta), rng_(rng), u(0.0, 1.0) {}

        size_t operator()(Field &field, size_t site) {

            auto corona = field.corona(site);
            double p_up = std::exp(beta_ * corona) / (std::exp(beta_ * corona) + std::exp(-beta_ * corona));

            auto r = u(rng_);
            if (r < p_up)
                field[site] = 1;
            else
                field[site] = -1;
            return 1;
        }

        void set_beta(double beta) { beta_ = beta; }

        double beta() const { return beta_; }

    private:
        double beta_;
        rng_t &rng_;
        std::uniform_real_distribution<double> u;
    };


    template<typename L, typename R>
    void hot_start(IsingField<L> &field, R &rng) {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (auto site = 0; site < field.n_elements; ++site) {
            if (u(rng) < 0.5)
                field[site] = 1;
            else
                field[site] = -1;
        }
    }

    template<typename F>
    int64_t M(const F &f) {
        return sum<int64_t>(f);
    }

    template<typename Float, typename F>
    Float magnetisation(const F &f) {
        return mean<Float>(f);
    }

    template<typename Float, typename F>
    Float energy(const F &f) {
        Float e = 0.0;
        for (int i = 0; i < f.n_elements; ++i)
            e += f.up_corona(i) * f[i];
        return -e / f.n_elements;
    }


    template<typename F>
    int64_t E(const F &f) {
        int64_t e = 0;
        for (int i = 0; i < f.n_elements; ++i)
            e += f.up_corona(i) * f[i];
        return -e;
    }


}