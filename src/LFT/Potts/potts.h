//
// Created by pbialas on 20.01.2022.
//

#pragma once

#include <cmath>
#include <random>
#include <cstdint>


#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace potts {

    template<uint8_t N, typename L>
    class PottsField : public Field<uint8_t, L> {
    public:
        static const uint8_t N_STATES = N;
    };

    template<typename L> using IsingField = Field<int8_t, L>;

    template<uint8_t N, typename L, typename RNG=std::mt19937_64>
    struct Metropolis {
        using size_t = typename L::size_t;
        using Field = PottsField<N, L>;
        using field_t = typename Field::field_t;
        using rng_t = RNG;

        Metropolis(double beta, rng_t &rng) : beta_(beta), rng_(rng), u_(0, N - 1), r_(0, 1) {}

        size_t operator()(Field &field, size_t site) {
            size_t accepted = 0;
            int energy = 0;
            field_t at_site = field[site];
            for (auto d = 0; d < field.DIM; ++d) {
                energy += at_site == field[field.lat.up(site, d)] ? 1 : -1;
            }

            for (auto d = 0; d < field.DIM; ++d) {
                energy += at_site == field[field.lat.dn(site, d)] ? 1 : -1;
            }
            const size_t n_hits = 8;
            for (size_t i = 0; i < n_hits; ++i) {
                field_t prop = (at_site + u(rng_)) % N;
                field_t prop_energy;
                for (auto d = 0; d < field.DIM; ++d) {
                    prop_energy += prop == field[field.lat.up(site, d)] ? 1 : -1;
                }

                for (auto d = 0; d < field.DIM; ++d) {
                    prop_energy += prop == field[field.lat.dn(site, d)] ? 1 : -1;
                }

                if (prop_energy <= energy)
                    goto accept;
                else {
                    if (energy - prop_energy < std::log(1 - r_(rng_)) / beta_)
                        goto accept;
                }
                continue;

                accept:
                field[site] = prop;
                accepted += 1;
                energy = prop_energy;
            }

            return accepted;
        }

        void set_beta(double beta) { beta_ = beta; }

        double beta() const { return beta_; }

    private:
    private:
        double beta_;
        rng_t &rng_;
        std::uniform_int_distribution<field_t> u_;
        std::uniform_real_distribution<field_t> r_;
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
    int32_t M(const F &f) {
        return sum<int32_t>(f);
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
        return e / f.n_elements;
    }


}