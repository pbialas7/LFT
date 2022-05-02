// Created by pbialas on 29.04.2022.
//

#pragma once

#include <cmath>
#include <random>

#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace phi4 {
    template<typename L> using ScalarField = Field<float, L>;

    template<typename L, typename RNG=std::mt19937_64>
    struct Metropolis {
        using size_t = typename L::size_t;
        using Field = ScalarField<L>;
        using field_t = typename Field::field_t;
        using rng_t = RNG;
        using proposal_dist = std::normal_distribution<field_t>;
        using param_type = typename proposal_dist::param_type;

        Metropolis(field_t kappa, field_t m2, field_t lambda, rng_t &rng) : kappa_(kappa), m2_(m2), lambda_(lambda),
                                                                            rng_(rng), proposal_{},
                                                                            uni_{0.0, 1.0} {}

        void set_eps(field_t eps) {
            proposal_.param(param_type{0.0, eps});
        }

        field_t point_action(field_t phi, field_t corona) {
            return -kappa_ * phi * corona + phi * phi * (0.5 * (m2_ + 4.0 * kappa_) + 0.25 * lambda_ * phi * phi);
        }

        size_t operator()(Field &field, size_t site) {

            auto corona = field.corona(site);

            size_t accepted = 0;
            auto phi = field[site];
            auto action = point_action(phi, corona);
            for (size_t i = 0; i < 4; i++) {
                auto phi_proposed = phi + proposal_(rng_);
                auto action_proposed = point_action(phi_proposed, corona);


                auto action_diff = action_proposed - action;

                if (-action_diff >= 0) {
                    goto accept;
                } else {
                    auto r = uni_(rng_);
                    if (std::exp(-action_diff) > r)
                        goto accept;
                }
                continue;

                accept:
                action = action_proposed;
                field[site] = phi = phi_proposed;
                accepted += 1;

            }
            return accepted;
        }

    private:
        field_t kappa_;
        field_t m2_;
        field_t lambda_;
        rng_t &rng_;
        proposal_dist proposal_;
        std::uniform_real_distribution<field_t> uni_;
    };

    template<typename Float, typename F>
    Float magnetisation(const F &f) {
        return mean<Float>(f);
    }


    template<typename L, typename R>
    void hot_start(ScalarField<L> &field, R &rng) {
        using field_t = typename ScalarField<L>::field_t;
        std::normal_distribution<field_t> r(0.0, 1.0);
        for (auto site = 0; site < field.n_elements; ++site) {
            field[site] = r(rng);
        }
    }

    template<typename L, typename R>
    void cold_start(ScalarField<L> &field) {
        for (auto site = 0; site < field.n_elements; ++site) {
            field[site] = 1.0;
        }
    }
}