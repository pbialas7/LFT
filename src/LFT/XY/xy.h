// Created by pbialas on 29.04.2022.
//

#pragma once

#include <cmath>
#include <random>

#include "spdlog/spdlog.h"

#include "Field/Field.h"

namespace xy {

    constexpr double PI = 3.1415926535897932384626433832795;
    constexpr double M_2PI = 2*PI;


    template<typename L> using ScalarField = Field<float, L>;

    template<typename L, typename RNG=std::mt19937_64>
    struct Metropolis {
        using size_t = typename L::size_t;
        using Field = ScalarField<L>;
        using field_t = typename Field::field_t;
        using rng_t = RNG;
        using proposal_dist = std::normal_distribution<field_t>;
        using param_type = typename proposal_dist::param_type;
        static const int DIM = L::DIM;


        Metropolis(field_t beta, int n_hits, rng_t &rng) : beta_(beta), n_hits_(n_hits), rng_(rng), proposal_{},
                                                           uni_{0.0, 1.0} {}

        void set_eps(field_t eps) {
            proposal_.param(param_type{0.0, eps});
        }

        field_t point_action(const Field &field, size_t site, field_t phi) {
            field_t action = 0.0;
            for (auto d = 0; d < DIM; ++d) {
                action += std::cos(field[field.lat.up(site, d)] - phi);
                action += std::cos(field[field.lat.dn(site, d)] - phi);
            }
            return action;
        }

        size_t operator()(Field &field, size_t site) {


            size_t accepted = 0;
            auto phi = field[site];
            auto action = point_action(field, site, phi);
            for (size_t i = 0; i < n_hits_; i++) {
                auto phi_proposed = phi + proposal_(rng_);
                if (phi_proposed > M_2PI)
                    phi_proposed -= M_2PI;
                else if (phi_proposed < 0.0)
                    phi_proposed += M_2PI;
                auto action_proposed = point_action(field, site, phi_proposed);


                auto action_diff = action_proposed - action;

                if (beta_ * action_diff >= 0) {
                    goto accept;
                } else {
                    auto r = uni_(rng_);
                    if (std::exp(beta_ * action_diff) > r)
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
        field_t beta_;
        rng_t &rng_;
        int n_hits_;
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
        std::uniform_real_distribution<field_t> r(0.0, 2 * PI);
        for (auto site = 0; site < field.n_elements; ++site) {
            field[site] = r(rng);
        }
    }

    template<typename L>
    void cold_start(ScalarField<L> &field) {
        for (auto site = 0; site < field.n_elements; ++site) {
            field[site] = 0.0;
        }
    }

    template<typename F, typename R = typename F::field_t>
    auto action(const F &field) {
        using L = typename F::lattice_t;
        R action = 0.0;
        R specific_heat = 0.0;
        for (typename L::size_t site = 0; site < field.lat.n_elements; ++site) {
            auto phi = field[site];
            R c = 0.0;
            for (auto d = 0; d < L::DIM; ++d) {
                c -= std::cos(field[field.lat.up(site, d)] - phi);
            }
            action += c;
            specific_heat += c * c;
        }

        return std::make_pair(action, specific_heat);
    }

    template<typename F, typename R = typename F::field_t>
    auto magnetisation(const F &field) {
        R mag_x = 0.0;
        R mag_y = 0.0;
        R mag2_x = 0.0;
        R mag2_y = 0.0;
        for (typename F::lattice_t::size_t site = 0; site < field.lat.n_elements; ++site) {
            auto c = std::cos(field[site]);
            auto s = std::sin(field[site]);
            mag_x += c;
            mag_y += s;
            mag2_x += c * c;
            mag2_y += s * s;
        }
        return std::make_tuple(mag_x, mag_y, mag2_x + mag2_y);
    }

}