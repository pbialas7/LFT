//
// Created by pbialas on 19.03.2026.
//

#pragma once

#include "ea.h"
#include <spdlog/spdlog.h>


namespace lft::ea {
    template<typename L>
    struct Replicas {
        Replicas(int q) : q(q), replica(q, nullptr) {
        }

        ~Replicas() {
            for (int i = 0; i < q; ++i) {
                if (replica[i])
                    delete replica[i];
            }
        }

        const int q;
        std::vector<SpinField<L> *> replica;

        SpinField<L> *&operator[](int i) {
            return replica[i];
        }


        template<typename RNG>
        size_t sweep_t1(int n, HeathBath<float, L> &heath_bath, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < q; ++j) {
                    acceptance += heath_bath.sweep(*replica[j], rng);
                }
            }
            return acceptance;
        }

        template<typename RNG>
        size_t sweep_mt(int n, HeathBath<float, L> &heath_bath, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < q; ++j) {
                    acceptance += heath_bath.sweep_mt(*replica[j], rng);
                }
            }
            return acceptance;
        }
    };

    template<typename L>
    struct ParallelTempering {
        ParallelTempering(int q, int n, const std::vector<float> &betas,
                          const JField<float, L> &J_a) : q(q), replicas(n, Replicas<L>(q)), betas(betas), J(J_a),
                                                         heath_bath(n, J_a) {
            assert(betas.size() == n);
            for (int i = 0; i < n; ++i) {
                heath_bath[i].set_beta(betas[i]);
            }
        }

        template<typename RNG>
        size_t sweep_mt(int n_sweeps, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n_sweeps; ++i) {
                for (int j = 0; j < betas.size(); j++)
                    replicas[j].sweep_mt(n_sweeps, heath_bath[j], rng);
            }
            return acceptance;
        }

        template<typename RNG>
        size_t sweep_t1(int n_sweeps, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n_sweeps; ++i) {
                for (int j = 0; j < betas.size(); j++)
                    replicas[j].sweep_t1(n_sweeps, heath_bath[j], rng);
            }
            return acceptance;
        }


        template<typename RNG>
        size_t exchange(int i_b1, int i_b2, int i_r, RNG &rng) {
            spdlog::trace("exchange {} {} {}", i_b1, i_b2, i_r);
            float b1 = betas[i_b1];
            float b2 = betas[i_b2];
            spdlog::trace("Exchange {} {}", b1, b2);
            spdlog::trace("Exchange {} {}", (const void *) replicas[i_b1][i_r], (const void *) replicas[i_b2][i_r]);
            float E1 = energy<float>(*replicas[i_b1][i_r], J);
            float E2 = energy<float>(*replicas[i_b2][i_r], J);
            float delta = (E1 - E2) * (b1 - b2);
            auto r = u(rng);

            if (r < std::exp(delta)) {
                std::swap(replicas[i_b1][i_r], replicas[i_b2][i_r]);
                return 1;
            }

            return 0;
        }

        template<typename RNG>
        size_t exchange(int i_r, RNG &rng) {
            size_t acceptance = 0;
            for (int j = 0; j < betas.size() - 1; j += 2)
                acceptance += exchange(j, j + 1, i_r, rng);
            for (int j = 1; j < betas.size() - 1; j += 2)
                acceptance += exchange(j, j + 1, i_r, rng);
            return acceptance;
        }

        template<typename RNG>
        size_t exchange(RNG &rng) {
            size_t acceptance = 0;
            for (int j = 0; j < q; j++)
                acceptance += exchange(j, rng);
        }

        int q;
        std::uniform_real_distribution<float> u;
        std::vector<Replicas<L> > replicas;
        std::vector<float> betas;
        JField<float, L> J;
        std::vector<HeathBath<float, L> > heath_bath;
    };
}
