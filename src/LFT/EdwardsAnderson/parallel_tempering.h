//
// Created by pbialas on 19.03.2026.
//

#pragma once

#include "ea.h"

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
                          const JField<float, L> &J_a) : replicas(n, Replicas<L>(q)), betas(betas), J(J),
                                                       heath_bath(n, J_a) {
            for (int i = 0; i < n; ++i) {
                heath_bath[i].set_beta(betas[i]);
            }
        }

        template<typename RNG>
        size_t sweep_mt(int n, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n; ++i) {
                replicas[i].sweep_mt(n, heath_bath[i], rng);
            }
            return acceptance;
        }

        template<typename RNG>
        size_t sweep_t1(int n, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n; ++i) {
                replicas[i].sweep_t1(n, heath_bath[i], rng);
            }
            return acceptance;
        }

        std::vector<Replicas<L> > replicas;
        std::vector<float> betas;
        JField<float, L> J;
        std::vector<HeathBath<float, L> > heath_bath;
    };
}
