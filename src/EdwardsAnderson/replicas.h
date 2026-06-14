//
// Created by pbialas on 12.06.2026.
//

#pragma once

#include "ea.h"

namespace lft::ea {
    template<typename L>
    struct Replicas {
        Replicas(int n_replicas) : n_replicas(n_replicas), replica(n_replicas, nullptr) {
        }

        ~Replicas() {
            for (int i = 0; i < n_replicas; ++i) {
                if (replica[i])
                    delete replica[i];
            }
        }

        const int n_replicas;
        std::vector<SpinField<L> *> replica;

        SpinField<L> *&operator[](int i) {
            return replica[i];
        }

        SpinField<L> *operator[](int i) const {
            return replica[i];
        }

        auto get_lattice()const  {
            return replica[0]->lat;
        }

        template<typename RNG>
        size_t sweep_t1(int n_sweeps, HeathBath<float, L> &heath_bath, RNG &rng) {
            size_t acceptance = 0;
            for (int i = 0; i < n_sweeps; ++i) {
                for (int j = 0; j < n_replicas; ++j) {
                    acceptance += heath_bath.sweep(*replica[j], rng);
                }
            }
            return acceptance;
        }
    };
}
