//
// Created by pbialas on 19.03.2026.
//

#pragma once

#include <iostream>
#include <spdlog/spdlog.h>

#include "ea.h"
#include "replicas.h"
#include "houdayer.h"


namespace lft::ea {
    template <typename L>
    struct ParallelTemperingMT {
        ParallelTemperingMT(int q, int n, const std::vector<float>& betas,
                            const JField<float, L>& J_a) : q(q), chains(n, Replicas<L>(q)), betas(betas), J(J_a),
                                                           heath_bath(n, J_a), accepted_v(n, 0), exchange_v(n, 0),
                                                           u(0.0, 1.0) {
            assert(betas.size() == n);
            for (int i = 0; i < n; ++i) {
                heath_bath[i].set_beta(betas[i]);
            }
        }


        template <typename RNG>
        size_t sweep_t1(int n_sweeps, RNG& rng) {
            size_t acceptance = 0;
            for (int j = 0; j < betas.size(); j++)
                chains[j].sweep_t1(n_sweeps, heath_bath[j], rng);
            return acceptance;
        }

        template <typename RNG>
        size_t sweep_mt(int n_sweeps, RNG& rng) {
            size_t acceptance = 0;
#pragma omp parallel for shared(rng) reduction(+: acceptance) schedule(static)
            for (int j = 0; j < betas.size(); j++) {
                auto t = omp_get_thread_num();
                chains[j].sweep_t1(n_sweeps, heath_bath[j], rng[t]);
            }
            return acceptance;
        }


        template <typename RNG>
        size_t exchange(int i_b1, int i_b2, int i_r, RNG& rng) {
            spdlog::trace("exchange {} {} {}", i_b1, i_b2, i_r);
            float b1 = betas[i_b1];
            float b2 = betas[i_b2];
            spdlog::trace("Exchange {} {}", b1, b2);
            spdlog::trace("Exchange {} {}", (const void*)chains[i_b1][i_r], (const void*)chains[i_b2][i_r]);
            float E1 = energy<float>(*chains[i_b1][i_r], J);
            float E2 = energy<float>(*chains[i_b2][i_r], J);
            float delta = (E1 - E2) * (b1 - b2);
            auto r = u(rng);

            if (r < std::exp(delta)) {
                std::swap(chains[i_b1][i_r], chains[i_b2][i_r]);
                return 1;
            }

            return 0;
        }

        template <typename RNG>
        size_t exchange(int i_r, RNG& rng) {
            size_t accepted = 0;
            for (int j = 0; j < betas.size() - 1; j += 2) {
                auto a = exchange(j, j + 1, i_r, rng);
                accepted += a;
                accepted_v[j] += a;
                exchange_v[j]++;
            }
            for (int j = 1; j < betas.size() - 1; j += 2) {
                auto a = exchange(j, j + 1, i_r, rng);
                accepted += a;
                accepted_v[j] += a;
                exchange_v[j]++;
            }
            return accepted;
        }


        template <typename RNG>
        size_t exchange(RNG& rng) {
            size_t accepted = 0;
            for (int j = 0; j < q; j++)
                accepted += exchange(j, rng);
            return accepted;
        }

        template <typename RNG>
        size_t exchange_mt(RNG& rng) {
            size_t accepted = 0;
#pragma omp parallel for shared(rng,accepted_v, exchange_v)  schedule(static)
            for (int j = 0; j < betas.size() - 1; j += 2) {
                auto t = omp_get_thread_num();
                for (int i_r = 0; i_r < q; i_r++) {
                    auto a = exchange(j, j + 1, i_r, rng[t]);
                    accepted += a;
                    accepted_v[j] += a;
                    exchange_v[j]++;
                }
            }
#pragma omp parallel for shared(rng,accepted_v, exchange_v)  schedule(static)
            for (int j = 1; j < betas.size() - 1; j += 2) {
                auto t = omp_get_thread_num();
                for (int i_r = 0; i_r < q; i_r++) {
                    auto a = exchange(j, j + 1, i_r, rng[t]);
                    accepted += a;
                    accepted_v[j] += a;
                    exchange_v[j]++;
                }
            }
            return accepted;
        }


        template <typename RNG>
        size_t houdayer(Replicas<L>& chain, int r1, int r2, RNG& rng) {
            int i_start = choose_starting(chain, r1, r2, rng);

            auto E1 = energy<float>(*chain[r1], J);
            auto E2 = energy<float>(*chain[r2], J);
            size_t cluster_size = 0;
            if (i_start >= 0) {
                cluster_size = flip_cluster(chain, r1, r2, i_start);
            }
            assert(cluster_size > 0);
            auto E1_after = energy<float>(*chain[r1], J);
            auto E2_after = energy<float>(*chain[r2], J);

            spdlog::trace("houdayer: E1 {} E2 {} E1_after {} E2_after {}", E1, E2, E1_after, E2_after);

            assert(E1_after + E2_after == E1+E2);

            return cluster_size;
        }

        template <typename RNG>
        size_t houdayer(RNG& rng) {
            size_t cluster_size = 0;
            for (int j = 0; j < betas.size(); j++) {
                cluster_size += houdayer(chains[j], 0, 1, rng);
            }
            return cluster_size;
        }

        void print_stats(std::ostream& os) {
            for (auto i = 0; i < betas.size() - 1; i++) {
                os << std::format("{:.3f} {:.3f} {:.2f}", betas[i], betas[i + 1],
                                  (double)accepted_v[i] / exchange_v[i]) << std::endl;
            }
        }

        void reset_stats() {
            std::fill(accepted_v.begin(), accepted_v.end(), 0);
            std::fill(exchange_v.begin(), exchange_v.end(), 0);
        }

        int q;
        std::uniform_real_distribution<float> u;
        std::vector<Replicas<L>> chains;
        std::vector<float> betas;
        JField<float, L> J;
        std::vector<HeathBath<float, L>> heath_bath;
        std::vector<size_t> accepted_v;
        std::vector<size_t> exchange_v;
    };
}
