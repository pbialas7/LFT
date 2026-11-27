//
// Created by pbialas on 19.03.2026.
//

#pragma once

#include <algorithm>
#include <cassert>
#include <format>
#include <iostream>
#include <vector>

#include <omp.h>
#include <spdlog/spdlog.h>

#include "ea.h"
#include "replicas.h"
#include "houdayer.h"


namespace lft::ea {
    template <typename L>
    struct ParallelTemperingMT {
        ParallelTemperingMT(int n_replicas, int n, const std::vector<float>& betas,
                            const JField<float, L>& J_a) : n_replicas(n_replicas), chains(n, Replicas<L>(n_replicas)), betas(betas), J(J_a),
                                                           heath_bath(n, J_a),
                                                           accepted_v(n, 0), exchange_v(n, 0), cluster_size_v(n, 0) {
            assert(betas.size() == n);
            for (int i = 0; i < n; ++i) {
                heath_bath[i].set_beta(betas[i]);
            }
        }


        /*
         * Random streams: every chain j has its own stream rngs[j] (rngs must provide at least betas.size()
         * streams via rngs[j]). Chain j uses it for its sweeps and Houdayer moves, and the exchange of the pair
         * (j, j+1) uses rngs[j]. The streams therefore do not depend on which OpenMP thread runs a chain, and the
         * results are the same for any number of threads. Chains, and pairs of the same parity, are processed in
         * parallel with `omp parallel for`; with one thread the loops run in order.
         */

        // Sweeps all replicas of every chain n_sweeps times. Returns the number of accepted updates.
        template <typename RNGS>
        size_t sweep(int n_sweeps, RNGS& rngs) {
            size_t acceptance = 0;
            const int n_chains = betas.size();
#pragma omp parallel for shared(rngs) reduction(+: acceptance) schedule(static)
            for (int j = 0; j < n_chains; j++) {
                acceptance += chains[j].sweep_t1(n_sweeps, heath_bath[j], rngs[j]);
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
            std::uniform_real_distribution<float> u(0.0f, 1.0f);
            auto r = u(rng);

            if (r < std::exp(delta)) {
                std::swap(chains[i_b1][i_r], chains[i_b2][i_r]);
                return 1;
            }

            return 0;
        }

        // Attempts exchanges between neighbouring betas for all replicas: first the pairs (0,1), (2,3), ..., then
        // (1,2), (3,4), .... Returns the number of accepted exchanges.
        template <typename RNGS>
        size_t exchange(RNGS& rngs) {
            size_t accepted = 0;
            const int n_chains = betas.size();
            // Pairs of the same parity do not overlap, so they can be processed in parallel.
            for (int parity = 0; parity < 2; ++parity) {
#pragma omp parallel for shared(rngs, accepted_v, exchange_v) reduction(+: accepted) schedule(static)
                for (int j = parity; j < n_chains - 1; j += 2) {
                    for (int i_r = 0; i_r < n_replicas; i_r++) {
                        auto a = exchange(j, j + 1, i_r, rngs[j]);
                        accepted += a;
                        accepted_v[j] += a;
                        exchange_v[j]++;
                    }
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

            auto E1_after = energy<float>(*chain[r1], J);
            auto E2_after = energy<float>(*chain[r2], J);

            spdlog::trace("houdayer: E1 {} E2 {} E1_after {} E2_after {}", E1, E2, E1_after, E2_after);

            assert(E1_after + E2_after == E1+E2);

            return cluster_size;
        }

        // Houdayer cluster moves between replicas i and i + n_replicas/2 of every chain. Returns the total size of
        // the flipped clusters; the per-chain totals are accumulated in cluster_size_v.
        template <typename RNGS>
        size_t houdayer(RNGS& rngs) {
            size_t cluster_size = 0;
            const int n_chains = betas.size();
#pragma omp parallel for shared(rngs, cluster_size_v) reduction(+: cluster_size) schedule(static)
            for (int j = 0; j < n_chains; j++) {
                size_t chain_size = 0;
                for (int i_r = 0; i_r < n_replicas / 2; i_r++)
                    chain_size += houdayer(chains[j], i_r, i_r + n_replicas / 2, rngs[j]);
                cluster_size_v[j] += chain_size;
                cluster_size += chain_size;
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

        int n_replicas;
        std::vector<Replicas<L>> chains;
        std::vector<float> betas;
        JField<float, L> J;
        std::vector<HeathBath<float, L>> heath_bath;
        std::vector<size_t> accepted_v;
        std::vector<size_t> exchange_v;
        std::vector<size_t> cluster_size_v;
    };
}
