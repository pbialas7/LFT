// 
// Created by pbialas on 13.06.2026.
//

#pragma once

#include <random>
#include <vector>

#include "replicas.h"

namespace lft::ea {
    template <typename L, typename RNG>
    int choose_starting(const Replicas<L>& replicas, int r1, int r2, RNG& rng) {
        int chosen = -1;
        std::size_t count = 0;
        auto lattice = replicas.get_lattice();

        for (int i = 0; i < lattice.n_elements; ++i) {
            if ((*replicas[r1])[i] != (*replicas[r2])[i]) {
                count++;
                // Replace with prob 1/count
                std::uniform_int_distribution<std::size_t> dist(1, count);
                if (dist(rng) == 1) {
                    chosen = i;
                }
            }
        }
        return chosen; // -1 if none satisfy
    }

    template <typename L>
    int flip_cluster(Replicas<L>& replicas, int r1, int r2, int start) {
        assert((*replicas[r1])[start] != (*replicas[r2])[start]);

        auto lattice = replicas.get_lattice();
        std::vector<int> visited(lattice.n_elements, 0);
        std::vector<int> queue(std::max(lattice.dims[0], lattice.dims[1]), 0);

        queue.push_back(start);
        visited[start] = 1;
        auto count = 0;
        while (!queue.empty()) {
            auto element = queue.back();
            queue.pop_back();

            std::swap((*replicas[r1])[element], (*replicas[r2])[element]);
            count++;

            for (auto d = 0; d < L::DIM; ++d) {
                auto up = lattice.up(element, d);
                if (!visited[up] && ((*replicas[r1])[up] != (*replicas[r2])[up])) {
                    queue.push_back(up);
                    visited[up] = 1;
                }
                auto dn = lattice.dn(element, d);
                if (!visited[dn] && ((*replicas[r1])[dn] != (*replicas[r2])[dn])) {
                    queue.push_back(dn);
                    visited[dn] = 1;
                }
            }
        }
        return count;
    }
}
