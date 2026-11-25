//
// Regression tests for the Houdayer cluster move (lft-90f).
//

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <random>
#include <vector>

#include "Field/Lattice.h"
#include "EdwardsAnderson/ea.h"
#include "EdwardsAnderson/replicas.h"
#include "EdwardsAnderson/houdayer.h"

namespace {
    // Independent BFS: the connected component of sites with s1 != s2 that contains `start`.
    template <typename L>
    std::vector<char> q_cluster(const lft::ea::Replicas<L>& r, int start) {
        const auto& f1 = *r[0];
        const auto& f2 = *r[1];
        const auto& lat = f1.lat;
        std::vector<char> in(lat.n_elements, 0);
        std::vector<std::size_t> stack{static_cast<std::size_t>(start)};
        in[start] = 1;
        while (!stack.empty()) {
            auto s = stack.back();
            stack.pop_back();
            for (int d = 0; d < L::DIM; ++d) {
                for (auto n : {lat.up(s, d), lat.dn(s, d)}) {
                    if (!in[n] && f1[n] != f2[n]) {
                        in[n] = 1;
                        stack.push_back(n);
                    }
                }
            }
        }
        return in;
    }

    template <typename L>
    void check_houdayer(const typename L::dim_t& dims, int n_moves, std::uint64_t seed) {
        using namespace lft::ea;
        L lat(dims, 'C');
        typename JLattice<L>::dim_t j_dims;
        j_dims[0] = L::DIM;
        for (int d = 0; d < L::DIM; ++d) j_dims[d + 1] = dims[d];
        JLattice<L> j_lat(j_dims, 'C');

        std::mt19937_64 rng(seed);
        auto J = lft::make_field(j_lat, 1.0f);
        init_gaussian(J, rng);

        Replicas<L> replicas(2);
        replicas[0] = new SpinField<L>(lat, 1);
        replicas[1] = new SpinField<L>(lat, 1);

        for (int m = 0; m < n_moves; ++m) {
            // Fresh random configurations each time so many cluster shapes are exercised.
            init_bernoulli(*replicas[0], rng);
            init_bernoulli(*replicas[1], rng);

            int start = choose_starting(replicas, 0, 1, rng);
            if (start < 0) continue;

            auto before1 = *replicas[0];
            auto before2 = *replicas[1];
            double E_before = energy<double>(*replicas[0], J) + energy<double>(*replicas[1], J);
            auto cluster = q_cluster(replicas, start);
            std::size_t cluster_size = 0;
            for (auto c : cluster) cluster_size += c;

            int count = flip_cluster(replicas, 0, 1, start);

            REQUIRE(count == static_cast<int>(cluster_size));

            // Exactly the cluster sites are swapped between the replicas; nothing else changes.
            for (std::size_t i = 0; i < lat.n_elements; ++i) {
                if (cluster[i]) {
                    REQUIRE((*replicas[0])[i] == before2[i]);
                    REQUIRE((*replicas[1])[i] == before1[i]);
                } else {
                    REQUIRE((*replicas[0])[i] == before1[i]);
                    REQUIRE((*replicas[1])[i] == before2[i]);
                }
            }

            double E_after = energy<double>(*replicas[0], J) + energy<double>(*replicas[1], J);
            REQUIRE_THAT(E_after, Catch::Matchers::WithinAbs(E_before, 1e-9));
        }
    }
}

TEST_CASE("Houdayer flip_cluster flips exactly the q=-1 cluster", "[ea][houdayer]") {
    using L2 = lft::Lattice<uint32_t, 2>;
    using L3 = lft::Lattice<uint32_t, 3>;

    SECTION("2D odd 3x3") { check_houdayer<L2>({3, 3}, 200, 11); }
    SECTION("2D 5x4") { check_houdayer<L2>({5, 4}, 200, 12); }
    SECTION("2D 4x5") { check_houdayer<L2>({4, 5}, 200, 13); }
    SECTION("2D even 6x6") { check_houdayer<L2>({6, 6}, 200, 14); }
    SECTION("3D odd 3x3x3") { check_houdayer<L3>({3, 3, 3}, 100, 15); }
    SECTION("3D even 4x4x4") { check_houdayer<L3>({4, 4, 4}, 100, 16); }
}
