//
// Tests for the EA building blocks in ea.h: local field, energy and the heat-bath sweep.
//

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <random>

#include "Field/Lattice.h"
#include "EdwardsAnderson/ea.h"
#include "EdwardsAnderson/parallel_tempering_mt.h"

namespace {
    template <typename L>
    struct EaSystem {
        L lat;
        lft::ea::JLattice<L> j_lat;
        lft::ea::SpinField<L> spins;
        lft::Field<float, lft::ea::JLattice<L>> J;

        static typename lft::ea::JLattice<L>::dim_t j_dims(const typename L::dim_t& dims) {
            typename lft::ea::JLattice<L>::dim_t jd;
            jd[0] = L::DIM;
            for (int d = 0; d < L::DIM; ++d) jd[d + 1] = dims[d];
            return jd;
        }

        EaSystem(const typename L::dim_t& dims, bool gaussian, std::uint64_t seed)
            : lat(dims, 'C'), j_lat(j_dims(dims), 'C'), spins(lat, 1), J(lft::make_field(j_lat, 1.0f)) {
            std::mt19937_64 rng(seed);
            if (gaussian)
                lft::ea::init_gaussian(J, rng);
            else
                lft::ea::init_bernoulli(J, rng);
            lft::ea::init_bernoulli(spins, rng);
        }
    };

    // The energy computed from the bonds to the upper neighbours, from the bonds to the lower neighbours and
    // from the full local field must all agree (each bond is counted once, once, and twice respectively).
    template <typename L>
    void check_energy_consistency(const typename L::dim_t& dims, bool gaussian) {
        using namespace lft::ea;
        EaSystem<L> sys(dims, gaussian, 123);
        double e_up = energy<double>(sys.spins, sys.J);
        double e_dn = 0.0, e_full = 0.0;
        for (std::size_t i = 0; i < sys.spins.n_elements; ++i) {
            e_dn -= sys.spins[i] * local_field_dn<double>(sys.spins, sys.J, i);
            e_full -= 0.5 * sys.spins[i] * local_field<double>(sys.spins, sys.J, i);
        }
        REQUIRE_THAT(e_dn, Catch::Matchers::WithinAbs(e_up, 1e-9 * (1 + std::abs(e_up))));
        REQUIRE_THAT(e_full, Catch::Matchers::WithinAbs(e_up, 1e-9 * (1 + std::abs(e_up))));
    }
}

TEST_CASE("EA energy equals the sum over upper, lower and all bonds", "[ea][energy]") {
    using L2 = lft::Lattice<uint32_t, 2>;
    using L3 = lft::Lattice<uint32_t, 3>;
    SECTION("2D +-J") { check_energy_consistency<L2>({6, 4}, false); }
    SECTION("2D gaussian") { check_energy_consistency<L2>({6, 4}, true); }
    SECTION("3D +-J") { check_energy_consistency<L3>({4, 2, 4}, false); }
    SECTION("3D gaussian") { check_energy_consistency<L3>({4, 2, 4}, true); }
}

TEST_CASE("EA local field of a ferromagnet", "[ea][energy]") {
    using L2 = lft::Lattice<uint32_t, 2>;
    lft::Lattice<uint32_t, 2> lat({4, 4}, 'C');
    L2::dim_t dims{4, 4};
    lft::ea::JLattice<L2> j_lat({2, 4, 4}, 'C');
    auto J = lft::make_field(j_lat, 1.0f);
    lft::ea::SpinField<L2> spins(lat, 1);
    for (std::size_t i = 0; i < spins.n_elements; ++i) {
        REQUIRE(lft::ea::local_field<float>(spins, J, i) == 4.0f);
    }
    REQUIRE(lft::ea::energy<double>(spins, J) == -2.0 * spins.n_elements);
}

TEST_CASE("HeathBath sweep updates every site once", "[ea][sweep]") {
    using L2 = lft::Lattice<uint32_t, 2>;
    EaSystem<L2> sys({8, 6}, true, 7);
    lft::ea::HeathBath<float, L2> hb(0.7f, sys.J);
    std::mt19937_64 rng(5);
    REQUIRE(hb.sweep(sys.spins, rng) == sys.spins.n_elements);

    EaSystem<L2> odd({3, 4}, true, 7);
    lft::ea::HeathBath<float, L2> hb_odd(0.7f, odd.J);
    REQUIRE_THROWS_AS(hb_odd.sweep(odd.spins, rng), std::invalid_argument);
}

TEST_CASE("ParallelTemperingMT Houdayer moves conserve E1 + E2", "[ea][houdayer]") {
    using L2 = lft::Lattice<uint32_t, 2>;
    EaSystem<L2> sys({8, 6}, true, 11);  // Gaussian J: the energies are not integers
    const std::vector<float> betas{0.4f, 0.9f};
    const int n_replicas = 4;
    lft::ea::ParallelTemperingMT<L2> pt(n_replicas, betas.size(), betas, sys.J);
    std::mt19937_64 init_rng(3);
    for (std::size_t b = 0; b < betas.size(); ++b)
        for (int r = 0; r < n_replicas; ++r)
            pt.chains[b][r] = new lft::ea::SpinField<L2>(sys.lat, 1);

    std::vector<std::mt19937_64> rngs{std::mt19937_64(5), std::mt19937_64(6)};  // one stream per chain
    std::size_t total_cluster_size = 0;
    for (int move = 0; move < 200; ++move) {
        double before = 0.0;
        for (std::size_t b = 0; b < betas.size(); ++b)
            for (int r = 0; r < n_replicas; ++r) {
                lft::ea::init_bernoulli(*pt.chains[b][r], init_rng);
                before += lft::ea::energy<double>(*pt.chains[b][r], sys.J);
            }
        total_cluster_size += pt.houdayer(rngs);  // in Debug builds this also runs the internal assert
        double after = 0.0;
        for (std::size_t b = 0; b < betas.size(); ++b)
            for (int r = 0; r < n_replicas; ++r)
                after += lft::ea::energy<double>(*pt.chains[b][r], sys.J);
        REQUIRE_THAT(after, Catch::Matchers::WithinAbs(before, 1e-9 * (std::abs(before) + 1)));
    }
    REQUIRE(total_cluster_size > 0);
}
