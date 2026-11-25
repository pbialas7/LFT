//
// Created by pbialas on 25.11.2025.
//

#pragma once


#include <cmath>
#include <random>

#include <omp.h>
#include <spdlog/spdlog.h>

#include "Field/Field.h"
#include "MonteCarlo/sweep.h"

namespace lft::ea {
    template <typename L>
    using SpinField = lft::Field<int8_t, L>;
    template <typename L>
    using JLattice = lft::Lattice<typename L::index_t, L::DIM + 1>;
    template <typename F, typename L>
    using JField = lft::Field<F, JLattice<L>>;


    template <typename F, typename RNG>
    void init_bernoulli(F& j_field, double p, RNG& rng) {
        std::bernoulli_distribution bern(p);
        for (int i = 0; i < j_field.n_elements; ++i) {
            if (bern(rng)) {
                j_field[i] = 1;
            }
            else {
                j_field[i] = -1;
            }
        }
    }

    template <typename F, typename RNG>
    void init_bernoulli(F& j_field, RNG& rng) {
        init_bernoulli(j_field, 0.5, rng);
    }


    template <typename F, typename RNG>
    void init_gaussian(F& j_field,
                       typename F::field_t mu,
                       typename F::field_t sigma,
                       RNG& rng) {
        using f_type = F::field_t;
        std::normal_distribution<f_type> normal(mu, sigma);
        for (int i = 0; i < j_field.n_elements; ++i) {
            j_field[i] = normal(rng);
        }
    }

    template <typename F, typename RNG>
    void init_gaussian(F& j_field, RNG& rng) {
        using f_type = F::field_t;
        init_gaussian(j_field, f_type(0.0), f_type(1.0), rng);
    }


    /*
     * Sums of s_j * J_ij over the bonds attached to site i. The coupling of the bond from site i in direction d is
     * stored at j[i + n_elements*d] (see JField), so the bond to the lower neighbour dn is stored at dn's index.
     * Terms are added to `corona` in order, so the up and down sums can be chained into one accumulator.
     */
    template <typename Acc, typename SF, typename JF>
    Acc local_field_up(const SF& field, const JF& j, std::size_t i, Acc corona = Acc(0)) {
        for (auto d = 0; d < SF::DIM; ++d) {
            auto idx = i + field.lat.n_elements * d;
            corona += field[field.lat.up(i, d)] * j[idx];
        }
        return corona;
    }

    template <typename Acc, typename SF, typename JF>
    Acc local_field_dn(const SF& field, const JF& j, std::size_t i, Acc corona = Acc(0)) {
        for (auto d = 0; d < SF::DIM; ++d) {
            auto dn_site = field.lat.dn(i, d);
            auto idx = dn_site + field.lat.n_elements * d;
            corona += field[dn_site] * j[idx];
        }
        return corona;
    }

    // Local field on site i: the sum over all 2*DIM bonds attached to it.
    template <typename Acc, typename SF, typename JF>
    Acc local_field(const SF& field, const JF& j, std::size_t i) {
        return local_field_dn(field, j, i, local_field_up<Acc>(field, j, i));
    }


    template <typename F, typename L>
    struct HeathBath {
        using size_t = L::size_t;
        using index_t = L::index_t;

        using field_class = SpinField<L>;
        using field_t = field_class::field_t;


        HeathBath(const JField<F, L>& j) : j_(j), j_lat_(j.lat), beta_() {}

        HeathBath(F beta, const JField<F, L>& j) : beta_(beta), j_(j), j_lat_(j.lat) {}

        void set_beta(F beta) { beta_ = beta; }


        template <typename R>
        size_t update(field_class& field, size_t i, R& rng) {
            return update_with_dist(field, i, rng, u_);
        }

        template <typename R, typename UDist>
        size_t update_with_dist(field_class& field, size_t i, R& rng, UDist& u) {
            auto r = std::log(F(1.0) / u(rng) - F(1.0));
            F corona = local_field<F>(field, j_, i);

            F r_p_up = -F(2.0) * beta_ * corona;

            if (r > r_p_up)
                field[i] = 1;
            else
                field[i] = -1;
            return 1;
        }

        template <typename R>
        size_t operator()(field_class& field, size_t i, R& rng) {
            return update(field, i, rng);
        }


        // One checkerboard sweep (even sites, then odd sites) using the generic lft::sweep.
        template <typename RNG>
        size_t sweep(SpinField<L>& field, RNG& rng) {
            auto site_update = [&](field_class& f, size_t site) { return update(f, site, rng); };
            return ::sweep(field, site_update);
        }

    private:
        F beta_;
        std::uniform_real_distribution<F> u_;
        JField<F, L> j_;
        JLattice<L> j_lat_;
    };

    template <typename Float, typename F>
    Float magnetisation(const F& f) {
        return sum<Float>(f);
    }

    template <typename Float, typename F, typename JF>
    Float energy(const F& field, const JF& j_) {
        Float e = 0.0;
        for (int i = 0; i < field.n_elements; ++i) {
            double corona = local_field_up<double>(field, j_, i);
            e += field[i] * corona;
        }

        return -e;
    }

    template <typename Float, typename F>
    Float overlap(const F& f1, const F& f2) {
        Float overlap_ = 0.0;
        for (int i = 0; i < f1.n_elements; ++i) {
            overlap_ += f1[i] * f2[i];
        }
        return overlap_ / f1.n_elements;
    }

    template <typename Float, typename F>
    Float link_overlap(const F& field1, const F& field2) {
        Float link_overlap_ = 0.0;
        for (int i = 0; i < field1.n_elements; ++i) {
            for (auto d = 0; d < F::DIM; ++d) {
                auto i_up = field1.lat.up(i, d);
                auto s_up1 = field1[i_up];
                auto s_up2 = field2[i_up];
                link_overlap_ += field1[i] * s_up1 * field2[i] * s_up2;
            }
        }

        return -link_overlap_ / (field1.n_elements * F::DIM);
    }
}
