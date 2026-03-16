//
// Created by pbialas on 20.01.2022.
//
#pragma once
#include <algorithm>
#include <vector>
#include <iostream>

#include <omp.h>

template<typename F, typename U, typename RNG>
F::lattice_t::size_t sweep(F &field, U update) {
    auto lat = field.lat;
    typename F::lattice_t::size_t accepted = 0;
    for (size_t i = 0; i < lat.n_elements / 2; i++) {
        auto site = lat.even(i);
        accepted += update(field, site);
    }

    for (auto i = 0; i < lat.n_elements / 2; i++) {
        auto site = lat.odd(i);
        accepted += update(field, site);
    }

    return accepted;
}

template<typename F, typename U, typename RNG>
F::lattice_t::size_t sweep_mt(F &field, U update, RNG &rng) {
    auto lat = field.lat;
    typename F::lattice_t::size_t accepted = 0;
#pragma omp parallel for reduction(+:accepted)
    for (size_t i = 0; i < lat.n_elements / 2; i++) {
        auto t = omp_get_thread_num();
        auto site = lat.even(i);
        accepted += update(field, site, rng[t]);
    }

#pragma omp parallel for reduction(+:accepted)
    for (auto i = 0; i < lat.n_elements / 2; i++) {
        auto t = omp_get_thread_num();
        auto site = lat.odd(i);
        accepted += update(field, site, rng[t]);
    }
    return accepted;
}


template<typename F, typename U>
typename F::lattice_t::size_t sweep(F &field, U update, std::vector<typename F::lattice_t::size_t> &fixed_sites) {
    auto lat = field.lat;
    typename F::lattice_t::size_t accepted = 0;
    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.even(i);
        if (!std::binary_search(fixed_sites.begin(), fixed_sites.end(), site))
            accepted += update(field, site);
    }

    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.odd(i);
        if (!std::binary_search(fixed_sites.begin(), fixed_sites.end(), site))
            accepted += update(field, site);
    }

    return accepted;
}
