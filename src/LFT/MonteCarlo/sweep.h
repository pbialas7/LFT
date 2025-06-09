//
// Created by pbialas on 20.01.2022.
//
#pragma once

template<typename F, typename U>
typename F::lattice_t::size_t sweep(F &field, U update) {
    auto lat = field.lat;
    typename F::lattice_t::size_t accepted = 0;
    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.even(i);
        accepted += update(field, site);
    }

    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.odd(i);
        accepted += update(field, site);
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