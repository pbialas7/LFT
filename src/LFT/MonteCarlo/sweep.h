//
// Created by pbialas on 20.01.2022.
//
#pragma once

template<typename F, typename U>
typename F::lattice_t::size_t sweep(F &field, U update) {
    auto lat = field.lat;
    typename F::lattice_t::size_t acceptance = 0;
    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.even(i);
        acceptance += update(field, site);
    }

    for (auto i = 0; i < lat.n_elements / 2; ++i) {
        auto site = lat.odd(i);
        acceptance += update(field, site);
    }

    return acceptance;
}