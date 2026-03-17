//
// Created by pbialas on 20.01.2022.
//
#pragma once


template<typename F, typename U>
F::lattice_t::size_t sweep(F &field, U &update) {
    typename F::lattice_t::size_t accepted = 0;

    for (size_t i = 0; i < field.lat.n_elements / 2; i++) {
        auto site = field.lat.even(i);
        accepted += update(field, site);
    }


    for (auto i = 0; i < field.lat.n_elements / 2; i++) {
        auto site = field.lat.odd(i);
        accepted += update(field, site);
    }

    return accepted;
}
