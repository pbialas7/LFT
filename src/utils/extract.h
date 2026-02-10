//
// Created by pbialas on 11.09.2025.
//

#pragma once

#include <vector>

#include "Field/Field.h"
#include "Field/Lattice.h"

template<typename F, typename L>
void extract_edges(const lft::Field<F, L> &field, std::vector<F> &edge) {
    using lattice_t = L;
    const lattice_t &lat = field.lat;
    static_assert(L::DIM == 2);
    const int d0 = field.lat.dims[0];
    const int d1 = field.lat.dims[1];
    const int len = d0 + d1 - 1;
    assert(edge.size()>=len);
    int i = 0;
    for (int j = 0; j < d0; j++) {
        edge[i] = field[lat.idx({i, 0})];
        ++i;
    }
    for (int j = 1; j < d1; ++j) {
        edge[i] = field[lat.idx({0, j})];
        ++i;
    }
}
