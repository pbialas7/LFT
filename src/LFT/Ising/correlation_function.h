//
// Created by pbialas on 05.04.23.
//


#pragma once

#include <iostream>

#include "ising.h"

template<typename Float, typename F, typename A>
void correlation(const F &f, A &out) {
    auto lat = f.lat;
    if (lat.dims[0] != lat.dims[1]) {
        std::cerr << "Lattice should be square for correlation function calculation.\n";
        exit(-1);
    }

    auto mag = ising::magnetisation<double>(f);

    auto n = lat.dims[0];
    auto L = lat.dims[0];
    for (int i = 0; i < lat.dims[0]; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; k++) {
                auto r = j + k;
                r = r < L ? r : r - L;
                out[k] += (f[{i, r}]) * (f[{i, j}]);
            }

    for (int i = 0; i < lat.dims[0]; ++i)
        for (int j = 0; j < L; ++j)
            for (int k = 0; k < n; k++) {
                auto r = j + k;
                r = r < L ? r : r - L;
                out[k] += (f[{r, i}]) * (f[{j, i}]);
            }

    for (int i = 0; i < n; ++i) {
        out[i] /= (2*lat.dims[0] * lat.dims[0]);
    }
}