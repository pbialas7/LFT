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

    auto n = lat.dims[0] / 2;
    for (int i = 0; i < lat.dims[0]; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; k++) {
                out[k] += (f[{i, j + k}] - mag) * (f[{i, j}] - mag);
            }

    for (int i = 0; i < lat.dims[0]; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; k++) {
                out[k] += (f[{j + k, i}] - mag) * (f[{j, i}] - mag);
            }

    for (int i = 0; i < n; ++i) {
        out[i] /= (0.5 * lat.dims[0] * lat.dims[0]);
        out[i] /= out[0];
    }
}