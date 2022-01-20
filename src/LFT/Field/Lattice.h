//
// Created by pbialas on 20.01.2022.
//


#pragma once

#include <array>
#include <vector>
#include <cstdint>

template<typename I=uint32_t>

class Lattice {
public:
    static const int DIM = 2;
    using dim_t = std::array<int, DIM>;
    using size_t = I;
private:

    static size_t n_elem(dim_t dims) {
        return dims[0] * dims[1];
    }

public:
    Lattice(std::array<int, DIM> dims) : dims(dims), n_elements(n_elem(dims)), up_(DIM * n_elements),
                                         dn_(DIM * n_elements) {
        for (auto iy = 0; iy < dims[1]; ++iy) {
            for (auto ix = 0; ix < dims[0]; ++ix) {

                auto i = idx({ix, iy});

                if ((ix + iy) % 2 == 0)
                    even_.push_back(i);
                else
                    odd_.push_back(i);

                auto u_ = ix + 1;
                if (u_ >= dims[0])
                    u_ -= dims[0];
                up_[i * DIM + 0] = idx({u_, iy});

                auto d_ = ix - 1;
                if (d_ < 0)
                    d_ += dims[0];
                dn_[i * DIM + 0] = idx({d_, iy});

                u_ = iy + 1;
                if (u_ >= dims[1])
                    u_ -= dims[1];
                up_[i * DIM + 1] = idx({ix, u_});

                d_ = iy - 1;
                if (d_ < 0)
                    d_ += dims[1];
                dn_[i * DIM + 1] = idx({ix, d_});
            }
        }
    }

    size_t idx(dim_t coords) const { return dims[0] * coords[1] + coords[0]; }

    size_t up(size_t i, int dim) const { return up_[DIM * i + dim]; }

    size_t dn(size_t i, int dim) const { return dn_[DIM * i + dim]; }


    size_t even(size_t i) const { return even_[i]; }

    size_t odd(size_t i) const { return odd_[i]; }

    const dim_t dims;
    const size_t n_elements;
private:
    std::vector<size_t> up_;
    std::vector<size_t> dn_;

    std::vector<size_t> even_;
    std::vector<size_t> odd_;
};



