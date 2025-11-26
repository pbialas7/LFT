//
// Created by pbialas on 20.01.2022.
//


#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include <spdlog/spdlog.h>


template<typename I=uint32_t, int D = 2>
class MultiIndex {
public:
    static const int DIM = D;
    using index_t = I;
    using dims_t = std::array<index_t, DIM>;

    MultiIndex(dims_t dims) : dims(dims) {
        std::fill_n(coords_.begin(), DIM, 0);
        set_volume();
    }

    MultiIndex(dims_t dims, dims_t idx) : dims(dims), coords_(idx) {
        set_volume();
    }

    MultiIndex &operator++() {
        coords_[0]++;
        for (auto i = 0; i < DIM - 1; i++) {
            if (coords_[i] >= dims[i]) {
                coords_[i + 1]++;
                coords_[i] = 0;
            }
        }
        return *this;
    }

    dims_t coords() const { return coords_; }

    size_t idx() const {
        index_t idx = 0;
        for (auto i = 0; i < DIM; i++) {
            idx += volumes_[i] * coords_[i];
        }
        return idx;
    }

    index_t volume() const { return volumes_[DIM]; }

    const std::array<index_t, DIM> dims;

private:
    void set_volume() {
        volumes_[0] = 1;
        for (auto i = 1; i < DIM + 1; ++i) {
            volumes_[i] = volumes_[i - 1] * dims[i - 1];
        }
    }

    std::array<index_t, DIM> coords_;

public:
    std::array<index_t, DIM + 1> volumes_;
};

template<typename I=uint32_t, int D = 2>
class Lattice {
public:
    static const int DIM = D;
    using dim_t = std::array<int, DIM>;
    using size_t = I;

private:
    static std::size_t n_elem(dim_t dims) {
        std::size_t n = dims[0];
        for (int i = 1; i < DIM; i++)
            n *= dims[i];
        return n;
    }

public:
    Lattice(std::array<int, DIM> dims) : dims(dims), n_elements(n_elem(dims)),
                                         nn_(2 * DIM * n_elements) {
        MultiIndex<int, DIM> m_index(dims);
        for (auto i = 0; i < m_index.volume(); i++, ++m_index) {
            auto idx_ = i;
            auto coords = m_index.coords();
            size_t i_sum = 0;
            for (int j = 0; j < DIM; j++)
                i_sum += coords[j];

            if (i_sum % 2 == 0) {
                even_.push_back(idx_);
            } else {
                odd_.push_back(idx_);
            }

            for (auto j = 0; j < DIM; j++) {
                auto u_ = coords[j] + 1;
                if (u_ >= dims[j])
                    u_ -= dims[j];
                auto u_coords_ = coords;
                u_coords_[j] = u_;


                nn_[2 * idx_ * DIM + DIM + j] = idx(u_coords_);
            }

            for (auto j = 0; j < DIM; j++) {
                auto d_ = coords[j] - 1;
                if (d_ < 0)
                    d_ += dims[j];
                auto d_coords_ = coords;
                d_coords_[j] = d_;

                nn_[2 * idx_ * DIM + j] = idx(d_coords_);
            }
        }

        spdlog::info("Lattice");
        spdlog::info("{}", (void *) odd_.data());

        set_volume();
    }


    size_t idx(dim_t coords) const {
        std::size_t idx = coords[DIM - 1];
        for (int i = DIM - 1; i > 0; i--) {
            idx = idx * dims[i - 1] + coords[i - 1];
        }
        return idx;
    }


    size_t nn(size_t i, int dir, int dim) const { return nn_[2 * DIM * i + (dir + 1) / 2 * DIM + dim]; }
    size_t nn(size_t i, int mu) const { return nn_[2 * DIM * i + mu]; }

    size_t up(size_t i, int dim) const { return nn(i, 1, dim); }
    size_t dn(size_t i, int dim) const { return nn(i, -1, dim); }


    size_t even(size_t i) const { return even_[i]; }

    size_t odd(size_t i) const { return odd_[i]; }

    const dim_t dims;
    std::array<int, DIM + 1> volumes;
    const size_t n_elements;

private:
    void set_volume() {
        volumes[0] = 1;
        for (auto i = 1; i < DIM + 1; ++i) {
            volumes[i] = volumes[i - 1] * dims[i - 1];
        }
    }

    std::vector<size_t> nn_;

    std::vector<size_t> even_;
    std::vector<size_t> odd_;
};
