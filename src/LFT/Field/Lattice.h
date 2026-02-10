//
// Created by pbialas on 20.01.2022.
//


#pragma once

#include <array>
#include <vector>

#include <spdlog/spdlog.h>


namespace {

    template<typename I, std::size_t DIM>
    std::size_t n_elem_(const std::array<I,DIM> &dims) {
        std::size_t n = dims[0];
        for (int i = 1; i < DIM; i++)
            n *= dims[i];
        return n;
    }

    template <typename I, std::size_t DIM>
    void set_strides_(const std::array<I, DIM>& dims, std::array<I, DIM>& strides, char order = 'F') {
        switch (order) {
        case 'F':
            strides[0] = 1;
            for (auto i = 1; i < DIM; ++i) {
                strides[i] = strides[i - 1] * dims[i - 1];
            }
            break;
        case 'C':
            strides[DIM] = 1;
            for (auto i = DIM - 1; i >= 0; --i) {
                strides[i] = strides[i + 1] * dims[i];
            }
            break;
        default:
            spdlog::error("Invalid order: {}, expected 'C' or 'F'", order);
        }
    }
}


namespace lft {
    /*
     * MultiIndex class is used to iterate over the lattice sites and to convert between coordinates and linear index.
     * It is used in the Lattice class to initialize the nearest neighbor indices and to separate even and odd sites.
     *
     *  i_0, i_1, i_2 -> i_0 +dims[0]* i_1 + dims[0]*dims[1]*i_2
     */
    template <typename I, int D>
    class MultiIndex {
    public:
        static const int DIM = D;
        using index_t = I;
        using dims_t = std::array<index_t, DIM>;

        MultiIndex(dims_t dims) : dims(dims) {
            std::fill_n(coords_.begin(), DIM, 0);
            set_volume();
            set_strides_(this->dims, strides_, 'F');
        }

        MultiIndex(dims_t dims, dims_t idx) : dims(dims), coords_(idx) {
            set_volume();
            set_strides_(this->dims, strides_, 'F');
        }

        MultiIndex& operator++() {
            coords_[0]++;
            for (auto i = 0; i < DIM - 1; i++) {
                if (coords_[i] >= dims[i]) {
                    ++coords_[i + 1];
                    coords_[i] = 0;
                }
            }
            return *this;
        }

        dims_t coords() const { return coords_; }

        size_t idx() const {
            index_t idx = 0;
            for (auto i = 0; i < DIM; i++) {
                idx += strides_[i] * coords_[i];
            }
            return idx;
        }

        index_t n_elements() const { return n_elements_; }

        const std::array<index_t, DIM> dims;


        auto strides() const { return strides_; }

    private:
        void set_volume() {
            n_elements_ = 1;
            for (auto i = 0; i < DIM; ++i) {
                n_elements_ *= dims[i];
            }
        }


        std::array<index_t, DIM> coords_;
        std::array<index_t, DIM> strides_;
        uint64_t n_elements_;
    };

    /*
     * Lattice class represents a D-dimensional lattice with periodic boundary conditions. It stores the dimensions of the lattice,
     * the number of elements, the nearest neighbor indices, and the even and odd sites. The nearest neighbor indices are stored in a flat vector,
     * where the first half corresponds to the negative direction and the second half corresponds to the positive direction
     * The even and odd sites are separated based on the sum of their coordinates, which is used for certain algorithms that require a checkerboard pattern.
     * The class provides methods to get the nearest neighbor indices, the even and odd sites, and
     * to convert between coordinates and linear index.
     */
    template <typename I=uint32_t, int D = 2>
    class Lattice {
    public:
        static const int DIM = D;
        using dim_t = std::array<I, DIM>;
        using index_t = I;
        using size_t = std::size_t;

    private:
        static std::size_t n_elem(dim_t dims) {
            std::size_t n = dims[0];
            for (int i = 1; i < DIM; i++)
                n *= dims[i];
            return n;
        }

    public:
        Lattice(std::array<I, DIM> dims) : dims(dims), n_elements(n_elem(dims)),
                                           nn_(2 * DIM * n_elements) {
            MultiIndex<I, DIM> m_index(dims);

            set_strides_(dims, strides_, 'F');

            for (auto i = 0; i < m_index.n_elements(); i++, ++m_index) {
                auto idx_ = i;
                auto coords = m_index.coords();
                size_t i_sum = 0;
                for (int j = 0; j < DIM; j++)
                    i_sum += coords[j];

                if (i_sum % 2 == 0) {
                    even_.push_back(idx_);
                }
                else {
                    odd_.push_back(idx_);
                }

                for (auto j = 0; j < DIM; j++) {
                    I u_;
                    if (coords[j] == dims[j] - 1)
                        u_ = 0;
                    else
                        u_ = coords[j] + 1;

                    auto u_coords_ = coords;
                    u_coords_[j] = u_;

                    nn_[2 * idx_ * DIM + DIM + j] = idx(u_coords_);
                }

                for (auto j = 0; j < DIM; j++) {
                    I d_;
                    if (coords[j] == 0)
                        d_ = dims[j] - 1;
                    else
                        d_ = coords[j] - 1;
                    auto d_coords_ = coords;
                    d_coords_[j] = d_;

                    nn_[2 * idx_ * DIM + j] = idx(d_coords_);
                }
            }
        }


        size_t idx(dim_t coords) const {
            std::size_t idx = 0;
            for (int i = 0; i < DIM; ++i) {
                idx += coords[i] * strides_[i];
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
        const size_t n_elements;

    private:
        dim_t strides_;

        std::vector<size_t> nn_;

        std::vector<size_t> even_;
        std::vector<size_t> odd_;
    };
}
