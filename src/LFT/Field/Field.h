//
// Created by pbialas on 20.01.2022.
//

#pragma once


#include <vector>
#include <ostream>

#include "Field/Lattice.h"


template <typename F, typename L>
class Field {
public:
    using field_t = F;
    using lattice_t = L;
    static const int DIM = lattice_t::DIM;

    const size_t n_elements;
    const lattice_t lat;

    explicit Field(const lattice_t& lat) : lat(lat), n_elements(lat.n_elements), field_(lat.n_elements) {
    }

    Field(const lattice_t& lat, field_t val) : lat(lat), n_elements(lat.n_elements), field_(lat.n_elements, val) {
    }

    field_t& operator[](size_t i) { return field_[i]; }

    field_t operator[](size_t i) const { return field_[i]; }

    field_t& operator[](typename lattice_t::dim_t dims) {
        return field_[lat.idx(dims)];
    }

    field_t operator[](typename lattice_t::dim_t dims) const {
        return field_[lat.idx(dims)];
    }


    field_t up_corona(typename lattice_t::index_t i) const {
        field_t corona_ = field_t(0);
        for (auto d = 0; d < DIM; ++d) {
            corona_ += field_[lat.up(i, d)];
        }
        return corona_;
    }

    field_t dn_corona(typename lattice_t::index_t i) const {
        field_t corona_ = field_t(0);
        for (auto d = 0; d < DIM; ++d) {
            corona_ += field_[lat.dn(i, d)];
        }
        return corona_;
    }

    field_t corona(typename lattice_t::index_t i) const {
        return up_corona(i) + dn_corona(i);
    }

    const field_t* data() const {
        return field_.data();
    }

    void fill(field_t val) {
        std::fill(field_.begin(), field_.end(), val);
    }

    void write(std::ostream& out) const {
        out.write((const char*)field_.data(), sizeof(field_t) * n_elements);
    }


    friend std::ostream& operator<<(std::ostream& out, const Field& field) {
        for (auto d = 0; d < field.n_elements; ++d) {
            out << field[d] << " ";
        }
        return out;
    }

private:
    std::vector<field_t> field_;
};

template <typename F, typename L>
Field<F, L> make_field(const L& lat) {
    return Field<F, L>(lat);
}

template <typename F, typename L>
Field<F, L> make_field(const L& lat, F val) {
    return Field<F, L>(lat, val);
}


template <typename Res, typename F>
Res sum(const F& field) {
    Res sum = Res(0);
    for (auto i = 0; i < field.n_elements; ++i)
        sum += field[i];
    return sum;
}

template <typename Float, typename F>
Float mean(const F& field) {
    return sum<Float, F>(field) / field.n_elements;
}
