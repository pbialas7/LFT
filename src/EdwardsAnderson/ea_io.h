//
// Helpers shared by the EA drivers: initialisation of the couplings and spins, and the per-sweep measurement.
//

#pragma once

#include <format>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

#include <spdlog/spdlog.h>

#include "ea.h"
#include "replicas.h"

namespace lft::ea {
    /**
     * Initializes the field with either random values or from a file. If file_path is empty, it initializes randomly.
     * If ising is true, it initializes with 1. If binary is true, it initializes with +/-1.
     * Otherwise, it initializes with Gaussian random numbers. A file must contain exactly DIM * N values;
     * otherwise std::runtime_error is thrown.
     */
    template <typename Field, typename RNG>
    void init_j_field(Field& field, const std::string& file_path, bool binary, bool ising, RNG& rng) {
        if (file_path.empty()) {
            if (!ising) {
                if (binary)
                    init_bernoulli(field, rng);
                else
                    init_gaussian(field, rng);
            }
        }
        else {
            std::ifstream ifs(file_path, std::ios::in);
            if (!ifs)
                throw std::runtime_error(std::format("Cannot open the coupling file {}", file_path));
            spdlog::info("Reading link variables from {}", file_path);
            const std::size_t n = field.lat.n_elements;
            for (std::size_t i = 0; i < n; ++i) {
                float val;
                if (!(ifs >> val)) {
                    throw std::runtime_error(
                        ifs.eof()
                            ? std::format("Coupling file {} has only {} values, but this lattice needs {}",
                                          file_path, i, n)
                            : std::format("Coupling file {}: value {} is not a number", file_path, i));
                }
                field[i] = val;
            }
            // The file must not contain anything else: extra values mean it was written for another lattice.
            ifs >> std::ws;
            if (!ifs.eof())
                throw std::runtime_error(std::format(
                    "Coupling file {} has more than the {} values this lattice needs (written for a larger "
                    "lattice or another dimension?)", file_path, n));
        }
    }

    /**
     * Writes the couplings as text, one value per line, in the order init_j_field reads them. Uses enough digits
     * (max_digits10) that reading the file back gives exactly the same float values.
     */
    template <typename Field>
    void write_j_field(std::ostream& out, const Field& field) {
        out << std::setprecision(std::numeric_limits<typename Field::field_t>::max_digits10) << field;
    }

    template <typename L, typename RNG>
    void init_spin_field(SpinField<L>& field, bool cold_start, RNG& rng) {
        if (cold_start) {
            field.fill(1);
        }
        else {
            init_bernoulli(field, rng);
        }
    }

    /**
     * Writes energy and magnetisation of every replica and, for two or more replicas, the overlap and link overlap
     * of replicas 0 and 1. Does nothing if em_stream_ptr is null.
     */
    template <typename L>
    void measure_em(std::fstream* em_stream_ptr, Replicas<L>& replicas, const JField<float, L>& j_field) {
        if (em_stream_ptr) {
            for (int j = 0; j < replicas.n_replicas; ++j) {
                *em_stream_ptr << energy<double>(*replicas[j], j_field) << " ";
                *em_stream_ptr << magnetisation<double>(*replicas[j]) << " ";
            }
            if (replicas.n_replicas > 1) {
                *em_stream_ptr << overlap<double>(*replicas[0], *replicas[1]) << " ";
                *em_stream_ptr << link_overlap<double>(*replicas[0], *replicas[1]) << "\n";
            }
            else
                *em_stream_ptr << "\n";
            em_stream_ptr->flush();
        }
    }
}
