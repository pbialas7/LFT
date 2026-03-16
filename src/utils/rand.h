//
// Created by pbialas on 16.03.2026.
//

#pragma once

#include <limits>
#include <iostream>
#include <sstream>
#include <vector>

using Float = double;

namespace lft::rand {
    template<typename S, int SEEDS_PER_GEN, int ALIGNMENT_IN_BYTES>
    class random_array_state {
    public:
        typedef S seed_t;

        static const int seeds_per_gen = SEEDS_PER_GEN;

        random_array_state(int n_gen)
            : n_gen_(n_gen),
              SEEDS_BLOCK_SIZE(ALIGNMENT_IN_BYTES *
                               ((sizeof(S) * SEEDS_PER_GEN + ALIGNMENT_IN_BYTES - 1) /
                                ALIGNMENT_IN_BYTES)),
              SEEDS_BLOCK(SEEDS_BLOCK_SIZE / sizeof(S)) {
            posix_memalign(
                (void **) &seeds_, ALIGNMENT_IN_BYTES, SEEDS_BLOCK_SIZE * n_gen_);
        }

        int n_gen() const { return n_gen_; }

        void set_seeds(const S *seeds) {
            for (int i = 0; i < n_gen_; ++i) {
                for (int j = 0; j < SEEDS_PER_GEN; ++j)
                    seeds_[i * SEEDS_BLOCK_SIZE / sizeof(S) + j] =
                            seeds[i * SEEDS_PER_GEN + j];
            }
        }

        void get_seeds(S *seeds) const {
            for (int i = 0; i < n_gen_; ++i) {
                for (int j = 0; j < SEEDS_PER_GEN; ++j)
                    seeds[i * SEEDS_PER_GEN + j] =
                            seeds_[i * SEEDS_BLOCK_SIZE / sizeof(S) + j];
            }
        }

        void gen_seeds(long seed) {
            srand48(seed);
            for (int i = 0; i < n_gen_; ++i) {
                for (int j = 0; j < SEEDS_PER_GEN; ++j)
                    seeds_[i * SEEDS_BLOCK_SIZE / sizeof(S) + j] =
                            lrand48() % (std::numeric_limits<S>::max());
            }
        }

        int fwrite_state(FILE *fout) {
            S *seeds = new S[SEEDS_PER_GEN * n_gen_];
            get_seeds(seeds);
            int bytes;
            const int n_params = 3;
            int params[n_params] = {n_gen_, SEEDS_PER_GEN, sizeof(S)};
            bytes = std::fwrite(params, sizeof(int), n_params, fout);
            bytes += std::fwrite(seeds, sizeof(S), n_gen_ * SEEDS_PER_GEN, fout);
            delete[] seeds;
            return bytes;
        }

        int fread_state(FILE *fin) {
            const int n_params = 3;
            int params[n_params];
            int bytes = std::fread(params, sizeof(int), n_params, fin);
            if (params[0] != n_gen_) {
                std::ostringstream msg;
                msg << "number of generators differ : program " << n_gen_ << " read "
                        << params[0];
                throw(msg);
            }
            if (params[1] != SEEDS_PER_GEN) {
                std::ostringstream msg;
                msg << "SEES_PER_GEN differ : program " << SEEDS_PER_GEN << " read "
                        << params[1];
                throw(msg);
            }
            if (params[2] != sizeof(S)) {
                std::ostringstream msg;
                std::cerr << "sizeof seeds differ : program " << sizeof(S) << " read "
                        << params[2];
                throw(msg);
            }
            S *seeds = new S[SEEDS_PER_GEN * n_gen_];
            bytes += std::fread(seeds, sizeof(S), n_gen_ * SEEDS_PER_GEN, fin);
            set_seeds(seeds);
            delete[] seeds;
            return bytes;
        }

    protected:
        const int SEEDS_BLOCK_SIZE;
        const int SEEDS_BLOCK;
        int n_gen_;
        S *seeds_;
    };


    class taus_array : public random_array_state<unsigned int, 4, 256> {
    public:
        using result_type = unsigned int;

        static constexpr result_type min() {
            return std::numeric_limits<result_type>::min();
        }

        static constexpr result_type max() {
            return std::numeric_limits<result_type>::max();
        }

        taus_array(int n) : random_array_state<unsigned int, 4, 256>(n) {
            gen_seeds(121245);
            taus_.reserve(n);
            for (int i = 0; i < n; ++i) {
                taus_.emplace_back(i, this);
            }
        }

        result_type operator()(int i) {
            taus_step(seeds_[i * SEEDS_BLOCK], 13, 19, 12, 4294967294u);
            taus_step(seeds_[i * SEEDS_BLOCK + 1], 2, 25, 4, 4294967288u);
            taus_step(seeds_[i * SEEDS_BLOCK + 2], 3, 11, 17, 4294967280u);
            LCG_step(seeds_[i * SEEDS_BLOCK + 3], 1664525u, 1013904223u);
            return (seeds_[i * SEEDS_BLOCK] ^ seeds_[i * SEEDS_BLOCK + 1] ^
                    seeds_[i * SEEDS_BLOCK + 2] ^ seeds_[i * SEEDS_BLOCK + 3]);
        }


        static taus_array *generator() { return generator_; }

        class taus {
        public:
            taus(int i, taus_array *array) : i_(i), array_(array) {
            }

            using result_type = taus_array::result_type;

            static constexpr result_type min() {
                return taus_array::min();
            }

            static constexpr result_type max() {
                return taus_array::max();
            }

            result_type operator()() {
                return array_->operator()(i_);
            }

        private:
            int i_;
            taus_array *array_;
        };

        static void init(int n, long int seed) {
            generator_ = new taus_array(n);
            generator_->gen_seeds(seed);
        }


        taus &operator[](int i) {
            return taus_[i];
        }

    private:
        void taus_step(unsigned &z, int S1, int S2, int S3, unsigned M) {
            unsigned b = (((z << S1) ^ z) >> S2);
            z = (((z & M) << S3) ^ b);
        }

        void LCG_step(unsigned &z, unsigned A, unsigned C) { z = (A * z + C); }

        static taus_array *generator_;
        std::pmr::vector<taus> taus_;
    };
}
