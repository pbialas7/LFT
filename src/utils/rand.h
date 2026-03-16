//
// Created by pbialas on 16.03.2026.
//

#pragma once

#include <random>

namespace lft::rand {
    template<int N, typename RNG>
    class Generators {
    public:
        std::array<RNG, N> generators;
    };
}
