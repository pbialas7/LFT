//
// Created by pbialas on 16.06.2026.
//

#pragma once

#include <filesystem>
namespace fs = std::filesystem;
#include <fstream>

namespace lft::ea {
    template<typename Float, typename F>
    void write_gsg_file(const F &field1, const F &field2, fs::path &path) {
        auto lat = field1.lat;
        auto out = std::fstream(path, std::fstream::out);
        for (int i = 0; i < lat.n_elem; i++) {
        }
    }
}
