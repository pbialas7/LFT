//
// Created by pbialas on 06.06.25.
//

#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
namespace fs = std::filesystem;
#include <ostream>
#include <fstream>


fs::path make_file_path(const fs::path& data_dir, const std::string& prefix, int Lx, int Ly, const std::string& name,
                        const std::string& extension);

fs::path make_file_path(const fs::path& data_dir, const std::string& prefix, int Lx, int Ly, int Lz,
                        const std::string& name,
                        const std::string& extension);

fs::path make_file_path(const fs::path& data_dir, const std::string& prefix, const std::string& name,
                        const std::string& extension);


// File name for a run on a lattice with extents `dims`: {prefix}_{L0:03}x{L1:03}_{name}.{ext} in 2D and
// {prefix}_{L0:02}x{L1:02}x{L2:02}_{name}.{ext} in 3D.
template <std::size_t N>
fs::path make_file_path(const fs::path& data_dir, const std::string& prefix, const std::array<std::uint32_t, N>& dims,
                        const std::string& name, const std::string& extension) {
    static_assert(N == 2 || N == 3, "make_file_path supports 2D and 3D lattices");
    if constexpr (N == 2)
        return make_file_path(data_dir, prefix, static_cast<int>(dims[0]), static_cast<int>(dims[1]), name, extension);
    else
        return make_file_path(data_dir, prefix, static_cast<int>(dims[0]), static_cast<int>(dims[1]),
                              static_cast<int>(dims[2]), name, extension);
}


std::fstream* optional_fstream_ptr(const fs::path& path, bool create, std::ios::openmode mode);
