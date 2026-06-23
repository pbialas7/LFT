//
// Created by pbialas on 06.06.25.
//

#include "fs.h"

#include <filesystem>
#include <format>

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, int Lx, int Ly, const std::string &name,
                        const std::string &extension) {
    auto file_name = std::format("{}_{:03d}x{:03d}_{}.{}", prefix, Lx, Ly, name, extension);
    return data_dir / file_name;
}

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, int Lx, int Ly, int Lz, const std::string &name,
                        const std::string &extension) {
    auto file_name = std::format("{}_{:02d}x{:02d}x{:02d}_{}.{}", prefix, Lx, Ly, Lz, name, extension);
    return data_dir / file_name;
}

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, const std::string &name, const std::string &extension) {
    auto file_name = prefix + "_" + name + "."+extension;
    return data_dir / file_name;
}

std::fstream *optional_fstream_ptr(const fs::path &path, bool create, std::ios::openmode mode) {
    if (create)
        return new std::fstream(path, mode);

    return nullptr;
}
