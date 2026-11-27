//
// Created by pbialas on 06.06.25.
//

#include "fs.h"

#include <filesystem>
#include <cerrno>
#include <format>
#include <stdexcept>
#include <system_error>

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

std::fstream open_checked(const fs::path &path, std::ios::openmode mode) {
    errno = 0;
    std::fstream stream(path, mode);
    if (!stream.is_open()) {
        auto reason = errno ? std::generic_category().message(errno) : std::string("unknown error");
        throw std::runtime_error(std::format("Cannot open {}: {}", path.string(), reason));
    }
    return stream;
}

void check_stream(const std::ios &stream, const fs::path &path) {
    if (!stream)
        throw std::runtime_error(std::format("Writing to {} failed", path.string()));
}

void close_checked(std::fstream &stream, const fs::path &path) {
    check_stream(stream, path);
    stream.close();
    if (stream.fail())
        throw std::runtime_error(std::format("Closing {} failed (data may be incomplete)", path.string()));
}

std::unique_ptr<std::fstream> optional_fstream_ptr(const fs::path &path, bool create, std::ios::openmode mode) {
    if (create)
        return std::make_unique<std::fstream>(open_checked(path, mode));

    return nullptr;
}
