//
// Created by pbialas on 06.06.25.
//

#pragma once

#include <filesystem>
namespace fs = std::filesystem;
#include <ostream>
#include <fstream>


fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, int Lx, int Ly, const std::string &name,
                        const std::string &extension);

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, const std::string &name,
                        const std::string &extension);


std::fstream *optional_fstream_ptr(const fs::path &path, bool create, std::ios::openmode mode);
