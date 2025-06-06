//
// Created by pbialas on 06.06.25.
//

#pragma once

#include <filesystem>

namespace fs = std::filesystem;

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, const std::string &name,
                        const std::string &extension);
