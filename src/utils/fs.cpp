//
// Created by pbialas on 06.06.25.
//

#include "fs.h"

#include <filesystem>

fs::path make_file_path(const fs::path &data_dir, const std::string &prefix, const std::string &name, const std::string &extension) {
    auto file_name = prefix + "_" + name + "."+extension;
    return data_dir / file_name;
}