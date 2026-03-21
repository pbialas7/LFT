cmake_minimum_required(VERSION 3.15)
project(EdwardsAnderson LANGUAGES CXX)

add_executable(ea_q ea_q.cpp ea.h
        ../utils/fs.h
        ../utils/fs.cpp
        ../utils/hardware.cpp
        ../utils/hardware.h
        options.cpp
        options.h
)
target_link_libraries(ea_q PRIVATE yaml-cpp::yaml-cpp lft lyra spdlog::spdlog OpenMP::OpenMP_CXX)
target_include_directories(ea_q PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_options(ea_q PRIVATE -g -fno-omit-frame-pointer)
