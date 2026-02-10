cmake_minimum_required(VERSION 3.15)
project(test LANGUAGES CXX)

add_executable(test test.cpp
        test_Lattice.cpp
        test_multiindex.cpp
        test_Field.cpp
        test_ising.cpp
        test_sweep.cpp
)
target_link_libraries(test PRIVATE lft spdlog::spdlog Catch2::Catch2)