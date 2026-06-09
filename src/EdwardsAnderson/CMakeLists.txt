cmake_minimum_required(VERSION 3.15)
project(EdwardsAnderson LANGUAGES CXX)


set(SRC_COMMON
        ea.h
        ../utils/fs.h
        ../utils/fs.cpp
        ../utils/hardware.cpp
        ../utils/hardware.h
        ../utils/logging.cpp
        ../utils/logging.h)

add_executable(ea ea_main.cpp ea.h
        parallel_tempering.h
        parallel_tempering.cpp
        options.cpp
        ${SRC_COMMON}
)
target_link_libraries(ea PRIVATE yaml-cpp::yaml-cpp lft lyra spdlog::spdlog OpenMP::OpenMP_CXX)
target_include_directories(ea PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_options(ea PRIVATE -g -fno-omit-frame-pointer)


add_executable(ea_pt_mt ea_main_pt_mt.cpp
        parallel_tempering_mt.h
        parallel_tempering_mt.cpp
        options_cli11.cpp
        ${SRC_COMMON}
)
target_link_libraries(ea_pt_mt PRIVATE yaml-cpp::yaml-cpp lft spdlog::spdlog OpenMP::OpenMP_CXX CLI11::CLI11)
target_include_directories(ea_pt_mt PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_options(ea_pt_mt PRIVATE -g -fno-omit-frame-pointer)