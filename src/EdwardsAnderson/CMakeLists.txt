cmake_minimum_required(VERSION 3.15)
project(EdwardsAnderson LANGUAGES CXX)



add_executable(ea ea_main.cpp ea.h
        ../utils/fs.h
        ../utils/fs.cpp
        ../utils/hardware.cpp
        ../utils/hardware.h
        options.cpp
        options.h
        ../utils/logging.cpp
        ../utils/logging.h
)
target_link_libraries(ea PRIVATE yaml-cpp::yaml-cpp lft lyra spdlog::spdlog OpenMP::OpenMP_CXX)
target_include_directories(ea PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_options(ea PRIVATE -g -fno-omit-frame-pointer)
message(${CMAKE_INSTALL_PREFIX})


