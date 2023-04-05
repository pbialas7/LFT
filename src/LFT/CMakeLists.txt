cmake_minimum_required(VERSION 3.15)
project(LFT LANGUAGES CXX)

add_library(lft INTERFACE Ising/correlation_function.h)
target_include_directories(lft INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)


