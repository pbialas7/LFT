cmake_minimum_required(VERSION 3.15)
project(main LANGUAGES CXX)

include(FetchContent)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
)
FetchContent_MakeAvailable(spdlog)


FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/pbialas7/fmt.git
        GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG devel
)
FetchContent_MakeAvailable(Catch2)

FetchContent_Declare(
        lyra
        GIT_REPOSITORY https://github.com/bfgroup/Lyra.git
        GIT_TAG develop

)
FetchContent_MakeAvailable(lyra)

add_subdirectory(src/LFT)
add_subdirectory(tests)

add_executable(ising src/Drivers/ising.cpp)
target_link_libraries(ising PRIVATE lft lyra spdlog::spdlog)
#set_target_properties(ising PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(ising1d src/Drivers/ising1d.cpp)
target_link_libraries(ising1d PRIVATE lft lyra spdlog::spdlog)

add_executable(phi4 src/Drivers/phi4.cpp)
target_link_libraries(phi4 PRIVATE lft lyra spdlog::spdlog fmt::fmt)