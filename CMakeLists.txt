cmake_minimum_required(VERSION 3.15)
project(LFT LANGUAGES CXX)

include(FetchContent)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 20)

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(spdlog)


FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG devel
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

FetchContent_Declare(
        lyra
        GIT_REPOSITORY https://github.com/bfgroup/Lyra.git
        GIT_TAG develop
        GIT_SHALLOW TRUE

)
FetchContent_MakeAvailable(lyra)

FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm.git
        GIT_TAG master
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(glm)


add_subdirectory(src/LFT)
add_subdirectory(tests)

add_executable(ising src/Drivers/ising.cpp src/LFT/Ising/correlation_function.h
        src/LFT/Ising/Wolff.h
        src/utils/fs.h
        src/utils/fs.cpp
        src/utils/extract.cpp
        src/utils/extract.h
        src/Drivers/ising_base_options.cpp
        src/Drivers/ising_base_options.h)

target_link_libraries(ising PRIVATE lft lyra spdlog::spdlog)
target_include_directories(ising PRIVATE src)
target_compile_options(ising PRIVATE
        $<$<CONFIG:Release>:-O3 -DNDEBUG>
)

#set_target_properties(ising PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(ising1d src/Drivers/ising1d.cpp)
target_link_libraries(ising1d PRIVATE lft lyra spdlog::spdlog)

add_executable(ising3d src/Drivers/ising3d.cpp)
target_link_libraries(ising3d PRIVATE lft lyra spdlog::spdlog)

add_executable(phi4 src/Drivers/phi4.cpp)
target_link_libraries(phi4 PRIVATE lft lyra spdlog::spdlog)

add_executable(xy src/Drivers/xy.cpp src/LFT/XY/xy.h)
target_link_libraries(xy PRIVATE lft lyra spdlog::spdlog)

add_executable(ea src/Drivers/ea.cpp
        src/LFT/EdwardsAnderson/ea.h
        src/Drivers/ising_base_options.cpp
        src/Drivers/ising_base_options.h
        src/utils/fs.h
        src/utils/fs.cpp
)
target_link_libraries(ea PRIVATE lft lyra spdlog::spdlog)
target_include_directories(ea PRIVATE src)

add_custom_target(all_targets DEPENDS ising ising1d ising3d phi4 xy ea)