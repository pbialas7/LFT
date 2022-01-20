//
// Created by pbialas on 20.01.2022.
//

#include "catch2/catch_session.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"


int main( int argc, char* argv[] ) {
    //    your setup ...
    //    auto test_console = spdlog::stdout_color_mt("test_logger");
    //    test_console->set_level(spdlog::level::debug);
    int result = Catch::Session().run( argc, argv );

    // your clean-up...

    return result;
}