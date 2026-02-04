// RAND/include/rand_logger.hpp
// Logging wrapper for RAND module using spdlog
#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <memory>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX  // Prevent Windows.h from defining min/max macros
#endif
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#endif

namespace randflash {

// Initialize console for UTF-8 output (Windows only)
inline void initConsole() {
#ifdef _WIN32
    // Set console output to UTF-8
    SetConsoleOutputCP(CP_UTF8);
    // Set console input to UTF-8
    SetConsoleCP(CP_UTF8);
    // Enable buffering for better performance
    setvbuf(stdout, nullptr, _IOFBF, 1000);
    setvbuf(stderr, nullptr, _IOFBF, 1000);
#endif
}

// Auto-initialize console on module load
namespace detail {
    struct ConsoleInitializer {
        ConsoleInitializer() { initConsole(); }
    };
    static ConsoleInitializer console_init;
}

// Get the RAND module logger
inline std::shared_ptr<spdlog::logger>& getLogger() {
    static std::shared_ptr<spdlog::logger> logger = spdlog::default_logger();
    return logger;
}

// Set custom logger (optional)
inline void setLogger(std::shared_ptr<spdlog::logger> logger) {
    getLogger() = std::move(logger);
}

// Set log level
inline void setLogLevel(spdlog::level::level_enum level) {
    getLogger()->set_level(level);
}

// Convenience macros
#define RAND_TRACE(...)    randflash::getLogger()->trace(__VA_ARGS__)
#define RAND_DEBUG(...)    randflash::getLogger()->debug(__VA_ARGS__)
#define RAND_INFO(...)     randflash::getLogger()->info(__VA_ARGS__)
#define RAND_WARN(...)     randflash::getLogger()->warn(__VA_ARGS__)
#define RAND_ERROR(...)    randflash::getLogger()->error(__VA_ARGS__)

} // namespace randflash
