// RAND/include/rand_logger.hpp
// Logging wrapper for RAND module using spdlog
#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <memory>

namespace randflash {

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
