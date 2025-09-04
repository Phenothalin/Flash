cmake_minimum_required(VERSION 3.18)
project(rand_flash_tp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- 找 Eigen (你已经在代码里使用了) ---
find_package(Eigen3 REQUIRED)

# --- 找 ThermoPack（假设已安装并导出 thermopack target）---
# 典型做法：cmake 安装后会生成 thermopack-config.cmake
# 你可以通过 -DThermoPack_DIR=/path/to/thermopack/cmake 指定
find_package(ThermoPack REQUIRED CONFIG)

# --- GTest（用于 tests）---
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# --- 头文件路径 ---
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# --- 库目标：rand_flash + ThermoPackAdapter ---
add_library(rand_flash
  src/rand_flash.cpp
  src/ThermoPackAdapter.cpp
)
target_link_libraries(rand_flash
  PUBLIC Eigen3::Eigen
         thermopack   # 来自 ThermoPack 的 imported target
)
target_include_directories(rand_flash
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# --- 测试 ---
add_executable(rand_flash_test tests/rand_flash_test.cpp)
target_link_libraries(rand_flash_test
  PRIVATE rand_flash
          GTest::gtest GTest::gtest_main
)
enable_testing()
add_test(NAME RandFlashTest COMMAND rand_flash_test)
