cmake_minimum_required(VERSION 3.18)

message(STATUS "Configuring embedded ThermoPack...")

set(THERMOPACK_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(THERMOPACK_INCLUDE ${THERMOPACK_ROOT}/include)

# 平台检测
if(APPLE)
    set(THERMOPACK_LIB_DIR ${THERMOPACK_ROOT}/lib/macos)
    set(THERMOPACK_SHARED_LIB ${THERMOPACK_LIB_DIR}/libthermopack.dylib)
    set(THERMOPACK_STATIC_LIB ${THERMOPACK_LIB_DIR}/libthermopack.a)
elseif(WIN32)
    set(THERMOPACK_LIB_DIR ${THERMOPACK_ROOT}/lib/windows)
    set(THERMOPACK_SHARED_LIB ${THERMOPACK_LIB_DIR}/libthermopack.dll)
    set(THERMOPACK_STATIC_LIB ${THERMOPACK_LIB_DIR}/libthermopack.a)
    set(THERMOPACK_IMPLIB ${THERMOPACK_LIB_DIR}/libthermopack.dll.a)
elseif(UNIX)
    set(THERMOPACK_LIB_DIR ${THERMOPACK_ROOT}/lib/linux)
    set(THERMOPACK_SHARED_LIB ${THERMOPACK_LIB_DIR}/libthermopack.so)
    set(THERMOPACK_STATIC_LIB ${THERMOPACK_LIB_DIR}/libthermopack.a)
endif()

# 静态/动态链接选项
option(THERMOPACK_USE_STATIC "Link ThermoPack statically" OFF)

if(THERMOPACK_USE_STATIC)
    set(THERMOPACK_LIB ${THERMOPACK_STATIC_LIB})
    set(THERMOPACK_LIB_TYPE STATIC)
else()
    set(THERMOPACK_LIB ${THERMOPACK_SHARED_LIB})
    set(THERMOPACK_LIB_TYPE SHARED)
endif()

# 创建 imported target
add_library(thermopack ${THERMOPACK_LIB_TYPE} IMPORTED GLOBAL)
set_target_properties(thermopack PROPERTIES
    IMPORTED_LOCATION ${THERMOPACK_LIB}
    INTERFACE_INCLUDE_DIRECTORIES ${THERMOPACK_INCLUDE}
)

# Windows: 设置导入库
if(WIN32 AND NOT THERMOPACK_USE_STATIC)
    set_target_properties(thermopack PROPERTIES
        IMPORTED_IMPLIB ${THERMOPACK_IMPLIB}
    )
endif()

# macOS: 链接 Accelerate framework (LAPACK/BLAS)
if(APPLE)
    find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    set_property(TARGET thermopack APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES ${ACCELERATE_FRAMEWORK})
endif()

message(STATUS "ThermoPack library: ${THERMOPACK_LIB}")
message(STATUS "ThermoPack include: ${THERMOPACK_INCLUDE}")
