cmake_minimum_required(VERSION 3.16)
project(graphics CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 源码目录: matrix/dual/grad_graph 等都是头文件模板, 测试只需 include
set(GRAPHICS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ---------------------------------------------------------------------------
# 第三方依赖: googletest
# ---------------------------------------------------------------------------
# **不在构建期联网**。FetchContent 会让 cmake 的成败取决于网络,
# 且同一份源码在不同机器上可能编出不同东西 (见 graphics/deps.manifest 头注释)。
#
#   ./graphics/scripts/fetch_deps.sh        # 先跑这个, 把源码放到 third_party/googletest
#
# 找不到就 **直接报错**, 不静默跳过 —— 静默跳过会让闸门悄悄消失
# (本仓有过一次: 整个 so3.hpp 没登记进 smoke, 用例全绿但一次都没测到)。
option(GRAPHICS_WITH_GTEST "启用 googletest 用例 (需先跑 fetch_deps.sh)" ON)

if(GRAPHICS_WITH_GTEST)
    set(GTEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest)
    if(NOT EXISTS ${GTEST_DIR}/CMakeLists.txt)
        message(FATAL_ERROR
            "\n找不到 googletest 源码: ${GTEST_DIR}/CMakeLists.txt\n"
            "\n先跑一次性依赖脚本 (它会下载并校验 sha256):\n"
            "    ./graphics/scripts/fetch_deps.sh\n"
            "\n确实要绕过 (例如只跑老的手写用例):\n"
            "    cmake -S graphics -B graphics/build -DGRAPHICS_WITH_GTEST=OFF\n")
    endif()

    # gtest 自身的构建选项: 不装、不建 gmock、不要它自己的用例
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(gtest_build_samples OFF CACHE BOOL "" FORCE)
    set(gtest_build_tests OFF CACHE BOOL "" FORCE)
    add_subdirectory(${GTEST_DIR} EXCLUDE_FROM_ALL)
    message(STATUS "googletest: ${GTEST_DIR}")
endif()

# demo 目录存放测试程序, 通过 ctest 运行
enable_testing()
add_subdirectory(demo)
