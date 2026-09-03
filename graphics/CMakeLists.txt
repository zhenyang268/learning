cmake_minimum_required(VERSION 3.16)
project(graphics CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 源码目录: matrix/dual/grad_graph 等都是头文件模板, 测试只需 include
set(GRAPHICS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# demo 目录存放测试程序, 通过 ctest 运行
enable_testing()
add_subdirectory(demo)
