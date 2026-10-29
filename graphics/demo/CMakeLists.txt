# demo 目录: 以 ctest 验收用例为主 (每个子目录是一个独立测试程序)。
# apps/ 是附加区 —— 现象演示, 不注册 ctest (见 apps/CMakeLists.txt)。

# ---------------------------------------------------------------------------
# 两个注册函数 —— 迁移期两种框架并存, 但**一个用例只准用其中一种**
# ---------------------------------------------------------------------------
# 命名规则保持 1 个可执行文件 = 1 个 ctest 用例, ctest 名 = 文件名。
# 这样 graphics/test-map.md 的「接口 -> 权威用例」映射不用改,
# 也保证 `./run_tests.sh -l` 的输出与文档一致。
# (不用 gtest_discover_tests: 它会把 17 项拆成几百项, 映射表立刻失效。)

# demo/ 自身也要进 include 路径, 供 #include "gtest_common.hpp"
set(DEMO_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# 老写法: 自己带 main 和 Checker 计分板
function(graphics_add_legacy_test name)
    add_executable(${name} ${name}.cpp)
    target_include_directories(${name} PRIVATE ${GRAPHICS_SRC_DIR})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# gtest 写法: 不写 main (GTest::gtest_main 提供), 判据用 TEST()/EXPECT_*
function(graphics_add_gtest name)
    if(NOT GRAPHICS_WITH_GTEST)
        message(FATAL_ERROR
            "${name}.cpp 已经是 gtest 用例, 但 GRAPHICS_WITH_GTEST=OFF。"
            " 去掉该选项或改回 graphics_add_legacy_test。")
    endif()
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE GTest::gtest_main)
    # ${DEMO_DIR} 让用例能写 #include "gtest_common.hpp" 而不必 ../ 相对跳
    target_include_directories(${name} PRIVATE ${GRAPHICS_SRC_DIR} ${DEMO_DIR})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_subdirectory(core)
add_subdirectory(autodiff)
add_subdirectory(linalg)
add_subdirectory(quat)
add_subdirectory(apps)
