# demo/apps —— 附加区: 现象演示 / 实验程序。不属于 ctest 验收体系。
#
# 与 demo/<模块>/ 的分工 (一条判定规则):
#   要出现在 ctest 的绿灯里  ->  demo/<模块>/, 四行 (add_executable + add_test)
#   只是看现象、不判对错    ->  这里,     三行 (add_executable + include, 无 add_test)
#
# 约定: 文件名 = CMake 目标名 —— run_apps.sh 靠扫目录发现, 与这里的 target 必须一致。
#
# 一律 EXCLUDE_FROM_ALL 的理由:
#   run_tests.sh 的构建是 `cmake --build $BUILD_DIR -j` (不带 target, 走 ALL)。
#   不用 EXCLUDE_FROM_ALL 的话, 本目录里的半成品/实验代码会**编译失败**,
#   进而把整个 run_tests.sh 卡在构建阶段 (exit 1), ctest 一条都跑不到。
#   加上之后: 默认构建不碰它们, 验收绿灯不受附加区影响。
#
# 构建 / 运行:
#   cmake --build graphics/build --target <name> -j     (单独构建)
#   ./graphics/run_apps.sh <name>                       (构建 + 运行, 推荐)

add_executable(string_modes EXCLUDE_FROM_ALL string_modes.cpp)
target_include_directories(string_modes PRIVATE
    ${GRAPHICS_SRC_DIR}
)
