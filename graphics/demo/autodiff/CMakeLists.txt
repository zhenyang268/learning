# 计算图自动微分测试
add_executable(compute_graph compute_graph.cpp)

# 依赖 src/autodiff 与 src/core 的头文件
target_include_directories(compute_graph PRIVATE
    ${GRAPHICS_SRC_DIR}
)

# 注册为 ctest 测试用例
add_test(NAME compute_graph COMMAND compute_graph)
