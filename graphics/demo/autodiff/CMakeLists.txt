# 计算图自动微分测试
add_executable(compute_graph compute_graph.cpp)

# 依赖 src/autodiff 与 src/core 的头文件
target_include_directories(compute_graph PRIVATE
    ${GRAPHICS_SRC_DIR}
)

# 注册为 ctest 测试用例
add_test(NAME compute_graph COMMAND compute_graph)

# Apply 的域分裂验收 (2026-10-16): 4 个重载 x 4 个域 + Function/MatrixFunction 节点
# + 嵌套对偶 Dual<Dual<T>> 的二阶导 (修 dual.hpp 的 log ADL 断链)
add_executable(apply apply.cpp)
target_include_directories(apply PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME apply COMMAND apply)

# 反向模式最小可用验收 (2026-10-17): 最小闭环三方对拍 + 12 个可用算子的伴随公式
# + adjoint 自动初始化 + 幂等; 占位算子如实断言"仍不可用"
add_executable(backward backward.cpp)
target_include_directories(backward PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME backward COMMAND backward)

# 分子布局雅可比 (D19 第二次修订, 2026-10-19): 形状 m×n + 链式法则自然序
# + 与图反向逐元对拍 + float/double 双域 + 非法输入
add_executable(jacobian jacobian.cpp)
target_include_directories(jacobian PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME jacobian COMMAND jacobian)
