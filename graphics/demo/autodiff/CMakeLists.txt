# 计算图自动微分测试
#
# 框架: gtest (2026-10-30 迁移) —— 判据名走 SCOPED_TRACE, 迁移前 36 项断言全保留
graphics_add_gtest(compute_graph)

# Apply 的域分裂验收 (2026-10-16): 4 个重载 x 4 个域 + Function/MatrixFunction 节点
# + 嵌套对偶 Dual<Dual<T>> 的二阶导 (修 dual.hpp 的 log ADL 断链)
# 框架: gtest (2026-10-30 迁移) —— 判据编号 [A1]–[F3] 保留
graphics_add_gtest(apply)

# 反向模式最小可用验收 (2026-10-17): 最小闭环三方对拍 + 12 个可用算子的伴随公式
# + adjoint 自动初始化 + 幂等; 占位算子如实断言"仍不可用"
# 框架: gtest (2026-10-30 迁移) —— [F] 分域由 TYPED_TEST_SUITE 保证
graphics_add_gtest(backward)

# 分子布局雅可比 (D19 第二次修订, 2026-10-19): 形状 m×n + 链式法则自然序
# + 与图反向逐元对拍 + float/double 双域 + 非法输入
# 框架: gtest (2026-10-30 迁移) —— double/float 双域由 TYPED_TEST_SUITE 保证
graphics_add_gtest(jacobian)