# linalg 层验收 —— 每个用例对应一个交付判据
# 公共测试工具在 demo/common.hpp
# 框架: gtest (2026-10-30 迁移) —— 判据名保留原样, 走 SCOPED_TRACE
graphics_add_gtest(solve)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[H] 走 SCOPED_TRACE
graphics_add_gtest(inverse)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[I] 走 SCOPED_TRACE
graphics_add_gtest(cholesky)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[E] 走 SCOPED_TRACE
graphics_add_gtest(eigen)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[F] 走 SCOPED_TRACE
graphics_add_gtest(svd)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[E] 走 SCOPED_TRACE
graphics_add_gtest(qr)

# 框架: gtest (2026-10-30 迁移) —— 段号 [A]–[F] 走 SCOPED_TRACE
graphics_add_gtest(dual_complex)

# 指数映射 (expm: exp 非单射 / SU(2) 双覆盖 / BCH 测阶 三刀 + 超强判据)
# 框架: gtest (2026-10-30 迁移) —— 判据编号 [A1]–[E4] 保留, 走 SCOPED_TRACE
graphics_add_gtest(expm)

# so(3) 三维向量算子 (Cross / Hat / Vee) —— A4 李群线载体, 只依赖 core/matrix.hpp
# 框架: gtest (2026-10-30 迁移) —— double/float 双域由 TYPED_TEST_SUITE 保证
graphics_add_gtest(so3)

# 显式闸门: 同一份 linalg 源码 x {double, Complex} 的**全入口实例化**烟测
# (不是数值验收 —— 数值判据在各算法用例; 这里判"每个公开入口都真的进过编译器")
# 框架: gtest (2026-10-30 迁移) —— [A] 的取函数指针手法原样保留, [B] 走 TYPED_TEST
graphics_add_gtest(smoke)
