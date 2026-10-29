# linalg 层验收 —— 每个用例对应一个交付判据
# 公共测试工具在 demo/common.hpp
# 框架: gtest (2026-10-30 迁移) —— 判据名保留原样, 走 SCOPED_TRACE
graphics_add_gtest(solve)

add_executable(inverse inverse.cpp)
target_include_directories(inverse PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME inverse COMMAND inverse)

add_executable(cholesky cholesky.cpp)
target_include_directories(cholesky PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME cholesky COMMAND cholesky)

add_executable(eigen eigen.cpp)
target_include_directories(eigen PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME eigen COMMAND eigen)

add_executable(svd svd.cpp)
target_include_directories(svd PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME svd COMMAND svd)

add_executable(qr qr.cpp)
target_include_directories(qr PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME qr COMMAND qr)

add_executable(dual_complex dual_complex.cpp)
target_include_directories(dual_complex PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME dual_complex COMMAND dual_complex)

# 指数映射 (expm: exp 非单射 / SU(2) 双覆盖 / BCH 测阶 三刀 + 超强判据)
# 框架: gtest (2026-10-30 迁移) —— 判据编号 [A1]–[E4] 保留, 走 SCOPED_TRACE
graphics_add_gtest(expm)

# so(3) 三维向量算子 (Cross / Hat / Vee) —— A4 李群线载体, 只依赖 core/matrix.hpp
add_executable(so3 so3.cpp)
target_include_directories(so3 PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME so3 COMMAND so3)

# 显式闸门: 同一份 linalg 源码 x {double, Complex} 的**全入口实例化**烟测
# (不是数值验收 —— 数值判据在各算法用例; 这里判"每个公开入口都真的进过编译器")
add_executable(smoke smoke.cpp)
target_include_directories(smoke PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME smoke COMMAND smoke)
