# core 层验收
#
# 框架标注 (迁移期):
#   gtest  : scalar_traits  (2026-10-30 迁移)
#   legacy : tensor         (骨架期, 24 条 [P] 占位判据断言"仍 throw";
#                            迁移与实现合并做, 避免把占位判据写两遍)

# 标量层契约验收 (Scalar / CommutativeScalar / fromReal / conj / abs2)
graphics_add_gtest(scalar_traits)

# Tensor / TensorView (契约真源: acceptance.md 红线 #4)
# 实现序 ① 几何工具 → ② Tensor → ③ 视图 → ④ 深拷贝 → ⑤ 逐元素
# (与 src/core/tensor.hpp 头注释的实现顺序一致)
graphics_add_legacy_test(tensor)
