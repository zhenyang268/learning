# quat 层验收 —— 判据真源: graphics/acceptance.md
# 已全部落地: 9 个入口无占位, 160 项判据; 剩 5 处 throw 全是契约性
#   (2 处 domain_error + 3 处 invalid_argument)
#
# 框架: gtest + TYPED_TEST (2026-10-30 迁移) —— double / float 双域由
# TYPED_TEST_SUITE 保证, 不再手写两遍调用。
graphics_add_gtest(quaternion)
