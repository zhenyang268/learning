# quat 层验收 (plan.md 阶段 2 第 5 条; V2.5 / V2.6)
# 现在是**骨架**: [A]/[B] 是真判据, [P] 是「断言仍 throw」的占位判据
# (实现掉一条, 对应 [P] 就变红 —— 红就是提醒把占位换成真判据)
add_executable(quaternion quaternion.cpp)
target_include_directories(quaternion PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME quaternion COMMAND quaternion)
