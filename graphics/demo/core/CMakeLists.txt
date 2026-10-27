# 标量层契约验收 (Scalar / CommutativeScalar / fromReal / conj / abs2)
add_executable(scalar_traits scalar_traits.cpp)

target_include_directories(scalar_traits PRIVATE
    ${GRAPHICS_SRC_DIR}
)

add_test(NAME scalar_traits COMMAND scalar_traits)

# Tensor / TensorView (matrix.md §6 第三次订正 / 决策表 #21, 2026-10-23)
# 容器只准哑; 运算一律自由函数 (plan.md 红线 #4 修订后仍如此)
add_executable(tensor tensor.cpp)
target_include_directories(tensor PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME tensor COMMAND tensor)
