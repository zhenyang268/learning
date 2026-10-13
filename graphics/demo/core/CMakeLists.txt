# 标量层契约验收 (Scalar / CommutativeScalar / fromReal / conj / abs2)
add_executable(scalar_traits scalar_traits.cpp)

target_include_directories(scalar_traits PRIVATE
    ${GRAPHICS_SRC_DIR}
)

add_test(NAME scalar_traits COMMAND scalar_traits)
