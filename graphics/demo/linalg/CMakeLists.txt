# linalg 层验收 (plan.md V1.1–V1.5) —— 每个用例对应一个交付判据
# 公共测试工具在 demo/common.hpp
add_executable(solve solve.cpp)
target_include_directories(solve PRIVATE
    ${GRAPHICS_SRC_DIR}
)
add_test(NAME solve COMMAND solve)

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
