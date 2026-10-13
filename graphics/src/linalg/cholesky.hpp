#pragma once

#include <stdexcept>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"

/**
 * cholesky 分解
 * A = L * L^H, L 为下三角矩阵
 * 
 * // TODO: 补上数学介绍
 * 1. 可以使用递归推导计算
 * 2. 可以直接使用元素公式填写
 * 3. 后续可以优化分块板和LDLT非正定版本
 * 
 * @param A 矩阵
 * @return 分解结果
 */

template <CommutativeScalar T>
Matrix<T> cholesky(const Matrix<T>& A) 
{
    Matrix<T> L(A.Rows(), A.Cols());
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < i; j++) {
            L(i, j) = A(i, j);
            for (int k = 0; k < j; k++) {
                L(i, j) -= L(i, k) * L(j, k);
            }
            L(i, j) /= L(j, j);
        }
        L(i, i) = A(i, i);
        for (int j = 0; j < i; j++) {
            L(i, i) -= L(i, j) * L(i, j);
        }
        if (L(i, i) <= 0) {
            return Matrix<T>();
        }
        L(i, i) = sqrt(L(i, i));
    }
    return L;
}

template <CommutativeScalar T>
Matrix<T> cholesky_r(const Matrix<T>& A)
{
    Matrix<T> L(A.Rows(), A.Cols());
    Matrix<T> B = A;
    if (cholesky_recursive(B, L, 0)) {
        return L;
    }
    return Matrix<T>();
}

template <CommutativeScalar T>
bool cholesky_recursive(Matrix<T>& A, Matrix<T>& L, int i)
{
    if (i >= A.Rows()) {
        return true;
    }
    if (A(i, i) <= 0) {
        return false;
    }
    L(i, i) = sqrt(A(i, i));
    // TODO: 这块还没做完, 只处理了对角线, 还得继续更新A以及非对角线矩阵
    return cholesky_recursive(A, L, i + 1);
}