#pragma once

#include "core/scalar.hpp"
#include "core/matrix.hpp"

template <CommutativeScalar T> Matrix<T> Inverse(const Matrix<T>& A);

template <CommutativeScalar T> Matrix<T>
Matrix Inverse(const Matrix<T>& A) const;

bool is_orthogonal();           // 正交
bool is_symmetric();            // 对称
bool is_positive_definite();    // 正定

Matrix lu_inverse();            // 求解梯度矩阵lu
Matrix diagnal_inverse();       // 对角矩阵
Matrix cholesky_inverse();      // 对称正定 LLT, 下三角
Matrix qr_inverse();            // QR分解, 正交*上三角