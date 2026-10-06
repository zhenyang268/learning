#pragma once

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"

// 对矩阵所有元素应用同一函数
template <typename T>
Matrix<T> Apply(Matrix<T> &A, FunctionType f)
{
    Matrix<T> B(A.rows, A.cols);
    for (int i = 0; i < A.rows * A.cols; i++) {
        B.data[i] = f(A.data[i]);
    }

    return B;
}

// 逐元素应用函数矩阵
template <typename T>
Matrix<T> Apply(Matrix<T> &A, Matrix<DualFunc<T>> &matrixf)
{
    if (matrixf.rows == 1 && matrixf.cols == 1) {
        return Apply(A, matrixf[0][0]); // 单函数退化为统一应用
    }

    // 应该只有一维, 但是后续染色节点打包可能涉及多行
    if (A.cols != matrixf.cols) {
        throw std::invalid_argument("Matrix apply dimension mismatch!");
    }

    Matrix<T> B(A.rows, A.cols);
    for (int i = 0; i < A.rows * A.cols; i++) {
        B.data[i] = dualFunc(A.data[i], matrixf.data[i]);
    }

    return B;
}

// 把value赋值为derive用于backward
template <typename T>
Matrix<T>& Derive(Matrix<T> &A)
{
    for (int i = 0; i < A.rows * A.cols; i++) {
        A.data[i].v = A.data[i].d;
    }

    return A;
}

// 累加: 
template <typename T>
Matrix<T> SumAll(Matrix<T> &A)
{
    Matrix<T> B = A;
    sumLine(B[0], B.Size());
    return B;
}

// 累乘
template <typename T>
Matrix<T> MultAll(Matrix<T> &A)
{
    Matrix<T> B = A;
    multLine(B.GetData(), B.Size());
    return B;
}