#pragma once

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include <vector>

template <typename T>
struct LUResult
{
    Matrix<T> L;
    Matrix<T> U;
    std::vector<int> piv;
};

template <CommutativeScalar T>
void SwapRow(Matrix<T> &A, int row)
{
    // 交换矩阵A的row行和max_row行
    for (int j = 0; j < A.Cols(); j++) {
        std::swap(A(row, j), A(max_row, j));
    }
}

// piv记录换行, 本质上就是在说, 怎么重排A使得PA可以无痛分解成LU, 当然B也得重排
// eg 1 / 1e-20 = 1e20数据爆炸, 没办法用, 必须用小数
template <CommutativeScalar T>
LUResult<T> LU(const Matrix<T>& A)
{
    LUResult<T> result;
    result.L = Matrix<T>::Identity(A.Rows());
    result.U = A;
    result.piv.resize(A.Cols());
    for (int j = 0; j < A.Cols(); j++) {
        int max_row = j;
        T max_val = result.U(j, j);
        for (int i = j; i < A.Rows(); i++) {
            if (result.U(i, j) > max_val) {
                max_row = i;
                max_val = result.U(i, j);
            }
        }
        if (ScalarTraits<T>::abs2(max_val) <= Tol) {
            throw std::runtime_error("matrix is singular");
            return result;
        }

        if (max_row != j) {
            SwapRow(result.U, j, max_row);
        }
        for (int i = j+1; i < A.Rows(); i++) {
            if (max_row != i) {
                result.L(i, j) = result.U(i, j) / max_val;
                // 记录m构造下三角L
                // 更新U
                for (int k = j; k < A.Cols(); k++) {
                    result.U(i, k) -= result.L(i, j) * result.U(i, k);
                }
            }
        }

        piv[j] = max_row;
    }
    return result;
}

template <CommutativeScalar T>
Matrix<T> Solve(const Matrix<T>& A, Matrix<T> B)
{
    LUResult<T> result = LU(A);
    // TODO: 实现LU分解后的求解逻辑
    // 1. 把piv应用到B上
    for (int i = 0; i < B.Rows(); i++) {
        if (result.piv[i] != i) {
            SwapRow(B, i, result.piv[i]);
        }
    }
    // 2. 把L-1应用到B上, LUX = PB
    for (int i = 0; i < B.Rows(); i++) {
        for (int k = 0; k < i; k++) {
            for (int j = 0; j < B.Cols(); j++) {
                B(i, j) -= result.L(i, k) * B(k, j);
            }
        }
    }
    // 3. 解UX = B
    Matrix<T> X(A.Rows(), B.Cols());
    for (int i = A.Rows() - 1; i >= 0; i--) {
        for (int j = 0; j < B.Cols(); j++) {
            for (int k = i+1; k < A.Cols(); k++) {
                B(i, j) -= result.U(i, k) * X(k, j);
            }
        }
        for (int j = 0; j < B.Cols(); j++) {
            X(i,j) = B(i, j) / result.U(i, i);
        }
    }
    return X;
}