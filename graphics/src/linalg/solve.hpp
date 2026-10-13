#pragma once

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include <limits>
#include <stdexcept>
#include <vector>

template <typename T>
struct LUResult
{
    Matrix<T> L;
    Matrix<T> U;
    std::vector<int> piv;
};

template <CommutativeScalar T>
void SwapRow(Matrix<T> &A, int r1, int r2)
{
    // 交换矩阵A的 r1 行和 r2 行
    for (int j = 0; j < A.Cols(); j++) {
        std::swap(A(r1, j), A(r2, j));
    }
}

// 换 U 的行时, L 里"已经算好的乘子"必须跟着搬 —— L(i,k) 的语义是
// "当时坐在位置 i 的那一行"的乘子, 行换了位置, 这段历史也得跟着走。
// 只换 [0, upto) 这一段: 更右的乘子还没算出来, 而单位下三角的对角 1 与行号绑定,
// 整行换会把 1 搬到上三角, 结构当场就坏。
template <CommutativeScalar T>
void SwapMultipliers(Matrix<T> &L, int r1, int r2, int upto)
{
    for (int k = 0; k < upto; k++) {
        std::swap(L(r1, k), L(r2, k));
    }
}

// DefaultTol<T>(n) 已上移到 core/scalar.hpp 标量层 (2026-10-14, S2a):
// 它是 D12 的容差策略, solve / inverse / cholesky / eigen / svd 共用。

// piv记录换行, 本质上就是在说, 怎么重排A使得PA可以无痛分解成LU, 当然B也得重排
// eg 1 / 1e-20 = 1e20数据爆炸, 没办法用, 必须用小数
template <CommutativeScalar T>
LUResult<T> LU(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols())
        throw std::invalid_argument("LU: A 必须是方阵");

    const typename ScalarTraits<T>::real_t tol = DefaultTol<T>(A.Rows());

    LUResult<T> result;
    result.L = Matrix<T>::Identity(A.Rows());
    result.U = A;
    result.piv.resize(A.Cols());
    for (int j = 0; j < A.Cols(); j++) {
        int max_row = j;
        T max_val = result.U(j, j);
        // 选主元只能比 |.|^2 (D7): Complex / Dual 没有 operator>, 泛型代码里不能写裸 >
        typename ScalarTraits<T>::real_t best = ScalarTraits<T>::abs2(max_val);
        for (int i = j; i < A.Rows(); i++) {
            const typename ScalarTraits<T>::real_t cur = ScalarTraits<T>::abs2(result.U(i, j));
            if (cur > best) {
                max_row = i;
                max_val = result.U(i, j);
                best = cur;
            }
        }
        if (ScalarTraits<T>::abs2(max_val) <= tol * tol) {
            // D12: 算法失败 -> throw 带上下文(容差), 消息类型中立
            throw std::runtime_error("LU: 主元近奇异 (判据 |p| <= tol, tol = 8*n*eps)");
            return result;
        }

        if (max_row != j) {
            SwapRow(result.U, j, max_row);
            SwapMultipliers(result.L, j, max_row, j);
        }

        // i > j 的每一行都要消: 换行后主元行已经坐在 j 上, "不消主元行"由 i > j 天然保证。
        for (int i = j+1; i < A.Rows(); i++) {
            result.L(i, j) = result.U(i, j) / max_val;
            // 记录m构造下三角L
            // 更新U
            for (int k = j; k < A.Cols(); k++) {
                result.U(i, k) = result.U(i, k) - result.L(i, j) * result.U(j, k);
            }
        }

        result.piv[j] = max_row;
    }
    return result;
}

template <CommutativeScalar T>
Matrix<T> Solve(const Matrix<T>& A, Matrix<T> B)
{
    if (B.Rows() != A.Rows())
        throw std::invalid_argument("Solve: B 的行数必须等于 A 的行数");

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
                B(i, j) = B(i, j) - result.L(i, k) * B(k, j);
            }
        }
    }
    // 3. 解UX = B
    Matrix<T> X(A.Rows(), B.Cols());
    for (int i = A.Rows() - 1; i >= 0; i--) {
        for (int j = 0; j < B.Cols(); j++) {
            for (int k = i+1; k < A.Cols(); k++) {
                B(i, j) = B(i, j) - result.U(i, k) * X(k, j);
            }
        }
        for (int j = 0; j < B.Cols(); j++) {
            X(i,j) = B(i, j) / result.U(i, i);
        }
    }
    return X;
}