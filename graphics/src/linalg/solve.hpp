#pragma once

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

// 运算符纪律: 只用 Scalar concept 保证的 + - * / 。
// Complex / Dual 都没有复合赋值 (也没有 operator>), 泛型 linalg 代码里写 `x -= y`
// 会在 Matrix<Complex<double>> 上直接编不过 —— 这不是 Complex 的 bug, 是 concept 的口径。
//
// LU 分解 (带部分选主元): P·A = L·U
//   piv[j] = 第 j 步选出的主元行号 (在"当前部分矩阵"里的行号);
//   按 j 递增顺序依次交换 (j, piv[j]), 就是置换 P 本身。
// 换行时乘子列必须跟着换 —— L 的行序 = 置换后 A 的行序, 见 SwapMultipliers。
template <typename T>
struct LUResult
{
    Matrix<T> L;
    Matrix<T> U;
    std::vector<int> piv;
};

// D12 / matrix.md §2.4: 容差 = eps 的倍数, 且随规模缩放 —— 不写魔数。
// 不能用全局 Tol: 那是 float 的 epsilon, 对 Matrix<double> 量级完全不对。
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t DefaultTol(int n)
{
    using real_t = typename ScalarTraits<T>::real_t;
    return real_t(8) * real_t(n) * std::numeric_limits<real_t>::epsilon();
}

// 整行交换 (r1, r2 都是在当前矩阵里的行号)
template <CommutativeScalar T>
void SwapRow(Matrix<T>& A, int r1, int r2)
{
    for (int j = 0; j < A.Cols(); j++) {
        std::swap(A(r1, j), A(r2, j));
    }
}

// L 的乘子列 [0, upto) 互换。
// 不能整行换: L 是单位下三角, 对角线上的 1 与行序号绑定, 整行换会把结构换坏
// (换完 P·A != L·U, 而结果看着"像那么回事")。
template <CommutativeScalar T>
void SwapMultipliers(Matrix<T>& L, int r1, int r2, int upto)
{
    for (int k = 0; k < upto; k++) {
        std::swap(L(r1, k), L(r2, k));
    }
}

// eg 1 / 1e-20 = 1e20 数据爆炸 -> 必须选主元, 否则小主元把小误差放大成大错
template <CommutativeScalar T>
LUResult<T> LU(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    if (A.Rows() != A.Cols())
        throw std::invalid_argument("LU: A 必须是方阵");

    const int n = A.Rows();
    const real_t tol = DefaultTol<T>(n);
    const real_t tol2 = tol * tol;   // 下面比 abs2 <= tol^2, 等价于 |p| <= tol 且免开方 (D7)

    LUResult<T> result;
    result.L = Matrix<T>::Identity(n);
    result.U = A;
    result.piv.assign(n, 0);

    for (int j = 0; j < n; j++) {
        // --- 选主元: 比较只走 abs2 —— Complex / Dual 没有 operator> (D7) ---
        int max_row = j;
        real_t best = ScalarTraits<T>::abs2(result.U(j, j));
        for (int i = j + 1; i < n; i++) {
            const real_t v = ScalarTraits<T>::abs2(result.U(i, j));
            if (v > best) {
                best = v;
                max_row = i;
            }
        }
        if (best <= tol2) {
            throw std::runtime_error("LU: 主元 |p| <= tol, 矩阵奇异或近奇异");
        }
        result.piv[j] = max_row;

        if (max_row != j) {
            SwapRow(result.U, j, max_row);
            SwapMultipliers(result.L, j, max_row, j);
        }

        // --- 消元: i > j 的每一行都要做 ---
        // 原先写成 `if (max_row != i) { ... }` 是错的: 主元行恰好落在 i 上时整段被跳过,
        // U 少消一行、L 少一个乘子, 编译能过而结果静默出错。
        const T pivot = result.U(j, j);
        for (int i = j + 1; i < n; i++) {
            const T mul = result.U(i, j) / pivot;
            result.L(i, j) = mul;
            for (int k = j; k < n; k++) {
                result.U(i, k) = result.U(i, k) - mul * result.U(j, k);
            }
        }
    }
    return result;
}

// 解 A·X = B (B 单列时即教科书的 A·x = b, D19 修订版: 列向量)
//   1) 换行作用到 B 上: P·A·X = P·B
//   2) 前代解 L·Y = P·B
//   3) 回代解 U·X = Y
template <CommutativeScalar T>
Matrix<T> Solve(const Matrix<T>& A, Matrix<T> B)
{
    if (A.Rows() != A.Cols())
        throw std::invalid_argument("Solve: A 必须是方阵");
    if (B.Rows() != A.Rows())
        throw std::invalid_argument("Solve: B 的行数必须等于 A 的阶数");

    const int n = A.Rows();
    const int nCols = B.Cols();

    const LUResult<T> lu = LU(A);

    // 1) 置换: 顺序与记录顺序一致 (j 递增), 依次交换 (j, piv[j])
    for (int j = 0; j < n; j++) {
        if (lu.piv[j] != j) {
            SwapRow(B, j, lu.piv[j]);
        }
    }

    // 2) 前代: L 单位下三角, 乘子已在 LU 里跟着换过行
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < i; k++) {
            const T mul = lu.L(i, k);
            for (int j = 0; j < nCols; j++) {
                B(i, j) = B(i, j) - mul * B(k, j);
            }
        }
    }

    // 3) 回代 (B 此时就是 Y; 不原地改 Y, 回代里要读它)
    Matrix<T> X(n, nCols);
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < nCols; j++) {
            T acc = B(i, j);
            for (int k = i + 1; k < n; k++) {
                acc = acc - lu.U(i, k) * X(k, j);
            }
            X(i, j) = acc / lu.U(i, i);
        }
    }
    return X;
}
