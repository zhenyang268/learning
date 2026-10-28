#pragma once

// ===========================================================================
// linalg/solve.hpp —— LU 分解与线性方程组求解
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <typename T> struct LUResult { Matrix<T> L, U; std::vector<int> piv; };
//
//   template <CommutativeScalar T> LUResult<T> LU(const Matrix<T>& A);
//   template <CommutativeScalar T> Matrix<T>   Solve(const Matrix<T>& A, Matrix<T> B);
//   template <CommutativeScalar T> void        SwapRow(Matrix<T>& A, int r1, int r2);
//   template <CommutativeScalar T> void        SwapMultipliers(Matrix<T>& L, int r1, int r2, int upto);
//
// 形状: A 必须 n×n 方阵。
//   · 非方阵                    -> throw std::invalid_argument
//   · 主元 |p| <= DefaultTol<T>(n) -> throw std::runtime_error (消息带容差与位置)
//   · Solve: A.Rows() != B.Rows() -> throw std::invalid_argument
//
// 返回: P·A = L·U, 其中 P 由 `piv` 记录 (piv[i] = 第 i 步与谁交换)。
//       L 为单位下三角, U 为上三角。
//       Solve(A,B) 返回 X 使 A·X = B, **X 与 B 同形状** (B 可以是矩阵 = 多右端项)。
//
// 关键口径: 全库**列向量** (D19)。所以 Solve(A, b) 的 b 是 n×1 列向量,
//   矩阵形式 A·X = B 是自然序, **不是** X·A = B。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   LU: 逐列消元, **每列都选主元** (列主元 = 部分选主元)。
//       piv[j] 记录第 j 步换了哪一行; L 存乘子, U 存消元结果。
//   Solve: 先用同一套消元把 B 一起做行变换, 再回代。不显式求 A^{-1}
//          —— 求逆的 flops 是解方程的数倍, 且精度更差。
//   SwapMultipliers: 换 U 的行时, L 里**已经算好的乘子**必须跟着搬。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么必须选主元 —— 不选会"数据爆炸"】
//   反例: 主元 1e-20, 下一行除以它得 1e20。这不是精度问题, 是**完全不可用**。
//   选主元 (把当前列绝对值最大的行换上来) 让乘子 |l_ij| <= 1, 误差可控。
//   判据: 不选主元的 LU 在 [[ε,1],[1,1]] 这类矩阵上会给出垃圾 —— 手算就能看出。
//
// 【为什么 piv 记录换行就够 —— P 不必显式存】
//   "怎么重排 A 使得 PA 可以无痛分解成 LU" 这句话就是 piv 的全部语义。
//   需要 P 时按 piv 从单位阵重放一遍即可; 存一个大矩阵反而多一份不同步的风险。
//
// 【乘子必须跟着换行 —— 修过的 bug】
//   L(i,k) 的语义是「**当时坐在位置 i 的那一行**的乘子」。行换了位置,
//   这段历史也得跟着走, 否则 L·U 复现的不是 PA 而是某个错位的东西。
//   ⚠ 只换 [0, upto) 这一段: 更右的乘子还没算出来; 而单位下三角的对角 1
//     与行号绑定, 整行换会把 1 搬到上三角 —— 结构当场就坏。
//   这类错**不会**让 L·U 变成明显的垃圾, 而是给出"数字看着合理但解错"的结果。
//
// 【容差为什么不在这里定义】
//   DefaultTol<T>(n) 在 core/scalar.hpp —— solve / inverse / cholesky / eigen /
//   svd 共用同一份。放在这里会让其它文件为了取容差而 include solve.hpp,
//   形成不必要的依赖链。见 numeric-contract.md §4。
// ===========================================================================

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"

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

// piv 记录换行 (为什么 P 不必显式存, 见文件头「三」)
// 反例: 1 / 1e-20 = 1e20, 数据爆炸 —— 这就是必须选主元的原因。
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
            // 算法失败 -> throw, 消息带上下文 (容差 + 位置)
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