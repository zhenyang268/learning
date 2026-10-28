#pragma once

#include <stdexcept>
#include <string>
#include "core/matrix.hpp"
#include "core/scalar.hpp"

/**
 * cholesky 分解 (Hermitian 正定): A = L * L^H, L 为下三角, 对角元为正实数
 *
 * 元素公式 (matrix.md 3.1; plan.md V2.2):
 *     d      = A(j,j) - sum_{k<j} L(j,k) * conj(L(j,k))            <- 实数
 *     L(j,j) = sqrt(d)
 *     L(i,j) = ( A(i,j) - sum_{k<j} L(i,k) * conj(L(j,k)) ) / L(j,j)    (i > j)
 *
 * 与 LU (solve.hpp) 的分工: LU 通用但必须选主元; 这里只吃正定矩阵, 换来
 *   1. 不选主元 —— 正定性保证每个主元严格为正, 于是没有换行, 也就没有
 *      "乘子要跟着换行"那类簿记 (S1 在 LU 里修掉的 bug 在这里不存在)
 *   2. flops 约一半 (n^3/3 vs LU 的 2n^3/3)
 *   3. L 唯一 (对角元取正号), 对拍干净
 *
 * 存在性定理: A 正定 <=> Cholesky 走通。所以本文件同时是 SPD 的判据 ——
 * inverse.hpp 的 is_positive_definite 直接调它 (matrix.md 3.2)。
 *
 * D7 : 共轭一律走 conj, 不把 x*x 当模方 (复数下必须 conj)
 * D12: 非方阵 / 空 / 非正定 -> throw, 不拿空矩阵当错误码
 *
 * @param A Hermitian 正定方阵 (Hermitian 由调用方保证; 本函数只读下三角)
 * @return L, 满足 A = L * L^H
 */

template <CommutativeScalar T>
Matrix<T> cholesky(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    if (A.Rows() != A.Cols())
        throw std::invalid_argument("cholesky: A 必须是方阵");
    const int n = A.Rows();
    if (n == 0)
        throw std::invalid_argument("cholesky: A 不能为空矩阵");

    // 容差与 solve.hpp 同一口径 (D12): 主元实部 <= tol 即判非正定/近奇异
    const real_t tol = DefaultTol<T>(n);

    Matrix<T> L(n, n);
    for (int j = 0; j < n; j++) {
        // 对角元: d = A(j,j) - sum_{k<j} |L(j,k)|^2
        // 写成 L(j,k)*conj(L(j,k)) 而不是 abs2: ScalarTraits<Dual<T>>::abs2 只看 value,
        // 会丢掉 epsilon 分量; 而 conj 对 Dual 是逐分量的, 于是 Dual 矩阵也能传播导数。
        T d = A(j, j);
        for (int k = 0; k < j; k++) {
            d = d - L(j, k) * ScalarTraits<T>::conj(L(j, k));
        }
        // 判据走 real(d): Complex / Dual 没有 operator<=, 而 Hermitian 下 d 本就是实数。
        // 写成 !(... > tol) 而不是 (... <= tol): NaN 也会落进 throw, 不留一路静默传播的口子。
        if (!(ScalarTraits<T>::real(d) > tol)) {
            throw std::runtime_error("cholesky: A 非正定或近奇异 (第 " + std::to_string(j)
                                     + " 个主元实部 <= tol, tol = 8*n*eps)");
        }
        using std::sqrt;          // ADL 两步法: 内建类型走 std, Complex/Dual 走类内 friend
        L(j, j) = sqrt(d);

        // 该列下方: L(i,j) = ( A(i,j) - sum_{k<j} L(i,k)*conj(L(j,k)) ) / L(j,j)
        for (int i = j + 1; i < n; i++) {
            T acc = A(i, j);
            for (int k = 0; k < j; k++) {
                acc = acc - L(i, k) * ScalarTraits<T>::conj(L(j, k));
            }
            L(i, j) = acc / L(j, j);
        }
    }
    return L;
}

// ---------------------------------------------------------------------------
// 解 A * X = B (A Hermitian 正定): A = L*L^H => L*(L^H*X) = B, 分两步
//   1. 前代 L * Y = B    (L 下三角, 对角不是 1, 所以每一步要除)
//   2. 回代 L^H * X = Y  (L^H 上三角)
// 语义与 LU 的 Solve 一致 (教科书列向量约定), 用来打 V2.2 的"解同一方程组"对拍。
// ---------------------------------------------------------------------------
template <CommutativeScalar T>
Matrix<T> cholesky_solve(const Matrix<T>& A, Matrix<T> B)
{
    if (B.Rows() != A.Rows())
        throw std::invalid_argument("cholesky_solve: B 的行数必须等于 A 的行数");

    const int n = A.Rows();
    const int m = B.Cols();
    const Matrix<T> L = cholesky(A);          // A 非正定在这里就 throw 了

    // 1) 前代: Y(i,:) = ( B(i,:) - sum_{k<i} L(i,k)*Y(k,:) ) / L(i,i)
    //    原地写到 B 上 —— 第 i 步读的 B(k,:) (k<i) 已经是 Y(k,:)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            T acc = B(i, j);
            for (int k = 0; k < i; k++) {
                acc = acc - L(i, k) * B(k, j);
            }
            B(i, j) = acc / L(i, i);
        }
    }

    // 2) 回代: X(i,:) = ( Y(i,:) - sum_{k>i} conj(L(k,i))*X(k,:) ) / conj(L(i,i))
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < m; j++) {
            T acc = B(i, j);
            for (int k = i + 1; k < n; k++) {
                acc = acc - ScalarTraits<T>::conj(L(k, i)) * B(k, j);
            }
            B(i, j) = acc / ScalarTraits<T>::conj(L(i, i));
        }
    }
    return B;
}
