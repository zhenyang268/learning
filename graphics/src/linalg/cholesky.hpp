#pragma once

#include <stdexcept>
#include <string>
#include "core/matrix.hpp"
#include "core/scalar.hpp"

// ===========================================================================
// linalg/cholesky.hpp —— Cholesky 分解 A = L·L^H
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <CommutativeScalar T>
//   Matrix<T> cholesky(const Matrix<T>& A)                      // 返回 L
//   template <CommutativeScalar T>
//   Matrix<T> cholesky_solve(const Matrix<T>& A, Matrix<T> B)   // 返回 x
//
// 形状: A 必须 n×n 方阵且 **Hermitian 正定** (n >= 1)。
//   · 非方阵 / 空矩阵                        -> throw std::invalid_argument
//   · 非正定 (主元实部 <= DefaultTol<T>(n))  -> throw std::runtime_error
//   · cholesky_solve: B.Rows() != n          -> throw std::invalid_argument
// 返回: L 为 n×n **下三角**, 对角元取正号 (唯一)。A = L·L^H。
//
// ⚠ 异常类型口径不统一 (历史遗留): 「非正定」抛的是 `std::runtime_error`,
//   而 `similar` 的其它入口 (expm / FromAxisAngle) 用的是 `domain_error`。
//   两者语义上都算"输入不在定义域", **是 runtime_error 用错了**。
//   改它会动到已有用例的 catch 分支, 故此处**如实记录不改**, 新代码请用
//   `domain_error` (见 numeric-contract.md §5 三分)。
//
// 前置条件: Hermitian **由调用方保证** —— 本函数只读下三角, 不检查
//   A(i,j) == conj(A(j,i))。传一个非 Hermitian 矩阵进来不会报错,
//   只会得到"下三角部分定义的"那个矩阵的分解。
//
// 不变量:
//   · L(j,j) > 0 对所有 j (正定性的直接结果)。
//   · L 唯一 —— 因为对角元被钉成正号。这一点让对拍非常干净。
//
// 域: CommutativeScalar (Dual 可用 —— 见末尾 conduit 说明)。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   · 自上而下按列推进, 不选主元。**不选主元是正确性的关键**: 正定性保证
//     每个主元严格为正, 于是不存在"需要换行"的情形; 一旦真的需要换,
//     说明矩阵不是正定, 这时候该报错, 而不是救它。
//   · 容差与 solve.hpp 同一口径: 主元实部 <= DefaultTol<T>(n) 判非正定。
//   · 复数域一律写 `L(j,k) * conj(L(j,k))` 而**不用** abs2 —— 原因见末尾。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【元素公式】
//     d      = A(j,j) - sum_{k<j} L(j,k)·conj(L(j,k))          ← 实数
//     L(j,j) = sqrt(d)
//     L(i,j) = ( A(i,j) - sum_{k<j} L(i,k)·conj(L(j,k)) ) / L(j,j)      (i > j)
//
// 【为什么它比 LU 便宜】
//   1. 不选主元 —— 没有换行, 也就没有"乘子要跟着换行搬"那类簿记
//      (LU 侧修过的 bug 在这里根本不存在)。
//   2. flops 约一半: n³/3 vs LU 的 2n³/3。
//   3. L 唯一 (对角元取正号), 对拍干净。
//
// 【存在性定理: A 正定 ⟺ Cholesky 走通】
//   所以本函数同时是 **SPD 的判据** —— inverse.hpp 的 is_positive_definite
//   直接调它, 不另写一个特征值判据。这是"用算法本身当判据"的少数正当场合:
//   因为等价性是双向的 (不是"能用它近似判断")。
//
// 【为什么写 L(j,k)*conj(L(j,k)) 而不是 abs2 —— 踩点】
//   `ScalarTraits<Dual<T>>::abs2` 只看 value 槽, **会丢掉 epsilon 分量**;
//   而 `conj` 对 Dual 是逐分量的, 于是 Dual 矩阵也能正确传播导数。
//   ⇒ 凡是**要传播导数**的地方, 用 x*conj(x) 而不是 abs2;
//     凡是**纯数值判据** (比主元、判收敛), 才用 abs2。
//   这两者容易混 —— 判据用 abs2、传播用 conj, 分工见 numeric-contract.md D7。
// ===========================================================================

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
