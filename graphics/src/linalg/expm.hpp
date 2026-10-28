#pragma once

// ===========================================================================
// linalg/expm.hpp —— 矩阵指数 exp(A)
// ===========================================================================
//
// 【状态】已实现两条路线: Taylor + scaling-and-squaring (expm) 与
//   Padé [13/13] + scaling-and-squaring (expm_pade)。互为对照基线。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <CommutativeScalar T> Matrix<T> expm(const Matrix<T>& A);
//   template <CommutativeScalar T> Matrix<T> expm_pade(const Matrix<T>& A);
//
// 形状: A 必须 n×n 方阵, n >= 1。
//   · 非方阵 / 空矩阵 -> throw std::invalid_argument
//
// 返回: exp(A), n×n。对任意 T (double / Complex / Dual) 成立。
//
// 精度边界 (判据设计必须知道):
//   · 误差是**相对**于元素量级 e^‖A‖ 的, 不是绝对的。判据必须写成
//     ‖X - X_ref‖ / ‖X_ref‖, 用绝对差会在 ‖A‖ 大时得出假 FAIL。
//   · 实测 (double, 见 demo/linalg/expm.cpp): 3×3 相对差 5.3e-16,
//     5×5 2.5e-15, ‖A‖≤100 的 5×5 1.5e-13, 复 4×4 1.1e-15。
//   · Padé 与 Taylor 之差是**最强判据** (两条独立路径)。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   Taylor 版: expm(A) = exp(A/2^s)^(2^s)
//     ① 求 s 使 ‖A/2^s‖_F <= 0.5
//     ② 对折半后的矩阵展开 Taylor 级数 (此时是递减型, 不会灾难性相消)
//     ③ 平方 s 次还原
//   Padé 版: 同 scaling, 但中段用 [13/13] 有理逼近
//     N = Σ c_k P^k, D = Σ (-1)^k c_k P^k, X = D^{-1}·N
//     **解线性系统而不是显式求逆** (沿用 solve.hpp 的纪律)
//   两条路共用同一个 scaling 与 squaring 骨架, 只有中段不同。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么必须 scaling-and-squaring —— 精度地板】
//   直接对 A 展开 Taylor: 中间项量级 ~ e^‖A‖/√(2π‖A‖), 而最终和可能远小于
//   最大中间项 (元素有正有负时灾难性相消)。于是精度被卡在
//       ~ (e^‖A‖/√(2π‖A‖)) · eps
//   这个地板与算法无关, 只与"直接在原尺度展开"有关。
//   先折半到 ‖As‖_F <= 0.5, 级数变成递减型, 精度回到 ~K·eps。
//   ⇒ 任何通用 exp 实现都不能省这一步。
//
// 【为什么 exp(−A)·exp(A) 不能当唯一判据】
//   条件数 ≈ e^{2‖A‖}。‖A‖ = 100 时那是 e^200 ≈ 7e86 —— 任何实现的残渣都会
//   被放大到这个量级。所以 `exp(A)exp(−A) ≈ I` 在 ‖A‖ 大时**必然**看起来很差,
//   这是**判据病态, 不是算法错**。
//   超强独立判据用 det(exp X) = exp(tr X) —— 两个标量, 都良态。
//
// 【为什么 Padé 解线性系统而不求逆】
//   求逆的 flops 更大 (2n³ vs 解一座的 n³/3), 而且精度更差。
//   判据: 只要不需要显式 D^{-1}, 就解方程。这条是全库纪律。
//
// 【E_pade_vs_taylor 那段曾被插两遍】
//   2026-10-17 一次「读—改—写」脚本不幂等, 同一段被插两次, 触发 redefinition
//   构建失败。教训: 改脚本必须 assert 目标串不存在 (见 docs/pitfalls.md)。
// ===========================================================================

#include <cmath>
#include <stdexcept>

#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "linalg/operator.hpp"   // Abs2: ‖A‖_F²
#include "linalg/solve.hpp"      // Padé 版的 D^{-1}·N (解线性系统, 不显式求逆)

template <CommutativeScalar T>
Matrix<T> expm(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    const int n = A.Rows();
    if (A.Rows() != A.Cols()) {
        throw std::invalid_argument("expm: A must be square");
    }
    if (n == 0) {
        throw std::invalid_argument("expm: A must be non-empty");
    }

    // ① 定 s: 最小的 s 使 ‖A‖_F / 2^s <= 0.5。用平方量比, 省一次开方 (0.25 = 0.5²)。
    //    范数齐次 ⇒ 不必真的把矩阵除一遍再量, 只在实数上反复 /4。
    real_t n2 = Abs2(A);                  // ‖A‖_F²
    int s = 0;
    while (n2 > real_t(0.25)) {
        n2 = n2 / real_t(4);              // (‖A‖/2)² = ‖A‖²/4
        s++;
    }

    // ② 缩小: As = A / 2^s (ldexp 取 2^-s; 不用 1/(1<<s), s 大时移位溢出 UB)
    const Matrix<T> As = A * ScalarTraits<T>::fromReal(std::ldexp(real_t(1), -s));

    // ③ Taylor 级数。相对判据提前退出 (项数 12~15 自然停), k<=30 只是防 NaN/Inf 的保险丝。
    //    标量除法必须走 fromReal —— Complex 无单参构造, term * 0.5 编译不过 (D7/D15)。
    const real_t tol = DefaultTol<T>(n);
    Matrix<T> result = Matrix<T>::Identity(n);
    Matrix<T> term   = Matrix<T>::Identity(n);
    for (int k = 1; k <= 30; k++) {
        term = term * As * ScalarTraits<T>::fromReal(real_t(1) / real_t(k));   // 矩阵乘 + 标量乘
        result += term;
        if (Abs2(term) <= tol * tol * Abs2(result)) break;   // 一切比较走 Abs2 (D7)
    }

    // ④ 缩回: 平方 s 次
    for (int i = 0; i < s; i++) {
        result = result * result;         // 矩阵乘, 不是 cwiseProduct
    }
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Padé 版 —— [13/13] 有理逼近 + scaling-and-squaring
//
//   exp(A) ≈ D^{-1}·N,  N = Σ_{k=0}^{m} c_k A^k,  D = Σ_{k=0}^{m} c_k (−A)^k
//   c_k = (2m−k)!·m! / ((2m)!·k!·(m−k)!),  m = 13 (Higham 的标准阶)
//
// 与上面 Taylor 版的关系: **互为对照基线** (判据 = 两条路之差, **相对** < 1e-8)。
// 为什么值得多这一条: 级数版靠"折半到 ‖As‖ ≤ 0.5"换精度; Padé 靠有理逼近,
// 同样项数下 ‖As‖ 的容许上限大得多 (θ₁₃ ≈ 5.37) ⇒ **s 更小、平方更少**。
// 两条路走的数学完全不同, 所以对拍才有信息量 (判据独立于实现路径)。
//
// 未做的优化: 标准实现会把 N/D 按 A²、A⁴、A⁶ 分块做 Horner, 把矩阵乘从 ~13 次压到 ~6 次;
// 本版直接求和 (清晰优先), 效率只由 θ₁₃ 带来的 s 缩减贡献一部分。
// ─────────────────────────────────────────────────────────────────────────────
template <CommutativeScalar T>
Matrix<T> expm_pade(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    const int n = A.Rows();
    if (A.Rows() != A.Cols()) {
        throw std::invalid_argument("expm_pade: A must be square");
    }
    if (n == 0) {
        throw std::invalid_argument("expm_pade: A must be non-empty");
    }

    constexpr int m = 13;
    const real_t theta13 = real_t(5.371920351148152);   // [13/13] 的平衡常数 (Higham 2009)

    // ① 定 s: 最小的 s 使 ‖A‖_F / 2^s <= θ₁₃。
    //    注意与 Taylor 版的差别 —— 那边阈值是 0.5, 这里 5.37, 所以 s 通常小 3 个数。
    real_t w = std::sqrt(Abs2(A));
    int s = 0;
    while (w > theta13) {
        w = w / real_t(2);
        s++;
    }

    // ② 缩小
    const Matrix<T> As = A * ScalarTraits<T>::fromReal(std::ldexp(real_t(1), -s));

    // ③ 系数 c_k: c_0 = 1, c_k = c_{k-1}·(m−k+1) / (k·(2m−k+1))
    //    (k=1 时给 13/26 = 1/2, 与 [13/13] 的标准系数表一致)
    real_t c[m + 1];
    c[0] = real_t(1);
    for (int k = 1; k <= m; k++) {
        c[k] = c[k - 1] * real_t(m - k + 1) / (real_t(k) * real_t(2 * m - k + 1));
    }

    // ④ N = Σ c_k As^k, D = Σ c_k (−As)^k  (对称: D 就是把奇次项取负)
    Matrix<T> N = Matrix<T>::Identity(n);
    Matrix<T> D = Matrix<T>::Identity(n);
    Matrix<T> P = Matrix<T>::Identity(n);            // P 依次是 As^1, As^2, ...
    for (int k = 1; k <= m; k++) {
        P = P * As;
        N += P * ScalarTraits<T>::fromReal(c[k]);
        D += P * ScalarTraits<T>::fromReal((k % 2 == 0) ? c[k] : -c[k]);
    }

    // ⑤ X = D^{-1}·N —— 解线性系统, 不显式求逆 (沿用 Solve 的既有纪律)
    Matrix<T> X = Solve(D, N);

    // ⑥ 缩回
    for (int i = 0; i < s; i++) {
        X = X * X;
    }
    return X;
}
