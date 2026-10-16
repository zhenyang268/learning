// linalg/expm.hpp
#pragma once
#include <cmath>
#include <stdexcept>
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "linalg/operator.hpp"   // Abs2: ‖A‖_F²
#include "linalg/solve.hpp"     // Padé 版的 D^{-1}·N (解线性系统, 不显式求逆)

// expm(A) = exp(A/2^s)^(2^s) —— Taylor 级数 + scaling-and-squaring.
// 直接对 A 展开会撞精度地板 (中间项 ~e^‖A‖/√(2π‖A‖), 与最终和灾难性相消);
// 先折半到 ‖As‖_F <= 0.5, 级数变成递减型, 精度回到 ~K·eps.
// Padé 版留给阶段 2 / V2.4, 本版是它的对照基线.
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
// Padé 版 (阶段 2 / V2.4) —— [13/13] 有理逼近 + scaling-and-squaring
//
//   exp(A) ≈ D^{-1}·N,  N = Σ_{k=0}^{m} c_k A^k,  D = Σ_{k=0}^{m} c_k (−A)^k
//   c_k = (2m−k)!·m! / ((2m)!·k!·(m−k)!),  m = 13 (Higham 的标准阶)
//
// 与上面 Taylor 版的关系: **互为对照基线** (V2.4 判据 = 两条路之差 < 1e-8)。
// 为什么值得多这一条: 级数版靠"折半到 ‖As‖ ≤ 0.5"换精度; Padé 靠有理逼近,
// 同样项数下 ‖As‖ 的容许上限大得多 (θ₁₃ ≈ 5.37) ⇒ **s 更小、平方更少**。
// 两条路走的数学完全不同, 所以对拍才有信息量 (matrix.md §10.9: 判据独立于实现路径)。
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
