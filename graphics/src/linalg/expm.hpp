// linalg/expm.hpp
#pragma once
#include <cmath>
#include <stdexcept>
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "linalg/operator.hpp"   // Abs2: ‖A‖_F²

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
