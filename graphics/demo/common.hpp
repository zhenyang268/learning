#pragma once
// ---------------------------------------------------------------------------
// demo/common.hpp -- C++ 测试侧公共工具 (plan.md 阶段 1 要求)
//
// 只放"怎么算对"的东西, 不放被测实现:
//   MaxAbs2 / NormInf   元素级量度 (比较一律走 abs2, D7; 需要真实 ∞ 范数时才开方)
//   IsApprox            ‖A - B‖∞ <= tol
//   RandMatrix<T>       实 / 复 / Dual 随机矩阵 (定 seed, 可复现)
//   Checker             [OK]/[FAIL] 记分板 + 末行汇总 (run_tests.sh -d 依赖这两个标记)
//
// 纪律 (matrix.md §10.9 第 2 条): 判据先行 —— 用例先写"多少算对", 代码对着判据写。
// ---------------------------------------------------------------------------

#include "core/common.hpp"
#include "core/complex.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "autodiff/dual.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <type_traits>

// ---- 元素级量度 -----------------------------------------------------------

// max |a_ij|^2: 比较用这个 (免开方, 阈值取 tol^2 即等价于 ‖·‖∞ <= tol)
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t MaxAbs2(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;
    // real_t 是具体实数类型(float/double), 这里的 0 是实数域的 0, 不是"T 域里的字面量"
    // (不能用 ScalarTraits<T>::zero(): 对 Complex/Dual 它返回的是 T, 不是 real_t)
    real_t best = real_t(0);
    for (int i = 0; i < A.Size(); i++) {
        const real_t v = ScalarTraits<T>::abs2(A.Data()[i]);
        if (v > best) best = v;
    }
    return best;
}

// 真正的 ∞ 范数 (只用在打印/打分, 不进算法分支 —— D7)
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t NormInf(const Matrix<T>& A)
{
    return std::sqrt(MaxAbs2(A));
}

template <CommutativeScalar T>
bool IsApprox(const Matrix<T>& A, const Matrix<T>& B,
              typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != B.Rows() || A.Cols() != B.Cols()) return false;
    return MaxAbs2(A - B) <= tol * tol;
}

// ---- 随机矩阵: 实 / 复 / Dual (V1.5 的 Dual<Complex> 也在内) ----------------

template <typename T> struct RandFill;

template <> struct RandFill<float> {
    static float Make(double a, double) { return static_cast<float>(a); }
};
template <> struct RandFill<double> {
    static double Make(double a, double) { return a; }
};
template <typename U> struct RandFill<Complex<U>> {
    static Complex<U> Make(double a, double b)
    {
        return Complex<U>(static_cast<U>(a), static_cast<U>(b));
    }
};
template <typename U> struct RandFill<Dual<U>> {
    static Dual<U> Make(double a, double b)
    {
        return Dual<U>(static_cast<U>(a), static_cast<U>(b));
    }
};

template <CommutativeScalar T>
Matrix<T> RandMatrix(int rows, int cols, unsigned seed = 1u,
                     double lo = -1.0, double hi = 1.0)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(lo, hi);
    Matrix<T> A(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // 复数填实部+虚部; 实数取第一个; Dual 填 value+derive (derive 小, 近似当作单位方向)
            A(i, j) = RandFill<T>::Make(dist(rng), dist(rng) * 0.125);
        }
    }
    return A;
}

// ---- 记分板 ---------------------------------------------------------------

struct Checker
{
    int ok = 0;
    int fail = 0;

    void check(bool pass, const std::string& name)
    {
        if (pass) {
            ok++;
            std::printf("  [OK]   %s\n", name.c_str());
        } else {
            fail++;
            std::printf("  [FAIL] %s\n", name.c_str());
        }
    }
};

inline int CheckReport(int ok, int fail, const char* title)
{
    std::printf("\n%s: 检查 %d 项 (FAIL %d)\n", title, ok + fail, fail);
    if (fail == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("ALL TESTS FAILED (%d failures)\n", fail);
    return 1;
}
