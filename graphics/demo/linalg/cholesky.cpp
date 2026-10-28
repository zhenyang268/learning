// cholesky.cpp -- linalg/cholesky.hpp 验收 (S2a / V2.2)
//
// 判据 (判据先行) —— 先写"多少算对", 再对着判据写实现:
//   ① 重构: ‖L·Lᴴ − A‖∞ < 1e-8  (double 与 Complex<double> 共用同一份 linalg 源码, V1.1)
//   ② V2.2 硬指标: SPD 矩阵上解同一方程组, Cholesky 解 与 LU(Solve) 解 之差 < 1e-8,
//      且 Cholesky 各自的残差 ‖A·X − B‖∞ < 1e-10
//   ③ 结构: L 的严格上三角全 0; 对角元为正实数 (复数的虚部可忽略)
//   ④ 解析对拍: 手算 2x2 / 3x3 / Hermitian 2x2 的 L 逐元一致
//   ⑤ 异常: 非方阵 / 0x0 / 非正定 / 半正定 / 秩亏边界 全部 throw (D12, 不拿空矩阵当错误码)
//   ⑥ SPD 判据语义: cholesky 走通 <=> A 正定 (inverse.hpp 的 is_positive_definite 判据本体)
//
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh                  # 全量
//   ./graphics/build/demo/linalg/cholesky    # 单跑

#include "../common.hpp"
#include "linalg/cholesky.hpp"
#include "linalg/solve.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

static Checker ck;

using D = double;
using C = Complex<double>;

constexpr double V22_TOL = 1e-8;      // V2.2 解之差硬阈值
constexpr double CHOL_RES_TOL = 1e-10;

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

// SPD 随机矩阵: A = M·Mᴴ + n·I —— 对任意 M 都正定 (特征值 >= n), 且实/复/Dual 通用
template <CommutativeScalar T>
static Matrix<T> RandSPD(int n, unsigned seed)
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T> M = RandMatrix<T>(n, n, seed);
    Matrix<T> A = M * M.Adjoint();
    for (int i = 0; i < n; i++) {
        A(i, i) = A(i, i) + ScalarTraits<T>::fromReal(real_t(n));
    }
    return A;
}

// L 的严格上三角是否精确为 0 (我们从没写过那里, 应当是精确 0)
template <CommutativeScalar T>
static bool UpperStrictZero(const Matrix<T>& L)
{
    using real_t = typename ScalarTraits<T>::real_t;
    for (int i = 0; i < L.Rows(); i++) {
        for (int j = i + 1; j < L.Cols(); j++) {
            if (ScalarTraits<T>::abs2(L(i, j)) != real_t(0)) return false;
        }
    }
    return true;
}

// 对角元是否为正实数: 实部 > 0, 且虚部可忽略
// 注意 "虚部" 必须写成 d - fromReal(Re d) 再取模方, 不能写 abs2(d) - re*re:
//   clang 默认 -ffp-contract=fast, 会把后者融成 fma(-re, re, abs2(d)),
//   结果变成 fl(d^2) - d^2 —— 即 d^2 的舍入残渣(~1e-16), 在 double 上恒大于 tol^2。
//   写成减法则是精确抵消(d - d = 0), 不受 FP contract 影响, 对复数也正好是 Im(d)^2。
template <CommutativeScalar T>
static bool DiagPositiveReal(const Matrix<T>& L, typename ScalarTraits<T>::real_t tol)
{
    using real_t = typename ScalarTraits<T>::real_t;
    for (int i = 0; i < L.Rows(); i++) {
        const T d = L(i, i);
        const real_t re = ScalarTraits<T>::real(d);
        if (!(re > real_t(0))) return false;
        if (ScalarTraits<T>::abs2(d - ScalarTraits<T>::fromReal(re)) > tol * tol) return false;
    }
    return true;
}

// ① 重构判据
// 用 try 包住: SPD 输入本不该抛异常, 但一旦实现写坏就可能抛 —— 那时要把它记成
// 一条 FAIL, 而不是让未捕获异常 abort 掉整个用例 (否则后面的判据一条都跑不到)。
template <CommutativeScalar T>
static void CheckReconstruct(const char* tag, int n, unsigned seed, double& worst)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const std::string name = std::string(tag) + " n=" + std::to_string(n)
                             + " seed=" + std::to_string(seed);
    try {
        Matrix<T> A = RandSPD<T>(n, seed);
        Matrix<T> L = cholesky(A);
        const real_t r = NormInf(L * L.Adjoint() - A);
        if (double(r) > worst) worst = double(r);
        ck.check(r < real_t(V22_TOL), name + ": ‖L·Lᴴ − A‖∞ = " + Fmt(double(r)) + " < 1e-8");
    } catch (const std::exception& e) {
        ck.check(false, name + ": 意外抛出 — " + e.what());
    }
}

// ② V2.2 解方程对拍
template <CommutativeScalar T>
static void CheckAgainstLU(const char* tag, int n, unsigned seed, double& worst_diff, double& worst_res)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const std::string name = std::string(tag) + " n=" + std::to_string(n);
    try {
        Matrix<T> A = RandSPD<T>(n, seed);
        Matrix<T> B = RandMatrix<T>(n, 1, seed + 100u);

        Matrix<T> Xc = cholesky_solve(A, B);
        Matrix<T> Xl = Solve(A, B);

        const real_t diff = NormInf(Xc - Xl);
        const real_t res  = NormInf(A * Xc - B);
        if (double(diff) > worst_diff) worst_diff = double(diff);
        if (double(res)  > worst_res)  worst_res  = double(res);

        ck.check(diff < real_t(V22_TOL),
                 name + ": ‖X_chol − X_lu‖∞ = " + Fmt(double(diff)) + " < 1e-8");
        ck.check(res < real_t(CHOL_RES_TOL),
                 name + ": ‖A·X − B‖∞ = " + Fmt(double(res)) + " < 1e-10");
    } catch (const std::exception& e) {
        ck.check(false, name + ": 解同一方程组意外抛出 — " + e.what());
        ck.check(false, name + ": 残差判据未跑 (同上)");
    }
}

int main()
{
    std::printf("== linalg/cholesky.hpp 验收 (S2a / V2.2) ==\n");

    // ---- A. 重构 (double) ----
    std::printf("\n[A] 重构 ‖L·Lᴴ − A‖∞ (double 随机 SPD)\n");
    double wd = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 4; sd++)
            CheckReconstruct<D>("double", n, sd, wd);

    // ---- B. 重构 (Complex, 同一份源码) ----
    std::printf("\n[B] 重构 ‖L·Lᴴ − A‖∞ (Complex<double>, 同一份 linalg 源码 V1.1)\n");
    double wc = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckReconstruct<C>("complex", n, sd, wc);

    // ---- C. V2.2: 与 LU 解同一方程组 (double) ----
    std::printf("\n[C] V2.2 解同一方程组 vs LU (double)\n");
    double wd_d = 0.0, wd_r = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckAgainstLU<D>("double", n, sd, wd_d, wd_r);

    // ---- D. V2.2: 与 LU 解同一方程组 (Complex) ----
    std::printf("\n[D] V2.2 解同一方程组 vs LU (Complex<double>)\n");
    double wc_d = 0.0, wc_r = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 2; sd++)
            CheckAgainstLU<C>("complex", n, sd, wc_d, wc_r);

    // ---- E. 结构 ----
    std::printf("\n[E] 结构: 严格上三角全 0; 对角元为正实数\n");
    for (int n = 2; n <= 5; n++) {
        Matrix<D> A = RandSPD<D>(n, 3u);
        Matrix<D> L = cholesky(A);
        const D tol = DefaultTol<D>(n);
        ck.check(UpperStrictZero(L), std::string("double n=") + std::to_string(n) + ": 严格上三角全 0");
        ck.check(DiagPositiveReal(L, tol), std::string("double n=") + std::to_string(n) + ": 对角元为正实数");
    }
    {
        Matrix<C> A = RandSPD<C>(3, 3u);
        Matrix<C> L = cholesky(A);
        const D tol = DefaultTol<C>(3);
        ck.check(UpperStrictZero(L), "complex n=3: 严格上三角全 0");
        ck.check(DiagPositiveReal(L, tol), "complex n=3: 对角元为正实数 (虚部可忽略)");
    }

    // ---- F. 解析对拍 ----
    std::printf("\n[F] 解析对拍: 手算 L 逐元一致\n");
    {
        Matrix<D> A2 = {{4.0, 12.0}, {12.0, 37.0}};
        Matrix<D> L2 = cholesky(A2);
        ck.check(std::abs(L2(0, 0) - 2.0) < 1e-13 && std::abs(L2(1, 0) - 6.0) < 1e-13
                 && std::abs(L2(1, 1) - 1.0) < 1e-13 && L2(0, 1) == 0.0,
                 "2x2 [[4,12],[12,37]]: L = [[2,0],[6,1]]");

        Matrix<D> A3 = {{4.0, 12.0, -16.0}, {12.0, 37.0, -43.0}, {-16.0, -43.0, 98.0}};
        Matrix<D> L3 = cholesky(A3);
        const double exp3[3][3] = {{2.0, 0.0, 0.0}, {6.0, 1.0, 0.0}, {-8.0, 5.0, 3.0}};
        bool same = true;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                if (std::abs(L3(i, j) - exp3[i][j]) > 1e-12) same = false;
        ck.check(same, "3x3 经典例: L = [[2,0,0],[6,1,0],[-8,5,3]]");

        // Hermitian 2x2: A = [[4, 2+2i],[2-2i, 6]]  ->  L = [[2,0],[1-i,2]]
        Matrix<C> Ah = {{C(4, 0), C(2, 2)}, {C(2, -2), C(6, 0)}};
        Matrix<C> Lh = cholesky(Ah);
        ck.check(ScalarTraits<C>::abs2(Lh(0, 0) - C(2, 0)) < 1e-26
                 && ScalarTraits<C>::abs2(Lh(1, 0) - C(1, -1)) < 1e-26
                 && ScalarTraits<C>::abs2(Lh(1, 1) - C(2, 0)) < 1e-26
                 && ScalarTraits<C>::abs2(Lh(0, 1)) == 0.0,   // Complex 无 operator==, 判模方为 0
                 "Hermitian 2x2 复: L = [[2,0],[1-i,2]]");
    }

    // ---- G. 接口契约: 只读下三角 ----
    std::printf("\n[G] 接口契约 (Hermitian 由调用方保证)\n");
    {
        Matrix<D> A = RandSPD<D>(4, 41u);
        Matrix<D> G = A;
        for (int i = 0; i < 4; i++)
            for (int j = i + 1; j < 4; j++)
                G(i, j) = 1e9;                     // 上三角填垃圾
        ck.check(IsApprox(cholesky(A), cholesky(G), 1e-14),
                 "只读下三角: 上三角填 1e9 不影响 L (故非 Hermitian 输入不会 throw)");
    }

    // ---- H. 异常路径 ----
    std::printf("\n[H] 异常路径 (D12: throw, 不拿空矩阵当错误码)\n");
    {
        auto throws = [](auto fn) {
            try { fn(); return false; } catch (const std::exception&) { return true; }
        };

        ck.check(throws([&] { (void)cholesky(Matrix<D>(2, 3)); }), "非方阵 -> throw");
        ck.check(throws([&] { (void)cholesky(Matrix<D>(0, 0)); }), "0x0 -> throw");
        ck.check(throws([&] { (void)cholesky(Matrix<D>{{1.0, 2.0}, {2.0, 1.0}}); }),
                 "对称非正定 [[1,2],[2,1]] -> throw");
        ck.check(throws([&] { (void)cholesky(Matrix<D>{{1.0, 1.0}, {1.0, 1.0}}); }),
                 "半正定 [[1,1],[1,1]] -> throw");
        ck.check(throws([&] { (void)cholesky(Matrix<D>{{-1.0, 0.0}, {0.0, -1.0}}); }),
                 "负定 -> throw");

        // 秩亏边界: A = M·Mᵀ (M 3x2) 半正定秩 2 —— 前几步都正常, 最后一个主元才塌成 0
        Matrix<D> M = RandMatrix<D>(3, 2, 33u);
        ck.check(throws([&] { (void)cholesky(M * M.Transpose()); }),
                 "秩亏边界 (A = M·Mᵀ, M 3x2): 最后一个主元才塌 -> throw");

        ck.check(throws([&] { (void)cholesky_solve(Matrix<D>{{1.0, 2.0}, {2.0, 1.0}}, Matrix<D>::Identity(2)); }),
                 "cholesky_solve 非正定 -> throw");
        ck.check(throws([&] { (void)cholesky_solve(Matrix<D>::Identity(3), Matrix<D>::Identity(2)); }),
                 "cholesky_solve B 行数不符 -> throw");
    }

    // ---- I. SPD 判据语义 ----
    std::printf("\n[I] SPD 判据语义: 走通 <=> 正定 (is_positive_definite 的判据本体)\n");
    {
        auto pd = [](auto A) { try { (void)cholesky(A); return true; } catch (...) { return false; } };
        ck.check(pd(RandSPD<D>(4, 21u)),                        "随机 SPD n=4 -> 走通");
        ck.check(pd(Matrix<D>::Identity(5)),                    "单位阵 -> 走通");
        ck.check(!pd(Matrix<D>{{2.0, 0.0}, {0.0, -3.0}}),       "对角含负数 -> 不走通");
        ck.check(!pd(Matrix<D>{{1.0, 3.0}, {3.0, 1.0}}),        "对称非正定 [[1,3],[3,1]] -> 不走通");
        ck.check(pd(RandSPD<C>(3, 22u)),                        "随机 Hermitian 正定 n=3 -> 走通");
    }

    // ---- J. 计时参考 (V2.2 附注; 仅打印, 不作判据) ----
    std::printf("\n[J] 计时参考 (n=160, 仅信息不判分)\n");
    {
        const int n = 160;
        Matrix<D> A = RandSPD<D>(n, 99u);
        Matrix<D> B = RandMatrix<D>(n, 1, 5u);

        auto t0 = std::chrono::steady_clock::now();
        const Matrix<D> Xc = cholesky_solve(A, B);
        auto t1 = std::chrono::steady_clock::now();
        const Matrix<D> Xl = Solve(A, B);
        auto t2 = std::chrono::steady_clock::now();

        const double tc = std::chrono::duration<double, std::milli>(t1 - t0).count();
        const double tl = std::chrono::duration<double, std::milli>(t2 - t1).count();
        std::printf("    cholesky_solve %.2f ms | LU Solve %.2f ms | 比值 %.2f\n", tc, tl, tc / tl);
        std::printf("    理论 flops: cholesky n^3/3, LU 2n^3/3 -> 比约 0.50 (Debug 构建, 只作参考)\n");
        std::printf("    两者解之差 ‖Xc − Xl‖∞ = %s\n", Fmt(double(NormInf(Xc - Xl))).c_str());
    }

    std::printf("\n最差: 重构 double %s | complex %s ; 解之差 double %s | complex %s\n",
                Fmt(wd).c_str(), Fmt(wc).c_str(), Fmt(wd_d).c_str(), Fmt(wc_d).c_str());
    return CheckReport(ck.ok, ck.fail, "cholesky");
}
