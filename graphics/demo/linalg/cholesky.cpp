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
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 迁移前 98 项断言 -> 迁移后 98 条 EXPECT/ADD_FAILURE (逐条对应, 一条没合并)。
//   · 段号 [A]–[I] 原样保留 (J 段是计时参考, 无判据), 通过 SCOPED_TRACE 携带。
//   · 老 helper 里的 `try { ... } catch (std::exception&) { ck.check(false, ...) }`
//     改成 `ADD_FAILURE()`: gtest 会捕获未处理异常并终止**该 TEST 剩余部分**,
//     所以这里保留 try/catch, 让"意外抛出"记为**一条**失败而不吞掉同段其它 seed。
//   · 异常判据走 EXPECT_EXACT_THROW (老写法 catch (std::exception&) 分不出抛哪一种)。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                    | 边界 / 病态输入                        | 专抓什么错
//  ----|-------------------------|----------------------------------------|--------------------------
//   A  | 1e-8 绝对 (相对 ‖A‖≈n)  | n=2..6 × 4 seed, A = M·Mᵀ + n·I         | 用上三角 L 重构 (ᴴ 写成 ᵀ)
//   B  | 1e-8 绝对               | 同上但 Complex<double>                 | 复域漏 conj (实域永远绿)
//   C  | 1e-8 解之差; 1e-10 残差 | n=2..6 × 3 seed, B 单列                | 前代/回代方向错
//   D  | 同上                    | Complex, n=2..4 × 2 seed               | 复域求解漏 conj
//   E  | 严格上三角**精确 0**    | 随机 SPD n=2..5 + complex n=3          | 写出上三角 (重构仍近似对)
//      |                         | 对角元正实部, 虚部 <= DefaultTol        | 对角元取负号/带相位
//   F  | **1e-13 / 1e-12 逐元**  | 2×2 / 3×3 经典例 / Hermitian 2×2        | 算法在"整数可手算"的例子上错
//   G  | 1e-14                   | **上三角填 1e9 垃圾**                   | 偷偷读上三角 (多半顺带 throw)
//   H  | **异常类型精确到子类**  | 非方阵 / 0×0 / 非正定 / 半正定 / 负定    | 用 logic_error 判基类假绿
//      |                         | 秩亏边界 M·Mᵀ (M 3×2) / B 行数不符       | 半正定被当成"能走通"
//   I  | 只判走通/不走通         | 单位阵 / 对角含负 / 对称非正定 / 复 Hermitian| SPD 判据退化成"看对角"
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **L·Lᴴ 写成 L·Lᵀ** (漏 conj): 实域**完全抓不到** ([A]/[C]/[E]/[F] 全绿),
//     只有 [B] (复重构) 与 [D] (复求解) 会红。这就是"复域必须单测"的实证。
//  ② **对角元不取正号** (sqrt 前不加 |·|, 用带符号的余量): [E] 的"对角元为正实数"必红,
//     而 [A] 的重构仍绿 (L 与 diag(±1)L 给出同一个 L·Lᴴ) —— 结构判据的独立价值。
//  ③ **秩亏边界不判 |pivot|** (只判 pivot < 0): [H] 的"秩亏边界"与"半正定"必红,
//     其余段全绿。
//  ④ 计时段 [J] 无论文无判据 (Debug 构建, 只作参考), 不参与注入对照。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh                     # 全量 (注册名: cholesky)
//   ./graphics/build/demo/linalg/cholesky       # 单跑
//   ./graphics/build/demo/linalg/cholesky --gtest_filter='*Contract*'
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/cholesky.hpp"
#include "linalg/solve.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace
{

using D = double;
using C = Complex<double>;

constexpr double V22_TOL = 1e-8;      // V2.2 解之差硬阈值
constexpr double CHOL_RES_TOL = 1e-10;

std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

// SPD 随机矩阵: A = M·Mᴴ + n·I —— 对任意 M 都正定 (特征值 >= n), 且实/复/Dual 通用
template <CommutativeScalar T>
Matrix<T> RandSPD(int n, unsigned seed)
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
bool UpperStrictZero(const Matrix<T>& L)
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
bool DiagPositiveReal(const Matrix<T>& L, typename ScalarTraits<T>::real_t tol)
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
void CheckReconstruct(const char* sec, const char* tag, int n, unsigned seed, double& worst)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const std::string name = std::string("[") + sec + "] " + tag + " n=" + std::to_string(n)
                             + " seed=" + std::to_string(seed);
    try {
        Matrix<T> A = RandSPD<T>(n, seed);
        Matrix<T> L = cholesky(A);
        const real_t r = NormInf(L * L.Adjoint() - A);
        if (double(r) > worst) worst = double(r);
        SCOPED_TRACE(name + ": ‖L·Lᴴ − A‖∞ < 1e-8");
        EXPECT_LT(double(r), V22_TOL) << Fmt(double(r));
    } catch (const std::exception& e) {
        ADD_FAILURE() << name << ": 意外抛出 — " << e.what();
    }
}

// ② V2.2 解方程对拍
template <CommutativeScalar T>
void CheckAgainstLU(const char* sec, const char* tag, int n, unsigned seed,
                    double& worst_diff, double& worst_res)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const std::string name = std::string("[") + sec + "] " + tag + " n=" + std::to_string(n);
    try {
        Matrix<T> A = RandSPD<T>(n, seed);
        Matrix<T> B = RandMatrix<T>(n, 1, seed + 100u);

        Matrix<T> Xc = cholesky_solve(A, B);
        Matrix<T> Xl = Solve(A, B);

        const real_t diff = NormInf(Xc - Xl);
        const real_t res  = NormInf(A * Xc - B);
        if (double(diff) > worst_diff) worst_diff = double(diff);
        if (double(res)  > worst_res)  worst_res  = double(res);

        { SCOPED_TRACE(name + ": ‖X_chol − X_lu‖∞ < 1e-8");
          EXPECT_LT(double(diff), V22_TOL) << Fmt(double(diff)); }
        { SCOPED_TRACE(name + ": ‖A·X − B‖∞ < 1e-10");
          EXPECT_LT(double(res), CHOL_RES_TOL) << Fmt(double(res)); }
    } catch (const std::exception& e) {
        ADD_FAILURE() << name << ": 解同一方程组意外抛出 — " << e.what();
        ADD_FAILURE() << name << ": 残差判据未跑 (同上)";
    }
}

}  // namespace

#define STEP(label) SCOPED_TRACE(label)

// ---- [A] 重构 (double) ----
TEST(Cholesky, ReconstructDouble)
{
    double wd = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 4; sd++)
            CheckReconstruct<D>("A", "double", n, sd, wd);
    std::printf("[cholesky] 重构 double 最差 %s\n", Fmt(wd).c_str());
}

// ---- [B] 重构 (Complex, 同一份源码) ----
TEST(Cholesky, ReconstructComplex)
{
    double wc = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckReconstruct<C>("B", "complex", n, sd, wc);
    std::printf("[cholesky] 重构 complex 最差 %s\n", Fmt(wc).c_str());
}

// ---- [C] V2.2: 与 LU 解同一方程组 (double) ----
TEST(Cholesky, SolveAgainstLU)
{
    double wd_d = 0.0, wd_r = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckAgainstLU<D>("C", "double", n, sd, wd_d, wd_r);
    std::printf("[cholesky] 解之差 double %s | 残差 %s\n", Fmt(wd_d).c_str(), Fmt(wd_r).c_str());
}

// ---- [D] V2.2: 与 LU 解同一方程组 (Complex) ----
TEST(Cholesky, SolveAgainstLUComplex)
{
    double wc_d = 0.0, wc_r = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 2; sd++)
            CheckAgainstLU<C>("D", "complex", n, sd, wc_d, wc_r);
    std::printf("[cholesky] 解之差 complex %s | 残差 %s\n", Fmt(wc_d).c_str(), Fmt(wc_r).c_str());
}

// ---- [E] 结构 ----
TEST(Cholesky, Structure)
{
    for (int n = 2; n <= 5; n++) {
        Matrix<D> A = RandSPD<D>(n, 3u);
        Matrix<D> L = cholesky(A);
        const D tol = DefaultTol<D>(n);
        const std::string tag = "[E] double n=" + std::to_string(n);
        { STEP(tag + ": 严格上三角全 0"); EXPECT_TRUE(UpperStrictZero(L)); }
        { STEP(tag + ": 对角元为正实数"); EXPECT_TRUE(DiagPositiveReal(L, tol)); }
    }
    {
        Matrix<C> A = RandSPD<C>(3, 3u);
        Matrix<C> L = cholesky(A);
        const D tol = DefaultTol<C>(3);
        { STEP("[E] complex n=3: 严格上三角全 0"); EXPECT_TRUE(UpperStrictZero(L)); }
        { STEP("[E] complex n=3: 对角元为正实数 (虚部可忽略)");
          EXPECT_TRUE(DiagPositiveReal(L, tol)); }
    }
}

// ---- [F] 解析对拍 ----
TEST(Cholesky, AnalyticalCompare)
{
    { STEP("[F] 2x2 [[4,12],[12,37]]: L = [[2,0],[6,1]]");
      Matrix<D> A2 = {{4.0, 12.0}, {12.0, 37.0}};
      Matrix<D> L2 = cholesky(A2);
      EXPECT_TRUE(std::abs(L2(0, 0) - 2.0) < 1e-13 && std::abs(L2(1, 0) - 6.0) < 1e-13
                  && std::abs(L2(1, 1) - 1.0) < 1e-13 && L2(0, 1) == 0.0)
          << "L = [[" << L2(0, 0) << ", " << L2(0, 1) << "], ["
          << L2(1, 0) << ", " << L2(1, 1) << "]]"; }

    { STEP("[F] 3x3 经典例: L = [[2,0,0],[6,1,0],[-8,5,3]]");
      Matrix<D> A3 = {{4.0, 12.0, -16.0}, {12.0, 37.0, -43.0}, {-16.0, -43.0, 98.0}};
      Matrix<D> L3 = cholesky(A3);
      const double exp3[3][3] = {{2.0, 0.0, 0.0}, {6.0, 1.0, 0.0}, {-8.0, 5.0, 3.0}};
      bool same = true;
      for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
              if (std::abs(L3(i, j) - exp3[i][j]) > 1e-12) same = false;
      EXPECT_TRUE(same); }

    // Hermitian 2x2: A = [[4, 2+2i],[2-2i, 6]]  ->  L = [[2,0],[1-i,2]]
    { STEP("[F] Hermitian 2x2 复: L = [[2,0],[1-i,2]]");
      Matrix<C> Ah = {{C(4, 0), C(2, 2)}, {C(2, -2), C(6, 0)}};
      Matrix<C> Lh = cholesky(Ah);
      EXPECT_TRUE(ScalarTraits<C>::abs2(Lh(0, 0) - C(2, 0)) < 1e-26
                  && ScalarTraits<C>::abs2(Lh(1, 0) - C(1, -1)) < 1e-26
                  && ScalarTraits<C>::abs2(Lh(1, 1) - C(2, 0)) < 1e-26
                  && ScalarTraits<C>::abs2(Lh(0, 1)) == 0.0); }   // Complex 无 ==, 判模方为 0
}

// ---- [G] 接口契约: 只读下三角 ----
TEST(Cholesky, ReadLowerOnly)
{
    Matrix<D> A = RandSPD<D>(4, 41u);
    Matrix<D> G = A;
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
            G(i, j) = 1e9;                     // 上三角填垃圾
    { STEP("[G] 只读下三角: 上三角填 1e9 不影响 L (故非 Hermitian 输入不会 throw)");
      EXPECT_TRUE(IsApprox(cholesky(A), cholesky(G), 1e-14))
          << "‖L(A) − L(G)‖∞ = " << NormInf(cholesky(A) - cholesky(G)); }
}

// ---- [H] 异常路径 (D12: throw, 不拿空矩阵当错误码) ----
TEST(Cholesky, Contract)
{
    { STEP("[H] 非方阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)cholesky(Matrix<D>(2, 3)), std::invalid_argument); }

    { STEP("[H] 0x0 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)cholesky(Matrix<D>(0, 0)), std::invalid_argument); }

    { STEP("[H] 对称非正定 [[1,2],[2,1]] -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)cholesky(Matrix<D>{{1.0, 2.0}, {2.0, 1.0}}), std::runtime_error); }

    { STEP("[H] 半正定 [[1,1],[1,1]] -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)cholesky(Matrix<D>{{1.0, 1.0}, {1.0, 1.0}}), std::runtime_error); }

    { STEP("[H] 负定 -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)cholesky(Matrix<D>{{-1.0, 0.0}, {0.0, -1.0}}), std::runtime_error); }

    // 秩亏边界: A = M·Mᵀ (M 3x2) 半正定秩 2 —— 前几步都正常, 最后一个主元才塌成 0
    { STEP("[H] 秩亏边界 (A = M·Mᵀ, M 3x2): 最后一个主元才塌 -> throw std::runtime_error");
      Matrix<D> M = RandMatrix<D>(3, 2, 33u);
      EXPECT_EXACT_THROW((void)cholesky(M * M.Transpose()), std::runtime_error); }

    { STEP("[H] cholesky_solve 非正定 -> throw std::runtime_error");
      EXPECT_EXACT_THROW(
          (void)cholesky_solve(Matrix<D>{{1.0, 2.0}, {2.0, 1.0}}, Matrix<D>::Identity(2)),
          std::runtime_error); }

    { STEP("[H] cholesky_solve B 行数不符 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW(
          (void)cholesky_solve(Matrix<D>::Identity(3), Matrix<D>::Identity(2)),
          std::invalid_argument); }
}

// ---- [I] SPD 判据语义 ----
TEST(Cholesky, PositiveDefiniteSemantics)
{
    auto pd = [](auto A) { try { (void)cholesky(A); return true; } catch (...) { return false; } };

    { STEP("[I] 随机 SPD n=4 -> 走通"); EXPECT_TRUE(pd(RandSPD<D>(4, 21u))); }
    { STEP("[I] 单位阵 -> 走通"); EXPECT_TRUE(pd(Matrix<D>::Identity(5))); }
    { STEP("[I] 对角含负数 -> 不走通"); EXPECT_FALSE(pd(Matrix<D>{{2.0, 0.0}, {0.0, -3.0}})); }
    { STEP("[I] 对称非正定 [[1,3],[3,1]] -> 不走通");
      EXPECT_FALSE(pd(Matrix<D>{{1.0, 3.0}, {3.0, 1.0}})); }
    { STEP("[I] 随机 Hermitian 正定 n=3 -> 走通"); EXPECT_TRUE(pd(RandSPD<C>(3, 22u))); }
}

// ---- [J] 计时参考 (V2.2 附注; 仅打印, 不作判据) ----
TEST(Cholesky, TimingReference)
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
    std::printf("[cholesky][J] cholesky_solve %.2f ms | LU Solve %.2f ms | 比值 %.2f\n", tc, tl, tc / tl);
    std::printf("[cholesky][J] 理论 flops: cholesky n^3/3, LU 2n^3/3 -> 比约 0.50 (Debug 构建, 只作参考)\n");
    std::printf("[cholesky][J] 两者解之差 ‖Xc − Xl‖∞ = %s\n", Fmt(double(NormInf(Xc - Xl))).c_str());
}