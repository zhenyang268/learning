// dual_complex.cpp -- V1.5 验收: linalg × autodiff 交叉 (Dual<Complex<double>>)
//
// V1.5 原文 :
//   `Dual<Complex<double>>` 作为 Matrix 元素编译通过, 2x2 复矩阵 LU 的 value
//   与直接用 `Complex<double>` 计算一致 (为波包演化铺路)。
//   配套用例 (Complex 上的求导图景):
//     z^2 单种子(全纯) / |z|^2 双种子(非全纯);
//     阈值: value 一致 1e-12, derive 一致 1e-10 (double)。
//
// 这一刀在解决什么 (三层, 不只是"复数求导"):
//   ① 工程: 泛型域覆盖的最后一格。double / Complex / Dual<double> 都验过,
//      Dual<Complex<double>> 从没被编译器碰过 —— "文件在"不等于"能编译"的同族风险。
//   ② 数学: Dual 的规则 (fg)' = f'g + fg' 隐含"导数唯一"。实数域可导即唯一;
//      复数域只有**全纯函数**才有唯一复导数。非全纯函数上 Dual 给的是"沿某个方向的
//      数", 随方向变化 —— 那是错的答案。
//   ③ 用法: 非全纯场合的补救 = 双种子。**一个种子只给出实雅可比 J 的一列**
//      (复值 derive 的实/虚部就是那一列的 2 个元素); 两个线性无关种子 (1 与 i)
//      才凑齐整张 2x2 实 J, 从而能分离出 Wirtinger 的两个分量:
//        ∂f/∂z  = (d0 − i·d1) / 2   —— 恒可定义, 全纯时它就是 f'(z)
//        ∂f/∂z̄ = (d0 + i·d1) / 2   —— **非全纯判据**: 这个 0 不 0
//      全纯时 J 的第 2 列由第 1 列唯一决定 (转 90°: d1 = i·d0), 所以一个种子就够 ——
//      这正是 C-R 方程的雅可比说法。
//
// 判据 (判据先行 —— 先写"多少算对", 再对着判据写实现):
//   [A] 泛型域   Matrix<Dual<Complex<double>>> 构造/随机生成/矩阵乘 零修改可用;
//                d(A^2) = A·A' + A'·A (矩阵层链式法则)
//   [B] 跨域一致 derive 置零时, Dual 矩阵解方程的结果 value 与纯 Complex 逐元一致 < 1e-12
//   [C] 全纯     z^2 的复导数唯一: 两种子 (s=1, s=i) 的 derive/s 相同, 且 = 2z0
//   [D] 非全纯   |z|^2 无复导数: 两方向 derive/s 必然不同 (**本段的真判据**);
//                Wirtinger 对拍: ∂f/∂z̄ = (d0+i·d1)/2 = z0 ≠ 0 (非全纯标志),
//                                ∂f/∂z  = (d0−i·d1)/2 = z̄0  (恒可定义, 不是判据)
//   [E] 实用形态 解 A(theta)x = b 时 dx/dtheta 自动传播, 与解析式 + 中心差分三方对拍
//   [F] 边界     1x1 通过; 非方阵 / 主元近奇异 全部 throw (D12)
//
// 运行 (由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh                   # 全量
//   ./graphics/build/demo/linalg/dual_complex  # 单跑

#include "../common.hpp"
#include "linalg/solve.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

static Checker ck;

using C  = Complex<double>;
using DC = Dual<C>;
using MC = Matrix<C>;
using MD = Matrix<DC>;

// V1.5 原文阈值
constexpr double VALUE_TOL  = 1e-12;
constexpr double DERIVE_TOL = 1e-10;

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

static double Mod(const C& z) { return std::sqrt(z.real * z.real + z.imag * z.imag); }

// 纯复数矩阵 -> Dual 矩阵, derive 全置 0
// (value 应与原矩阵逐位一致; derive 为 0 表示"这不是求导, 只是换个域算同一个答案")
static MD Lift(const MC& A)
{
    MD B(A.Rows(), A.Cols());
    const C z0(0.0, 0.0);
    for (int i = 0; i < A.Rows(); i++)
        for (int j = 0; j < A.Cols(); j++)
            B(i, j) = DC(A(i, j), z0);
    return B;
}

// Dual 矩阵的 value 部分 与 纯复数矩阵 的最大逐元差 (模)
static double ValueErr(const MD& A, const MC& B)
{
    double e = 0.0;
    for (int i = 0; i < A.Rows(); i++)
        for (int j = 0; j < A.Cols(); j++)
            e = std::max(e, Mod(A(i, j).value - B(i, j)));
    return e;
}

// ---- [A] 泛型域: Matrix<Dual<Complex<double>>> 零修改可用 --------------------
// "零修改"的证据: 下面是与 solve.cpp / cholesky.cpp 完全相同的 Matrix 运算,
// 没有任何为 Dual 写的特化分支。
static void RunGeneric()
{
    std::printf("\n[A] 泛型域: Matrix<Dual<Complex<double>>> 基本可用\n");

    // ① 随机生成 —— 顺带验收 common.hpp 的 RandFill<Dual<Complex<U>>> 特化
    MD R = RandMatrix<DC>(3, 3, 41u);
    bool finite = true;
    for (int i = 0; i < R.Size(); i++) {
        const DC& z = R.Data()[i];
        if (!std::isfinite(z.value.real) || !std::isfinite(z.value.imag) ||
            !std::isfinite(z.derive.real) || !std::isfinite(z.derive.imag)) finite = false;
    }
    ck.check(finite, "RandMatrix<Dual<Complex<double>>> 生成 3x3");

    // ② 构造 + 矩阵乘: value 与同 value 的纯 Complex 矩阵一致
    MD A(2, 2);
    A(0, 0) = DC(C(1.0,  2.0), C(0.0,  1.0));
    A(0, 1) = DC(C(0.5, -0.5), C(1.0,  0.0));
    A(1, 0) = DC(C(-1.0, 0.0), C(0.0,  0.0));
    A(1, 1) = DC(C(3.0,  2.0), C(0.0, -1.0));

    MC Av(2, 2), Dv(2, 2);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) { Av(i, j) = A(i, j).value; Dv(i, j) = A(i, j).derive; }

    MD P = A * A;
    MC p = Av * Av;
    const double emul = ValueErr(P, p);
    ck.check(emul < VALUE_TOL, "A·A 的 value 与纯 Complex 一致 (" + Fmt(emul) + ")");

    // ③ 矩阵层链式法则: d(A^2) = A·A' + A'·A  (A' = 逐元素 derive 矩阵)
    MC wantD = Av * Dv + Dv * Av;
    double eder = 0.0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            eder = std::max(eder, Mod(P(i, j).derive - wantD(i, j)));
    ck.check(eder < DERIVE_TOL, "d(A^2) = A·A' + A'·A (" + Fmt(eder) + ")");
}

// ---- [B] 跨域一致: derive 置零时, Dual 解与纯复数解逐元一致 ----------------
static void RunSolveCross(int n, unsigned seed)
{
    MC A = RandMatrix<C>(n, n, seed);
    for (int i = 0; i < n; i++) A(i, i) = A(i, i) + C(double(n) + 1.0, 0.0);  // 压小条件数
    MC b = RandMatrix<C>(n, 1, seed + 100u);

    MC xc = Solve(A, b);
    MD xd = Solve(Lift(A), Lift(b));

    const double e   = ValueErr(xd, xc);
    const double res = NormInf(A * xc - b);

    const std::string tag = "n=" + std::to_string(n) + " ";
    ck.check(e < VALUE_TOL,   tag + "Solve: Dual<Complex> 的 value = 纯 Complex 解 (" + Fmt(e) + ")");
    ck.check(res < 1e-10,     tag + "残差 ‖A·x − b‖∞ (" + Fmt(res) + ")");
}

// ---- [C] 全纯: z^2 的复导数唯一 (两种子必须给同一个 derive/s) ----------------
static void RunHolomorphic()
{
    std::printf("\n[C] 全纯 z^2: 复导数唯一, 与种子方向无关\n");

    const C z0(0.6, -0.8);
    const C seeds[2] = { C(1.0, 0.0), C(0.0, 1.0) };   // 两个线性无关的复方向
    C val[2], ratio[2];

    for (int k = 0; k < 2; k++) {
        DC z(z0, seeds[k]);
        DC w = z * z;
        val[k]   = w.value;
        ratio[k] = w.derive / seeds[k];      // 商: 把种子方向剥掉
    }

    const double eval = Mod(val[0] - z0 * z0);
    ck.check(eval < VALUE_TOL, "z^2 的 value = z0^2 (" + Fmt(eval) + ")");

    const C analytic(2.0 * z0.real, 2.0 * z0.imag);
    const double d1 = Mod(ratio[0] - ratio[1]);
    const double d2 = Mod(ratio[0] - analytic);
    ck.check(d1 < DERIVE_TOL, "两种子 derive/s 一致 — 复导数唯一 (" + Fmt(d1) + ")");
    ck.check(d2 < DERIVE_TOL, "derive/s = 2·z0 解析对拍 (" + Fmt(d2) + ")");
}

// ---- [D] 非全纯: |z|^2 没有复导数, 必须双种子 -------------------------------
static void RunNonHolomorphic()
{
    std::printf("\n[D] 非全纯 |z|^2: 复导数不存在, 双种子凑齐 2x2 实雅可比\n");

    const C z0(0.6, -0.8);
    const C seeds[2] = { C(1.0, 0.0), C(0.0, 1.0) };
    C val[2], ratio[2];

    for (int k = 0; k < 2; k++) {
        DC z(z0, seeds[k]);
        DC f = ScalarTraits<DC>::conj(z) * z;     // |z|^2 = zbar·z
        val[k]   = f.value;
        ratio[k] = f.derive / seeds[k];           // 必须剥掉种子方向才能比 —— 见下
    }

    // value 与种子无关, 恒为 |z0|^2 = 1
    const double v1 = Mod(val[0] - C(1.0, 0.0));
    const double v2 = Mod(val[1] - C(1.0, 0.0));
    ck.check(v1 < VALUE_TOL && v2 < VALUE_TOL, "|z|^2 的 value = 1 (与种子无关)");

    // 归一化后的 "derive/s" 随方向变化 —— 这才是"复导数不存在"的判据。
    // 注意不能用原始 derive 相减: z^2 (全纯) 的两种子态 derive 也恰好不同
    // (2z0 vs 2z0*i), 拿原始值比会漏判 —— 实测把本段换成 z*z 验证过该漏洞。
    const double spread = Mod(ratio[0] - ratio[1]);
    ck.check(spread > 1e-3,
             "两方向 derive/s 不同 (|Δ|=" + Fmt(spread) + ") — 复导数不存在");

    // 双种子重组: comb = d0 − i·d1，按 Wirtinger 它是 2·∂f/∂z。
    // 对 |z|^2 恰好 = 2·z̄0，所以数值上就是那个 2·zbar —— 但**别把它当成非全纯判据**:
    // ∂f/∂z 对任何实可微的 f 都定义, 全纯的 z^2 上它同样非零 (那里 = 2z0)。
    // 「非全纯」只体现在 ∂f/∂z̄ = (d0 + i·d1)/2 ≠ 0 (本段 = z0 ≠ 0),
    // 以及上面那条 spread —— 方向不同, 归一化后的导数就不同。
    const C comb = ratio[0] * seeds[0] - C(0.0, 1.0) * (ratio[1] * seeds[1]);
    const C want(2.0 * z0.real, -2.0 * z0.imag);
    const double dc = Mod(comb - want);
    ck.check(dc < DERIVE_TOL, "d0 − i·d1 = 2·∂f/∂z = 2·z̄0 (Wirtinger 重组) (" + Fmt(dc) + ")");
}

// ---- [E] 实用形态: A(theta)x = b 的 dx/dtheta 自动传播 ----------------------
// 这是 V1.5 "为波包演化铺路"的真实形状: 线性系统对参数的灵敏度由 Dual 直接给出,
// 不需要手推求导公式。
static void RunPractical()
{
    std::printf("\n[E] 实用形态: 解 A(theta)·x = b, 求 dx/dtheta (三方对拍)\n");

    const int n = 3;
    MC A0 = RandMatrix<C>(n, n, 21u);
    for (int i = 0; i < n; i++) A0(i, i) = A0(i, i) + C(double(n) + 2.0, 0.0);
    MC A1 = RandMatrix<C>(n, n, 22u);
    MC b  = RandMatrix<C>(n, 1, 23u);

    const double th0 = 0.37, h = 1e-6;

    // ① 解析: x = A^-1 b,  dx/dtheta = -A^-1·A1·x
    const MC Ath = A0 + A1 * C(th0, 0.0);
    const MC x0  = Solve(Ath, b);
    const MC dxa = Solve(Ath, (A1 * x0) * C(-1.0, 0.0));

    // ② Dual 自动传播: A(theta) 每格写成 Dual, theta 的 derive = 1
    MD Ad(n, n);
    const DC theta(C(th0, 0.0), C(1.0, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Ad(i, j) = DC(A0(i, j), C(0.0, 0.0)) + theta * DC(A1(i, j), C(0.0, 0.0));
    const MD xd = Solve(Ad, Lift(b));

    // ③ 中心差分参考
    const MC xp = Solve(A0 + A1 * C(th0 + h, 0.0), b);
    const MC xm = Solve(A0 + A1 * C(th0 - h, 0.0), b);
    const MC dxf = (xp - xm) * C(1.0 / (2.0 * h), 0.0);

    double e_val = 0.0, e_ana = 0.0, e_fd = 0.0;
    for (int i = 0; i < n; i++) {
        e_val = std::max(e_val, Mod(xd(i, 0).value  - x0(i, 0)));
        e_ana = std::max(e_ana, Mod(xd(i, 0).derive - dxa(i, 0)));
        e_fd  = std::max(e_fd,  Mod(xd(i, 0).derive - dxf(i, 0)));
    }

    ck.check(e_val < VALUE_TOL,  "solve 的 value = x(theta0) (" + Fmt(e_val) + ")");
    ck.check(e_ana < DERIVE_TOL, "自动传播的 dx/dtheta = 解析 −A^-1·A1·x (" + Fmt(e_ana) + ")");
    ck.check(e_fd  < 1e-7,       "自动传播的 dx/dtheta = 中心差分 (" + Fmt(e_fd) + ")");
}

// ---- [F] 边界 --------------------------------------------------------------
static void RunBoundary()
{
    std::printf("\n[F] 边界\n");

    // 1x1 复数
    MC a1(1, 1); a1(0, 0) = C(2.0, 1.0);
    MC b1(1, 1); b1(0, 0) = C(3.0, -1.0);
    const MC x1 = Solve(a1, b1);
    const double e1 = Mod(x1(0, 0) - b1(0, 0) / a1(0, 0));
    ck.check(e1 < VALUE_TOL, "1x1 复数 Solve (" + Fmt(e1) + ")");

    const MD x1d = Solve(Lift(a1), Lift(b1));
    ck.check(Mod(x1d(0, 0).value - x1(0, 0)) < VALUE_TOL, "1x1 Dual<Complex> Solve: value 一致");

    // 非方阵 -> throw
    bool threw = false;
    try {
        MC bad(2, 3);
        MC rhs(2, 1);
        (void)Solve(bad, rhs);
    } catch (const std::exception&) { threw = true; }
    ck.check(threw, "非方阵 Solve 抛异常 (D12)");

    // 主元近奇异 -> throw
    threw = false;
    try {
        MC sing(2, 2);
        sing(0, 0) = C(1.0, 0.0); sing(0, 1) = C(2.0, 0.0);
        sing(1, 0) = C(2.0, 0.0); sing(1, 1) = C(4.0, 0.0);   // 第二行 = 2 × 第一行
        MC rhs(2, 1);
        (void)Solve(sing, rhs);
    } catch (const std::exception&) { threw = true; }
    ck.check(threw, "奇异矩阵 Solve 抛异常 (D12)");
}

int main()
{
    std::printf("=== V1.5: Dual<Complex<double>> 泛型域 + 复数求导 ===\n");
    std::printf("阈值: value < %.0e, derive < %.0e\n", VALUE_TOL, DERIVE_TOL);

    // 每段独立包 try —— 实现抛异常不能 terminate 掉整个用例 (S2a/S3 各踩过一次 SIGABRT)
    auto guard = [](const char* name, void (*fn)()) {
        try {
            fn();
        } catch (const std::exception& e) {
            ck.fail++;
            std::printf("  [FAIL] %s: 抛出异常 \"%s\"\n", name, e.what());
        }
    };

    guard("A 泛型域",   RunGeneric);
    guard("B 跨域一致", [] { RunSolveCross(2, 31u); RunSolveCross(4, 32u); RunSolveCross(6, 33u); });
    guard("C 全纯",     RunHolomorphic);
    guard("D 非全纯",   RunNonHolomorphic);
    guard("E 实用形态", RunPractical);
    guard("F 边界",     RunBoundary);

    return CheckReport(ck.ok, ck.fail, "dual_complex");
}
