// apply.cpp -- autodiff `Apply` 的域分裂验收 + Function 节点数值 + 嵌套对偶修复
//
// ---------------------------------------------------------------------------
// 为什么单独要这个用例 (2026-10-16)
// ---------------------------------------------------------------------------
// `Apply` 此前只有一个 `Apply(Matrix<T>, ...)`, 体内却假定 T = Dual<U> —— 好处是
// compute_graph (Dual<float> 域) 一直能过, 代价是**标量域从未被实例化过**:
// `Apply(Matrix<double>, Sin)` 报 `error: use of undeclared identifier 'dualFunc'`
// 而**没有任何用例会去调它**。这正是 AGENTS.md §3.1 那类问题 ("文件在、编译器没看过")
// 的又一变体 —— 只不过这次的盲区是"域", 不是"文件"。
//
// 现 `Apply` 按域拆成 4 个重载 (Dual 专用 x2 + 标量专用 x2, 见 jacobian.hpp),
// 靠 [temp.deduct.partial] 部分排序分流。本用例的职责是:
//   [A] 4 个重载 x 4 个域全部**被实例化**且数值对
//   [B] Dual 域 derive 与**解析导数**逐条对拍 (权威判据)
//   [C] Dual 域 derive 与**中心差分**对拍 (判据独立于实现路径)
//   [D] grad_graph 真实路径: Function / MatrixFunction 两个 GradOp (此前从未跑过)
//   [E] 嵌套对偶 Dual<Dual<double>> 的 ln / pow —— 修 dual.hpp 的 log ADL 断链
//   [F] 维度不匹配 -> throw
//
// ---------------------------------------------------------------------------
// 本用例覆盖不到的: 偏序守卫的反向证据
// ---------------------------------------------------------------------------
// jacobian.hpp 的标量版体内有 `static_assert(!IsDual<T>::value, ...)`。要证明它
// **真的在拦**, 必须让 Dual 域落到标量版 —— 那一定会编译失败, 所以属**注入对照**,
// 不能放进正常用例。照抄 smoke.cpp 的两条注入纪律:
//   ① 行为型: 把 (2) 的 static_assert 改成 `!!IsDual<T>::value` -> 本用例编不过
//   ② 编译型: 把 (1) 的签名 `Matrix<Dual<T>>` 误写成 `Matrix<T>` -> 双域调用即报
//      "static_assert failed: Apply: Dual 域应由 Dual 专用重载接管"
// (2026-10-16 两条均已人工验证, 见 dev-progress.md §6)
//
// 运行: ./graphics/run_tests.sh      (注册名: apply)
//       ./graphics/build/demo/autodiff/apply

#include "../common.hpp"
#include "autodiff/dual.hpp"
#include "autodiff/grad_graph.hpp"
#include "autodiff/jacobian.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using D  = Dual<double>;
using CD = Dual<Complex<double>>;

static Checker ck;

static const FunctionType kAll[5] = {
    FunctionType::Sin, FunctionType::Cos, FunctionType::Exp,
    FunctionType::Log, FunctionType::Sqrt
};

// 取样点: 3 个**正数** —— 5 个函数 (含 Log / Sqrt) 在此都有定义
static const double kXs[3] = {0.7, 2.0, 3.5};

static const char* FName(FunctionType f)
{
    switch (f) {
        case FunctionType::Sin:  return "sin";
        case FunctionType::Cos:  return "cos";
        case FunctionType::Exp:  return "exp";
        case FunctionType::Log:  return "log";
        case FunctionType::Sqrt: return "sqrt";
    }
    return "?";
}

// 解析值 / 解析导数 —— 判据来源, 与被测实现**完全独立** (走 libc 的 std::*)
static double AnalyticValue(FunctionType f, double x)
{
    switch (f) {
        case FunctionType::Sin:  return std::sin(x);
        case FunctionType::Cos:  return std::cos(x);
        case FunctionType::Exp:  return std::exp(x);
        case FunctionType::Log:  return std::log(x);
        case FunctionType::Sqrt: return std::sqrt(x);
    }
    return x;
}
static double AnalyticDeriv(FunctionType f, double x)
{
    switch (f) {
        case FunctionType::Sin:  return  std::cos(x);
        case FunctionType::Cos:  return -std::sin(x);
        case FunctionType::Exp:  return  std::exp(x);
        case FunctionType::Log:  return  1.0 / x;
        case FunctionType::Sqrt: return  0.5 / std::sqrt(x);
    }
    return 0.0;
}

// ===========================================================================
// [A] 4 个重载 x 4 个域: 实例化 + 数值
// ===========================================================================
static void TestA_AllOverloadsAllDomains()
{
    std::printf("\n[A] 4 个重载 x {double, Complex, Dual, Dual<Complex>}\n");
    const double tol = 1e-15;

    // --- (2) 标量专用: double ---
    {
        Matrix<double> A(3, 1);
        for (int j = 0; j < 3; j++) A(j, 0) = kXs[j];
        bool ok = true;
        for (FunctionType f : kAll) {
            auto B = Apply(A, f);
            for (int j = 0; j < 3; j++) {
                if (std::fabs(B(j, 0) - AnalyticValue(f, kXs[j])) > tol) ok = false;
            }
        }
        ck.check(ok, std::string("[A1] (2) Matrix<double> x FunctionType, 5 函数 x 3 点"));
    }

    // --- (4) 标量专用: double + Matrix<FunctionType> (每格不同函数) ---
    {
        Matrix<double> A(3, 1);
        for (int j = 0; j < 3; j++) A(j, 0) = kXs[j];
        Matrix<FunctionType> F(3, 1, std::vector<FunctionType>{
            FunctionType::Sin, FunctionType::Exp, FunctionType::Log});
        auto B = Apply(A, F);
        bool ok = std::fabs(B(0, 0) - std::sin(0.7)) < tol
               && std::fabs(B(1, 0) - std::exp(2.0)) < tol
               && std::fabs(B(2, 0) - std::log(3.5)) < tol;
        ck.check(ok, std::string("[A2] (4) Matrix<double> x Matrix<FunctionType> 逐格不同"));
    }

    // --- (1) Dual 专用: value 与 derive 都要对 ---
    {
        Matrix<D> A(3, 1);
        for (int j = 0; j < 3; j++) A(j, 0) = D(kXs[j], 1.0);   // seed = 1
        bool ok = true;
        for (FunctionType f : kAll) {
            auto B = Apply(A, f);
            for (int j = 0; j < 3; j++) {
                if (std::fabs(B(j, 0).value  - AnalyticValue(f, kXs[j])) > tol) ok = false;
                if (std::fabs(B(j, 0).derive - AnalyticDeriv(f, kXs[j])) > 1e-14) ok = false;
            }
        }
        ck.check(ok, std::string("[A3] (1) Matrix<Dual<double>> x FunctionType, value+derive"));
    }

    // --- (3) Dual 专用 + Matrix<FunctionType> ---
    {
        Matrix<D> A(2, 1);
        A(0, 0) = D(0.7, 1.0); A(1, 0) = D(1.2, 1.0);
        Matrix<FunctionType> F(2, 1, std::vector<FunctionType>{
            FunctionType::Sin, FunctionType::Exp});
        auto B = Apply(A, F);
        bool ok = std::fabs(B(0, 0).value  - std::sin(0.7)) < tol
               && std::fabs(B(0, 0).derive - std::cos(0.7)) < 1e-14
               && std::fabs(B(1, 0).value  - std::exp(1.2)) < tol
               && std::fabs(B(1, 0).derive - std::exp(1.2)) < 1e-14;
        ck.check(ok, std::string("[A4] (3) Matrix<Dual<double>> x Matrix<FunctionType>"));
    }

    // --- (2) 标量专用: Complex< double>, 与 std::complex 对拍 (独立实现) ---
    {
        const std::complex<double> z(0.5, 0.2);
        Matrix<Complex<double>> A(1, 1);
        A(0, 0) = Complex<double>(z.real(), z.imag());
        bool ok = true;
        for (FunctionType f : kAll) {
            std::complex<double> want;
            switch (f) {
                case FunctionType::Sin:  want = std::sin(z);  break;
                case FunctionType::Cos:  want = std::cos(z);  break;
                case FunctionType::Exp:  want = std::exp(z);  break;
                case FunctionType::Log:  want = std::log(z);  break;
                case FunctionType::Sqrt: want = std::sqrt(z); break;
            }
            auto B = Apply(A, f);
            if (std::fabs(B(0, 0).real - want.real()) > 1e-15) ok = false;
            if (std::fabs(B(0, 0).imag - want.imag()) > 1e-15) ok = false;
        }
        ck.check(ok, std::string("[A5] (2) Matrix<Complex<double>>, 对拍 std::complex"));
    }

    // --- (1) Dual 专用: Dual<Complex<double>> ---
    {
        const std::complex<double> z(0.5, 0.2);
        Matrix<CD> A(1, 1);
        A(0, 0) = CD(Complex<double>(z.real(), z.imag()), Complex<double>(1.0, 0.0));
        auto B = Apply(A, FunctionType::Sin);
        std::complex<double> want = std::sin(z);
        bool ok = std::fabs(B(0, 0).value.real - want.real()) < 1e-15
               && std::fabs(B(0, 0).value.imag - want.imag()) < 1e-15;
        ck.check(ok, std::string("[A6] (1) Matrix<Dual<Complex<double>>> 走通"));
    }
}

// ===========================================================================
// [B] Dual 域 derive vs 解析导数: 权威判据 (上面 [A3] 已含, 这里单列宽覆盖)
// ===========================================================================
static void TestB_AnalyticDerivWide()
{
    std::printf("\n[B] derive vs 解析导数 (10 点 x 5 函数)\n");
    const int nx = 10;
    bool ok = true;
    double worst = 0.0;
    double perF[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (int j = 0; j < nx; j++) {
        const double x = 0.25 + 0.3 * j;          // 0.25 .. 2.95, 全正
        Matrix<D> A(1, 1);
        A(0, 0) = D(x, 1.0);
        for (int k = 0; k < 5; k++) {
            const FunctionType f = kAll[k];
            auto B = Apply(A, f);
            const double err = std::fabs(B(0, 0).derive - AnalyticDeriv(f, x));
            if (err > perF[k]) perF[k] = err;
            if (err > worst)   worst   = err;
            if (err > 1e-14)   ok      = false;
        }
    }
    for (int k = 0; k < 5; k++)
        std::printf("     %-5s 最差 |derive - 解析| = %.3e\n", FName(kAll[k]), perF[k]);
    std::printf("     合计最差 = %.3e\n", worst);
    ck.check(ok, std::string("[B1] 50 组 (10 点 x 5 函数) 全部 <= 1e-14"));
}

// ===========================================================================
// [C] Dual 域 derive vs 中心差分: 判据独立于实现路径
// ===========================================================================
static void TestC_CentralDifference()
{
    std::printf("\n[C] derive vs 中心差分 (独立判据)\n");
    const double h = 1e-6;
    bool ok = true;
    double worst = 0.0;
    for (double x : {0.7, 1.4, 2.3, 3.1}) {
        for (FunctionType f : kAll) {
            Matrix<D> A(1, 1); A(0, 0) = D(x, 1.0);
            Matrix<D> P(1, 1); P(0, 0) = D(x + h, 1.0);
            Matrix<D> M(1, 1); M(0, 0) = D(x - h, 1.0);
            const double auto_d = Apply(A, f)(0, 0).derive;
            const double num_d  = (Apply(P, f)(0, 0).value - Apply(M, f)(0, 0).value) / (2.0 * h);
            const double err = std::fabs(auto_d - num_d);
            // 阈值取 1e-8: 中心差分截断 O(h^2)=1e-12, 舍入 ~eps/h ≈ 1e-10, 都在其下
            if (err > worst) worst = err;
            if (err > 1e-8) ok = false;
        }
    }
    std::printf("     最差 |自动微分 - 中心差分| = %.3e\n", worst);
    ck.check(ok, std::string("[C1] 20 组 (4 点 x 5 函数) 全部 <= 1e-8"));
}

// ===========================================================================
// [D] grad_graph 真实路径: Function / MatrixFunction (改用例里从未跑过)
// ===========================================================================
static void TestD_GradGraphFunctionNodes()
{
    std::printf("\n[D] grad_graph: Function / MatrixFunction 节点 (真实调用路径)\n");

    // Function 节点: y = sin(x), x = 0.7
    {
        Matrix<FunctionType> F(1, 1, std::vector<FunctionType>{FunctionType::Sin});
        auto fn = std::make_shared<GradNode<D>>(GradOp::Function, std::move(F));
        Matrix<D> m(1, 1); m(0, 0) = D(0.7, 1.0);
        auto in = std::make_shared<GradNode<D>>(GradOp::Input, std::move(m));
        fn->FillInput(in);
        fn->Forward();
        bool ok = std::fabs(fn->matrix(0, 0).value  - std::sin(0.7)) < 1e-15
               && std::fabs(fn->matrix(0, 0).derive - std::cos(0.7)) < 1e-14;
        ck.check(ok, std::string("[D1] Function(sin) 前向 value + derive"));
    }

    // MatrixFunction 节点: y = [sin(x0), exp(x1)]
    {
        Matrix<FunctionType> F(2, 1, std::vector<FunctionType>{
            FunctionType::Sin, FunctionType::Exp});
        auto fn = std::make_shared<GradNode<D>>(GradOp::MatrixFunction, std::move(F));
        Matrix<D> m(2, 1);
        m(0, 0) = D(0.7, 1.0); m(1, 0) = D(1.2, 1.0);
        auto in = std::make_shared<GradNode<D>>(GradOp::Input, std::move(m));
        fn->FillInput(in);
        fn->Forward();
        bool ok = std::fabs(fn->matrix(0, 0).value  - std::sin(0.7)) < 1e-15
               && std::fabs(fn->matrix(0, 0).derive - std::cos(0.7)) < 1e-14
               && std::fabs(fn->matrix(1, 0).value  - std::exp(1.2)) < 1e-15
               && std::fabs(fn->matrix(1, 0).derive - std::exp(1.2)) < 1e-14;
        ck.check(ok, std::string("[D2] MatrixFunction([sin, exp]) 逐格不同 + derive"));
    }

    // MatrixFunction 的退化路径 (1x1 函数矩阵 -> 落回统一应用重载)
    {
        Matrix<FunctionType> F(1, 1, std::vector<FunctionType>{FunctionType::Sqrt});
        auto fn = std::make_shared<GradNode<D>>(GradOp::MatrixFunction, std::move(F));
        Matrix<D> m(1, 1); m(0, 0) = D(2.25, 1.0);
        auto in = std::make_shared<GradNode<D>>(GradOp::Input, std::move(m));
        fn->FillInput(in);
        fn->Forward();
        bool ok = std::fabs(fn->matrix(0, 0).value  - 1.5) < 1e-15
               && std::fabs(fn->matrix(0, 0).derive - (0.5 / 1.5)) < 1e-14;
        ck.check(ok, std::string("[D3] MatrixFunction 退化 (1x1 -> 统一应用) 路径"));
    }
}

// ===========================================================================
// [E] 嵌套对偶 Dual<Dual<double>>: 二阶导 (修 dual.hpp 的 log ADL 断链)
// ===========================================================================
static void TestE_NestedDualSecondOrder()
{
    std::printf("\n[E] 嵌套对偶: ln / pow 的二阶导 (dual.hpp log 断链回归)\n");
    using DD = Dual<Dual<double>>;

    // seed: value.derive = 1 (第一方向), derive.value = 1 (第二方向)
    // 标准 hyper-dual 布局 => derive.derive 即 f''
    const double x0 = 2.0;
    DD x(Dual<double>(x0, 1.0), Dual<double>(1.0, 0.0));

    // ln: f' = 1/x, f'' = -1/x^2
    {
        auto y = ln(x);
        const double v  = y.value.value;
        const double d1 = y.value.derive;
        const double d2 = y.derive.derive;
        std::printf("     ln : v=%.15g d1=%.15g d2=%.15g\n", v, d1, d2);
        bool ok = std::fabs(v  - std::log(x0)) < 1e-15
               && std::fabs(d1 - 1.0 / x0)     < 1e-15
               && std::fabs(d2 - (-1.0 / (x0 * x0))) < 1e-14;
        // 这段在修 log 断链之前**根本编不过** (no matching function for call to 'log')
        ck.check(ok, std::string("[E1] ln(Dual<Dual<double>>): v / f' / f'' 三阶全对"));
    }

    // pow: f = x^3, f' = 3x^2 = 12, f'' = 6x = 12  (x0=2)
    {
        DD y3(Dual<double>(3.0, 0.0), Dual<double>(0.0, 0.0));
        auto z = pow(x, y3);
        const double v  = z.value.value;
        const double d1 = z.value.derive;
        const double d2 = z.derive.derive;
        std::printf("     x^3: v=%.15g d1=%.15g d2=%.15g\n", v, d1, d2);
        bool ok = std::fabs(v  - 8.0)  < 1e-13
               && std::fabs(d1 - 12.0) < 1e-13
               && std::fabs(d2 - 12.0) < 1e-13;
        ck.check(ok, std::string("[E2] pow(Dual<Dual<double>>, 3): v / f' / f'' 三阶全对"));
    }

    // exp / sin 也走一遍 (它们的嵌套链本来就通, 作为对照)
    {
        auto ye = exp(x);
        bool ok_e = std::fabs(ye.value.value - std::exp(x0)) < 1e-13
                 && std::fabs(ye.derive.derive - std::exp(x0)) < 1e-13;
        auto ys = sin(x);
        bool ok_s = std::fabs(ys.value.value - std::sin(x0)) < 1e-15
                 && std::fabs(ys.derive.derive - (-std::sin(x0))) < 1e-14;
        ck.check(ok_e, std::string("[E3] exp 嵌套对照 (链本就通)"));
        ck.check(ok_s, std::string("[E4] sin 嵌套对照 (链本就通)"));
    }
}

// ===========================================================================
// [F] 维度不匹配 -> throw
// ===========================================================================
static void TestF_DimensionMismatch()
{
    std::printf("\n[F] 维度不匹配 -> throw\n");
    // 列向量约定: A 是 3x1, 函数矩阵误给成 2x1 (长度不符)
    Matrix<double> A(3, 1);
    Matrix<FunctionType> F(2, 1, std::vector<FunctionType>{
        FunctionType::Sin, FunctionType::Cos});

    bool threw = false;
    try { (void)Apply(A, F); } catch (const std::invalid_argument&) { threw = true; }
    ck.check(threw, std::string("[F1] 标量版: A 3x1 vs f 2x1 -> invalid_argument"));

    Matrix<D> Ad(3, 1);
    bool threw_u = false;
    try { (void)Apply(Ad, F); } catch (const std::invalid_argument&) { threw_u = true; }
    ck.check(threw_u, std::string("[F2] Dual 版: 同样 throw"));

    // 反向对照: 1x1 函数矩阵应当**广播**, 不抛 (否则上面的 throw 可能只是"见到非同形就抛")
    Matrix<FunctionType> F1(1, 1, std::vector<FunctionType>{FunctionType::Sin});
    bool ok_bcast = false;
    try { (void)Apply(A, F1); ok_bcast = true; } catch (const std::exception&) { ok_bcast = false; }
    ck.check(ok_bcast, std::string("[F3] 1x1 函数矩阵 -> 广播, 不抛"));
}

// ===========================================================================

int main()
{
    TestA_AllOverloadsAllDomains();
    TestB_AnalyticDerivWide();
    TestC_CentralDifference();
    TestD_GradGraphFunctionNodes();
    TestE_NestedDualSecondOrder();
    TestF_DimensionMismatch();
    return CheckReport(ck.ok, ck.fail, "autodiff/apply");
}
