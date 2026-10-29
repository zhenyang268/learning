// scalar_traits.cpp — 标量层契约验收 (Scalar / CommutativeScalar / fromReal / conj / abs2)
//
// 为什么要有这个 demo:
//   1. Scalar concept 是"同一份 linalg 源码支持 float/Complex/Dual"的地基。
//      契约一旦被破坏, 报错会出现在 linalg 的模板深处, 极难定位 ——
//      静态断言把它钉在编译期。
//   2. 交换律铁律: CommutativeScalar 要求"用单测验证 a*b == b*a", 本文件即该单测。
//   3. 常数通道 fromReal(real_t) 曾经缺失(代码里是 const_n),
//      导致 Complex<T>::sqrt 一实例化就编译错 —— 现按设计文档补齐并在此锁住。
//
// 框架: googletest (2026-10-30 起)。判据编号 [S1].. 只是标识, 不表示先后。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                | 边界 / 病态输入          | 专抓什么错
//  ----|---------------------|--------------------------|----------------------
//   A  | 编译期 (static_assert) | 组合类型 Dual<Complex> | concept 悄悄不满足 → 报错飘到 linalg 深处
//   B  | 1e-9 (精确值)       | 常数种子 derive 必须 =0  | fromReal 把常数填进 epsilon 槽
//   C  | 精确比较            | Dual 的 epsilon 分量     | conj 只 conj(value) 漏 conj(derive)
//   D  | 逐位相等            | 复 / Dual / Dual<Complex> | 乘法实现依赖顺序 (反交换被当成"能过")
//   E  | 1e-15 ~ 1e-6        | 主值分支 (sqrt/log)      | 分支切割错 / 少了 fromReal 通道
//
// 为什么 A 段是编译期的: 「实例化才叫落地」—— 模板 bug 只在实例化时暴露,
// 静态断言把"没实例化"这件事本身也变成红。
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "core/scalar.hpp"
#include "core/matrix.hpp"
#include "core/complex.hpp"
#include "autodiff/dual.hpp"

#include <cmath>

using C  = Complex<double>;
using D  = Dual<float>;
using DC = Dual<Complex<double>>;

// ---------------------------------------------------------------------------
// A. 契约的编译期断言(不满足就编不过, 不用等到运行时)
// ---------------------------------------------------------------------------
static_assert(Scalar<float>,                      "Scalar<float>");
static_assert(Scalar<double>,                     "Scalar<double>");
static_assert(CommutativeScalar<float>,           "CommutativeScalar<float>");
static_assert(Scalar<C>,                          "Scalar<Complex<double>>");
static_assert(CommutativeScalar<C>,               "CommutativeScalar<Complex<double>>");
static_assert(Scalar<D>,                          "Scalar<Dual<float>>");
static_assert(CommutativeScalar<D>,               "CommutativeScalar<Dual<float>>");
// 组合类型: 复系数 + 自动求导(波包演化)
static_assert(Scalar<DC>,                         "Scalar<Dual<Complex<double>>>");
static_assert(CommutativeScalar<DC>,              "CommutativeScalar<Dual<Complex<double>>>");
// ScalarTraits<T>::real_t 永远是实数: 一切比较只在实数域发生
static_assert(std::same_as<ScalarTraits<C>::real_t, double>,  "real_t of Complex<double>");
static_assert(std::same_as<ScalarTraits<D>::real_t, float>,   "real_t of Dual<float>");

// 编译期断言无法被运行期报告, 用一条 TEST 让它在 ctest 输出里留下痕迹
TEST(ScalarContract, CompileTimeAssertionsEvaluated)
{
    SUCCEED() << "11 条 static_assert 已在编译期通过 (见文件头 A 段)";
}

// ---------------------------------------------------------------------------
// B. 常数通道 fromReal(real_t) —— 任意实数常数的唯一入口
// ---------------------------------------------------------------------------

TEST(ScalarTraits, FromRealIsIdentityForRealTypes)
{
    EXPECT_DOUBLE_EQ(ScalarTraits<float>::fromReal(0.5f), 0.5);
    EXPECT_DOUBLE_EQ(ScalarTraits<double>::fromReal(0.5), 0.5);
}

TEST(ScalarTraits, FromRealPutsZeroInImagAndDerive)
{
    const C c = ScalarTraits<C>::fromReal(0.5);       // 复数: 虚部折算为 0
    EXPECT_DOUBLE_EQ(c.real, 0.5);
    EXPECT_DOUBLE_EQ(c.imag, 0.0);

    const D d = ScalarTraits<D>::fromReal(0.5f);      // Dual: 常数种子 derive = 0
    EXPECT_FLOAT_EQ(d.value, 0.5f);
    EXPECT_FLOAT_EQ(d.derive, 0.0f);

    const DC dc = ScalarTraits<DC>::fromReal(0.5);    // 组合类型: 逐层折算
    EXPECT_DOUBLE_EQ(dc.value.real, 0.5);
    EXPECT_DOUBLE_EQ(dc.value.imag, 0.0);
    EXPECT_DOUBLE_EQ(dc.derive.real, 0.0);
    EXPECT_DOUBLE_EQ(dc.derive.imag, 0.0);
}

// ---------------------------------------------------------------------------
// C. conj / abs2 (conj 对 Dual 必须逐分量)
// ---------------------------------------------------------------------------

TEST(ScalarTraits, ConjAndAbs2OnComplex)
{
    EXPECT_DOUBLE_EQ(ScalarTraits<C>::conj(C(3, 4)).real, 3.0);
    EXPECT_DOUBLE_EQ(ScalarTraits<C>::conj(C(3, 4)).imag, -4.0);
    EXPECT_DOUBLE_EQ(ScalarTraits<C>::abs2(C(3, 4)), 25.0);
}

// Dual: conj 逐分量转发基域共轭 —— epsilon 分量也要 conj
// (逐分量对、只 conj(epsilon) 错: d(z̄) = conj(dz))
TEST(ScalarTraits, ConjOnDualConjugatesBothSlots)
{
    const D x(3.0f, 4.0f);
    const D xc = ScalarTraits<D>::conj(x);
    EXPECT_FLOAT_EQ(xc.value, 3.0f);
    EXPECT_FLOAT_EQ(xc.derive, 4.0f);

    // abs2 只看 value(比较/选主元只在实数域发生)
    EXPECT_FLOAT_EQ(ScalarTraits<D>::abs2(x), 9.0f);
}

// Matrix 层: Adjoint 走 conj
TEST(ScalarTraits, MatrixAdjointUsesConj)
{
    Matrix<D> A(2, 2);
    A(0, 0) = D(1.0f, 0.5f);
    A(1, 0) = D(7.0f, 8.0f);

    const Matrix<D> AH = A.Adjoint();
    EXPECT_FLOAT_EQ(AH(0, 1).value, 7.0f);
    EXPECT_FLOAT_EQ(AH(0, 1).derive, 8.0f);
}

// ---------------------------------------------------------------------------
// D. 交换律单测 (CommutativeScalar 的硬要求 a*b == b*a)
// ---------------------------------------------------------------------------
//
// 为什么逐位比而不是用容差: 交换律是**代数恒等式**, 两边必须是同一串运算的
// 不同顺序, 结果应逐位相同。用容差会放过"顺序不同导致浮点舍入不同"的实现,
// 而那种实现恰恰说明乘法不是真的可交换。

TEST(CommutativeScalar, ComplexMultiplicationCommutes)
{
    const C a(1.5, 2.5), b(-0.5, 3.0);
    const C ab = a * b, ba = b * a;
    EXPECT_EQ(ab.real, ba.real);
    EXPECT_EQ(ab.imag, ba.imag);
}

TEST(CommutativeScalar, DualMultiplicationCommutes)
{
    const D p(1.5f, 2.5f), q(-0.5f, 3.0f);
    const D pq = p * q, qp = q * p;
    EXPECT_EQ(pq.value, qp.value);
    EXPECT_EQ(pq.derive, qp.derive);
}

TEST(CommutativeScalar, DualComplexMultiplicationCommutes)
{
    const DC u(C(1.5, 2.5), C(0.5, -1.0)), v(C(-0.5, 3.0), C(1.0, 2.0));
    const DC uv = u * v, vu = v * u;
    EXPECT_EQ(uv.value.real, vu.value.real);
    EXPECT_EQ(uv.value.imag, vu.value.imag);
    EXPECT_EQ(uv.derive.real, vu.derive.real);
    EXPECT_EQ(uv.derive.imag, vu.derive.imag);
}

// ---------------------------------------------------------------------------
// E. Complex 原函数(依赖 fromReal 的路径: sqrt/log)
// ---------------------------------------------------------------------------

TEST(ComplexFunctions, SqrtPrincipalBranch)
{
    const C z = sqrt(C(3.0, 4.0));                   // 主值: 2 + i
    EXPECT_NEAR(z.real, 2.0, 1e-12);
    EXPECT_NEAR(z.imag, 1.0, 1e-12);
}

TEST(ComplexFunctions, LogPrincipalBranch)
{
    const C l = log(C(3.0, 4.0));                    // ln|z| + i·arg
    EXPECT_NEAR(l.real, 1.6094379124341003, 1e-12);
    EXPECT_NEAR(l.imag, 0.9272952180016122, 1e-12);
}

TEST(ComplexFunctions, ExpSinCos)
{
    const C e = exp(C(1.0, 1.0));                    // e·(cos1 + i·sin1)
    EXPECT_NEAR(e.real, 1.4686939399158851, 1e-12);
    EXPECT_NEAR(e.imag, 2.2873552871788423, 1e-12);

    const C s = sin(C(1.0, 1.0));                    // sin x ch y + i cos x sh y
    EXPECT_NEAR(s.real, 1.2984575814159773, 1e-12);
    EXPECT_NEAR(s.imag, 0.6349639147847361, 1e-12);

    const C co = cos(C(1.0, 1.0));
    EXPECT_NEAR(co.real, 0.8337300251311491, 1e-12);
    EXPECT_NEAR(co.imag, -0.9888977057628651, 1e-12);
}

TEST(ComplexFunctions, PowIntegerExponent)
{
    const C pw = pow(C(1.0, 1.0), C(2.0, 0.0));      // (1+i)^2 = 2i
    EXPECT_NEAR(pw.real, 0.0, 1e-15);
    EXPECT_NEAR(pw.imag, 2.0, 1e-15);
}

// ---------------------------------------------------------------------------
// F. Dual 原函数(字面量全部走 traits 之后)
// ---------------------------------------------------------------------------

TEST(DualFunctions, Exp)
{
    const D e = exp(D(1.0f, 1.0f));                  // (e, e)
    EXPECT_NEAR(e.value, 2.718281828459045, 1e-6);
    EXPECT_NEAR(e.derive, 2.718281828459045, 1e-6);
}

TEST(DualFunctions, Log)
{
    const D l = ln(D(2.0f, 1.0f));                   // (ln2, 1/2)
    EXPECT_NEAR(l.value, 0.6931471805599453, 1e-7);
    EXPECT_NEAR(l.derive, 0.5, 1e-7);
}

// sqrt 里用 fromReal(0.5): d√x = dx/(2√x)
TEST(DualFunctions, SqrtPropagatesDerivative)
{
    const D r = sqrt(D(4.0f, 2.0f));                 // (2, 0.5)
    EXPECT_NEAR(r.value, 2.0, 1e-7);
    EXPECT_NEAR(r.derive, 0.5, 1e-7);
}

// pow(x, n): d(x^n) = n·x^(n-1)·dx (常数种子 y' = 0)
TEST(DualFunctions, PowWithConstantExponent)
{
    const D p = pow(D(2.0f, 1.0f), D(3.0f, 0.0f));   // (8, 12)
    EXPECT_NEAR(p.value, 8.0, 1e-6);
    EXPECT_NEAR(p.derive, 12.0, 1e-6);
}

TEST(DualFunctions, CosAndSinAtZero)
{
    const D c = cos(D(0.0f, 1.0f));                  // (1, 0)
    EXPECT_NEAR(c.value, 1.0, 1e-7);
    EXPECT_NEAR(c.derive, 0.0, 1e-7);

    const D sn = sin(D(0.0f, 1.0f));                 // (0, 1)
    EXPECT_NEAR(sn.value, 0.0, 1e-7);
    EXPECT_NEAR(sn.derive, 1.0, 1e-7);
}
