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
// 编译与运行(由 CMake 注册为 ctest 用例):
//   cmake -S graphics -B graphics/build -DCMAKE_BUILD_TYPE=Debug
//   cmake --build graphics/build -j
//   ctest --test-dir graphics/build --output-on-failure

#include "core/scalar.hpp"
#include "core/matrix.hpp"
#include "core/complex.hpp"
#include "autodiff/dual.hpp"

#include <cmath>
#include <cstdio>

using C  = Complex<double>;
using D  = Dual<float>;
using DC = Dual<Complex<double>>;

// ---------------------------------------------------------------------------
// 1. 契约的编译期断言(不满足就编不过, 不用等到运行时)
// ---------------------------------------------------------------------------
static_assert(Scalar<float>,                      "Scalar<float>");
static_assert(Scalar<double>,                     "Scalar<double>");
static_assert(CommutativeScalar<float>,           "CommutativeScalar<float>");
static_assert(Scalar<C>,                          "Scalar<Complex<double>>");
static_assert(CommutativeScalar<C>,               "CommutativeScalar<Complex<double>>");
static_assert(Scalar<D>,                          "Scalar<Dual<float>>");
static_assert(CommutativeScalar<D>,               "CommutativeScalar<Dual<float>>");
// V1.5 的组合类型: 复系数 + 自动求导(波包演化)
static_assert(Scalar<DC>,                         "Scalar<Dual<Complex<double>>>");
static_assert(CommutativeScalar<DC>,              "CommutativeScalar<Dual<Complex<double>>>");
// ScalarTraits<T>::real_t 永远是实数: 一切比较只在实数域发生
static_assert(std::same_as<ScalarTraits<C>::real_t, double>,  "real_t of Complex<double>");
static_assert(std::same_as<ScalarTraits<D>::real_t, float>,   "real_t of Dual<float>");

// ---------------------------------------------------------------------------
// 测试辅助
// ---------------------------------------------------------------------------
static int failures = 0;

void check(const char* name, double got, double want, double eps = 1e-9)
{
    bool ok = std::fabs(got - want) <= eps;
    if (!ok) failures++;
    std::printf("  %-34s got=%15.9f  want=%15.9f  %s\n",
                name, got, want, ok ? "[OK]" : "[FAIL]");
}

void check_true(const char* name, bool ok)
{
    if (!ok) failures++;
    std::printf("  %-34s %s\n", name, ok ? "[OK]" : "[FAIL]");
}

// ---------------------------------------------------------------------------
// 2. 常数通道 fromReal(real_t)
// ---------------------------------------------------------------------------
void test_from_real()
{
    std::printf("\n[1] fromReal(real_t): 任意实数常数的唯一入口\n");

    check("float 恒等", ScalarTraits<float>::fromReal(0.5f), 0.5);
    check("double 恒等", ScalarTraits<double>::fromReal(0.5), 0.5);

    C c = ScalarTraits<C>::fromReal(0.5);            // 复数: 虚部折算为 0
    check("Complex 实部", c.real, 0.5);
    check("Complex 虚部 (=0)", c.imag, 0.0);

    D d = ScalarTraits<D>::fromReal(0.5f);           // Dual: 常数种子 derive = 0
    check("Dual value", d.value, 0.5);
    check("Dual derive (=0)", d.derive, 0.0);

    DC dc = ScalarTraits<DC>::fromReal(0.5);         // 组合类型: 逐层折算
    check("Dual<Complex> value.real", dc.value.real, 0.5);
    check("Dual<Complex> value.imag (=0)", dc.value.imag, 0.0);
    check("Dual<Complex> derive (=0)", dc.derive.real, 0.0);
}

// ---------------------------------------------------------------------------
// 3. conj / abs2 (conj 对 Dual 必须逐分量 —— 见 dual.hpp 末尾)
// ---------------------------------------------------------------------------
void test_conj_abs2()
{
    std::printf("\n[2] conj / abs2\n");

    check("Complex conj 实部", ScalarTraits<C>::conj(C(3, 4)).real, 3.0);
    check("Complex conj 虚部", ScalarTraits<C>::conj(C(3, 4)).imag, -4.0);
    check("Complex abs2", ScalarTraits<C>::abs2(C(3, 4)), 25.0);

    // Dual: conj 逐分量转发基域共轭 —— epsilon 分量也要 conj
    // (逐分量对、只 conj(epsilon) 错: d(z̄) = conj(dz))
    D x(3.0f, 4.0f);
    D xc = ScalarTraits<D>::conj(x);
    check("Dual conj value", xc.value, 3.0);
    check("Dual conj derive", xc.derive, 4.0);
    // abs2 只看 value(比较/选主元只在实数域发生)
    check("Dual abs2 (=value^2)", ScalarTraits<D>::abs2(x), 9.0);

    // Matrix 层: Adjoint 走 conj
    Matrix<D> A(2, 2);
    A(0, 0) = D(1.0f, 0.5f);
    A(1, 0) = D(7.0f, 8.0f);
    Matrix<D> AH = A.Adjoint();
    check("Adjoint 转置位置", AH(0, 1).value, 7.0);
    check("Adjoint conj derive", AH(0, 1).derive, 8.0);
}

// ---------------------------------------------------------------------------
// 4. 交换律单测 (CommutativeScalar 必须 a*b == b*a)
// ---------------------------------------------------------------------------
void test_commutative()
{
    std::printf("\n[3] 交换律 a*b == b*a (CommutativeScalar 的硬要求)\n");

    C a(1.5, 2.5), b(-0.5, 3.0);
    C ab = a * b, ba = b * a;
    check_true("Complex: a*b == b*a",
               ab.real == ba.real && ab.imag == ba.imag);

    D p(1.5f, 2.5f), q(-0.5f, 3.0f);
    D pq = p * q, qp = q * p;
    check_true("Dual: a*b == b*a",
               pq.value == qp.value && pq.derive == qp.derive);

    DC u(C(1.5, 2.5), C(0.5, -1.0)), v(C(-0.5, 3.0), C(1.0, 2.0));
    DC uv = u * v, vu = v * u;
    check_true("Dual<Complex>: a*b == b*a",
               uv.value.real == vu.value.real && uv.value.imag == vu.value.imag &&
               uv.derive.real == vu.derive.real && uv.derive.imag == vu.derive.imag);
}

// ---------------------------------------------------------------------------
// 5. Complex 原函数(依赖 fromReal 的路径: sqrt/log)
// ---------------------------------------------------------------------------
void test_complex_functions()
{
    std::printf("\n[4] Complex<double> 原函数对拍\n");

    C z = sqrt(C(3.0, 4.0));                       // 主值: 2 + i
    check("sqrt(3+4i) 实部", z.real, 2.0);
    check("sqrt(3+4i) 虚部", z.imag, 1.0);

    C l = log(C(3.0, 4.0));                        // ln|z| + i·arg
    check("log(3+4i) 实部 (=ln5)", l.real, 1.6094379124341003);
    check("log(3+4i) 虚部 (=atan2)", l.imag, 0.9272952180016122);

    C e = exp(C(1.0, 1.0));                        // e·(cos1 + i·sin1)
    check("exp(1+i) 实部", e.real, 1.4686939399158851);
    check("exp(1+i) 虚部", e.imag, 2.2873552871788423);

    C s = sin(C(1.0, 1.0));                        // sin x ch y + i cos x sh y
    check("sin(1+i) 实部", s.real, 1.2984575814159773);
    check("sin(1+i) 虚部", s.imag, 0.6349639147847361);

    C co = cos(C(1.0, 1.0));
    check("cos(1+i) 实部", co.real, 0.8337300251311491);
    check("cos(1+i) 虚部", co.imag, -0.9888977057628651);

    C pw = pow(C(1.0, 1.0), C(2.0, 0.0));          // (1+i)^2 = 2i
    check("pow(1+i, 2) 实部 (=0)", pw.real, 0.0, 1e-15);
    check("pow(1+i, 2) 虚部 (=2)", pw.imag, 2.0, 1e-15);
}

// ---------------------------------------------------------------------------
// 6. Dual 原函数(字面量全部走 traits 之后)
// ---------------------------------------------------------------------------
void test_dual_functions()
{
    std::printf("\n[5] Dual<float> 原函数对拍 (value 与 derive 分开比)\n");

    D e = exp(D(1.0f, 1.0f));                      // (e, e)
    check("exp(1) value", e.value, 2.718281828459045, 1e-6);
    check("exp(1) derive", e.derive, 2.718281828459045, 1e-6);

    D l = ln(D(2.0f, 1.0f));                       // (ln2, 1/2)
    check("ln(2) value", l.value, 0.6931471805599453, 1e-7);
    check("ln(2) derive", l.derive, 0.5, 1e-7);

    // sqrt 里用 fromReal(0.5): d√x = dx/(2√x)
    D r = sqrt(D(4.0f, 2.0f));                     // (2, 0.5)
    check("sqrt(4) value", r.value, 2.0, 1e-7);
    check("sqrt(4) derive", r.derive, 0.5, 1e-7);

    // pow(x, n): d(x^n) = n·x^(n-1)·dx (常数种子 y' = 0)
    D p = pow(D(2.0f, 1.0f), D(3.0f, 0.0f));       // (8, 12)
    check("pow(2,3) value", p.value, 8.0, 1e-6);
    check("pow(2,3) derive", p.derive, 12.0, 1e-6);

    D c = cos(D(0.0f, 1.0f));                      // (1, 0)
    check("cos(0) value", c.value, 1.0, 1e-7);
    check("cos(0) derive (=0)", c.derive, 0.0, 1e-7);

    D sn = sin(D(0.0f, 1.0f));                     // (0, 1)
    check("sin(0) value", sn.value, 0.0, 1e-7);
    check("sin(0) derive (=1)", sn.derive, 1.0, 1e-7);
}

// ---------------------------------------------------------------------------

int main()
{
    test_from_real();
    test_conj_abs2();
    test_commutative();
    test_complex_functions();
    test_dual_functions();

    std::printf("\n%s (%d failures)\n",
                failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED", failures);
    return failures ? 1 : 0;
}
