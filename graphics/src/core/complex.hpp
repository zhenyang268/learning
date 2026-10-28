#pragma once

// ===========================================================================
// core/complex.hpp —— 复数 Complex<T>
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
// 表示: 两个实槽 real / imag。T 是实浮点 (float / double / long double)。
//
// 构造: Complex() = (0,0); Complex(re, im)。拷贝/移动全部 = default。
//
// 运算 (实参全在类内):
//   +  -  一元-   *   /     复数四则; **除法用 Smith 算法** (防溢出/下溢)
//   成员 real / imag 公开可读写 (不封装 —— 哑数据, 同 Matrix)
//
// 原函数 (友元, 靠 ADL 找到):
//   exp(z)    = e^a (cos b + i sin b)                       欧拉公式
//   log(z)    = ln|z| + i·arg(z)                            主值, 分支切割在负实轴
//   sqrt(z)   = |z| 半角公式, Re >= 0                       主值, 分支切割同上
//   pow(z, w) = exp(w·log z)                                主值
//
// **没有**复合赋值 (+= 等) —— 全库纪律: 容器与标量类型尽量只给值语义,
//   复合赋值留给调用方写成 a = a + b。`Dual` 同此。
// **没有**全序: 不提供 < <= > >=。要比较只能逐槽比, 或比 abs2
//   (见 numeric-contract.md D7)。
//
// 域约束:
//   · **不用 std::complex**。理由见末尾。
//   · `Complex<T>` 要进 linalg, 必须由 `core/scalar.hpp` 特化
//     `ScalarTraits<Complex<T>>` (zero/one/conj/abs2/real/fromReal)。
//
// ---------------------------------------------------------------------------
// 二、内部机制
// ---------------------------------------------------------------------------
//   · 除法走 Smith: 先比 |c| 与 |d|, 按大的那个归一化, 避免 c²+d² 溢出
//     或下溢到 0。朴素公式 (a+bi)/(c+di) = (ac+bd + (bc-ad)i)/(c²+d²)
//     在 |c| 或 |d| 接近 1e±150 时会出 inf/inf 或 0/0。
//   · sqrt 用 hypot 而非 sqrt(x²+y²) —— 后者在 |x|,|y| ~ 1e200 时平方溢出。
//     半角公式里 (r±x) 也用「不抵消」的写法 (见实现处), 避免 x 接近 -r 时相消。
//   · exp / log / pow / sqrt 全是**主值**, 分支切割在负实轴: 两侧分别给 +i/-i,
//     靠 copysign 取符号。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么自己写而不 std::complex】
//   ① `std::complex<T>` 的 `operator<` 在 C++20 起被移除但历史版本语义混乱,
//      而本项目要求「算法层零实复分支」—— 需要一个**行为完全可控**的 T。
//   ② 要能特化 ScalarTraits (zero/one/fromReal), 并在 Dual 域做嵌套
//      (Dual<Complex>)。std::complex 的成员函数集是封闭的, 加不了这些。
//   ③ 可测: 自己的实现能逐条对拍解析式; 库实现的黑盒行为在不同编译器上
//      可能有差异 (尤其 pow 的分支切割)。
//
// 【除法为什么必须 Smith】
//   朴素算法在 |c|,|d| ~ 1e200 时 c²+d² 直接溢成 inf ⇒ 结果 0;
//   在 |c|,|d| ~ 1e-200 时 c²+d² 下溢成 0 ⇒ 结果 inf 或 nan。
//   而正确的商其实完全在正常范围内 (例如 (1e200+i)/(1e200+i) = 1)。
//   Smith 算法先除以模大的那个分量, 全程量级受控。
//   判据: 用 1e200 量级的输入, 结果必须有限且正确 —— 朴素版会给 inf/nan。
//
// 【exp 的实现口径】
//   exp(a+bi) = e^a (cos b + i sin b) —— 欧拉公式, 直接算, 不做级数。
//   ⚠ 这里不是「复数指数映射」(那是 expm 对矩阵做的), 别混。
// ===========================================================================

#include <cmath>

#include "core/scalar_traits.hpp"

template <typename T>
class Complex
{
public:
    T real;
    T imag;

    Complex() = default;
    Complex(T real, T imag) : real(real), imag(imag) {}

    Complex(const Complex&) = default;
    Complex& operator=(const Complex&) = default;

    Complex(Complex&&) = default;
    Complex& operator=(Complex&&) = default;

    ~Complex() = default;

    Complex operator+(const Complex& other) const
    {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const
    {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator-() const
    {
        return Complex(-real, -imag);
    }

    Complex operator*(const Complex& other) const
    {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    // Smith (1962) 算法：不做 c^2 + d^2，按 |c| 与 |d| 分流，避免上溢/下溢。
    // 代数上与原式恒等；分流条件保证 |r| <= 1，中间量不超过 ~2*max(|c|,|d|)。
    // 除零不判、不抛，让 IEEE 的 inf/NaN 沿运算传播（D12）。
    Complex operator/(const Complex& other) const
    {
        using std::abs;                 // 两步法：T 是 float/double 时走 std::abs
        const T c = other.real;
        const T d = other.imag;

        if (abs(c) >= abs(d)) {
            const T r = d / c;          // |r| <= 1
            const T denom = c + d * r;  // = (c^2 + d^2) / c
            return Complex((real + imag * r) / denom,
                           (imag - real * r) / denom);
        } else {
            const T r = c / d;          // |r| < 1
            const T denom = d + c * r;  // = (c^2 + d^2) / d
            return Complex((real * r + imag) / denom,
                           (imag * r - real) / denom);
        }
    }

    // complex functions
    friend Complex sin(const Complex& z)
    {
        using std::sin;
        using std::cos;
        using std::sinh;
        using std::cosh;
        return Complex(sin(z.real) * cosh(z.imag), cos(z.real) * sinh(z.imag));
    }

    friend Complex cos(const Complex& z)
    {
        using std::sin;
        using std::cos;
        using std::sinh;
        using std::cosh;
        return Complex(cos(z.real) * cosh(z.imag), -sin(z.real) * sinh(z.imag));
    }

    friend Complex tan(const Complex& z)
    {
        return sin(z) / cos(z);
    }

    friend Complex exp(const Complex& z)
    {
        using std::exp;
        return Complex(exp(z.real) * cos(z.imag), exp(z.real) * sin(z.imag));
    }

    friend Complex log(const Complex& z)
    {
        using std::log;
        using std::atan2;
        return Complex(ScalarTraits<T>::fromReal(0.5) * log(z.real * z.real + z.imag * z.imag), atan2(z.imag, z.real));
    }

    friend Complex pow(const Complex& z, const Complex& w)
    {
        return exp(w * log(z));
    }

    // sqrt 主值: Re >= 0; 分支切割在负实轴（两侧 -> ±i·sqrt|x|, 待核标注同 ln）
    // 两支分工避开抵消（Smith 精神）; r 用 hypot 防溢出。
    friend Complex sqrt(const Complex& z)
    {
        using std::sqrt;  using std::hypot;  using std::copysign;
        const T x = z.real;
        const T y = z.imag;
        const T zero = ScalarTraits<T>::zero();

        if (x == zero && y == zero)                       // z = 0: 否则 b = 0/0 = NaN
            return Complex(zero, zero);

        const T r = hypot(x, y);                          // |z|, 不走 sqrt(x^2+y^2)

        if (x >= zero) {
            const T a = sqrt((r + x) * ScalarTraits<T>::fromReal(0.5));  // r+x 不抵消
            return Complex(a, y / (a + a));               // b = y/(2a), 避开 r-x; a+a 避开字面量 2
        } else {
            const T b = copysign(sqrt((r - x) * ScalarTraits<T>::fromReal(0.5)), y);  // r-x 不抵消
            return Complex(y / (b + b), b);               // a = y/(2b), 避开 r+x
        }
    }

};
// ===========================================================================
// 数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里, 实现处只留"为什么这么写"的工程理由。**
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   z = a + bi     复数 (代码里 real = a, imag = b)
//   a              实部 (real 槽)
//   b              虚部 (imag 槽)
//   i              虚数单位, i² = −1
//   conj(z)        共轭 = a − bi
//   |z|            模 = √(a² + b²)
//   arg(z)         辐角 = atan2(b, a)
//   T              基域类型 (float / double / long double)
//   e              自然常数 2.718...
//   φ              相位角
//
// ---------------------------------------------------------------------------
// (Z1) 基本运算
// ---------------------------------------------------------------------------
//        加:   (a + bi) + (c + dii) = (a+c) + (b+d)i
//        乘:   (a + bi)(c + di) = (ac − bd) + (ad + bc)i
//        共轭: conj(a + bi) = a − bi
//        模方: |z|² = a² + b²   (代码里 abs2)
//
// ---------------------------------------------------------------------------
// (Z2) 除法 (Smith 算法)
// ---------------------------------------------------------------------------
//        朴素算法:
//              (a + bi)/(c + di) = [ (ac + bd) + (bc − ad)i ] / (c² + d²)
//
//        ⚠ 朴素算法在 |c|,|d| 极大或极小时**溢出/下溢**:
//          · |c|,|d| ~ 1e200 时 c²+d² 直接溢成 inf ⇒ 结果 0
//          · |c|,|d| ~ 1e-200 时 c²+d² 下溢成 0 ⇒ 结果 inf/nan
//          而正确的商其实完全在正常范围内 —— 例如 (1e200+i)/(1e200+i) = 1。
//
//        Smith 算法: 先比 |c| 与 |d|, 按**大的那个**归一化
//              if |c| >= |d|:
//                  r = d/c
//                  denom = c + d·r                    ← 不溢出
//                  商 = (a + b·r)/denom + (b − a·r)/denom · i
//              else:   (对称地按 d 归一化)
//
//        全程量级受控, 不会溢出。
//        判据: 用 1e200 量级的输入, 结果必须有限且正确 —— 朴素版会给 inf/nan。
//
// ---------------------------------------------------------------------------
// (Z3) 指数 (欧拉公式)
// ---------------------------------------------------------------------------
//              exp(a + bi) = e^a · (cos b + i sin b)
//
//        ⚠ 这**不是**矩阵指数 expm(A) —— 那是另一个东西 (linalg/expm.hpp)。
//          两者只在"1×1 矩阵"时重合: exp([a]) = [e^a]。
//
// ---------------------------------------------------------------------------
// (Z4) 对数 (主值)
// ---------------------------------------------------------------------------
//              log(z) = ln|z| + i·arg(z),      arg(z) ∈ (−π, π]
//
//        **分支切割在负实轴**: 负实轴两侧分别给 +π 与 −π。
//        这意味着 log 在负实轴上**不连续** —— 这是复变函数多值性的必然结果,
//        不是实现缺陷。需要连续版本的场合 (如复幂) 必须自己做相位展开。
//
// ---------------------------------------------------------------------------
// (Z5) 开方 (主值, 半角公式)
// ---------------------------------------------------------------------------
//              √z = √((|z| + a)/2) + sign(b)·√((|z| − a)/2) · i
//
//        · 实部取正号 ⇒ 结果在右半平面 (Re ≥ 0), 这是主值的定义
//        · 符号取 b 的符号 ⇒ 虚部与原复数同侧 (负实轴上按 copysign 约定, 见实现)
//        ⚠ 主值也**不连续**: 跨负实轴时虚部符号翻转。
//
//        实现里用 **hypot(a,b)** 而不是 sqrt(a²+b²): 后者在 |a|,|b| ~ 1e200 时
//        平方溢出。而且 (|z| ± a) 用"不抵消"的写法, 避免 a 接近 −|z| 时相消。
//
// ---------------------------------------------------------------------------
// (Z6) 幂 (主值)
// ---------------------------------------------------------------------------
//              z^w = exp(w · log z)
//
//        因为 log 是多值的, z^w 一般也是多值的; 取主值即取 log 的主值。
//        ⚠ 这与"实数幂展开"不等价 —— 例如 (−8)^(1/3) 的主值**不是** −2,
//          而是 1 + 1.732i (因为 log(−8) 的主值辐角是 +π)。
//          要实数立方根请自己判符号, 不要用 pow。
//
// ---------------------------------------------------------------------------
// (Z7) 为什么不用 std::complex (三条理由)
// ---------------------------------------------------------------------------
//        ① `std::complex<T>` 的 `operator<` 语义在不同标准版本里变过,
//           而本项目要求「算法层零实复分支」—— 需要一个**行为完全可控**的 T。
//        ② 要能特化 ScalarTraits (zero/one/conj/abs2/real/fromReal),
//           并在 Dual 域做嵌套 (Dual<Complex>)。std::complex 是封闭的, 加不了。
//        ③ 可测: 自己的实现能逐条对拍解析式; 库实现的黑盒行为
//           在不同编译器上可能有差异 (尤其 pow 的分支切割)。
//
// ---------------------------------------------------------------------------
// (Z8) 为什么没有全序
// ---------------------------------------------------------------------------
//        复数域**不是有序域** —— 不存在与加法/乘法相容的全序。
//        所以不提供 < <= > >=。
//
//        要比较只能:
//        · 比模长  |z₁| vs |z₂|        (常用)
//        · 逐分量比 a 与 b             (判相等)
//        ⇒ 这就是为什么全库的数值判据一律走 `abs2` (D7)。
// ===========================================================================
