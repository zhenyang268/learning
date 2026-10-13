#pragma once

/**
 * 通过template实现的标量类
 * 通过traits屏蔽不同类型之间的运算
 * 通过concepts实现类型检查
 */

#include <utility>
#include <cstdint>
#include <concepts>
#include <type_traits>
#include <cmath>
#include <limits>
#include "core/complex.hpp"
#include "autodiff/dual.hpp"
#include "core/scalar_traits.hpp"


template <typename T>
concept Scalar = std::semiregular<T> && requires (T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
    { a * b } -> std::same_as<T>;
    { a / b } -> std::same_as<T>;
    { -a } -> std::same_as<T>;
    { ScalarTraits<T>::zero() } -> std::same_as<T>;
    { ScalarTraits<T>::one() } -> std::same_as<T>;
    { ScalarTraits<T>::conj(std::declval<T>()) } -> std::same_as<T>;
    { ScalarTraits<T>::abs2(std::declval<T>()) } -> std::same_as<typename ScalarTraits<T>::real_t>;
    // 任意实数常数的唯一入口(0/1 走 zero/one 快路); matrix.md 7.1 第 2 条。
    // 写进 concept 是为了堵住"泛型代码裸写字面量"这个最隐蔽的实数假设。
    { ScalarTraits<T>::fromReal(std::declval<typename ScalarTraits<T>::real_t>()) } -> std::same_as<T>;
    // 取实部: 一切数值判据(正定性 / 收敛 / 主元)只发生在实数域 —— 与 abs2 同口径,
    // Dual 只看 value 不带 epsilon 分量。matrix.md 7.1。
    { ScalarTraits<T>::real(std::declval<T>()) } -> std::same_as<typename ScalarTraits<T>::real_t>;
};

template <>
struct ScalarTraits<float> {
    using real_t = float;
    static constexpr bool commutative = true;

    static real_t zero() { return 0.0f; }
    static real_t one() { return 1.0f; }

    static real_t conj(real_t x) { return x; }
    static real_t abs2(real_t x) { return x * x; }
    static real_t real(real_t x) { return x; }        // 实数域: 实部就是自己

    static real_t fromReal(real_t x) { return x; }   // float/double: 恒等折算
};

template <>
struct ScalarTraits<double> {
    using real_t = double;
    static constexpr bool commutative = true;

    static real_t zero() { return 0.0; }
    static real_t one() { return 1.0; }

    static real_t conj(real_t x) { return x; }
    static real_t abs2(real_t x) { return x * x; }
    static real_t real(real_t x) { return x; }        // 实数域: 实部就是自己

    static real_t fromReal(real_t x) { return x; }   // float/double: 恒等折算
};

template <typename T> requires Scalar<T>
struct ScalarTraits<Complex<T>> {
    using real_t = typename ScalarTraits<T>::real_t;
    static constexpr bool commutative = ScalarTraits<T>::commutative;

    static Complex<T> zero() { return Complex<T>(ScalarTraits<T>::zero(), ScalarTraits<T>::zero()); }
    static Complex<T> one() { return Complex<T>(ScalarTraits<T>::one(), ScalarTraits<T>::zero()); }

    static Complex<T> conj(const Complex<T>& x) { return Complex<T>(x.real, -x.imag); }
    static real_t abs2(const Complex<T>& x) { return x.real * x.real + x.imag * x.imag; }
    // 取实部: 递归下钻, 兼容 Complex<Dual<...>> 这类嵌套 (real_t 由基域决定)
    static real_t real(const Complex<T>& x) { return ScalarTraits<T>::real(x.real); }

    // 复数常量: 虚部折算为 0 (D15 的 traits 化)
    static Complex<T> fromReal(real_t x)
    {
        return Complex<T>(ScalarTraits<T>::fromReal(x), ScalarTraits<T>::zero());
    }
};

template <typename T> requires Scalar<T>
struct ScalarTraits<Dual<T>> {
    using real_t = typename ScalarTraits<T>::real_t;
    static constexpr bool commutative = ScalarTraits<T>::commutative;

    static Dual<T> zero() { return Dual<T>(ScalarTraits<T>::zero(), ScalarTraits<T>::zero()); }
    static Dual<T> one() { return Dual<T>(ScalarTraits<T>::one(), ScalarTraits<T>::zero()); }

    // conj 逐分量转发基域共轭: conj(Dual) = Dual(conj(value), conj(derive))
    // (epsilon 共轭是坑 —— 逐分量对、只 conj(epsilon) 错, matrix.md 7.1 第 6 条)
    static Dual<T> conj(const Dual<T>& x) { return Dual<T>(ScalarTraits<T>::conj(x.value),
                                                        ScalarTraits<T>::conj(x.derive)); }
    // abs2 只看 value: 比较/选主元只在实数域发生, 不带 epsilon 分量
    static real_t abs2(const Dual<T>& x) { return ScalarTraits<T>::abs2(x.value); }
    // 取实部: 与 abs2 同口径 —— 判据是纯数值决策, 不带 epsilon 分量
    static real_t real(const Dual<T>& x) { return ScalarTraits<T>::real(x.value); }

    // 常数种子: value = 该实数, derive = 0
    static Dual<T> fromReal(real_t x)
    {
        return Dual<T>(ScalarTraits<T>::fromReal(x), ScalarTraits<T>::zero());
    }
};

// 细化 concept：在 requires 里查这个 flag。
// 没登记的类型 → 成员不存在 → 软失败判 false（报错干净，不是两百行天书）
template <typename T>
concept CommutativeScalar = Scalar<T> && requires {
    requires ScalarTraits<T>::commutative;      // 必须是编译期 true
};

// 容差: eps 的倍数, 且随矩阵规模缩放 —— D12(主元 |p| <= tol 显式抛) / §10.8 ①(不写魔数 1e-9)。
// 不能用全局 Tol: 那是 float 的 epsilon(1.19e-7), 对 double(2.2e-16) 差 9 个数量级,
// 会按错误尺度判奇异性(§2.4: 容差不随离散尺度缩放 -> 把良定矩阵误判成奇异)。
// 放标量层(§1.4 数值细则): solve / inverse / cholesky / eigen / svd 共用,
// 避免 linalg 各文件为了取一个容差而互相 include。
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t DefaultTol(int n)
{
    using real_t = typename ScalarTraits<T>::real_t;
    return real_t(8) * real_t(n) * std::numeric_limits<real_t>::epsilon();
}