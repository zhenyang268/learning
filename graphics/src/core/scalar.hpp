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
    { ScalarTraits<T>::const_n(std::declval<T>()) } -> std::same_as<T>;
};

template <>
struct ScalarTraits<float> {
    using real_t = float;
    static constexpr bool commutative = true;

    static real_t zero() { return 0.0f; }
    static real_t one() { return 1.0f; }

    static real_t conj(real_t x) { return x; }
    static real_t abs2(real_t x) { return x * x; }

    static real_t const_n(real_t x) { return x; }
};

template <>
struct ScalarTraits<double> {
    using real_t = double;
    static constexpr bool commutative = true;

    static real_t zero() { return 0.0; }
    static real_t one() { return 1.0; }

    static real_t conj(real_t x) { return x; }
    static real_t abs2(real_t x) { return x * x; }

    static real_t const_n(real_t x) { return x; }
};

template <typename T> requires Scalar<T>
struct ScalarTraits<Complex<T>> {
    using real_t = typename ScalarTraits<T>::real_t;
    static constexpr bool commutative = ScalarTraits<T>::commutative;

    static Complex<T> zero() { return Complex<T>(ScalarTraits<T>::zero(), ScalarTraits<T>::zero()); }
    static Complex<T> one() { return Complex<T>(ScalarTraits<T>::one(), ScalarTraits<T>::zero()); }

    static Complex<T> conj(const Complex<T>& x) { return Complex<T>(x.real, -x.imag); }
    static real_t abs2(const Complex<T>& x) { return x.real * x.real + x.imag * x.imag; }

    static Complex<T> const_n(const T& x) { return Complex<T>(x, ScalarTraits<T>::zero()); }
};

template <typename T> requires Scalar<T>
struct ScalarTraits<Dual<T>> {
    using real_t = typename ScalarTraits<T>::real_t;
    static constexpr bool commutative = ScalarTraits<T>::commutative;

    static Dual<T> zero() { return Dual<T>(ScalarTraits<T>::zero(), ScalarTraits<T>::zero()); }
    static Dual<T> one() { return Dual<T>(ScalarTraits<T>::one(), ScalarTraits<T>::zero()); }

    static Dual<T> conj(const Dual<T>& x) { return Dual<T>(ScalarTraits<T>::conj(x.real),
                                                        ScalarTraits<T>::conj(x.dual)); }
    static real_t abs2(const Dual<T>& x) { return ScalarTraits<T>::abs2(x.real); }

    static Dual<T> const_n(const T& x) { return Dual<T>(x, ScalarTraits<T>::zero()); }
};

// 细化 concept：在 requires 里查这个 flag。
// 没登记的类型 → 成员不存在 → 软失败判 false（报错干净，不是两百行天书）
template <typename T>
concept CommutativeScalar = Scalar<T> && requires {
    requires ScalarTraits<T>::commutative;      // 必须是编译期 true
};