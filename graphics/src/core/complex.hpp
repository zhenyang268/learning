#pragma once

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
        return Complex(ScalarTraits<T>::const_n(T(0.5)) * log(z.real * z.real + z.imag * z.imag), atan2(z.imag, z.real));
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