#pragma once

#include <cmath>
#include <functional>
#include <vector>
#include <map>
#include "core/common.hpp"
#include "core/scalar_traits.hpp"

template <typename T>
class Dual
{
public:
    T value;    // 函数值
    T derive;   // 微分值

    // 必须有默认初始化, 0.0f
    Dual(T v = T(), T d = T()) : value(v), derive(d) {}

    Dual operator+(const Dual& other) const
    {
        return Dual(value + other.value, derive + other.derive);
    }

    Dual operator-(const Dual& other) const
    {
        return Dual(value - other.value, derive - other.derive);
    }

    Dual operator-() const
    {
        return Dual(-value, -derive);
    }

    Dual operator*(const Dual& other) const
    {
        return Dual(value * other.value,
                value * other.derive + derive * other.value);
    }

    Dual operator/(const Dual& other) const
    {
        // 放弃除0判断, 让用户自己处理
        T v = value / other.value;
        T d = (derive * other.value - other.derive * value) / (other.value * other.value);
        return Dual(v, d);
    }

    Dual& operator=(Dual dual)
    {
        value = dual.value;
        derive = dual.derive;
        return *this;
    }

    // 数学函数
    // 常数用于函数矩阵拓展
    friend Dual constant(const Dual& x)
    {
        return Dual(x.value, ScalarTraits<T>::zero());
    }

    friend Dual sin(const Dual& x)
    {
        using std::sin;
        using std::cos;
        return Dual(sin(x.value), cos(x.value) * x.derive);
    }

    friend Dual cos(const Dual& x)
    {
        using std::sin;
        using std::cos;
        return Dual(cos(x.value), (-ScalarTraits<T>::one()) * sin(x.value) * x.derive);
    }

    friend Dual exp(const Dual& x)
    {
        using std::exp;
        T e = exp(x.value);
        return Dual(e, e * x.derive);
    }

    friend Dual ln(const Dual& x)
    {
        using std::log;
        return Dual(log(x.value), x.derive / x.value);
    }

    // x ^ y 也支持 x ^ n
    friend Dual pow(const Dual &x, const Dual &y)
    {
        using std::pow;
        using std::log;
        T v = pow(x.value, y.value);
        T d = v * (y.derive * log(x.value) + x.derive * y.value / x.value);
        return Dual(v, d);
    }

    friend Dual sqrt(const Dual& x)
    {
        using std::sqrt;
        T v = sqrt(x.value);
        return Dual(v, x.derive * ScalarTraits<T>::const_n(0.5) / v);
    }

    friend Dual sumLine(const std::vector<Dual>& arr)
    {
        Dual dual(ScalarTraits<T>::zero(), ScalarTraits<T>::zero());
        for (const Dual& x : arr) {
            dual = dual + x;
        }
        return dual;
    }

    // 对一行 Dual 累乘: Y = Π xi
    // forward: dY = Σ_i (Π_{j≠i} v_j)·dxi —— 前后缀积 O(n), 零因子安全
    // (修原实现 total/arr[i].value 在零因子处 0/0=NaN 污染整行的 bug, matrix.md §7.1 第 4 条)
    friend Dual multLine(const std::vector<Dual>& arr)
    {
        const int n = static_cast<int>(arr.size());
        if (n == 0) return Dual(ScalarTraits<T>::one(), ScalarTraits<T>::zero());

        std::vector<T> pre(n), suf(n);   // pre[i] = Π_{j<i} v_j, suf[i] = Π_{j>i} v_j
        pre[0] = ScalarTraits<T>::one();
        for (int i = 1; i < n; i++) pre[i] = pre[i - 1] * arr[i - 1].value;
        suf[n - 1] = ScalarTraits<T>::one();
        for (int i = n - 2; i >= 0; i--) suf[i] = suf[i + 1] * arr[i + 1].value;

        T total = pre[n - 1] * arr[n - 1].value;
        T d = ScalarTraits<T>::zero();
        for (int i = 0; i < n; i++) d = d + pre[i] * suf[i] * arr[i].derive;
        return Dual(total, d);
    }

    friend Dual dualFunc(const Dual &x, FunctionType type)
    {
        switch (type)
        {
        case FunctionType::Sin:
            return sin(x);
        case FunctionType::Cos:
            return cos(x);
        case FunctionType::Exp:
            return exp(x);
        case FunctionType::Log:
            return ln(x);
        case FunctionType::Sqrt:
            return sqrt(x);
        default:
            return x;
        }
    }
};