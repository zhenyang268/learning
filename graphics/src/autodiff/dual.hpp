#pragma once

#include <cmath>
#include <functional>
#include <vector>
#include "common.hpp"

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
        return Dual(x.value, 0);
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
        return Dual(cos(x.value), -1 * sin(x.value) * x.derive);
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

    // 对一行 Dual 累加: Y = Σ xi
    // forward 模式: dY/dxi = 1
    // 输出每个位置value都等于总和, 微分为1
    friend Dual sumLine(Dual* arr, int n)
    {
        Dual dual;
        for (int i = 0; i < n; i++) {
            dual.value  = dual.value + arr[i].value;
        }
        for (int i = 0; i < n; i++) {
            arr[i].value  = dual.value;
            arr[i].derive = ScalarTraits<T>::one();
        }
        return dual;
    }

    // 对一行 Dual 累乘: Y = Π xi
    // forward 模式: dY = Σ (total/xi)·dxi, 输出每个位置都等于总积及其全微分
    friend Dual multLine(Dual* arr, int n)
    {
        Dual dual;
        dual.value = T(1);
        for (int i = 0; i < n; i++) {
            dual.value = dual.value * arr[i].value;
        }

        T d = T(0);
        for (int i = 0; i < n; i++) {
            d += (dual.value / arr[i].value) * arr[i].derive;
        }
        for (int i = 0; i < n; i++) {
            arr[i].value  = dual.value;
            arr[i].derive = d;
        }
        return dual;
    }

    friend Dual dualFunc(Dual &x, FuncType type)
    {
        switch (type)
        {
        case FuncType::SIN:
            return sin(x);
        case FuncType::COS:
            return cos(x);
        case FuncType::EXP:
            return exp(x);
        case FuncType::LN:
            return ln(x);
        case FuncType::SQRT:
            return sqrt(x);
        default:
            return x;
        }
    }
};

template <typename T>
using DualFunc = std::function<T(T&)>;

// DualFunc<Dual<float>>