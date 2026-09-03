#pragma once

#include <cmath>
#include <functional>
#include <vector>

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
        if (other.value == 0) {
            return Dual(0, 0);
        }
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
        return Dual(std::sin(x.value), std::cos(x.value) * x.derive);
    }

    friend Dual cos(const Dual& x)
    {
        return Dual(std::cos(x.value), -1 * std::sin(x.value) * x.derive);
    }

    friend Dual exp(const Dual& x)
    {
        T e = std::exp(x.value);
        return Dual(e, e * x.derive);
    }

    friend Dual ln(const Dual& x)
    {
        return Dual(std::log(x.value), x.derive / x.value);
    }

    // x ^ y 也支持 x ^ n
    friend Dual pow(const Dual &x, const Dual &y)
    {
        T v = std::pow(x.value, y.value);
        T d = v * (y.derive * std::log(x.value) + x.derive * y.value / x.value);
        return Dual(v, d);
    }

    friend Dual sqrt(const Dual& x)
    {
        T v = std::sqrt(x.value);
        return Dual(v, x.derive / (2.0f * v));
    }

    friend Dual setDeriveOne(Dual& dual)
    {
        dual.derive = T(1);
        return dual;
    }

    friend Dual setDeriveZero(Dual& dual)
    {
        dual.derive = T(0);
        return dual;
    }

    // 对一行 Dual 累加: Y = Σ xi
    // forward 模式: dY = Σ dxi, 输出每个位置都等于总和及其全微分
    friend Dual sumLine(Dual* arr, int n)
    {
        Dual dual;
        for (int i = 0; i < n; i++) {
            dual.value  += arr[i].value;
            dual.derive += arr[i].derive;
        }
        for (int i = 0; i < n; i++) {
            arr[i].value  = dual.value;
            arr[i].derive = dual.derive;
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
            dual.value *= arr[i].value;
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
};

template <typename T>
using DualFunc = std::function<T(T&)>;

// DualFunc<Dual<float>>