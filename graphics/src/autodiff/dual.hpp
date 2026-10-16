#pragma once

#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>
#include <map>
#include "core/common.hpp"
#include "core/scalar_traits.hpp"

// 域判别: IsDual<Dual<U>> = true, 其余 (double / float / Complex<...>) false。
// 用途 —— `Apply` 的重载分裂 (见 jacobian.hpp): Dual 专用版与标量版共存时, 靠
// [temp.deduct.partial] 部分排序把 Dual 域分流到 Dual 版; 本 trait 是标量版体内的
// 守卫 `static_assert(!IsDual<T>::value, ...)` —— 一旦偏序选错, 立刻编不过而不是静默走错路。
template <typename T> class Dual;                       // 前置声明, 供下面的特化用
template <typename T> struct IsDual : std::false_type {};
template <typename U> struct IsDual<Dual<U>> : std::true_type {};

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

    // 拷贝赋值不手写 —— 原 operator= 与隐式版本逐字段等价, 却触发了 -Wdeprecated-copy:
    // 触发条件是"有 user-provided 拷贝赋值 + 没有 user-declared 拷贝构造",
    // 把参数从按值改成 const& 并不解除(实测)。交给编译器合成即可,
    // 拷贝构造/拷贝赋值/移动构造/移动赋值全都在, 且 `d = 0.5f` 仍经非 explicit
    // 构造转成临时 Dual 后走拷贝赋值。

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

    // ln 与 log 在 Dual 域绑成同一个函数 —— 这不是便利别名, 是**修 ADL 断链**:
    // 标量域与 Complex 域的函数名都是 log, 本类却只有 ln, 而上面 ln 的体内写的是
    // `log(x.value)` (依赖调用 ⇒ 名字查找 = 定义点普通查找 + 实例化点 ADL)。
    // 当 T = Dual<U> (嵌套对偶 = 二阶导) 时 `x.value` 的类型是 Dual<U>, ADL 只会去找
    // Dual<U> 的友元 —— 它只有 ln、没有 log ⇒ `error: no matching function for call to
    // 'log'` (2026-10-16 实测, dual.hpp:86)。pow 的体内同样用了 log, 一并受益。
    friend Dual log(const Dual& x) { return ln(x); }

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
        return Dual(v, x.derive * ScalarTraits<T>::fromReal(0.5) / v);
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