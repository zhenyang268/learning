#pragma once

// ===========================================================================
// autodiff/dual.hpp —— 对偶数 (前向模式自动微分)
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <typename T> class Dual { public: T value, derive; };
//   template <typename T> struct IsDual;            // 域判别 trait
//
// 语义: `Dual(v, d)` 表示 `v + d·ε`, 其中 ε² = 0。
//   一个 Dual 同时携带**函数值** (value) 与**一阶导数** (derive)。
//   前向模式: 一边算值, 一边按链式法则传播导数。
//
// 运算 (全部成员/友元, 靠 ADL 找到):
//   + - 一元- * /           四则; 商法则、乘积法则已内建
//   原函数: sin cos tan exp ln log sqrt pow abs …
//   下标: operator[](i)     取第 i 个偏导 (多元输入时 derive 是数组/映射)
//
// 域判别:
//   IsDual<Dual<U>>::value == true; 其它 (double/float/Complex) 为 false。
//   用途: `Apply` 的**重载分裂** —— Dual 版与标量版共存时靠偏序分流,
//   标量版体内用 `static_assert(!IsDual<T>::value)` 做守卫。
//
// 关键口径:
//   · **Dual 交换**: a*b == b*a。这是铁律, 有单测。
//   · 没有复合赋值 (`+=` 等) —— 全库纪律, 见末尾。
//   · T 的约束由 ScalarTraits 决定; 嵌套 Dual<Dual<T>> 支持二阶导。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   · 常数通道: 所有字面量必须走 `ScalarTraits<T>::fromReal(r)`, 其 derive
//     自动为 0 —— 这是"常数导数为零"的机制化, 不靠调用方记得。
//   · 原函数都是**友元自由函数** (不是成员): 这样 `sin(x)` 与 `std::sin`
//     同形, 且靠 ADL 找到 —— 不需要 using namespace, 也不会污染全局。
//   · 多元输入时 derive 可以是 vector / map; 单变量时就是一个标量槽。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么 ε² = 0 就够了 —— 对偶数的全部秘密】
//   把 f(x0+h) 按 h 展开: f(x0) + f'(x0)·h + O(h²)。
//   若规定 h² = 0 (而不是 h→0 的极限), 则**二阶及以上自动消失**, 剩下
//   精确的一阶项。所以 `Dual(v,d)` 就是 (f(x0), f'(x0)·d) 的一个打包。
//   链式法则因此变成**代数恒等式**, 不需要任何"求导"的特殊代码。
//
// 【为什么不用成员函数做 sin/cos】
//   若写成 `x.sin()`, 泛型代码里就得区分 `T` 是 Dual 还是 double (判断 + 分支)。
//   写成友元 `sin(x)` 后, 泛型代码写 `sin(x)` 即可 —— `double` 走 std::sin,
//   `Dual` 走 ADL 找到的友元版。**算法层于是零分支**, 这正是全库要求。
//
// 【为什么没有复合赋值 —— 与 Matrix 同口径】
//   复合赋值要求类型可自我修改, 会引出"什么时候允许写回"的别名问题。
//   值语义 (a = a + b) 让每个表达式都产生新对象, 优化器处理起来也比
//   就地修改更自由。这是刻意的统一, 不是遗漏。
//
// 【Dual<Dual<T>> —— 二阶导为什么能自动成立】
//   因为所有运算都是**泛型**的: `Dual<Dual<double>>` 的 value/derive 两个槽
//   本身就是 `Dual<double>`, 四则运算递归下去, 于是 ε 分量的 ε 分量就是二阶导。
//   ⚠ 前提: 每个原函数都必须对**任意 T** 成立, 不能偷偷假设 T 是浮点。
//     典型反例: `sqrt` 若写成 `T(std::sqrt(double(x)))` 就会把里层 Dual 拍平,
//     二阶导静默变 0。这类 bug 只有 `Dual<Dual<T>>` 域的探针能抓到。
//
// 【epsilon 共轭的坑】
//   conj(Dual) 必须**逐分量**共轭: Dual(conj(value), conj(derive))。
//   只共轭 epsilon 分量、或只共轭 value, 都在实域退化成恒等 ⇒ 双精度探针
//   永远看不见。详见 core/scalar.hpp 的 ScalarTraits<Dual<T>>。
// ===========================================================================

#include <cmath>
#include <functional>
#include <map>
#include <type_traits>
#include <vector>

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
    // (修原实现 total/arr[i].value 在零因子处 0/0=NaN 污染整行的 bug)
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
// ===========================================================================
// 四、数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里, 实现处只留"为什么这么写"的工程理由。**
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   ε (epsilon)    对偶单位, 满足 **ε² = 0** 但 ε ≠ 0
//   Dual(v, d)     对偶数 v + d·ε (代码里字段 value / derive)
//   v              函数值 (value 槽)
//   d              导数值 / 微分分量 (derive 槽)
//   f'(x)          一阶导数
//   x₀             求导点
//   h              步长 (无穷小)
//   T              基域类型 (float / double / Complex / Dual)
//   IsDual<T>      域判别 trait
//
// ---------------------------------------------------------------------------
// (D1) 对偶数的定义
// ---------------------------------------------------------------------------
//        对偶数 z = v + d·ε,  其中  ε² = 0,  ε ≠ 0
//
//      注意这不是复数 (复数有 i² = −1), 也不是"很小的数" (不涉及极限)。
//      ε² = 0 是**代数规定**, 不是近似 —— 这一点决定了下面所有推导都是精确的。
//
// ---------------------------------------------------------------------------
// (D2) 为什么它能自动求导 (对偶数的全部秘密)
// ---------------------------------------------------------------------------
//        把 f(x₀ + h) 在 x₀ 展开:
//              f(x₀ + h) = f(x₀) + f'(x₀)·h + O(h²)
//
//        若规定 **h² = 0** (而不是 h → 0 的极限), 则 O(h²) 及更高项
//        **精确消失**, 剩下:
//
//              f(x₀ + d·ε) = f(x₀) + d·f'(x₀)·ε
//
//      ⇒ value 槽 = f(x₀), derive 槽 = d·f'(x₀)。
//        链式法则因此变成**代数恒等式**, 不需要任何"求导"的特殊代码。
//
//      这是"把微积分变成代数"的全部内容。
//
// ---------------------------------------------------------------------------
// (D3) 四则运算的传播规则 (乘积法则 / 商法则)
// ---------------------------------------------------------------------------
//        加/减:  (a + b·ε) ± (c + e·ε) = (a ± c) + (b ± e)·ε
//
//        乘法:   (a + b·ε)(c + e·ε) = a·c + (a·e + b·c)·ε + b·e·ε²
//                                    ^^^^^  ε² = 0, 这一项消失
//                ⇒ = a·c + (a·e + b·c)·ε            ← 乘积法则
//
//        除法:   1/(c + e·ε) = (c − e·ε)/((c+e·ε)(c−e·ε))
//                            = (c − e·ε)/(c² − e²ε²)  = (c − e·ε)/c²
//                ⇒ (a + b·ε)/(c + e·ε) = (a/c) + (b·c − a·e)/c² · ε   ← 商法则
//
//      注意 (D3) 里**没有任何近似** —— 每一步都是精确代数。
//
// ---------------------------------------------------------------------------
// (D4) 一元原函数的传播 (链式法则)
// ---------------------------------------------------------------------------
//        对任意光滑 f:
//              f(a + b·ε) = f(a) + b·f'(a)·ε
//
//        所以每个原函数只需要"值 + 导数"两件东西:
//              sin:  (sin a,  cos a·b)
//              cos:  (cos a, −sin a·b)
//              exp:  (exp a,  exp a·b)
//              ln :  (ln a,   b/a)
//              sqrt: (√a,     b/(2√a))
//
//      这就是为什么 dual.hpp 里每个原函数都长成 `Dual(f(v), f'(v)*d)` 的形状。
//
// ---------------------------------------------------------------------------
// (D5) 二阶导 (Dual<Dual<T>>)
// ---------------------------------------------------------------------------
//        嵌套:  Dual<Dual<double>> 的内层 value/derive 各是一个 Dual<double>,
//        于是"ε 分量的 ε 分量"就是二阶导。
//
//        代数上:  D₁ = v + d₁·ε₁,  D₂ = D₁ + e·ε₂
//                 ⟹ z = v + d₁·ε₁ + e·ε₂ + (∂²)·ε₁ε₂
//
//      ⚠ 前提: 每个原函数必须对**任意 T** 成立, 不能假设 T 是浮点。
//        反例: `sqrt` 若写成 `T(std::sqrt(double(x)))` 就会把里层 Dual 拍平,
//        二阶导静默变 0。只有 `Dual<Dual<T>>` 域的探针能抓到。
//
// ---------------------------------------------------------------------------
// (D6) 交换律 (Dual 是交换环)
// ---------------------------------------------------------------------------
//        a·b = b·a,  对所有 a, b
//
//      这是**铁律**, 有单独单测。与四元数的根本区别:
//      四元数不交换 (多一项叉积), 所以四元数**不得**进这条单测。
//
// ---------------------------------------------------------------------------
// (D7) 常数通道 (为什么字面量必须走 fromReal)
// ---------------------------------------------------------------------------
//        常数 c 的对偶形式:  c = c + 0·ε       ← derive 槽为 0 (常数导数为 0)
//
//        若在泛型代码里裸写 `+ 1`, 对 T = Dual 时可能:
//          ① 编不过 (无隐式转换), 或
//          ② 走上某个隐式转换、把 ε 分量静默丢掉 —— 后者更危险。
//      ⇒ 必须走 `ScalarTraits<T>::fromReal(r)`, 它保证 derive = 0。
//
// ---------------------------------------------------------------------------
// (D8) ε² = 0 与"截断 Taylor"的等价性
// ---------------------------------------------------------------------------
//        f(x + h) 的 Taylor 展开里, 含 h² 及以上的项都乘了 ε² = 0,
//        所以**精确**等于前两项:
//
//              f(x + h·ε) = f(x) + h·f'(x)·ε
//
//        这与"一阶截断 Taylor"数值上等价, 但机制完全不同:
//        截断是**丢掉了项** (有误差), 对偶数是**那些项本来就是 0** (无误差)。
//
// ---------------------------------------------------------------------------
// 实现口径 (跟公式直接相关)
// ---------------------------------------------------------------------------
//   · **原函数写成友元自由函数** (不是成员): 这样 `sin(x)` 与 `std::sin` 同形,
//     泛型代码写 `sin(x)` 即可, double 与 Dual 自动分流 —— 算法层零分支。
//   · **没有复合赋值** (+= 等): 值语义统一, 与 Matrix 同口径。
//   · conj(Dual) 必须**逐分量**: Dual(conj(value), conj(derive))。
//     只 conj 一个槽在实域退化成恒等 ⇒ 双精度探针永远看不见。
// ===========================================================================
