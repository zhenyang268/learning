#pragma once

// ===========================================================================
// core/scalar.hpp —— 标量约束与 traits
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
// 两个 concept + 一个工厂, 全部是**编译期**的东西, 没有运行期开销。
//
//   Scalar<T>             能做 + - * / 一元-, 且 traits 提供 zero/one/conj/
//                         abs2/fromReal/real。**不要求**乘法交换。
//   CommutativeScalar<T>  在 Scalar 基础上再要求 traits::commutative == true。
//                         linalg 的一切算法入口用它 (见末尾「为什么」)。
//   DefaultTol<T>(n)      容差 = 8·n·eps(real_t)。n 传矩阵规模。
//
// **允许的 T**: float / double / Complex<任意 Scalar> / Dual<任意 Scalar>。
//   嵌套也成立: Dual<Complex<double>>。
//
// 要接入新类型 T, 你只需特化 `ScalarTraits<T>` 六个静态函数, **不用改算法**:
//   zero() / one() / conj(x) / abs2(x) / real(x) / fromReal(r)
//   外加一个 `static constexpr bool commutative`。
//
// 错误用法: 在泛型代码里裸写 `0` / `1` / `1.0`。必须走
//   ScalarTraits<T>::zero()  (或 fromReal(0))
// 裸字面量在 T = Complex / Dual 时**编不过**; 更糟的是在 T = Dual 时可能
//   走上隐式转换、把 epsilon 分量静默丢掉。
//
// ---------------------------------------------------------------------------
// 二、内部机制
// ---------------------------------------------------------------------------
//   · `Scalar` concept 里 `{ ScalarTraits<T>::fromReal(...) } -> same_as<T>`
//     这一条是**故意的**: 把 fromReal 摆进 concept 而不是只写在文档里,
//     是为了让「泛型代码裸写字面量」在**概念约束处**就失败, 而不是等到某个
//     深处的模板实例化 —— 后者的报错是两百行天书。
//   · ScalarTraits 的偏序: Complex<T> / Dual<T> 的 `requires Scalar<T>` 保证了
//     递归下钻是良基的 (real_t 最终落到 float/double)。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么需要 commutative 这个 flag】
//   linalg 的算法 (LU / QR / Jacobi …) 都隐含假设 a·b = b·a —— 例如行交换
//   与消元的顺序可以互换、点积可以随意重排。四元数满足除环但不交换, 一旦
//   传进来, 结果是「编得过、算得错」。所以用 concept 在**入口**挡掉, 而不是
//   在算法里到处加 assert。
//   ⇒ 这正是 D16「四元数不做矩阵元素」的机制化落实。
//
// 【conj 对 Dual 为什么必须逐分量】
//   conj(Dual) = Dual(conj(value), conj(derive))。
//   看似「对」的写法是只共轭 epsilon 分量、或只共轭 value —— 两者都在
//   实数域退化成恒等, 于是**双精度探针永远看不见**。只有 Dual<Complex>
//   才暴露。这是一类典型事故: **实域测试免疫、复域测试才现形**。
//
// 【abs2 / real 为什么只看 value 槽】
//   一切数值判据 (选主元、判收敛、比正定性) 只在**实数域**发生, 不带
//   epsilon 分量 —— epsilon 是"导数方向的标记", 拿它参与大小比较没有意义。
//   见 numeric-contract.md D7。
//
// 【容差为什么放这里而不放 core/common.hpp】
//   common.hpp 里那个 `Tol = EPSILON` 是 **float** 的 eps (1.19e-7),
//   对 double (2.2e-16) 差 9 个数量级 —— 拿它判奇异性会把良定矩阵误判成奇异。
//   而 DefaultTol<T> 按 T 的 real_t 取 eps, 才没这个问题。
//   放标量层的另一个理由: solve / inverse / cholesky / eigen / svd 都要用它,
//   放 common 会让 linalg 各文件为了取一个容差而互相 include。
//   ⇒ **新代码一律用 DefaultTol<T>(n), 不要碰全局 Tol。**
// ===========================================================================

#include <cmath>
#include <concepts>
#include <limits>
#include <utility>

#include "autodiff/dual.hpp"
#include "core/complex.hpp"
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
    // 任意实数常数的唯一入口 (0/1 走 zero/one 快路)。
    // 写进 concept 是为了堵住"泛型代码裸写字面量"这个最隐蔽的实数假设。
    { ScalarTraits<T>::fromReal(std::declval<typename ScalarTraits<T>::real_t>()) } -> std::same_as<T>;
    // 取实部: 一切数值判据(正定性 / 收敛 / 主元)只发生在实数域 —— 与 abs2 同口径,
    // Dual 只看 value 不带 epsilon 分量。
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
    // (epsilon 共轭是坑 —— 逐分量对、只 conj(epsilon) 错; 见文件头「三」)
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

// 容差 = 8·n·eps(real_t)。为什么不用全局 Tol、为什么放标量层
// —— 见文件头「三、容差为什么放这里」。
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t DefaultTol(int n)
{
    using real_t = typename ScalarTraits<T>::real_t;
    return real_t(8) * real_t(n) * std::numeric_limits<real_t>::epsilon();
}
// ===========================================================================
// 四、数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里。**
// 本文件是**约束层** (concept + traits), 没有算法, 但它定义的六个函数
// 各自对应一个数学运算 —— 这一节说明"它们数学上是什么"。
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   T              标量域类型 (float / double / Complex / Dual)
//   real_t         T 对应的**实数**类型 (由 ScalarTraits<T>::real_t 给出)
//   x, y, z        标量
//   conj(x)        x 的共轭
//   |x|            x 的模 / 绝对值
//   ‖·‖           范数
//   ε             机器 epsilon
//   n             矩阵规模 (给容差用)
//   commutative    交换律标记
//
// ---------------------------------------------------------------------------
// (S1) traits 六个函数各自对应的数学运算
// ---------------------------------------------------------------------------
//        zero()         加法单位元 0        (满足 x + 0 = x)
//        one()          乘法单位元 1        (满足 x · 1 = x)
//        conj(x)        共轭                (实域退化为恒等, 复域翻虚部符号)
//        abs2(x)        |x|²                (返回 real_t, 避免开方)
//        real(x)        Re(x)               (投影到实轴)
//        fromReal(r)    实数 → T 的嵌入     (保证导数/虚部为 0)
//
//        这六个是 linalg 算法**唯一的**标量原语。算法层只用它们,
//        不直接碰 T 的算术细节 —— 于是实/复/Dual 的差异全部集中在 traits 里。
//
// ---------------------------------------------------------------------------
// (S2) 为什么 zero/one 必须是静态函数而不是 constexpr 常数
// ---------------------------------------------------------------------------
//        对 `double`, 0 可以是 constexpr。但对 `Complex<T>` 与 `Dual<T>`,
//        "0" 是**构造出来的对象**, 且构造过程可能要调到基域的 zero()。
//        统一成静态函数, 四种域写法一致, concept 也只需检查一种形式。
//
// ---------------------------------------------------------------------------
// (S3) conj 的三域含义
// ---------------------------------------------------------------------------
//        实数:   conj(x) = x                        ← 恒等
//        复数:   conj(a + bi) = a − bi              ← 翻虚部
//        对偶:   conj(Dual(v,d)) = Dual(conj v, conj d)   ← **逐分量**
//
//        ⚠ 对偶的 conj 必须逐分量。只 conj 一个槽的写法在实域**完全等价于正确版**
//          ⇒ 双精度探针永远看不见, 只有 Dual<Complex> 才暴露。
//          这是一类典型事故: 实域免疫、复域才现形。
//
// ---------------------------------------------------------------------------
// (S4) abs2 而不是 abs
// ---------------------------------------------------------------------------
//              abs2(x) = |x|² = conj(x) · x
//
//        为什么提供平方版:
//        ① 省一次开方 (算法里比较大小, 平方单调所以不影响判定)
//        ② 阈值可以同步取平方: |x| <= tol  ⟺  |x|² <= tol²
//        ③ 开方引入额外误差, 平方没有
//
//        ⚠ abs2 会**扔掉相位信息**: 两个仅差相位的复数被判为相等。
//          所以不能用它判"复向量是否相同", 只能判"模长是否相同"。
//          判复向量相同要逐分量比 real 与 imag 两个槽。
//
// ---------------------------------------------------------------------------
// (S5) real(x) 的两层用途
// ---------------------------------------------------------------------------
//        ① 数学: 取实部 Re(x)
//        ② 工程: 一切**数值判据** (比主元、判收敛、判正定) 都在实数域发生
//
//        对 Dual, real(x) = real(x.value) —— **只看 value 槽, 不带 epsilon**。
//        因为 epsilon 是"导数方向的标记", 拿它参与大小比较没有数学意义。
//
// ---------------------------------------------------------------------------
// (S6) fromReal 的语义 (常数通道)
// ---------------------------------------------------------------------------
//        fromReal(r):  ℝ → T,   且保证 r 被嵌入后"导数/虚部为 0"
//
//        实数:   fromReal(r) = r
//        复数:   fromReal(r) = r + 0i          ← 虚部为 0
//        对偶:   fromReal(r) = Dual(r, 0)       ← 导数为 0 (常数导数为零)
//
//        这就是为什么泛型代码里**禁止裸写字面量**:
//        裸写 `1` 在 Dual 域可能走上某个隐式转换、静默丢掉 epsilon 分量,
//        而 `fromReal(1)` 保证 derive = 0。
//
// ---------------------------------------------------------------------------
// (S7) CommutativeScalar 的数学含义
// ---------------------------------------------------------------------------
//        CommutativeScalar<T>  ⟺  Scalar<T> 且  ∀ a,b:  a·b = b·a
//
//        **为什么 linalg 需要它**: LU / QR / Jacobi 等算法都隐含假设
//        乘法可交换 —— 例如"行交换与消元顺序可互换"、"点积可随意重排"。
//        四元数满足除环但不交换, 一旦传进来, 结果是**编得过、算得错**。
//
//        ⇒ 用 concept 在**入口**挡掉, 而不是在算法里到处加 assert。
//          这是 D16 (四元数不做矩阵元素) 的机制化落实。
//
// ---------------------------------------------------------------------------
// (S8) 容差公式
// ---------------------------------------------------------------------------
//        DefaultTol<T>(n) = 8 · n · ε(real_t)
//
//        · 按**实域**取 ε: 对 Complex<T> / Dual<T> 用的是基域的 ε, 不是 T 的。
//        · 随 n 缩放: 累计误差随运算量增长 (n 步消元 ⇒ n 倍 ε 量级)。
//        · 系数 8: 留一点余量, 不是理论推导出来的常数。
//
//        ⚠ **不要**用 core/common.hpp 的 Tol: 那是 float 的 ε (1.19e-7),
//          对 double 差 9 个数量级, 会把良定矩阵误判成奇异。
//
// ---------------------------------------------------------------------------
// (S9) 嵌套域的递归下钻
// ---------------------------------------------------------------------------
//        Complex<Dual<T>> 与 Dual<Complex<T>> 都是合法域。
//        递归关系:  real_t(Complex<U>) = real_t(U),   real_t(Dual<U>) = real_t(U)
//
//        所以任何嵌套最终都落到 float / double —— 这就是为什么
//        `DefaultTol<T>(n)` 总能拿到一个具体的 ε。
//        前提: 每一层的 `requires Scalar<T>` 保证递归是良基的 (不会无限展开)。
// ===========================================================================
