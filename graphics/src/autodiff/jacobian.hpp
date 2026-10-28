#pragma once

#include <vector>
#include <stdexcept>
#include <functional>
#include <string>       // CheckFuncMatrixShape 的报错信息
#include <cmath>        // 标量域的 sin/cos/exp/log/sqrt (elemFunc)
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "autodiff/dual.hpp"

// DualFunc: 逐元素函数类型 (目标签名 T(const T&))
// 当前函数节点走 FunctionType 枚举, DualFunc 暂无使用者 —— 别名占位
template <typename T>
using DualFunc = std::function<T(const T&)>;

// ---- 标量域逐元素分派 ----
// `dualFunc` (dual.hpp) 是 Dual 域专用: 它是 Dual<T> 的 in-class friend, 只能经 ADL 找到,
// 所以标量域根本看不见这个名字。标量域 (double / float / Complex<...>) 走这里。
// 函数名一律用 log/sin/... (标量域口径); Complex<T> 的数学函数同样是 in-class friend,
// 靠 `using std::xxx;` 引入 + ADL 命中。Dual 域不会走到本函数 (由下面的 Dual 专用重载接管)。
template <typename T>
T elemFunc(T x, FunctionType f)
{
    switch (f) {
        case FunctionType::Sin:  { using std::sin;  return sin(x);  }
        case FunctionType::Cos:  { using std::cos;  return cos(x);  }
        case FunctionType::Exp:  { using std::exp;  return exp(x);  }
        case FunctionType::Log:  { using std::log;  return log(x);  }
        case FunctionType::Sqrt: { using std::sqrt; return sqrt(x); }
    }
    return x;   // FunctionType 已全覆盖; 此行只为满足返回路径 (将来加枚举时 -Wswitch 会报警)
}

// ---- 函数矩阵的形状契约 (2026-10-19, 列向量翻案时修的真 bug) ----
//
// 旧版比较的是 `A.Cols() != f.Cols()` —— 那是**行向量**约定下的"长度"。列向量下
// A 是 n×1、f 也是 n×1, Cols() 恒为 1, 这个检查退化成 "1 == 1 永远通过";
// 于是 n 与 f 的长度不一致时**不抛**, 一路走到 `f.Data()[i]` 越界读 (实测确认)。
//
// 现契约与向量方向无关:
//   · f 是 1x1  -> 广播到 A 的每一格
//   · 否则      -> 必须与 A **同形** (Rows 与 Cols 都比)
// 这样行向量/列向量/批量矩阵三种输入都被同一个判据覆盖。
template <typename T, typename U>
void CheckFuncMatrixShape(const Matrix<T>& A, const Matrix<U>& f, const char* who)
{
    if (f.Rows() == 1 && f.Cols() == 1) return;          // 单函数: 广播
    if (f.Rows() != A.Rows() || f.Cols() != A.Cols()) {
        throw std::invalid_argument(
            std::string(who) + ": 函数矩阵形状与 A 不一致 (A " +
            std::to_string(A.Rows()) + "x" + std::to_string(A.Cols()) + ", f " +
            std::to_string(f.Rows()) + "x" + std::to_string(f.Cols()) +
            "; 只允许 1x1 广播或同形)");
    }
}

// ---- Apply 的重载分裂 (2026-10-16) ----
// 原先只有一个 `Apply(Matrix<T>, ...)`, 体内却假定 T = Dual<U> ⇒ 标量域一调用就报
// `error: use of undeclared identifier 'dualFunc'` (错误信息把人引向"是不是漏 include", 误导)。
// 现按域拆成两个重载: Dual 专用 + 标量专用。共存不需要任何 hack —— `Matrix<Dual<U>>`
// 比 `Matrix<T>` 更特化, [temp.deduct.partial] 的部分排序必然把 Dual 域判给 Dual 版。
// 标量版体内的 static_assert 是守卫: 万一偏序选错, 立刻编不过, 而不是静默走错路径。

// (1) Dual 专用: 走 dualFunc, 传播导数
template <typename T>
Matrix<Dual<T>> Apply(const Matrix<Dual<T>>& A, FunctionType f)
{
    Matrix<Dual<T>> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = dualFunc(A.Data()[i], f);
    }
    return B;
}

// (2) 标量专用: 走 elemFunc, 只算值 (返回类型里没有 Dual)
template <typename T>
Matrix<T> Apply(const Matrix<T>& A, FunctionType f)
{
    static_assert(!IsDual<T>::value,
                  "Apply: Dual 域应由 Dual 专用重载接管 —— 偏序选错了");
    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = elemFunc(A.Data()[i], f);
    }
    return B;
}

// 逐元素应用函数矩阵 (MatrixFunction 节点前向): 每格用自己的函数
// (3) Dual 专用
template <typename T>
Matrix<Dual<T>> Apply(const Matrix<Dual<T>>& A, const Matrix<FunctionType>& f)
{
    if (f.Rows() == 1 && f.Cols() == 1) {
        return Apply(A, f(0, 0)); // 单函数退化为统一应用 (落回 (1))
    }

    CheckFuncMatrixShape(A, f, "Apply");

    Matrix<Dual<T>> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = dualFunc(A.Data()[i], f.Data()[i]);
    }
    return B;
}

// (4) 标量专用
template <typename T>
Matrix<T> Apply(const Matrix<T>& A, const Matrix<FunctionType>& f)
{
    static_assert(!IsDual<T>::value,
                  "Apply: Dual 域应由 Dual 专用重载接管 —— 偏序选错了");
    if (f.Rows() == 1 && f.Cols() == 1) {
        return Apply(A, f(0, 0)); // 单函数退化为统一应用 (落回 (2))
    }

    CheckFuncMatrixShape(A, f, "Apply");

    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = elemFunc(A.Data()[i], f.Data()[i]);
    }
    return B;
}

// 把 value 赋值为 derive, 用于 backward (Function 节点: 局部导数暂存 derive)
// 按值收参并返回: 调用点写的是 Derive(Apply(...)), 原先收 Matrix<T>& 接不了临时量
// (Backward 一实例化就报 no matching function for call to 'Derive')
template <typename T>
Matrix<T> Derive(Matrix<T> A)
{
    for (int i = 0; i < A.Size(); i++) {
        A.Data()[i].value = A.Data()[i].derive;
    }
    return A;
}

// 局部导数 f'(x): **直接在内层标量上算**, 完全不碰 Dual 的 derive 分量。
//
// 为什么不能用 `Derive(Apply(A, f))`: 那个结果是 f'(x)·A.derive, **依赖外部种子**
// (实测: 种子 0 ⇒ 全 0; 种子 3 ⇒ 3 倍)。反向模式不是"沿种子传播", derive 在这里没有位置。
// 反向时重算 f'(x.value) 与"前向时存下来"**等价**(纯函数), 只多花一次 O(n) 标量计算 ——
// 省掉那次计算是性能优化, 不是正确性要求。
template <typename U>
U localFunc(U x, FunctionType f)
{
    switch (f) {
        case FunctionType::Sin:  { using std::cos;  return  cos(x);  }
        case FunctionType::Cos:  { using std::sin;  return -sin(x);  }
        case FunctionType::Exp:  { using std::exp;  return  exp(x);  }
        case FunctionType::Log:  return ScalarTraits<U>::one() / x;
        case FunctionType::Sqrt: { using std::sqrt; return ScalarTraits<U>::fromReal(0.5) / sqrt(x); }
    }
    return x;
}

template <typename T>
Matrix<Dual<T>> LocalDeriv(const Matrix<Dual<T>>& A, FunctionType f)
{
    Matrix<Dual<T>> D(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++)
        D.Data()[i] = Dual<T>(localFunc(A.Data()[i].value, f), ScalarTraits<T>::zero());
    return D;
}

// 逐格不同函数的版本; 1x1 函数矩阵退化为统一应用 (与 Forward 侧的 Apply 同规则)
template <typename T>
Matrix<Dual<T>> LocalDeriv(const Matrix<Dual<T>>& A, const Matrix<FunctionType>& f)
{
    if (f.Rows() == 1 && f.Cols() == 1) return LocalDeriv(A, f(0, 0));
    Matrix<Dual<T>> D(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++)
        D.Data()[i] = Dual<T>(localFunc(A.Data()[i].value, f.Data()[i]), ScalarTraits<T>::zero());
    return D;
}

// ---- Matrix 级约简包装 ----
// Dual 级原语 sumLine/multLine 在 dual.hpp (Dual 域约简, dual.hpp 不能
// 反向 include matrix.hpp, 故 Matrix 级包装只能放在本文件)

// 全元素和: 返回 1 个 Dual(和, 全微分)
template <typename T>
Dual<T> SumAll(const Matrix<Dual<T>>& A)
{
    std::vector<Dual<T>> flat(A.Data(), A.Data() + A.Size());
    return sumLine(flat);
}

// 全元素积: 返回 1 个 Dual(总积, 全微分)
template <typename T>
Dual<T> MultAll(const Matrix<Dual<T>>& A)
{
    std::vector<Dual<T>> flat(A.Data(), A.Data() + A.Size());
    return multLine(flat);
}

// MultAll 的反向: 每格 = Dual(总积, ∂总积/∂x_i = Π_{j≠i} x_j)
// 同一套前后缀积 (O(n) 空间), 零因子安全
template <typename T>
Matrix<Dual<T>> MultAllBackward(const Matrix<Dual<T>>& A)
{
    const int n = A.Size();
    Matrix<Dual<T>> out(A.Rows(), A.Cols());
    if (n == 0) return out;

    std::vector<Dual<T>> flat(A.Data(), A.Data() + A.Size());
    std::vector<T> pre(n), suf(n);   // pre[i] = Π_{j<i} v_j, suf[i] = Π_{j>i} v_j
    pre[0] = ScalarTraits<T>::one();
    for (int i = 1; i < n; i++) pre[i] = pre[i - 1] * flat[i - 1].value;
    suf[n - 1] = ScalarTraits<T>::one();
    for (int i = n - 2; i >= 0; i--) suf[i] = suf[i + 1] * flat[i + 1].value;

    T total = pre[n - 1] * flat[n - 1].value;
    for (int i = 0; i < n; i++) {
        out.Data()[i] = Dual<T>(total, pre[i] * suf[i]);
    }
    return out;
}

// ===========================================================================
// 分子布局雅可比 (D19 第二次修订, 2026-10-19)
// ===========================================================================
//
//   f: R^n -> R^m, x 是 n×1 **列**向量, y = f(x) 也是 m×1 列向量
//   J_ij = ∂y_i/∂x_j, 形状 **m×n** (行 = 输出分量, 列 = 输入分量)
//
// 一次前向 sweep (种子 e_j) 填 J 的**第 j 列** —— 这是与旧口径的唯一实质差别
// (旧: 分母布局 n×m, 一次 sweep 填第 i 行)。
//
// 为什么翻案 (D19 第二次修订的三条理由):
//   1) 链式法则回到自然序 J_{f∘g} = J_f·J_g (旧口径是 J_g·J_f, 镜像);
//   2) 标量输出 f: R^n -> R 时 J 是 1×n 行, 梯度 = J^T 是 n×1, 与被微对象 x 同形;
//   3) autodiff 与 linalg 从此同一套口径, 交界处不再需要显式转置。
//
// 判据 (demo/autodiff/jacobian.cpp): 解析 / 中心差分 / 图反向 三方对拍。

template <typename T, typename F>
Matrix<T> Jacobian(F f, const Matrix<Dual<T>>& x, int m)
{
    if (x.Cols() != 1) {
        throw std::invalid_argument(
            "Jacobian: x 必须是 n×1 列向量 (收到 " + std::to_string(x.Rows()) +
            "x" + std::to_string(x.Cols()) + ")");
    }
    if (m < 0) {
        throw std::invalid_argument("Jacobian: m 不能为负");
    }

    const int n = x.Rows();
    Matrix<T> J(m, n);

    for (int j = 0; j < n; j++) {
        Matrix<Dual<T>> seed = x;                        // 拷贝: 不改调用方的 x
        for (int i = 0; i < seed.Size(); i++) {
            seed.Data()[i].derive = ScalarTraits<T>::zero();
        }
        seed.Data()[j].derive = ScalarTraits<T>::one();  // 种子 e_j

        const Matrix<Dual<T>> y = f(seed);
        if (y.Cols() != 1 || y.Rows() != m) {
            throw std::invalid_argument(
                "Jacobian: 输出必须是 m×1 列向量 (声明 m=" + std::to_string(m) +
                ", 实得 " + std::to_string(y.Rows()) + "x" +
                std::to_string(y.Cols()) + ")");
        }
        for (int i = 0; i < m; i++) {
            J(i, j) = y.Data()[i].derive;                // 填**第 j 列**
        }
    }
    return J;
}

// 不给 m 的版本: 先用全零种子跑一次拿输出形状 (顺手校验输出是列向量)
template <typename T, typename F>
Matrix<T> Jacobian(F f, const Matrix<Dual<T>>& x)
{
    if (x.Cols() != 1) {
        throw std::invalid_argument(
            "Jacobian: x 必须是 n×1 列向量 (收到 " + std::to_string(x.Rows()) +
            "x" + std::to_string(x.Cols()) + ")");
    }
    Matrix<Dual<T>> probe = x;
    for (int i = 0; i < probe.Size(); i++) {
        probe.Data()[i].derive = ScalarTraits<T>::zero();
    }
    const Matrix<Dual<T>> y0 = f(probe);
    if (y0.Cols() != 1) {
        throw std::invalid_argument(
            "Jacobian: 输出必须是列向量 (实得 " + std::to_string(y0.Rows()) +
            "x" + std::to_string(y0.Cols()) + ")");
    }
    return Jacobian(f, x, y0.Rows());
}

// ===========================================================================
// 四、数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里, 实现处只留"为什么这么写"的工程理由。**
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   f: Rⁿ → Rᵐ     向量值函数 (代码里是模板参数 F)
//   x              n×1 输入列向量
//   y = f(x)       m×1 输出列向量
//   J              雅可比矩阵, m×n
//   J_ij           ∂y_i / ∂x_j
//   m, n           输出的分量个数 / 输入的分量个数
//   e_j            第 j 个标准基向量 (n×1)
//   Dual<T>        对偶数 (见 autodiff/dual.hpp)
//   A              元素级运算的输入矩阵
//   SumAll/MultAll 矩阵级约简 (求和 / 连乘)
//
// ---------------------------------------------------------------------------
// (J1) 雅可比矩阵的定义 (分子布局, D19)
// ---------------------------------------------------------------------------
//        J ∈ R^{m×n},     J_ij = ∂y_i / ∂x_j
//
//        **行 = 输出分量, 列 = 输入分量** —— 这是分子布局 (numerator layout)。
//
//        旧口径 (分母布局) 是 J_ij = ∂y_j/∂x_i, 形状 n×m, 是这里的转置。
//        全库已于 2026-10-19 统一为分子布局 (D19 第二次修订)。
//
// ---------------------------------------------------------------------------
// (J2) 为什么分子布局 (链式法则的自然序)
// ---------------------------------------------------------------------------
//        分子布局:  J_{f∘g} = J_f · J_g          ← 自然序, 与复合书写同向
//        分母布局:  J_{f∘g} = J_g · J_f          ← 镜像
//
//        上面那条"与书写同向"的性质在大复合链里非常重要 ——
//        分母布局下每写一层复合就要把乘法顺序倒过来, 人工推导极易出错。
//
// ---------------------------------------------------------------------------
// (J3) 前向模式如何算 J 的第 j 列 (本文件的算法)
// ---------------------------------------------------------------------------
//        给 x 播种第 j 个分量:  x = x₀ + e_j·ε
//        则 f(x) = f(x₀) + (∂f/∂x_j)·ε
//
//      ⇒ **一次前向求值就得到 J 的整列** (m 个数)。
//        要整个 J 就跑 n 次 (每个 j 一次) —— 这是前向模式的固有代价 O(n)。
//
//        对偶: 反向模式一次得到**整行**, 代价 O(m) (见 grad_graph.hpp)。
//        ⇒ 输入少输出多 (m >> n) 用前向; 输入多输出少 (n >> m) 用反向。
//
// ---------------------------------------------------------------------------
// (J4) 标量函数的特殊情况
// ---------------------------------------------------------------------------
//        若 m = 1 (f: Rⁿ → R), 则 J 是 **1×n 行向量**:
//              J = [ ∂f/∂x_1  ∂f/∂x_2  ...  ∂f/∂x_n ]
//
//        而**梯度**定义为 J 的转置 (n×1 列向量):
//              ∇f = Jᵀ
//
//        ⚠ 梯度和雅可比**不是同一个东西**: 一个是列向量一个是行向量。
//          混用会在与其它矩阵相乘时得到形状不匹配或静默转置的错误。
//
// ---------------------------------------------------------------------------
// (J5) 元素级运算的雅可比 (Apply 的语义)
// ---------------------------------------------------------------------------
//        若 f 逐元素作用 (y_i = φ(x_i)), 则
//              J = diag( φ'(x_1), ..., φ'(x_n) )      ← **对角阵**
//
//        这就是为什么 Apply 的雅可比"看起来像向量" —— 非对角元恒为 0。
//        实现里不构造这个对角阵 (浪费 n² 内存), 只存对角。
//
//        广播情形 (f 是 1×1 而 A 是 m×n):  J = φ'(a) · 全 1 矩阵。
//
// ---------------------------------------------------------------------------
// (J6) 链式法则的矩阵形式 (为什么形状必须严格检查)
// ---------------------------------------------------------------------------
//        y = f(g(x)):   J = J_f(g(x)) · J_g(x)
//
//        形状: (m×k)·(k×n) = (m×n)   ← 中间那个 k 必须对得上
//
//      ⇒ 矩阵乘的形状检查不是"防御性编程", 而是**链式法则成立的必要条件**。
//        这里形状错 = 数学上无意义, 不是"碰巧能算"。
//
// ---------------------------------------------------------------------------
// (J7) 有限差分的对拍判据 (独立于实现路径)
// ---------------------------------------------------------------------------
//        J_ij ≈ ( f_i(x + h·e_j) − f_i(x − h·e_j) ) / (2h)      ← 中心差分
//
//        误差 ~ O(h²)·f''' + O(eps/h)。
//        h = 1e-6 (double) 时两项都在 1e-10 量级, 与解析雅可比差 < 1e-8 可判过。
//
//        ⚠ 差分**只读 .value**, 不碰 derive 槽 —— 它才能当独立判据。
//          如果用 Dual 的 derive 去验 Dual 的 derive, 就只是自证。
//
// ---------------------------------------------------------------------------
// 实现口径 (跟公式直接相关)
// ---------------------------------------------------------------------------
//   · **x 必须是 n×1 列向量** (D19); 输出必须是 m×1。形状不对 throw。
//   · 播种时**拷贝**输入 (不改调用方的 x), 每次只把第 j 个 derive 置 1。
//   · 无 m 版本: 先用全零种子探一次拿输出形状, 再按 m 跑 —— 多一次求值,
//     换接口不必让调用方手动给 m。
// ===========================================================================
