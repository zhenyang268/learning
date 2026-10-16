#pragma once

#include <vector>
#include <stdexcept>
#include <functional>
#include <cmath>        // 标量域的 sin/cos/exp/log/sqrt (elemFunc)
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "autodiff/dual.hpp"

// DualFunc: 计划中的逐元素函数类型 (matrix.md §7.1 第 5 条, 目标签名 T(const T&))
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

    // 应该只有一维, 但是后续计算图节点打包可能涉及多行
    if (A.Cols() != f.Cols()) {
        throw std::invalid_argument("Matrix apply dimension mismatch!");
    }

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

    // 应该只有一维, 但是后续计算图节点打包可能涉及多行
    if (A.Cols() != f.Cols()) {
        throw std::invalid_argument("Matrix apply dimension mismatch!");
    }

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

// ---- Matrix 级约简包装 (matrix.md §7.1 第 4 条) ----
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
