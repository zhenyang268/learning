#pragma once

#include <cstring>
#include <cassert>

/*
#include <numbers>
#include <limits>

constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
constexpr double NaN_D = std::numeric_limits<double>::quiet_NaN();
constexpr float Inf = std::numeric_limits<float>::infinity();
constexpr double Inf_D = std::numeric_limits<double>::infinity();

constexpr float PI = std::numbers::pi_v<float>;
constexpr double PI_D = std::numbers::pi_v<double>;

#define DEG2RAD(degree) (degree) * PI / 180.0f
#define RAD2DEG(rad) (rad) * 180.0f / PI

constexpr float EPSILON = std::numeric_limits<float>::epsilon();
constexpr double EPSILON_D = std::numeric_limits<double>::epsilon();

constexpr float Tol = EPSILON;
constexpr double Tol_D = EPSILON_D;
*/

enum class MatrixKind 
{ 
    Generic,
    Diagonal,
    Symmetric, // 酉对称
    Orthogonal, // 正交矩阵
    RigidTransform,
    SymmetricPositiveDefinite, /*Hermitian正定*/
    Projection 
};

// 单元操作直接定义在Dual类的友元函数中, Dual.sin
// 矩阵级约简/逐元素应用在 autodiff/jacobian.hpp (SumAll/MultAll/Apply)
enum class GradOp
{
    Input,          // 输入节点
    Constant,       // 常数/常量矩阵节点
    Output,         // 输出透传节点
    Add,
    Sub,
    Mult,           // 逐元素乘 (标量/同形)
    CWiseMult,      // 逐元素乘法
    CWisePow,       // 逐元素幂运算
    Divide,
    MatrixMult,     // 可以实现mask和choose, 比如x 10维,可以只取0 1 2三维
    Transpose,      // 转置是线性变换
    Inverse,        // 求逆用伴随内积配对不变
    Sigmoid,
    Relu,
    Tanh,
    Softmax,
    CrossEntropy,    // 交叉熵, 损失函数, 这种就是计算图优化, 我们拆开
    SumAll,          // 全元素求和
    MultAll,         // 全元素累乘
    LnMultAll,       // ln版累乘函数
    Function,        // 对所有的input进行相同的操作
    MatrixFunction,  // 使用不同的函数应用到不同的input上
};

// 一元函数表, 只用在 Matrix<FunctionType>(函数节点) 这种情况里面
enum class FunctionType
{
    Sin,
    Cos,
    Log,
    Exp,
    Sqrt
};
// ===========================================================================
// 数学公式与符号约定
// ===========================================================================
//
// 本文件提供的是**常数与宏**, 不是算法, 所以公式很少但都很关键。
//
// ---------------------------------------------------------------------------
// 符号表
// ---------------------------------------------------------------------------
//   π (PI)         圆周率 3.14159265358979323846
//   ε (EPSILON)    机器 epsilon: 1.0 与"比它大的最小浮点数"之差
//   Tol            本文件的历史容差常数 (= EPSILON, **不建议新代码使用**)
//   degree, rad    角度 / 弧度
//
// ---------------------------------------------------------------------------
// (M1) 角度与弧度
// ---------------------------------------------------------------------------
//        rad = degree · π / 180
//        degree = rad · 180 / π
//
//      代码里对应 DEG2RAD / RAD2DEG 两个宏。
//      ⚠ 这两个宏**没有加括号保护** (见实现处), 传 `a+b` 进去会算错。
//        新代码建议直接用 `x * PI_D / 180.0` 手写, 或自己包一个 constexpr 函数。
//
// ---------------------------------------------------------------------------
// (M2) 机器 epsilon 的定义
// ---------------------------------------------------------------------------
//        ε = min{ δ > 0 : 1.0 + δ ≠ 1.0 }        (在给定浮点类型下)
//
//        含义: 相对精度极限。float 的 ε ≈ 1.19e-7, double 的 ε ≈ 2.22e-16。
//        **两者差 9 个数量级** —— 这是本文件 Tol 不能被误用的根本原因。
//
// ---------------------------------------------------------------------------
// (M3) 为什么本文件的 Tol 不建议使用 (踩过的坑)
// ---------------------------------------------------------------------------
//        Tol = EPSILON = float 的 ε (1.19e-7)
//
//        如果拿它判 double 矩阵的奇异性, 会把**良定矩阵误判成奇异**
//        (double 的残渣在 1e-16 量级, 离 1e-7 还很远)。
//
//        ⇒ 正确做法是用 `DefaultTol<T>(n)` (定义在 core/scalar.hpp):
//              DefaultTol<T>(n) = 8 · n · ε(real_t)
//          它是**按 T 的实域取 ε**, 且随规模 n 缩放。
//
//        本文件的 Tol / Tol_D 是历史遗留, 仅为兼容老代码保留。
//
// ---------------------------------------------------------------------------
// (M4) 为什么不用 M_PI (踩过的坑)
// ---------------------------------------------------------------------------
//        M_PI 不是标准 C++ —— 它在 <cmath> 里由 POSIX 扩展提供,
//        在严格标准模式下 (`-std=c++20` 不带 gnu 扩展) 可能不存在。
//        本文件用 `std::numbers::pi_v<T>` (C++20 标准), 保证可移植。
//
// ---------------------------------------------------------------------------
// 实现口径
// ---------------------------------------------------------------------------
//   · PI 是 float 精度, PI_D 是 double 精度 —— **别混用**。
//     用 float 精度算 double 的三角函数会引入 1e-8 级误差。
//   · EPSILON / EPSILON_D 分别对应两种精度, 名字上带 _D 的是 double 版。
// ===========================================================================
