#pragma once

#include <cmath>
#include <cstring>
#include <iostream>
#include <initializer_list>
#include <numbers>
#include <cassert>
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


// 不能按用途去分, 怎么用是上层的责任, 可以提供静态工厂辅助生成
enum class MatrixType
{
    General,                // 通用矩阵
    Rotation,               // 旋转矩阵
    Scale,                  // 缩放矩阵
    Translation,            // 平移矩阵
    RigidTransform,         // 刚体变换, 旋转平移
    Projection,             // 投影矩阵
    IllConditioned,         // 病态矩阵
    SymmetricPositive,      // 对称正定矩阵
};

enum class MatrixKind 
{ 
    Generic,
    Diagonal,
    Orthogonal, /*酉*/
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