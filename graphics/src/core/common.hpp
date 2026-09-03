#pragma once

#include <cmath>
#include <cstring>
#include <iostream>
#include <initializer_list>

constexpr float PI = 3.14159265358979323846f;
#define DEG2RAD(degree) #degree * PI / 180.0f
#define RAD2DEG(rad) #rad * 180.0f / PI
constexpr float EPSILON = 1e-6f;

enum MatrixType
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

// 单元操作直接定义在Dual类的友元函数中, Dual.sin
// 多元函数操作定义在Matrix类的友元函数中, Matrix<T>.SumAll
enum GradOp
{
    Input,          // 输入行
    Output,         // 输出行
    Constant,       // 线性关系中的常数, 按参数处理
    Add,
    Sub,
    Mult,           // 逐元素乘法
    Devide,
    MatrixMult,     // 可以实现mask和choose, 比如x 10维,可以只取0 1 2三维
    Transpose,
    Inverse,
    Sigmoid,
    Relu,
    Tanh,
    Softmax,
    CrossEntropy,   // 交叉熵, 损失函数, 这种就是计算图优化, 我们拆开
    SumAll,         // 求和函数
    MultAll,        // 累乘函数
    LnMultAll,      // ln版累乘函数
    Function,       // 对所有的input进行相同的操作
    MatrixFunction, // 使用不同的函数应用到不同的input上
};
