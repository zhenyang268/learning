#pragma once

// ===========================================================================
// linalg/operator.hpp —— 矩阵的 Frobenius 内积与范数
// ===========================================================================
//
// 【状态】已实现 (暂驻: 按设计本文件将来是 LinearOperator 的家, 见末尾)。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <CommutativeScalar T> T Dot(const Matrix<T>& a, const Matrix<T>& b);
//   template <CommutativeScalar T>
//   typename ScalarTraits<T>::real_t Abs2(const Matrix<T>& m);
//
//   Dot     Frobenius 内积 <a,b> = Σ conj(a_i)·b_i, 返回**标量域 T**
//   Abs2    ‖m‖_F² = Σ |m_i|², 返回**实数域 real_t**
//
// 形状: Dot 要求 a 与 b **元素个数与形状都相同**, 否则 throw std::invalid_argument。
//       Abs2 对任意形状都成立。
//
// 关键口径:
//   · Dot **共轭第一个参数** —— 与 abs2 同口径 (numeric-contract.md D7)。
//   · 一切比较走 Abs2、阈值也取平方 (abs2 <= tol*tol), 避免开方引入误差。
//     要真实 Frobenius 范数时才 sqrt, 且只在打印/打分路径。
//
// ---------------------------------------------------------------------------
// 二、内部机制
// ---------------------------------------------------------------------------
//   两者都是单层循环, 按 Data() 线性扫描 (不按 (i,j) 走, 免去乘法)。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【Dot 必须共轭第一个参数 —— "double 免疫"老坑】
//   漏掉 conj 在实矩阵下完全看不出 (double 的 conj 恒等), 只有复数域才暴露。
//   后果不是"少了点精度", 而是内积失去正定性 (<x,x> 不再是实数范数),
//   依赖它的判据 (正交性、投影、最小二乘) 会静默给出错误结论。
//   ⇒ 复数域的判据必须**至少跑一次 Complex 域**; 只跑 double 等于没测。
//
// 【为什么 Abs2 返 real_t 而 Dot 返 T】
//   Dot 的结果在标量域 (复数内积是复数), 而 Abs2 按定义是实数。
//   写成同一返回类型会让其中一个被迫做无意义转换, 并在 Dual 域丢掉 epsilon。
//
// 【为什么这两个暂住在这里】
//   按分层, 它们属于 linalg 的通用工具, 而 operator.hpp 的正式职责是
//   LinearOperator (线性算子 = 基 + 边界条件 + 截断矩阵)。
//   在 LinearOperator 落地前, 先放这里供 qr.hpp 使用; 届时再定去留。
//   判据: 加 LinearOperator 时若发现这俩更该在别处, 搬家即可 —— 它们没有
//   反向依赖 (qr.hpp 是唯一消费者)。
// ===========================================================================

#include <cmath>
#include <stdexcept>
#include <vector>

#include "core/common.hpp"
#include "core/complex.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"

// Dot: Frobenius 内积。共轭**第一个**参数 —— 理由见文件头「三」。
template <CommutativeScalar T>
T Dot(const Matrix<T>& a, const Matrix<T>& b)
{
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) {
        throw std::invalid_argument("Dot: dimension mismatch");
    }
    T result = ScalarTraits<T>::zero();
    for (int i = 0; i < a.Rows() * a.Cols(); i++) {
        result = result + ScalarTraits<T>::conj(a.Data()[i]) * b.Data()[i];   // 无复合赋值 (D7)
    }
    return result;
}

// Abs2: ‖M‖_F² = Σ |m_i|² (元素级平方和)。要真实 Frobenius 范数再开方;
// 一切比较走 abs2、阈值也用平方 (abs2 <= tol*tol), 避免开方误差。
// 为什么返 real_t 而不是 T —— 见文件头「三」。
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t Abs2(const Matrix<T>& m)
{
    typename ScalarTraits<T>::real_t result = 0;
    for (int i = 0; i < m.Rows() * m.Cols(); i++) {
        result = result + ScalarTraits<T>::abs2(m.Data()[i]);
    }
    return result;
}