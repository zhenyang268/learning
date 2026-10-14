#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"

// 落点说明: 按 plan.md §2, operator.hpp 的正式职责是 LinearOperator(线性算子)
// + 基/边界条件工厂; Dot/Abs2 是元素级矩阵工具, 暂放此处供 qr.hpp 使用
// (S6 派支线处置时再定去留 —— 见 plan.md S6)。

// Dot: Frobenius 内积 <a, b> = Σ conj(a_i)·b_i
// 复数域必须共轭**第一个**参数 —— 与 abs2 同口径 (D7);
// 漏 conj 是"double 免疫"老坑 (实矩阵下 conj 恒等, 双精度探针永远看不到)。
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
// 一切比较走 abs2、阈值也用平方 (D12: abs2 <= tol*tol), 避免开方误差。
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t Abs2(const Matrix<T>& m)
{
    typename ScalarTraits<T>::real_t result = 0;
    for (int i = 0; i < m.Rows() * m.Cols(); i++) {
        result = result + ScalarTraits<T>::abs2(m.Data()[i]);
    }
    return result;
}