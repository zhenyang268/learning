#pragma once

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"

/** 列向量约束
 * 1. 左乘矩阵代表坐标系变换 Vy Py = Vx Px, Vy = Vx G, Px = G Py
 * 2. 所以GTAG就代表, 换坐标系, 运算, 再换回去
 */

 

// ---- tag (D14): 结构由调用方声明, 不进类型系统 ----
struct Hermitian {};
struct General {};              // 阶段 2 才实现, 此处只放声明

// ---- 结果 (特征值恒为实数 —— Hermitian 的结构决定, 与标量无关) ----
template <typename T>
struct Eigendecomposition
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T>           Q;       // 特征向量按列存放, n×n, 酉/正交
    std::vector<real_t> values;  // 特征值, 实数
};

// ---- ① 主入口: Q 与 Λ ----
template <CommutativeScalar T>
Eigendecomposition<T> eig(const Matrix<T>& A, Hermitian);

// ---- ② 只要谱 ----
template <CommutativeScalar T>
std::vector<typename ScalarTraits<T>::real_t>
eigenvalues(const Matrix<T>& A, Hermitian);

// ---- ③ 预留: 阶段 2 (非对称可能有复特征值), 只声明不实现 ----
// Eigendecomposition<T> eig(const Matrix<T>& A, General);
// std::vector<Complex<real_t>> eigenvalues(const Matrix<T>& A, General);