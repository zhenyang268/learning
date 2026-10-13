#pragma once

#include <cmath>
#include <vector>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"

/** 列向量约束
 * 1. 左乘矩阵代表坐标系变换 Vy Py = Vx Px, Vy = Vx G, Px = G Py
 * 2. 所以GTAG就代表, 换坐标系, 运算, 再换回去
 */

// ---- 结果 (特征值恒为实数 —— Hermitian 的结构决定, 与标量无关) ----
template <typename T>
struct Eigendecomposition
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T>           Q;       // 特征向量按列存放, n×n, 酉/正交
    Matrix<T>           Lamda;   // 特征值, 实数
};

// ---- ① 主入口: Q 与 Λ ----
template <CommutativeScalar T>
Eigendecomposition<T> eig(const Matrix<T>& A)
{
    if (is_symmetric(A)) {
        return eig_symmetric(A);
    }
    else {
        return eig_general(A);
    }
}


template <CommutativeScalar T>
static void rotate(Matrix<T> &A, Matrix<T> &R)
{
    Matrix<T> RT = R.Transpose();
    A = R * A * RT;
}

template <CommutativeScalar T>
Eigendecomposition<T> eig_symmetric(const Matrix<T>& A)
{
    Eigendecomposition<T> res;
    res.Q = Matrix<T>::Identity(A.Rows(), A.Cols());
    res.Lamda = A;
    Matrix<T>& L = res.Lamda;
    Matrix<T> R = Matrix<T>::Identity(A.Rows(), A.Cols());

    for (int i = 0; i < A.Rows(); i++) {
        for (int j = i+1; j < A.Cols(); j++) {
            double l1, l2, theta;
            jacobi(L(i, i), L(i, j), L(j, j), l1, l2, theta);
            R(i, i) = R(j, j) = std::cos(theta);
            R(i, j) = -std::sin(theta);
            R(j, i) = std::sin(theta);
            // 应用RTLR, 影响两行两列
            rotate(L, R);
            // 已经把i, j变成0了, R1R2... Rn就是要找的正交矩阵, 也就是特征向量
            res.Q = res.Q * R;

            R(i, i) = R(j, j) = ScalarTraits()::one();
            R(i, j) = R(j, i) = ScalarTraits()::zero();
        }
    }
    // 最后res.lamda就是转换完的对角矩阵, 特征值
    // res.Q存的就是特征向量
    return res;
}

// A = [[a, b], [b, c]]，2x2 实对称，正交对角化
// tan(2*theta) = 2b / (a - c)
void jacobi(double a, double b, double c,
              double& l1, double& l2, double& theta)
{
    theta = 0.5 * std::atan2(2.0 * b, a - c);        // 旋转角，无特判
    double half = 0.5 * (a + c);
    double d    = std::hypot(0.5 * (a - c), b);      // 比手写 sqrt 更稳，防溢出
    l1 = half + d;
    l2 = half - d;
    // 校验：l1+l2 == a+c,  l1*l2 == a*c-b*b
}


// ---- ③ 预留: 阶段 2 (非对称可能有复特征值), 只声明不实现 ----
// Eigendecomposition<T> eig(const Matrix<T>& A, General);
// std::vector<Complex<real_t>> eigenvalues(const Matrix<T>& A, General);