#pragma once

#include <stdexcept>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "linalg/operator.hpp"

/**
 * QR分解
 * 
 * 1. householder QR
 * 2. givens QR
 * 3. shimit QR
 * 
 * @param A 矩阵
 */

template <CommutativeScalar T>
struct QRcomposition
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T>           Q;       // 特征向量按列存放, m*m, 酉/正交
    Matrix<T>           R;      // 上三角矩阵
};

/**
 * householder QR分解
 * Hx = alpha * e1, e1 = [1, 0, 0, ..., 0]T
 * H = I - 2 * v * vT / vT * v
 * Hv = -v  分解 x = v + v_ 所以 x - Hx // v
 * 取v = x - alpha * e1 并且 alpha = -||x||即可
 */
template <CommutativeScalar T>
QRcomposition<T> QR_Householder(const Matrix<T>& A)
{
    Matrix<T> X; // 记录每一个列向量
    QRcomposition<T> res;
    res.Q = Matrix<T>::Identity(A.Rows(), A.Rows());
    res.R = A;
    real_t alpha;

    for (int j = 0; j < A.Cols(); j++) {
        X = Matrix<T>(A.Rows(), 1);
        alpha = 0; // 存VHV, ||V||
        for (int i = j; i < A.Rows(); i++) {
            X(i, 0) = res.R(i, j);
        }
        aplha = Abs2(X);
        
        X(j, 0) = X(j, 0) + sqrt(alpha);
        Matrix<T> H = Matrix<T>::Identity(A.Rows(), A.Rows());
        for (int i = j; i < H.Rows(); i++) {
            for (int k = j; k < H.Cols(); k++) {
                H(i, k) -= 2 * X(i, 0) * ScalarTraits<T>::conj(X(k, 0)) / alpha;
            }
        }
        res.R = H * res.R;
        // 这里j,j下面应该全是0了
        res.Q = res.Q * H.Transpose();
    }

}

template <CommutativeScalar T>
QRcomposition<T> QR_givens(const Matrix<T>& A)
{
    QRcomposition<T> res;
    res.Q = Matrix<T>::Identity(A.Rows(), A.Rows());
    res.R = A;
    real_t thelta;
    Matrix<T> QT = Matrix<T>::Identity(A.Rows(), A.Rows());

    for (int i = 0; i < A.Rows(); i++) {
        for (int j = i+1; j < i; j++) {
            // 左乘 [cos, sin] [-sin, cos]这样把左下角变成0
            thelta = atan2(res.R(i, j)/res.R(j, j));
            QT(i, i) = cos(thelta);
            QT(i, j) = -sin(thelta);
            QT(j, i) = sin(thelta);
            QT(j, j) = cos(thelta);
            res.R = QT * res.R;

            QT(i, j) = -QT(i, j);
            QT(j, i) = -QT(j, i);
            res.Q = res.Q * QT;

            QT(i, i) = ScalarTraits<T>::One();
            QT(j, j) = ScalarTraits<T>::One();
            QT(i, j) = ScalarTraits<T>::Zero();
            QT(j, i) = ScalarTraits<T>::Zero();
        }
    }


    return res;
}

/**
 * shimit QR分解
 * 通过上三角矩阵Q对A进行施密特正交化
 * 这里就是ak减去自己在前面正交基上的投影
 */

template <CommutativeScalar T>
QRcomposition<T> QR_Shimit(const Matrix<T>& A)
{
    QRComposition<T> res;
    res.Q = Matrix<T>::Identity(A.Cows(), A.Cows());
    res.R = A;
    Matrix<T> Q = Matrix<T>::Identity(A.Cows(), A.Cows());
    for (int j = 0; j < A.Cols(); j++) {
        Matrix<T> vk = Matrix<T>(A.Rows(), 1);
        Matrix<T> vj = Matrix<T>(1, A.Rows());
        for (int k = 0; k < j; k++) {
            for (int i = 0; i < j; i++) {
                vk(i, 0) = res.R(i, k);
                vj(i, 0) = res.R(i, j);
            }
            Q(k, j) = Dot(vk, vj) / Dot(vj, vj);

            for (int i = 0; i < j; i++) {
                res.R(i, j) = res.R(i, j) - Q(k, j) * vk(i, 0);
            }
        }
        
        res.Q = Q * res.Q;

        for (int k = 0; k < j; k++) {
            Q(k, j) = ScalarTraits<T>::Zero();
        }
    }
    return res;
}