#pragma once

#include <vector>
#include <stdexcept>
#include <functional>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "autodiff/dual.hpp"

// DualFunc: 计划中的逐元素函数类型 (matrix.md §7.1 第 5 条, 目标签名 T(const T&))
// 当前函数节点走 FunctionType 枚举, DualFunc 暂无使用者 —— 别名占位
template <typename T>
using DualFunc = std::function<T(const T&)>;

// 对矩阵所有元素应用同一函数 (Function 节点前向)
template <typename T>
Matrix<T> Apply(const Matrix<T>& A, FunctionType f)
{
    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = dualFunc(A.Data()[i], f);
    }
    return B;
}

// 逐元素应用函数矩阵 (MatrixFunction 节点前向): 每格用自己的函数
template <typename T>
Matrix<T> Apply(const Matrix<T>& A, const Matrix<FunctionType>& f)
{
    if (f.Rows() == 1 && f.Cols() == 1) {
        return Apply(A, f(0, 0)); // 单函数退化为统一应用
    }

    // 应该只有一维, 但是后续计算图节点打包可能涉及多行
    if (A.Cols() != f.Cols()) {
        throw std::invalid_argument("Matrix apply dimension mismatch!");
    }

    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Size(); i++) {
        B.Data()[i] = dualFunc(A.Data()[i], f.Data()[i]);
    }
    return B;
}

// 把 value 赋值为 derive, 用于 backward (Function 节点: 局部导数暂存 derive)
template <typename T>
Matrix<T>& Derive(Matrix<T>& A)
{
    for (int i = 0; i < A.Size(); i++) {
        A.Data()[i].value = A.Data()[i].derive;
    }
    return A;
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
