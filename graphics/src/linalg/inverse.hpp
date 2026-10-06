#pragma once

#include <stdexcept>
#include "core/common.hpp"      // MatrixKind 已定义于此, 直接复用不重复定义
#include "core/matrix.hpp"
#include "core/scalar.hpp"

// #include "linalg/solve.hpp"  // TODO: LU 落地后打开 (matrix.md §3.1:
//   LUResult{L,U,piv} / LU(A) / Solve(A,b) 解 X·A = b, 行向量约定 D19)

// ---- 结构探测: 每个判断是一次数学练习 (matrix.md §3.2) ----
// D7: 比较用 abs2, 共轭用 conj, 算法层零分支; D12: 容差 tol 由调用方给
template <CommutativeScalar T>
bool is_diagonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol);

template <CommutativeScalar T>
bool is_symmetric(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol);
// A ≈ A^H: 实数域=对称, 复数域=Hermitian (traits 写法免费升级)

template <CommutativeScalar T>
bool is_orthogonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol);
// A·A^H ≈ I: 实数域=正交, 复数域=酉 (traits 写法免费升级)

template <CommutativeScalar T>
bool is_positive_definite(const Matrix<T>& A);
// 判据: Cholesky 不选主元能否走通

template <CommutativeScalar T>
MatrixKind Classify(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol);
// 探测顺序: Diagonal -> Orthogonal -> SymmetricPositiveDefinite -> Generic

// ---- 各分支专用求逆 (实现体 TODO, 用户填写) ----
template <CommutativeScalar T>
Matrix<T> lu_inverse(const Matrix<T>& A);
// 阶段 1 主路径: LU 部分选主元 + 逐"行"回代 (D19 镜像!)
// A·A^{-1} = I => A^{-1} 第 i 行 r_i 满足 r_i·A = e_i => r_i = Solve(A, e_i)
// 教科书"逐列回代"是列向量语境, 本库行向量语境整体转置 = 逐行回代
// 奇异判定: 选主元 |pivot|^2 <= tol 则抛 (D12)

template <CommutativeScalar T>
Matrix<T> diagonal_inverse(const Matrix<T>& A);   // 对角元倒数, 零/近零抛 (D12)

template <CommutativeScalar T>
Matrix<T> cholesky_inverse(const Matrix<T>& A);   // A = L·L^H, A^{-1} = (L^{-1})^H·L^{-1}

template <CommutativeScalar T>
Matrix<T> qr_inverse(const Matrix<T>& A);         // A = Q·R, A^{-1} = R^{-1}·Q^H

// ---- 对外入口 (matrix.md §3.2 枚举分发) ----
template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A);            // Classify -> switch 分发
// Orthogonal 分支直接 A.Adjoint() (R^{-1} = R^H, 无需专用函数)

template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A, MatrixKind hint);  // 已知结构跳过探测