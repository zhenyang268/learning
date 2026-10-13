#pragma once

#include <stdexcept>
#include "core/common.hpp"      // MatrixKind 已定义于此, 直接复用不重复定义
#include "core/matrix.hpp"
#include "core/scalar.hpp"

#include "linalg/solve.hpp"   // LU / Solve / DefaultTol (S1 已落地, 2026-10-14)
//   LUResult{L,U,piv} / LU(A) / Solve(A,b) 解 A·x = b, 教科书列向量约定, D19 修订版

// ---- 结构探测: 每个判断是一次数学练习 (matrix.md §3.2) ----
// D7: 比较用 abs2, 共轭用 conj, 算法层零分支; D12: 容差 tol 由调用方给
template <CommutativeScalar T>
bool is_diagonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < A.Cols(); j++) {
            if (i == j) continue;
            if (ScalarTraits<T>::abs2(A(i, j)) > tol * tol) return false;
        }
    }
    return true;
}

template <CommutativeScalar T>
bool is_symmetric(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    Matrix<T> AH = A.Adjoint();          // AH(i,j) = conj(A(j,i))
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < i; j++) {
            // 对称判据 A ≈ A^H 即 A(i,j) ≈ AH(i,j) (= conj(A(j,i)))。
            // 原写 AH(j,i) = conj(A(i,j)), 于是比的是"A(i,j) 的虚部是否为 0",
            // 实数矩阵恒真 -> 所有实方阵都会被判成 Symmetric。
            if (ScalarTraits<T>::abs2(A(i, j) - AH(i, j)) > tol * tol) return false;
        }
    }
    return true;
}
// 酉对称, Hermitian 对称, 实对称是免费自带的, 不用if
// A ≈ A^H: 实数域=对称, 复数域=Hermitian (traits 写法免费升级)

template <CommutativeScalar T>
bool is_orthogonal(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) return false;
    Matrix<T> G = A.Adjoint() * A;       // 正交/酉 <=> A^H·A = I
    for (int i = 0; i < A.Rows(); i++) {
        for (int j = 0; j < A.Cols(); j++) {
            // 只查对角线会漏检非对角元: [[1,x],[x,1]] 也会被判成正交
            const T expect = (i == j) ? ScalarTraits<T>::one() : ScalarTraits<T>::zero();
            if (ScalarTraits<T>::abs2(G(i, j) - expect) > tol * tol) return false;
        }
    }
    return true;
}
// A·A^H ≈ I: 实数域=正交, 复数域=酉 (traits 写法免费升级)

template <CommutativeScalar T>
bool is_positive_definite(const Matrix<T>& A)
{
    // 判据本应是"Cholesky 不选主元能否走通", 依赖 cholesky.hpp (plan.md S2a)
    (void)A;
    throw std::logic_error("is_positive_definite: 未实现 (依赖 cholesky.hpp, 见 plan.md S2a)");
}

template <CommutativeScalar T>
MatrixKind Classify(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (is_diagonal(A, tol)) return MatrixKind::Diagonal;
    if (is_orthogonal(A, tol)) return MatrixKind::Orthogonal;
    if (is_symmetric(A, tol)) return MatrixKind::Symmetric;
    // SPD 档待 cholesky 接入 (plan.md S2a) 后再打开。暂缺该档时对称矩阵落到 Generic,
    // 走 LU 求逆, 结果依然正确, 只是没用上 Cholesky 的优化。
    // if (is_positive_definite(A)) return MatrixKind::SymmetricPositiveDefinite;
    return MatrixKind::Generic;
}
// 探测顺序: Diagonal -> Orthogonal -> Symmetric -> (SPD, 待 cholesky 接入) -> Generic

// ---- 各分支专用求逆 (未接入 cholesky/qr 的分支显式 throw, 不留静默 no-op) ----
template <CommutativeScalar T>
Matrix<T> lu_inverse(const Matrix<T>& A)
{
    Matrix<T> I = Matrix<T>::Identity(A.Rows());
    return Solve(A, I);
}
// 奇异判定: 选主元 |pivot|^2 <= tol 则抛 (D12)


template <CommutativeScalar T>
Matrix<T> diagonal_inverse(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (A.Rows() != A.Cols()) throw std::invalid_argument("diagonal_inverse: A 必须是方阵");
    Matrix<T> B(A.Rows(), A.Cols());
    for (int i = 0; i < A.Rows(); i++) {
        if (ScalarTraits<T>::abs2(A(i, i)) <= tol * tol)   // 判据 |d| <= tol (D12), 比较走 abs2
            throw std::runtime_error("diagonal_inverse: 对角元近奇异 (判据 |d| <= tol)");
        B(i, i) = ScalarTraits<T>::one() / A(i, i);
    }
    return B;
}
// 对角元倒数, 零/近零抛 (D12)


// A = L·L^H, A^{-1} = (L^{-1})^H·L^{-1}   —— 依赖 cholesky.hpp (plan.md S2a)
template <CommutativeScalar T>
Matrix<T> cholesky_inverse(const Matrix<T>& A)
{
    (void)A;
    throw std::logic_error("cholesky_inverse: 未实现 (依赖 cholesky.hpp, 见 plan.md S2a)");
}

// A = Q·R, A^{-1} = R^{-1}·Q^H   —— 依赖 qr.hpp (plan.md S6)
template <CommutativeScalar T>
Matrix<T> qr_inverse(const Matrix<T>& A)
{
    (void)A;
    throw std::logic_error("qr_inverse: 未实现 (依赖 qr.hpp, 见 plan.md S6)");
}

// ---- 对外入口 (matrix.md §3.2 枚举分发) ----

// 已知结构: 跳过探测直接走该结构对应的实现 (Orthogonal 直接 A.Adjoint(), A^{-1} = A^H)。
template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A, MatrixKind hint)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("Inverse: A 必须是非空方阵");

    switch (hint) {
        case MatrixKind::Diagonal:                  return diagonal_inverse(A, DefaultTol<T>(A.Rows()));
        case MatrixKind::Orthogonal:                return A.Adjoint();
        case MatrixKind::Symmetric:
        case MatrixKind::Generic:                   return lu_inverse(A);
        case MatrixKind::SymmetricPositiveDefinite: return cholesky_inverse(A);   // 未实现 -> throw
        case MatrixKind::RigidTransform:            throw std::logic_error("Inverse: RigidTransform 分支未实现");
        case MatrixKind::Projection:                throw std::logic_error("Inverse: 投影矩阵不可逆 (A^2 = A != I)");
    }
    throw std::logic_error("Inverse: 未覆盖的 MatrixKind");
}

// 通用入口: 先 Classify 探测结构, 再分发。
template <CommutativeScalar T>
Matrix<T> Inverse(const Matrix<T>& A)
{
    return Inverse(A, Classify(A, DefaultTol<T>(A.Rows())));
}