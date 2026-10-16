#pragma once

#include <stdexcept>
#include "core/common.hpp"      // MatrixKind 已定义于此, 直接复用不重复定义
#include "core/matrix.hpp"
#include "core/scalar.hpp"

#include "linalg/solve.hpp"   // LU / Solve / DefaultTol (S1 已落地, 2026-10-14)
//   LUResult{L,U,piv} / LU(A) / Solve(A,b) 解 A·x = b, 教科书列向量约定, D19 修订版
#include "linalg/cholesky.hpp" // cholesky / cholesky_solve (S2a 已落地, 2026-10-14)
//   SPD 判据 = "Cholesky 不选主元能否走通" (matrix.md §3.2) —— 本文件只做分发, 不重复实现分解
#include "linalg/qr.hpp"       // QR_Householder (S6b 已落地) —— qr_inverse 走 A^{-1} = R^{-1}·Q^H
//   qr.hpp 只 include core/* + operator.hpp, **不反向依赖本文件**, 故无循环。
//   (eigen.hpp 反过来 include 本文件, 于是顺序是 core -> solve/cholesky/qr -> inverse -> eigen。)

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

// SPD 判据 (matrix.md §3.2): "Cholesky 不选主元能否走通" —— 存在性定理 A 正定 <=> 分解走通。
// 只把 cholesky 抛出的 runtime_error (非正定 / 近奇异) 解读为"否"; 别的异常照常外传, 不吞。
template <CommutativeScalar T>
bool is_positive_definite(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        return false;   // SPD 蕴含方阵; 非方阵拿去问 cholesky 只会换回一个入参错, 直接答"否"
    try {
        (void)cholesky(A);
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}
// 口径: 近奇异(主元实部 <= 8n·eps)也算"否" —— 与 cholesky 同一把尺子, 不放近奇异矩阵过去。

template <CommutativeScalar T>
MatrixKind Classify(const Matrix<T>& A, typename ScalarTraits<T>::real_t tol)
{
    if (is_diagonal(A, tol)) return MatrixKind::Diagonal;
    if (is_orthogonal(A, tol)) return MatrixKind::Orthogonal;
    if (is_symmetric(A, tol)) {
        // SPD ⊂ Symmetric, 必须在这条**对称分支内部**筛子类。平铺成同级 if 的话,
        // 上一行先 return Symmetric 就把它吃掉了, 本档永远不可达 —— S2a-2 的坑。
        if (is_positive_definite(A)) return MatrixKind::SymmetricPositiveDefinite;
        return MatrixKind::Symmetric;
    }
    return MatrixKind::Generic;
}
// 探测顺序: Diagonal -> Orthogonal -> Symmetric{ -> SPD } -> Generic
// Diagonal 排在 SPD 之前是有意的: 对角矩阵即使正定也走 Diagonal(逐元倒数更省)。

// ---- 各分支专用求逆 (五个结构各一条路径; 近奇异/非正定在各自内部 throw, 不留静默 no-op) ----
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


// A = L·L^H, A^{-1} = (L^{-1})^H·L^{-1}   —— S2a-2 接入
// cholesky_solve(A, I) 展开正是 "L·Y = I 前代 + L^H·X = Y 回代", 与上式同一条路径,
// 于是直接复用那份压过 98 项验收的代码, 不另写一遍三角求逆 (少一份要维护的正确性)。
template <CommutativeScalar T>
Matrix<T> cholesky_inverse(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("cholesky_inverse: A 必须是非空方阵");
    return cholesky_solve(A, Matrix<T>::Identity(A.Rows()));   // A 非正定在此 throw (D12)
}

// A = Q·R (thin; 方阵时 Q、R 都是 n×n), A^{-1} = R^{-1}·Q^H。
// 只做两件**已有**的事: ① Householder QR (qr.hpp) ② 对上三角 R 的一次回代 (R·X = I ⇒ X = R^{-1})。
// 不写成 Solve(R, I): R 本身已经上三角, 再走一遍 LU 等于把 O(n³) 又花一遍;
// 也不显式求 R^{-1} 再去乘 —— 回代逐列解出 X 就是同一件事的紧凑写法 (与 qr.hpp 自带注释同口径)。
// 奇异判据沿用 lu_inverse / diagonal_inverse 的口径: |R_ii| <= tol ⇒ 抛 (D12)。
// 注意它**不接矩阵分发**: MatrixKind 里没有 QR 档 (也不该有 —— QR 是算法不是结构),
// 所以本函数与 lu_inverse / cholesky_inverse 并列, 是一条"显式点名"的独立路径。
template <CommutativeScalar T>
Matrix<T> qr_inverse(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("qr_inverse: A 必须是非空方阵");

    const int n = A.Rows();
    const QRcomposition<T> qr = QR_Householder(A);      // m >= n 由方阵保证
    const typename ScalarTraits<T>::real_t tol = DefaultTol<T>(n);

    Matrix<T> X = Matrix<T>::Identity(n);               // 右端项 B = I, 原地被解成 R^{-1}
    for (int i = n - 1; i >= 0; i--) {
        if (ScalarTraits<T>::abs2(qr.R(i, i)) <= tol * tol)
            throw std::runtime_error("qr_inverse: R 对角元近奇异 (秩亏或近奇异, 判据 |R_ii| <= tol)");
        const T inv_diag = ScalarTraits<T>::one() / qr.R(i, i);
        for (int j = 0; j < n; j++) {
            T acc = X(i, j);                            // 此刻第 i 行还是单位阵的 δ_ij, 当右端项用
            for (int k = i + 1; k < n; k++)
                acc = acc - qr.R(i, k) * X(k, j);       // k > i 的行已回代完, 可直接引用
            X(i, j) = acc * inv_diag;
        }
    }
    return X * qr.Q.Adjoint();                          // R^{-1}·Q^H
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
        case MatrixKind::SymmetricPositiveDefinite: return cholesky_inverse(A);
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