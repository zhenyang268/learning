#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"
#include "linalg/eigen.hpp"   // eig_symmetric: A^H·A 是 Hermitian, 直接走这一支

// ---------------------------------------------------------------------------
// SVD 奇异值分解 (plan.md 阶段 1 §3.3 / matrix.md §3.3)
//
//   A = U·Σ·V^H,   A 为 m×n,  r = min(m, n)   —— 经济型 (只保留 r 个奇异值)
//     U: m×r  左奇异向量按列 (列正交归一 / 酉)
//     Σ: r×r  对角, 对角元 σ_i >= 0 且降序
//     V: n×r  右奇异向量按列 (列正交归一 / 酉)
//
// 路线: 经 A^H·A 的 Jacobi 特征分解组装 (复用 S3 验完的 eig_symmetric)。
//   B = A^H·A (n×n Hermitian 半正定) = V·Λ·V^H,  σ_i = √λ_i (非负根 + 降序),
//   u_i = A·v_i / σ_i (σ_i > tol);  σ_i ≈ 0 时左奇异向量在正交补里任取。
//
// 为什么 σ 的精度不拖累重构: U·Σ·V^H = Σ_i (A·v_i/σ_i)·σ_i·v_i^H = A·(V·V^H),
//   σ 与 1/σ 精确相消 —— 重构只取决于 V 的列正交归一、且张成 A 的行空间。
//   (单独要 σ 的高精度则不行: A^H·A 把条件数平方, σ 只有 √eps 级 —— 那需阶段 2 的
//    Golub–Kahan 双对角化, plan.md 明确本阶段走 A^H·A 路线。)
// ---------------------------------------------------------------------------

template <typename T>
struct SVDcomposition
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T> U;   // m×r, 左奇异向量按列
    Matrix<T> V;   // n×r, 右奇异向量按列
    Matrix<T> S;   // r×r, 对角矩阵, 对角元为实 (奇异值非负且降序)
};

namespace linalg_detail {

// 给 U 的第 col 列补一个与前面 col 列正交归一的单位列。
// 只在 σ ≈ 0 时用到: u = A·v/σ 数值不可用, 而零奇异值对应的左奇异向量本就在
// "与 A 值域正交补"里任取。取一个确定性结果: 依次扫标准基 e_0..e_{m-1},
// 对已有列做 Gram-Schmidt, 挑投影后仍够长的那个 (即与已有列线性无关)。
template <CommutativeScalar T>
inline void FillOrthonormalColumn(Matrix<T>& U, int col)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int m = U.Rows();
    for (int cand = 0; cand < m; cand++) {
        std::vector<T> w(m, ScalarTraits<T>::zero());
        w[cand] = ScalarTraits<T>::one();
        for (int c = 0; c < col; c++) {
            T dot = ScalarTraits<T>::zero();
            for (int p = 0; p < m; p++) dot = dot + ScalarTraits<T>::conj(U(p, c)) * w[p];
            for (int p = 0; p < m; p++) w[p] = w[p] - dot * U(p, c);
        }
        real_t nrm2 = real_t(0);
        for (int p = 0; p < m; p++) nrm2 = nrm2 + ScalarTraits<T>::abs2(w[p]);
        if (nrm2 > real_t(0.25)) {                       // 投影后仍够长 -> 与已有列无关
            const T scale = ScalarTraits<T>::fromReal(real_t(1) / std::sqrt(nrm2));
            for (int p = 0; p < m; p++) U(p, col) = w[p] * scale;
            return;
        }
    }
    throw std::runtime_error("svd: 无法为正交补找到基向量 (内部错误)");
}

} // namespace linalg_detail

// ---- 主入口: 经济型 SVD (A = U·Σ·V^H) ----
template <CommutativeScalar T>
SVDcomposition<T> svd(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    const int m = A.Rows();
    const int n = A.Cols();
    if (m == 0 || n == 0)
        throw std::invalid_argument("svd: A 不能为空矩阵");

    const int r = std::min(m, n);

    // ---- ① B = A^H·A (n×n Hermitian PSD), Jacobi 特征分解 ----
    // B 按构造即 Hermitian, 直接走 eig_symmetric (不走 eig 的 is_symmetric 探测)。
    const Matrix<T> B = A.Adjoint() * A;
    const Eigendecomposition<T> e = eig_symmetric(B);

    // 取对角元 λ_i (理论上 >= 0; 浮点下可能出 -1e-17) 与对应特征向量列
    std::vector<real_t> lam(n);
    for (int i = 0; i < n; i++) lam[i] = ScalarTraits<T>::real(e.Lamda(i, i));

    // 降序排序 (V1.4 要求"奇异值非负有序"), 同步记录特征向量的列号置换
    std::vector<int> ord(n);
    for (int i = 0; i < n; i++) ord[i] = i;
    std::sort(ord.begin(), ord.end(),
              [&](int a, int b) { return lam[a] > lam[b]; });

    // ---- ② σ_i = √max(λ_i, 0), 只取前 r 个 ----
    std::vector<real_t> sig(r);
    for (int i = 0; i < r; i++) {
        const real_t li = lam[ord[i]];
        sig[i] = std::sqrt(li > real_t(0) ? li : real_t(0));
    }
    const real_t tol_sig = DefaultTol<T>(r) * sig[0];    // σ 尺度上的零判据 (D12)

    // ---- ③ 组装 U / Σ / V ----
    SVDcomposition<T> res;
    res.U = Matrix<T>::Zero(m, r);
    res.S = Matrix<T>::Zero(r, r);
    res.V = Matrix<T>::Zero(n, r);

    for (int i = 0; i < r; i++) {
        const int col = ord[i];

        for (int k = 0; k < n; k++) res.V(k, i) = e.Q(k, col);   // 右奇异向量 = B 的特征向量
        res.S(i, i) = ScalarTraits<T>::fromReal(sig[i]);

        if (sig[i] > tol_sig) {
            // u_i = A·v_i / σ_i (逐行做 A 的第 p 行与 v_i 的内积)
            const T inv = ScalarTraits<T>::fromReal(real_t(1) / sig[i]);
            for (int p = 0; p < m; p++) {
                T acc = ScalarTraits<T>::zero();
                for (int k = 0; k < n; k++) acc = acc + A(p, k) * res.V(k, i);
                res.U(p, i) = acc * inv;
            }
        } else {
            linalg_detail::FillOrthonormalColumn(res.U, i);      // σ ≈ 0: 左向量在正交补里取
        }
    }

    return res;
}


