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

// ===========================================================================
// linalg/svd.hpp —— SVD 奇异值分解 (经济型)
// ===========================================================================
//
// 【状态】已实现 (经 A^H·A 的 Jacobi 特征分解路线)。
//
// ---------------------------------------------------------------------------
// 一、API 契约
// ---------------------------------------------------------------------------
//   template <typename T> struct SVDcomposition { Matrix<T> U, V, S; };
//   template <CommutativeScalar T>
//   SVDcomposition<T> svd(const Matrix<T>& A);
//
//   形状: A 为 m×n, r = min(m, n) —— **经济型**, 只保留 r 个奇异值。
//     U: m×r  左奇异向量按列 (列正交归一 / 酉)
//     S: r×r  对角, 对角元 σ_i >= 0 且**降序**
//     V: n×r  右奇异向量按列 (列正交归一 / 酉)
//   满足 A = U·S·V^H。
//
// 异常: 空矩阵 -> throw std::invalid_argument。
//
// 不变量:
//   · σ_i 非负且降序 (降序是**约定**, 不是数学必然 —— 标准库不保证, 这里做了)。
//   · U^H·U = I_r, V^H·V = I_r。
//
// 已知精度边界 (**必读, 否则会得出"实现错了"的错误结论**):
//   · **重构 A = U·S·V^H 精度良好** (≈ eps 级)。
//   · **单独取 σ 的精度只有 √eps 级** —— 因为走 A^H·A 把条件数平方了。
//     κ(A)=1e8 时 σ 的相对误差约 1e-8, 这是**方法决定的**, 不是 bug。
//   · 需要 σ 本身高精度 -> 必须换 Golub–Kahan 双对角化 (未做)。
//
// 秩判定: 不要用 "σ_i > 0"。用 σ_i > κ_max·σ_0 (相对阈值), 或
//   σ_i > DefaultTol<T>(n)·σ_0。见 numeric-contract.md §4.3。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   B = A^H·A  (n×n Hermitian 半正定)  = V·Λ·V^H   ← 复用 eig_symmetric
//   σ_i = √λ_i        (非负根, 再降序排序)
//   u_i = A·v_i / σ_i  (σ_i > tol)
//   σ_i ≈ 0 时: 左奇异向量在**正交补**里任取 (见 FillOrthonormal 一类做法),
//                因为此时 A 在这个方向上没有信息, 任何正交单位向量都合法。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【为什么 σ 的精度不拖累重构 —— 这是最容易被误判为 bug 的一条】
//     U·S·V^H = Σ_i (A·v_i/σ_i)·σ_i·v_i^H = A·(V·V^H)
//   σ 与 1/σ **精确相消** —— 重构只取决于 V 的列正交归一、且张成 A 的行空间。
//   所以 σ 只有 √eps 也不影响 A ≈ U·S·V^H 的精度。
//   ⚠ 但如果你单独拿 σ 去和解析值比, 会看到 √eps 级的差 —— 那不是实现错。
//
// 【A^H·A 把条件数平方 —— 为什么还是选了它】
//   这是"能用 + 便宜"对"精度好"的取舍。Golub–Kahan 双对角化精度好但实现
//   复杂得多。当前阶段的消费点 (秩判定、最小二乘、伪逆) 都能容忍 √eps,
//   而它们全都**依赖重构而非 σ 本身**。判据: 如果将来出现一个真的要 σ 高精度
//   的场景 (例如病态矩阵的谱条件数估计), 那时候再换, 现在不提前优化。
//
// 【降序排序是约定, 不是数学】
//   特征分解的返回顺序取决于算法内部细节 (Jacobi 的旋转顺序)。所以排序是
//   本文件**必须自己做的**一步, 且要与 V 的列同步置换 —— 只排 σ 不排 V
//   会让 U/S/V 三者错位, 得到"数值都对但配错了"的结果。
// ===========================================================================

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


