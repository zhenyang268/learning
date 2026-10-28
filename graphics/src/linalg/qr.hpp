#pragma once

// ===========================================================================
// linalg/qr.hpp —— QR 分解 (三种实现, 并列互不复用)
// ===========================================================================
//
// 【状态】已实现三条路径。
//
// ---------------------------------------------------------------------------
// 一、API 契约 (三条路径共用同一口径)
// ---------------------------------------------------------------------------
//   template <typename T> struct QRcomposition { Matrix<T> Q, R; };
//
//   QRcomposition<T> QR_Householder(const Matrix<T>& A);   // 反射 (Householder)
//   QRcomposition<T> QR_givens(const Matrix<T>& A);        // 旋转 (Givens)
//   QRcomposition<T> QR_Shimit(const Matrix<T>& A);        // 修正 Gram–Schmidt
//
//   (**thin / reduced** QR; 要求 m >= n)
//     Q: m×n, 列正交/酉 —— Q^H·Q = I_n
//     R: n×n, 上三角
//     满足 A = Q·R
//
// 异常:
//   · 空矩阵 / m < n  -> throw std::invalid_argument
//     (宽矩阵在 thin 约定下退化 —— **统一抛错, 不去猜另一种语义**)
//
// 不变量:
//   · Q^H·Q = I_n (注意是 **H 不是 T**, 见末尾)。
//   · R 的上三角之外恒为 0; R 对角元可正可负 (三条路径的符号约定各自不同)。
//
// ---------------------------------------------------------------------------
// 二、调用顺序 / 内部机制
// ---------------------------------------------------------------------------
//   三条路径**互不复用实现** —— 各自独立算, 这样才能互相做**对拍判据**。
//   (如果共用内核, 对拍就只是"同一个 bug 出现两次"。)
//   共同的小工具放在 `qr_detail`:
//     LeftCols / TopRows   取左上 n 列 / 上 n 行 (从完整 QR 里截 thin 部分)
//     FillOrthonormalColumn  线性相关时在正交补里补一个合法方向
//
//   三条路径的取舍:
//     Householder  最稳, 单位置消整列; 是默认首选
//     Givens       一次消一个元素, 适合稀疏/带状 (这里只作对照)
//     Shimit       (修正 Gram–Schmidt) 正交性最差, 但能显式给出每一步的投影
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑
// ---------------------------------------------------------------------------
// 【"shimit" 是 Shimit 的拼法 —— 指的其实是经典 Gram–Schmidt 的修正版】
//   CGS (classical) 与 MGS (modified) 的区别只在**用哪个 q 做投影**:
//     CGS: r_ij = q_i^H·a_j  (i,j 都用原向量) —— 数值正交性会丢
//     MGS: 边算边扣 (用刚更新过的 a_j)        —— 正交性好得多
//   这里走 MGS。判据: Q^H·Q 与 I 的偏差, CGS 在病态下能到 1e-8 量级。
//
// 【Q^H 而不是 Q^T —— 复数域陷阱, 已两次踩】
//   复数域下 "列正交" 说的是 Q^H·Q = I, 不是 Q^T·Q = I。
//   写 Transpose 在 double 下**完全看不出问题** (conj 恒等), 只有 Complex
//   域的用例会红。已踩两次: givens 旋转、qr_inverse。
//   ⇒ 本文件任何"正交性"相关的写法, 必查用的是 Adjoint。
//
// 【线性相关列怎么处理 —— 不 throw, 补一个正交补方向】
//   第 j 列与前 j-1 列线性相关时, 数学上 R(j,j) = 0, Q 的第 j 列**不唯一**
//   (正交补里任何单位向量都合法)。此时:
//     · R(j,j) = 0 (如实记录"这一列没有新信息")
//     · Q 的第 j 列由 FillOrthonormalColumn 造一个与已有列正交的单位向量
//   ⇒ 秩亏矩阵也能走通, 由调用方通过 R 的对角元判断秩。这比 throw 更有用。
// ===========================================================================

#include <cmath>
#include <stdexcept>

#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "linalg/operator.hpp"

// 统一口径 (三条路径共用, thin / reduced QR; 要求 m >= n):
//   Q: m×n, 列正交/酉 (Q^H·Q = I_n)
//   R: n×n, 上三角, A = Q·R
// 宽矩阵 (m < n) 在 thin 约定下退化 —— 三个入口统一抛错, 不去猜另一种语义。
template <CommutativeScalar T>
struct QRcomposition
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T>           Q;       // 正交/酉基, m×n, 按列存放
    Matrix<T>           R;       // 上三角, n×n
};

namespace qr_detail {

// 取前 c 列 (m×k -> m×c)
template <CommutativeScalar T>
Matrix<T> LeftCols(const Matrix<T>& A, int c)
{
    Matrix<T> res(A.Rows(), c);
    for (int i = 0; i < A.Rows(); i++)
        for (int j = 0; j < c; j++) res(i, j) = A(i, j);
    return res;
}

// 取前 r 行 (k×n -> r×n)
template <CommutativeScalar T>
Matrix<T> TopRows(const Matrix<T>& A, int r)
{
    Matrix<T> res(r, A.Cols());
    for (int i = 0; i < r; i++)
        for (int j = 0; j < A.Cols(); j++) res(i, j) = A(i, j);
    return res;
}

// 给 X 的第 col 列补一个"与前面 col 列都正交的单位向量"(秩亏时的兜底)。
// 做法就是拿标准基 e_c 逐个试, 减掉在已有列上的投影, 剩下还够长的那个。
// 不参与重构 (该列对应 R 的对角元为 0), 只为让 Q 保持列正交。
template <CommutativeScalar T>
void FillOrthonormalColumn(Matrix<T>& X, int col)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int m = X.Rows();
    const real_t tol = DefaultTol<T>(m);

    for (int cand = 0; cand < m; cand++) {
        Matrix<T> v = Matrix<T>::Zero(m, 1);
        v(cand, 0) = ScalarTraits<T>::one();
        for (int k = 0; k < col; k++) {
            T acc = ScalarTraits<T>::zero();
            for (int i = 0; i < m; i++)
                acc = acc + ScalarTraits<T>::conj(X(i, k)) * v(i, 0);
            for (int i = 0; i < m; i++)
                v(i, 0) = v(i, 0) - acc * X(i, k);
        }
        real_t n2 = real_t(0);
        for (int i = 0; i < m; i++) n2 = n2 + ScalarTraits<T>::abs2(v(i, 0));
        if (n2 > tol * tol) {
            const real_t nrm = std::sqrt(n2);
            for (int i = 0; i < m; i++)
                X(i, col) = v(i, 0) * ScalarTraits<T>::fromReal(real_t(1) / nrm);
            return;
        }
    }
    throw std::runtime_error("QR: 无法补出正交方向 (列数超过空间维数)");
}

} // namespace qr_detail

/**
 * householder QR分解
 * Hx = alpha * e1, e1 = [1, 0, 0, ..., 0]T
 * H = I - 2 * v * v^H / (v^H * v)
 * 取 v = x - alpha * e1, alpha = -phase(x1) * ||x||, 则 Hx = alpha * e1, 且
 *   v1 = x1 + phase(x1)*||x|| = phase(x1) * (|x1| + ||x||)   —— 是"加", 不相消
 * 每步只作用在 [j, m) 子块上 (x 的前 j 个分量本来就是 0)。
 */
template <CommutativeScalar T>
QRcomposition<T> QR_Householder(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int m = A.Rows();
    const int n = A.Cols();
    if (m == 0 || n == 0)  throw std::invalid_argument("QR_Householder: 空矩阵");
    if (n > m)             throw std::invalid_argument("QR_Householder: 需要 m >= n (thin QR)");

    Matrix<T> Qf = Matrix<T>::Identity(m);   // 完整 m×m 累积, 末尾截前 n 列
    Matrix<T> Rf = A;

    for (int j = 0; j < n; j++) {
        // x = Rf(j:m-1, j)
        real_t x2 = real_t(0);
        for (int i = j; i < m; i++) x2 = x2 + ScalarTraits<T>::abs2(Rf(i, j));
        const real_t xnorm = std::sqrt(x2);
        if (!(xnorm > DefaultTol<T>(m))) continue;   // 该列已在(数值)对角线上, 跳过

        Matrix<T> v = Matrix<T>::Zero(m, 1);
        for (int i = j; i < m; i++) v(i, 0) = Rf(i, j);

        // phase = x1/|x1| (x1 = 0 时取 1); alpha = -phase * ||x||
        const real_t a1 = std::sqrt(ScalarTraits<T>::abs2(v(j, 0)));
        const T phase = (a1 > real_t(0))
                      ? v(j, 0) * ScalarTraits<T>::fromReal(real_t(1) / a1)
                      : ScalarTraits<T>::one();
        v(j, 0) = v(j, 0) + phase * ScalarTraits<T>::fromReal(xnorm);

        // 分母必须是"改后"的 v^H·v = 2||x||(||x|| + |x1|)。
        // 原实现用改前的 ||x||² 当分母 —— H 不再是反射 (不对合), 消不掉 x。
        const real_t v2 = Abs2(v);
        if (!(v2 > real_t(0))) continue;
        const T coef = ScalarTraits<T>::fromReal(real_t(2) / v2);

        Matrix<T> H = Matrix<T>::Identity(m);
        for (int i = j; i < m; i++)
            for (int k = j; k < m; k++)
                H(i, k) = H(i, k) - coef * v(i, 0) * ScalarTraits<T>::conj(v(k, 0));

        Rf = H * Rf;
        // Q = Q·H^H。H 是 Hermitian (vv^H 自伴), 所以 H^H = H; 但复数域下
        // 必须写 Adjoint() —— Transpose() 只在实矩阵下恰好等价。
        Qf = Qf * H.Adjoint();
    }

    QRcomposition<T> res;
    res.Q = qr_detail::LeftCols(Qf, n);
    res.R = qr_detail::TopRows(Rf, n);
    return res;
}

/**
 * givens QR分解
 * 用一个只作用在第 j、i 两行上的平面旋转, 把 R(i, j) (i > j) 逐点消成 0。
 * G = [[ c,  s], [-conj(s), c]]  (c 为实数, |c|² + |s|² = 1) 是酉阵。
 * 对 (f, g) = (R(j,j), R(i,j)) 取
 *   r = hypot(|f|, |g|),  c = |f|/r,  phase = f/|f|,  s = phase * conj(g)/r
 * 则 G·[f; g] = [(f/|f|)·r; 0] —— 第 i 行归零, 第 j 行取得模长 r。
 * 复数域靠 phase 补自由度: 一个实数角 θ 只能消一个实数条件, 而 g 有实虚两部
 * (与 eigen 里 ph = z/|z| 同一件事)。
 */
template <CommutativeScalar T>
QRcomposition<T> QR_givens(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int m = A.Rows();
    const int n = A.Cols();
    if (m == 0 || n == 0) throw std::invalid_argument("QR_givens: 空矩阵");
    if (n > m)            throw std::invalid_argument("QR_givens: 需要 m >= n (thin QR)");

    Matrix<T> Qf = Matrix<T>::Identity(m);
    Matrix<T> Rf = A;
    const real_t tol2 = DefaultTol<T>(m) * DefaultTol<T>(m);

    // 外层列 j, 内层行 i > j —— 消的是下三角。
    // 原实现写成"外层行 i、内层列 j > i", 是它的镜像 (去消上三角), 而且 j 的上界用了
    // A.Rows(): R 只有 n 列, j 一超过 n 就越界 (ASan 报的 heap-buffer-overflow)。
    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < m; i++) {
            const T f = Rf(j, j);
            const T g = Rf(i, j);
            const real_t f2 = ScalarTraits<T>::abs2(f);
            const real_t g2 = ScalarTraits<T>::abs2(g);
            if (!(g2 > tol2)) continue;              // 已经是 0, 不必旋转

            const real_t r = std::sqrt(f2 + g2);
            const real_t af = std::sqrt(f2);
            const real_t c = af / r;                 // cos, 实数非负
            const T phase = (af > real_t(0))
                          ? f * ScalarTraits<T>::fromReal(real_t(1) / af)    // f/|f|
                          : ScalarTraits<T>::one();
            const T s = phase * ScalarTraits<T>::conj(g)
                      * ScalarTraits<T>::fromReal(real_t(1) / r);

            Matrix<T> G = Matrix<T>::Identity(m);
            G(j, j) = ScalarTraits<T>::fromReal(c);
            G(j, i) = s;
            G(i, j) = -ScalarTraits<T>::conj(s);
            G(i, i) = ScalarTraits<T>::fromReal(c);

            Rf = G * Rf;
            Qf = Qf * G.Adjoint();
        }
    }

    QRcomposition<T> res;
    res.Q = qr_detail::LeftCols(Qf, n);
    res.R = qr_detail::TopRows(Rf, n);
    return res;
}

/**
 * Gram-Schmidt QR分解 (原注称 "shimit", 即 Schmidt)
 * 对 A 的每一列做正交化: v_j = a_j - Σ_{k<j} r_kj · q_k, 然后 q_j = v_j / r_jj。
 * 采用**修正** GS (MGS): 每减掉一次投影, 就用**当前**的 v_j 去算下一个系数
 * (经典 GS 用原始 a_j 一次算完全部系数 —— 数值上更易失正交; MGS 是它的顺序改进版)。
 */
template <CommutativeScalar T>
QRcomposition<T> QR_Shimit(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int m = A.Rows();
    const int n = A.Cols();
    if (m == 0 || n == 0) throw std::invalid_argument("QR_Shimit: 空矩阵");
    if (n > m)            throw std::invalid_argument("QR_Shimit: 需要 m >= n (thin QR)");

    QRcomposition<T> res;
    res.Q = Matrix<T>::Zero(m, n);     // m×n: A 各列正交化后的单位基 (原实现填的是 n×n 系数阵)
    res.R = Matrix<T>::Zero(n, n);     // n×n: 投影系数 (上三角) + 对角范数
    Matrix<T> V = A;                   // 工作副本; 第 j 列在循环里被逐步正交化, 原地更新

    const real_t tol = DefaultTol<T>(m);

    for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
            // r_kj = <q_k, v_j> = Σ conj(q_ik)·v_ij   (内积共轭第一个参数, D7)
            // 这里不除 <q_k, q_k> —— Q 的每一列在**它自己那一轮**末尾已经归一化了,
            // 所以被引用的基向量恒为单位向量, <q_k,q_k> ≡ 1, 分母被归一化吸收掉了。
            // (若不归一化, 就得写 acc / <v_k, v_k>, 且每减一次投影都要现算一遍分母。)
            T acc = ScalarTraits<T>::zero();
            for (int i = 0; i < m; i++)
                acc = acc + ScalarTraits<T>::conj(res.Q(i, k)) * V(i, j);
            res.R(k, j) = acc;
            for (int i = 0; i < m; i++)
                V(i, j) = V(i, j) - acc * res.Q(i, k);
        }

        real_t n2 = real_t(0);
        for (int i = 0; i < m; i++) n2 = n2 + ScalarTraits<T>::abs2(V(i, j));
        const real_t nrm = std::sqrt(n2);

        if (nrm > tol) {
            res.R(j, j) = ScalarTraits<T>::fromReal(nrm);          // 归一化因子落在 R 对角
            for (int i = 0; i < m; i++)
                res.Q(i, j) = V(i, j) * ScalarTraits<T>::fromReal(real_t(1) / nrm);
        } else {
            // 第 j 列与前面线性相关: R(j,j) = 0, Q 的第 j 列在正交补里任取一个方向
            res.R(j, j) = ScalarTraits<T>::zero();
            qr_detail::FillOrthonormalColumn(res.Q, j);
        }
    }
    return res;
}