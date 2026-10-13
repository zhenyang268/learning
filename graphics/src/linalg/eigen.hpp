#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"
#include "linalg/inverse.hpp"   // is_symmetric(A, tol) + DefaultTol<T>(n) (S2/S2a-2 已落地)

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
    Matrix<T>           Lamda;   // 对角矩阵, 对角元为实 (Hermitian 的结构决定, 与标量无关)
};

// ---- ① 主入口: Q 与 Λ ----
// 仅限对称/Hermitian (plan.md 阶段 1 范围); 非对称矩阵特征值可能为复数, 属阶段 2 (matrix.md §3.3)。
template <CommutativeScalar T>
Eigendecomposition<T> eig(const Matrix<T>& A)
{
    if (A.Rows() != A.Cols() || A.Rows() == 0)
        throw std::invalid_argument("eig: A 必须是非空方阵");
    if (is_symmetric(A, DefaultTol<T>(A.Rows())))   // 实对称与 Hermitian 统一走这一支
        return eig_symmetric(A);
    throw std::logic_error(
        "eig: 非对称/Hermitian 矩阵的特征值可能为复数, 需阶段 2 的 Hessenberg+QR (matrix.md §3.3)");
}

// ---- ② 2x2 Hermitian 子块的一次 Jacobi 旋转 ----
// 子块 A2 = [[a, z], [conj(z), b]] (Hermitian 对角元必为实, 故 a, b 走 real_t)。
// 取酉阵 U = [[c, -s·ph], [s·conj(ph), c]], c = cosθ, s = sinθ, ph = z/|z| = e^{iφ}。
// 直接验算 U^H·A2·U 的 (p,q) 元 = e^{iφ}·(ρ·cos2θ − (a−b)·sin2θ/2) (ρ = |z|), 故取零当且仅当
//     tan(2θ) = 2|z| / (a − b)      即     θ = 0.5·atan2(2|z|, a − b)
// 复数域比实对称只多一个相位: 实对称下 z 为实数, ph = ±1, U 退化成普通实旋转。
// 旋转后的两个对角元正是该子块的特征值, 写成解析式 (比 c²a+s²b+2cs|z| 少一次相消):
//     a' = (a+b)/2 + D/2,   b' = (a+b)/2 − D/2,   D = hypot(a−b, 2|z|)
// 注意 a', b' 与上面 θ 的取法**必须配套**: 取 θ+π/2 同样能消零, 但会把两个特征值对调。
template <CommutativeScalar T>
void jacobi_rotation(typename ScalarTraits<T>::real_t a,
                     typename ScalarTraits<T>::real_t b,
                     const T& z,
                     typename ScalarTraits<T>::real_t& c,
                     typename ScalarTraits<T>::real_t& s,
                     T& ph,
                     typename ScalarTraits<T>::real_t& a_new,
                     typename ScalarTraits<T>::real_t& b_new)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const real_t rho = std::sqrt(ScalarTraits<T>::abs2(z));       // |z| >= 0
    const real_t th  = real_t(0.5) * std::atan2(real_t(2) * rho, a - b);
    c = std::cos(th);
    s = std::sin(th);
    ph = (rho > real_t(0)) ? z / ScalarTraits<T>::fromReal(rho)
                           : ScalarTraits<T>::one();              // e^{iφ}
    const real_t half = (a + b) * real_t(0.5);
    const real_t d    = std::hypot(a - b, real_t(2) * rho) * real_t(0.5);
    a_new = half + d;
    b_new = half - d;
}

// ---- ③ 对称/Hermitian 的特征分解: cyclic Jacobi 迭代 ----
// 为什么必须**循环**: 一次旋转只在精确算术下把 A(p,q) 置零, 且后续旋转会往 (p,q) 位置
// 重新灌入非零 (fill-in) —— 所以得一轮轮扫 (sweep) 全部 i<j, 直到整个非对角块足够小。
// 收敛性: 每轮 off(A)² = Σ_{i≠j}|A(i,j)|² 单调下降, 渐近二次收敛;
// 实对称与 Hermitian 共用这一份代码 (V1.1)。
// 容差标尺: 酉相似变换保持 Frobenius 范数, 所以 ‖A‖_F 只需开头算一次 —— 容差相对它取 (D12)。
template <CommutativeScalar T>
Eigendecomposition<T> eig_symmetric(const Matrix<T>& A)
{
    using real_t = typename ScalarTraits<T>::real_t;

    if (A.Rows() != A.Cols())
        throw std::invalid_argument("eig_symmetric: A 必须是方阵");
    const int n = A.Rows();
    if (n == 0)
        throw std::invalid_argument("eig_symmetric: A 不能为空矩阵");

    Eigendecomposition<T> res;
    res.Q     = Matrix<T>::Identity(n);
    res.Lamda = A;                        // 原地相似变换成对角形
    Matrix<T>& L = res.Lamda;

    real_t frob2 = real_t(0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            frob2 = frob2 + ScalarTraits<T>::abs2(A(i, j));
    if (frob2 == real_t(0)) return res;                  // 零矩阵: 已经是对角形

    const real_t frob  = std::sqrt(frob2);
    const real_t tol   = DefaultTol<T>(n) * frob;        // 总容差: off(A) <= tol
    const real_t tol2  = tol * tol;
    const int    nd    = n * (n - 1);                    // 非对角元个数
    const real_t skip2 = (nd > 0) ? tol2 / real_t(nd) : real_t(0);
    // 单元素跳过阈值取 tol²/nd: 即使剩余元素全被跳过, 总和也不超过 tol²。

    const int kMaxSweep = 100;
    for (int sweep = 0; sweep < kMaxSweep; sweep++) {
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                if (ScalarTraits<T>::abs2(L(p, q)) <= skip2) continue;   // 已够小, 转了只是噪声

                real_t c, s, na, nb;
                T ph;
                jacobi_rotation<T>(ScalarTraits<T>::real(L(p, p)),
                                   ScalarTraits<T>::real(L(q, q)),
                                   L(p, q), c, s, ph, na, nb);
                // real_t 与 T 不能直接相乘: Complex 没有单参构造, real_t*Complex 不会隐式转换,
                // 一律先 fromReal() 折成 T 再乘 (D15 的 traits 化)。
                const T tc = ScalarTraits<T>::fromReal(c);
                const T sc = ScalarTraits<T>::fromReal(s) * ScalarTraits<T>::conj(ph);  // s·e^{-iφ}
                const T sp = ScalarTraits<T>::fromReal(s) * ph;                         // s·e^{iφ}

                // (a) 两行两列之外: (k,p)/(k,q) 用**旧值**同步更新, 再按 Hermitian 补齐对角侧
                for (int k = 0; k < n; k++) {
                    if (k == p || k == q) continue;
                    const T kp = L(k, p);
                    const T kq = L(k, q);
                    const T np = tc * kp + sc * kq;
                    const T nq = tc * kq - sp * kp;
                    L(k, p) = np;
                    L(k, q) = nq;
                    L(p, k) = ScalarTraits<T>::conj(np);
                    L(q, k) = ScalarTraits<T>::conj(nq);
                }

                // (b) 对角 2x2 用解析式, 并把 (p,q) 精确置零 (理论上已为 0)
                L(p, p) = ScalarTraits<T>::fromReal(na);
                L(q, q) = ScalarTraits<T>::fromReal(nb);
                L(p, q) = ScalarTraits<T>::zero();
                L(q, p) = ScalarTraits<T>::zero();

                // (c) 特征向量累积 Q <- Q·U, 只动第 p、q 两列 (U 仅在 (p,q) 子块非平凡)
                for (int k = 0; k < n; k++) {
                    const T kp = res.Q(k, p);
                    const T kq = res.Q(k, q);
                    res.Q(k, p) = tc * kp + sc * kq;
                    res.Q(k, q) = tc * kq - sp * kp;
                }
            }
        }

        // 收敛检查: off(A)² <= tol²
        real_t off2 = real_t(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j) off2 = off2 + ScalarTraits<T>::abs2(L(i, j));
        if (off2 <= tol2) return res;
    }

    throw std::runtime_error("eig_symmetric: Jacobi 迭代 " + std::to_string(kMaxSweep)
                             + " 轮仍未收敛 (A 可能不是 Hermitian?)");
}


// ---- ④ 预留: 阶段 2 (非对称可能有复特征值), 只声明不实现 ----
// Eigendecomposition<T> eig(const Matrix<T>& A, General);
// std::vector<Complex<real_t>> eigenvalues(const Matrix<T>& A, General);