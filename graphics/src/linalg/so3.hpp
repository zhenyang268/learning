#pragma once

// ===========================================================================
// linalg/so3.hpp —— 三维向量算子与 so(3) 李代数 (fan-plan A4 李群线的载体)
// ===========================================================================
//
// 为什么单开一个文件 (三层否决):
//   · 不进 core/  —— core 是「标量代数 + 哑容器」(complex / scalar / traits /
//                    matrix)。叉积是**几何**运算, 不是代数基础设施。且 §3.7 刚
//                    立过规矩: 不为一次性用途扩 core 接口 (与 Rotate 就地写
//                    trace 同理)。
//   · 不进 quat/  —— D16 把四元数封死在本目录且不依赖 linalg; 放这儿会让
//                    so(3) 那条线反向依赖 quat, 层次就乱了。
//   · 不塞 operator.hpp —— 那文件自己写了「Dot/Abs2 是暂放, S6 派支线再定
//                    去留」, 别往一个待搬迁的抽屉里再塞东西。
//
// 本文件只 include core/matrix.hpp, **不碰任何分解/求解** —— 所以它是 linalg
// 里最底的一层, 谁都能用。
//
// ---------------------------------------------------------------------------
// 三个入口, 一条链
// ---------------------------------------------------------------------------
//   Cross(a, b) : R^3 x R^3 -> R^3   叉积 (三个分量直算, 不建矩阵)
//   Hat(u)      : R^3      -> so(3)  3x1 列向量 -> 3x3 反对称
//   Vee(K)      : so(3)    -> R^3    3x3 反对称   -> 3x1 列向量
//
//   Vee(Hat(u)) = u  (恒等);  Cross(a, b) = Hat(a) * b = Vee(Hat(a)Hat(b))
//
// ---------------------------------------------------------------------------
// 口径
// ---------------------------------------------------------------------------
//   · 入参出参一律 **3x1 列向量** (D19 第二次修订, 2026-10-19: 全库列向量 +
//     分子布局)。形状不对 -> throw std::invalid_argument, 不静默接受
//     (Apply 那次越界读的教训)。
//   · **不造固定 3 维的 Vector3 类** —— 3x1 就是 Matrix<T>(3,1)。造个
//     Vector3 等于重新引入双表示。
//     口径澄清 (2026-10-23): matrix.md 决策 #21 立项的是**运行期长度**的
//     Vector<T>, 不是固定 3 的 Vector3 —— 那条与本节不冲突。
//     原则同源: **容器保持哑的, 能力做成薄层自由函数**
//     (matrix.md §6 第三次订正: 不造"胖" Tensor 类; 容器放行, 运算仍进自由函数)。
//   · Hat **不归一化**: hat(λu) = λ·hat(u) 是线性的, 归一化会破坏线性,
//       而 Rodrigues 需要 θ 单独乘在外面。要归一化的是
//       Quaternion::FromAxisAngle (那里是半角, 归一化是语义的一部分)。
//   · Vee 只读三个位置 (2,1) (0,2) (1,0) —— 与 Python hat/vee 同口径
//     (demo/pyapps/rotation/_rot_core.py:71 / :79), 保证数值可对拍。
//     **不做**"入参是否反对称"的检查: 那要算 K+Kᵀ 的范数, 是 O(n²) 还带
//     阈值选取问题 (§4 的 √eps 纪律)。契约写清即可。
// ===========================================================================

#include <stdexcept>

#include "core/matrix.hpp"

// ---- Cross: 叉积 ----------------------------------------------------------
//
//   a x b = (a_y b_z - a_z b_y,  a_z b_x - a_x b_z,  a_x b_y - a_y b_x)
//
// 右手系, e_x x e_y = e_z —— 与 Hamilton 约定 ij = k 一致 (quaternion.hpp:110)。
// 写成自由函数而不是 Matrix 成员: 叉积是 **R^3 专属** 的, 挂到通用的
// Matrix 上等于给所有形状开一个只在 3 维有意义的门。
template <typename T>
Matrix<T> Cross(const Matrix<T>& a, const Matrix<T>& b)
{
    if (a.Rows() != 3 || a.Cols() != 1 || b.Rows() != 3 || b.Cols() != 1) {
        throw std::invalid_argument(
            "Cross: 入参必须是两个 3x1 列向量 (D19 全库列向量 + 分子布局)");
    }
    const T ax = a(0, 0), ay = a(1, 0), az = a(2, 0);
    const T bx = b(0, 0), by = b(1, 0), bz = b(2, 0);

    Matrix<T> r(3, 1);
    r(0, 0) = ay * bz - az * by;
    r(1, 0) = az * bx - ax * bz;
    r(2, 0) = ax * by - ay * bx;
    return r;
}

// ---- Hat: R^3 -> so(3) ----------------------------------------------------
//
//   Hat(u) = [[0, -u_z, u_y],
//             [u_z, 0, -u_x],
//             [-u_y, u_x, 0]]        满足  Hat(u) b = u x b
//
// so(3) = 3x3 实反对称矩阵, 是 SO(3) 在单位元处的切空间 —— 李代数的具体载体。
// Rodrigues 闭式 R = exp(θ Hat(u)) 成立的**唯一关键**是 Hat(u)^3 = -|u|^2 Hat(u)
// (单位轴下退化为 -Hat(u)), 它把幂次降回 1/2 次: 见 report.md §1, 实测 2.2e-16。
template <typename T>
Matrix<T> Hat(const Matrix<T>& u)
{
    if (u.Rows() != 3 || u.Cols() != 1) {
        throw std::invalid_argument(
            "Hat: u 必须是 3x1 列向量 (D19 全库列向量 + 分子布局)");
    }
    const T ux = u(0, 0), uy = u(1, 0), uz = u(2, 0);

    Matrix<T> K(3, 3);
    K(0, 0) = T(0);  K(0, 1) = -uz;    K(0, 2) = uy;
    K(1, 0) = uz;    K(1, 1) = T(0);   K(1, 2) = -ux;
    K(2, 0) = -uy;   K(2, 1) = ux;     K(2, 2) = T(0);
    return K;
}

// ---- Vee: so(3) -> R^3 ----------------------------------------------------
//
//   Vee(K) = (K_21, K_02, K_10)      与 Hat 互逆 (在反对称矩阵上)
//
// 用途: 从旋转矩阵反解轴角时的原料 ——
//   θ = atan2( ‖vee(R - Rᵀ)‖/2 , (tr R - 1)/2 )
// 这个 atan2 口径是良态的; 换成 acos((tr-1)/2) 会在 θ→0 处导数发散
// (report.md §7 第 1 条, Python 侧踩过两次假误差)。见 Quaternion::AngleOf。
template <typename T>
Matrix<T> Vee(const Matrix<T>& K)
{
    if (K.Rows() != 3 || K.Cols() != 3) {
        throw std::invalid_argument("Vee: K 必须是 3x3 矩阵 (so(3) 元素)");
    }
    Matrix<T> u(3, 1);
    u(0, 0) = K(2, 1);
    u(1, 0) = K(0, 2);
    u(2, 0) = K(1, 0);
    return u;
}
