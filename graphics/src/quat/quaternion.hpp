#pragma once

// ===========================================================================
// quat/quaternion.hpp —— 四元数与三维旋转 (plan.md 阶段 2 第 5 条, B1 / V2.5 / V2.6)
// ===========================================================================
//
// 定位 (D16 定案): 四元数**封死在本目录**, 不做矩阵元素 ——
//   不进 Matrix<Quaternion<T>>, 因此**不需要** ScalarTraits<Quaternion<T>>
//   (zero / one / fromReal / conj / abs2 一整套都不欠)。
//
// 与 dual.hpp 的关系: **只共用"类的形态"** —— 模板类 + 公开成员 + 无复合赋值 +
//   无比较运算符。代数是另一套, 别类推:
//     · Dual 交换 (a*b == b*a 是 matrix.md §1.3 的铁律单测)
//     · 四元数**不交换** (多一项 v1 x v2) ⇒ 四元数**不得**进那条单测
//     · Dual 的 sin/cos/exp/ln/pow/sqrt 在这里**一个都没有**: 那些是标量语义;
//       四元数只有 Exp/Log (流形指数映射), 且已被 FromAxisAngle / AngleOf 覆盖。
//
// 默认构造给**单位元**, 不给零四元数 —— 这是与 Dual 唯一"故意不一致"的地方:
//   Dual() = (0,0) 是合法的零对偶数, 但 Quaternion() = (0,0,0,0) 没有逆、不是旋转。
//   默认建造出一个合法旋转, 比默认建造出一个数学废物安全。
//
// ---------------------------------------------------------------------------
// 实现进度 (四步, 依赖序; 与 dev-progress.md「下一版基线」B1 对齐)
// ---------------------------------------------------------------------------
//   ① 基本代数     + / - / 一元 - / 标量乘 ............ 已实现
//                  Hamilton 积 operator* ............. 待实现 (本文件的魂)
//   ② 共轭与范数   Conjugate / SquaredNorm / Norm ..... 已实现
//                  Inverse / Normalize ................ 待实现
//   ③ 互转 + 旋转  FromAxisAngle / ToMatrix /
//                  FromMatrix / Rotate ................ 待实现
//   ④ 插值         Slerp / Nlerp / AngleOf ............ 待实现
//
// 「待实现」的体一律 throw std::logic_error, **不留静默 no-op**
// (本项目纪律: 宁可明确报错, 也不要静默算错 —— grad_graph 那 4 类 bug 的教训)。
// demo/quat/quaternion.cpp 里对每条占位都有 `断言它 throw` 的用例;
// 你实现掉一条, 对应用例就会变红 —— 那就是提醒你把占位判据换成真判据。
//
// ---------------------------------------------------------------------------
// 预期用法 (写完后的样子)
// ---------------------------------------------------------------------------
//   using Q = Quaternion<double>;
//   Q a = Q::FromAxisAngle(0.0, 1.0, 0.0, 1.2);      // 绕 +y 转 1.2 rad
//   Matrix<double> R = a.ToMatrix();                  // 3x3 旋转矩阵
//   Matrix<double> v(3, 1, {1.0, 0.0, 0.0});
//   Matrix<double> w = a.Rotate(v);                   // = q v q*
//   Q b = Q::FromMatrix(R);                           // Shepperd 反解
//   Q c = Slerp(a, b, 0.5);                           // 球面插值 (匀速)
//   double th = a.AngleOf();                          // 2*atan2(|v|, w)
//
// ---------------------------------------------------------------------------
// 判据来源
// ---------------------------------------------------------------------------
//   graphics/demo/pyapps/rotation/  (A5.1–A5.6, report.md 223 行)
//   先在 Python 侧拿数, 再在 C++ 侧对拍 —— 不要凭公式推数值。
// ===========================================================================

#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "core/matrix.hpp"

template <typename T>
class Quaternion
{
    // 把域写进签名 (2026-10-16 那条教训): Apply 的签名贪婪写成 Matrix<T>,
    // 结果标量域报的是 "undeclared identifier 'dualFunc'", 把人引向"是不是漏了 include"。
    // 这里一行 static_assert 把 T 钉在实数浮点上 —— 越界调用直接给出准确原因。
    static_assert(std::is_floating_point_v<T>,
                  "Quaternion<T>: T 必须是实数浮点 (float / double / long double)。"
                  "四元数的共轭与范数只在实数域有意义。");

public:
    T w;   // 实部
    T x;   // i
    T y;   // j
    T z;   // k

    Quaternion(T w_ = T(1), T x_ = T(0), T y_ = T(0), T z_ = T(0))
        : w(w_), x(x_), y(y_), z(z_) {}

    static Quaternion Identity() { return Quaternion(T(1), T(0), T(0), T(0)); }

    // ---- ① 基本代数 (已实现的部分) ----------------------------------------

    Quaternion operator+(const Quaternion& o) const
    {
        return Quaternion(w + o.w, x + o.x, y + o.y, z + o.z);
    }

    Quaternion operator-(const Quaternion& o) const
    {
        return Quaternion(w - o.w, x - o.x, y - o.y, z - o.z);
    }

    Quaternion operator-() const { return Quaternion(-w, -x, -y, -z); }

    // 标量乘: 两侧都要 —— 插值的加权和靠它, 不靠"两个四元数相加"
    Quaternion operator*(T s) const { return Quaternion(w * s, x * s, y * s, z * s); }
    friend Quaternion operator*(T s, const Quaternion& q) { return q * s; }

    // ---- ① 的魂: Hamilton 积 (待实现) ------------------------------------
    //
    //   (w1 + v1)(w2 + v2) = (w1 w2 - v1.v2) + (w1 v2 + w2 v1 + v1 x v2)
    //
    // 最后一项 v1 x v2 是**唯一的非交换来源**: 漏掉它, 四元数就退化成
    // "一个复数 + 两个死分量", 复合旋转全错。
    //
    // 三个容易错的点:
    //   1) 叉积符号: Hamilton 约定 ij = k, ji = -k (即 v1 x v2 用右手系叉积)。
    //      写反 => 复合旋转方向反。**注意这个错误往返测试抓不到** (正反抵消),
    //      只有「qx*qz != qz*qx」和「与矩阵乘同态 M(a*b) = M(a)M(b)」能抓住。
    //   2) 顺序: 左乘 = 绕世界轴, 右乘 = 绕物体轴 (A5.6); 别把 a/b 对调。
    //   3) 叉积三个分量手写, 别指望 Matrix 帮忙 (quat/ 不依赖 linalg, 见 D16)。
    Quaternion operator*(const Quaternion& other) const
    {
        return Quaternion(w * other.w - x * other.x - y * other.y - z * other.z,
                           w * other.x + x * other.w + y * other.z - z * other.y,
                           w * other.y + y * other.w + z * other.x - x * other.z,
                           w * other.z + z * other.w + x * other.y - y * other.x);
    }

    // ---- ② 共轭与范数 ----------------------------------------------------

    Quaternion Conjugate() const { return Quaternion(w, -x, -y, -z); }

    // 平方范数单列: Inverse 用 q*/|q|^2 就不必开方; 只有 Normalize 才需要 Norm
    T SquaredNorm() const { return w * w + x * x + y * y + z * z; }
    T Norm() const { return std::sqrt(SquaredNorm()); }

    // 待实现 —— 决策点 (D12「除零两层」): 零四元数没有逆。
    //   建议**只对精确 0 抛**, 不要带容差: 近零四元数仍然可逆 (只是逆很大),
    //   用 DefaultTol 判会把合法输入拒掉。
    //   单位四元数下有捷径 q^{-1} = q* —— 但别写死成"只支持单位四元数",
    //   Slerp 的中间量与 FromMatrix 的原值都可能暂时不是单位长度。
    Quaternion Inverse() const
    {
        throw std::logic_error("Quaternion::Inverse: 待实现 (B1 第 ② 步)");
    }

    Quaternion Normalize() const
    {
        throw std::logic_error("Quaternion::Normalize: 待实现 (B1 第 ② 步)");
    }

    // ---- ③ 与旋转矩阵互转 + 旋转向量 (全部待实现) --------------------------

    // q = cos(θ/2) + sin(θ/2)·û —— **半角**是四元数的全部秘密:
    // 四元数乘法 = 旋转复合, 所以"转 θ"要两次半角。
    // 轴 (ux,uy,uz) 先归一化; 传零轴是非法输入 (θ=0 时轴无定义, 调用方自己判)。
    static Quaternion FromAxisAngle(T ux, T uy, T uz, T theta)
    {
        (void)ux; (void)uy; (void)uz; (void)theta;
        throw std::logic_error("Quaternion::FromAxisAngle: 待实现 (B1 第 ③ 步)");
    }

    // 单位四元数 -> 3x3 旋转矩阵。
    // 注意: 内部**自己归一化** (与 Python 版 quat_to_matrix 一致), 别假设调用方已归好。
    Matrix<T> ToMatrix() const
    {
        throw std::logic_error("Quaternion::ToMatrix: 待实现 (B1 第 ③ 步)");
    }

    // 3x3 旋转矩阵 -> 单位四元数。**必须四分支 Shepperd**:
    //   按 tr(R) 与三个对角元里最大的那个选分支, 否则在 θ→π (tr→-1) 处
    //   除以小量、精度崩。单分支 w = sqrt(1+tr)/2 是教科书上最常见的错写法。
    // 返回的代表元符号是任意的 (M(q) = M(-q)) —— 想和别的 q 比, 请比 |q0·q1|。
    static Quaternion FromMatrix(const Matrix<T>& R)
    {
        (void)R;
        throw std::logic_error("Quaternion::FromMatrix: 待实现 (B1 第 ③ 步)");
    }

    // 直接算 v' = q v q* (标量+向量形式, 不构造 4x4):
    //   v' = v + 2w(û x v) + 2û x (û x v)
    // 入参约定: 3x1 列向量; 其他形状 throw (契约写清楚, 别静默接受)。
    Matrix<T> Rotate(const Matrix<T>& v) const
    {
        (void)v;
        throw std::logic_error("Quaternion::Rotate: 待实现 (B1 第 ③ 步)");
    }

    // ---- ④ 插值 (待实现) -------------------------------------------------

    // 旋转角, 取 [0, π] 的最短表示 (w < 0 时先翻号)。
    // **必须用 2*atan2(|v|, w)**, 不要 2*arccos(w) —— arccos 在 θ→0 处导数发散,
    // Python 侧已经因此踩过两次假误差 (report.md §7 第 1 条)。
    T AngleOf() const
    {
        throw std::logic_error("Quaternion::AngleOf: 待实现 (B1 第 ④ 步)");
    }
};

// ---- ④ Slerp / Nlerp (自由函数, 与 Slerp(a,b,t) 的读法一致) ----------------
//
// Slerp 的三个**必需**件 (缺一条就是 bug, 不是优化):
//   1) d = q0·q1 < 0 时翻 q1 的符号取短弧 —— 不翻会走 350° 长弧 (姿态上只差 10°)。
//   2) d > 0.9995 退化成 nlerp —— 否则 sin(Ω)→0, 权重除零。
//   3) 权重用 sin 不用线性 —— 线性权重会让 |q| 掉到 1 以下 (实测最低 0.857)。
template <typename T>
Quaternion<T> Slerp(const Quaternion<T>& q0, const Quaternion<T>& q1, T t)
{
    (void)q0; (void)q1; (void)t;
    throw std::logic_error("Slerp: 待实现 (B1 第 ④ 步)");
}

template <typename T>
Quaternion<T> Nlerp(const Quaternion<T>& q0, const Quaternion<T>& q1, T t)
{
    (void)q0; (void)q1; (void)t;
    throw std::logic_error("Nlerp: 待实现 (B1 第 ④ 步)");
}
