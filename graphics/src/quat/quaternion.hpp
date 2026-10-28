#pragma once

// ===========================================================================
// quat/quaternion.hpp —— 四元数与三维旋转 (B1 / V2.5 / V2.6)
// ===========================================================================
//
// 定位 (D16 定案): 四元数**封死在本目录**, 不做矩阵元素 ——
//   不进 Matrix<Quaternion<T>>, 因此**不需要** ScalarTraits<Quaternion<T>>
//   (zero / one / fromReal / conj / abs2 一整套都不欠)。
//
// 与 dual.hpp 的关系: **只共用"类的形态"** —— 模板类 + 公开成员 + 无复合赋值 +
//   无比较运算符。代数是另一套, 别类推:
//     · Dual 交换 (a*b == b*a 是铁律单测)
//     · 四元数**不交换** (多一项 v1 x v2) ⇒ 四元数**不得**进那条单测
//     · Dual 的 sin/cos/exp/ln/pow/sqrt 在这里**一个都没有**: 那些是标量语义;
//       四元数只有 Exp/Log (流形指数映射), 且已被 FromAxisAngle / AngleOf 覆盖。
//
// 默认构造给**单位元**, 不给零四元数 —— 这是与 Dual 唯一"故意不一致"的地方:
//   Dual() = (0,0) 是合法的零对偶数, 但 Quaternion() = (0,0,0,0) 没有逆、不是旋转。
//   默认建造出一个合法旋转, 比默认建造出一个数学废物安全。
//
// ---------------------------------------------------------------------------
// 实现进度: 四步全部落地
// ---------------------------------------------------------------------------
//   ① 基本代数     + / - / 一元 - / 标量乘 / Hamilton 积 ...... 已实现
//   ② 共轭与范数   Conjugate / SquaredNorm / Norm /
//                  IsZero / Inverse / Normalize ............... 已实现
//   ③ 互转 + 旋转  FromAxisAngle / ToMatrix / FromMatrix /
//                  Rotate ................................... 已实现
//   ④ 插值         AngleOf / Dot / Slerp / Nlerp ............. 已实现
//
// T 域约束: static_assert(is_floating_point_v<T>) —— 只有 float / double /
//   long double。所以**不用** ScalarTraits<T>::zero()/one(): 直接 T(0)/T(1)
//   就是精确的, 走 traits 只多一层间接。
//
// ---------------------------------------------------------------------------
// 错误类型三分 (与 numeric-contract.md §5 一致)
// ---------------------------------------------------------------------------
//   · 形状/参数非法         -> std::invalid_argument  (v 不是 3x1、轴是零向量)
//   · 数学上无定义 (定义域)  -> std::domain_error      (零四元数求逆/归一化)
//   · 「还没实现」           -> std::logic_error       (本文件已无)
//
// ⚠ **别用 std::runtime_error**: 它不派生自 std::logic_error, 语义上也错
//   (那是"运行环境出问题", 不是"输入不在定义域")。曾经的坑。
//
// ---------------------------------------------------------------------------
// 用法
// ---------------------------------------------------------------------------
//   using Q = Quaternion<double>;
//   Q a = Q::FromAxisAngle(0.0, 1.0, 0.0, 1.2);      // 绕 +y 转 1.2 rad
//   Matrix<double> R = a.ToMatrix();                  // 3x3 旋转矩阵
//   Matrix<double> v(3, 1, {1.0, 0.0, 0.0});
//   Matrix<double> w = a.Rotate(v);                   // = q v q*, 3x1
//   Q b = Q::FromMatrix(R);                           // Shepperd 反解
//   Q c = Slerp(a, b, 0.5);                           // 球面插值 (匀速)
//   double th = a.AngleOf();                          // 2*atan2(|v|, w), [0,pi]
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

// Nlerp 在 Slerp 里被调用 (退化分支), 所以要先声明。
template <typename T> class Quaternion;
template <typename T> Quaternion<T> Nlerp(const Quaternion<T>&, const Quaternion<T>&, T);

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

    // ---- ① 基本代数 -------------------------------------------------------

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
    //
    // 刻意**不提供** operator/(T): 四元数是可除代数, 但 q1/q2 有左右之分
    // (q1*q2^-1 还是 q2^-1*q1), 一提供 / 就有二义性。要除就显式写 * (1/s)。
    Quaternion operator*(T s) const { return Quaternion(w * s, x * s, y * s, z * s); }
    friend Quaternion operator*(T s, const Quaternion& q) { return q * s; }

    // ---- ① 的魂: Hamilton 积 ----------------------------------------------
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

    // ---- ② 共轭与范数 -----------------------------------------------------

    Quaternion Conjugate() const { return Quaternion(w, -x, -y, -z); }

    // 平方范数单列: Inverse 用 q*/|q|^2 就不必开方; 只有 Normalize 才需要 Norm
    T SquaredNorm() const { return w * w + x * x + y * y + z * z; }
    T Norm() const { return std::sqrt(SquaredNorm()); }

    // 是否精确为零四元数 (四个分量全为 0)。零四元数没有逆、不能归一化。
    //
    // 为什么**不带容差**: 近零四元数仍然可逆 (只是逆很大), 用容差判会把
    // 合法输入拒掉 —— D12 的两层口径里, 这属于"数学上无定义"那一层,
    // 不是"数值上不可靠"那一层。判据只有一条: 精确全零。
    bool IsZero() const
    {
        return w == T(0) && x == T(0) && y == T(0) && z == T(0);
    }

    // q^-1 = q* / |q|^2
    //
    // 单位四元数下有捷径 q^-1 = q* —— 但**别写死成"只支持单位四元数"**:
    // Slerp 的中间量与 FromMatrix 的原值都可能暂时不是单位长度。
    Quaternion Inverse() const
    {
        if (IsZero()) {
            throw std::domain_error(
                "Quaternion::Inverse: 零四元数没有逆 (判据: 四个分量精确全为 0)");
        }
        return Conjugate() * (T(1) / SquaredNorm());
    }

    Quaternion Normalize() const
    {
        if (IsZero()) {
            throw std::domain_error(
                "Quaternion::Normalize: 零四元数无法归一化 (判据: 四个分量精确全为 0)");
        }
        return *this * (T(1) / Norm());
    }

    // ---- ③ 与旋转矩阵互转 + 旋转向量 ---------------------------------------

    // q = cos(θ/2) + sin(θ/2)·û —— **半角**是四元数的全部秘密:
    // 四元数乘法 = 旋转复合, 所以"转 θ"要两次半角。
    // 轴 (ux,uy,uz) 先归一化; 传零轴是非法输入 (θ=0 时轴无定义, 调用方自己判)。
    //
    // 零轴判据同样**不带容差**: 单位化是尺度不变的除法, 极小的非零轴归一化后
    // 方向完全正确, 数值上没问题。只有精确零轴才真的无定义。
    // (曾经写成 norm <= DefaultTol<T>() —— 一是签名要求传 n 编不过, 二是容差
    //  在这里没有正当性: 它会把合法的微小轴误判成非法。)
    static Quaternion FromAxisAngle(T ux, T uy, T uz, T theta)
    {
        const T norm = std::sqrt(ux * ux + uy * uy + uz * uz);
        if (norm == T(0)) {
            throw std::invalid_argument(
                "Quaternion::FromAxisAngle: 旋转轴不能是零向量 (方向无定义)");
        }

        const T h = theta / T(2);
        const T c = std::cos(h);
        const T s = std::sin(h) / norm;      // 一次除法; 别写成 ux/norm*sin(h) 三遍
        return Quaternion(c, ux * s, uy * s, uz * s);
    }

    /*
        v' = cosθ * v + sinθ * u x v + (1 - cosθ) * u x u x v = R * v
        u x v = [u]v = (uu^T - I) * v
        1. w^2 - |v|^2 = cosθ
        2. 2vv^T = (1 - cosθ)uu^T, v = sin(θ/2)*u
        3. 2w[v] = sinθ * [u]
    */

    // 单位四元数 -> 3x3 旋转矩阵。
    // 注意: 内部**自己归一化** (与 Python 版 quat_to_matrix 一致), 别假设调用方已归好。
    //
    // 逐格写, 别用对称性推: R(0,1)=2(xy-wz) 与 R(1,0)=2(xy+wz) **不是转置关系**
    // (差在 w 项的符号)。用"对称"推导必错。
    Matrix<T> ToMatrix() const
    {
        const Quaternion q = Normalize();
        const T W = q.w, X = q.x, Y = q.y, Z = q.z;

        Matrix<T> R(3, 3);
        R(0, 0) = T(1) - T(2) * (Y * Y + Z * Z);
        R(0, 1) = T(2) * (X * Y - W * Z);
        R(0, 2) = T(2) * (X * Z + W * Y);
        R(1, 0) = T(2) * (X * Y + W * Z);
        R(1, 1) = T(1) - T(2) * (X * X + Z * Z);
        R(1, 2) = T(2) * (Y * Z - W * X);
        R(2, 0) = T(2) * (X * Z - W * Y);
        R(2, 1) = T(2) * (Y * Z + W * X);
        R(2, 2) = T(1) - T(2) * (X * X + Y * Y);
        return R;
    }

    /*
        1. tr(R) = 1 + 2cosθ     ← 注意符号! 本文件 R 的写法下 (1-2(y²+z²) 构造),
                                    对角和 = 4w²-1 = 1+2cosθ
        2.  R(2,1) - R(1,2) = 4w·x
            R(0,2) - R(2,0) = 4w·y
            R(1,0) - R(0,1) = 4w·z

        R = [[1-2(y²+z²), 2(xy−wz), 2(xz+wy)],
             [2(xy+wz), 1-2(x²+z²), 2(yz−wx)],
             [2(xz−wy), 2(yz+wx), 1-2(x²+y²)]]

        分支判据只看 t > 0 (即 cosθ > -1/2, θ < 120°), 与 t 的具体形式无关,
        所以上面那条符号注释不影响代码 —— 但写错了会误导后来的人。
    */

    // 3x3 旋转矩阵 -> 单位四元数。**必须四分支 Shepperd**:
    //   按 tr(R) 与三个对角元里最大的那个选分支, 否则在 θ→π (tr→-1) 处
    //   除以小量、精度崩。单分支 w = sqrt(1+tr)/2 是教科书上最常见的错写法。
    // 返回的代表元符号是任意的 (M(q) = M(-q)) —— 想和别的 q 比, 请比 |q0·q1|。
    //
    // 契约: R 必须是 SO(3) 元素。**不做正交性检查** (要算 RᵀR−I, O(n³) 还带
    //   阈值选取问题)。传非旋转矩阵进来, 结果是 NaN 或不代表任何旋转。
    static Quaternion FromMatrix(const Matrix<T>& R)
    {
        if (R.Rows() != 3 || R.Cols() != 3)
            throw std::invalid_argument("Quaternion::FromMatrix: R 必须是 3x3");

        const T t = R(0, 0) + R(1, 1) + R(2, 2);   // Matrix 无 Trace(), 直接加
        T w, x, y, z;

        if (t > T(0)) {                                        // 分支 1: w 当主角
            const T s = std::sqrt(t + T(1)) * T(2);
            w = T(0.25) * s;
            x = (R(2, 1) - R(1, 2)) / s;
            y = (R(0, 2) - R(2, 0)) / s;
            z = (R(1, 0) - R(0, 1)) / s;
        } else if (R(0, 0) > R(1, 1) && R(0, 0) > R(2, 2)) {    // 分支 2: x 当主角
            const T s = std::sqrt(T(1) + R(0, 0) - R(1, 1) - R(2, 2)) * T(2);
            w = (R(2, 1) - R(1, 2)) / s;
            x = T(0.25) * s;
            y = (R(0, 1) + R(1, 0)) / s;
            z = (R(0, 2) + R(2, 0)) / s;
        } else if (R(1, 1) > R(2, 2)) {                         // 分支 3: y 当主角
            const T s = std::sqrt(T(1) + R(1, 1) - R(0, 0) - R(2, 2)) * T(2);
            w = (R(0, 2) - R(2, 0)) / s;
            x = (R(0, 1) + R(1, 0)) / s;
            y = T(0.25) * s;
            z = (R(1, 2) + R(2, 1)) / s;
        } else {                                                // 分支 4: z 当主角
            const T s = std::sqrt(T(1) + R(2, 2) - R(0, 0) - R(1, 1)) * T(2);
            w = (R(1, 0) - R(0, 1)) / s;
            x = (R(0, 2) + R(2, 0)) / s;
            y = (R(1, 2) + R(2, 1)) / s;
            z = T(0.25) * s;
        }
        return Quaternion(w, x, y, z).Normalize();
    }

    // 直接算 v' = q v q*, 不走矩阵:
    //   v' = v + 2w(û x v) + 2û x (û x v)        (与 q v q* 等价)
    //
    // 这里走"纯四元数积"这条更短的等价路径: 把 v 当成实部为 0 的四元数
    // p = (0, v), 则 q p q* 仍是纯四元数, 其向量部分就是旋转后的 v。
    //
    // 这条路径与 ToMatrix()*v **完全独立** (一条四元数乘, 一条矩阵乘) ——
    // 所以 demo 用它们互相对拍, 这是 Rotate 最靠得住的判据。
    //
    // 入参约定: 3x1 列向量 (D19); 其他形状 throw (契约写清楚, 别静默接受)。
    // 返回 3x1 列向量。
    //
    // ⚠ 曾经错写成 `return v_prime.ToMatrix();` —— 返回 3x3 是**量纲错**
    //   (把"旋转后的向量"当成了"旋转矩阵"), 而且 v_prime 是纯四元数,
    //   归一化后 w=0 那一维直接丢了。这类错编译器不报, 只有形状判据能抓。
    Matrix<T> Rotate(const Matrix<T>& v) const
    {
        if (v.Rows() != 3 || v.Cols() != 1)
            throw std::invalid_argument(
                "Quaternion::Rotate: v 必须是 3x1 列向量 (D19 全库列向量)");

        const Quaternion n = Normalize();
        const Quaternion p(T(0), v(0, 0), v(1, 0), v(2, 0));
        const Quaternion r = n * p * n.Conjugate();    // q v q* (逆=共轭, n 已单位化)

        Matrix<T> out(3, 1);
        out(0, 0) = r.x;
        out(1, 0) = r.y;
        out(2, 0) = r.z;
        return out;
    }

    // ---- ④ 角度 -----------------------------------------------------------

    // 旋转角, 取 [0, π] 的最短表示。
    //
    // **必须用 2*atan2(|v|, w)**, 不要 2*arccos(w) —— arccos 在 θ→0 处导数发散
    // (相对误差 ~ sqrt(eps)/θ), Python 侧已经因此踩过两次假误差 (report.md §7 第 1 条)。
    // 例: θ=1e-8 时 arccos 直接返回 0, atan2 仍精确。
    //
    // w<0 时先取 |w|: q 与 -q 是同一个旋转, 但 w<0 时 atan2(|v|, w) 会落到
    // (π/2, π] 再乘 2 得 (π, 2π], 不是"最短表示"。取 |w| 等价于整体翻号
    // (|v| 不变, 因为每个分量都翻、模长不变), 只是省掉四次乘法。
    //
    // 注: 零四元数返回 0 (它不是旋转, 但比抛异常更少扰乱调用方);
    //     非单位四元数也正确 (atan2(|v|,w) 是尺度不变的, 不需要先归一化)。
    T AngleOf() const
    {
        const T ww = (w < T(0)) ? -w : w;
        return T(2) * std::atan2(std::sqrt(x * x + y * y + z * z), ww);
    }
};

// ---- ④ 点积 / Slerp / Nlerp (自由函数) ------------------------------------
//
// Dot 用自由函数而不是成员: 插值里要写成对称的 Dot(q0,q1), 写成 q0.Dot(q1) 会
// 让人误以为有方向性 (其实没有)。与 linalg::operator.hpp 的 Dot(Matrix) 是
// 同名不同参的重载, 不冲突 —— 本文件不 include linalg/, 更不可能误调 (D16)。

/*
    cos Ω = q0 · q1 = w0w1 + x0x1 + y0y1 + z0z1
    q(t) = [sin((1-t)Ω)/sinΩ]·q0 + [sin(tΩ)/sinΩ]·q1        t ∈ [0,1]
    r = (q0^-1 · q1), Re(r) = q0·q1
    q0r^t = q(t), q(0) = q0, q(1) = q1
    q(t) = q0·exp(log(q0r^t))
*/
template <typename T>
T Dot(const Quaternion<T>& a, const Quaternion<T>& b)
{
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

// 归一化线性插值: 线性插值后归一化。比 Slerp 快, 但**不匀速** (两端慢中间快)。
//
// 短弧翻号**同样必需**: 不翻时 q1 ≈ -q0 会让 (1-t)q0 + t·q1 穿过零向量附近,
// 归一化把噪声放大; q1 = -q0 时更直接得到零向量 -> Normalize 抛异常。
// (这是 Slerp 必需件 1 在 Nlerp 上的同一件事, 别以为是 Slerp 专属。)
//
// 定义放在 Slerp **之前**: Slerp 的退化分支要调它, 而模板的隐式实例化
// 要求声明可见 (只靠文件开头那行前向声明也能编过, 但放前面更直白)。
template <typename T>
Quaternion<T> Nlerp(const Quaternion<T>& q0, const Quaternion<T>& q1, T t)
{
    Quaternion<T> b = q1;
    if (Dot(q0, q1) < T(0)) {
        b = -q1;
    }
    return (q0 * (T(1) - t) + b * t).Normalize();
}

// Slerp 的三个**必需**件 (缺一条就是 bug, 不是优化):
//   1) d = q0·q1 < 0 时翻 q1 的符号取短弧 —— 不翻会走 350° 长弧
//      (姿态上只差 10°)。这是**必需步骤**, 不是优化。
//   2) d > 0.9995 退化成 nlerp —— 否则 sin(Ω)→0, 权重除零。
//   3) 权重用 sin 不用线性 —— 线性权重会让 |q| 掉到 1 以下 (实测最低 0.857)。
//
// 一个反直觉的点: 这里的 acos(d) **是安全的**, 因为 d>0.9995 已经先走 nlerp 了,
// 剩下的 Ω >= 0.0316 rad, 1/sin(Ω) <= 31.6, 不发散。这跟 AngleOf 那里
// θ 可以任意小是两回事 —— 别看到 acos 就以为也要换成 atan2。
//
// 端点: t=0 精确返回 q0; t=1 返回 ±q1 (翻过号就是 -q1, 代表同一旋转)。
//   所以测端点**不能逐分量比**, 要比 |Dot(result, q1)| = 1。
//
// ⚠ 权重必须乘倒数, 不能写 `expr / std::sin(theta)` —— 四元数没有
//   operator/(T) (见 ① 里"刻意不提供"那段), 那样写编不过。
//   也不该为这一处去加 operator/, 否则 q1/q2 的二义性就漏进来了。
template <typename T>
Quaternion<T> Slerp(const Quaternion<T>& q0, const Quaternion<T>& q1, T t)
{
    T d = Dot(q0, q1);
    Quaternion<T> b = q1;
    if (d < T(0)) {          // 必需件 1: 取短弧
        b = -q1;
        d = -d;
    }

    if (d > T(0.9995)) {     // 必需件 2: 几乎同向, 退化避免除零
        return Nlerp(q0, b, t);
    }

    const T omega = std::acos(d);
    const T inv = T(1) / std::sin(omega);
    return q0 * (std::sin((T(1) - t) * omega) * inv)
         + b  * (std::sin(t * omega) * inv);
}