// quaternion.cpp -- 四元数验收 (B1)
//
// ---------------------------------------------------------------------------
// 判据的设计原则 (graphics 3.3): **判据独立于实现路径**
// ---------------------------------------------------------------------------
//   [A] ① 基本代数    构造 / 加减 / 一元负 / 标量乘          双重精确 + 容差
//   [B] ② 共轭与范数  Conjugate / SquaredNorm / Norm
//   [C] ① 的魂        Hamilton 积 —— 乘法表 + 反自同构 + 非交换 + 双覆盖
//   [D] ② 逆与归一化  **左右逆都测** + 非单位四元数 + 与共轭区分
//   [E] ③ 轴角与矩阵  FromAxisAngle / ToMatrix / FromMatrix / Rotate
//   [F] ④ 角度与插值  AngleOf / Slerp / Nlerp
//   [G] 契约          零四元数 / 零轴 / 形状非法 必须 throw
//
// 每一节都尽量用**两条独立代码路径对拍**, 而不是拿手算常数比对:
//   · Rotate(v) vs ToMatrix()*v            —— 四元数乘 vs 矩阵乘
//   · FromMatrix(ToMatrix(q)) vs q         —— 往返
//   · Slerp 的匀速性用 atan2 口径量角       —— 避开 arccos 病态
//
// ---------------------------------------------------------------------------
// 分域纪律 (3.3 第 1 条)
// ---------------------------------------------------------------------------
// T = double 与 T = float **各跑一遍** (模板 bug 只在实例化时暴露;
// 本项目单域探针连续误判过 5 次, 正反两个方向都栽过)。
// float 的阈值放宽到 ~1e-5 (float eps ≈ 1.2e-7, sin/cos 相对误差 ~1e-7)。
//
// 框架: googletest + TYPED_TEST (2026-10-30 迁移) —— 双域由 TYPED_TEST_SUITE 保证,
// 不再靠手写两遍调用; 判据编号 [X#] 通过 SCOPED_TRACE 保留, 与 test-map.md 对得上。
//
// ---------------------------------------------------------------------------
// 已知盲区 (写下来防误判)
// ---------------------------------------------------------------------------
// ① **叉积符号写反抓不到往返测试** (正反抵消) —— 只有 [C1]-[C5] 基元乘法表
//    和 [C11] 反自同构能抓住。
// ② **Hat 式的整体翻号**同理: Quaternion 取负等于同一旋转, 所以
//    "q 与 -q 等价" 这件事让**所有"比 |Dot|"的判据都看不见符号错**。
//    符号必须靠基元乘法表钉。
// ③ FromMatrix 的**四个分支**必须各自走到, 否则分支 2/3/4 的符号错误全漏。
//    只测 R=I 是无效判据 (那时 tr=3, 永远走分支 1)。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                    | 边界 / 病态输入            | 专抓什么错
//  ----|-------------------------|----------------------------|----------------------
//   A  | 1e-6 abs; 单位元用 ==   | 零四元数、极值分量          | 构造/四则的顺序与符号
//   B  | |q|²=30 逐位精确        | 1e-300 极小非零 (IsZero)   | 共轭漏翻号 / 容差误判近零
//   C  | 1e-6, 结合律放宽 ×10    | 基元 i/j/k、q(θ+2π)        | **叉积符号**、乘法顺序
//   D  | 1e-6; q⁻¹≠q* 差 > 1e-3  | **非单位 q,|q|=√30**       | 写死"只支持单位四元数"
//   E  | 1e-5 (float 同阈)       | θ=0 / θ=π / 非单位轴 / 4 分支| Rotate 量纲、Shepperd 分支漏
//   F  | 匀速 ptp ≤ 1e-4         | 0°→350° 短弧、d>0.9995     | 漏短弧翻号、sin 权重写成线性
//   G  | 异常**类型**精确到子类   | 零四元数 / 零轴 / 3x3 vs 3x1| runtime_error 不派生 logic_error
//
// 关键阈值的来历:
//   · **1e-6 而非 1e-12**: 多处经过 sin/cos (float 相对误差 ~1e-7), 卡太紧会假红。
//   · **float 与 double 同用 1e-5 (E 段)**: float eps ≈ 1.2e-7, 再紧就测的是舍入。
//   · **D2 的 1e-3 是"差足够大"**: 非单位 q 下 q⁻¹ 与 q* 差 O(1), 不是比精度。
//   · **F9 的 1e-4 是 ptp (峰峰值)**: 均匀时应为 0, 只留浮点噪声余量。
//   · **F10 期望 5° 而不是 175°**: 两者差 O(1), 一抓一个准; 阈值 1e-4 够宽。
//   · **G 段用 EXPECT_EXACT_THROW**: std::invalid_argument 与 std::out_of_range
//     都派生自 logic_error, 判基类会假绿 (见 demo/common.hpp 的说明)。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh                          # 全量 (注册名: quaternion)
//   ./graphics/build/demo/quat/quaternion            # 单跑
//   ./graphics/build/demo/quat/quaternion --gtest_filter='*Hamilton*'   # 只跑一节
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "quat/quaternion.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// ---- 小工具 ---------------------------------------------------------------

template <typename T>
std::string FmtQ(const Quaternion<T>& q)
{
    char buf[160];
    std::snprintf(buf, sizeof(buf), "(%.9g, %.9g, %.9g, %.9g)",
                  static_cast<double>(q.w), static_cast<double>(q.x),
                  static_cast<double>(q.y), static_cast<double>(q.z));
    return std::string(buf);
}

// 逐分量比较。**故意不写 operator==**: 四元数的"相等"有两个意思 ——
// 逐分量相等 (代数) vs |q0.q1| = 1 (同一旋转, 因为 M(q) = M(-q))。
// 把 == 定义成其中任何一个都会误导另一半。比较留在用例里, 写清是哪一种。
template <typename T>
bool SameQ(const Quaternion<T>& a, const Quaternion<T>& b, T tol)
{
    return std::fabs(a.w - b.w) <= tol && std::fabs(a.x - b.x) <= tol
        && std::fabs(a.y - b.y) <= tol && std::fabs(a.z - b.z) <= tol;
}

// 同一旋转 (不管代表元符号): |q0·q1| = 1
template <typename T>
bool SameRotation(const Quaternion<T>& a, const Quaternion<T>& b, T tol)
{
    return std::fabs(std::fabs(Dot(a, b)) - T(1)) <= tol;
}

template <typename T>
Matrix<T> V3(T a, T b, T c)
{
    return Matrix<T>(3, 1, std::vector<T>{a, b, c});
}

// 矩阵∞范数差 (元素级最大绝对差, 与全库口径一致: D7 只比实数域)
template <typename T>
double MatDiff(const Matrix<T>& A, const Matrix<T>& B)
{
    if (A.Rows() != B.Rows() || A.Cols() != B.Cols()) return 1e9;
    double best = 0.0;
    for (int i = 0; i < A.Rows() * A.Cols(); i++) {
        const double d = std::fabs(static_cast<double>(A.Data()[i]) -
                                   static_cast<double>(B.Data()[i]));
        if (d > best) best = d;
    }
    return best;
}

// 绕 +y 转 theta (半角) —— 多处判据要用
template <typename T>
Quaternion<T> QY(T theta)
{
    const T h = theta / T(2);
    return Quaternion<T>(std::cos(h), T(0), std::sin(h), T(0));
}

// 判据编号保留: [X#] 标签进 SCOPED_TRACE, 失败时与 test-map.md 直接对上
#define STEP(label) SCOPED_TRACE(label)
#define CHECK_TRUE(label, cond)      \
    do {                             \
        STEP(label);                 \
        EXPECT_TRUE(cond);           \
    } while (0)

template <typename T>
class QuaternionTest : public ::testing::Test
{
protected:
    static constexpr T kTol = T(1e-6);      // A~D 段
    static constexpr T kTolE = T(1e-5);     // E~F 段 (经 sin/cos, float 更松)
    static constexpr T kPi = T(3.14159265358979323846);
};

using QuaternionTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(QuaternionTest, QuaternionTypes);

// ---- [A] ① 基本代数 -------------------------------------------------------

TYPED_TEST(QuaternionTest, Algebra)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = TestFixture::kTol;

    const Q q0;
{ STEP("[A1] 默认构造 = 单位元 (1,0,0,0)");
    EXPECT_TRUE(q0.w == T(1) && q0.x == T(0) && q0.y == T(0) && q0.z == T(0))
    << FmtQ(q0); }


    CHECK_TRUE("[A2] Identity() 与默认构造一致", SameQ(q0, Q::Identity(), T(0)));

    const Q a(T(1), T(2), T(3), T(4));
{ STEP("[A3] 四参构造逐位就位");
    EXPECT_TRUE(a.w == T(1) && a.x == T(2) && a.y == T(3) && a.z == T(4)) << FmtQ(a); }


    const Q b(T(0.5), T(-1.5), T(2), T(-0.25));
{ STEP("[A4] 加法逐位");
    EXPECT_TRUE(SameQ(a + b, Q(T(1.5), T(0.5), T(5), T(3.75)), tol)) << FmtQ(a + b); }

{ STEP("[A5] 减法逐位");
    EXPECT_TRUE(SameQ(a - b, Q(T(0.5), T(3.5), T(1), T(4.25)), tol)) << FmtQ(a - b); }


    CHECK_TRUE("[A6] 一元负: a + (-a) = 0 (逐位精确)",
               SameQ(a + (-a), Q(T(0), T(0), T(0), T(0)), T(0)));
{ STEP("[A7] 标量乘 q*2 == 2*q, 且逐位翻倍");
    EXPECT_TRUE(SameQ(a * T(2), T(2) * a, T(0)) && (a * T(2)).x == T(4)) << FmtQ(a * T(2)); }


    CHECK_TRUE("[A8] 标量乘对加法分配",
               SameQ((a + b) * T(2), a * T(2) + b * T(2), tol));
}

// ---- [B] ② 共轭 / 范数 / IsZero -------------------------------------------

TYPED_TEST(QuaternionTest, ConjugateNorm)
{
    using T = TypeParam;
    using Q = Quaternion<T>;

    const Q a(T(1), T(2), T(3), T(4));
    const Q c = a.Conjugate();
{ STEP("[B1] 共轭只翻向量部分");
    EXPECT_TRUE(c.w == T(1) && c.x == T(-2) && c.y == T(-3) && c.z == T(-4)) << FmtQ(c); }


    // 1 + 4 + 9 + 16 = 30, 整数运算在 float/double 下都精确
    CHECK_TRUE("[B2] |q|^2 = 30 (逐位精确)", a.SquaredNorm() == T(30));
{ STEP("[B3] |q| = sqrt(30)");
    EXPECT_NEAR(static_cast<double>(a.Norm()), std::sqrt(30.0), 1e-6); }


    CHECK_TRUE("[B4] 共轭对合: (q*)* = q (逐位精确)",
               SameQ(a.Conjugate().Conjugate(), a, T(0)));

    CHECK_TRUE("[B5] |q*|^2 = |q|^2 (逐位精确)", c.SquaredNorm() == a.SquaredNorm());

    CHECK_TRUE("[B6] IsZero: 全零为 true, 非零为 false",
               Q(T(0), T(0), T(0), T(0)).IsZero() && !a.IsZero());

    // 极小但非零 **不是** zero —— 这是"不带容差"的正面判据
    CHECK_TRUE("[B7] 极小非零四元数 IsZero() == false (不带容差)",
               !Q(T(1e-300), T(0), T(0), T(0)).IsZero()
                   || sizeof(T) == sizeof(float));   // float 下 1e-300 就是 0
}

// ---- [C] ① 的魂: Hamilton 积 ----------------------------------------------
//
// 为什么这一节的判据要这么密: 叉积**符号写反**, 往返测试 (q * q^-1) 和姿态往返
// 都**抓不到** (正反抵消), 只有
//   ① 基元乘法表 (i*j=+k 等 5 条)  —— 钉死 ij=k 还是 ji=k
//   ② 共轭反自同构 (ab)* = b* a*   —— 钉死乘法**顺序**
//   ③ 非交换 [i,k] ≠ 0             —— 钉死"确实不是复数"
// 三者合起来才能唯一确定 Hamilton 约定。

TYPED_TEST(QuaternionTest, HamiltonProduct)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = TestFixture::kTol;

    const Q i(T(0), T(1), T(0), T(0));
    const Q j(T(0), T(0), T(1), T(0));
    const Q k(T(0), T(0), T(0), T(1));
    const Q one(T(1), T(0), T(0), T(0));

    // 乘法表 (Hamilton 约定 ij = k, ji = -k)
{ STEP("[C1] i*j = +k");
    EXPECT_TRUE(SameQ(i * j, k, tol)) << FmtQ(i * j); }

{ STEP("[C2] j*k = +i");
    EXPECT_TRUE(SameQ(j * k, i, tol)) << FmtQ(j * k); }

{ STEP("[C3] k*i = +j");
    EXPECT_TRUE(SameQ(k * i, j, tol)) << FmtQ(k * i); }

{ STEP("[C4] j*i = -k  (反向)");
    EXPECT_TRUE(SameQ(j * i, -k, tol)) << FmtQ(j * i); }

{ STEP("[C5] i*k = -j  (反向)");
    EXPECT_TRUE(SameQ(i * k, -j, tol)) << FmtQ(i * k); }


    // 非交换: ik = -j 而 ki = +j, 差 2j —— 量级 O(1), 不是舍入残渣
    CHECK_TRUE("[C6] 不同基元反交换: ik + ki = 0",
               SameQ(i * k + k * i, Q(T(0), T(0), T(0), T(0)), tol));

    CHECK_TRUE("[C7] 非交换量 |ik - ki|^2 = 4 (O(1), 不是舍入)",
               (i * k - k * i).SquaredNorm() == T(4));

    const Q a(T(1), T(2), T(3), T(4));
    CHECK_TRUE("[C8] 单位元: a*1 = 1*a = a (逐位精确)",
               SameQ(a * one, a, T(0)) && SameQ(one * a, a, T(0)));

    // 结合律 / 分配律
    const Q b(T(0.5), T(-1.5), T(2), T(-0.25));
    const Q c(T(-2), T(0.75), T(1), T(3));
    CHECK_TRUE("[C9] 结合律 (ab)c = a(bc)",
               SameQ((a * b) * c, a * (b * c), tol * T(10)));
    CHECK_TRUE("[C10] 分配律 a(b+c) = ab + ac",
               SameQ(a * (b + c), a * b + a * c, tol * T(10)));

    // 共轭是**反**自同构: (ab)* = b* a* —— 顺序反了就不成立, 专抓乘法写反
    CHECK_TRUE("[C11] 共轭反自同构 (ab)* = b* a*",
               SameQ((a * b).Conjugate(), b.Conjugate() * a.Conjugate(), tol * T(10)));

    // 范数乘性 (Lagrange 恒等式 |ab| = |a||b|): 判据独立于四个分量各自的值
{ STEP("[C12] |ab|^2 = |a|^2 |b|^2");
    EXPECT_LE(static_cast<double>(std::fabs((a * b).SquaredNorm()
    - a.SquaredNorm() * b.SquaredNorm())),
    static_cast<double>(tol * a.SquaredNorm() * b.SquaredNorm() * T(10))); }


    // 同轴旋转相加: q(θ1)·q(θ2) = q(θ1+θ2)  —— 半角逻辑的正确性
    CHECK_TRUE("[C13] q(0.8)*q(0.6) = q(1.4) (同轴相加)",
               SameQ(QY(T(0.8)) * QY(T(0.6)), QY(T(1.4)), tol));

    // 双覆盖: 转 2π 得 -q (转 4π 才回到 q)
    CHECK_TRUE("[C14] q(θ+2π) = -q(θ)  (SU(2) 双覆盖)",
               SameQ(QY(T(0.5) + T(2) * TestFixture::kPi), -QY(T(0.5)), tol));
}

// ---- [D] ② 逆与归一化 -----------------------------------------------------
//
// 这一节的核心是**左右逆都测**: 四元数不交换, 所以 q·q⁻¹ 与 q⁻¹·q 是两次
// 独立的乘法, 只测一边漏掉一半代数。

TYPED_TEST(QuaternionTest, InverseNormalize)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = TestFixture::kTol;

    // 关键: 用**非单位**四元数 (|q| = sqrt(30) ≈ 5.48)。
    // 单位四元数下 q⁻¹ = q*, 看不出来是不是真的除了 |q|²。
    const Q a(T(1), T(2), T(3), T(4));
    const Q inv = a.Inverse();

    CHECK_TRUE("[D1] 左右逆: q·q⁻¹ = q⁻¹·q = 1 (非单位 q)",
               SameQ(a * inv, Q::Identity(), tol) && SameQ(inv * a, Q::Identity(), tol));

    // 非单位四元数下 q⁻¹ != q* —— 差是 O(1), 专抓"写死成只支持单位四元数"
{ STEP("[D2] 非单位 q: q⁻¹ != q*  (差 O(1))");
    EXPECT_FALSE(SameQ(inv, a.Conjugate(), T(1e-3)))
    << "|q⁻¹-w - q*-w| = " << std::fabs(inv.w - a.Conjugate().w); }


    // 精确值: q*/|q|² = (1,-2,-3,-4)/30
{ STEP("[D3] q⁻¹ = q*/|q|² 逐位");
    EXPECT_TRUE(SameQ(inv, Q(T(1) / T(30), T(-2) / T(30), T(-3) / T(30), T(-4) / T(30)), tol))
    << FmtQ(inv); }


    CHECK_TRUE("[D4] (q⁻¹)⁻¹ = q", SameQ(inv.Inverse(), a, tol * T(10)));

    // Normalize: |q/|q|| = 1
    const Q n = a.Normalize();
{ STEP("[D5] |Normalize(q)| = 1");
    EXPECT_NEAR(static_cast<double>(n.Norm()), 1.0, 1e-6); }

    EXPECT_TRUE(SameQ(n, n, T(0))) << FmtQ(n);

    // 归一化不改方向: n 与 a 逐分量成比例
    CHECK_TRUE("[D6] Normalize = q * (1/|q|)",
               SameQ(n, a * (T(1) / a.Norm()), tol));

    // 单位四元数: Normalize 是恒等 (幂等)
    CHECK_TRUE("[D7] Normalize 幂等 (已单位化再归不变)", SameQ(n.Normalize(), n, tol));

    // 单位四元数下 q⁻¹ == q* —— 这条与 [D2] 一起把"两种情形"都钉住
    CHECK_TRUE("[D8] 单位 q: q⁻¹ = q*", SameQ(n.Inverse(), n.Conjugate(), tol));
}

// ---- [E] ③ 轴角与矩阵互转 + Rotate ----------------------------------------
//
// Rotate 最靠得住的判据是**与 ToMatrix()*v 对拍** —— 两条完全独立的路径
// (一条走四元数乘 q v q*, 一条走 3x3 矩阵乘)。

TYPED_TEST(QuaternionTest, AxisAngleAndMatrix)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = TestFixture::kTolE;
    const T pi = TestFixture::kPi;

    // ---- FromAxisAngle: 半角 ----
    const Q r = Q::FromAxisAngle(T(0), T(1), T(0), T(1.2));
    const T h = T(1.2) / T(2);
{ STEP("[E1] FromAxisAngle 用**半角**");
    EXPECT_TRUE(SameQ(r, Q(std::cos(h), T(0), std::sin(h), T(0)), tol)) << FmtQ(r); }


    // 单位轴: |q| = 1
{ STEP("[E2] FromAxisAngle 结果单位长");
    EXPECT_NEAR(static_cast<double>(r.Norm()), 1.0, static_cast<double>(tol)); }


    // 非单位轴输入: 内部归一化, 结果与单位轴版本一致
    const Q r2 = Q::FromAxisAngle(T(0), T(5), T(0), T(1.2));      // 轴长 5
    CHECK_TRUE("[E3] 非单位轴 (0,5,0) 归一化后与 (0,1,0) 一致", SameQ(r2, r, tol));

    // θ=0 => 单位元
    CHECK_TRUE("[E4] θ=0 => 单位元 (任意非零轴)",
               SameQ(Q::FromAxisAngle(T(1), T(2), T(3), T(0)), Q::Identity(), tol));

    // 绕 z 转 π => (0,0,0,1), |w| = 0 (cos(π/2)=0)
    CHECK_TRUE("[E5] 绕 z 转 π => (0,0,0,1)",
               SameQ(Q::FromAxisAngle(T(0), T(0), T(1), pi), Q(T(0), T(0), T(0), T(1)), tol));

    // ---- ToMatrix: 逐格对拍解析值 ----
    // 绕 +y 转 θ 的解析矩阵 (标准右手系):
    //   [[ cosθ, 0, sinθ], [0,1,0], [-sinθ, 0, cosθ]]
    const T th = T(0.7);
    const Matrix<T> R = QY(th).ToMatrix();
    const Matrix<T> want(3, 3, std::vector<T>{
        std::cos(th), T(0),        std::sin(th),
        T(0),        T(1),        T(0),
        -std::sin(th), T(0),       std::cos(th)});
{ STEP("[E6] ToMatrix 绕 y: 逐格对拍解析式");
    EXPECT_LE(MatDiff(R, want), static_cast<double>(tol)); }


    // 是旋转矩阵: R^T R = I, det = 1
    const Matrix<T> RtR = R.Transpose() * R;
{ STEP("[E7] ToMatrix: RᵀR = I");
    EXPECT_LE(MatDiff(RtR, Matrix<T>::Identity(3)), static_cast<double>(tol)); }


    const T det = R(0,0)*(R(1,1)*R(2,2) - R(1,2)*R(2,1))
                - R(0,1)*(R(1,0)*R(2,2) - R(1,2)*R(2,0))
                + R(0,2)*(R(1,0)*R(2,1) - R(1,1)*R(2,0));
{ STEP("[E8] ToMatrix: det R = 1");
    EXPECT_NEAR(static_cast<double>(det), 1.0, static_cast<double>(tol)); }


    // ToMatrix 内部归一化: 非单位 q 与它的归一化给同一个 R
    const Q scaled = QY(th) * T(3.7);
{ STEP("[E9] ToMatrix 内部归一化 (q*3.7 与 q 同 R)");
    EXPECT_LE(MatDiff(scaled.ToMatrix(), R), static_cast<double>(tol)); }


    // ---- FromMatrix: 四分支**各自走到** ----
    // 分支判据: t = tr(R) > 0 走分支 1; 否则看 R00/R11/R22 谁最大。
    // 绕 x/y/z 转 π 时 tr = -1, 且对应对角元为 1 => 分别命中分支 2/3/4。
    // 只测 R=I (tr=3) 会让分支 2/3/4 全漏 —— 那是无效判据。
    auto branchCheck = [&](const Q& q, const char* name) {
        const Q back = Q::FromMatrix(q.ToMatrix());
        STEP(name);
        EXPECT_TRUE(SameRotation(back, q, tol)) << "back = " << FmtQ(back);
    };
    branchCheck(QY(T(0.6)), "[E10] 分支1 (θ=0.6, tr>0)");
    branchCheck(Q::FromAxisAngle(T(1), T(0), T(0), pi), "[E10] 分支2 (绕x转π, R00=1)");
    branchCheck(Q::FromAxisAngle(T(0), T(1), T(0), pi), "[E10] 分支3 (绕y转π, R11=1)");
    branchCheck(Q::FromAxisAngle(T(0), T(0), T(1), pi), "[E10] 分支4 (绕z转π, R22=1)");
    // 再加一条一般轴: 混合轴 120° (三个对角元相等, 考验并列)
    branchCheck(Q::FromAxisAngle(T(1), T(1), T(1), T(2) * pi / T(3)), "[E10] 混合轴 120°");

    // ---- Rotate: 与 ToMatrix()*v 对拍 (两条独立路径) ----
    const Q q = Q::FromAxisAngle(T(0.3), T(-0.7), T(0.5), T(1.1));
    const Matrix<T> v = V3<T>(T(1.2), T(-0.4), T(2.3));
    const Matrix<T> viaQ = q.Rotate(v);
    const Matrix<T> viaM = q.ToMatrix() * v;
{ STEP("[E11] Rotate 返回 3x1 列向量 (不是矩阵!)");
    EXPECT_EQ(viaQ.Rows(), 3); }

    EXPECT_EQ(viaQ.Cols(), 1);
{ STEP("[E12] Rotate(v) == ToMatrix()*v");
    EXPECT_LE(MatDiff(viaQ, viaM), static_cast<double>(tol)); }


    // 保长
    auto dot3 = [](const Matrix<T>& a) {
        return a(0,0)*a(0,0) + a(1,0)*a(1,0) + a(2,0)*a(2,0);
    };
{ STEP("[E13] Rotate 保长: |v'| = |v|");
    EXPECT_LE(static_cast<double>(std::fabs(dot3(viaQ) - dot3(v))),
    static_cast<double>(tol) * static_cast<double>(dot3(v)) * 10.0); }


    // 轴上不动点: 绕 y 转任意角, v = (0,1,0) 不变
    const Matrix<T> axisY = V3<T>(T(0), T(1), T(0));
{ STEP("[E14] Rotate: 轴上的向量不动");
    EXPECT_LE(MatDiff(QY(T(1.1)).Rotate(axisY), axisY), static_cast<double>(tol)); }

{ STEP("[E15] Rotate(单位四元数) = 恒等");
    EXPECT_LE(MatDiff(Q::Identity().Rotate(v), v), static_cast<double>(tol)); }

{ STEP("[E16] Rotate: q 与 -q 结果相同 (双覆盖)");
    EXPECT_LE(MatDiff(q.Rotate(v), (-q).Rotate(v)), static_cast<double>(tol)); }

}

// ---- [F] ④ 角度与插值 -----------------------------------------------------

TYPED_TEST(QuaternionTest, AngleAndInterpolation)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = TestFixture::kTolE;
    const T pi = TestFixture::kPi;

    // ---- AngleOf ----
{ STEP("[F1] AngleOf(单位元) = 0");
    EXPECT_NEAR(static_cast<double>(Q::Identity().AngleOf()), 0.0, static_cast<double>(tol)); }

{ STEP("[F2] AngleOf 还原 θ");
    EXPECT_NEAR(static_cast<double>(QY(T(1.2)).AngleOf()), 1.2, static_cast<double>(tol)); }


    // **最短表示**: 转 1.5π 应返回 0.5π (而不是 1.5π)
    // 这是 w<0 翻号那条分支的唯一判据 —— 不翻号会得到 1.5π
    const T got = QY(T(1.5) * pi).AngleOf();
{ STEP("[F3] AngleOf 取最短: q(1.5π) => 0.5π  (不翻号会得 1.5π=4.712)");
    EXPECT_NEAR(static_cast<double>(got), static_cast<double>(T(0.5) * pi), 1e-4); }


    // 小角极限: θ=1e-6 时仍准确 (arccos 在这里已经掉精度)
{ STEP("[F4] AngleOf 小角 θ=1e-6 仍准 (atan2 优于 arccos)");
    EXPECT_LE(static_cast<double>(std::fabs(QY(T(1e-6)).AngleOf() - T(1e-6))), 1e-8); }


    // 非单位四元数也正确 (atan2 尺度不变)
{ STEP("[F5] AngleOf 对非单位 q 也正确");
    EXPECT_NEAR(static_cast<double>((QY(T(0.8)) * T(4.2)).AngleOf()), 0.8,
    static_cast<double>(tol)); }


    // ---- Slerp 端点 ----
    const Q a = QY(T(0.3));
    const Q b = Q::FromAxisAngle(T(1), T(0), T(0), T(1.4));

    CHECK_TRUE("[F6] Slerp t=0 => q0 (逐位)", SameQ(Slerp(a, b, T(0)), a, tol));
    CHECK_TRUE("[F7] Slerp t=1 => q1 (比 |Dot|=1, 因可能翻号)",
               SameRotation(Slerp(a, b, T(1)), b, tol));
{ STEP("[F8] Slerp 中点 |q|=1 (权重用 sin 的判据)");
    EXPECT_NEAR(static_cast<double>(Slerp(a, b, T(0.5)).Norm()), 1.0,
    static_cast<double>(tol)); }


    // ---- Slerp 匀速性: 相邻增量等长 ----
    {
        const int N = 32;
        T dmin = T(1e9), dmax = T(0);
        Q prev = Slerp(a, b, T(0));
        for (int i = 1; i <= N; i++) {
            const Q cur = Slerp(a, b, T(i) / T(N));
            const T d = Dot(prev, cur);              // 相邻步的四元数内积
            const T ang = T(2) * std::acos(std::fabs(d) > T(1) ? T(1) : std::fabs(d));
            if (ang < dmin) dmin = ang;
            if (ang > dmax) dmax = ang;
            prev = cur;
        }
{ STEP("[F9] Slerp 匀速: 逐帧角增量 ptp (阈值 1e-4)");
    EXPECT_LE(static_cast<double>(dmax - dmin), 1e-4); }

    }

    // ---- Slerp 短弧: 这是**必需件**, 不是优化 ----
    // q0 = 0°, q1 = 绕 z 转 350°。短弧应走 -10° (中点约 -5°)。
    // 不翻号会走 +175° 长弧, 中点差 O(1)。
    {
        const T deg = pi / T(180);
        const Q z0 = Q::FromAxisAngle(T(0), T(0), T(1), T(0));
        const Q z1 = Q::FromAxisAngle(T(0), T(0), T(1), T(350) * deg);
        const T mid = Slerp(z0, z1, T(0.5)).AngleOf();
        const T want = T(5) * deg;                  // 最短表示 => 5°
{ STEP("[F10] Slerp 短弧: 0°→350° 中点 = 5° (不翻号会得 175°)");
    EXPECT_LE(static_cast<double>(std::fabs(mid - want)), 1e-4); }

    }

    // 端点同向 (d>0.9995) 走退化分支也不出 NaN
    {
        const Q q = QY(T(1.0));
        const Q r = Slerp(q, q * T(1.0000001), T(0.5));
{ STEP("[F11] Slerp 近同向 (d>0.9995) 退化 nlerp, 无 NaN");
    EXPECT_TRUE(std::isfinite(static_cast<double>(r.w))); }

        EXPECT_NEAR(static_cast<double>(r.Norm()), 1.0, static_cast<double>(tol));
    }

    // ---- Nlerp ----
    CHECK_TRUE("[F12] Nlerp t=0 => q0", SameQ(Nlerp(a, b, T(0)), a.Normalize(), tol));
    CHECK_TRUE("[F13] Nlerp t=1 => q1", SameRotation(Nlerp(a, b, T(1)), b, tol));
{ STEP("[F14] Nlerp 结果单位长");
    EXPECT_NEAR(static_cast<double>(Nlerp(a, b, T(0.37)).Norm()), 1.0,
    static_cast<double>(tol)); }


    // Nlerp 短弧翻号: q1 = -q0 时若不翻号, (1-t)q0 + t·q1 在 t=0.5 处为零向量
    // => Normalize 抛 domain_error。这条专抓 Nlerp 漏翻号。
    {
        const Q q0 = QY(T(0.9));
        const Q neg = -q0;                          // 同一旋转的另一个代表元
        bool ok = true;
        Q mid;
        try {
            mid = Nlerp(q0, neg, T(0.5));
        } catch (...) {
            ok = false;
        }
{ STEP("[F15] Nlerp(q, -q) 不爆 (短弧翻号必需)");
    EXPECT_TRUE(ok) << "抛异常了: 说明 Nlerp 没做短弧翻号"; }

        if (ok) EXPECT_TRUE(SameRotation(mid, q0, T(1e-4)));
    }

    // Nlerp **不匀速** —— 断言它不匀, 防止以后有人"优化"掉 Slerp。
    // 这是反向判据: 若 Nlerp 变成匀速, 说明它偷偷改成了 Slerp 的实现。
    {
        const int N = 32;
        T dmin = T(1e9), dmax = T(0);
        Q prev = Nlerp(a, b, T(0));
        for (int i = 1; i <= N; i++) {
            const Q cur = Nlerp(a, b, T(i) / T(N));
            const T d = std::fabs(Dot(prev, cur));
            const T ang = T(2) * std::acos(d > T(1) ? T(1) : d);
            if (ang < dmin) dmin = ang;
            if (ang > dmax) dmax = ang;
            prev = cur;
        }
{ STEP("[F16] Nlerp 确实不匀速 (ptp > 1e-6)");
    EXPECT_GT(static_cast<double>(dmax - dmin), 1e-6); }

    }
}

// ---- [G] 契约: 三类错误必须 throw -----------------------------------------
//
// 这里比老写法更严: 老写法 catch (const domain_error&) 会放过派生类,
// 现在判**精确类型** (见 demo/common.hpp 头注释)。

TYPED_TEST(QuaternionTest, Contract)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const Q zero(T(0), T(0), T(0), T(0));

    // 零四元数求逆 / 归一化 -> domain_error (数学上无定义, 不是参数错)
{ STEP("[G1] 零四元数 Inverse => domain_error");
    EXPECT_EXACT_THROW((void)zero.Inverse(), std::domain_error); }

{ STEP("[G2] 零四元数 Normalize => domain_error");
    EXPECT_EXACT_THROW((void)zero.Normalize(), std::domain_error); }


    // 零轴 -> invalid_argument (参数非法)
{ STEP("[G3] FromAxisAngle 零轴 => invalid_argument");
    EXPECT_EXACT_THROW((void)Q::FromAxisAngle(T(0), T(0), T(0), T(1)),
    std::invalid_argument); }


    // 非 3x3 -> invalid_argument
{ STEP("[G4] FromMatrix 非 3x3 => invalid_argument");
    EXPECT_EXACT_THROW((void)Q::FromMatrix(Matrix<T>(2, 2)), std::invalid_argument); }

{ STEP("[G5] FromMatrix 3x1 => invalid_argument");
    EXPECT_EXACT_THROW((void)Q::FromMatrix(Matrix<T>(3, 1)), std::invalid_argument); }


    // Rotate 非 3x1 -> invalid_argument (列向量口径, D19)
{ STEP("[G6] Rotate 传 3x3 => invalid_argument");
    EXPECT_EXACT_THROW((void)Q::Identity().Rotate(Matrix<T>(3, 3)), std::invalid_argument); }

{ STEP("[G7] Rotate 传 1x3 行向量 => invalid_argument");
    EXPECT_EXACT_THROW((void)Q::Identity().Rotate(Matrix<T>(1, 3)), std::invalid_argument); }

}
