// so3.cpp -- so(3) 三维向量算子验收 (fan-plan A4 李群线的载体)
//
// ---------------------------------------------------------------------------
// 判据的设计原则 (graphics 3.3): **判据独立于实现路径**
// ---------------------------------------------------------------------------
// 这里的每一条都用**两条不同的代码路径**对拍, 而不是拿手算常数比对:
//   [A] Hat 的结构      —— 反对称 K + Kᵀ = 0, 对角精确 0 (结构性, 不受舍入影响)
//   [B] Hat/Vee 互逆    —— Vee∘Hat = id 且 Hat∘Vee = id (两个方向都测)
//   [C] Cross 的代数    —— 基元乘法表 + 反交换 + 正交 (良态口径)
//   [D] Cross vs Hat·b  —— 手写三分量 vs 走 Matrix 乘法, 两条独立路径
//   [E] Hat³ = -|u|²Hat —— Rodrigues 闭式成立的**唯一关键** (report.md §1)
//   [F] Cross vs 四元数 —— 纯四元数积的向量部分, 把 R³ / so(3) / 四元数钉在一起
//   [G] 形状契约        —— 传错形状必须 throw, 不静默越界读 (Apply 的教训)
//
// [F] 放在这里而不是 demo/quat/ 的理由: 它是唯一一条**同时**穿过三条线的恒等式,
//   放在能同时看见三者的地方才有意义。(demo 层不受 D16「quat/ 不依赖 linalg」
//   的约束 —— 那条管的是 src 内部依赖, 不是 demo。)
//
// ---------------------------------------------------------------------------
// 分域纪律 (3.3 第 1 条) + 框架 (googletest, 2026-10-30 迁移)
// ---------------------------------------------------------------------------
// 模板 bug 只在实例化时暴露: double 过不代表 float 过 (本项目连续误判过 5 次,
// 正反两个方向都栽过)。所以每条数值判据都跑 double + float, float 阈值放宽到
// 1e-5 (float eps ≈ 1.2e-7, sin/cos 相对误差 ~1e-7)。
//
// 迁移说明: 老写法是 `TestXxx<double>(ck, "double")` + `TestXxx<float>(ck, "float ")`
//   各手写一遍 —— 现在由 **TYPED_TEST_SUITE(double, float)** 保证两个域都跑。
//   判据编号 [A1]–[G5] 全部保留 (SCOPED_TRACE 携带), 迁移前 51 项 -> 迁移后 51 条。
//   [G] 段与类型无关, 仍是普通 TEST。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                   | 边界 / 病态输入          | 专抓什么错
//  ----|------------------------|--------------------------|------------------------
//   A  | 反对称用 == (逐位精确)  | 一般轴 (1,2,3)           | Hat 六位置符号、归一化
//   B  | 1e-6                   | 零向量、**任意**反对称 K  | Vee 只写对一半 (如 u(1)两次)
//   C  | 1e-6                   | 基元 e_x/e_y/e_z、平行向量| 叉积符号写反、右手系
//   D  | 1e-5 (李括号 ×10)       | 一般向量 a,b             | 李括号因子写成 2 (应为 1)
//   E  | 1e-5 (一般轴 ×10)       | **单位轴** vs 一般轴      | **抓不到整体翻号!** 见盲区
//   F  | 1e-6                   | 平行纯四元数             | 四元数向量部分 ≠ 叉积
//   G  | 异常类型 invalid_argument| 2x1 / 1x3 / 3x1          | 静默越界读 (Apply 的教训)
//
// **E 段的边界必须写清楚: 它是结构判据, 不是符号判据。**
//   注入实测 (2026-10-20): 把 Hat 整体翻号, [E1][E2][E3] 仍全绿 ——
//   Hat³ = -|u|²Hat 两边都是奇次, 整体翻号互相抵消。
//   ⇒ 符号由 [A3] (六位置逐位) 与 [D1] (Cross 对拍) 钉, **不由 [E] 钉**。
//   (2026-10-30 迁移后重跑了这个注入, 结论不变 —— 见本次日志。)
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh                 # 全量 (本用例注册名: so3)
//   ./graphics/build/demo/linalg/so3        # 单跑
//   ./graphics/build/demo/linalg/so3 --gtest_filter='*double*'   # 只跑一个域
//
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/so3.hpp"
#include "quat/quaternion.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#define STEP(label) SCOPED_TRACE(label)

// 判据编号保留: [X#] 标签进 SCOPED_TRACE, 失败时与 test-map.md 直接对上
#define CHECK_TRUE(label, cond)   \
    do {                          \
        STEP(label);              \
        EXPECT_TRUE(cond);        \
    } while (0)

namespace
{

// 3x1 **列**向量 (D19 第二次修订, 2026-10-19: 全库列向量 + 分子布局)
template <typename T>
Matrix<T> V3(T a, T b, T c)
{
    return Matrix<T>(3, 1, std::vector<T>{a, b, c});
}

template <typename T>
std::string FmtV(const Matrix<T>& v)
{
    char buf[160];
    std::snprintf(buf, sizeof(buf), "(%.9g, %.9g, %.9g)",
                  static_cast<double>(v(0, 0)), static_cast<double>(v(1, 0)),
                  static_cast<double>(v(2, 0)));
    return std::string(buf);
}

template <typename T>
bool SameV(const Matrix<T>& a, const Matrix<T>& b, T tol)
{
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) return false;
    for (int i = 0; i < a.Rows(); i++) {
        if (std::fabs(a(i, 0) - b(i, 0)) > tol) return false;
    }
    return true;
}

// 逐分量比较。故意不给 Quaternion 写 operator== —— 四元数的"相等"有两个意思
// (逐分量相等 vs |q0·q1| = 1), 定义成任一个都会误导另一半。
template <typename T>
bool SameQ(const Quaternion<T>& a, const Quaternion<T>& b, T tol)
{
    return std::fabs(a.w - b.w) <= tol && std::fabs(a.x - b.x) <= tol
        && std::fabs(a.y - b.y) <= tol && std::fabs(a.z - b.z) <= tol;
}

// 固定取值, 不用随机源 —— 判据要可复现
template <typename T> Matrix<T> SampleA() { return V3<T>(T(0.3), T(-1.7), T(2.4)); }
template <typename T> Matrix<T> SampleB() { return V3<T>(T(1.9), T(0.6), T(-1.2)); }

}  // namespace

template <typename T>
class So3Test : public ::testing::Test
{
protected:
    static constexpr T kTol = T(1e-6);
};

using So3Types = ::testing::Types<double, float>;
TYPED_TEST_SUITE(So3Test, So3Types);

// ---- [A] Hat 的结构: 反对称 -----------------------------------------------

TYPED_TEST(So3Test, HatStructure)
{
    using T = TypeParam;
    const Matrix<T> u = V3<T>(T(1), T(2), T(3));
    const Matrix<T> K = Hat(u);

    // K + Kᵀ = 0 —— 结构判据, 逐位精确
    bool antisym = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (K(i, j) + K(j, i) != T(0)) antisym = false;
        }
    }
    CHECK_TRUE("[A1] Hat(u) 反对称: K + Kᵀ = 0 (逐位精确)", antisym);

    CHECK_TRUE("[A2] Hat(u) 对角元精确 0",
               K(0, 0) == T(0) && K(1, 1) == T(0) && K(2, 2) == T(0));

    // 六个位置的**符号**: Hamilton 约定 ij = k 的落点
    CHECK_TRUE("[A3] Hat 六位置逐位就位 (右手系)",
               K(2, 1) == T(1) && K(0, 2) == T(2) && K(1, 0) == T(3)
                   && K(1, 2) == T(-1) && K(2, 0) == T(-2) && K(0, 1) == T(-3));

    // 线性: hat(λu) = λ·hat(u) —— 归一化会破坏这条, 所以 Hat **绝不**归一化
    const T lam = T(2.5);
    Matrix<T> su(3, 1);
    for (int i = 0; i < 3; i++) su(i, 0) = u(i, 0) * lam;
    { STEP("[A4] Hat 线性: hat(λu) = λ·hat(u)");
      EXPECT_TRUE(IsApprox(Hat(su), Hat(u) * lam, T(1e-6)))
          << "‖hat(λu) − λ·hat(u)‖∞ = " << NormInf(Hat(su) - Hat(u) * lam); }
}

// ---- [B] Hat / Vee 互逆 ---------------------------------------------------

TYPED_TEST(So3Test, HatVeeRoundTrip)
{
    using T = TypeParam;
    const T tol = T(1e-6);
    const Matrix<T> u = SampleA<T>();

    // 方向一: Vee(Hat(u)) = u
    { STEP("[B1] Vee(Hat(u)) = u");
      EXPECT_TRUE(SameV(Vee(Hat(u)), u, tol)) << "实测 " << FmtV(Vee(Hat(u))); }

    // 方向二: Hat(Vee(K)) = K —— K 取**任意**反对称阵, 不取 Hat(u) 以免同义反复
    const Matrix<T> K(3, 3, std::vector<T>{T(0), T(-3.1), T(2.2),
                                           T(3.1), T(0), T(-0.7),
                                           T(-2.2), T(0.7), T(0)});
    { STEP("[B2] Hat(Vee(K)) = K (任意反对称 K)");
      EXPECT_TRUE(IsApprox(Hat(Vee(K)), K, tol))
          << "‖Hat(Vee(K)) − K‖∞ = " << NormInf(Hat(Vee(K)) - K); }

    // 零向量
    { STEP("[B3] Cross(0, u) = 0 (逐位精确)");
      EXPECT_TRUE(SameV(Cross(V3<T>(T(0), T(0), T(0)), u), V3<T>(T(0), T(0), T(0)), T(0))); }
}

// ---- [C] Cross 的代数性质 -------------------------------------------------

TYPED_TEST(So3Test, CrossAlgebra)
{
    using T = TypeParam;
    const T tol = T(1e-6);
    const Matrix<T> ex = V3<T>(T(1), T(0), T(0));
    const Matrix<T> ey = V3<T>(T(0), T(1), T(0));
    const Matrix<T> ez = V3<T>(T(0), T(0), T(1));
    const Matrix<T> zero = V3<T>(T(0), T(0), T(0));

    auto chk = [&](const Matrix<T>& got, const Matrix<T>& want, const std::string& name) {
        STEP(name);
        EXPECT_TRUE(SameV(got, want, tol)) << "实测 " << FmtV(got);
    };

    // 右手系基元乘法表 —— 必须与 quaternion.cpp [C1]-[C5] 是同一张表
    chk(Cross(ex, ey), ez, "[C1] e_x × e_y = +e_z");
    chk(Cross(ey, ez), ex, "[C2] e_y × e_z = +e_x");
    chk(Cross(ez, ex), ey, "[C3] e_z × e_x = +e_y");
    // 注意: Matrix 没有一元负号 (core/matrix.hpp 只有二元 operator-),
    // 所以取负走标量乘 —— 也不为此去扩 core 的接口。
    chk(Cross(ey, ex), ez * T(-1), "[C4] e_y × e_x = -e_z (反向)");
    chk(Cross(ex, ez), ey * T(-1), "[C5] e_x × e_z = -e_y (反向)");

    const Matrix<T> a = SampleA<T>();
    const Matrix<T> b = SampleB<T>();

    // 反交换
    { STEP("[C6] 反交换: a×b + b×a = 0");
      EXPECT_TRUE(SameV(Cross(a, b) + Cross(b, a), zero, tol)); }

    // 与两个因子都正交 (手算点积, 不引 linalg::Dot —— 本文件保持只依赖 core)
    const Matrix<T> c = Cross(a, b);
    const T da = a(0, 0) * c(0, 0) + a(1, 0) * c(1, 0) + a(2, 0) * c(2, 0);
    const T db = b(0, 0) * c(0, 0) + b(1, 0) * c(1, 0) + b(2, 0) * c(2, 0);
    { STEP("[C7] (a×b) ⊥ a, b");
      EXPECT_TRUE(std::fabs(da) <= tol && std::fabs(db) <= tol)
          << "实测 " << static_cast<double>(da) << " / " << static_cast<double>(db); }

    // 平行 => 0 (整数运算在 float/double 下精确)
    const Matrix<T> p = V3<T>(T(1), T(2), T(3));
    const Matrix<T> q = V3<T>(T(2), T(4), T(6));
    { STEP("[C8] 平行向量叉积逐位精确 0");
      EXPECT_TRUE(SameV(Cross(p, q), zero, T(0))); }
}

// ---- [D] Cross 与 Hat 的一致性 (两条独立路径) -----------------------------

TYPED_TEST(So3Test, CrossVsHat)
{
    using T = TypeParam;
    const T tol = T(1e-6);
    const Matrix<T> a = SampleA<T>();
    const Matrix<T> b = SampleB<T>();

    // Cross (手写三分量) vs Hat(a) * b (走 Matrix 乘法) —— 两条完全不同的代码路径
    { STEP("[D1] Cross(a,b) = Hat(a)·b");
      EXPECT_TRUE(IsApprox(Cross(a, b), Hat(a) * b, tol))
          << "实测 " << FmtV(Cross(a, b)) << " vs " << FmtV(Hat(a) * b); }

    // 李括号: [Hat(a), Hat(b)] = Hat(a×b) —— **因子是 1, 不是 2**。
    //   (a=e_x, b=e_y 手验: 左边 [[0,-1,0],[1,0,0],[0,0,0]] = Hat(e_z)。)
    // 易混的点: vee(R - Rᵀ) 那个口径里出现的 1/2 是另一回事
    //   (θ = atan2(‖vee(R-Rᵀ)‖/2, (tr-1)/2)), 别把两个 2 记混。
    // 这条的意义: 它证明 Hat 是**李代数同构** (保括号), 不只是个线性双射。
    const Matrix<T> comm = Hat(a) * Hat(b) - Hat(b) * Hat(a);   // [Hat(a), Hat(b)]
    const Matrix<T> want = Hat(Cross(a, b));
    { STEP("[D2] 李括号 [Hat(a),Hat(b)] = Hat(a×b)");
      EXPECT_TRUE(IsApprox(comm, want, tol * T(10)))
          << "‖[Hat(a),Hat(b)] − Hat(a×b)‖∞ = " << NormInf(comm - want); }
}

// ---- [E] Hat 立方恒等式 (Rodrigues 闭式的唯一关键) -----------------------

// ⚠ 盲区 (2026-10-20 注入实测): 这一节**抓不到 Hat 的整体符号翻转**。
//   把 Hat 改成 -Hat (等价左手系), [E1][E2][E3] 仍然全绿 ——
//   因为两边都是奇次 (左边 K³, 右边 |u|²K), 整体翻号互相抵消。
//   ⇒ **符号由 [A3] (六位置逐位) 和 [D1] (Cross 对拍) 钉, 不由 [E] 钉。**
//   同理, 往返类判据 (Vee∘Hat = id) 也不免疫: 见 quaternion.hpp:110 那条
//   「叉积符号写反, 往返测试抓不到」。
TYPED_TEST(So3Test, HatCube)
{
    using T = TypeParam;
    const T tol = T(1e-6);

    // 单位轴: Hat(u)³ = -Hat(u)
    const Matrix<T> ey = V3<T>(T(0), T(1), T(0));
    const Matrix<T> K = Hat(ey);
    { STEP("[E1] 单位轴: Hat(u)³ = -Hat(u)");
      EXPECT_TRUE(IsApprox((K * K) * K, K * T(-1), tol))
          << "‖K³ + K‖∞ = " << NormInf((K * K) * K + K); }

    // 一般轴: Hat(u)³ = -|u|²·Hat(u) —— 把尺度带上, 验证不是碰巧
    const Matrix<T> u = V3<T>(T(1), T(2), T(3));        // |u|² = 14
    const Matrix<T> M = Hat(u);
    const T n2 = T(1) * T(1) + T(2) * T(2) + T(3) * T(3);
    { STEP("[E2] 一般轴: Hat(u)³ = -|u|²·Hat(u)  (|u|²=14)");
      EXPECT_TRUE(IsApprox((M * M) * M, M * (-n2), tol * T(10)))
          << "‖M³ + |u|²M‖∞ = " << NormInf((M * M) * M + M * n2); }

    // 零轴退化: Hat(0)³ = 0 = -0·Hat(0)
    const Matrix<T> z = Hat(V3<T>(T(0), T(0), T(0)));
    { STEP("[E3] 零轴: Hat(0)³ = 0 (逐位精确)");
      EXPECT_TRUE(IsApprox((z * z) * z, z, T(0))); }
}

// ---- [F] 与四元数对拍: 纯四元数积的向量部分 = 叉积 ------------------------

TYPED_TEST(So3Test, CrossVsQuaternion)
{
    using T = TypeParam;
    using Q = Quaternion<T>;
    const T tol = T(1e-6);
    const Matrix<T> a = SampleA<T>();
    const Matrix<T> b = SampleB<T>();

    // 纯四元数 (0, a) · (0, b) = ( -a·b , a×b )
    const Q p = Q(T(0), a(0, 0), a(1, 0), a(2, 0))
              * Q(T(0), b(0, 0), b(1, 0), b(2, 0));

    const T dot = a(0, 0) * b(0, 0) + a(1, 0) * b(1, 0) + a(2, 0) * b(2, 0);
    { STEP("[F1] 纯四元数积的实部 = -(a·b)");
      EXPECT_TRUE(std::fabs(p.w + dot) <= tol)
          << "实测 " << static_cast<double>(p.w) << " vs "
          << static_cast<double>(-dot); }

    const Matrix<T> cv = Cross(a, b);
    { STEP("[F2] 纯四元数积的向量部分 = a×b");
      EXPECT_TRUE(SameV(V3<T>(p.x, p.y, p.z), cv, tol))
          << "实测 " << FmtV(V3<T>(p.x, p.y, p.z)) << " vs " << FmtV(cv); }

    // 推论: 两个平行向量的纯四元数**交换** (叉积为 0) —— 四元数非交换性的来源
    // 正是叉积, 这条把 quaternion.hpp 的注释 (第 110 行) 变成可执行判据
    const Matrix<T> u = V3<T>(T(2), T(-1), T(0.5));
    const Matrix<T> u2 = V3<T>(T(-6), T(3), T(-1.5));      // = -3u
    const Q a1 = Q(T(0), u(0, 0), u(1, 0), u(2, 0));
    const Q b1 = Q(T(0), u2(0, 0), u2(1, 0), u2(2, 0));
    { STEP("[F3] 平行纯四元数交换 (叉积=0 ⇒ 交换)");
      EXPECT_TRUE(SameQ(a1 * b1, b1 * a1, tol)); }
}

// ---- [G] 形状契约: 传错形状必须 throw --------------------------------------

TEST(So3, ShapeContract)
{
    using D = double;

    { STEP("[G1] Hat: 2x1 非法 (要 3x1 列向量)");
      EXPECT_EXACT_THROW((void)Hat(Matrix<D>(2, 1)), std::invalid_argument); }
    { STEP("[G2] Hat: 1x3 非法 (行向量不接受, D19)");
      EXPECT_EXACT_THROW((void)Hat(Matrix<D>(1, 3)), std::invalid_argument); }
    { STEP("[G3] Vee: 3x1 非法 (要 3x3)");
      EXPECT_EXACT_THROW((void)Vee(Matrix<D>(3, 1)), std::invalid_argument); }
    { STEP("[G4] Cross: 第二个入参 2x1 非法");
      EXPECT_EXACT_THROW((void)Cross(V3<D>(1, 0, 0), Matrix<D>(2, 1)), std::invalid_argument); }
    { STEP("[G5] Cross: 第一个入参 2x1 非法");
      EXPECT_EXACT_THROW((void)Cross(Matrix<D>(2, 1), V3<D>(1, 0, 0)), std::invalid_argument); }
}