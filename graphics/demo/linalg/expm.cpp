// expm.cpp -- expm 验收: 指数映射 (两条路线 + 对易关系)
//
// 判据: 每段的独立判据见本文件各段 [A]–[E] 的注释。
//   对易关系数值验证: so(3) 生成元 exp(Ât) 作用向量 = 旋转;  [Jx, Jy] = Jz。
//
// ---------------------------------------------------------------------------
// 这份文件的状态
// ---------------------------------------------------------------------------
// 已落地: 22 项编号判据 (A1–A5 / B1–B4 / C1–C3 / D1–D6 / E1–E4) 全部为真判据;
//   断言数是 25 条 —— [E1] 对 4 个规模各断一次 (同一编号, 4 条断言)。
//   [E] 是 V2.4: expm_pade (阶段 2 第 3 条) 与 Taylor 版互为对照基线。
// 纪律 (判据先行): 阈值与"多少算对"先写在注释里, 再写代码。
//
// 五段的分工:
//   [A] 刀 1  expm 本体      幂零恒等式 (判精确) + expm(0) + 守卫 throw
//   [B] 刀 2  so(3) 对易     [Ji,Jj] = ε_ijk Jk  (判精确 0)
//   [C] 刀 3  三线对拍       级数 / 谱路线 / Rodrigues 闭式 互相对拍
//   [D] 几何 + 超强判据      R·u=u, tr R=1+2cosθ, det R=1, RᵀR=I; det(exp X)=exp(tr X)
//   [E] Padé 与 Taylor 互为对照基线 (必须相对判据)
//
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 判据编号 [A1]…[E4] 原样保留, 通过 SCOPED_TRACE 携带; 每段一个 TEST。
//   · 迁移前 25 项断言 -> 迁移后 25 条 EXPECT (一条没合并)。
//   · [A4]/[A5] 老写法 catch (const std::exception&), 现在判**精确类型**
//     std::invalid_argument (契约见 linalg/expm.hpp「一、API 契约」)。
//     判精确类型更严: 派生类冒充会假绿, 见 gtest_common.hpp 文件头。
//
// ---------------------------------------------------------------------------
// 独立标尺 (另一份实现: /tmp/probe_v16/probe.cpp, 由 agent 逐条校核)
// ---------------------------------------------------------------------------
//   [Jx,Jy]−Jz = 0 (精确)      | 幂零 expm(N) = I+N (精确)   | 非方阵/空矩阵 throw
//   级数 vs Rodrigues ≤ 1.1e-15 | 级数 vs 谱 ≤ 1.2e-15        | 谱 vs Rodrigues ≤ 6.1e-16
//   R·u−u ≤ 1.3e-15            | tr R −(1+2cosθ) ≤ 1.7e-15   | det R − 1 ≤ 1.2e-15
//   det(exp X)/exp(tr X) − 1 ≤ 1.0e-15 (随机 3×3)            | exp(2π[u]×) − I = 6.7e-16
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                    | 边界 / 病态输入                    | 专抓什么错
//  ----|-------------------------|------------------------------------|--------------------------
//   A  | **精确 0** (不加容差)   | 幂零指数 1 的上三角 (只留一个非对角)| 幂零判精确却被人加容差
//      | throw 类型精确          | 2×3 非方阵 / 0×0 空矩阵            | 静默返回空矩阵
//   B  | **精确 0**              | Jx/Jy/Jz 三基 (整型系数)           | 生成元符号 / 循环次序写反
//   C  | 1e-12 绝对 (实测 ≤1.2e-15)| θ ∈ {0,0.3,1,2,3,−0.9,2π}        | 级数项数不足 / 谱路线 e^{−iλ} 翻号
//   D  | 1e-12 绝对              | 绕 z 轴 θ=1; 随机 3×3 |a|<=1     | 不是旋转却仍往"像旋转"的方向凑
//   E  | **相对** 1e-8; E3/E4 见注 | 5×5 |a|<=100 (元素量级 e^100)、复 4×4 | **Padé 的 D 奇次项翻号**
//      |                         | 2×2 幂零 (判逐位精确)              | (X = I, 后段平方仍是 I)
//
// 关键阈值的来历:
//   · **[A]/[B] 判精确 0**: 代数恒等式 (幂零 / 对易) 在浮点下也能精确 ——
//     这是少见的"不该加容差"的场合, 加了容差就抓不到整体翻号。
//   · **[C]/[D] 的 1e-12**: 实测 ≤1.2e-15, 别往 1e-16 卡。
//   · **[E1]/[E2] 必须相对**: 元素量级 e^{‖A‖}, ‖A‖≤100 时绝对差是 1e+29 (假 FAIL),
//     相对差只有 1.5e-13 (numeric-contract.md §4.3④)。
//   · **[E4] 判 `d == 0.0`**: 与刀 1 同一条纪律 (幂零 ⇒ 逐位精确), 别加容差。
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **Padé 的 D 奇次项不取负** (X3: D = Σ c_k(−A)^k 变成 Σ c_k A^k):
//     N == D ⇒ X = D⁻¹N = I, 后面的平方仍是 I ⇒ expm_pade 恒返回 I。
//     **只有 [E] 段会红** ([E1] 4 条 + [E2] + [E3] + [E4]), 刀 1–4 全绿 ——
//     这正是"Padé 必须有独立判据"的原因, 也是本文件最容易漏掉的整段错误。
//  ② **Taylor 版去掉 scaling** (式 X1 的 s 不再除): ‖A‖ 小时看不出来,
//     只有 [E1] 的 5×5 |a|<=100 会红; 刀 1–4 的 ‖A‖ 都 <= 2π, 可能仍全绿。
//  ③ **[C] 的三线对拍抓不到"三线一起错"**: 级数 / 谱 / Rodrigues 都是本文件写的,
//     若 Rodrigues 与级数同错一个符号, [C3] 反而更绿 —— 所以 [B] 的基元对易子必须保留,
//     它是不依赖这三条路线的第四条尺 (整型系数, 判精确)。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh -d    # 全量 + 逐用例明细
//   ./graphics/build/demo/linalg/expm
//   ./graphics/build/demo/linalg/expm --gtest_filter='*Pade*'
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/eigen.hpp"
#include "linalg/expm.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace
{

using D  = double;
using C  = Complex<double>;
using MD = Matrix<D>;
using MC = Matrix<C>;

// core/common.hpp 里已有一个 `constexpr float PI`(精度只有 float) —— 这里要 double 版,
// 名字避开它, 否则重定义报错 (实测)。
constexpr D K_PI = 3.14159265358979323846;

// 超越函数对拍阈值
// 代数恒等式 (幂零/对易) 不用它, 判精确 0。
constexpr D TOL = 1e-12;

std::string Fmt(D v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

// ---- 手工小工具 (库不提供 Det/Trace; 只用到 3x3) --------------------------

MD Skew(D x, D y, D z)
{
    MD S = MD::Zero(3, 3);
    S(0, 1) = -z; S(0, 2) = y;
    S(1, 0) = z;  S(1, 2) = -x;
    S(2, 0) = -y; S(2, 1) = x;
    return S;
}
D Tr(const MD& A) { D s = 0; for (int i = 0; i < A.Rows(); i++) s += A(i, i); return s; }
D Det3(const MD& A)
{
    return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
         - A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
         + A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
}

// Rodrigues 闭式 (标尺, 不是被判对象): R = cosθ·I + (1−cosθ)·u uᵀ + sinθ·[u]×
// (原先是 float 旧草稿 src/view/rotate.hpp 里的注释；该文件与整个 src/view/ 已于 2026-10-26 删除)
MD Rodrigues(const D u[3], D th)
{
    const D c = std::cos(th), s = std::sin(th);
    MD R = MD::Identity(3) * c;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R(i, j) = R(i, j) + (1.0 - c) * u[i] * u[j];
    return R + Skew(u[0], u[1], u[2]) * s;
}

// 谱路线 (第三条路, 靠 eig_symmetric 而不是级数):
//   A = θ[u]× 实反对称 ⇒ H = iA 是 Hermitian ⇒ Q Λ Qᴴ (Λ 实) ⇒ exp(A) = Q·diag(e^{−iλ})·Qᴴ
//   本库只有 Hermitian 特征分解, 所以必须绕这一步 —— 不必写复数通用特征值分解。
MC SpectralExp(const MD& A)
{
    const int n = A.Rows();
    MC H(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) H(i, j) = C(0.0, A(i, j));    // i * a
    const Eigendecomposition<C> e = eig_symmetric(H);
    MC Dg = MC::Zero(n, n);
    for (int k = 0; k < n; k++) {
        const D lam = e.Lamda(k, k).real;
        Dg(k, k) = C(std::cos(lam), -std::sin(lam));             // exp(−i·λ)
    }
    return e.Q * Dg * e.Q.Adjoint();
}

}  // namespace

// 判据编号进 SCOPED_TRACE; 标签与断言必须在同一个块内 (见 docs/pitfalls.md §十二)
#define STEP(label) SCOPED_TRACE(label)

// ===========================================================================
// [A] 刀 1: expm 本体
// ===========================================================================

TEST(Expm, BasicNilpotent)
{
    // 幂零前提: 严格上三角 + 只留一个非对角元 ⇒ N² 逐位为 0。
    // 坑: [[0,a,b],[0,0,c],[0,0,0]] 的 A²(0,2) = a·c ≠ 0, 幂零指数是 3, 前提不成立。
    MD N = MD::Zero(3, 3);
    N(0, 1) = 0.7;
    const MD Id = MD::Identity(3);
    const MD eN = expm(N);

    { STEP("[A1] 前提核对 MaxAbs2(N²) == 0.0 (精确)");
      EXPECT_EQ(MaxAbs2(N * N), 0.0) << "实测 " << Fmt(MaxAbs2(N * N)); }

    { STEP("[A2] 幂零 expm(N) 逐位 == I+N (精确)");
      EXPECT_EQ(MaxAbs2(eN - (Id + N)), 0.0) << "实测 " << Fmt(MaxAbs2(eN - (Id + N))); }

    { STEP("[A3] expm(0) == I (精确)");
      EXPECT_EQ(MaxAbs2(expm(MD::Zero(3, 3)) - Id), 0.0)
          << "实测 " << Fmt(MaxAbs2(expm(MD::Zero(3, 3)) - Id)); }

    { STEP("[A4] 非方阵 throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)expm(MD::Zero(2, 3)), std::invalid_argument); }

    { STEP("[A5] 空矩阵 throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)expm(MD(0, 0)), std::invalid_argument); }
}

// ===========================================================================
// [B] 刀 2: so(3) 三基与对易关系
// ===========================================================================

TEST(Expm, So3Commutator)
{
    MD Jx = MD::Zero(3, 3); Jx(1, 2) = -1; Jx(2, 1) = 1;
    MD Jy = MD::Zero(3, 3); Jy(0, 2) =  1; Jy(2, 0) = -1;
    MD Jz = MD::Zero(3, 3); Jz(0, 1) = -1; Jz(1, 0) =  1;

    const MD cyc[3] = { Jx * Jy - Jy * Jx, Jy * Jz - Jz * Jy, Jz * Jx - Jx * Jz };
    const MD tgt[3] = { Jz, Jx, Jy };
    const char* nm[3] = { "[Jx,Jy] - Jz", "[Jy,Jz] - Jx", "[Jz,Jx] - Jy" };

    for (int k = 0; k < 3; k++) {
        STEP(std::string("[B") + char('1' + k) + "] " + nm[k] + " == 0 (精确)");
        EXPECT_EQ(MaxAbs2(cyc[k] - tgt[k]), 0.0) << "实测 " << Fmt(MaxAbs2(cyc[k] - tgt[k]));
    }

    { STEP("[B4] Jx 反对称 Jx + Jxᵀ == 0 (精确)");
      EXPECT_EQ(MaxAbs2(Jx + Jx.Transpose()), 0.0)
          << "实测 " << Fmt(MaxAbs2(Jx + Jx.Transpose())); }
}

// ===========================================================================
// [C] 刀 3: 三线对拍 (级数 / 谱 / Rodrigues 闭式)
// ===========================================================================

TEST(Expm, ThreeWayCrossCheck)
{
    D u[3] = { 1.0, 2.0, 3.0 };
    const D nn = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    for (int i = 0; i < 3; i++) u[i] /= nn;
    const MD K = Skew(u[0], u[1], u[2]);

    // 含负角与整圈: 0 检查恒等, 2π 检查"绕一圈回到 I"
    const D thetas[7] = { 0.0, 0.3, 1.0, 2.0, 3.0, -0.9, 2.0 * K_PI };

    D m_ser_rod = 0.0, m_ser_spec = 0.0, m_spec_rod = 0.0;
    for (D th : thetas) {
        const MD series = expm(K * th);
        const MD rod    = Rodrigues(u, th);
        const MC spec   = SpectralExp(K * th);
        MD spec_r(3, 3);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) spec_r(i, j) = spec(i, j).real;
        m_ser_rod  = std::max(m_ser_rod,  NormInf(series - rod));
        m_ser_spec = std::max(m_ser_spec, NormInf(series - spec_r));
        m_spec_rod = std::max(m_spec_rod, NormInf(spec_r - rod));
    }

    { STEP("[C1] max‖级数 − Rodrigues‖∞ < 1e-12");
      EXPECT_LT(m_ser_rod, TOL) << "实测 " << Fmt(m_ser_rod); }

    { STEP("[C2] max‖级数 − 谱路线‖∞ < 1e-12");
      EXPECT_LT(m_ser_spec, TOL) << "实测 " << Fmt(m_ser_spec); }

    { STEP("[C3] max‖谱 − Rodrigues‖∞ < 1e-12");
      EXPECT_LT(m_spec_rod, TOL) << "实测 " << Fmt(m_spec_rod); }
}

// ===========================================================================
// [D] 几何判据 + 超强判据
// ===========================================================================

TEST(Expm, GeometryAndDetIdentity)
{
    D u[3] = { 0.0, 0.0, 1.0 };
    const D th = 1.0;
    const MD R = expm(Skew(u[0], u[1], u[2]) * th);
    const MD ucol = MD({ { u[0] }, { u[1] }, { u[2] } });

    { STEP("[D1] R·u − u < 1e-12");
      EXPECT_LT(NormInf(R * ucol - ucol), TOL) << "实测 " << Fmt(NormInf(R * ucol - ucol)); }

    { STEP("[D2] tr R − (1+2cosθ) < 1e-12");
      EXPECT_LT(std::fabs(Tr(R) - (1.0 + 2.0 * std::cos(th))), TOL)
          << "实测 " << Fmt(std::fabs(Tr(R) - (1.0 + 2.0 * std::cos(th)))); }

    { STEP("[D3] det R − 1 < 1e-12");
      EXPECT_LT(std::fabs(Det3(R) - 1.0), TOL)
          << "实测 " << Fmt(std::fabs(Det3(R) - 1.0)); }

    { STEP("[D4] RᵀR − I < 1e-12");
      EXPECT_LT(NormInf(R * R.Transpose() - MD::Identity(3)), TOL)
          << "实测 " << Fmt(NormInf(R * R.Transpose() - MD::Identity(3))); }

    // 超强判据: det(exp X) = exp(tr X) (Jacobi 公式), 对任意方阵成立、与算法无关
    D worst = 0.0;
    for (int s = 1; s <= 4; s++) {
        const MD X = RandMatrix<D>(3, 3, 100u + s, -1.0, 1.0);
        const D lhs = Det3(expm(X));
        const D rhs = std::exp(Tr(X));
        worst = std::max(worst, std::fabs(lhs - rhs) / std::fabs(rhs));
    }
    { STEP("[D5] max|det(exp X)/exp(tr X) − 1| < 1e-12");
      EXPECT_LT(worst, TOL) << "实测 " << Fmt(worst); }

    // 刀 4 前哨: exp(2π[u]×) = I  (exp 不单射)
    const MD Rp = expm(Skew(u[0], u[1], u[2]) * (2.0 * K_PI));
    { STEP("[D6] exp(2π[u]×) − I < 1e-12");
      EXPECT_LT(NormInf(Rp - MD::Identity(3)), TOL)
          << "实测 " << Fmt(NormInf(Rp - MD::Identity(3))); }
}

// ===========================================================================
// [E] V2.4: Padé 版 (expm_pade) —— 与 Taylor 版互为对照基线
// ===========================================================================
// 判据 (阈值): expm(A) 与 Taylor 级数版差 < 1e-8。
// **必须用相对判据**: ‖A‖ 大时元素量级是 e^{‖A‖}, 绝对差没有意义 ——
// 实测 |a|<=100 的 5x5, 绝对差 1.1e+29 会假报 FAIL, 相对差只有 1.5e-13。
// (老坑, 见 numeric-contract.md §4.3④「expm 判据必须相对」。)

TEST(Expm, PadeVsTaylor)
{
    struct Case { int n; double scale; const char* tag; };
    const Case cases[] = {{3, 1.0, "实 3x3 |a|<=1  "}, {3, 10.0, "实 3x3 |a|<=10 "},
                          {5, 1.0, "实 5x5 |a|<=1  "}, {5, 100.0, "实 5x5 |a|<=100"}};
    for (const Case& c : cases) {
        const MD A  = RandMatrix<D>(c.n, c.n, 7u, -c.scale, c.scale);
        const MD E1 = expm(A), E2 = expm_pade(A);
        STEP(std::string("[E1] ") + c.tag + " ||Δtaylor-pade||/||E|| < 1e-8");
        EXPECT_LT(NormInf(E1 - E2) / NormInf(E1), 1e-8)
            << "实测 " << Fmt(NormInf(E1 - E2) / NormInf(E1));
    }

    {   // 复数域
        const Matrix<C> Ac = RandMatrix<C>(4, 4, 11u, -2.0, 2.0);
        const Matrix<C> E1 = expm(Ac), E2 = expm_pade(Ac);
        const double rel = NormInf(E1 - E2) / NormInf(E1);
        STEP("[E2] 复 4x4 ||Δtaylor-pade||/||E|| < 1e-8");
        EXPECT_LT(rel, 1e-8) << "实测 " << Fmt(rel);
    }

    {   // so(3): Padé 也必须落回 Rodrigues 这条独立闭式
        const double th = 1.3;
        const MD K = Skew(0.6, 0.0, 0.8);
        const MD R = MD::Identity(3) + K * std::sin(th) + (K * K) * (1.0 - std::cos(th));
        const double d = NormInf(expm_pade(K * th) - R);
        STEP("[E3] Padé 落回 Rodrigues 闭式 < 1e-12");
        EXPECT_LT(d, TOL) << "实测 " << Fmt(d);
    }

    {   // 幂零: 逐位精确 (与刀 1 同一条纪律 —— 别加容差)
        MD N(2, 2); N(0, 1) = 7.0;
        const double d = NormInf(expm_pade(N) - (MD::Identity(2) + N));
        STEP("[E4] 幂零 expm_pade(N) = I+N **逐位精确**");
        EXPECT_EQ(d, 0.0) << "实测 " << Fmt(d);
    }
}