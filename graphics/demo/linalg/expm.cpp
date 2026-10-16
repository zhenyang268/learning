// expm.cpp -- V1.6 验收: 指数映射 expm  (16 项判据全绿, ctest 10/10)
//
// plan.md V1.6 原文:
//   对易关系数值验证: so(3) 生成元 exp(Ât) 作用向量 = 旋转;  [Jx, Jy] = Jz。
// 刀序与每刀的独立判据见 learning-todo/lie-plan.md (刀1 本体 / 刀2 对易 / 刀3 三线对拍)。
//
// ---------------------------------------------------------------------------
// 这份文件的状态
// ---------------------------------------------------------------------------
// 已落地: 20 项判据 (A1–A5 / B1–B4 / C1–C3 / D1–D6 / E1–E4) 全部为真判据。
//   [E] 是 V2.4: expm_pade (阶段 2 第 3 条) 与 Taylor 版互为对照基线。
// 纪律 (matrix.md §10.9 判据先行): 阈值与"多少算对"先写在注释里, 再写代码。
//
// 四段的分工:
//   [A] 刀 1  expm 本体      幂零恒等式 (判精确) + expm(0) + 守卫 throw
//   [B] 刀 2  so(3) 对易     [Ji,Jj] = ε_ijk Jk  (判精确 0)
//   [C] 刀 3  三线对拍       级数 / 谱路线 / Rodrigues 闭式 互相对拍
//   [D] 几何 + 超强判据      R·u=u, tr R=1+2cosθ, det R=1, RᵀR=I; det(exp X)=exp(tr X)
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
// 阈值纪律 (踩过的坑, 照做)
// ---------------------------------------------------------------------------
//   · **代数恒等式判"精确"**: 对易子、幂零 — `MaxAbs2(...) == 0.0`, 别加容差
//     (这是少见的"浮点也能要求精确"的场合, 见 lie-plan §4 刀1 补充)。
//   · **超越函数对拍判 1e-12**: 实测 ≤1.2e-15, 别往 1e-16 卡。
//   · 一切比较走 `Abs2` / `NormInf`, 不写元素 `==`。
//   · 复数域的标量乘**必须**给 `Complex`: `M * C(t, 0.0)` —— `Complex` 无单参构造,
//     `M * 2.0` 直接编译不过 (探针实测)。
//
// 运行:
//   ./graphics/run_tests.sh -d    # 全量 + 逐用例明细
//   ./graphics/build/demo/linalg/expm

#include "../common.hpp"
#include "linalg/eigen.hpp"
#include "linalg/expm.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

static Checker ck;

using D  = double;
using C  = Complex<double>;
using MD = Matrix<D>;
using MC = Matrix<C>;

// core/common.hpp 里已有一个 `constexpr float PI`(精度只有 float) —— 这里要 double 版,
// 名字避开它, 否则重定义报错 (实测)。
constexpr D K_PI = 3.14159265358979323846;

// 超越函数对拍阈值 (plan.md V1.6 未给数, 按实测 ≤1.2e-15 取 1e-12 留三个量级余量)
// 代数恒等式 (幂零/对易) 不用它, 判精确 0。
constexpr D TOL = 1e-12;

static std::string Fmt(D v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}
static std::string Bool(bool b) { return b ? "true" : "false"; }

// ---- 手工小工具 (库不提供 Det/Trace; 只用到 3x3) --------------------------

static MD Skew(D x, D y, D z)
{
    MD S = MD::Zero(3, 3);
    S(0, 1) = -z; S(0, 2) = y;
    S(1, 0) = z;  S(1, 2) = -x;
    S(2, 0) = -y; S(2, 1) = x;
    return S;
}
static D Tr(const MD& A) { D s = 0; for (int i = 0; i < A.Rows(); i++) s += A(i, i); return s; }
static D Det3(const MD& A)
{
    return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
         - A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
         + A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
}

// Rodrigues 闭式 (标尺, 不是被判对象): R = cosθ·I + (1−cosθ)·u uᵀ + sinθ·[u]×
// (源起 src/view/rotate.hpp 的注释; 那里是 float 旧草稿, 不要 include)
static MD Rodrigues(const D u[3], D th)
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
static MC SpectralExp(const MD& A)
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

// ===========================================================================
// [A] 刀 1: expm 本体
// ===========================================================================
static void A_expm_basic()
{
    std::printf("\n[A] 刀1 expm 本体\n");

    // 幂零前提: 严格上三角 + 只留一个非对角元 ⇒ N² 逐位为 0。
    // 坑: [[0,a,b],[0,0,c],[0,0,0]] 的 A²(0,2) = a·c ≠ 0, 幂零指数是 3, 前提不成立。
    MD N = MD::Zero(3, 3);
    N(0, 1) = 0.7;
    const MD Id = MD::Identity(3);
    const MD eN = expm(N);

    ck.check(MaxAbs2(N * N) == 0.0,
             "[A1] 前提核对 MaxAbs2(N²) == 0.0 (精确)      实测 " + Fmt(MaxAbs2(N * N)));
    ck.check(MaxAbs2(eN - (Id + N)) == 0.0,
             "[A2] 幂零 expm(N) 逐位 == I+N (精确)         实测 " + Fmt(MaxAbs2(eN - (Id + N))));
    ck.check(MaxAbs2(expm(MD::Zero(3, 3)) - Id) == 0.0,
             "[A3] expm(0) == I (精确)                     实测 " + Fmt(MaxAbs2(expm(MD::Zero(3, 3)) - Id)));

    bool throw_sq = false;
    try { expm(MD::Zero(2, 3)); } catch (const std::exception&) { throw_sq = true; }
    ck.check(throw_sq, "[A4] 非方阵 throw                            实测 " + Bool(throw_sq));

    bool throw_empty = false;
    try { expm(MD(0, 0)); } catch (const std::exception&) { throw_empty = true; }
    ck.check(throw_empty, "[A5] 空矩阵 throw                            实测 " + Bool(throw_empty));
}

// ===========================================================================
// [B] 刀 2: so(3) 三基与对易关系
// ===========================================================================
static void B_so3_commutator()
{
    std::printf("\n[B] 刀2 so(3) 对易关系\n");

    MD Jx = MD::Zero(3, 3); Jx(1, 2) = -1; Jx(2, 1) = 1;
    MD Jy = MD::Zero(3, 3); Jy(0, 2) =  1; Jy(2, 0) = -1;
    MD Jz = MD::Zero(3, 3); Jz(0, 1) = -1; Jz(1, 0) =  1;

    const MD cyc[3] = { Jx * Jy - Jy * Jx, Jy * Jz - Jz * Jy, Jz * Jx - Jx * Jz };
    const MD tgt[3] = { Jz, Jx, Jy };
    const char* nm[3] = { "[Jx,Jy] - Jz", "[Jy,Jz] - Jx", "[Jz,Jx] - Jy" };

    for (int k = 0; k < 3; k++) {
        ck.check(MaxAbs2(cyc[k] - tgt[k]) == 0.0,
                 std::string("[B") + char('1' + k) + "] " + nm[k] +
                 " == 0 (精确)         实测 " + Fmt(MaxAbs2(cyc[k] - tgt[k])));
    }
    ck.check(MaxAbs2(Jx + Jx.Transpose()) == 0.0,
             "[B4] Jx 反对称 Jx + Jxᵀ == 0 (精确)        实测 " + Fmt(MaxAbs2(Jx + Jx.Transpose())));
}

// ===========================================================================
// [C] 刀 3: 三线对拍 (级数 / 谱 / Rodrigues 闭式)
// ===========================================================================
static void C_three_way()
{
    std::printf("\n[C] 刀3 三线对拍 (u = (1,2,3)/√14)\n");

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

    ck.check(m_ser_rod < TOL, "[C1] max‖级数 − Rodrigues‖∞ < 1e-12    实测 " + Fmt(m_ser_rod));
    ck.check(m_ser_spec < TOL, "[C2] max‖级数 − 谱路线‖∞   < 1e-12    实测 " + Fmt(m_ser_spec));
    ck.check(m_spec_rod < TOL, "[C3] max‖谱   − Rodrigues‖∞ < 1e-12    实测 " + Fmt(m_spec_rod));
}

// ===========================================================================
// [D] 几何判据 + 超强判据
// ===========================================================================
static void D_geometry_and_det()
{
    std::printf("\n[D] 几何判据 + det(exp X) = exp(tr X)\n");

    D u[3] = { 0.0, 0.0, 1.0 };
    const D th = 1.0;
    const MD R = expm(Skew(u[0], u[1], u[2]) * th);
    const MD ucol = MD({ { u[0] }, { u[1] }, { u[2] } });

    ck.check(NormInf(R * ucol - ucol) < TOL,
             "[D1] R·u − u < 1e-12                        实测 " + Fmt(NormInf(R * ucol - ucol)));
    ck.check(std::fabs(Tr(R) - (1.0 + 2.0 * std::cos(th))) < TOL,
             "[D2] tr R − (1+2cosθ) < 1e-12               实测 " + Fmt(std::fabs(Tr(R) - (1.0 + 2.0 * std::cos(th)))));
    ck.check(std::fabs(Det3(R) - 1.0) < TOL,
             "[D3] det R − 1 < 1e-12                      实测 " + Fmt(std::fabs(Det3(R) - 1.0)));
    ck.check(NormInf(R * R.Transpose() - MD::Identity(3)) < TOL,
             "[D4] RᵀR − I < 1e-12                        实测 " + Fmt(NormInf(R * R.Transpose() - MD::Identity(3))));

    // 超强判据: det(exp X) = exp(tr X) (Jacobi 公式), 对任意方阵成立、与算法无关
    D worst = 0.0;
    for (int s = 1; s <= 4; s++) {
        const MD X = RandMatrix<D>(3, 3, 100u + s, -1.0, 1.0);
        const D lhs = Det3(expm(X));
        const D rhs = std::exp(Tr(X));
        worst = std::max(worst, std::fabs(lhs - rhs) / std::fabs(rhs));
    }
    ck.check(worst < TOL, "[D5] max|det(exp X)/exp(tr X) − 1| < 1e-12   实测 " + Fmt(worst));

    // 刀 4 前哨: exp(2π[u]×) = I  (exp 不单射)
    const MD Rp = expm(Skew(u[0], u[1], u[2]) * (2.0 * K_PI));
    ck.check(NormInf(Rp - MD::Identity(3)) < TOL,
             "[D6] exp(2π[u]×) − I < 1e-12                实测 " + Fmt(NormInf(Rp - MD::Identity(3))));
}

// ---------------------------------------------------------------------------
// [E] V2.4: Padé 版 (expm_pade) —— 与 Taylor 版互为对照基线
// ---------------------------------------------------------------------------
// 判据 (V2.4 原文): expm(A) 与 Taylor 级数版差 < 1e-8。
// **必须用相对判据**: ‖A‖ 大时元素量级是 e^{‖A‖}, 绝对差没有意义 ——
// 实测 |a|<=100 的 5x5, 绝对差 1.1e+29 会假报 FAIL, 相对差只有 1.5e-13。
// (老坑, 见 matrix.md §3.6「expm 判据必须相对」。)
static void E_pade_vs_taylor()
{
    std::printf("\n-- [E] V2.4: Padé 版 vs Taylor 版 (相对判据) --\n");

    struct Case { int n; double scale; const char* tag; };
    const Case cases[] = {{3, 1.0, "实 3x3 |a|<=1  "}, {3, 10.0, "实 3x3 |a|<=10 "},
                          {5, 1.0, "实 5x5 |a|<=1  "}, {5, 100.0, "实 5x5 |a|<=100"}};
    for (const Case& c : cases) {
        const MD A  = RandMatrix<D>(c.n, c.n, 7u, -c.scale, c.scale);
        const MD E1 = expm(A), E2 = expm_pade(A);
        const D rel = NormInf(E1 - E2) / NormInf(E1);
        ck.check(rel < 1e-8,
                 std::string("[E1] ") + c.tag + " ||Δtaylor-pade||/||E|| < 1e-8  实测 " + Fmt(rel));
    }

    {   // 复数域
        const Matrix<C> Ac = RandMatrix<C>(4, 4, 11u, -2.0, 2.0);
        const Matrix<C> E1 = expm(Ac), E2 = expm_pade(Ac);
        const double rel = NormInf(E1 - E2) / NormInf(E1);
        ck.check(rel < 1e-8, "[E2] 复 4x4 ||Δtaylor-pade||/||E|| < 1e-8          实测 " + Fmt(rel));
    }

    {   // so(3): Padé 也必须落回 Rodrigues 这条独立闭式
        const double th = 1.3;
        const MD K = Skew(0.6, 0.0, 0.8);
        const MD R = MD::Identity(3) + K * std::sin(th) + (K * K) * (1.0 - std::cos(th));
        const double d = NormInf(expm_pade(K * th) - R);
        ck.check(d < TOL, "[E3] Padé 落回 Rodrigues 闭式 < 1e-12            实测 " + Fmt(d));
    }

    {   // 幂零: 逐位精确 (与刀 1 同一条纪律 —— 别加容差)
        MD N(2, 2); N(0, 1) = 7.0;
        const double d = NormInf(expm_pade(N) - (MD::Identity(2) + N));
        ck.check(d == 0.0, "[E4] 幂零 expm_pade(N) = I+N **逐位精确**        实测 " + Fmt(d));
    }
}

int main()
{
    std::printf("== V1.6 expm 验收 (刀1 本体 / 刀2 对易 / 刀3 三线对拍 / 几何+超强判据 / V2.4 Padé) ==\n");
    A_expm_basic();
    B_so3_commutator();
    C_three_way();
    D_geometry_and_det();
    E_pade_vs_taylor();
    return CheckReport(ck.ok, ck.fail, "V1.6 expm");
}
