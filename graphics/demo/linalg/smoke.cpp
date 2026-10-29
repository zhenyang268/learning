// smoke.cpp -- V1.1 显式闸门 (S6 收口第 ① 项)
//
//   Matrix<double> 与 Matrix<Complex<double>> 使用**同一份** linalg 源码
//   零修改编译通过 (traits 设计正确性的硬验收)。
//
// ---------------------------------------------------------------------------
// 为什么 V1.2-V1.6 全过了, 还要单独补这一件
// ---------------------------------------------------------------------------
// 2026-10-14 的体检发现: 当时 528 行 linalg/*.hpp **从未被任何编译单元 include 过** ——
// "文件存在"与"编译器看过"是两件事, 而模板代码的 bug 只在**实例化**时暴露
// (只 #include 抓不到; 见 AGENTS.md §3.1「已落地」判定)。
// V1.2-V1.6 各自只覆盖了**自己用到的**那几个入口 (solve/inverse/cholesky/eigen/svd/qr/expm)。
// V1.1 要的证据不是"某个入口算得对", 而是"**每一个公开入口**都在 double 与 Complex 下
// 各实例化过一次"。两者的区别在: 将来某人给某入口加了新的模板分支、或某个入口失去全部调用方
// 之后再被改坏, 专项用例都发现不了, 本文件能。
//
// ---------------------------------------------------------------------------
// 两段, 各有独立作用 (不要合并)
// ---------------------------------------------------------------------------
//   [A] 实例化闸门 —— 对每个公开入口**取函数指针** (odr-use 强制实例化)。
//       刻意不调用: 目的是让"实例化"这件事独立于"数值对不对"。
//       新增公开入口时**必须**补一行 (EntryTable + GateCount 两处), 否则闸门静默失效。
//   [B] 最小数值闸门 —— 每个入口真的跑一次, 给一条轻量判据。
//       阈值比专项用例**松** (专项用例才是权威阈值); 这里只判"跑得通 + 方向对"。
//
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 迁移前 56 项断言 -> 迁移后 56 条 EXPECT (逐条对应, 一条没合并)。
//   · **[A] 段的手法原样保留**: `EntryTable` + `GateCount` + "取函数指针、刻意不调用"
//     一行没动 —— 它不是"某个模板跑两遍", 而是一张**登记表**。
//     [A1]/[A2] 仍是两条独立断言 (double / Complex 各数一次), **不改成 TYPED_TEST**
//     (改成 TYPED_TEST 会让 kEntryCount 与 EntryTable 的一一对应关系散在两个实例里)。
//   · [B] 段才是"同一模板跑两域" ⇒ NumericGate<T> 由 TYPED_TEST_SUITE 跑 double + Complex。
//   · 异常判据走 EXPECT_EXACT_THROW。
//
// ---------------------------------------------------------------------------
// 支线定案 (S6 遗留三条, 2026-10-16 全部结算)
// ---------------------------------------------------------------------------
//   · `Dot` (operator.hpp) 零调用方 —— 本文件成为它的**第一个调用方** ([B12])。
//     判据顺带把它和 `Abs2` 钉在一起: <v,v> = Σ|v_i|², 复数域必须共轭第一个参数。
//   · `MatrixType` —— **用户定案**: 它按用途划分 ("怎么用是上层的责任"), 已从 core/common.hpp
//     挪到应用层 `demo/common.hpp`; core 层只留按代数结构分的 `MatrixKind` (带 tol 判定)。
//   · `qr_inverse` (inverse.hpp) —— 已补 `A^{-1} = R^{-1}·Q^H` (S-App6), [B6d] 随之改为数值判据。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                  | 边界 / 病态输入               | 专抓什么错
//  ----|-----------------------|-------------------------------|--------------------------
//   A  | **入口计数 == 40**    | double / Complex 各一次        | 新增入口没登记 (闸门静默变小)
//   B  | 1e-10 (smoke 级)      | 4×4 良态; SPD = M·Mᴴ + 4I      | 入口跑不通 / 方向错
//      | 1e-8 (V1.2/V1.4 借用) | 负定 Hermitian (Gershgorin 判据) | SPD 档误判
//      | 精确 0 (B11a/B13)     | expm(0) / 容器出口              | 精确性退化
//
// **本文件的阈值一律比专项用例松** —— 它是闸门, 不是权威判据。
// 改这里的阈值不会提高可信度, 只会让闸门更晚报警。
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **闸门只保证"登记过的入口"被实例化, 不保证"全部入口"被登记** ——
//     历史事故两次: `expm_pade` 漏登记 (expm.hpp)、`so3.hpp` 整个文件漏登记。
//     注入对照: 从 EntryTable 删掉一行 + 不同步改 kEntryCount ⇒ [A1]/[A2] **仍绿**
//     (计数是 GateCount 数的, 不是 EntryTable 数的)。这条盲区只能用
//     test-map.md §4.1 的"接口全图"兜, 本文件兜不住。
//  ② 删掉一行 **且** 同步把 kEntryCount 减 1 ⇒ 全绿且闸门变小。**这是设计使然**。
//  ③ [B] 段阈值 1e-10 抓不到"数值略退化" (专项用例抓 1e-8/1e-12 那些才抓得到)。
//
// ---------------------------------------------------------------------------
// 范围说明
// ---------------------------------------------------------------------------
// 本闸门覆盖 core + linalg。autodiff 层 (jacobian.hpp 的 Apply/Derive/SumAll/MultAll)
// 由 `scalar_traits` / `compute_graph` / `dual_complex` 三个用例覆盖, 不在此重复。
//
// 运行:
//   ./graphics/run_tests.sh            # 全量 (本用例注册名: smoke)
//   ./graphics/build/demo/linalg/smoke # 单跑
//   ./graphics/run_tests.sh -w         # 顺带开 -Wall -Wextra 看实例化告警

#include "gtest_common.hpp"

#include "linalg/cholesky.hpp"
#include "linalg/eigen.hpp"
#include "linalg/expm.hpp"
#include "linalg/inverse.hpp"
#include "linalg/operator.hpp"
#include "linalg/qr.hpp"
#include "linalg/solve.hpp"
#include "linalg/svd.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using D = double;
using C = Complex<double>;

// [B] 段阈值: 一律比专项用例松。专项判据 (V1.2 1e-8 求逆残差、V1.4 1e-8 重构 ...) 才是权威。
constexpr D TOL_SMOKE = 1e-10;

std::string Fmt(D v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

}  // namespace

// ===========================================================================
// [A] 实例化闸门: 每个公开入口一行
//     `&f<T>` 取地址 => odr-use => 强制实例化该模板特化 (含函数体)。
//     static constexpr 数据成员在此处只声明, 真正触发实例化的是 GateCount 里的取用。
// ===========================================================================
template <typename T>
struct EntryTable
{
    using Mat    = Matrix<T>;
    using real_t = typename ScalarTraits<T>::real_t;

    // ---- 自由函数: solve.hpp ----
    static constexpr auto lu          = &LU<T>;
    static constexpr auto solve       = &Solve<T>;
    static constexpr auto swap_row    = &SwapRow<T>;
    static constexpr auto swap_mult   = &SwapMultipliers<T>;
    // ---- inverse.hpp (Inverse 有两个重载, 取地址必须显式消歧) ----
    static constexpr auto inv_auto    = static_cast<Mat (*)(const Mat&)>(&Inverse<T>);
    static constexpr auto inv_hint    = static_cast<Mat (*)(const Mat&, MatrixKind)>(&Inverse<T>);
    static constexpr auto classify    = &Classify<T>;
    static constexpr auto is_diag     = &is_diagonal<T>;
    static constexpr auto is_sym      = &is_symmetric<T>;
    static constexpr auto is_orth     = &is_orthogonal<T>;
    static constexpr auto is_spd      = &is_positive_definite<T>;
    static constexpr auto lu_inv      = &lu_inverse<T>;
    static constexpr auto diag_inv    = &diagonal_inverse<T>;
    static constexpr auto chol_inv    = &cholesky_inverse<T>;
    static constexpr auto qr_inv      = &qr_inverse<T>;
    // ---- cholesky.hpp ----
    static constexpr auto chol        = &cholesky<T>;
    static constexpr auto chol_solve  = &cholesky_solve<T>;
    // ---- eigen.hpp ----
    static constexpr auto eig_all     = &eig<T>;
    static constexpr auto eig_sym     = &eig_symmetric<T>;
    static constexpr auto jacobi_rot  = &jacobi_rotation<T>;
    // ---- svd.hpp ----
    static constexpr auto svd_all     = &svd<T>;
    // ---- qr.hpp (三个入口, 一个都不能漏 —— S6b 三个函数都出过算法层 bug) ----
    static constexpr auto qr_house    = &QR_Householder<T>;
    static constexpr auto qr_givens   = &QR_givens<T>;
    static constexpr auto qr_mgs      = &QR_Shimit<T>;
    // ---- expm.hpp ----
    static constexpr auto expm_all    = &expm<T>;
    // ---- operator.hpp ----
    static constexpr auto dot         = &Dot<T>;
    static constexpr auto abs2        = &Abs2<T>;

    // ---- Matrix<T> 成员 (静态工厂 + 基础出口; 重载成员如 operator()/at/Data 略, 它们每个用例都在用) ----
    static constexpr auto m_zero      = &Mat::Zero;
    static constexpr auto m_ones      = &Mat::Ones;
    static constexpr auto m_identity  = &Mat::Identity;
    static constexpr auto m_diagonal  = &Mat::Diagonal;
    static constexpr auto m_transpose = &Mat::Transpose;
    static constexpr auto m_adjoint   = &Mat::Adjoint;
    static constexpr auto m_cwise_mul = &Mat::cwiseProduct;
    static constexpr auto m_cwise_div = &Mat::cwiseDivision;
    static constexpr auto m_plus      = &Mat::operator+;
    static constexpr auto m_minus     = &Mat::operator-;
    static constexpr auto m_rows      = &Mat::Rows;
    static constexpr auto m_cols      = &Mat::Cols;
    static constexpr auto m_size      = &Mat::Size;
};

// 期望覆盖的入口条数 —— 与 EntryTable 一一对应, 改动其中之一必须同时改另一处。
// 这条计数本身就是闸门的一部分: 漏抄一行 (闸门静默变小) 会被它抓到。
static constexpr int kEntryCount = 40;   // 27 自由函数 + 13 成员

template <typename T>
int GateCount()
{
    using E = EntryTable<T>;
    int n = 0;
    const auto touch = [&n](auto fn) { if (fn != nullptr) n++; };

    touch(E::lu);          touch(E::solve);       touch(E::swap_row);    touch(E::swap_mult);
    touch(E::inv_auto);    touch(E::inv_hint);    touch(E::classify);
    touch(E::is_diag);     touch(E::is_sym);      touch(E::is_orth);     touch(E::is_spd);
    touch(E::lu_inv);      touch(E::diag_inv);    touch(E::chol_inv);    touch(E::qr_inv);
    touch(E::chol);        touch(E::chol_solve);
    touch(E::eig_all);     touch(E::eig_sym);     touch(E::jacobi_rot);
    touch(E::svd_all);
    touch(E::qr_house);    touch(E::qr_givens);   touch(E::qr_mgs);
    touch(E::expm_all);
    touch(E::dot);         touch(E::abs2);

    touch(E::m_zero);      touch(E::m_ones);      touch(E::m_identity);  touch(E::m_diagonal);
    touch(E::m_transpose); touch(E::m_adjoint);
    touch(E::m_cwise_mul); touch(E::m_cwise_div);
    touch(E::m_plus);      touch(E::m_minus);
    touch(E::m_rows);      touch(E::m_cols);      touch(E::m_size);

    return n;
}

namespace
{

// ===========================================================================
// [B] 最小数值闸门: 每个入口跑一次 + 一条轻量判据
// ===========================================================================
template <typename T>
void NumericGate(const char* tag)
{
    using real_t = typename ScalarTraits<T>::real_t;
    using Mat    = Matrix<T>;

    const real_t TB   = real_t(TOL_SMOKE);
    const real_t TB2  = TB * TB;
    const auto   R    = [](real_t v) { return ScalarTraits<T>::fromReal(v); };
    const auto   suit = [&tag](const char* id, const std::string& body) {
        return std::string("[") + id + "] " + tag + " " + body;
    };

    // ---- [B1] Solve: 残差 ‖A·x − b‖∞ / (1 + ‖b‖∞) ----
    const Mat A = RandMatrix<T>(4, 4, 1001u, -1.0, 1.0) + Mat::Identity(4) * R(real_t(4.0));
    const Mat b = RandMatrix<T>(4, 1, 1002u, -1.0, 1.0);
    const Mat x = Solve(A, b);
    const real_t res_solve = NormInf(A * x - b) / (real_t(1) + NormInf(b));
    { SCOPED_TRACE(suit("B1", "Solve 残差 < 1e-10"));
      EXPECT_LT(res_solve, TB) << Fmt(double(res_solve)); }

    // ---- [B2] LU: L 单位下三角, U 上三角 ----
    const LUResult<T> lu = LU(A);
    real_t unit_dev = real_t(0), low_dev = real_t(0);
    for (int i = 0; i < 4; i++)
        unit_dev = std::max(unit_dev, ScalarTraits<T>::abs2(lu.L(i, i) - ScalarTraits<T>::one()));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < i; j++)
            low_dev = std::max(low_dev, ScalarTraits<T>::abs2(lu.U(i, j)));
    { SCOPED_TRACE(suit("B2", "LU 结构 L对角=I / U下三角=0"));
      EXPECT_TRUE(std::sqrt(double(unit_dev)) < TOL_SMOKE && std::sqrt(double(low_dev)) < TOL_SMOKE)
          << "偏差 " << Fmt(std::sqrt(double(unit_dev))) << " / " << Fmt(std::sqrt(double(low_dev))); }

    // ---- [B3] Inverse 自动分发: ‖A·A⁻¹ − I‖∞ (V1.2 判据 1e-8) ----
    const Mat Ai = Inverse(A);
    const real_t res_inv = NormInf(A * Ai - Mat::Identity(4));
    { SCOPED_TRACE(suit("B3", "Inverse 残差 < 1e-8"));
      EXPECT_LT(res_inv, real_t(1e-8)) << Fmt(double(res_inv)); }

    // ---- [B4] Inverse(A, hint) 与自动分发同结果 ----
    const Mat Ai_h = Inverse(A, MatrixKind::Generic);
    { SCOPED_TRACE(suit("B4", "Inverse(hint=Generic) 与自动分发一致"));
      EXPECT_LE(MaxAbs2(Ai - Ai_h), TB2) << Fmt(double(NormInf(Ai - Ai_h))); }

    // ---- [B5] 结构探测: is_diagonal / is_symmetric / is_orthogonal + Classify ----
    const Mat M4 = RandMatrix<T>(4, 4, 1003u, -1.0, 1.0);
    const Mat SPD = M4 * M4.Adjoint() + Mat::Identity(4) * R(real_t(4.0));   // 必然 SPD/Hermitian 正定
    const Mat HERM = A + A.Adjoint();                                        // Hermitian
    // 负定 Hermitian: 对角 ~8−20 = −12, 非对角 ≤ 2, Gershgorin 圆全在左半轴 => 必非正定
    const Mat NEG = HERM - Mat::Identity(4) * R(real_t(20.0));
    const real_t tol4 = DefaultTol<T>(4);
    { SCOPED_TRACE(suit("B5a", "is_diagonal / is_symmetric / is_orthogonal / !is_positive_definite"));
      EXPECT_TRUE(is_diagonal(Mat::Identity(3), DefaultTol<T>(3)) &&
                  is_symmetric(HERM, tol4) && is_orthogonal(Mat::Identity(3), DefaultTol<T>(3)) &&
                  !is_positive_definite(NEG)); }
    { SCOPED_TRACE(suit("B5b", "Classify(I) == Diagonal"));
      EXPECT_EQ(Classify(Mat::Identity(3), DefaultTol<T>(3)), MatrixKind::Diagonal); }
    { SCOPED_TRACE(suit("B5c", "Classify(SPD) == SymmetricPositiveDefinite"));
      EXPECT_EQ(Classify(SPD, tol4), MatrixKind::SymmetricPositiveDefinite); }
    { SCOPED_TRACE(suit("B5d", "Classify(负定 Hermitian) == Symmetric (SPD 档不误判)"));
      EXPECT_EQ(Classify(NEG, tol4), MatrixKind::Symmetric); }

    // ---- [B6] 各分支求逆 ----
    const Mat Dg = Mat::Diagonal(std::vector<T>{ R(real_t(2)), R(real_t(-4)),
                                                 R(real_t(8)), R(real_t(16)) });
    { SCOPED_TRACE(suit("B6a", "diagonal_inverse 残差"));
      EXPECT_LT(NormInf(Dg * diagonal_inverse(Dg, DefaultTol<T>(4)) - Mat::Identity(4)), TB); }
    { SCOPED_TRACE(suit("B6b", "lu_inverse 与 Inverse 一致"));
      EXPECT_LE(MaxAbs2(lu_inverse(A) - Ai), TB2); }
    { SCOPED_TRACE(suit("B6c", "cholesky_inverse 残差"));
      EXPECT_LT(NormInf(SPD * cholesky_inverse(SPD) - Mat::Identity(4)), TB); }

    // ---- [B6d] qr_inverse (S-App6 补回代): 数值判据 + 与 lu_inverse 一致 + 奇异抛 ----
    const Mat Aqi = qr_inverse(A);
    const real_t res_qri = NormInf(A * Aqi - Mat::Identity(4));
    { SCOPED_TRACE(suit("B6d", "qr_inverse 残差 < 1e-8 且与 lu_inverse 一致"));
      EXPECT_TRUE(res_qri < real_t(1e-8) && MaxAbs2(Aqi - Ai) <= TB2) << Fmt(double(res_qri)); }
    { SCOPED_TRACE(suit("B6e", "qr_inverse(零矩阵, R 对角奇异) -> throw std::runtime_error"));
      EXPECT_EXACT_THROW((void)qr_inverse(Mat::Zero(3, 3)), std::runtime_error); }

    // ---- [B7] Cholesky: 重构 + 解方程 ----
    const Mat L = cholesky(SPD);
    { SCOPED_TRACE(suit("B7a", "cholesky 重构 L·Lᴴ − A"));
      EXPECT_LT(NormInf(L * L.Adjoint() - SPD), TB); }
    const Mat xs = cholesky_solve(SPD, b);
    { SCOPED_TRACE(suit("B7b", "cholesky_solve 残差"));
      EXPECT_LT(NormInf(SPD * xs - b) / (real_t(1) + NormInf(b)), TB); }

    // ---- [B8] QR 三入口: 重构 + 列正交 (thin 口径 Q: m×n) ----
    const Mat A43 = RandMatrix<T>(4, 3, 1004u, -1.0, 1.0);
    const auto qr_check = [&](const char* nm, const QRcomposition<T>& q) {
        const real_t rec  = NormInf(q.Q * q.R - A43);
        const real_t orth = NormInf(q.Q.Adjoint() * q.Q - Mat::Identity(3));
        SCOPED_TRACE(suit("B8", std::string(nm) + " 重构 / 列正交"));
        EXPECT_TRUE(rec < TB && orth < TB)
            << Fmt(double(rec)) << " / " << Fmt(double(orth));
    };
    qr_check("QR_Householder", QR_Householder(A43));
    qr_check("QR_givens", QR_givens(A43));
    qr_check("QR_Shimit", QR_Shimit(A43));

    // ---- [B9] 特征分解: 两条入口重构 + jacobi_rotation 的守恒量 ----
    const Eigendecomposition<T> ed = eig_symmetric(HERM);
    { SCOPED_TRACE(suit("B9a", "eig_symmetric 重构"));
      EXPECT_LT(NormInf(ed.Q * ed.Lamda * ed.Q.Adjoint() - HERM), TB); }
    const Eigendecomposition<T> e2 = eig(HERM);
    { SCOPED_TRACE(suit("B9b", "eig 重构"));
      EXPECT_LT(NormInf(e2.Q * e2.Lamda * e2.Q.Adjoint() - HERM), TB); }
    {
        // 2x2 Hermitian 子块 [[a, z], [conj z], b]] 一次旋转后: 特征值满足 迹 = a+b, 积 = ab − |z|²
        const real_t a = real_t(3), bb = real_t(-1);
        const T z = R(real_t(0.5));
        real_t c = 0, s = 0, a_new = 0, b_new = 0;
        T ph = ScalarTraits<T>::zero();
        jacobi_rotation<T>(a, bb, z, c, s, ph, a_new, b_new);
        const real_t rho2    = ScalarTraits<T>::abs2(z);
        const real_t tr_dev  = std::fabs(a_new + b_new - (a + bb));
        const real_t det_dev = std::fabs(a_new * b_new - (a * bb - rho2));
        SCOPED_TRACE(suit("B9", "jacobi_rotation 迹/行列式守恒"));
        EXPECT_TRUE(tr_dev < TB && det_dev < TB)
            << "偏差 " << Fmt(double(tr_dev)) << " / " << Fmt(double(det_dev));
    }

    // ---- [B10] SVD: 重构 + 奇异值非负降序 ----
    const SVDcomposition<T> sd = svd(A43);
    { SCOPED_TRACE(suit("B10a", "svd 重构"));
      EXPECT_LT(NormInf(sd.U * sd.S * sd.V.Adjoint() - A43), TB); }
    bool sv_ok = true;
    for (int i = 0; i < sd.S.Rows(); i++) {
        const real_t sv = ScalarTraits<T>::real(sd.S(i, i));
        if (sv < real_t(0)) sv_ok = false;
        if (i > 0 && sv > ScalarTraits<T>::real(sd.S(i - 1, i - 1))) sv_ok = false;
    }
    { SCOPED_TRACE(suit("B10b", "svd 奇异值非负且降序"));
      EXPECT_TRUE(sv_ok); }

    // ---- [B11] expm: expm(0) = I (精确) + 小范数下 expm(X)·expm(−X) = I ----
    { SCOPED_TRACE(suit("B11a", "expm(0) == I (精确)"));
      EXPECT_EQ(MaxAbs2(expm(Mat::Zero(4, 4)) - Mat::Identity(4)), real_t(0)); }
    const Mat X = RandMatrix<T>(3, 3, 1005u, -1.0, 1.0) * R(real_t(0.25));
    const real_t res_exp = NormInf(expm(X) * expm(X * R(real_t(-1))) - Mat::Identity(3));
    { SCOPED_TRACE(suit("B11b", "expm(X)·expm(−X) 残差"));
      EXPECT_LT(res_exp, TB) << Fmt(double(res_exp)); }

    // ---- [B12] Dot / Abs2: <v,v> = Σ|v_i|² (Dot 的第一个调用方) ----
    const Mat v = RandMatrix<T>(3, 1, 1006u, -1.0, 1.0);
    const real_t dot_dev = ScalarTraits<T>::abs2(Dot(v, v) - ScalarTraits<T>::fromReal(Abs2(v)));
    { SCOPED_TRACE(suit("B12", "Dot(v,v) − Abs2(v) = 0"));
      EXPECT_LE(dot_dev, TB2) << Fmt(std::sqrt(double(dot_dev))); }

    // ---- [B13] Matrix 容器入口 ----
    {
        const Mat O  = Mat::Ones(2, 2);
        const Mat I2 = Mat::Identity(2);
        bool okc = (O.Rows() == 2 && O.Cols() == 2 && O.Size() == 4 && O.Data() != nullptr);
        okc = okc && MaxAbs2(O.cwiseProduct(O) - O) == real_t(0);
        okc = okc && MaxAbs2(O.cwiseDivision(O) - O) == real_t(0);
        okc = okc && MaxAbs2(Mat::Diagonal(std::vector<T>{ R(real_t(2)), R(real_t(3)) }).Transpose()
                             - Mat::Diagonal(std::vector<T>{ R(real_t(2)), R(real_t(3)) })) == real_t(0);
        okc = okc && MaxAbs2(O + O - O * R(real_t(2))) == real_t(0);
        okc = okc && MaxAbs2(I2.Adjoint() - I2) == real_t(0);
        // 比较一律走 abs2 (D7): Complex/Dual 没有 operator==, 写 `x == one()` 编不过 (实测)
        okc = okc && ScalarTraits<T>::abs2(I2.at(1, 1) - ScalarTraits<T>::one()) == real_t(0);
        okc = okc && ScalarTraits<T>::abs2(O(0, 0) - ScalarTraits<T>::one()) == real_t(0);
        SCOPED_TRACE(suit("B13", "Zero/Ones/Identity/Diagonal/Transpose/Adjoint/cwise/+/-/at 出口"));
        EXPECT_TRUE(okc);
    }
}

}  // namespace

// ---- [A] 实例化闸门: 两域各数一次 (不是 TYPED_TEST —— 见文件头) ----
TEST(Smoke, InstantiationGate)
{
    const int nd = GateCount<D>();
    const int nc = GateCount<C>();
    { SCOPED_TRACE("[A1] double 实例化入口数 / 期望 " + std::to_string(kEntryCount));
      EXPECT_EQ(nd, kEntryCount); }
    { SCOPED_TRACE("[A2] Complex<double> 实例化入口数 / 期望 " + std::to_string(kEntryCount));
      EXPECT_EQ(nc, kEntryCount); }
}

// ---- [B] 最小数值闸门: 同一模板跑两域 ⇒ TYPED_TEST ----
template <typename T>
class SmokeNumericTest : public ::testing::Test
{
protected:
    static const char* Tag() { return std::is_same<T, double>::value ? "double" : "complex"; }
};

using SmokeNumericTypes = ::testing::Types<D, C>;
TYPED_TEST_SUITE(SmokeNumericTest, SmokeNumericTypes);

TYPED_TEST(SmokeNumericTest, NumericGate)
{
    NumericGate<TypeParam>(TestFixture::Tag());
}