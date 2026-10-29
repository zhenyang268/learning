// inverse.cpp -- linalg/inverse.hpp 验收 (S2: Inverse 分发; S2a-2: SPD 档接入; S-App6: qr_inverse)
//
// 判据清单 (判据先行: 先写「多少算对」再写实现):
//   ① V1.2 硬指标: 随机矩阵 ‖A·A⁻¹ − I‖∞ < 1e-8, 且 ‖A⁻¹·A − I‖∞ < 1e-8
//   ② double 与 Complex<double> 共用同一份 linalg 源码 (V1.1)
//   ③ 结构分支: 对角 -> 逐元倒数; 正交(旋转) -> A^T; 单位阵 -> I (解析对拍)
//   ④ Classify 语义: 对角/旋转/随机/对称 分别落到对应 kind
//   ⑤ Inverse(A, hint) 与 Inverse(A) 一致 (hint 只跳过探测, 不改结果)
//   ⑥ 异常: 奇异 / 非方阵 / 0x0
//   ⑦ SPD 档 (S2a-2): Classify 分出 SymmetricPositiveDefinite 走 cholesky_inverse,
//      与 LU 路径对拍 + 正定性判据 (结构探测口径)
//   ⑧ qr_inverse (S-App6, 2026-10-16): A = Q·R ⇒ A⁻¹ = R⁻¹·Qᴴ (对上三角 R 回代);
//      残差 + 与 lu_inverse 逐元一致 + 奇异/非方阵/空矩阵 throw; double 与 Complex 各一遍
//
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 迁移前 85 项断言 -> 迁移后 85 条 EXPECT (逐条对应, 一条没合并)。
//   · 判据名原样保留并加段号 [A]–[H], 通过 SCOPED_TRACE 携带; 段 ↔ 上面的清单:
//       A/B -> ① ②, C -> ③, D -> ④, E -> ⑤, F -> ⑥, G -> ⑦, H -> ⑧。
//   · 异常判据走 EXPECT_EXACT_THROW (老写法 `catch (const std::exception&)` 连
//     `logic_error` 都收, 判不出"抛哪一种"; 契约见 inverse.hpp「一、API 契约」)。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                       | 边界 / 病态输入                     | 专抓什么错
//  ----|----------------------------|-------------------------------------|--------------------------
//   A  | 1e-8 (V1.2 硬指标)         | n=2..6 × 4 seed, 随机良态           | LU 路径的残差退化
//   B  | 1e-8                       | Complex<double>, n=2..4 × 3 seed    | 实域能过、复域编不过/算错
//   C  | 对角 1e-15 / 正交 1e-12    | diag(2,-4,0.5), 2×2 旋转, I3, Hermitian| 快路径退化成通用 LU
//   D  | DefaultTol(3)              | **A^T·A 对角全 1 但非对角 0.5**      | is_orthogonal 只看对角线
//      |                            | 对称不定 / 上三角非对称              | SPD 档吞掉 Symmetric 档
//   E  | 1e-14                      | n=3..5 随机                          | hint 改了结果 (只该跳过探测)
//   F  | **异常类型精确到子类**     | 奇异 / 1×1? 非方阵 / 0×0 / 假 SPD hint| 用 logic_error 判基类假绿
//   G  | 1e-8 残差; 1e-10 两路差    | 真/假两侧都测 (半正定奇异、全负对角)  | **平铺 if 让 SPD 档永不可达**
//   H  | 1e-8 残差; 1e-10 两路差    | **复数 2×2** (Qᴴ 写成 Qᵀ 只有复域能抓)| 漏 conj / 与 LU 不是同一个逆
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **Classify 把 is_symmetric 排到 is_diagonal 前面**: 对角阵被对称档吃掉 ⇒
//     [D] 的「对角 -> Diagonal」与 [G] 的「单位阵 -> Diagonal」必红。
//     ⚠ 残差类判据**全绿** (对称路径也能算对) —— 这正是「必须专门证明分支被命中」的理由。
//  ② **qr_inverse 里 Qᴴ 写成 Qᵀ**: 实域**完全抓不到** (实矩阵 Qᵀ = Qᴴ),
//     只有 [H] ② 的复数 2×2 会红。别以为 double 那段能兜住。
//  ③ **is_orthogonal 只查对角线**: [D] 的 `A^T·A 对角全 1 但非对角非零` 是唯一必红判据
//     (随机矩阵的非对角本来就不为 0, 抓不到)。
//  ④ **cholesky 路径与 LU 路径差 1e-10 的判据**抓不到"两条路同错",
//     两条路都错的场合只有 [C] 的解析对拍 (diag/旋转/单位阵) 能兜。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh                        # 全量 (注册名: inverse)
//   ./graphics/build/demo/linalg/inverse           # 单跑
//   ./graphics/build/demo/linalg/inverse --gtest_filter='*QrInverse*'
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/inverse.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace
{

using D = double;
using C = Complex<double>;

constexpr double V12_TOL = 1e-8;   // V1.2 硬阈值

std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

// ‖A·X − I‖∞ (元素最大模; demo/common.hpp 的 NormInf 就是这个口径)
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t InvResid(const Matrix<T>& A, const Matrix<T>& X)
{
    return NormInf(A * X - Matrix<T>::Identity(A.Rows()));
}

template <CommutativeScalar T>
void CheckRandomInv(const char* sec, const char* tag, int n, unsigned seed, double& worst)
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T> A = RandMatrix<T>(n, n, seed);
    Matrix<T> X = Inverse(A);
    const real_t r1 = InvResid(A, X);
    const real_t r2 = InvResid(X, A);
    const real_t r  = (r1 > r2) ? r1 : r2;
    if (double(r) > worst) worst = double(r);

    SCOPED_TRACE(std::string("[") + sec + "] " + tag + " n=" + std::to_string(n)
                 + " seed=" + std::to_string(seed) + ": ‖A·A⁻¹ − I‖∞ < 1e-8");
    EXPECT_LT(double(r), V12_TOL) << Fmt(double(r));
}

}  // namespace

// 判据名进 SCOPED_TRACE; 标签与断言必须在同一个块内 (见 docs/pitfalls.md §十二)
#define STEP(label) SCOPED_TRACE(label)

// ---- [A] double 随机矩阵 (V1.2 硬指标) ----
TEST(Inverse, RandomDouble)
{
    double worst_d = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 4; sd++)
            CheckRandomInv<D>("A", "double", n, sd, worst_d);

    std::printf("[inverse] double 最差残差 %s\n", Fmt(worst_d).c_str());
}

// ---- [B] Complex 随机矩阵 (同一份源码) ----
TEST(Inverse, RandomComplex)
{
    double worst_c = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckRandomInv<C>("B", "complex", n, sd, worst_c);

    std::printf("[inverse] complex 最差残差 %s\n", Fmt(worst_c).c_str());
}

// ---- [C] 结构分支: 解析对拍 ----
TEST(Inverse, StructuralBranches)
{
    { STEP("[C] 对角 [[2,-4,0.5]]: 逆 = 逐元倒数");
      Matrix<D> Dg = Matrix<D>::Diagonal({2.0, -4.0, 0.5});
      Matrix<D> Di = Inverse(Dg);
      EXPECT_TRUE(std::abs(Di(0, 0) - 0.5) < 1e-15
                  && std::abs(Di(1, 1) + 0.25) < 1e-15
                  && std::abs(Di(2, 2) - 2.0) < 1e-15)
          << "diag = (" << Di(0, 0) << ", " << Di(1, 1) << ", " << Di(2, 2) << ")"; }

    const double th = 0.7;
    Matrix<D> R = {{std::cos(th), -std::sin(th)}, {std::sin(th), std::cos(th)}};
    Matrix<D> Ri = Inverse(R);

    { STEP("[C] 旋转 2x2: 逆 = A^T");
      EXPECT_TRUE(IsApprox(Ri, R.Transpose(), 1e-12))
          << "‖A⁻¹ − Aᵀ‖∞ = " << NormInf(Ri - R.Transpose()); }

    { STEP("[C] 旋转 2x2: 残差 < 1e-12");
      EXPECT_LT(NormInf(R * Ri - Matrix<D>::Identity(2)), 1e-12); }

    { STEP("[C] 单位阵: 逆 = 自身");
      Matrix<D> I3 = Matrix<D>::Identity(3);
      EXPECT_TRUE(IsApprox(Inverse(I3), I3, 1e-15)); }

    // Hermitian 正定 (复): S2a-2 起探测落 SymmetricPositiveDefinite, 走 cholesky_inverse
    { STEP("[C] Hermitian 2x2 (复): 残差 < 1e-12");
      Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};
      Matrix<C> Hi = Inverse(H);
      EXPECT_LT(NormInf(H * Hi - Matrix<C>::Identity(2)), 1e-12); }
}

// ---- [D] Classify 语义 ----
TEST(Inverse, ClassifySemantics)
{
    const double tol = DefaultTol<D>(3);
    Matrix<D> Dg = Matrix<D>::Diagonal({1.0, 2.0, 3.0});
    Matrix<D> Ro = {{std::cos(0.3), -std::sin(0.3)}, {std::sin(0.3), std::cos(0.3)}};
    Matrix<D> Ge = RandMatrix<D>(3, 3, 7u);
    Matrix<D> Sy = {{2.0, 1.0}, {1.0, 3.0}};    // 特征值 3.618 / 1.382 -> 正定
    Matrix<D> Id = {{1.0, 2.0}, {2.0, 1.0}};    // 特征值 3 / -1 -> 不定
    Matrix<D> NS = {{1.0, 2.0}, {0.0, 1.0}};    // 上三角非对称

    { STEP("[D] 对角 -> Diagonal");
      EXPECT_EQ(Classify(Dg, tol), MatrixKind::Diagonal); }
    { STEP("[D] 旋转 -> Orthogonal");
      EXPECT_EQ(Classify(Ro, tol), MatrixKind::Orthogonal); }
    { STEP("[D] 随机 -> Generic");
      EXPECT_EQ(Classify(Ge, tol), MatrixKind::Generic); }
    { STEP("[D] 实对称正定 -> SymmetricPositiveDefinite (S2a-2 打开该档)");
      EXPECT_EQ(Classify(Sy, tol), MatrixKind::SymmetricPositiveDefinite); }
    { STEP("[D] 实对称非正定 -> Symmetric (没被 SPD 档吞掉)");
      EXPECT_EQ(Classify(Id, tol), MatrixKind::Symmetric); }
    { STEP("[D] 上三角非对称 -> 不判 Symmetric");
      EXPECT_NE(Classify(NS, tol), MatrixKind::Symmetric); }

    // A = [[1,c],[0,s]] 满足 A^T·A = [[1,c],[c,1]]: 对角线全 1, 非对角非零。
    // 只查对角线的 is_orthogonal 会把它误判成正交, 全元素检查则不会。
    { STEP("[D] A^T·A 对角全 1 但非对角非零 -> 不判 Orthogonal");
      Matrix<D> FO = {{1.0, 0.5}, {0.0, std::sqrt(1.0 - 0.25)}};
      EXPECT_NE(Classify(FO, tol), MatrixKind::Orthogonal); }
}

// ---- [E] hint 一致性 ----
TEST(Inverse, HintConsistency)
{
    for (int n = 3; n <= 5; n++) {
        Matrix<D> A  = RandMatrix<D>(n, n, 11u);
        Matrix<D> X1 = Inverse(A);
        Matrix<D> X2 = Inverse(A, MatrixKind::Generic);
        SCOPED_TRACE("[E] n=" + std::to_string(n) + ": hint=Generic 与自动探测一致");
        EXPECT_TRUE(IsApprox(X1, X2, 1e-14)) << "‖X1 − X2‖∞ = " << NormInf(X1 - X2);
    }
}

// ---- [F] 异常路径 ----
TEST(Inverse, Contract)
{
    Matrix<D> singular = {{1.0, 2.0}, {2.0, 4.0}};
    { STEP("[F] 奇异矩阵 -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)Inverse(singular), std::runtime_error); }

    Matrix<D> rect(2, 3);
    { STEP("[F] 非方阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)Inverse(rect), std::invalid_argument); }

    Matrix<D> empty(0, 0);
    { STEP("[F] 0x0 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)Inverse(empty), std::invalid_argument); }

    // SPD 分支已于 S2a-2 接入, 不再 throw —— 正确性验收在 G 段。
    // 这里只留"hint 声明 SPD 但实际不正定"这条失败路径: 入参错就显式抛, 不静默算错。
    Matrix<D> NotSpd = {{1.0, 2.0}, {2.0, 1.0}};
    { STEP("[F] hint=SymmetricPositiveDefinite 但 A 不正定 -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)Inverse(NotSpd, MatrixKind::SymmetricPositiveDefinite),
                         std::runtime_error); }
    // qr_inverse 已于 S-App6 实现 (2026-10-16) —— 它的正确性与异常验收全部挪到 [H] 段。
}

// ---- [G] SPD 档 (S2a-2: cholesky 接入; 结构探测口径 "SPD 检测 = Cholesky 不选主元能否走通") ----
TEST(Inverse, SpdPath)
{
    const double tol2 = DefaultTol<D>(2);

    // ① 正定性判据: 真 / 假 两侧都测 (只测"能过"的一侧, 等于没测)
    Matrix<D> Spd  = {{2.0, 1.0}, {1.0, 3.0}};          // 特征值 3.618 / 1.382
    Matrix<D> Ind  = {{1.0, 2.0}, {2.0, 1.0}};          // 特征值 3 / -1
    Matrix<D> Neg  = Matrix<D>::Diagonal({-1.0, -2.0}); // 全负
    Matrix<D> Sem  = {{1.0, 2.0}, {2.0, 4.0}};          // 特征值 5 / 0 (半正定, 奇异)
    Matrix<D> Rect(2, 3);

    { STEP("[G] [[2,1],[1,3]] (特征值 3.6/1.4) -> true");
      EXPECT_TRUE(is_positive_definite(Spd)); }
    { STEP("[G] [[1,2],[2,1]] (特征值 3/-1) -> false");
      EXPECT_FALSE(is_positive_definite(Ind)); }
    { STEP("[G] 全负对角 -> false");
      EXPECT_FALSE(is_positive_definite(Neg)); }
    { STEP("[G] 半正定奇异 [[1,2],[2,4]] -> false (主元塌到 0)");
      EXPECT_FALSE(is_positive_definite(Sem)); }
    { STEP("[G] 单位阵 -> true");
      EXPECT_TRUE(is_positive_definite(Matrix<D>::Identity(4))); }
    { STEP("[G] 非方阵 -> false (不问 cholesky, 直接答否)");
      EXPECT_FALSE(is_positive_definite(Rect)); }

    // ② Classify 真把 SPD 与 Symmetric 分开 —— 平铺 if 会让 SPD 档永远不可达, 这条是回归闸门
    { STEP("[G] Classify(正定对称) -> SymmetricPositiveDefinite");
      EXPECT_EQ(Classify(Spd, tol2), MatrixKind::SymmetricPositiveDefinite); }
    { STEP("[G] Classify(非正定对称) -> Symmetric (SPD 档没把它吞掉)");
      EXPECT_EQ(Classify(Ind, tol2), MatrixKind::Symmetric); }
    { STEP("[G] Classify(单位阵) -> Diagonal (对角档排在 SPD 之前, 有意为之)");
      EXPECT_EQ(Classify(Matrix<D>::Identity(3), DefaultTol<D>(3)), MatrixKind::Diagonal); }

    // ③ 随机 SPD 对拍: cholesky 路径 vs LU 路径 (hint 只换算法, 不换结果)
    //    SPD 生成: S = M·Mᵀ + n·I —— 前项半正定, 加 n·I 保证严格正定且条件数可控。
    //    M·Mᵀ 逐位对称 (a*b == b*a, 求和顺序也相同), 所以不会跟 is_symmetric 的容差打架。
    double worst_spd = 0.0;
    for (int n = 2; n <= 6; n++) {
        Matrix<D> M = RandMatrix<D>(n, n, 200u + unsigned(n));
        Matrix<D> S = M * M.Transpose();
        for (int i = 0; i < n; i++) S(i, i) = S(i, i) + double(n);

        Matrix<D> Xc = Inverse(S, MatrixKind::SymmetricPositiveDefinite);   // cholesky
        Matrix<D> Xl = Inverse(S, MatrixKind::Generic);                     // LU
        const double rc = double(InvResid(S, Xc));
        const double rl = double(InvResid(S, Xl));
        const double dx = double(NormInf(Xc - Xl));
        const double r  = (rc > rl) ? rc : rl;
        if (r > worst_spd) worst_spd = r;

        const std::string tag = "[G] 随机 SPD n=" + std::to_string(n);
        { STEP(tag + ": Classify -> SymmetricPositiveDefinite");
          EXPECT_EQ(Classify(S, DefaultTol<D>(n)), MatrixKind::SymmetricPositiveDefinite); }
        { STEP(tag + ": ‖A·A⁻¹ − I‖∞ < 1e-8");
          EXPECT_LT(r, V12_TOL) << Fmt(r); }
        { STEP(tag + ": cholesky 与 LU 逐元差 < 1e-10");
          EXPECT_LT(dx, 1e-10) << Fmt(dx); }
    }
    std::printf("[inverse] (随机 SPD 最差残差 %s)\n", Fmt(worst_spd).c_str());

    // ④ 复 Hermitian 正定走同一份源码 (别又只测 double)
    Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};   // 特征值 4 / 1
    { STEP("[G] 复 Hermitian [[2,1-i],[1+i,3]] -> true");
      EXPECT_TRUE(is_positive_definite(H)); }
    { STEP("[G] 复 Hermitian 正定 -> SymmetricPositiveDefinite");
      EXPECT_EQ(Classify(H, DefaultTol<C>(2)), MatrixKind::SymmetricPositiveDefinite); }

    Matrix<C> Hc = Inverse(H, MatrixKind::SymmetricPositiveDefinite);
    Matrix<C> Hl = Inverse(H, MatrixKind::Generic);
    { STEP("[G] 复 SPD: cholesky 残差 < 1e-12");
      EXPECT_LT(NormInf(H * Hc - Matrix<C>::Identity(2)), 1e-12); }
    { STEP("[G] 复 SPD: cholesky 与 LU 逐元差 < 1e-12");
      EXPECT_LT(NormInf(Hc - Hl), 1e-12); }

    // ⑤ 自动探测入口 (走 Classify) 与 hint=Generic 一致 —— 与 E 段同口径
    { STEP("[G] Inverse(SPD) 自动探测与 hint=Generic 一致");
      Matrix<D> Aut = {{2.0, 1.0}, {1.0, 3.0}};
      EXPECT_TRUE(IsApprox(Inverse(Aut), Inverse(Aut, MatrixKind::Generic), 1e-14)); }
}

// ---- [H] qr_inverse (S-App6): A = Q·R ⇒ A⁻¹ = R⁻¹·Qᴴ ----
TEST(Inverse, QrInverse)
{
    double worst_qr = 0.0;

    // ① 随机方阵: 残差 + 与 LU 路径逐元一致 —— 两条独立实现互为判据
    //    (只判"残差小"是不够的: A·X ≈ I 在 X 有系统性偏差时也可能近似成立,
    //     与 lu_inverse 对拍才把"算的是同一个逆"钉住。)
    for (int n = 2; n <= 5; n++) {
        Matrix<D> A  = RandMatrix<D>(n, n, 500u + unsigned(n));
        Matrix<D> Xq = qr_inverse(A);
        Matrix<D> Xl = lu_inverse(A);
        const double rq = double(InvResid(A, Xq));
        const double dx = double(NormInf(Xq - Xl));
        if (rq > worst_qr) worst_qr = rq;

        SCOPED_TRACE("[H] double n=" + std::to_string(n)
                     + ": qr_inverse 残差 < 1e-8, 与 lu_inverse 逐元差 < 1e-10");
        EXPECT_TRUE(rq < V12_TOL && dx < 1e-10)
            << "残差 " << Fmt(rq) << ", 逐元差 " << Fmt(dx);
    }

    // ② 复数域必须单测 (D7 老坑: 漏 conj 时实数域永远 0 FAIL —— QR 里 Qᴴ 极易写成 Qᵀ)
    {
        Matrix<C> A  = {{C(2, 1), C(1, -3)}, {C(0, 2), C(1, 1)}};
        Matrix<C> Xq = qr_inverse(A);
        Matrix<C> Xl = lu_inverse(A);
        const double rq = double(InvResid(A, Xq));
        const double dx = double(NormInf(Xq - Xl));
        if (rq > worst_qr) worst_qr = rq;

        SCOPED_TRACE("[H] complex 2x2: qr_inverse 残差 < 1e-8, 与 lu_inverse 逐元差 < 1e-10");
        EXPECT_TRUE(rq < V12_TOL && dx < 1e-10)
            << "残差 " << Fmt(rq) << ", 逐元差 " << Fmt(dx);
    }

    // ③ 异常路径: 奇异 (R 对角出 ~0) / 非方阵 / 空矩阵
    Matrix<D> singular = {{1.0, 2.0}, {2.0, 4.0}};
    { STEP("[H] qr_inverse 奇异矩阵 -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)qr_inverse(singular), std::runtime_error); }

    Matrix<D> rect(2, 3);
    { STEP("[H] qr_inverse 非方阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)qr_inverse(rect), std::invalid_argument); }

    Matrix<D> empty(0, 0);
    { STEP("[H] qr_inverse 0x0 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)qr_inverse(empty), std::invalid_argument); }

    std::printf("[inverse] qr_inverse 最差残差 %s\n", Fmt(worst_qr).c_str());
}