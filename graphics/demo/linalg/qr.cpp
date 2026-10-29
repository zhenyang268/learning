// qr.cpp -- linalg/qr.hpp 验收 (S6b: householder / givens / mgs 三个分解统一口径 / V2.1 硬判据)
//
// 判据清单 (判据先行: 先写「多少算对」再写实现):
//   ① V2.1 硬指标: 重构 ‖Q·R − A‖∞ < 1e-8 (实 / 复)
//   ② Q 列正交 / 酉: ‖Q^H·Q − I_n‖∞ < 1e-12
//   ③ R 严格下三角恒为 0 (上三角口径)
//   ④ 尺寸口径统一: Q 是 m×n, R 是 n×n (thin / reduced QR)
//   ⑤ 三个分解 (Householder / Givens / MGS) 在同一个矩阵上各自满足 ①–④
//   ⑥ 解析对拍: 对角阵 / 单位阵 / 2x2 闭式 |R11|,|R22| (含复)
//   ⑦ 边界: 1x1, m>n (3x2 / 5x3 / 4x1), 秩亏 (列相关 -> 对角为 0 但 Q 仍列正交);
//      异常: 空矩阵 / m < n (thin 约定下统一 throw)
//
// 为什么看 ‖Q^H·Q − I‖ 而不是直接比 Q: 同一个 A, 三个实现给出的 Q 差一个对角相位
// (R 对角符号约定各不相同), 所以 Q 本身不可比; 但 Q·R = A 与 Q^H·Q = I 与约定无关。
//
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 迁移前 189 项断言 -> 迁移后 189 条 EXPECT/ADD_FAILURE (逐条对应, 一条没合并)。
//   · 段号 [A]–[E] 原样保留; 每个 (矩阵, 实现) 组合都进 SCOPED_TRACE。
//   · 老 helper 的 catch 改成 `ADD_FAILURE()` (一条失败, 同段其它实现继续跑)。
//   · 异常判据走 EXPECT_EXACT_THROW: 空矩阵与 m<n 都是 `invalid_argument`。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                     | 边界 / 病态输入                     | 专抓什么错
//  ----|--------------------------|-------------------------------------|--------------------------
//   A  | 重构 1e-8; 正交 1e-12    | 2×2 / 3×2 / 5×3 / 4×1 / 1×1 / 5×5    | 三种实现各自的分支漏 (1×1 最易漏)
//      | 下三角 1e-14             | **x1 = -1 且尾巴 1e-6** (病态)       | Householder 符号选择 (不带 phase
//      |                          |                                     | 会灾难性相消, 精度掉到 1e-6)
//   B  | 同上                     | 复 2×2 / 4×3 / 2×1 / **纯相位列**    | 复域漏 conj; 纯相位逼出相位分支
//      |                          | 复 2×2 病态 x1<0                    | (实旋转对纯相位列无能为力)
//   C  | 1e-12 逐元 (|diag R|)    | diag(3,-4,5) / I3 / 2×2 闭式        | 排序假设 (QR 不排序) / 对角符号
//      |                          | 复 diag(i,2) / 复 2×2 闭式           | 复域 |R_ii| 取模而不是取实部
//   D  | 重构 1e-8; 正交 1e-12    | **列线性相关** / **4 个零矩阵块**    | 秩亏时补正交方向的那段 (qr.hpp:138)
//   E  | **异常类型精确到子类**   | 0×0 / 2×3 (m<n)                     | 用 logic_error 判基类假绿
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **R 的严格下三角判据 1e-14 看起来像"精确 0", 其实不是**:
//     它是浮点残渣阈值 (矩阵乘的舍入), 别在别的文件里照抄成"判精确"。
//  ② **Householder 的 phase 选择** (v1 = x1 + phase(x1)·‖x‖):
//     抓它的只有 [A] 的 `实 2x2 病态 x1<0` 与 [B] 的 `复 2x2 病态 x1<0` ——
//     其余矩阵上正确性无损, 只是精度差, 1e-8 的阈值完全看不见。
//  ③ **mgs 的列正交退化** (不做二次正交化) 在良态矩阵上抓不到,
//     只有 [D] 的秩亏 2×2 与 [A] 的病态例可能显形。
//  ④ 三个实现**共用同一份量度函数** —— 若量度本身写错 (如 OrthResid 忘了减 I),
//     三条判据会一起假绿。跨实现不可比 Q (差一个对角相位), 所以这是本文件的结构性局限。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh
//   ./graphics/build/demo/linalg/qr
//   ./graphics/build/demo/linalg/qr --gtest_filter='*Contract*'
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/qr.hpp"

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

constexpr double V21_TOL = 1e-8;    // V2.1 硬阈值
constexpr double ORTH_TOL = 1e-12;  // 列正交阈值
// 下三角残渣阈值: 消元的第 j 列在浮点下必留 ~eps 级残渣 (矩阵乘法的舍入),
// 判"== 0"是判错了 —— 判据要求的是"量级够小", 不是"逐位为 0"。
constexpr double LOWER_TOL = 1e-14;

std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

template <CommutativeScalar T>
double Abs1(const T& x) { return std::sqrt(double(ScalarTraits<T>::abs2(x))); }

// ---- 量度 ------------------------------------------------------------------

// ‖Q^H·Q − I_n‖∞ —— 列正交 / 酉
template <CommutativeScalar T>
double OrthResid(const Matrix<T>& Q)
{
    return double(NormInf(Q.Adjoint() * Q - Matrix<T>::Identity(Q.Cols())));
}

// ‖Q·R − A‖∞ —— V2.1 判据本体
template <CommutativeScalar T>
double ReconResid(const Matrix<T>& A, const QRcomposition<T>& qr)
{
    return double(NormInf(qr.Q * qr.R - A));
}

// R 的严格下三角部分最大模 (应该恒为 0)
template <CommutativeScalar T>
double LowerMax(const Matrix<T>& R)
{
    double best = 0.0;
    for (int i = 1; i < R.Rows(); i++)
        for (int j = 0; j < i && j < R.Cols(); j++)
            best = std::max(best, Abs1(R(i, j)));
    return best;
}

// ---- 三个实现登记表 (与 qr.hpp 的公开入口一一对应) -------------------------

template <CommutativeScalar T>
struct Entry
{
    const char* name;
    QRcomposition<T> (*fn)(const Matrix<T>&);
};

template <CommutativeScalar T>
const Entry<T> kEntries[3] = {
    {"householder", &QR_Householder<T>},
    {"givens",      &QR_givens<T>},
    {"mgs",         &QR_Shimit<T>},
};

// ---- 主判据: 三个实现跑同一套 ----------------------------------------------

template <CommutativeScalar T>
void RunCase(const char* sec, const Matrix<T>& A, const std::string& tag)
{
    const int m = A.Rows();
    const int n = A.Cols();

    for (const Entry<T>& e : kEntries<T>) {
        const std::string who = std::string("[") + sec + "] " + tag + " [" + e.name + "]";
        try {
            const QRcomposition<T> qr = e.fn(A);

            const bool shape = (qr.Q.Rows() == m && qr.Q.Cols() == n &&
                                qr.R.Rows() == n && qr.R.Cols() == n);
            const double rr = ReconResid(A, qr);
            const double ro = OrthResid(qr.Q);
            const double lo = LowerMax(qr.R);

            { SCOPED_TRACE(who + " 尺寸 Q=m×n, R=n×n");
              EXPECT_TRUE(shape) << "Q " << qr.Q.Rows() << "x" << qr.Q.Cols()
                                 << ", R " << qr.R.Rows() << "x" << qr.R.Cols(); }
            { SCOPED_TRACE(who + " 重构 ‖QR−A‖∞ < 1e-8");
              EXPECT_LT(rr, V21_TOL) << Fmt(rr); }
            { SCOPED_TRACE(who + " Q 列正交 ‖QᴴQ−I‖∞ < 1e-12");
              EXPECT_LT(ro, ORTH_TOL) << Fmt(ro); }
            { SCOPED_TRACE(who + " R 严格下三角 < 1e-14");
              EXPECT_LT(lo, LOWER_TOL) << Fmt(lo); }
        } catch (const std::exception& ex) {
            ADD_FAILURE() << who << " 抛异常: " << ex.what();
        }
    }
}

// ---- 解析对拍: 检查 |R| 的对角 (与实现无关的解析值) -------------------------

template <CommutativeScalar T>
void CheckDiagR(const Matrix<T>& A, const std::vector<double>& want, const std::string& tag)
{
    for (const Entry<T>& e : kEntries<T>) {
        try {
            const QRcomposition<T> qr = e.fn(A);
            bool all = true;
            for (int i = 0; i < (int)want.size(); i++)
                if (std::abs(Abs1(qr.R(i, i)) - want[i]) > 1e-12) all = false;
            std::string got;
            for (int i = 0; i < qr.R.Rows(); i++) {
                if (i) got += ", ";
                got += Fmt(Abs1(qr.R(i, i)));
            }
            SCOPED_TRACE("[C] " + tag + " [" + e.name + "] |diag R| = " + got);
            EXPECT_TRUE(all);
        } catch (const std::exception& ex) {
            ADD_FAILURE() << "[C] " << tag << " [" << e.name << "] 抛异常: " << ex.what();
        }
    }
}

}  // namespace

#define STEP(label) SCOPED_TRACE(label)

// ---- [A] 实数: 方阵 / 长方阵 / 退化形状 ----
TEST(Qr, RealCases)
{
    { Matrix<D> A2 = {{1.0, 2.0}, {3.0, 4.0}}; RunCase("A", A2, "实 2x2"); }

    { Matrix<D> A32 = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}; RunCase("A", A32, "实 3x2"); }

    { Matrix<D> A53 = RandMatrix<D>(5, 3, 11u); RunCase("A", A53, "实 5x3 随机"); }

    { Matrix<D> A41 = RandMatrix<D>(4, 1, 13u); RunCase("A", A41, "实 4x1"); }

    { Matrix<D> A11 = {{7.0}}; RunCase("A", A11, "实 1x1"); }

    { Matrix<D> A55 = RandMatrix<D>(5, 5, 17u); RunCase("A", A55, "实 5x5 随机"); }

    // 病态: 第 0 列 x1 = -1 为负、且 |x1| ≈ ||x|| (尾部只有 1e-6)。
    // Householder 的 v1 = x1 + phase(x1)*||x|| —— 带 phase 是"加"(不消), 不带就是
    // -1 + 1 = 0 的灾难性相消, 会丢掉 1e-6 那一项、并把误差放大到相对 1e-6。
    // 这是"符号选择"在实数域唯一能显形的场景 (正确性无损, 精度有损)。
    { Matrix<D> AIll = {{-1.0, 0.0}, {1e-6, 1.0}}; RunCase("A", AIll, "实 2x2 病态 x1<0"); }
}

// ---- [B] 复数: 方阵 / 长方阵 / 列向量 ----
TEST(Qr, ComplexCases)
{
    { Matrix<C> B2 = {{C(1.0, 1.0), C(2.0, -1.0)}, {C(3.0, 0.0), C(0.0, 2.0)}};
      RunCase("B", B2, "复 2x2"); }

    { Matrix<C> B43 = RandMatrix<C>(4, 3, 21u); RunCase("B", B43, "复 4x3 随机"); }

    { Matrix<C> B21 = Matrix<C>(2, 1);
      B21(0, 0) = C(1.0, 2.0);
      B21(1, 0) = C(3.0, -1.0);
      RunCase("B", B21, "复 2x1"); }

    // 纯相位列 (实部全 0) —— 实旋转在这里完全无能为力, 逼出相位分支
    { Matrix<C> Bph = Matrix<C>(2, 2);
      Bph(0, 0) = C(0.0, 1.0);  Bph(0, 1) = C(0.0, 0.0);
      Bph(1, 0) = C(0.0, -2.0); Bph(1, 1) = C(0.0, 3.0);
      RunCase("B", Bph, "复 2x2 纯相位"); }

    // 病态: x1 = -1 (负实数), 尾部 1e-6 —— 同实数域那条, 验证 phase 也吃这一路
    { Matrix<C> CIll = {{C(-1.0, 0.0), C(0.0, 1.0)}, {C(1e-6, 0.0), C(1.0, 0.0)}};
      RunCase("B", CIll, "复 2x2 病态 x1<0"); }
}

// ---- [C] 解析对拍: |diag R| 该是多少就是多少 ----
TEST(Qr, AnalyticDiagR)
{
    // 对角阵: 各列正交, |R_ii| = |a_ii|; 顺序不变 (QR 不排序)
    { Matrix<D> Dg = Matrix<D>::Diagonal({3.0, -4.0, 5.0});
      CheckDiagR(Dg, {3.0, 4.0, 5.0}, "实 diag(3,-4,5)"); }

    { Matrix<D> Id = Matrix<D>::Identity(3);
      CheckDiagR(Id, {1.0, 1.0, 1.0}, "实 I3"); }

    // 2x2 闭式: A = [[1,2],[3,4]] -> |R11| = sqrt(10), |R22| = |det|/|R11| = 2/sqrt(10)
    { Matrix<D> A2 = {{1.0, 2.0}, {3.0, 4.0}};
      CheckDiagR(A2, {std::sqrt(10.0), 2.0 / std::sqrt(10.0)}, "实 2x2 闭式"); }

    // 复对角: |R_ii| = |a_ii|
    { Matrix<C> Cd = Matrix<C>::Diagonal({C(0.0, 1.0), C(2.0, 0.0)});
      CheckDiagR(Cd, {1.0, 2.0}, "复 diag(i,2)"); }

    // 复 2x2: |R11| = ||col0|| = sqrt(2);  det = 2 - (1+i) = 1-i, |det| = sqrt(2),
    // 故 |R22| = |det|/|R11| = sqrt(2)/sqrt(2) = 1
    { Matrix<C> C2 = {{C(1.0, 0.0), C(1.0, 1.0)}, {C(1.0, 0.0), C(2.0, 0.0)}};
      CheckDiagR(C2, {std::sqrt(2.0), 1.0}, "复 2x2 闭式"); }
}

// ---- [D] 秩亏: 列线性相关 ----
TEST(Qr, RankDeficient)
{
    { Matrix<D> Rk = {{1.0, 2.0}, {2.0, 4.0}};    // col1 = 2 * col0
      RunCase("D", Rk, "实 2x2 秩亏"); }

    { Matrix<D> Rk0 = {{0.0, 0.0}, {0.0, 0.0}};   // 零矩阵: 全部列都秩亏
      RunCase("D", Rk0, "实 2x2 零矩阵"); }
}

// ---- [E] 异常口径 ----
TEST(Qr, Contract)
{
    { STEP("[E] householder 空矩阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_Householder(Matrix<D>()), std::invalid_argument); }
    { STEP("[E] givens 空矩阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_givens(Matrix<D>()), std::invalid_argument); }
    { STEP("[E] mgs 空矩阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_Shimit(Matrix<D>()), std::invalid_argument); }

    Matrix<D> Wide = RandMatrix<D>(2, 3, 31u);                 // m < n: thin 约定下退化
    { STEP("[E] householder m<n -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_Householder(Wide), std::invalid_argument); }
    { STEP("[E] givens m<n -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_givens(Wide), std::invalid_argument); }
    { STEP("[E] mgs m<n -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)QR_Shimit(Wide), std::invalid_argument); }
}