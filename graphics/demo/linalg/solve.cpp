// solve.cpp -- linalg/solve.hpp 验收 (S1: LU + Solve)
//
// 为什么先做这一步:
//   V1.2 (求逆残差) / V1.3 (Jacobi 重构) / V1.4 (SVD 重构) / V1.5 (Dual 复矩阵 LU)
//   的数值全部由 LU + Solve 提供 —— 它错, 后面每一步的残差都不值得信。
//
// 判据清单 (判据先行: 先写「多少算对」再写实现):
//   ① 结构: L 单位下三角, U 上三角
//   ② P·A = L·U —— P 由 piv 独立重建, 不复用 LU 内部逻辑
//   ③ ‖A·X − B‖∞ / (1+‖B‖∞) < 1e-10   (double 与 Complex<double> 共用同一份源码)
//   ④ 与独立实现的"朴素列主元高斯消元"逐元素对拍 (两条代码路径)
//   ⑤ 手算用例: 3×3 强制在第 1 步换行 (两行都已带乘子) 的 A^{-1} 解析值
//   ⑥ 异常路径: 奇异 / 非方阵 / B 行数不符
//
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 判据名原样保留, 通过 SCOPED_TRACE 携带; 上面的 ①–⑥ 对应下面的 A–E 五个 TEST:
//       ① ② ③ ④ + 换行分支 -> [A] RandomBatchDouble
//       ⑤ (3×3 / 2×2)       -> [B] Hand3x3Pivoting / [C] Hand2x2
//       ① ② ③ (Complex)    -> [D] ComplexDomain
//       ⑥                   -> [E] Contract
//       平凡结构            -> [F] TrivialStructure
//   · 迁移前 34 项断言 -> 迁移后 34 条 EXPECT (逐条对应, 一条没合并)。
//   · 异常判据走 EXPECT_EXACT_THROW (老写法 catch (const std::runtime_error&)
//     会放过派生类; 判精确类型才是契约判据)。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                        | 边界 / 病态输入                  | 专抓什么错
//  ----|-----------------------------|----------------------------------|--------------------------
//   A  | 结构 1e-30; P·A 1e-12 相对  | n=2..5 × 8 seed (随机, 良态)     | 消元顺序 / piv 应用错位
//      | 残差 1e-10 相对; 对拍 1e-10 |                                  | 用同一份逻辑验自己
//   B  | 乘子 1e-14; 解 1e-12        | **第 0、1 步都换行**, 两行都已带  | **换行时乘子不跟着搬**
//      |                             | 乘子 (L(1,0)=0.5 / L(2,0)=0.25)  | (L·U ≠ P·A, 解仍"看着合理")
//   C  | 1e-14                       | 2×2 A·x=b, 解为整数 (1,2)        | 回代方向 / 多右端项错列
//   D  | 结构 1e-30; 1e-12 / 1e-10   | Complex<double>, 6 seed × n=2..4 | 用 `>` 选主元 (复域编不过)
//      |                             |                                  | 或比较忘了 abs2
//   E  | **异常类型精确到子类**      | 奇异 2×2 / 2×3 非方阵 / B 行数不符| 用 logic_error 判基类假绿
//   F  | 1e-15 / 1e-14               | A=I, A=diag(2,-4,0.5)            | 把 I / 对角当特例走错分支
//
// ---------------------------------------------------------------------------
// 已知盲区 (写下来防误判, 也是注入对照的方案)
// ---------------------------------------------------------------------------
//  ① **换行时不搬乘子** (删掉 SwapMultipliers 调用): 必红的是 [B] 的
//     "乘子随换行" 与 "P·A = L·U"。随机批量 [A] 在 n<=5 时也常红, 但**不保证**
//     (某些 seed 恰好不换行) —— 所以 [B] 的手算用例是这条 bug 的**唯一必红判据**。
//  ② **只换 [0, upto) 写成整行换**: 单位下三角的对角 1 被搬到上三角 ->
//     [A]/[B] 的"L/U 结构"必红 (结构判据卡 1e-30, 是构造出来的, 不是算出来的)。
//  ③ **不选主元** (max_row 恒 = j): 本文件的随机矩阵都良态, **可能全绿** ——
//     这是已知盲区, 本文件没有"不选主元必红"的判据。要抓它得喂病态输入
//     (如 [[1e-20,1],[1,1]]), 现在没有, 记在这里而不是假装测过。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh                       # 全量 (注册名: solve)
//   ./graphics/build/demo/linalg/solve            # 单跑
//   ./graphics/build/demo/linalg/solve --gtest_filter='*Contract*'   # 只跑契约段
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/solve.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using C = Complex<double>;

// 结构零/单位一的判据 (结构是构造出来的, 不是算出来的, 可以卡得很死)
template <CommutativeScalar T>
typename ScalarTraits<T>::real_t StructZero()
{
    using real_t = typename ScalarTraits<T>::real_t;
    return real_t(1e-30);
}

// ---- 独立对拍实现: 朴素列主元高斯消元 (与 LU 两条代码路径, 仅在 SwapRow 上共用) ----
template <CommutativeScalar T>
Matrix<T> NaiveGaussSolve(Matrix<T> A, Matrix<T> B)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int n = A.Rows();
    const int nCols = B.Cols();

    for (int j = 0; j < n; j++) {
        int p = j;
        real_t best = ScalarTraits<T>::abs2(A(j, j));
        for (int i = j + 1; i < n; i++) {
            const real_t v = ScalarTraits<T>::abs2(A(i, j));
            if (v > best) {
                best = v;
                p = i;
            }
        }
        if (p != j) {
            SwapRow(A, j, p);
            SwapRow(B, j, p);
        }
        for (int i = j + 1; i < n; i++) {
            const T mul = A(i, j) / A(j, j);
            for (int k = j; k < n; k++) A(i, k) -= mul * A(j, k);
            for (int k = 0; k < nCols; k++) B(i, k) -= mul * B(j, k);
        }
    }

    Matrix<T> X(n, nCols);
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < nCols; j++) {
            T acc = B(i, j);
            for (int k = i + 1; k < n; k++) acc -= A(i, k) * X(k, j);
            X(i, j) = acc / A(i, i);
        }
    }
    return X;
}

// 由 piv 独立重建置换矩阵: 按 j 递增依次交换 (j, piv[j]) —— 与 Solve 里的应用方式一致,
// 但代码完全独立, 否则就是在"用同一个错误验证自己"
template <CommutativeScalar T>
Matrix<T> PermFromPiv(const std::vector<int>& piv, int n)
{
    Matrix<T> P = Matrix<T>::Identity(n);
    for (int j = 0; j < n; j++) {
        if (piv[j] != j) SwapRow(P, j, piv[j]);
    }
    return P;
}

template <CommutativeScalar T>
bool LUStructureOk(const LUResult<T>& lu, int n)
{
    const auto z = StructZero<T>();
    for (int i = 0; i < n; i++) {
        if (ScalarTraits<T>::abs2(lu.L(i, i) - ScalarTraits<T>::one()) > z) return false;
        for (int j = 0; j < n; j++) {
            if (j > i && ScalarTraits<T>::abs2(lu.L(i, j)) > z) return false;
            if (j < i && ScalarTraits<T>::abs2(lu.U(i, j)) > z) return false;
        }
    }
    return true;
}

// 手算用例: A = [[2,0,1],[1,0,2],[4,1,0]], det = -3
//   LU 第 0 步选到第 2 行 (|4| 最大) -> 第 1 步 |U(1,1)|=0.25 < |U(2,1)|=0.5 又要换行,
//   而这两行此刻都已经各带一个乘子 (0.25 与 0.5) —— 正是"乘子跟着换行"的考点
Matrix<double> HandA()
{
    return Matrix<double>{{2.0, 0.0, 1.0}, {1.0, 0.0, 2.0}, {4.0, 1.0, 0.0}};
}

Matrix<double> HandAInv()
{
    return Matrix<double>{{2.0 / 3.0, -1.0 / 3.0, 0.0},
                          {-8.0 / 3.0, 4.0 / 3.0, 1.0},
                          {-1.0 / 3.0, 2.0 / 3.0, 0.0}};
}

}  // namespace

// 判据名进 SCOPED_TRACE; 标签与断言必须在同一个块内 (见 docs/pitfalls.md §十二)
#define STEP(label) SCOPED_TRACE(label)
#define CHECK_TRUE(label, cond)   \
    do {                          \
        STEP(label);              \
        EXPECT_TRUE(cond);        \
    } while (0)

// =====================================================================
// A. 随机矩阵批量 (double): 结构 / P·A = L·U / 残差 / 与朴素高斯对拍
// =====================================================================

TEST(LinalgSolve, RandomBatchDouble)
{
    double worst_pa_lu = 0.0, worst_res = 0.0, worst_naive = 0.0;
    int pivot_taken = 0;

    for (int n : {2, 3, 4, 5}) {
        bool struct_ok = true, pa_lu_ok = true, res_ok = true, naive_ok = true;
        for (unsigned seed = 1; seed <= 8; seed++) {
            const auto A = RandMatrix<double>(n, n, seed * 31u + n);
            const auto B = RandMatrix<double>(n, 2, seed * 71u + n);

            const auto lu = LU(A);
            if (!LUStructureOk(lu, n)) struct_ok = false;
            for (int j = 0; j < n; j++) {
                if (lu.piv[j] != j) pivot_taken++;
            }

            const auto P = PermFromPiv<double>(lu.piv, n);
            const double e_pa = NormInf(P * A - lu.L * lu.U);
            worst_pa_lu = std::max(worst_pa_lu, e_pa / (1.0 + NormInf(A)));
            if (e_pa > 1e-12 * (1.0 + NormInf(A))) pa_lu_ok = false;

            const auto X = Solve(A, B);
            const double e_res = NormInf(A * X - B);
            worst_res = std::max(worst_res, e_res / (1.0 + NormInf(B)));
            if (e_res > 1e-10 * (1.0 + NormInf(B))) res_ok = false;

            const auto Xn = NaiveGaussSolve<double>(A, B);
            const double e_nv = NormInf(X - Xn);
            worst_naive = std::max(worst_naive, e_nv / (1.0 + NormInf(X)));
            if (e_nv > 1e-10 * (1.0 + NormInf(X))) naive_ok = false;
        }
        const std::string tag = "[A] n=" + std::to_string(n) + " (8 seed): ";
        { STEP(tag + "L 单位下三角 + U 上三角"); EXPECT_TRUE(struct_ok); }
        { STEP(tag + "P·A = L·U (相对误差 < 1e-12)");
          EXPECT_TRUE(pa_lu_ok) << "最差相对误差 = " << worst_pa_lu; }
        { STEP(tag + "‖A·X − B‖∞ 相对 < 1e-10");
          EXPECT_TRUE(res_ok) << "最差相对残差 = " << worst_res; }
        { STEP(tag + "与朴素列主元高斯消元对拍 < 1e-10");
          EXPECT_TRUE(naive_ok) << "最差对拍差 = " << worst_naive; }
    }
    { STEP("[A] 随机批量确实走到过换行分支 (pivot_taken > 0)");
      EXPECT_TRUE(pivot_taken > 0) << "pivot_taken = " << pivot_taken; }

    std::printf("[solve] 最差相对误差: P·A−L·U = %.3e | A·X−B = %.3e | 对拍 = %.3e\n",
                worst_pa_lu, worst_res, worst_naive);
}

// =====================================================================
// B. 手算用例: 第 1 步强制换行 + 两行都已带乘子
// =====================================================================

TEST(LinalgSolve, Hand3x3Pivoting)
{
    const auto A = HandA();
    const auto lu = LU(A);

    { STEP("[B] 手算 3×3: L/U 结构"); EXPECT_TRUE(LUStructureOk(lu, 3)); }

    { STEP("[B] 手算 3×3: piv = [2,2,?] (第 0、1 步都要换行)");
      EXPECT_TRUE(lu.piv[0] == 2 && lu.piv[1] == 2)
          << "piv = [" << lu.piv[0] << ", " << lu.piv[1] << ", " << lu.piv[2] << "]"; }

    // 乘子是否跟着换行了: L(1,0) 应为 0.5, L(2,0) 应为 0.25 (换行前是反的)
    { STEP("[B] 手算 3×3: 乘子随换行 (L(1,0)=0.5, L(2,0)=0.25)");
      EXPECT_TRUE(std::fabs(lu.L(1, 0) - 0.5) < 1e-14 && std::fabs(lu.L(2, 0) - 0.25) < 1e-14)
          << "L(1,0) = " << lu.L(1, 0) << ", L(2,0) = " << lu.L(2, 0); }

    { STEP("[B] 手算 3×3: L(2,1)=0.5"); EXPECT_NEAR(lu.L(2, 1), 0.5, 1e-14); }

    { STEP("[B] 手算 3×3: P·A = L·U");
      const auto P = PermFromPiv<double>(lu.piv, 3);
      EXPECT_TRUE(IsApprox(P * A, lu.L * lu.U, 1e-14))
          << "‖P·A − L·U‖∞ = " << NormInf(P * A - lu.L * lu.U); }

    { STEP("[B] 手算 3×3: Solve(A,I) = A^{-1} (解析值)");
      const auto X = Solve(A, Matrix<double>::Identity(3));
      EXPECT_TRUE(IsApprox(X, HandAInv(), 1e-12))
          << "‖X − A^{-1}‖∞ = " << NormInf(X - HandAInv()); }

    { STEP("[B] 手算 3×3: A·A^{-1} = I");
      const auto X = Solve(A, Matrix<double>::Identity(3));
      EXPECT_TRUE(IsApprox(A * X, Matrix<double>::Identity(3), 1e-12))
          << "‖A·A^{-1} − I‖∞ = " << NormInf(A * X - Matrix<double>::Identity(3)); }
}

// 手算 2×2: x 应为整数 (1, 2), 便于人眼核对
TEST(LinalgSolve, Hand2x2)
{
    Matrix<double> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double> b{{5.0}, {11.0}};
    const auto x = Solve(A, b);

    { STEP("[C] 手算 2×2 A·x = b: x = (1, 2)");
      EXPECT_TRUE(std::fabs(x(0, 0) - 1.0) < 1e-14 && std::fabs(x(1, 0) - 2.0) < 1e-14)
          << "x = (" << x(0, 0) << ", " << x(1, 0) << ")"; }
}

// =====================================================================
// D. 复数: 与 double 走同一份源码 (V1.1 的一半)
// =====================================================================

TEST(LinalgSolve, ComplexDomain)
{
    double worst_c = 0.0;
    bool c_struct = true, c_pa = true, c_res = true;

    for (unsigned seed = 1; seed <= 6; seed++) {
        for (int n : {2, 3, 4}) {
            const auto A = RandMatrix<C>(n, n, seed * 13u + n);
            const auto B = RandMatrix<C>(n, 2, seed * 17u + n);
            const auto lu = LU(A);
            if (!LUStructureOk(lu, n)) c_struct = false;

            const auto P = PermFromPiv<C>(lu.piv, n);
            const double e_pa = NormInf(P * A - lu.L * lu.U);
            if (e_pa > 1e-12 * (1.0 + NormInf(A))) c_pa = false;

            const auto X = Solve(A, B);
            const double e_res = NormInf(A * X - B);
            worst_c = std::max(worst_c, e_res / (1.0 + NormInf(B)));
            if (e_res > 1e-10 * (1.0 + NormInf(B))) c_res = false;
        }
    }

    { STEP("[D] Complex<double>: L/U 结构 (abs2 选主元, 无 operator>)"); EXPECT_TRUE(c_struct); }
    { STEP("[D] Complex<double>: P·A = L·U"); EXPECT_TRUE(c_pa); }
    { STEP("[D] Complex<double>: ‖A·X − B‖∞ 相对 < 1e-10");
      EXPECT_TRUE(c_res) << "最差相对残差 = " << worst_c; }

    std::printf("[solve] 复数最差相对残差 = %.3e\n", worst_c);
}

// =====================================================================
// E. 边界与异常路径
// =====================================================================

TEST(LinalgSolve, Contract)
{
    const Matrix<double> sing{{1.0, 2.0}, {2.0, 4.0}};   // 第 2 行 = 2 × 第 1 行

    { STEP("[E] 奇异矩阵 -> throw std::runtime_error");
      EXPECT_EXACT_THROW((void)LU(sing), std::runtime_error); }

    { STEP("[E] 奇异矩阵 Solve -> throw");
      EXPECT_EXACT_THROW((void)Solve(sing, Matrix<double>::Identity(2)), std::runtime_error); }

    { STEP("[E] 非方阵 LU -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)LU(Matrix<double>(2, 3)), std::invalid_argument); }

    { STEP("[E] B 行数不符 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)Solve(Matrix<double>::Identity(3), Matrix<double>(2, 1)),
                         std::invalid_argument); }
}

// =====================================================================
// F. 平凡结构 (不给优化, 但要保证不退化)
// =====================================================================

TEST(LinalgSolve, TrivialStructure)
{
    const auto B = RandMatrix<double>(3, 2, 99u);

    { STEP("[F] Solve(I, B) = B");
      EXPECT_TRUE(IsApprox(Solve(Matrix<double>::Identity(3), B), B, 1e-15)); }

    const double d[3] = {2.0, -4.0, 0.5};
    Matrix<double> D = Matrix<double>::Diagonal({d[0], d[1], d[2]});
    const auto X = Solve(D, B);
    bool diag_ok = true;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            if (std::fabs(X(i, j) - B(i, j) / d[i]) > 1e-14) diag_ok = false;
        }
    }
    { STEP("[F] Solve(Diagonal, B) 逐列相除"); EXPECT_TRUE(diag_ok); }
}