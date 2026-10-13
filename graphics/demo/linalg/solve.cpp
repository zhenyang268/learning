// solve.cpp -- linalg/solve.hpp 验收 (S1: LU + Solve)
//
// 为什么先做这一步:
//   V1.2 (求逆残差) / V1.3 (Jacobi 重构) / V1.4 (SVD 重构) / V1.5 (Dual 复矩阵 LU)
//   的数值全部由 LU + Solve 提供 —— 它错, 后面每一步的残差都不值得信。
//
// 判据 (matrix.md §10.9 第 2 条: 判据先行):
//   ① 结构: L 单位下三角, U 上三角
//   ② P·A = L·U —— P 由 piv 独立重建, 不复用 LU 内部逻辑
//   ③ ‖A·X − B‖∞ / (1+‖B‖∞) < 1e-10   (double 与 Complex<double> 共用同一份源码)
//   ④ 与独立实现的"朴素列主元高斯消元"逐元素对拍 (两条代码路径)
//   ⑤ 手算用例: 3×3 强制在第 1 步换行 (两行都已带乘子) 的 A^{-1} 解析值
//   ⑥ 异常路径: 奇异 / 非方阵 / B 行数不符
//
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh            # 全量
//   ./graphics/build/demo/linalg/solve # 单跑

#include "../common.hpp"
#include "linalg/solve.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using C = Complex<double>;

static Checker ck;

// 结构零/单位一的判据 (结构是构造出来的, 不是算出来的, 可以卡得很死)
template <CommutativeScalar T>
static typename ScalarTraits<T>::real_t StructZero()
{
    using real_t = typename ScalarTraits<T>::real_t;
    return real_t(1e-30);
}

// ---- 独立对拍实现: 朴素列主元高斯消元 (与 LU 两条代码路径, 仅在 SwapRow 上共用) ----
template <CommutativeScalar T>
static Matrix<T> NaiveGaussSolve(Matrix<T> A, Matrix<T> B)
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
static Matrix<T> PermFromPiv(const std::vector<int>& piv, int n)
{
    Matrix<T> P = Matrix<T>::Identity(n);
    for (int j = 0; j < n; j++) {
        if (piv[j] != j) SwapRow(P, j, piv[j]);
    }
    return P;
}

template <CommutativeScalar T>
static bool LUStructureOk(const LUResult<T>& lu, int n)
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
static Matrix<double> HandA()
{
    return Matrix<double>{{2.0, 0.0, 1.0}, {1.0, 0.0, 2.0}, {4.0, 1.0, 0.0}};
}

static Matrix<double> HandAInv()
{
    return Matrix<double>{{2.0 / 3.0, -1.0 / 3.0, 0.0},
                          {-8.0 / 3.0, 4.0 / 3.0, 1.0},
                          {-1.0 / 3.0, 2.0 / 3.0, 0.0}};
}

int main()
{
    std::printf("[solve] LU / Solve 验收\n");

    // =====================================================================
    // A. 随机矩阵批量 (double): 结构 / P·A = L·U / 残差 / 与朴素高斯对拍
    // =====================================================================
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
        char name[160];
        std::snprintf(name, sizeof(name), "n=%d (8 seed): L 单位下三角 + U 上三角", n);
        ck.check(struct_ok, name);
        std::snprintf(name, sizeof(name), "n=%d (8 seed): P·A = L·U (相对误差 < 1e-12)", n);
        ck.check(pa_lu_ok, name);
        std::snprintf(name, sizeof(name), "n=%d (8 seed): ‖A·X − B‖∞ 相对 < 1e-10", n);
        ck.check(res_ok, name);
        std::snprintf(name, sizeof(name), "n=%d (8 seed): 与朴素列主元高斯消元对拍 < 1e-10", n);
        ck.check(naive_ok, name);
    }
    ck.check(pivot_taken > 0, "随机批量确实走到过换行分支 (pivot_taken > 0)");

    std::printf("        最差相对误差: P·A−L·U = %.3e | A·X−B = %.3e | 对拍 = %.3e\n",
                worst_pa_lu, worst_res, worst_naive);

    // =====================================================================
    // B. 手算用例: 第 1 步强制换行 + 两行都已带乘子
    // =====================================================================
    {
        const auto A = HandA();
        const auto lu = LU(A);
        ck.check(LUStructureOk(lu, 3), "手算 3×3: L/U 结构");
        ck.check(lu.piv[0] == 2 && lu.piv[1] == 2,
                 "手算 3×3: piv = [2,2,?] (第 0、1 步都要换行)");

        // 乘子是否跟着换行了: L(1,0) 应为 0.5, L(2,0) 应为 0.25 (换行前是反的)
        ck.check(std::fabs(lu.L(1, 0) - 0.5) < 1e-14 && std::fabs(lu.L(2, 0) - 0.25) < 1e-14,
                 "手算 3×3: 乘子随换行 (L(1,0)=0.5, L(2,0)=0.25)");
        ck.check(std::fabs(lu.L(2, 1) - 0.5) < 1e-14, "手算 3×3: L(2,1)=0.5");

        const auto P = PermFromPiv<double>(lu.piv, 3);
        ck.check(IsApprox(P * A, lu.L * lu.U, 1e-14), "手算 3×3: P·A = L·U");

        const auto X = Solve(A, Matrix<double>::Identity(3));
        ck.check(IsApprox(X, HandAInv(), 1e-12), "手算 3×3: Solve(A,I) = A^{-1} (解析值)");
        ck.check(IsApprox(A * X, Matrix<double>::Identity(3), 1e-12), "手算 3×3: A·A^{-1} = I");
    }

    // 手算 2×2: x 应为整数 (1, 2), 便于人眼核对
    {
        Matrix<double> A{{1.0, 2.0}, {3.0, 4.0}};
        Matrix<double> b{{5.0}, {11.0}};
        const auto x = Solve(A, b);
        ck.check(std::fabs(x(0, 0) - 1.0) < 1e-14 && std::fabs(x(1, 0) - 2.0) < 1e-14,
                 "手算 2×2 A·x = b: x = (1, 2)");
    }

    // =====================================================================
    // C. 复数: 与 double 走同一份源码 (V1.1 的一半)
    // =====================================================================
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
        ck.check(c_struct, "Complex<double>: L/U 结构 (abs2 选主元, 无 operator>)");
        ck.check(c_pa, "Complex<double>: P·A = L·U");
        ck.check(c_res, "Complex<double>: ‖A·X − B‖∞ 相对 < 1e-10");
        std::printf("        复数最差相对残差 = %.3e\n", worst_c);
    }

    // =====================================================================
    // D. 边界与异常路径
    // =====================================================================
    {
        Matrix<double> sing{{1.0, 2.0}, {2.0, 4.0}};   // 第 2 行 = 2 × 第 1 行
        bool threw = false;
        try {
            (void)LU(sing);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        ck.check(threw, "奇异矩阵 -> throw std::runtime_error");

        threw = false;
        try {
            (void)Solve(sing, Matrix<double>::Identity(2));
        } catch (const std::runtime_error&) {
            threw = true;
        }
        ck.check(threw, "奇异矩阵 Solve -> throw");

        threw = false;
        try {
            (void)LU(Matrix<double>(2, 3));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        ck.check(threw, "非方阵 LU -> throw std::invalid_argument");

        threw = false;
        try {
            (void)Solve(Matrix<double>::Identity(3), Matrix<double>(2, 1));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        ck.check(threw, "B 行数不符 -> throw std::invalid_argument");
    }

    // =====================================================================
    // E. 平凡结构 (不给优化, 但要保证不退化)
    // =====================================================================
    {
        const auto B = RandMatrix<double>(3, 2, 99u);
        ck.check(IsApprox(Solve(Matrix<double>::Identity(3), B), B, 1e-15),
                 "Solve(I, B) = B");

        const double d[3] = {2.0, -4.0, 0.5};
        Matrix<double> D = Matrix<double>::Diagonal({d[0], d[1], d[2]});
        const auto X = Solve(D, B);
        bool diag_ok = true;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                if (std::fabs(X(i, j) - B(i, j) / d[i]) > 1e-14) diag_ok = false;
            }
        }
        ck.check(diag_ok, "Solve(Diagonal, B) 逐列相除");
    }

    return CheckReport(ck.ok, ck.fail, "[solve]");
}
