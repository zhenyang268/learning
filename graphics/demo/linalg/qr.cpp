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
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh
//   ./graphics/build/demo/linalg/qr

#include "../common.hpp"
#include "linalg/qr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

static Checker ck;

using D = double;
using C = Complex<double>;

constexpr double V21_TOL = 1e-8;    // V2.1 硬阈值
constexpr double ORTH_TOL = 1e-12;  // 列正交阈值
// 下三角残渣阈值: 消元的第 j 列在浮点下必留 ~eps 级残渣 (矩阵乘法的舍入),
// 判"== 0"是判错了 —— 判据要求的是"量级够小", 不是"逐位为 0"。
constexpr double LOWER_TOL = 1e-14;

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

template <CommutativeScalar T>
static double Abs1(const T& x) { return std::sqrt(double(ScalarTraits<T>::abs2(x))); }

// ---- 量度 ------------------------------------------------------------------

// ‖Q^H·Q − I_n‖∞ —— 列正交 / 酉
template <CommutativeScalar T>
static double OrthResid(const Matrix<T>& Q)
{
    return double(NormInf(Q.Adjoint() * Q - Matrix<T>::Identity(Q.Cols())));
}

// ‖Q·R − A‖∞ —— V2.1 判据本体
template <CommutativeScalar T>
static double ReconResid(const Matrix<T>& A, const QRcomposition<T>& qr)
{
    return double(NormInf(qr.Q * qr.R - A));
}

// R 的严格下三角部分最大模 (应该恒为 0)
template <CommutativeScalar T>
static double LowerMax(const Matrix<T>& R)
{
    double best = 0.0;
    for (int i = 1; i < R.Rows(); i++)
        for (int j = 0; j < i && j < R.Cols(); j++)
            best = std::max(best, Abs1(R(i, j)));
    return best;
}

template <typename F>
static bool Throws(F&& f)
{
    try { f(); } catch (...) { return true; }
    return false;
}

// ---- 主判据: 三个实现跑同一套 ----------------------------------------------

template <CommutativeScalar T>
static void RunCase(const Matrix<T>& A, const std::string& tag)
{
    const int m = A.Rows();
    const int n = A.Cols();

    struct Entry { const char* name; QRcomposition<T> (*fn)(const Matrix<T>&); };
    const Entry fns[3] = {
        {"householder", &QR_Householder<T>},
        {"givens",      &QR_givens<T>},
        {"mgs",         &QR_Shimit<T>},
    };

    for (const Entry& e : fns) {
        const std::string who = tag + " [" + e.name + "]";
        try {
            const QRcomposition<T> qr = e.fn(A);

            const bool shape = (qr.Q.Rows() == m && qr.Q.Cols() == n &&
                                qr.R.Rows() == n && qr.R.Cols() == n);
            const double rr = ReconResid(A, qr);
            const double ro = OrthResid(qr.Q);
            const double lo = LowerMax(qr.R);

            ck.check(shape, who + " 尺寸 Q=m×n, R=n×n");
            ck.check(rr < V21_TOL, who + " 重构 ‖QR−A‖∞ = " + Fmt(rr) + " < 1e-8");
            ck.check(ro < ORTH_TOL, who + " Q 列正交 ‖QᴴQ−I‖∞ = " + Fmt(ro) + " < 1e-12");
            ck.check(lo < LOWER_TOL, who + " R 严格下三角 < 1e-14 (实测 " + Fmt(lo) + ")");
        } catch (const std::exception& ex) {
            ck.check(false, who + " 抛异常: " + ex.what());
        }
    }
}

// ---- 解析对拍: 检查 |R| 的对角 (与实现无关的解析值) -------------------------

template <CommutativeScalar T>
static void CheckDiagR(const Matrix<T>& A, const std::vector<double>& want,
                       const std::string& tag)
{
    struct Entry { const char* name; QRcomposition<T> (*fn)(const Matrix<T>&); };
    const Entry fns[3] = {
        {"householder", &QR_Householder<T>},
        {"givens",      &QR_givens<T>},
        {"mgs",         &QR_Shimit<T>},
    };
    for (const Entry& e : fns) {
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
            ck.check(all, tag + " [" + e.name + "] |diag R| = " + got);
        } catch (const std::exception& ex) {
            ck.check(false, tag + " [" + e.name + "] 抛异常: " + ex.what());
        }
    }
}

int main()
{
    std::printf("==== QR 分解验收 (S6b: householder / givens / mgs) ====\n\n");

    std::printf("-- A. 实数: 方阵 / 长方阵 / 退化形状 --\n");
    {
        Matrix<D> A2 = {{1.0, 2.0}, {3.0, 4.0}};
        RunCase(A2, "实 2x2");

        Matrix<D> A32 = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        RunCase(A32, "实 3x2");

        Matrix<D> A53 = RandMatrix<D>(5, 3, 11u);
        RunCase(A53, "实 5x3 随机");

        Matrix<D> A41 = RandMatrix<D>(4, 1, 13u);
        RunCase(A41, "实 4x1");

        Matrix<D> A11 = {{7.0}};
        RunCase(A11, "实 1x1");

        Matrix<D> A55 = RandMatrix<D>(5, 5, 17u);
        RunCase(A55, "实 5x5 随机");

        // 病态: 第 0 列 x1 = -1 为负、且 |x1| ≈ ||x|| (尾部只有 1e-6)。
        // Householder 的 v1 = x1 + phase(x1)*||x|| —— 带 phase 是"加"(不消), 不带就是
        // -1 + 1 = 0 的灾难性相消, 会丢掉 1e-6 那一项、并把误差放大到相对 1e-6。
        // 这是"符号选择"在实数域唯一能显形的场景 (正确性无损, 精度有损)。
        Matrix<D> AIll = {{-1.0, 0.0}, {1e-6, 1.0}};
        RunCase(AIll, "实 2x2 病态 x1<0");
    }

    std::printf("\n-- B. 复数: 方阵 / 长方阵 / 列向量 --\n");
    {
        Matrix<C> B2 = {{C(1.0, 1.0), C(2.0, -1.0)}, {C(3.0, 0.0), C(0.0, 2.0)}};
        RunCase(B2, "复 2x2");

        Matrix<C> B43 = RandMatrix<C>(4, 3, 21u);
        RunCase(B43, "复 4x3 随机");

        Matrix<C> B21 = Matrix<C>(2, 1);
        B21(0, 0) = C(1.0, 2.0);
        B21(1, 0) = C(3.0, -1.0);
        RunCase(B21, "复 2x1");

        // 纯相位列 (实部全 0) —— 实旋转在这里完全无能为力, 逼出相位分支
        Matrix<C> Bph = Matrix<C>(2, 2);
        Bph(0, 0) = C(0.0, 1.0);  Bph(0, 1) = C(0.0, 0.0);
        Bph(1, 0) = C(0.0, -2.0); Bph(1, 1) = C(0.0, 3.0);
        RunCase(Bph, "复 2x2 纯相位");

        // 病态: x1 = -1 (负实数), 尾部 1e-6 —— 同实数域那条, 验证 phase 也吃这一路
        Matrix<C> CIll = {{C(-1.0, 0.0), C(0.0, 1.0)}, {C(1e-6, 0.0), C(1.0, 0.0)}};
        RunCase(CIll, "复 2x2 病态 x1<0");
    }

    std::printf("\n-- C. 解析对拍: |diag R| 该是多少就是多少 --\n");
    {
        // 对角阵: 各列正交, |R_ii| = |a_ii|; 顺序不变 (QR 不排序)
        Matrix<D> Dg = Matrix<D>::Diagonal({3.0, -4.0, 5.0});
        CheckDiagR(Dg, {3.0, 4.0, 5.0}, "实 diag(3,-4,5)");

        Matrix<D> Id = Matrix<D>::Identity(3);
        CheckDiagR(Id, {1.0, 1.0, 1.0}, "实 I3");

        // 2x2 闭式: A = [[1,2],[3,4]] -> |R11| = sqrt(10), |R22| = |det|/|R11| = 2/sqrt(10)
        Matrix<D> A2 = {{1.0, 2.0}, {3.0, 4.0}};
        CheckDiagR(A2, {std::sqrt(10.0), 2.0 / std::sqrt(10.0)}, "实 2x2 闭式");

        // 复对角: |R_ii| = |a_ii|
        Matrix<C> Cd = Matrix<C>::Diagonal({C(0.0, 1.0), C(2.0, 0.0)});
        CheckDiagR(Cd, {1.0, 2.0}, "复 diag(i,2)");

        // 复 2x2: |R11| = ||col0|| = sqrt(2);  det = 2 - (1+i) = 1-i, |det| = sqrt(2),
        // 故 |R22| = |det|/|R11| = sqrt(2)/sqrt(2) = 1
        Matrix<C> C2 = {{C(1.0, 0.0), C(1.0, 1.0)}, {C(1.0, 0.0), C(2.0, 0.0)}};
        CheckDiagR(C2, {std::sqrt(2.0), 1.0}, "复 2x2 闭式");
    }

    std::printf("\n-- D. 秩亏: 列线性相关 --\n");
    {
        Matrix<D> Rk = {{1.0, 2.0}, {2.0, 4.0}};    // col1 = 2 * col0
        RunCase(Rk, "实 2x2 秩亏");

        Matrix<D> Rk0 = {{0.0, 0.0}, {0.0, 0.0}};   // 零矩阵: 全部列都秩亏
        RunCase(Rk0, "实 2x2 零矩阵");
    }

    std::printf("\n-- E. 异常口径 --\n");
    {
        ck.check(Throws([] { (void)QR_Householder(Matrix<D>()); }), "householder 空矩阵 -> throw");
        ck.check(Throws([] { (void)QR_givens(Matrix<D>()); }),      "givens 空矩阵 -> throw");
        ck.check(Throws([] { (void)QR_Shimit(Matrix<D>()); }),      "mgs 空矩阵 -> throw");

        Matrix<D> Wide = RandMatrix<D>(2, 3, 31u);                 // m < n: thin 约定下退化
        ck.check(Throws([&] { (void)QR_Householder(Wide); }), "householder m<n -> throw");
        ck.check(Throws([&] { (void)QR_givens(Wide); }),      "givens m<n -> throw");
        ck.check(Throws([&] { (void)QR_Shimit(Wide); }),      "mgs m<n -> throw");
    }

    return CheckReport(ck.ok, ck.fail, "QR (S6b)");
}
