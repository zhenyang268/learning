// inverse.cpp -- linalg/inverse.hpp 验收 (S2: Inverse 分发 / V1.2 硬判据)
//
// 判据 (matrix.md §10.9 第 2 条: 判据先行):
//   ① V1.2 硬指标: 随机矩阵 ‖A·A⁻¹ − I‖∞ < 1e-8, 且 ‖A⁻¹·A − I‖∞ < 1e-8
//   ② double 与 Complex<double> 共用同一份 linalg 源码 (V1.1)
//   ③ 结构分支: 对角 -> 逐元倒数; 正交(旋转) -> A^T; 单位阵 -> I (解析对拍)
//   ④ Classify 语义: 对角/旋转/随机/实对称 分别落到对应 kind
//   ⑤ Inverse(A, hint) 与 Inverse(A) 一致 (hint 只跳过探测, 不改结果)
//   ⑥ 异常: 奇异 / 非方阵 / 0x0
//   ⑦ SPD / QR 分支当前显式 throw (未接入 cholesky/qr, plan.md S2a/S6) —— 不求静默算错
//
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh             # 全量
//   ./graphics/build/demo/linalg/inverse # 单跑

#include "../common.hpp"
#include "linalg/inverse.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

static Checker ck;

using D = double;
using C = Complex<double>;

constexpr double V12_TOL = 1e-8;   // V1.2 硬阈值

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

// ‖A·X − I‖∞ (元素最大模; demo/common.hpp 的 NormInf 就是这个口径)
template <CommutativeScalar T>
static typename ScalarTraits<T>::real_t InvResid(const Matrix<T>& A, const Matrix<T>& X)
{
    return NormInf(A * X - Matrix<T>::Identity(A.Rows()));
}

template <CommutativeScalar T>
static void CheckRandomInv(const char* tag, int n, unsigned seed, double& worst)
{
    using real_t = typename ScalarTraits<T>::real_t;
    Matrix<T> A = RandMatrix<T>(n, n, seed);
    Matrix<T> X = Inverse(A);
    const real_t r1 = InvResid(A, X);
    const real_t r2 = InvResid(X, A);
    const real_t r  = (r1 > r2) ? r1 : r2;
    if (double(r) > worst) worst = double(r);
    ck.check(r < real_t(V12_TOL),
             std::string(tag) + " n=" + std::to_string(n) + " seed=" + std::to_string(seed)
             + ": ‖A·A⁻¹ − I‖∞ = " + Fmt(double(r)) + " < 1e-8");
}

int main()
{
    std::printf("== linalg/inverse.hpp 验收 (S2 / V1.2) ==\n");

    // ---- A. double 随机矩阵 ----
    std::printf("\n[A] double 随机矩阵 (V1.2 硬指标)\n");
    double worst_d = 0.0;
    for (int n = 2; n <= 6; n++)
        for (unsigned sd = 1; sd <= 4; sd++)
            CheckRandomInv<D>("double", n, sd, worst_d);

    // ---- B. Complex 随机矩阵 (同一份源码) ----
    std::printf("\n[B] Complex<double> 随机矩阵 (同一份 linalg 源码, V1.1)\n");
    double worst_c = 0.0;
    for (int n = 2; n <= 4; n++)
        for (unsigned sd = 1; sd <= 3; sd++)
            CheckRandomInv<C>("complex", n, sd, worst_c);

    // ---- C. 结构分支: 解析对拍 ----
    std::printf("\n[C] 结构分支 (解析逆对拍)\n");
    {
        Matrix<D> Dg = Matrix<D>::Diagonal({2.0, -4.0, 0.5});
        Matrix<D> Di = Inverse(Dg);
        ck.check(std::abs(Di(0, 0) - 0.5) < 1e-15
                 && std::abs(Di(1, 1) + 0.25) < 1e-15
                 && std::abs(Di(2, 2) - 2.0) < 1e-15,
                 "对角 [[2,-4,0.5]]: 逆 = 逐元倒数");

        const double th = 0.7;
        Matrix<D> R = {{std::cos(th), -std::sin(th)}, {std::sin(th), std::cos(th)}};
        Matrix<D> Ri = Inverse(R);
        ck.check(IsApprox(Ri, R.Transpose(), 1e-12), "旋转 2x2: 逆 = A^T");
        ck.check(NormInf(R * Ri - Matrix<D>::Identity(2)) < 1e-12, "旋转 2x2: 残差 < 1e-12");

        Matrix<D> I3 = Matrix<D>::Identity(3);
        ck.check(IsApprox(Inverse(I3), I3, 1e-15), "单位阵: 逆 = 自身");

        // Hermitian 正定 (复): 探测应落 Symmetric, 走 Generic(LU)
        Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};
        Matrix<C> Hi = Inverse(H);
        ck.check(NormInf(H * Hi - Matrix<C>::Identity(2)) < 1e-12, "Hermitian 2x2 (复): 残差 < 1e-12");
    }

    // ---- D. Classify 语义 ----
    std::printf("\n[D] Classify 语义\n");
    {
        const double tol = DefaultTol<D>(3);
        Matrix<D> Dg = Matrix<D>::Diagonal({1.0, 2.0, 3.0});
        Matrix<D> Ro = {{std::cos(0.3), -std::sin(0.3)}, {std::sin(0.3), std::cos(0.3)}};
        Matrix<D> Ge = RandMatrix<D>(3, 3, 7u);
        Matrix<D> Sy = {{2.0, 1.0}, {1.0, 3.0}};
        Matrix<D> NS = {{1.0, 2.0}, {0.0, 1.0}};   // 上三角非对称

        ck.check(Classify(Dg, tol) == MatrixKind::Diagonal,   "对角 -> Diagonal");
        ck.check(Classify(Ro, tol) == MatrixKind::Orthogonal, "旋转 -> Orthogonal");
        ck.check(Classify(Ge, tol) == MatrixKind::Generic,    "随机 -> Generic");
        ck.check(Classify(Sy, tol) == MatrixKind::Symmetric,  "实对称 -> Symmetric");
        ck.check(Classify(NS, tol) != MatrixKind::Symmetric,  "上三角非对称 -> 不判 Symmetric");

        // A = [[1,c],[0,s]] 满足 A^T·A = [[1,c],[c,1]]: 对角线全 1, 非对角非零。
        // 只查对角线的 is_orthogonal 会把它误判成正交, 全元素检查则不会。
        Matrix<D> FO = {{1.0, 0.5}, {0.0, std::sqrt(1.0 - 0.25)}};
        ck.check(Classify(FO, tol) != MatrixKind::Orthogonal, "A^T·A 对角全 1 但非对角非零 -> 不判 Orthogonal");
    }

    // ---- E. hint 一致性 ----
    std::printf("\n[E] Inverse(A, hint) 与 Inverse(A) 一致\n");
    for (int n = 3; n <= 5; n++) {
        Matrix<D> A  = RandMatrix<D>(n, n, 11u);
        Matrix<D> X1 = Inverse(A);
        Matrix<D> X2 = Inverse(A, MatrixKind::Generic);
        ck.check(IsApprox(X1, X2, 1e-14),
                 std::string("n=") + std::to_string(n) + ": hint=Generic 与自动探测一致");
    }

    // ---- F. 异常路径 ----
    std::printf("\n[F] 异常路径\n");
    {
        auto throws = [](auto fn) {
            try { fn(); return false; } catch (const std::exception&) { return true; }
        };

        Matrix<D> singular = {{1.0, 2.0}, {2.0, 4.0}};
        ck.check(throws([&] { (void)Inverse(singular); }), "奇异矩阵 -> throw");

        Matrix<D> rect(2, 3);
        ck.check(throws([&] { (void)Inverse(rect); }), "非方阵 -> throw");

        Matrix<D> empty(0, 0);
        ck.check(throws([&] { (void)Inverse(empty); }), "0x0 -> throw");

        ck.check(throws([&] { (void)Inverse(Matrix<D>::Identity(2), MatrixKind::SymmetricPositiveDefinite); }),
                 "SPD 分支 -> throw (cholesky 未接入)");
        ck.check(throws([&] { (void)qr_inverse(Matrix<D>::Identity(2)); }),
                 "qr_inverse -> throw (qr 未接入)");
    }

    std::printf("\n最差残差: double %s | complex %s  (判据 < 1e-8)\n",
                Fmt(worst_d).c_str(), Fmt(worst_c).c_str());
    return CheckReport(ck.ok, ck.fail, "inverse");
}
