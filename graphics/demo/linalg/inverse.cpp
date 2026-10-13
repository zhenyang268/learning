// inverse.cpp -- linalg/inverse.hpp 验收 (S2: Inverse 分发; S2a-2: SPD 档接入)
//
// 判据 (matrix.md §10.9 第 2 条: 判据先行):
//   ① V1.2 硬指标: 随机矩阵 ‖A·A⁻¹ − I‖∞ < 1e-8, 且 ‖A⁻¹·A − I‖∞ < 1e-8
//   ② double 与 Complex<double> 共用同一份 linalg 源码 (V1.1)
//   ③ 结构分支: 对角 -> 逐元倒数; 正交(旋转) -> A^T; 单位阵 -> I (解析对拍)
//   ④ Classify 语义: 对角/旋转/随机/对称 分别落到对应 kind
//   ⑤ Inverse(A, hint) 与 Inverse(A) 一致 (hint 只跳过探测, 不改结果)
//   ⑥ 异常: 奇异 / 非方阵 / 0x0
//   ⑦ SPD 档 (S2a-2): Classify 分出 SymmetricPositiveDefinite 走 cholesky_inverse,
//      与 LU 路径对拍 + 正定性判据 (matrix.md §3.2); QR 分支仍显式 throw (待 plan.md S6)
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

        // Hermitian 正定 (复): S2a-2 起探测落 SymmetricPositiveDefinite, 走 cholesky_inverse
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
        Matrix<D> Sy = {{2.0, 1.0}, {1.0, 3.0}};    // 特征值 3.618 / 1.382 -> 正定
        Matrix<D> Id = {{1.0, 2.0}, {2.0, 1.0}};    // 特征值 3 / -1 -> 不定
        Matrix<D> NS = {{1.0, 2.0}, {0.0, 1.0}};   // 上三角非对称

        ck.check(Classify(Dg, tol) == MatrixKind::Diagonal,   "对角 -> Diagonal");
        ck.check(Classify(Ro, tol) == MatrixKind::Orthogonal, "旋转 -> Orthogonal");
        ck.check(Classify(Ge, tol) == MatrixKind::Generic,    "随机 -> Generic");
        ck.check(Classify(Sy, tol) == MatrixKind::SymmetricPositiveDefinite,
                 "实对称正定 -> SymmetricPositiveDefinite (S2a-2 打开该档)");
        ck.check(Classify(Id, tol) == MatrixKind::Symmetric,
                 "实对称非正定 -> Symmetric (没被 SPD 档吞掉)");
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

        // SPD 分支已于 S2a-2 接入, 不再 throw —— 正确性验收在 G 段。
        // 这里只留"hint 声明 SPD 但实际不正定"这条失败路径: 入参错就显式抛, 不静默算错。
        Matrix<D> NotSpd = {{1.0, 2.0}, {2.0, 1.0}};
        ck.check(throws([&] { (void)Inverse(NotSpd, MatrixKind::SymmetricPositiveDefinite); }),
                 "hint=SymmetricPositiveDefinite 但 A 不正定 -> throw");
        ck.check(throws([&] { (void)qr_inverse(Matrix<D>::Identity(2)); }),
                 "qr_inverse -> throw (qr 未接入)");
    }

    // ---- G. SPD 档 (S2a-2: cholesky 接入; matrix.md §3.2 "SPD 检测 = Cholesky 不选主元能否走通") ----
    std::printf("\n[G] SPD 档 (S2a-2)\n");
    {
        const double tol2 = DefaultTol<D>(2);

        // ① 正定性判据: 真 / 假 两侧都测 (只测"能过"的一侧, 等于没测)
        Matrix<D> Spd  = {{2.0, 1.0}, {1.0, 3.0}};          // 特征值 3.618 / 1.382
        Matrix<D> Ind  = {{1.0, 2.0}, {2.0, 1.0}};          // 特征值 3 / -1
        Matrix<D> Neg  = Matrix<D>::Diagonal({-1.0, -2.0}); // 全负
        Matrix<D> Sem  = {{1.0, 2.0}, {2.0, 4.0}};          // 特征值 5 / 0 (半正定, 奇异)
        Matrix<D> Rect(2, 3);

        ck.check(is_positive_definite(Spd),  "[[2,1],[1,3]] (特征值 3.6/1.4) -> true");
        ck.check(!is_positive_definite(Ind), "[[1,2],[2,1]] (特征值 3/-1) -> false");
        ck.check(!is_positive_definite(Neg), "全负对角 -> false");
        ck.check(!is_positive_definite(Sem), "半正定奇异 [[1,2],[2,4]] -> false (主元塌到 0)");
        ck.check(is_positive_definite(Matrix<D>::Identity(4)), "单位阵 -> true");
        ck.check(!is_positive_definite(Rect), "非方阵 -> false (不问 cholesky, 直接答否)");

        // ② Classify 真把 SPD 与 Symmetric 分开 —— 平铺 if 会让 SPD 档永远不可达, 这条是回归闸门
        ck.check(Classify(Spd, tol2) == MatrixKind::SymmetricPositiveDefinite,
                 "Classify(正定对称) -> SymmetricPositiveDefinite");
        ck.check(Classify(Ind, tol2) == MatrixKind::Symmetric,
                 "Classify(非正定对称) -> Symmetric (SPD 档没把它吞掉)");
        ck.check(Classify(Matrix<D>::Identity(3), DefaultTol<D>(3)) == MatrixKind::Diagonal,
                 "Classify(单位阵) -> Diagonal (对角档排在 SPD 之前, 有意为之)");

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

            const std::string tag = std::string("随机 SPD n=") + std::to_string(n);
            ck.check(Classify(S, DefaultTol<D>(n)) == MatrixKind::SymmetricPositiveDefinite,
                     tag + ": Classify -> SymmetricPositiveDefinite");
            ck.check(r < V12_TOL, tag + ": ‖A·A⁻¹ − I‖∞ = " + Fmt(r) + " < 1e-8");
            ck.check(dx < 1e-10, tag + ": cholesky 与 LU 逐元差 = " + Fmt(dx) + " < 1e-10");
        }
        std::printf("    (随机 SPD 最差残差 %s)\n", Fmt(worst_spd).c_str());

        // ④ 复 Hermitian 正定走同一份源码 (别又只测 double)
        Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};   // 特征值 4 / 1
        ck.check(is_positive_definite(H), "复 Hermitian [[2,1-i],[1+i,3]] -> true");
        ck.check(Classify(H, DefaultTol<C>(2)) == MatrixKind::SymmetricPositiveDefinite,
                 "复 Hermitian 正定 -> SymmetricPositiveDefinite");
        Matrix<C> Hc = Inverse(H, MatrixKind::SymmetricPositiveDefinite);
        Matrix<C> Hl = Inverse(H, MatrixKind::Generic);
        ck.check(NormInf(H * Hc - Matrix<C>::Identity(2)) < 1e-12, "复 SPD: cholesky 残差 < 1e-12");
        ck.check(NormInf(Hc - Hl) < 1e-12, "复 SPD: cholesky 与 LU 逐元差 < 1e-12");

        // ⑤ 自动探测入口 (走 Classify) 与 hint=Generic 一致 —— 与 E 段同口径
        Matrix<D> Aut = {{2.0, 1.0}, {1.0, 3.0}};
        ck.check(IsApprox(Inverse(Aut), Inverse(Aut, MatrixKind::Generic), 1e-14),
                 "Inverse(SPD) 自动探测与 hint=Generic 一致");
    }

    std::printf("\n最差残差: double %s | complex %s  (判据 < 1e-8)\n",
                Fmt(worst_d).c_str(), Fmt(worst_c).c_str());
    return CheckReport(ck.ok, ck.fail, "inverse");
}
