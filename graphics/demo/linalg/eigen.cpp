// eigen.cpp -- linalg/eigen.hpp 验收 (S3: cyclic Jacobi 迭代 / V1.3 硬判据)
//
// 判据清单 (判据先行: 先写「多少算对」再写实现):
//   ① V1.3 硬指标: 对称阵重构 max|A − Q·Λ·Q^H| < 1e-8; Hermitian 复矩阵过同一份源码
//   ② 特征向量定义: 逐列 A·q_j = λ_j·q_j
//   ③ Q 正交/酉: ‖Q^H·Q − I‖∞ < 1e-12
//   ④ 相似变换不变量: Σλ = tr(A)
//   ⑤ 解析对拍: 2x2 闭式特征值
//   ⑥ 边界: 1x1 / 已对角 / 零矩阵; 异常: 非方阵 / 空 / 非对称 (阶段 2 未实现)
//
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh
//   ./graphics/build/demo/linalg/eigen

#include "../common.hpp"
#include "linalg/eigen.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

static Checker ck;

using D = double;
using C = Complex<double>;

constexpr double V13_TOL = 1e-8;   // V1.3 硬阈值

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

// 取实部 (D7: 一切数值判据只在实数域发生)
template <CommutativeScalar T>
static double Re(const T& x) { return double(ScalarTraits<T>::real(x)); }

// ‖A − Q·Λ·Q^H‖∞ —— V1.3 判据本体
template <CommutativeScalar T>
static double ReconResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    return double(NormInf(A - e.Q * e.Lamda * e.Q.Adjoint()));
}

// ‖Q^H·Q − I‖∞ —— 正交 / 酉
template <CommutativeScalar T>
static double OrthResid(const Matrix<T>& Q)
{
    return double(NormInf(Q.Adjoint() * Q - Matrix<T>::Identity(Q.Rows())));
}

// ‖Q^H·A·Q − Λ‖∞ —— 真的把 A 转成对角了吗
template <CommutativeScalar T>
static double DiagResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    return double(NormInf(e.Q.Adjoint() * A * e.Q - e.Lamda));
}

// max_j ‖A·q_j − λ_j·q_j‖∞ —— 特征对的定义 (比"重构"更强: 逐列认证)
template <CommutativeScalar T>
static double EigenpairResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const int n = A.Rows();
    real_t worst = real_t(0);
    for (int j = 0; j < n; j++) {
        Matrix<T> q(n, 1);
        for (int i = 0; i < n; i++) q(i, 0) = e.Q(i, j);
        Matrix<T> r = A * q;
        const T lam = e.Lamda(j, j);
        for (int i = 0; i < n; i++) r(i, 0) = r(i, 0) - lam * q(i, 0);
        const real_t v = NormInf(r);
        if (v > worst) worst = v;
    }
    return double(worst);
}

// |Σλ − tr(A)| —— 相似变换保迹 (与特征值取值无关的独立校验)
template <CommutativeScalar T>
static double TraceGap(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    double ta = 0.0, tl = 0.0;
    for (int i = 0; i < A.Rows(); i++) {
        ta += Re(A(i, i));
        tl += Re(e.Lamda(i, i));
    }
    return std::abs(ta - tl);
}

// 跑 eig 并把异常记成一条 FAIL —— 实现抛异常不该让整个用例 abort
// (S2a 的教训: 未捕获异常会 terminate 掉用例, 其余判据一条都看不到)。
template <typename T>
static bool RunEig(const Matrix<T>& A, Eigendecomposition<T>& out, const std::string& tag)
{
    try {
        out = eig(A);
        return true;
    } catch (const std::exception& ex) {
        ck.check(false, tag + ": eig 抛异常 -> " + std::string(ex.what()));
        return false;
    }
}

int main()
{
    std::printf("== linalg/eigen.hpp 验收 (S3 / V1.3) ==\n");

    double worst_recon_d = 0.0, worst_recon_c = 0.0;

    // ---- A. 2x2 解析对拍 (闭式特征值) ----
    std::printf("\n[A] 2x2 解析对拍\n");
    {
        auto sorted = [](double x, double y) {
            return (x > y) ? std::make_pair(x, y) : std::make_pair(y, x);
        };
        const double s5 = std::sqrt(5.0);
        const double s2 = std::sqrt(2.0);

        Matrix<D> A1 = {{2.0, 1.0}, {1.0, 3.0}};          // λ = (5±√5)/2
        auto e1 = eig(A1);
        auto p1 = sorted(Re(e1.Lamda(0, 0)), Re(e1.Lamda(1, 1)));
        ck.check(std::abs(p1.first - (5.0 + s5) / 2) < 1e-12
                 && std::abs(p1.second - (5.0 - s5) / 2) < 1e-12,
                 "[[2,1],[1,3]]: λ = (5±√5)/2 = 3.618034 / 1.381966");
        ck.check(ReconResid(A1, e1) < 1e-14, "[[2,1],[1,3]]: 重构残差 < 1e-14");

        Matrix<D> A2 = {{1.0, -1.0}, {-1.0, 3.0}};        // λ = 2±√2
        auto e2 = eig(A2);
        auto p2 = sorted(Re(e2.Lamda(0, 0)), Re(e2.Lamda(1, 1)));
        ck.check(std::abs(p2.first - (2.0 + s2)) < 1e-12
                 && std::abs(p2.second - (2.0 - s2)) < 1e-12,
                 "[[1,-1],[-1,3]]: λ = 2±√2 = 3.414214 / 0.585786");

        Matrix<D> A3 = {{3.0, 0.0}, {0.0, -2.0}};         // 本来就是对角
        auto e3 = eig(A3);
        auto p3 = sorted(Re(e3.Lamda(0, 0)), Re(e3.Lamda(1, 1)));
        ck.check(std::abs(p3.first - 3.0) < 1e-15 && std::abs(p3.second + 2.0) < 1e-15,
                 "已对角 [[3,0],[0,-2]]: λ 不变");
        ck.check(IsApprox(e3.Q, Matrix<D>::Identity(2), 1e-15), "已对角: Q = I");
    }

    // ---- B. 实对称随机矩阵: V1.3 硬指标 ----
    std::printf("\n[B] 实对称随机矩阵 (V1.3)\n");
    for (int n = 2; n <= 6; n++) {
        for (unsigned sd = 1; sd <= 3; sd++) {
            Matrix<D> M = RandMatrix<D>(n, n, 300u + 7u * sd + unsigned(n));
            Matrix<D> S = M + M.Transpose();      // S(i,j) = M(i,j)+M(j,i): 逐位对称

            const std::string tag = std::string("对称 n=") + std::to_string(n)
                                  + " seed=" + std::to_string(sd);
            Eigendecomposition<D> e;
            if (!RunEig(S, e, tag)) continue;     // 异常已记成一条 FAIL

            const double rr = ReconResid(S, e);
            const double ro = OrthResid(e.Q);
            const double rd = DiagResid(S, e);
            const double rp = EigenpairResid(S, e);
            const double rg = TraceGap(S, e);
            if (rr > worst_recon_d) worst_recon_d = rr;
            ck.check(rr < V13_TOL, tag + ": 重构 ‖A − QΛQᵀ‖∞ = " + Fmt(rr) + " < 1e-8");
            ck.check(ro < 1e-12,   tag + ": Q 正交 ‖QᵀQ − I‖∞ = " + Fmt(ro) + " < 1e-12");
            ck.check(rd < 1e-12,   tag + ": 对角化 ‖QᵀAQ − Λ‖∞ = " + Fmt(rd) + " < 1e-12");
            ck.check(rp < 1e-12,   tag + ": 特征对 max‖A·q_j − λ_j·q_j‖∞ = " + Fmt(rp) + " < 1e-12");
            ck.check(rg < 1e-10,   tag + ": 迹不变 |Σλ − tr(A)| = " + Fmt(rg) + " < 1e-10");
        }
    }

    // ---- C. 复 Hermitian 随机矩阵 (同一份源码, V1.1) ----
    std::printf("\n[C] 复 Hermitian 随机矩阵 (同一份 linalg 源码, V1.1)\n");
    for (int n = 2; n <= 4; n++) {
        for (unsigned sd = 1; sd <= 2; sd++) {
            Matrix<C> M = RandMatrix<C>(n, n, 500u + 11u * sd + unsigned(n));
            Matrix<C> H = M + M.Adjoint();        // Hermitian (逐位精确)

            const std::string tag = std::string("Hermitian n=") + std::to_string(n)
                                  + " seed=" + std::to_string(sd);
            Eigendecomposition<C> e;
            if (!RunEig(H, e, tag)) continue;     // 异常已记成一条 FAIL

            const double rr = ReconResid(H, e);
            const double ro = OrthResid(e.Q);
            const double rd = DiagResid(H, e);
            const double rp = EigenpairResid(H, e);
            const double rg = TraceGap(H, e);
            if (rr > worst_recon_c) worst_recon_c = rr;
            ck.check(rr < V13_TOL, tag + ": 重构 ‖A − QΛQᴴ‖∞ = " + Fmt(rr) + " < 1e-8");
            ck.check(ro < 1e-12,   tag + ": Q 酉 ‖QᴴQ − I‖∞ = " + Fmt(ro) + " < 1e-12");
            ck.check(rd < 1e-12,   tag + ": 对角化 ‖QᴴAQ − Λ‖∞ = " + Fmt(rd) + " < 1e-12");
            ck.check(rp < 1e-12,   tag + ": 特征对 max‖A·q_j − λ_j·q_j‖∞ = " + Fmt(rp) + " < 1e-12");
            ck.check(rg < 1e-10,   tag + ": 迹不变 |Σλ − tr(A)| = " + Fmt(rg) + " < 1e-10");
        }
    }

    // ---- D. 复 Hermitian 的特征值必须是实数 (Hermitian 的结构决定) ----
    std::printf("\n[D] 复 Hermitian: 特征值为实数\n");
    {
        Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};   // λ = 4 / 1
        auto e = eig(H);
        double worst_imag = 0.0;
        for (int i = 0; i < H.Rows(); i++) {
            const double im = std::abs(e.Lamda(i, i).imag);
            if (im > worst_imag) worst_imag = im;
        }
        ck.check(worst_imag < 1e-15, "Hermitian 2x2: max|Im λ| = " + Fmt(worst_imag) + " < 1e-15");

        const double a = Re(e.Lamda(0, 0)), b = Re(e.Lamda(1, 1));
        const double hi = (a > b) ? a : b, lo = (a > b) ? b : a;
        ck.check(std::abs(hi - 4.0) < 1e-12 && std::abs(lo - 1.0) < 1e-12,
                 "[[2,1-i],[1+i,3]]: λ = 4 / 1");
    }

    // ---- E. 边界与异常 ----
    std::printf("\n[E] 边界与异常\n");
    {
        Matrix<D> A1 = {{7.0}};
        auto e1 = eig(A1);
        ck.check(std::abs(Re(e1.Lamda(0, 0)) - 7.0) < 1e-15, "1x1: λ = 7");

        Matrix<D> Z = Matrix<D>::Zero(3, 3);
        auto ez = eig(Z);
        ck.check(double(NormInf(ez.Lamda)) < 1e-15, "零矩阵: λ 全 0");
        ck.check(IsApprox(ez.Q, Matrix<D>::Identity(3), 1e-15), "零矩阵: Q = I");

        Matrix<D> Dg = Matrix<D>::Diagonal({3.0, -1.0, 2.0});
        auto ed = eig(Dg);
        ck.check(ReconResid(Dg, ed) < 1e-14, "对角 [[3,-1,2]]: 重构残差 < 1e-14");

        auto throws = [](auto fn) {
            try { fn(); return false; } catch (const std::exception&) { return true; }
        };

        Matrix<D> NS = {{1.0, 2.0}, {0.0, 1.0}};   // 上三角非对称
        ck.check(throws([&] { (void)eig(NS); }),   "非对称 -> throw (阶段 2 未实现)");
        Matrix<D> Rect(2, 3);
        ck.check(throws([&] { (void)eig(Rect); }), "非方阵 -> throw");
        Matrix<D> Emp(0, 0);
        ck.check(throws([&] { (void)eig(Emp); }),  "0x0 -> throw");
    }

    std::printf("\n最差重构残差: double %s | complex %s  (V1.3 判据 < 1e-8)\n",
                Fmt(worst_recon_d).c_str(), Fmt(worst_recon_c).c_str());
    return CheckReport(ck.ok, ck.fail, "eigen");
}
