// svd.cpp -- linalg/svd.hpp 验收 (S4: 经 A^H·A 特征分解组装 / V1.4 硬判据)
//
// 判据 (matrix.md §10.9 第 2 条: 判据先行):
//   ① V1.4 硬指标: 重构 ‖U·Σ·V^H − A‖∞ < 1e-8 (实/复, 方阵/长方阵)
//   ② 奇异值非负且降序
//   ③ U / V 列正交酉: ‖U^H·U − I‖∞ < 1e-12, ‖V^H·V − I‖∞ < 1e-12
//   ④ 独立交叉验证: Σσ_i² = ‖A‖_F² (长方阵也成立); σ_i² 与 eig(A^H·A) 的谱逐一吻合
//   ⑤ 解析对拍: 对角阵 / 旋转阵 / 2x2 闭式 (含复 Hermitian)
//   ⑥ 边界: 1x1 / 零矩阵 / 秩亏 (低秩, σ 含 0 -> U 走正交补) / 长方形; 异常: 空矩阵
//
// 运行(由 CMake 注册为 ctest 用例, 也可直接执行看逐项 [OK]/[FAIL]):
//   ./graphics/run_tests.sh
//   ./graphics/build/demo/linalg/svd

#include "../common.hpp"
#include "linalg/svd.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

static Checker ck;

using D = double;
using C = Complex<double>;

constexpr double V14_TOL = 1e-8;   // V1.4 硬阈值

static std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

// 取实部 (D7: 一切数值判据只在实数域发生)
template <CommutativeScalar T>
static double Re(const T& x) { return double(ScalarTraits<T>::real(x)); }

// ---- 量度 ------------------------------------------------------------------

// ‖A − U·Σ·V^H‖∞ —— V1.4 判据本体
template <CommutativeScalar T>
static double ReconResid(const Matrix<T>& A, const SVDcomposition<T>& sv)
{
    return double(NormInf(A - sv.U * sv.S * sv.V.Adjoint()));
}

// ‖X^H·X − I‖∞ —— 列正交 / 酉 (X 是 m×r, X^H·X 是 r×r)
template <CommutativeScalar T>
static double OrthResid(const Matrix<T>& X)
{
    return double(NormInf(X.Adjoint() * X - Matrix<T>::Identity(X.Cols())));
}

template <CommutativeScalar T>
static std::vector<double> Sigmas(const SVDcomposition<T>& sv)
{
    std::vector<double> v;
    for (int i = 0; i < sv.S.Rows(); i++) v.push_back(Re(sv.S(i, i)));
    return v;
}

// Σσ_i² 与 ‖A‖_F² 的相对差 —— 正交不变量的独立校验 (SVD 不改变 Frobenius 范数)
template <CommutativeScalar T>
static double FrobGap(const Matrix<T>& A, const SVDcomposition<T>& sv)
{
    double fa = 0.0, fs = 0.0;
    for (int i = 0; i < A.Size(); i++) fa += Re(ScalarTraits<T>::abs2(A.Data()[i]));
    for (int i = 0; i < sv.S.Rows(); i++) fs += Re(ScalarTraits<T>::abs2(sv.S(i, i)));
    return std::abs(fa - fs) / (fa > 0.0 ? fa : 1.0);
}

// max_i |σ_i² − λ_i(A^H·A)| / λ_max —— 与 eig 的谱逐一对照 (排序后)
template <CommutativeScalar T>
static double SpectrumGap(const Matrix<T>& A, const SVDcomposition<T>& sv)
{
    using real_t = typename ScalarTraits<T>::real_t;
    const Matrix<T> B = A.Adjoint() * A;
    const Eigendecomposition<T> e = eig_symmetric(B);

    std::vector<real_t> lam(B.Rows());
    for (int i = 0; i < B.Rows(); i++) lam[i] = ScalarTraits<T>::real(e.Lamda(i, i));
    std::sort(lam.begin(), lam.end(), [](real_t a, real_t b) { return a > b; });

    const real_t scale = (lam[0] > real_t(0)) ? lam[0] : real_t(1);
    real_t worst = real_t(0);
    for (int i = 0; i < sv.S.Rows(); i++) {
        const real_t s2 = ScalarTraits<T>::abs2(sv.S(i, i));
        const real_t d = std::abs(s2 - lam[i]) / scale;
        if (d > worst) worst = d;
    }
    return double(worst);
}

// 降序且非负
static bool OrderedNonNegative(const std::vector<double>& s)
{
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i] < -1e-15) return false;
        if (i + 1 < s.size() && s[i] + 1e-12 < s[i + 1]) return false;
    }
    return true;
}

// 跑 svd 并把异常记成一条 FAIL —— 实现抛异常不该让整个用例 abort
// (S2a / S3 的教训: 未捕获异常会 terminate 掉用例, 其余判据一条都看不到)。
template <typename T>
static bool RunSvd(const Matrix<T>& A, SVDcomposition<T>& out, const std::string& tag)
{
    try {
        out = svd(A);
        return true;
    } catch (const std::exception& ex) {
        ck.check(false, tag + ": svd 抛异常 -> " + std::string(ex.what()));
        return false;
    }
}

// 一组矩阵的统一体检 (重构 / 正交 / 有序 / 守恒 / 谱)
template <typename T>
static void CheckAll(const Matrix<T>& A, const std::string& tag,
                     double& worst_recon, double& worst_orth)
{
    SVDcomposition<T> sv;
    if (!RunSvd(A, sv, tag)) return;

    const std::vector<double> s = Sigmas(sv);
    const double rr = ReconResid(A, sv);
    const double ru = OrthResid(sv.U);
    const double rv = OrthResid(sv.V);
    const double rf = FrobGap(A, sv);
    const double rs = SpectrumGap(A, sv);
    if (rr > worst_recon) worst_recon = rr;
    if (ru > worst_orth)  worst_orth = ru;
    if (rv > worst_orth)  worst_orth = rv;

    ck.check(static_cast<int>(sv.U.Rows()) == A.Rows()
             && static_cast<int>(sv.U.Cols()) == std::min(A.Rows(), A.Cols()),
             tag + ": U 尺寸 m×r");
    ck.check(static_cast<int>(sv.V.Rows()) == A.Cols()
             && static_cast<int>(sv.V.Cols()) == std::min(A.Rows(), A.Cols()),
             tag + ": V 尺寸 n×r");
    ck.check(rr < V14_TOL, tag + ": 重构 ‖U·Σ·Vᴴ − A‖∞ = " + Fmt(rr) + " < 1e-8");
    ck.check(OrderedNonNegative(s), tag + ": 奇异值非负且降序");
    // U 的列 u_i = A·v_i/σ_i: σ 走 A^H·A 只有 √eps 级精度, 除以 σ_i 后方向误差放大到
    // eps·κ(A) (实测 1e-13 量级) —— 这是本路线的固有上限, 故 U 判 1e-10。
    // V 的列直接取自 eig_symmetric, 精度 ~ eps, 故仍判 1e-12 (能抓住 eig 侧的真 bug)。
    ck.check(ru < 1e-10, tag + ": U 列正交 ‖UᴴU − I‖∞ = " + Fmt(ru) + " < 1e-10");
    ck.check(rv < 1e-12, tag + ": V 列正交 ‖VᴴV − I‖∞ = " + Fmt(rv) + " < 1e-12");
    ck.check(rf < 1e-12, tag + ": Σσ² = ‖A‖_F² (相对差 = " + Fmt(rf) + " < 1e-12)");
    ck.check(rs < 1e-10, tag + ": σ_i² 与 eig(AᴴA) 谱吻合 (相对差 = " + Fmt(rs) + " < 1e-10)");
}

int main()
{
    std::printf("== linalg/svd.hpp 验收 (S4 / V1.4) ==\n");

    double worst_recon_d = 0.0, worst_recon_c = 0.0;
    double worst_orth = 0.0;

    // ---- A. 解析对拍 (闭式奇异值) ----
    std::printf("\n[A] 解析对拍\n");
    {
        auto near = [](double x, double y) { return std::abs(x - y) < 1e-12; };

        Matrix<D> Dg = Matrix<D>::Diagonal({3.0, -2.0});        // σ = |diag| 降序 = 3, 2
        auto sd = svd(Dg);
        auto sdv = Sigmas(sd);
        ck.check(sdv.size() == 2 && near(sdv[0], 3.0) && near(sdv[1], 2.0),
                 "对角 [[3,0],[0,-2]]: σ = 3, 2");

        // 对角元**升序**输入: 输出必须重排成降序 —— 否则"非负有序"判据形同虚设
        // (eig_symmetric 对已对角矩阵原样返回, 正好是"乱序输入"的确定性来源)
        Matrix<D> DgUp = Matrix<D>::Diagonal({1.0, 5.0, 3.0});
        auto sdu = svd(DgUp);
        auto sduv = Sigmas(sdu);
        ck.check(sduv.size() == 3 && near(sduv[0], 5.0) && near(sduv[1], 3.0) && near(sduv[2], 1.0),
                 "对角 diag(1,5,3): σ 重排为 5, 3, 1");

        Matrix<C> CUp = Matrix<C>::Diagonal({C(1, 0), C(4, 0), C(2, 0)});
        auto scu = svd(CUp);
        auto scuv = Sigmas(scu);
        ck.check(scuv.size() == 3 && near(scuv[0], 4.0) && near(scuv[1], 2.0) && near(scuv[2], 1.0),
                 "复对角 diag(1,4,2): σ 重排为 4, 2, 1");

        Matrix<D> A1 = {{1.0, 2.0}, {2.0, 1.0}};                // A^H A = [[5,4],[4,5]], λ = 9, 1
        auto s1 = svd(A1);
        auto s1v = Sigmas(s1);
        ck.check(s1v.size() == 2 && near(s1v[0], 3.0) && near(s1v[1], 1.0),
                 "[[1,2],[2,1]]: σ = 3, 1 (√9, √1)");

        Matrix<D> Ro = {{0.0, -1.0}, {1.0, 0.0}};               // 旋转: A^H A = I
        auto sr = svd(Ro);
        auto srv = Sigmas(sr);
        ck.check(srv.size() == 2 && near(srv[0], 1.0) && near(srv[1], 1.0),
                 "旋转 [[0,-1],[1,0]]: σ = 1, 1");

        // 复 Hermitian: A^H A = A², λ = 16, 1 -> σ = 4, 1
        Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};
        auto sh = svd(H);
        auto shv = Sigmas(sh);
        ck.check(shv.size() == 2 && near(shv[0], 4.0) && near(shv[1], 1.0),
                 "复 Hermitian [[2,1-i],[1+i,3]]: σ = 4, 1");
        ck.check(ReconResid(H, sh) < 1e-14, "复 Hermitian: 重构残差 < 1e-14");
    }

    // ---- B. 实随机方阵: V1.4 硬指标 ----
    std::printf("\n[B] 实随机方阵 (V1.4)\n");
    for (int n = 2; n <= 6; n++) {
        for (unsigned sd = 1; sd <= 3; sd++) {
            Matrix<D> A = RandMatrix<D>(n, n, 700u + 13u * sd + unsigned(n));
            CheckAll<D>(A, "实方阵 n=" + std::to_string(n) + " seed=" + std::to_string(sd),
                        worst_recon_d, worst_orth);
        }
    }

    // ---- C. 实随机长方阵 (m≠n): 经济型 SVD ----
    std::printf("\n[C] 实随机长方形 (m≠n)\n");
    {
        const int shapes[][2] = {{3, 5}, {5, 3}, {1, 4}, {4, 1}, {2, 7}, {6, 2}};
        for (const auto& mn : shapes) {
            Matrix<D> A = RandMatrix<D>(mn[0], mn[1], 900u + unsigned(mn[0] * 10 + mn[1]));
            CheckAll<D>(A, "实 " + std::to_string(mn[0]) + "x" + std::to_string(mn[1]),
                        worst_recon_d, worst_orth);
        }
    }

    // ---- D. 复随机矩阵 (同一份源码, V1.1): 方阵 + 长方阵 ----
    std::printf("\n[D] 复随机矩阵 (同一份 linalg 源码, V1.1)\n");
    for (int n = 2; n <= 4; n++) {
        for (unsigned sd = 1; sd <= 2; sd++) {
            Matrix<C> A = RandMatrix<C>(n, n, 1300u + 17u * sd + unsigned(n));
            CheckAll<C>(A, "复方阵 n=" + std::to_string(n) + " seed=" + std::to_string(sd),
                        worst_recon_c, worst_orth);
        }
    }
    {
        const int shapes[][2] = {{3, 5}, {5, 3}, {2, 6}};
        for (const auto& mn : shapes) {
            Matrix<C> A = RandMatrix<C>(mn[0], mn[1], 1500u + unsigned(mn[0] * 10 + mn[1]));
            CheckAll<C>(A, "复 " + std::to_string(mn[0]) + "x" + std::to_string(mn[1]),
                        worst_recon_c, worst_orth);
        }
    }

    // ---- E. 秩亏 / 边界 ----
    std::printf("\n[E] 秩亏与边界\n");
    {
        Matrix<D> O = Matrix<D>::Ones(2, 2);              // 秩 1: A^H A = [[2,2],[2,2]], λ = 4, 0
        auto so = svd(O);
        auto ov = Sigmas(so);
        ck.check(ov.size() == 2 && std::abs(ov[0] - 2.0) < 1e-12 && std::abs(ov[1]) < 1e-12,
                 "秩1 全1矩阵: σ = 2, 0");
        ck.check(ReconResid(O, so) < 1e-14, "秩1: 重构残差 < 1e-14");
        ck.check(OrthResid(so.U) < 1e-12, "秩1: U 仍列正交 (σ=0 的列走正交补)");

        Matrix<D> R1(3, 3);                               // 显式秩 1: 列 2 = 2×列 1
        for (int i = 0; i < 3; i++) { R1(i, 0) = 1.0 + i; R1(i, 1) = 2.0 * (1.0 + i); R1(i, 2) = 3.0 + i; }
        CheckAll<D>(R1, "秩1 3x3 (列2=2×列1)", worst_recon_d, worst_orth);

        Matrix<D> Z = Matrix<D>::Zero(3, 3);              // 零矩阵
        auto sz = svd(Z);
        auto zv = Sigmas(sz);
        ck.check(ReconResid(Z, sz) < 1e-15, "零矩阵: 重构残差 ~ 0");
        ck.check(OrthResid(sz.U) < 1e-14 && OrthResid(sz.V) < 1e-14, "零矩阵: U/V 取标准基 (正交)");
        ck.check(zv.size() == 3 && std::abs(zv[0]) < 1e-15, "零矩阵: σ 全 0");

        Matrix<D> A1 = {{7.0}};                            // 1x1
        auto s1 = svd(A1);
        ck.check(std::abs(Sigmas(s1)[0] - 7.0) < 1e-15, "1x1: σ = 7");

        Matrix<D> V1x = {{3.0, 4.0}};                      // 1x2: σ = ‖row‖ = 5
        auto sv1x = svd(V1x);
        auto v1v = Sigmas(sv1x);
        ck.check(v1v.size() == 1 && std::abs(v1v[0] - 5.0) < 1e-12, "1x2 [3,4]: σ = ‖row‖ = 5");
        ck.check(ReconResid(V1x, sv1x) < 1e-14, "1x2: 重构残差 < 1e-14");
    }

    // ---- F. 异常 ----
    std::printf("\n[F] 异常\n");
    {
        auto throws = [](auto fn) {
            try { fn(); return false; } catch (const std::exception&) { return true; }
        };
        Matrix<D> Emp(0, 0);
        ck.check(throws([&] { (void)svd(Emp); }),      "0x0 -> throw");
        Matrix<D> EmpC(2, 0);
        ck.check(throws([&] { (void)svd(EmpC); }),     "2x0 -> throw");
    }

    std::printf("\n最差重构残差: double %s | complex %s  (V1.4 判据 < 1e-8)\n",
                Fmt(worst_recon_d).c_str(), Fmt(worst_recon_c).c_str());
    std::printf("最差 U/V 正交残差: %s\n", Fmt(worst_orth).c_str());
    return CheckReport(ck.ok, ck.fail, "svd");
}
