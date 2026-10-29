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
// ---------------------------------------------------------------------------
// 框架: googletest (2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · 老写法的 `main()` + `Checker` 计分板**已删除** (GTest::gtest_main 提供 main)。
//   · 迁移前 119 项断言 -> 迁移后 119 条 EXPECT/ADD_FAILURE (逐条对应, 一条没合并)。
//   · 段号 [A]–[E] 原样保留, 通过 SCOPED_TRACE 携带。
//   · 老 helper `RunEig` 的 catch 改成 `ADD_FAILURE()`: 一条失败 + `continue`
//     (gtest 若让异常逃出测试体, 该 TEST 剩余部分全不跑, 与老口径不同)。
//   · 异常判据走 EXPECT_EXACT_THROW: 非对称是 `logic_error` (阶段 2 未实现),
//     非方阵 / 空是 `invalid_argument` —— 判基类分不出这两者。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值                     | 边界 / 病态输入                    | 专抓什么错
//  ----|--------------------------|------------------------------------|--------------------------
//   A  | 1e-12 (λ) / 1e-14 (重构)| 2×2 闭式三例: 一般 / 负非对角 / 已对角| 排序假设 (λ 顺序不该被假设)
//   B  | 1e-8 重构; 1e-12 其余   | n=2..6 × 3 seed, S = M + Mᵀ        | 只判重构 (不判逐列特征对)
//   C  | 同上                     | Complex Hermitian, n=2..4 × 2 seed | 复域漏 conj / ᴴ 写成 ᵀ
//   D  | 1e-15 (虚部) / 1e-12    | Hermitian 2×2 (λ = 4/1)            | 复特征值带虚部仍"重构对"
//   E  | 1e-15 / 1e-14           | 1×1 / 3×3 零矩阵 / 对角 / 非对称 / 非方阵| 空矩阵当错误码、对称探测漏判
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **重构 ‖A − QΛQᵀ‖ 抓不到"Q 的两列互换 + Λ 对角元没跟着换"**? —— 抓得到
//     (那样 QΛQᵀ 就不等于 A)。但**抓不到**「Q 被整体乘一个置换」这类只在
//     特征值排序上不同的错: 所以 [A] 只比**排序后**的 λ, 不比顺序本身。
//  ② **Q 正交判据看不见特征向量符号**: q_j 与 −q_j 同样正交、同样重构。
//     符号错误只有 [B]/[C] 的**逐列特征对**判据 (A·q_j = λ_j q_j) 能抓。
//  ③ **Jacobi 不收敛** (把 kMaxSweep 改成 1): [B]/[C] 的重构与对角化判据全红,
//     而 [A] 的 2×2 小例子可能仍绿 —— 小规模先收敛。
//  ④ **[D] 的"特征值为实数"是唯一能抓"复域用了 ᵀ 而不是 ᴴ"的**:
//     漏 conj 时重构判据在 Hermitian 上也可能凑得很近, 但虚部会立刻非 0。
//
// ---------------------------------------------------------------------------
// 运行:
//   ./graphics/run_tests.sh
//   ./graphics/build/demo/linalg/eigen
//   ./graphics/build/demo/linalg/eigen --gtest_filter='*Hermitian*'
// ---------------------------------------------------------------------------

#include "gtest_common.hpp"

#include "linalg/eigen.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

using D = double;
using C = Complex<double>;

constexpr double V13_TOL = 1e-8;   // V1.3 硬阈值

std::string Fmt(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return std::string(buf);
}

// 取实部 (D7: 一切数值判据只在实数域发生)
template <CommutativeScalar T>
double Re(const T& x) { return double(ScalarTraits<T>::real(x)); }

// ‖A − Q·Λ·Q^H‖∞ —— V1.3 判据本体
template <CommutativeScalar T>
double ReconResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    return double(NormInf(A - e.Q * e.Lamda * e.Q.Adjoint()));
}

// ‖Q^H·Q − I‖∞ —— 正交 / 酉
template <CommutativeScalar T>
double OrthResid(const Matrix<T>& Q)
{
    return double(NormInf(Q.Adjoint() * Q - Matrix<T>::Identity(Q.Rows())));
}

// ‖Q^H·A·Q − Λ‖∞ —— 真的把 A 转成对角了吗
template <CommutativeScalar T>
double DiagResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
{
    return double(NormInf(e.Q.Adjoint() * A * e.Q - e.Lamda));
}

// max_j ‖A·q_j − λ_j·q_j‖∞ —— 特征对的定义 (比"重构"更强: 逐列认证)
template <CommutativeScalar T>
double EigenpairResid(const Matrix<T>& A, const Eigendecomposition<T>& e)
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
double TraceGap(const Matrix<T>& A, const Eigendecomposition<T>& e)
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
bool RunEig(const Matrix<T>& A, Eigendecomposition<T>& out, const std::string& tag)
{
    try {
        out = eig(A);
        return true;
    } catch (const std::exception& ex) {
        ADD_FAILURE() << tag << ": eig 抛异常 -> " << ex.what();
        return false;
    }
}

// B/C 段共用的五条判据 (每 seed 一组)。dense = true 时用 ᵀ 的措辞, false 用 ᴴ。
template <typename T>
void CheckSymmetricCase(const Matrix<T>& A, const std::string& tag, bool hermitian, double& worst)
{
    Eigendecomposition<T> e;
    if (!RunEig(A, e, tag)) return;     // 异常已记成一条 FAIL

    const double rr = ReconResid(A, e);
    const double ro = OrthResid(e.Q);
    const double rd = DiagResid(A, e);
    const double rp = EigenpairResid(A, e);
    const double rg = TraceGap(A, e);
    if (rr > worst) worst = rr;

    { SCOPED_TRACE(tag + (hermitian ? ": 重构 < 1e-8" : ": 重构 < 1e-8"));
      EXPECT_LT(rr, V13_TOL) << Fmt(rr); }
    { SCOPED_TRACE(tag + (hermitian ? ": Q 酉 ‖QᴴQ − I‖∞ < 1e-12" : ": Q 正交 ‖QᵀQ − I‖∞ < 1e-12"));
      EXPECT_LT(ro, 1e-12) << Fmt(ro); }
    { SCOPED_TRACE(tag + (hermitian ? ": 对角化 ‖QᴴAQ − Λ‖∞ < 1e-12" : ": 对角化 ‖QᵀAQ − Λ‖∞ < 1e-12"));
      EXPECT_LT(rd, 1e-12) << Fmt(rd); }
    { SCOPED_TRACE(tag + ": 特征对 max‖A·q_j − λ_j·q_j‖∞ < 1e-12");
      EXPECT_LT(rp, 1e-12) << Fmt(rp); }
    { SCOPED_TRACE(tag + ": 迹不变 |Σλ − tr(A)| < 1e-10");
      EXPECT_LT(rg, 1e-10) << Fmt(rg); }
}

}  // namespace

#define STEP(label) SCOPED_TRACE(label)

// ---- [A] 2x2 解析对拍 (闭式特征值) ----
TEST(Eigen, Analytical2x2)
{
    auto sorted = [](double x, double y) {
        return (x > y) ? std::make_pair(x, y) : std::make_pair(y, x);
    };
    const double s5 = std::sqrt(5.0);
    const double s2 = std::sqrt(2.0);

    { STEP("[A] [[2,1],[1,3]]: λ = (5±√5)/2 = 3.618034 / 1.381966");
      Matrix<D> A1 = {{2.0, 1.0}, {1.0, 3.0}};          // λ = (5±√5)/2
      auto e1 = eig(A1);
      auto p1 = sorted(Re(e1.Lamda(0, 0)), Re(e1.Lamda(1, 1)));
      EXPECT_LT(std::abs(p1.first - (5.0 + s5) / 2), 1e-12) << "λ1 = " << p1.first;
      EXPECT_LT(std::abs(p1.second - (5.0 - s5) / 2), 1e-12) << "λ2 = " << p1.second; }

    { STEP("[A] [[2,1],[1,3]]: 重构残差 < 1e-14");
      Matrix<D> A1 = {{2.0, 1.0}, {1.0, 3.0}};
      auto e1 = eig(A1);
      EXPECT_LT(ReconResid(A1, e1), 1e-14) << Fmt(ReconResid(A1, e1)); }

    { STEP("[A] [[1,-1],[-1,3]]: λ = 2±√2 = 3.414214 / 0.585786");
      Matrix<D> A2 = {{1.0, -1.0}, {-1.0, 3.0}};        // λ = 2±√2
      auto e2 = eig(A2);
      auto p2 = sorted(Re(e2.Lamda(0, 0)), Re(e2.Lamda(1, 1)));
      EXPECT_LT(std::abs(p2.first - (2.0 + s2)), 1e-12) << "λ1 = " << p2.first;
      EXPECT_LT(std::abs(p2.second - (2.0 - s2)), 1e-12) << "λ2 = " << p2.second; }

    Matrix<D> A3 = {{3.0, 0.0}, {0.0, -2.0}};         // 本来就是对角
    auto e3 = eig(A3);
    auto p3 = sorted(Re(e3.Lamda(0, 0)), Re(e3.Lamda(1, 1)));

    { STEP("[A] 已对角 [[3,0],[0,-2]]: λ 不变");
      EXPECT_LT(std::abs(p3.first - 3.0), 1e-15);
      EXPECT_LT(std::abs(p3.second + 2.0), 1e-15); }

    { STEP("[A] 已对角: Q = I");
      EXPECT_TRUE(IsApprox(e3.Q, Matrix<D>::Identity(2), 1e-15)); }
}

// ---- [B] 实对称随机矩阵: V1.3 硬指标 ----
TEST(Eigen, SymmetricRandom)
{
    double worst_recon_d = 0.0;
    for (int n = 2; n <= 6; n++) {
        for (unsigned sd = 1; sd <= 3; sd++) {
            Matrix<D> M = RandMatrix<D>(n, n, 300u + 7u * sd + unsigned(n));
            Matrix<D> S = M + M.Transpose();      // S(i,j) = M(i,j)+M(j,i): 逐位对称
            const std::string tag = "[B] 对称 n=" + std::to_string(n)
                                  + " seed=" + std::to_string(sd);
            CheckSymmetricCase(S, tag, false, worst_recon_d);
        }
    }
    std::printf("[eigen] 重构 double 最差 %s\n", Fmt(worst_recon_d).c_str());
}

// ---- [C] 复 Hermitian 随机矩阵 (同一份源码, V1.1) ----
TEST(Eigen, HermitianRandom)
{
    double worst_recon_c = 0.0;
    for (int n = 2; n <= 4; n++) {
        for (unsigned sd = 1; sd <= 2; sd++) {
            Matrix<C> M = RandMatrix<C>(n, n, 500u + 11u * sd + unsigned(n));
            Matrix<C> H = M + M.Adjoint();        // Hermitian (逐位精确)
            const std::string tag = "[C] Hermitian n=" + std::to_string(n)
                                  + " seed=" + std::to_string(sd);
            CheckSymmetricCase(H, tag, true, worst_recon_c);
        }
    }
    std::printf("[eigen] 重构 complex 最差 %s\n", Fmt(worst_recon_c).c_str());
}

// ---- [D] 复 Hermitian 的特征值必须是实数 (Hermitian 的结构决定) ----
TEST(Eigen, HermitianEigenvaluesReal)
{
    Matrix<C> H = {{C(2, 0), C(1, -1)}, {C(1, 1), C(3, 0)}};   // λ = 4 / 1
    auto e = eig(H);

    double worst_imag = 0.0;
    for (int i = 0; i < H.Rows(); i++) {
        const double im = std::abs(e.Lamda(i, i).imag);
        if (im > worst_imag) worst_imag = im;
    }
    { STEP("[D] Hermitian 2x2: max|Im λ| < 1e-15");
      EXPECT_LT(worst_imag, 1e-15) << Fmt(worst_imag); }

    const double a = Re(e.Lamda(0, 0)), b = Re(e.Lamda(1, 1));
    const double hi = (a > b) ? a : b, lo = (a > b) ? b : a;
    { STEP("[D] [[2,1-i],[1+i,3]]: λ = 4 / 1");
      EXPECT_LT(std::abs(hi - 4.0), 1e-12) << "λmax = " << hi;
      EXPECT_LT(std::abs(lo - 1.0), 1e-12) << "λmin = " << lo; }
}

// ---- [E] 边界与异常 ----
TEST(Eigen, BoundaryAndContract)
{
    { STEP("[E] 1x1: λ = 7");
      Matrix<D> A1 = {{7.0}};
      auto e1 = eig(A1);
      EXPECT_LT(std::abs(Re(e1.Lamda(0, 0)) - 7.0), 1e-15); }

    Matrix<D> Z = Matrix<D>::Zero(3, 3);
    auto ez = eig(Z);
    { STEP("[E] 零矩阵: λ 全 0");
      EXPECT_LT(double(NormInf(ez.Lamda)), 1e-15); }
    { STEP("[E] 零矩阵: Q = I");
      EXPECT_TRUE(IsApprox(ez.Q, Matrix<D>::Identity(3), 1e-15)); }

    { STEP("[E] 对角 [[3,-1,2]]: 重构残差 < 1e-14");
      Matrix<D> Dg = Matrix<D>::Diagonal({3.0, -1.0, 2.0});
      auto ed = eig(Dg);
      EXPECT_LT(ReconResid(Dg, ed), 1e-14) << Fmt(ReconResid(Dg, ed)); }

    Matrix<D> NS = {{1.0, 2.0}, {0.0, 1.0}};   // 上三角非对称
    { STEP("[E] 非对称 -> throw std::logic_error (阶段 2 未实现)");
      EXPECT_EXACT_THROW((void)eig(NS), std::logic_error); }

    Matrix<D> Rect(2, 3);
    { STEP("[E] 非方阵 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)eig(Rect), std::invalid_argument); }

    Matrix<D> Emp(0, 0);
    { STEP("[E] 0x0 -> throw std::invalid_argument");
      EXPECT_EXACT_THROW((void)eig(Emp), std::invalid_argument); }
}