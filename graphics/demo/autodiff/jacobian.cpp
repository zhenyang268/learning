// jacobian.cpp -- 分子布局雅可比 Jacobian() 的验收 (D19 第二次修订, 2026-10-19)
//
// ---------------------------------------------------------------------------
// 这个用例在钉什么
// ---------------------------------------------------------------------------
// `jacobian.hpp` 这个文件名存在很久了, 但**里面从来没有 Jacobian() 函数** ——
// 只有 Derive / Apply / SumAll / MultAll 这些零件。所以"分子布局"长期是一句没有
// 代码承载的口号: 谁都可以声称自己是分子布局, 没有任何判据能证伪。
//
// 本用例把口径钉成可执行的三件事:
//   [A] 形状 = m×n (行 = 输出分量, 列 = 输入分量) —— 分子布局的定义本身
//   [D] 链式法则 = **自然序** J_{f∘g} = J_f · J_g  —— 这一条能把它和分母布局区分开
//   [E] 与图反向模式逐元对拍 —— 前向 sweep 与反向传播两条独立实现路径必须一致
//
// ---------------------------------------------------------------------------
// 判据纪律 (graphics 3.3) + 框架 (googletest, 2026-10-30 迁移)
// ---------------------------------------------------------------------------
//   · [B] 解析值: 与被测实现完全独立 (手写闭式, 走 std::*)
//   · [C] 中心差分: **只读 value 不碰 derive**, 是独立于 Dual 的第二条路径
//   · [E] 图反向: 与 Jacobian() 的 sweep 实现互不相干
//   · 每个数值判据都配了 **double 与 float 两域** (单域探针连续误判过 5 次)
//   · 浮点判据一律判区间, 不写 == 0; [F] 的 n=1 退化专门覆盖病态输入
//
// 迁移说明: 老写法 `RunXxx<double>(ck, "double", tol)` + `RunXxx<float>(...)` 各写一遍,
//   现在由 **TYPED_TEST_SUITE(double, float)** 保证两域; **阈值随域不同** (double 1e-15 /
//   float 1e-6 等), 所以塞进 `JacTol<T>` 特化, **一个数都没改**。
//   迁移前 18 项断言 -> 迁移后 18 条 (5 判据 x 2 域 + [G] 4 条)。
//   [G] 与类型无关, 仍是普通 TEST。
//
// ---------------------------------------------------------------------------
// 判据边界表 (格式见 demo/README.md §4.3.1)
// ---------------------------------------------------------------------------
//   段 | 阈值 (double / float)      | 边界 / 病态输入              | 专抓什么错
//  ----|----------------------------|------------------------------|--------------------------
//   A  | 形状比较, 无阈值           | f: R^3 -> R^2                | 分子/分母布局搞反 (3x2 vs 2x3)
//   B  | 1e-15 / 1e-6               | x = (2, 3, 0.5), 含 0 元素   | sweep 填错列 / 漏写 0
//   C  | 1e-7 (h=1e-6) / 1e-3 (h=1e-3)| 同上 6 元全查              | 只有 Dual 一条路算对
//   D  | 1e-14 / 1e-5               | f: R^3->R^2 ∘ g: R^2->R^3    | **分母布局** (J_g·J_f 乘不起来)
//   E  | 1e-14 / 1e-5               | y = W·x, W 2x3, 播种 e_i     | 前向与反向两条路径不一致
//   F  | 1e-14 / 1e-5               | **m=1** (标量输出) 与 **n=1** | J 退化成行/列时形状错
//   G  | **异常类型精确到子类**     | 1x3 行向量 / 1x2 输出 / 谎报 m / Apply 形状不符 | 契约只写在注释里没判据
//
// ---------------------------------------------------------------------------
// 已知盲区 (注入对照的方案写在这里)
// ---------------------------------------------------------------------------
//  ① **[D] 是唯一能把分子布局与分母布局区分开的判据**: 分母布局下要写 J_g·J_f,
//     而 3x2 · 2x3 根本乘不起来 —— 形状就先炸。其余 [A]/[B]/[C]/[E] 在两种布局下
//     都能"凑"出结果 (只是转置关系), 抓不到。
//  ② **[F1] 的 m=1 是唯一钉住"标量输出 ⇒ J 是 1×n 行, 梯度 = J^T"的判据**:
//     n=1 的 [F2] 与它形状不同, 不能互相代替。
//  ③ **[C] 中心差分只读 value** —— 它抓不到"derive 通道本身写错" (那条只有 [B] 抓得到);
//     反过来 [B] 与被测实现共享 dualFunc, 抓不到"sin 导数抄错" ([C] 抓得到)。
//
// 运行: ./graphics/run_tests.sh             (注册名: jacobian)
//       ./graphics/build/demo/autodiff/jacobian
//       ./graphics/build/demo/autodiff/jacobian --gtest_filter='*double*'

#include "gtest_common.hpp"

#include "autodiff/dual.hpp"
#include "autodiff/grad_graph.hpp"
#include "autodiff/jacobian.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

// 域 -> 阈值 (原文件 main() 里两个域各写一组, 数值原样搬过来)
template <typename T>
struct JacTol
{
    static constexpr double analytic = 1e-15;
    static constexpr double diff     = 1e-7;
    static constexpr double diff_h   = 1e-6;
    static constexpr double chain    = 1e-14;
    static constexpr double graph    = 1e-14;
    static constexpr double degen    = 1e-14;
};
template <>
struct JacTol<float>
{
    static constexpr double analytic = 1e-6;
    static constexpr double diff     = 1e-3;
    static constexpr double diff_h   = 1e-3;
    static constexpr double chain    = 1e-5;
    static constexpr double graph    = 1e-5;
    static constexpr double degen    = 1e-5;
};

const char* DomainName(double) { return "double"; }
const char* DomainName(float) { return "float "; }

// n×1 **列**向量 (列向量约定)
template <typename T>
Matrix<Dual<T>> ColD(const std::vector<T>& vs)
{
    std::vector<Dual<T>> d;
    d.reserve(vs.size());
    for (T v : vs) d.emplace_back(v, T(0));
    return Matrix<Dual<T>>(static_cast<int>(vs.size()), 1, d);
}

template <typename T>
void PrintM(const char* tag, const Matrix<T>& M)
{
    std::printf("     %-10s %dx%d\n", tag, M.Rows(), M.Cols());
    for (int i = 0; i < M.Rows(); i++) {
        std::printf("               [");
        for (int j = 0; j < M.Cols(); j++) std::printf("%10.6f ", static_cast<double>(M(i, j)));
        std::printf("]\n");
    }
}

// 被反复使用的测试函数 f: R^3 -> R^2, y = [x0*x1, sin(x2)]
template <typename T>
Matrix<Dual<T>> FAnalytic(const Matrix<Dual<T>>& x)
{
    using DQ = Dual<T>;
    Matrix<DQ> y(2, 1);
    y.Data()[0] = x.Data()[0] * x.Data()[1];
    y.Data()[1] = dualFunc(x.Data()[2], FunctionType::Sin);
    return y;
}

}  // namespace

#define STEP(label) SCOPED_TRACE(label)

template <typename T>
class JacobianTest : public ::testing::Test {};

using JacobianTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(JacobianTest, JacobianTypes);

// ===========================================================================
// [A][B] 形状 + 解析对拍
//   f: R^3 -> R^2,  y = [x0*x1, sin(x2)]
//   解析 J = [[x1, x0, 0], [0, 0, cos(x2)]]   (2x3)
// ===========================================================================
TYPED_TEST(JacobianTest, ShapeAndAnalytic)
{
    using T  = TypeParam;
    using DQ = Dual<T>;
    const char* tag = DomainName(T());
    const double tol = JacTol<T>::analytic;

    const T x0 = T(2), x1 = T(3), x2 = T(0.5);
    const Matrix<DQ> x = ColD<T>({x0, x1, x2});
    const Matrix<T> J = Jacobian(FAnalytic<T>, x);

    { STEP(std::string("[A] ") + tag + " J 形状 = 2x3 (m x n, 行=输出 列=输入)");
      EXPECT_TRUE(J.Rows() == 2 && J.Cols() == 3)
          << "实得 " << J.Rows() << "x" << J.Cols(); }

    const double c05 = std::cos(0.5);
    bool ok = std::fabs(static_cast<double>(J(0, 0)) - 3.0) < tol
           && std::fabs(static_cast<double>(J(0, 1)) - 2.0) < tol
           && std::fabs(static_cast<double>(J(0, 2)))        < tol
           && std::fabs(static_cast<double>(J(1, 0)))        < tol
           && std::fabs(static_cast<double>(J(1, 1)))        < tol
           && std::fabs(static_cast<double>(J(1, 2)) - c05)  < tol;
    PrintM("J", J);
    { STEP(std::string("[B] ") + tag + " 逐元对拍解析 J = [[x1,x0,0],[0,0,cos x2]]");
      EXPECT_TRUE(ok); }
}

// ===========================================================================
// [C] 中心差分: 只读 value, 独立于 derive 通道的第二条路径
// ===========================================================================
TYPED_TEST(JacobianTest, CentralDifference)
{
    using T  = TypeParam;
    using DQ = Dual<T>;
    const char* tag = DomainName(T());
    const double tol = JacTol<T>::diff;
    const double h   = JacTol<T>::diff_h;

    // 只读 value: 与 derive 通道完全无关的第二条路径
    auto val = [&](const std::vector<double>& xv) {
        Matrix<DQ> m = ColD<T>({T(xv[0]), T(xv[1]), T(xv[2])});
        const Matrix<DQ> y = FAnalytic<T>(m);
        return std::vector<double>{static_cast<double>(y.Data()[0].value),
                                   static_cast<double>(y.Data()[1].value)};
    };

    const std::vector<double> xv{2.0, 3.0, 0.5};
    const Matrix<DQ> x = ColD<T>({T(xv[0]), T(xv[1]), T(xv[2])});
    const Matrix<T> J = Jacobian(FAnalytic<T>, x);

    double worst = 0.0;
    for (int j = 0; j < 3; j++) {
        std::vector<double> p = xv, m = xv;
        p[j] += h;
        m[j] -= h;
        const std::vector<double> yp = val(p), ym = val(m);
        for (int i = 0; i < 2; i++) {
            const double num = (yp[i] - ym[i]) / (2.0 * h);
            const double err = std::fabs(static_cast<double>(J(i, j)) - num);
            if (err > worst) worst = err;
        }
    }
    std::printf("[jacobian] %s 最差 |J - 中心差分| = %.3e\n", tag, worst);

    { STEP(std::string("[C] ") + tag + " 6 元全部对拍中心差分");
      EXPECT_LT(worst, tol) << "最差 = " << worst; }
}

// ===========================================================================
// [D] 链式法则自然序: J_{f∘g} = J_f · J_g
//   这一条是分子布局与分母布局的**判别式**:
//   分母布局下同一条链要写成 J_g · J_f, 且形状 3x2 · 2x3 = 3x3, 根本乘不起来。
// ===========================================================================
TYPED_TEST(JacobianTest, ChainRuleNaturalOrder)
{
    using T  = TypeParam;
    using DQ = Dual<T>;
    const char* tag = DomainName(T());
    const double tol = JacTol<T>::chain;

    // g: R^2 -> R^3,  g(x) = [x0+x1, x0*x1, sin(x0)]
    auto g = [](const Matrix<DQ>& x) {
        Matrix<DQ> y(3, 1);
        y.Data()[0] = x.Data()[0] + x.Data()[1];
        y.Data()[1] = x.Data()[0] * x.Data()[1];
        y.Data()[2] = dualFunc(x.Data()[0], FunctionType::Sin);
        return y;
    };
    // f: R^3 -> R^2,  f(y) = [y0*y2, y1+y2]
    auto f = [](const Matrix<DQ>& y) {
        Matrix<DQ> z(2, 1);
        z.Data()[0] = y.Data()[0] * y.Data()[2];
        z.Data()[1] = y.Data()[1] + y.Data()[2];
        return z;
    };
    auto h = [&](const Matrix<DQ>& x) { return f(g(x)); };

    const Matrix<DQ> x = ColD<T>({T(0.7), T(1.3)});
    const Matrix<T> Jf = Jacobian(f, g(x));   // 2x3
    const Matrix<T> Jg = Jacobian(g, x);      // 3x2
    const Matrix<T> Jh = Jacobian(h, x);      // 2x2
    const Matrix<T> prod = Jf * Jg;           // 自然序: 2x3 · 3x2 = 2x2

    double worst = 0.0;
    if (prod.Rows() != Jh.Rows() || prod.Cols() != Jh.Cols()) {
        worst = 1e9;
    } else {
        for (int i = 0; i < prod.Rows(); i++)
            for (int j = 0; j < prod.Cols(); j++) {
                const double err = std::fabs(static_cast<double>(prod(i, j))
                                           - static_cast<double>(Jh(i, j)));
                if (err > worst) worst = err;
            }
    }
    PrintM("Jf (2x3)", Jf);
    PrintM("Jg (3x2)", Jg);
    PrintM("Jf*Jg", prod);
    PrintM("J(fog)", Jh);
    std::printf("[jacobian] %s 最差 |Jf*Jg - J(fog)| = %.3e   (自然序; 分母布局要写 Jg*Jf)\n",
                tag, worst);

    { STEP(std::string("[D] ") + tag + " J_{f∘g} = J_f · J_g (自然序)");
      EXPECT_LT(worst, tol) << "最差 = " << worst; }
}

// ===========================================================================
// [E] 与图反向模式对拍
//   y = W·x (W 2x3, x 3x1)。给第 i 个输出播种 e_i 后, x.adjoint 就是 J 的第 i 行。
//   这是与 Jacobian() 的"一次 sweep 填一列"完全独立的一条路径。
// ===========================================================================
TYPED_TEST(JacobianTest, AgainstGraphBackward)
{
    using T    = TypeParam;
    using DQ   = Dual<T>;
    using Node = GradNode<DQ>;
    const char* tag = DomainName(T());
    const double tol = JacTol<T>::graph;

    const double wv[6] = {1, 2, 3, 4, 5, 6};
    Matrix<DQ> Wm(2, 3);
    for (int i = 0; i < 6; i++) Wm.Data()[i] = DQ(T(wv[i]), T(0));
    Matrix<DQ> xm = ColD<T>({T(1), T(1), T(1)});

    auto W = std::make_shared<Node>(GradOp::Input, Wm);
    auto x = std::make_shared<Node>(GradOp::Input, xm);
    auto y = std::make_shared<Node>(GradOp::MatrixMult, Matrix<DQ>(2, 1));
    y->FillInput(W);
    y->FillInput(x);

    ComputationGraph<DQ> g(W, y);
    g.Forward();

    Matrix<T> Jrev(2, 3);
    for (int i = 0; i < 2; i++) {
        Matrix<DQ> seed(2, 1);
        seed.Data()[i] = DQ(T(1), T(0));      // 种子 e_i
        g.Backward(seed);
        for (int j = 0; j < 3; j++) Jrev(i, j) = x->adjoint.Data()[j].value;
    }

    auto lin = [&](const Matrix<DQ>& xin) { return Wm * xin; };
    const Matrix<T> Jfw = Jacobian(lin, xm);

    double worst = 0.0;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++) {
            const double err = std::fabs(static_cast<double>(Jrev(i, j))
                                       - static_cast<double>(Jfw(i, j)));
            if (err > worst) worst = err;
        }
    PrintM("J (反向)", Jrev);
    PrintM("J (前向)", Jfw);
    std::printf("[jacobian] %s 最差 |反向 - 前向| = %.3e\n", tag, worst);

    { STEP(std::string("[E] ") + tag + " 图反向 vs Jacobian() 逐元一致");
      EXPECT_LT(worst, tol) << "最差 = " << worst; }
}

// ===========================================================================
// [F] 退化情形: m=1（标量输出 -> J 是 1×n 行, 梯度 = J^T）与 n=1
// ===========================================================================
TYPED_TEST(JacobianTest, DegenerateShapes)
{
    using T  = TypeParam;
    using DQ = Dual<T>;
    const char* tag = DomainName(T());
    const double tol = JacTol<T>::degen;

    {   // m = 1: f(x) = x0 + x1 + x2 -> J = [1, 1, 1] (1x3)
        auto f = [](const Matrix<DQ>& x) {
            Matrix<DQ> y(1, 1);
            y.Data()[0] = x.Data()[0] + x.Data()[1] + x.Data()[2];
            return y;
        };
        const Matrix<DQ> x = ColD<T>({T(1), T(2), T(3)});
        const Matrix<T> J = Jacobian(f, x);
        const bool shp = (J.Rows() == 1 && J.Cols() == 3);
        const bool val = std::fabs(static_cast<double>(J(0, 0)) - 1.0) < tol
                      && std::fabs(static_cast<double>(J(0, 1)) - 1.0) < tol
                      && std::fabs(static_cast<double>(J(0, 2)) - 1.0) < tol;
        { STEP(std::string("[F1] ") + tag + " 标量输出: J 是 1x3 **行** (= 梯度的转置)");
          EXPECT_TRUE(shp && val) << "形状 " << J.Rows() << "x" << J.Cols(); }
    }

    {   // n = 1: f(x) = [sin x, x^2] -> J = [[cos x], [2x]] (2x1)
        auto f = [](const Matrix<DQ>& x) {
            Matrix<DQ> y(2, 1);
            y.Data()[0] = dualFunc(x.Data()[0], FunctionType::Sin);
            y.Data()[1] = x.Data()[0] * x.Data()[0];
            return y;
        };
        const Matrix<DQ> x = ColD<T>({T(1.3)});
        const Matrix<T> J = Jacobian(f, x);
        const bool shp = (J.Rows() == 2 && J.Cols() == 1);
        const bool val = std::fabs(static_cast<double>(J(0, 0)) - std::cos(1.3)) < 1e-6
                      && std::fabs(static_cast<double>(J(1, 0)) - 2.6) < 1e-6;
        { STEP(std::string("[F2] ") + tag + " 单输入: J 是 2x1 列");
          EXPECT_TRUE(shp && val) << "形状 " << J.Rows() << "x" << J.Cols(); }
    }
}

// ===========================================================================
// [G] 非法输入 -> throw (契约写进签名不够, 要有判据钉住)
// ===========================================================================
TEST(Jacobian, InvalidInputs)
{
    using DQ = Dual<double>;

    // x 是行向量 -> 拒绝
    { STEP("[G1] x 是 1x3 行向量 -> invalid_argument");
      Matrix<DQ> row(1, 3);                       // 1x3 行向量
      auto f = [](const Matrix<DQ>& v) { return v; };
      EXPECT_EXACT_THROW((void)Jacobian(f, row), std::invalid_argument); }

    // 输出是行向量 -> 拒绝
    { STEP("[G2] 输出是 1x2 行向量 -> invalid_argument");
      Matrix<DQ> x = ColD<double>({1.0, 2.0});
      auto f = [](const Matrix<DQ>& v) {
          Matrix<DQ> y(1, 2);                     // 1x2 行向量
          y.Data()[0] = v.Data()[0];
          y.Data()[1] = v.Data()[1];
          return y;
      };
      EXPECT_EXACT_THROW((void)Jacobian(f, x), std::invalid_argument); }

    // 声明 m 与实际不符 -> 拒绝
    { STEP("[G3] 声明 m 与实际输出不符 -> invalid_argument");
      Matrix<DQ> x = ColD<double>({1.0, 2.0});
      auto f = [](const Matrix<DQ>& v) {
          Matrix<DQ> y(2, 1);
          y.Data()[0] = v.Data()[0];
          y.Data()[1] = v.Data()[1];
          return y;
      };
      EXPECT_EXACT_THROW((void)Jacobian(f, x, 3), std::invalid_argument); }   // 谎报 m=3

    // Apply 的函数矩阵形状不符 -> 拒绝 (2026-10-19 修的那个越界读)
    { STEP("[G4] Apply: A 3x1 但函数矩阵 2x1 -> invalid_argument (旧版漏检)");
      Matrix<DQ> A(3, 1);
      Matrix<FunctionType> F(2, 1, std::vector<FunctionType>{
          FunctionType::Sin, FunctionType::Cos});
      EXPECT_EXACT_THROW((void)Apply(A, F), std::invalid_argument); }
}