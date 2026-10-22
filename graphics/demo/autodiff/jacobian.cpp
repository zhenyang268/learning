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
// 判据纪律 (graphics 3.3)
// ---------------------------------------------------------------------------
//   · [B] 解析值: 与被测实现完全独立 (手写闭式, 走 std::*)
//   · [C] 中心差分: **只读 value 不碰 derive**, 是独立于 Dual 的第二条路径
//   · [E] 图反向: 与 Jacobian() 的 sweep 实现互不相干
//   · 每个数值判据都配了 **double 与 float 两域** (单域探针连续误判过 5 次)
//   · 浮点判据一律判区间, 不写 == 0; [F] 的 n=1 退化专门覆盖病态输入
//
// 运行: ./graphics/run_tests.sh             (注册名: jacobian)
//       ./graphics/build/demo/autodiff/jacobian

#include "../common.hpp"
#include "autodiff/dual.hpp"
#include "autodiff/grad_graph.hpp"
#include "autodiff/jacobian.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

static Checker ck;

// n×1 **列**向量 (列向量约定)
template <typename T>
static Matrix<Dual<T>> ColD(const std::vector<T>& vs)
{
    std::vector<Dual<T>> d;
    d.reserve(vs.size());
    for (T v : vs) d.emplace_back(v, T(0));
    return Matrix<Dual<T>>(static_cast<int>(vs.size()), 1, d);
}

template <typename T>
static void PrintM(const char* tag, const Matrix<T>& M)
{
    std::printf("     %-10s %dx%d\n", tag, M.Rows(), M.Cols());
    for (int i = 0; i < M.Rows(); i++) {
        std::printf("               [");
        for (int j = 0; j < M.Cols(); j++) std::printf("%10.6f ", static_cast<double>(M(i, j)));
        std::printf("]\n");
    }
}

// ===========================================================================
// [A][B][C] 形状 + 解析 + 中心差分
//   f: R^3 -> R^2,  y = [x0*x1, sin(x2)]
//   解析 J = [[x1, x0, 0], [0, 0, cos(x2)]]   (2x3)
// ===========================================================================
template <typename T>
static void RunAnalytic(Checker& ck, const char* tag, double tol)
{
    using DQ = Dual<T>;
    auto f = [](const Matrix<DQ>& x) {
        Matrix<DQ> y(2, 1);
        y.Data()[0] = x.Data()[0] * x.Data()[1];
        y.Data()[1] = dualFunc(x.Data()[2], FunctionType::Sin);
        return y;
    };

    const T x0 = T(2), x1 = T(3), x2 = T(0.5);
    const Matrix<DQ> x = ColD<T>({x0, x1, x2});
    const Matrix<T> J = Jacobian(f, x);

    ck.check(J.Rows() == 2 && J.Cols() == 3,
             std::string("[A] ") + tag + " J 形状 = 2x3 (m x n, 行=输出 列=输入)  实得 " +
                 std::to_string(J.Rows()) + "x" + std::to_string(J.Cols()));

    const double c05 = std::cos(0.5);
    bool ok = std::fabs(static_cast<double>(J(0, 0)) - 3.0) < tol
           && std::fabs(static_cast<double>(J(0, 1)) - 2.0) < tol
           && std::fabs(static_cast<double>(J(0, 2)))        < tol
           && std::fabs(static_cast<double>(J(1, 0)))        < tol
           && std::fabs(static_cast<double>(J(1, 1)))        < tol
           && std::fabs(static_cast<double>(J(1, 2)) - c05)  < tol;
    PrintM("J", J);
    ck.check(ok, std::string("[B] ") + tag + " 逐元对拍解析 J = [[x1,x0,0],[0,0,cos x2]]");
}

template <typename T>
static void RunCentralDiff(Checker& ck, const char* tag, double tol, double h)
{
    using DQ = Dual<T>;
    auto f = [](const Matrix<DQ>& x) {
        Matrix<DQ> y(2, 1);
        y.Data()[0] = x.Data()[0] * x.Data()[1];
        y.Data()[1] = dualFunc(x.Data()[2], FunctionType::Sin);
        return y;
    };
    // 只读 value: 与 derive 通道完全无关的第二条路径
    auto val = [&](const std::vector<double>& xv) {
        Matrix<DQ> m = ColD<T>({T(xv[0]), T(xv[1]), T(xv[2])});
        const Matrix<DQ> y = f(m);
        return std::vector<double>{static_cast<double>(y.Data()[0].value),
                                   static_cast<double>(y.Data()[1].value)};
    };

    const std::vector<double> xv{2.0, 3.0, 0.5};
    const Matrix<DQ> x = ColD<T>({T(xv[0]), T(xv[1]), T(xv[2])});
    const Matrix<T> J = Jacobian(f, x);

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
    std::printf("     最差 |J - 中心差分| = %.3e\n", worst);
    ck.check(worst < tol, std::string("[C] ") + tag + " 6 元全部对拍中心差分");
}

// ===========================================================================
// [D] 链式法则自然序: J_{f∘g} = J_f · J_g
//   这一条是分子布局与分母布局的**判别式**:
//   分母布局下同一条链要写成 J_g · J_f, 且形状 3x2 · 2x3 = 3x3, 根本乘不起来。
// ===========================================================================
template <typename T>
static void RunChainRule(Checker& ck, const char* tag, double tol)
{
    using DQ = Dual<T>;
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
    std::printf("     最差 |Jf*Jg - J(fog)| = %.3e   (自然序; 分母布局要写 Jg*Jf)\n", worst);
    ck.check(worst < tol, std::string("[D] ") + tag + " J_{f∘g} = J_f · J_g (自然序)");
}

// ===========================================================================
// [E] 与图反向模式对拍
//   y = W·x (W 2x3, x 3x1)。给第 i 个输出播种 e_i 后, x.adjoint 就是 J 的第 i 行。
//   这是与 Jacobian() 的"一次 sweep 填一列"完全独立的一条路径。
// ===========================================================================
template <typename T>
static void RunAgainstGraph(Checker& ck, const char* tag, double tol)
{
    using DQ = Dual<T>;
    using Node = GradNode<DQ>;

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
    std::printf("     最差 |反向 - 前向| = %.3e\n", worst);
    ck.check(worst < tol, std::string("[E] ") + tag + " 图反向 vs Jacobian() 逐元一致");
}

// ===========================================================================
// [F] 退化情形: m=1（标量输出 -> J 是 1×n 行, 梯度 = J^T）与 n=1
// ===========================================================================
template <typename T>
static void RunDegenerate(Checker& ck, const char* tag, double tol)
{
    using DQ = Dual<T>;

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
        ck.check(shp && val,
                 std::string("[F1] ") + tag + " 标量输出: J 是 1x3 **行** (= 梯度的转置)");
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
        ck.check(shp && val, std::string("[F2] ") + tag + " 单输入: J 是 2x1 列");
    }
}

// ===========================================================================
// [G] 非法输入 -> throw (契约写进签名不够, 要有判据钉住)
// ===========================================================================
static void TestInvalidInputs()
{
    std::printf("\n[G] 非法输入 -> throw\n");
    using DQ = Dual<double>;

    // x 是行向量 -> 拒绝
    {
        bool threw = false;
        try {
            Matrix<DQ> row(1, 3);                       // 1x3 行向量
            auto f = [](const Matrix<DQ>& v) { return v; };
            (void)Jacobian(f, row);
        } catch (const std::invalid_argument&) { threw = true; }
        ck.check(threw, std::string("[G1] x 是 1x3 行向量 -> invalid_argument"));
    }
    // 输出是行向量 -> 拒绝
    {
        bool threw = false;
        try {
            Matrix<DQ> x = ColD<double>({1.0, 2.0});
            auto f = [](const Matrix<DQ>& v) {
                Matrix<DQ> y(1, 2);                     // 1x2 行向量
                y.Data()[0] = v.Data()[0];
                y.Data()[1] = v.Data()[1];
                return y;
            };
            (void)Jacobian(f, x);
        } catch (const std::invalid_argument&) { threw = true; }
        ck.check(threw, std::string("[G2] 输出是 1x2 行向量 -> invalid_argument"));
    }
    // 声明 m 与实际不符 -> 拒绝
    {
        bool threw = false;
        try {
            Matrix<DQ> x = ColD<double>({1.0, 2.0});
            auto f = [](const Matrix<DQ>& v) {
                Matrix<DQ> y(2, 1);
                y.Data()[0] = v.Data()[0];
                y.Data()[1] = v.Data()[1];
                return y;
            };
            (void)Jacobian(f, x, 3);                    // 谎报 m=3
        } catch (const std::invalid_argument&) { threw = true; }
        ck.check(threw, std::string("[G3] 声明 m 与实际输出不符 -> invalid_argument"));
    }
    // Apply 的函数矩阵形状不符 -> 拒绝 (2026-10-19 修的那个越界读)
    {
        bool threw = false;
        try {
            Matrix<DQ> A(3, 1);
            Matrix<FunctionType> F(2, 1, std::vector<FunctionType>{
                FunctionType::Sin, FunctionType::Cos});
            (void)Apply(A, F);
        } catch (const std::invalid_argument&) { threw = true; }
        ck.check(threw,
                 std::string("[G4] Apply: A 3x1 但函数矩阵 2x1 -> invalid_argument (旧版漏检)"));
    }
}

// ===========================================================================

int main()
{
    std::printf("== 分子布局雅可比 Jacobian() (D19 第二次修订) ==\n");
    std::printf("   布局: J_ij = d y_i / d x_j, 形状 m x n (行 = 输出, 列 = 输入)\n");

    std::printf("\n-- double 域 --\n");
    RunAnalytic<double>(ck, "double", 1e-15);
    RunCentralDiff<double>(ck, "double", 1e-7, 1e-6);
    RunChainRule<double>(ck, "double", 1e-14);
    RunAgainstGraph<double>(ck, "double", 1e-14);
    RunDegenerate<double>(ck, "double", 1e-14);

    std::printf("\n-- float 域 (分域纪律: 同一模板换域再实例化一次) --\n");
    RunAnalytic<float>(ck, "float ", 1e-6);
    RunCentralDiff<float>(ck, "float ", 1e-3, 1e-3);
    RunChainRule<float>(ck, "float ", 1e-5);
    RunAgainstGraph<float>(ck, "float ", 1e-5);
    RunDegenerate<float>(ck, "float ", 1e-5);

    TestInvalidInputs();

    return CheckReport(ck.ok, ck.fail, "autodiff/jacobian");
}
