// compute_graph.cpp -- 计算图前向模式基线 (列向量约定, 2026-10-19)
//
// ---------------------------------------------------------------------------
// 布局 (D19 第二次修订)
// ---------------------------------------------------------------------------
// 向量一律 **n×1 列向量**, 雅可比取**分子布局** J_ij = ∂y_i/∂x_j (m×n)。
// 本文件只跑前向 (Dual seed), 不涉及雅可比组装 —— 那是 `jacobian.cpp` 的职责。
//
// 本文件 2026-10-19 从行向量 (1×n) 改为列向量 (n×1), 同时把 float 统一成 double:
// 旧版因 float 的舍入 (~eps*|y|/h) 被迫把差分步长放大到 h=1e-2 才能不假报 FAIL;
// double 下 h=1e-6 的舍入只有 ~3e-9, 判据可以收紧到 1e-8。
//
// 四个用例 (语义未变, 只是索引从 [0][k] 改成 [k][0]):
//   [1] 单层线性   y = w*x + b
//   [2] 逐元素乘   y = x1 ⊙ x2
//   [3] 多层嵌套   y = w2*(w1*x + b1) + b2
//   [4] 损失求和   y = sum(x1 ⊙ x2)
//
// 运行: ./graphics/run_tests.sh            (注册名: compute_graph)
//       ./graphics/build/demo/autodiff/compute_graph

#include "autodiff/grad_graph.hpp"
#include "core/matrix.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using D       = Dual<double>;
using Node    = GradNode<D>;
using NodePtr = std::shared_ptr<Node>;
using Graph   = ComputationGraph<D>;

// ---------------------------------------------------------------------------
// 节点构造辅助
// ---------------------------------------------------------------------------

// 1x1 标量叶子节点 (Input / Constant)
NodePtr leaf(GradOp op, double v)
{
    Matrix<D> m(1, 1, std::vector<D>{D(v)});
    return std::make_shared<Node>(op, std::move(m));
}

// Nx1 **列**向量叶子节点
NodePtr leaf_vec(GradOp op, const std::vector<double>& vals)
{
    std::vector<D> d;
    d.reserve(vals.size());
    for (double v : vals) d.emplace_back(v);
    Matrix<D> m(static_cast<int>(vals.size()), 1, d);
    return std::make_shared<Node>(op, std::move(m));
}

// 二元运算节点, 形状与第一个输入一致 (Nx1)
NodePtr binop(GradOp op, const NodePtr& a, const NodePtr& b)
{
    const int n = a->matrix.Rows();
    Matrix<D> m(n, 1, std::vector<D>(n));
    auto node = std::make_shared<Node>(op, std::move(m));
    node->FillInput(a);
    node->FillInput(b);
    return node;
}

// Output 透传节点, 形状与输入一致 (Nx1)
NodePtr out(const NodePtr& src)
{
    const int n = src->matrix.Rows();
    Matrix<D> m(n, 1, std::vector<D>(n));
    auto node = std::make_shared<Node>(GradOp::Output, std::move(m));
    node->FillInput(src);
    return node;
}

// ---------------------------------------------------------------------------
// 测试辅助
// ---------------------------------------------------------------------------

static int failures = 0;

void check(const char* name, double got, double want, double eps = 1e-8)
{
    const bool ok = std::fabs(got - want) <= eps;
    if (!ok) failures++;
    std::printf("  %-22s got=%10.5f  want=%10.5f  %s\n",
                name, got, want, ok ? "[OK]" : "[FAIL]");
}

// 中心差分: 扰动 leaf 的第 leaf_row 个分量, 读取输出第 out_row 个分量
// h=1e-6: 本文件四个用例对每个被扰动变量都是线性/双线性 ⇒ 截断误差恒为 0,
// 只剩舍入 ~ eps*|y|/h ≈ 1e-16*30/1e-6 = 3e-9, 落在 eps=1e-8 之下。
double numeric_grad(Graph& g, NodePtr leaf, int leaf_row = 0, int out_row = 0,
                    double h = 1e-6)
{
    const double x0 = leaf->matrix[leaf_row][0].value;
    leaf->matrix[leaf_row][0].value = x0 + h;  g.Forward();
    const double yp = g.output->matrix[out_row][0].value;
    leaf->matrix[leaf_row][0].value = x0 - h;  g.Forward();
    const double ym = g.output->matrix[out_row][0].value;
    leaf->matrix[leaf_row][0].value = x0;      g.Forward();   // 恢复现场
    return (yp - ym) / (2.0 * h);
}

// ---------------------------------------------------------------------------
// 用例 1: 单层线性 y = w*x + b
// ---------------------------------------------------------------------------
void test_linear()
{
    std::printf("\n[1] 单层线性  y = w*x + b   (x=3, w=2, b=1)\n");

    auto x    = leaf(GradOp::Input, 3.0);
    auto w    = leaf(GradOp::Input, 2.0);
    auto b    = leaf(GradOp::Input, 1.0);
    auto mul  = binop(GradOp::Mult, x, w);    // a = x*w
    auto add  = binop(GradOp::Add,  mul, b);  // y = a + b
    auto outn = out(add);
    Graph g(x, outn);

    x->matrix[0][0].derive = 0;
    w->matrix[0][0].derive = 0;
    b->matrix[0][0].derive = 0;
    g.Forward();
    check("y 前向值 (=7)", outn->matrix[0][0].value, 7.0);

    // dy/dx = w = 2
    x->matrix[0][0].derive = 1; w->matrix[0][0].derive = 0; b->matrix[0][0].derive = 0;
    g.Forward();
    check("dy/dx (=2)", outn->matrix[0][0].derive, 2.0);
    check("dy/dx 数值", numeric_grad(g, x),        2.0);

    // dy/dw = x = 3
    x->matrix[0][0].derive = 0; w->matrix[0][0].derive = 1; b->matrix[0][0].derive = 0;
    g.Forward();
    check("dy/dw (=3)", outn->matrix[0][0].derive, 3.0);
    check("dy/dw 数值", numeric_grad(g, w),        3.0);

    // dy/db = 1
    x->matrix[0][0].derive = 0; w->matrix[0][0].derive = 0; b->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db (=1)", outn->matrix[0][0].derive, 1.0);
    check("dy/db 数值", numeric_grad(g, b),        1.0);
}

// ---------------------------------------------------------------------------
// 用例 2: 逐元素乘法(列向量) y = x1 ⊙ x2
// ---------------------------------------------------------------------------
void test_mul()
{
    std::printf("\n[2] 逐元素乘法(列向量)  y = x1 ⊙ x2  (x1=(2,3), x2=(5,7))\n");

    auto x1   = leaf_vec(GradOp::Input, {2.0, 3.0});
    auto x2   = leaf_vec(GradOp::Input, {5.0, 7.0});
    auto mul  = binop(GradOp::Mult, x1, x2);   // 2x1
    auto outn = out(mul);
    Graph g(x1, outn);

    auto zero = [&] {
        for (int i = 0; i < 2; i++) {
            x1->matrix[i][0].derive = 0;
            x2->matrix[i][0].derive = 0;
        }
    };

    // 前向值: y = (10, 21)
    zero();
    g.Forward();
    check("y[0] (=10)", outn->matrix[0][0].value, 10.0);
    check("y[1] (=21)", outn->matrix[1][0].value, 21.0);

    // ∂y0/∂x1[0] = x2[0] = 5
    zero(); x1->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y0/∂x1[0] (=5)", outn->matrix[0][0].derive, 5.0);
    check("∂y0/∂x1[0] 数值", numeric_grad(g, x1, 0, 0), 5.0);

    // ∂y1/∂x1[1] = x2[1] = 7
    zero(); x1->matrix[1][0].derive = 1;
    g.Forward();
    check("∂y1/∂x1[1] (=7)", outn->matrix[1][0].derive, 7.0);
    check("∂y1/∂x1[1] 数值", numeric_grad(g, x1, 1, 1), 7.0);

    // ∂y0/∂x2[0] = x1[0] = 2
    zero(); x2->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y0/∂x2[0] (=2)", outn->matrix[0][0].derive, 2.0);
    check("∂y0/∂x2[0] 数值", numeric_grad(g, x2, 0, 0), 2.0);

    // ∂y1/∂x2[1] = x1[1] = 3
    zero(); x2->matrix[1][0].derive = 1;
    g.Forward();
    check("∂y1/∂x2[1] (=3)", outn->matrix[1][0].derive, 3.0);
    check("∂y1/∂x2[1] 数值", numeric_grad(g, x2, 1, 1), 3.0);
}

// ---------------------------------------------------------------------------
// 用例 3: 多层嵌套 y = w2*(w1*x + b1) + b2
// ---------------------------------------------------------------------------
void test_mlp()
{
    std::printf("\n[3] 多层  y = w2*(w1*x + b1) + b2\n");
    std::printf("       (x=2, w1=3, b1=1, w2=4, b2=0.5)\n");

    auto x  = leaf(GradOp::Input, 2.0);
    auto w1 = leaf(GradOp::Input, 3.0);
    auto b1 = leaf(GradOp::Input, 1.0);
    auto w2 = leaf(GradOp::Input, 4.0);
    auto b2 = leaf(GradOp::Input, 0.5);

    auto m1   = binop(GradOp::Mult, x,  w1);   // x*w1
    auto h    = binop(GradOp::Add,  m1, b1);   // w1*x + b1
    auto m2   = binop(GradOp::Mult, h,  w2);   // h*w2
    auto y    = binop(GradOp::Add,  m2, b2);   // + b2
    auto outn = out(y);
    Graph g(x, outn);

    auto zero = [&] {
        for (NodePtr p : {x, w1, b1, w2, b2}) p->matrix[0][0].derive = 0;
    };

    // 前向: h = 3*2+1 = 7, y = 4*7+0.5 = 28.5
    zero();
    g.Forward();
    check("h (=7)",    h->matrix[0][0].value,    7.0);
    check("y (=28.5)", outn->matrix[0][0].value, 28.5);

    // dy/dx  = w1*w2 = 12
    zero(); x->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dx (=12)", outn->matrix[0][0].derive, 12.0);
    check("dy/dx 数值",  numeric_grad(g, x),        12.0);

    // dy/dw1 = x*w2 = 8
    zero(); w1->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dw1 (=8)", outn->matrix[0][0].derive, 8.0);
    check("dy/dw1 数值", numeric_grad(g, w1),       8.0);

    // dy/db1 = w2 = 4
    zero(); b1->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db1 (=4)", outn->matrix[0][0].derive, 4.0);
    check("dy/db1 数值", numeric_grad(g, b1),       4.0);

    // dy/dw2 = h = 7
    zero(); w2->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dw2 (=7)", outn->matrix[0][0].derive, 7.0);
    check("dy/dw2 数值", numeric_grad(g, w2),       7.0);

    // dy/db2 = 1
    zero(); b2->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db2 (=1)", outn->matrix[0][0].derive, 1.0);
    check("dy/db2 数值", numeric_grad(g, b2),       1.0);
}

// ---------------------------------------------------------------------------
// 用例 4(附加): 损失求和 y = sum(x1 ⊙ x2)  —— 验证 SumAll / SumLine 路径
// ---------------------------------------------------------------------------
void test_sumall()
{
    std::printf("\n[4] 损失求和  y = sum(x1 ⊙ x2)  (x1=(2,3), x2=(5,7))\n");

    auto x1   = leaf_vec(GradOp::Input, {2.0, 3.0});
    auto x2   = leaf_vec(GradOp::Input, {5.0, 7.0});
    auto mul  = binop(GradOp::Mult, x1, x2);          // (10, 21)
    // 约简节点输出恒为 1x1: 建图期就得给对形状。原先给的是 1x2(照抄 binop 的写法),
    // 于是下游 out() 按 1x2 建节点, 直到 forward 才被改成 1x1 —— 现在靠 Output 直接
    // 覆盖 matrix 没露馅, 但换成 y = sum(a)*2 这种就会在 cwiseProduct 抛尺寸不匹配。
    auto sum = std::make_shared<Node>(GradOp::SumAll,
                                      Matrix<D>(1, 1, std::vector<D>(1)));
    sum->FillInput(mul);
    auto outn = out(sum);
    Graph g(x1, outn);

    auto zero = [&] {
        for (int i = 0; i < 2; i++) {
            x1->matrix[i][0].derive = 0;
            x2->matrix[i][0].derive = 0;
        }
    };

    // 前向: sum = 10 + 21 = 31
    zero();
    g.Forward();
    check("sum 前向值 (=31)", outn->matrix[0][0].value, 31.0);

    // ∂y/∂x1[0] = x2[0] = 5
    zero(); x1->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y/∂x1[0] (=5)", outn->matrix[0][0].derive, 5.0);
    check("∂y/∂x1[0] 数值", numeric_grad(g, x1, 0, 0), 5.0);

    // ∂y/∂x1[1] = x2[1] = 7
    zero(); x1->matrix[1][0].derive = 1;
    g.Forward();
    check("∂y/∂x1[1] (=7)", outn->matrix[0][0].derive, 7.0);
    check("∂y/∂x1[1] 数值", numeric_grad(g, x1, 1, 0), 7.0);

    // ∂y/∂x2[0] = x1[0] = 2
    zero(); x2->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y/∂x2[0] (=2)", outn->matrix[0][0].derive, 2.0);
    check("∂y/∂x2[0] 数值", numeric_grad(g, x2, 0, 0), 2.0);
}

// ---------------------------------------------------------------------------

int main()
{
    test_linear();
    test_mul();
    test_mlp();
    test_sumall();

    std::printf("\n%s (%d failures)\n",
                failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED", failures);
    return failures ? 1 : 0;
}
