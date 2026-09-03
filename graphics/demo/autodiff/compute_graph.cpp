// compute_graph.cpp — 计算图 forward-mode 自动微分测试
//
// 用例(由易到难):
//   1. 单层线性   y = w*x + b
//   2. 逐元素乘法 y = x1 ⊙ x2        (向量)
//   3. 多层嵌套   y = w2*(w1*x + b1) + b2
//   4. 损失求和   y = sum(x1 ⊙ x2)   (SumAll 累加节点)
//
// 每个用例同时用「解析梯度」和「有限差分(数值梯度)」交叉验证。
//
// 编译:  g++ -std=c++17 -I ../../src compute_graph.cpp -o compute_graph
//
// forward-mode 求梯度的约定:
//   Input 节点是叶子, Forward() 不会覆盖其 derive 字段, 由测试预先设置 seed:
//   目标输入置 1, 其余输入置 0, 运行 Forward 后输出节点的 derive 即 ∂y/∂目标输入。

#include "core/matrix.hpp"
#include "autodiff/dual.hpp"
#include "autodiff/grad_graph.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using D       = Dual<float>;
using Node    = GradNode<D>;
using NodePtr = std::shared_ptr<Node>;
using Graph   = ComputationGraph<D>;

// ---------------------------------------------------------------------------
// 节点构造辅助
// ---------------------------------------------------------------------------

// 1x1 叶子节点 (Input / Constant)
NodePtr leaf(GradOp op, float v)
{
    Matrix<D> m(1, 1, std::vector<D>{D(v)});
    return std::make_shared<Node>(op, std::move(m));
}

// 1xN 叶子节点
NodePtr leaf_row(GradOp op, const std::vector<float>& vals)
{
    std::vector<D> d;
    d.reserve(vals.size());
    for (float v : vals) d.emplace_back(v);
    Matrix<D> m(1, (int)vals.size(), d);
    return std::make_shared<Node>(op, std::move(m));
}

// 二元运算节点, 形状与第一个输入一致 (1xN)
NodePtr binop(GradOp op, const NodePtr& a, const NodePtr& b)
{
    int n = a->matrix.Cols();
    Matrix<D> m(1, n, std::vector<D>(n));
    auto node = std::make_shared<Node>(op, std::move(m));
    node->FillInput(a);
    node->FillInput(b);
    return node;
}

// Output 透传节点, 形状与输入一致 (1xN)
NodePtr out(const NodePtr& src)
{
    int n = src->matrix.Cols();
    Matrix<D> m(1, n, std::vector<D>(n));
    auto node = std::make_shared<Node>(GradOp::Output, std::move(m));
    node->FillInput(src);
    return node;
}

// ---------------------------------------------------------------------------
// 测试辅助
// ---------------------------------------------------------------------------

static int failures = 0;

void check(const char* name, float got, float want, float eps = 1e-3f)
{
    bool ok = std::fabs(got - want) <= eps;
    if (!ok) failures++;
    std::printf("  %-22s got=%10.5f  want=%10.5f  %s\n",
                name, got, want, ok ? "[OK]" : "[FAIL]");
}

// 有限差分: 扰动 leaf 的第 leaf_col 列, 读取输出第 out_col 列, 得到数值梯度
float numeric_grad(Graph& g, NodePtr leaf, int leaf_col = 0, int out_col = 0,
                   float h = 1e-3f)
{
    float x0 = leaf->matrix[0][leaf_col].value;
    leaf->matrix[0][leaf_col].value = x0 + h;  g.Forward();
    float yp = g.output->matrix[0][out_col].value;
    leaf->matrix[0][leaf_col].value = x0 - h;  g.Forward();
    float ym = g.output->matrix[0][out_col].value;
    leaf->matrix[0][leaf_col].value = x0;      g.Forward(); // 恢复现场
    return (yp - ym) / (2.0f * h);
}

// ---------------------------------------------------------------------------
// 用例 1: 单层线性 y = w*x + b
// ---------------------------------------------------------------------------
void test_linear()
{
    std::printf("\n[1] 单层线性  y = w*x + b   (x=3, w=2, b=1)\n");

    auto x   = leaf(GradOp::Input, 3.0f);
    auto w   = leaf(GradOp::Input, 2.0f);
    auto b   = leaf(GradOp::Input, 1.0f);
    auto mul = binop(GradOp::Mult, x, w);   // a = x*w
    auto add = binop(GradOp::Add,  mul, b); // y = a + b
    auto outn = out(add);
    Graph g(x, outn);

    // 前向值
    x->matrix[0][0].derive = 0;
    w->matrix[0][0].derive = 0;
    b->matrix[0][0].derive = 0;
    g.Forward();
    check("y 前向值 (=7)", outn->matrix[0][0].value, 7.0f);

    // dy/dx = w = 2
    x->matrix[0][0].derive = 1; w->matrix[0][0].derive = 0; b->matrix[0][0].derive = 0;
    g.Forward();
    check("dy/dx (=2)",   outn->matrix[0][0].derive, 2.0f);
    check("dy/dx 数值",   numeric_grad(g, x),         2.0f);

    // dy/dw = x = 3
    x->matrix[0][0].derive = 0; w->matrix[0][0].derive = 1; b->matrix[0][0].derive = 0;
    g.Forward();
    check("dy/dw (=3)",   outn->matrix[0][0].derive, 3.0f);
    check("dy/dw 数值",   numeric_grad(g, w),         3.0f);

    // dy/db = 1
    x->matrix[0][0].derive = 0; w->matrix[0][0].derive = 0; b->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db (=1)",   outn->matrix[0][0].derive, 1.0f);
    check("dy/db 数值",   numeric_grad(g, b),         1.0f);
}

// ---------------------------------------------------------------------------
// 用例 2: 逐元素乘法(向量) y = x1 ⊙ x2
// ---------------------------------------------------------------------------
void test_mul()
{
    std::printf("\n[2] 逐元素乘法(向量)  y = x1 ⊙ x2  (x1=(2,3), x2=(5,7))\n");

    auto x1 = leaf_row(GradOp::Input, {2.0f, 3.0f});
    auto x2 = leaf_row(GradOp::Input, {5.0f, 7.0f});
    auto mul = binop(GradOp::Mult, x1, x2);   // 1x2
    auto outn = out(mul);
    Graph g(x1, outn);

    auto zero = [&] {
        for (int j = 0; j < 2; j++) { x1->matrix[0][j].derive = 0; x2->matrix[0][j].derive = 0; }
    };

    // 前向值: y = (10, 21)
    zero();
    g.Forward();
    check("y[0] (=10)", outn->matrix[0][0].value, 10.0f);
    check("y[1] (=21)", outn->matrix[0][1].value, 21.0f);

    // ∂y1/∂x1[0] = x2[0] = 5
    zero(); x1->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y1/∂x1[0] (=5)", outn->matrix[0][0].derive, 5.0f);
    check("∂y1/∂x1[0] 数值", numeric_grad(g, x1, 0, 0), 5.0f);

    // ∂y2/∂x1[1] = x2[1] = 7
    zero(); x1->matrix[0][1].derive = 1;
    g.Forward();
    check("∂y2/∂x1[1] (=7)", outn->matrix[0][1].derive, 7.0f);
    check("∂y2/∂x1[1] 数值", numeric_grad(g, x1, 1, 1), 7.0f);

    // ∂y1/∂x2[0] = x1[0] = 2
    zero(); x2->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y1/∂x2[0] (=2)", outn->matrix[0][0].derive, 2.0f);
    check("∂y1/∂x2[0] 数值", numeric_grad(g, x2, 0, 0), 2.0f);

    // ∂y2/∂x2[1] = x1[1] = 3
    zero(); x2->matrix[0][1].derive = 1;
    g.Forward();
    check("∂y2/∂x2[1] (=3)", outn->matrix[0][1].derive, 3.0f);
    check("∂y2/∂x2[1] 数值", numeric_grad(g, x2, 1, 1), 3.0f);
}

// ---------------------------------------------------------------------------
// 用例 3: 多层嵌套 y = w2*(w1*x + b1) + b2
// ---------------------------------------------------------------------------
void test_mlp()
{
    std::printf("\n[3] 多层  y = w2*(w1*x + b1) + b2\n");
    std::printf("       (x=2, w1=3, b1=1, w2=4, b2=0.5)\n");

    auto x  = leaf(GradOp::Input, 2.0f);
    auto w1 = leaf(GradOp::Input, 3.0f);
    auto b1 = leaf(GradOp::Input, 1.0f);
    auto w2 = leaf(GradOp::Input, 4.0f);
    auto b2 = leaf(GradOp::Input, 0.5f);

    auto m1  = binop(GradOp::Mult, x,  w1);  // x*w1
    auto h   = binop(GradOp::Add,  m1, b1);  // w1*x + b1
    auto m2  = binop(GradOp::Mult, h,  w2);  // h*w2
    auto y   = binop(GradOp::Add,  m2, b2);  // + b2
    auto outn = out(y);
    Graph g(x, outn);

    auto zero = [&] {
        for (NodePtr p : {x, w1, b1, w2, b2}) p->matrix[0][0].derive = 0;
    };

    // 前向: h = 3*2+1 = 7, y = 4*7+0.5 = 28.5
    zero();
    g.Forward();
    check("h (=7)",       h->matrix[0][0].value,     7.0f);
    check("y (=28.5)",    outn->matrix[0][0].value, 28.5f);

    // dy/dx  = w1*w2 = 12
    zero(); x->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dx (=12)",  outn->matrix[0][0].derive, 12.0f);
    check("dy/dx 数值",   numeric_grad(g, x),         12.0f);

    // dy/dw1 = x*w2 = 8
    zero(); w1->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dw1 (=8)",  outn->matrix[0][0].derive, 8.0f);
    check("dy/dw1 数值",  numeric_grad(g, w1),        8.0f);

    // dy/db1 = w2 = 4
    zero(); b1->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db1 (=4)",  outn->matrix[0][0].derive, 4.0f);
    check("dy/db1 数值",  numeric_grad(g, b1),        4.0f);

    // dy/dw2 = h = 7
    zero(); w2->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/dw2 (=7)",  outn->matrix[0][0].derive, 7.0f);
    check("dy/dw2 数值",  numeric_grad(g, w2),        7.0f);

    // dy/db2 = 1
    zero(); b2->matrix[0][0].derive = 1;
    g.Forward();
    check("dy/db2 (=1)",  outn->matrix[0][0].derive, 1.0f);
    check("dy/db2 数值",  numeric_grad(g, b2),        1.0f);
}

// ---------------------------------------------------------------------------
// 用例 4(附加): 损失求和 y = sum(x1 ⊙ x2)  —— 验证 SumAll / SumLine 路径
// ---------------------------------------------------------------------------
void test_sumall()
{
    std::printf("\n[4] 损失求和  y = sum(x1 ⊙ x2)  (x1=(2,3), x2=(5,7))\n");

    auto x1 = leaf_row(GradOp::Input, {2.0f, 3.0f});
    auto x2 = leaf_row(GradOp::Input, {5.0f, 7.0f});
    auto mul = binop(GradOp::Mult, x1, x2);          // (10, 21)
    auto sum = std::make_shared<Node>(GradOp::SumAll,
                                      Matrix<D>(1, 2, std::vector<D>(2)));
    sum->FillInput(mul);
    auto outn = out(sum);
    Graph g(x1, outn);

    auto zero = [&] {
        for (int j = 0; j < 2; j++) { x1->matrix[0][j].derive = 0; x2->matrix[0][j].derive = 0; }
    };

    // 前向: sum = 10 + 21 = 31
    zero();
    g.Forward();
    check("sum 前向值 (=31)", outn->matrix[0][0].value, 31.0f);

    // ∂y/∂x1[0] = x2[0] = 5
    zero(); x1->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y/∂x1[0] (=5)", outn->matrix[0][0].derive, 5.0f);
    check("∂y/∂x1[0] 数值", numeric_grad(g, x1, 0, 0), 5.0f);

    // ∂y/∂x1[1] = x2[1] = 7
    zero(); x1->matrix[0][1].derive = 1;
    g.Forward();
    check("∂y/∂x1[1] (=7)", outn->matrix[0][0].derive, 7.0f);
    check("∂y/∂x1[1] 数值", numeric_grad(g, x1, 1, 0), 7.0f);

    // ∂y/∂x2[0] = x1[0] = 2
    zero(); x2->matrix[0][0].derive = 1;
    g.Forward();
    check("∂y/∂x2[0] (=2)", outn->matrix[0][0].derive, 2.0f);
    check("∂y/∂x2[0] 数值", numeric_grad(g, x2, 0, 0), 2.0f);
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
