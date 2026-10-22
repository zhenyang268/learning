// backward.cpp -- grad_graph 反向模式的最小可用验收 (2026-10-17)
//
// ---------------------------------------------------------------------------
// 为什么会有这个用例
// ---------------------------------------------------------------------------
// 2026-10-17 体检发现: `GradNode::Backward` / `ComputationGraph::Backward` 这套反向模式
// **能编译、能跑，但基本是错的** —— 22 个 GradOp 里只有 5 个的伴随公式正确，其余不是抛异常
// 就是**静默出错**。而 `compute_graph.cpp` 只跑 Forward，所以这些错**一个都测不出来**。
// 本用例把修好之后的行为钉住。
//
// ---------------------------------------------------------------------------
// 原实现里的 bug 清单 (全部实测确认; 修法见 grad_graph.hpp 内联注释)
// ---------------------------------------------------------------------------
//   [架构] ComputationGraph::Backward() 直接 throw 拒绝执行
//          —— 理由写的是"adjoint 从没按形状初始化", 注释是对的, 但既然知道了就该初始化。
//   [架构] GradNode::adjoint 从不按形状初始化 (默认 0x0), 而所有分支都写 `node1->adjoint += ...`
//          ⇒ 第一个节点就抛 Matrix dimension mismatch。
//   [静默] `case Mult` **缺失** —— Forward 里 Mult ≡ CWiseMult (同一个 cwiseProduct),
//          但 Backward 只写了 CWiseMult ⇒ **乘法节点反向静默 no-op** (不抛错、不做事)。
//   [静默] `case Output` **缺失** ⇒ 透传节点不回传伴随, 静默 no-op。
//   [静默] `case Transpose`: `node1->adjoint = node1->matrix.Transpose();`
//          —— 两头都错: 用了**原值** (而非上游伴随), 又用 `=` 覆盖 (而非 `+=` 累加)。
//   [抛异常] `case Divide`: `adjoint * node2->matrix.Transpose()` —— 逐元素除法误写成矩阵乘。
//   [抛异常] `case MultAll`: 1×n 的 adjoint1 与 1×1 的上游伴随做 cwiseProduct, 维度不匹配;
//            且 MultAllBackward 返回的 Dual 把偏导放在 derive 分量, 需先 Derive 提取。
//   [抛异常] `case Function`/`MatrixFunction`: `adjoint * Derive(Apply(...))` 用了矩阵乘
//            (1×n · 1×n 不合法), 逐元素函数应为逐元素乘。
//   [语义错] `Derive(Apply(A, f))` 得到的是 **f'(x)·A.derive** (依赖外部种子), 不是 f'(x)。
//            种子为 0 ⇒ 局部导数全 0。故改用 `LocalDeriv` (直接算 f'(value), 与种子无关)。
//   [静默] `case Sigmoid/Relu/Tanh/Softmax` **case 写了但 body 是空的** —— 比缺 case 更隐蔽,
//          `-Wswitch` 抓不到。**本次未动** (Forward 侧本来就是 throw 占位)。
//   [注释过期] `case Inverse` 的 Forward 仍 throw "待 linalg::Inverse 落地 (阶段 1)",
//          而 linalg::Inverse 早在 V1.2 落地验收了。
//
// 本用例覆盖 [A]~[D]; [E] **如实断言占位算子仍不可用** (不假装它可用 —— 照抄 smoke.cpp 的规矩)。
//
// 运行: ./graphics/run_tests.sh          (注册名: backward)
//       ./graphics/build/demo/autodiff/backward

#include "../common.hpp"
#include "autodiff/dual.hpp"
#include "autodiff/grad_graph.hpp"
#include "autodiff/jacobian.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using D    = Dual<double>;
using Node = GradNode<D>;
using P    = std::shared_ptr<Node>;

static Checker ck;

// Nx1 **列**向量 (D19 第二次修订, 2026-10-19: 全库统一列向量 + 分子布局)
static Matrix<D> MV(const std::vector<double>& vs, double seed = 0.0)
{
    std::vector<D> d;
    for (double v : vs) d.emplace_back(v, seed);
    return Matrix<D>((int)vs.size(), 1, d);
}
// NxM
static Matrix<D> MM(const std::vector<std::vector<double>>& rows)
{
    const int r = (int)rows.size(), c = (int)rows[0].size();
    std::vector<D> d;
    for (auto& row : rows) for (double v : row) d.emplace_back(v, 0.0);
    return Matrix<D>(r, c, d);
}
static P Leaf(const Matrix<D>& m) { return std::make_shared<Node>(GradOp::Input, m); }
static P Op(GradOp op, const std::vector<P>& ins, int r = 1, int c = 1)
{
    auto n = std::make_shared<Node>(op, Matrix<D>(r, c));
    for (auto& x : ins) n->FillInput(x);
    return n;
}
static P OpF(GradOp op, const std::vector<P>& ins,
             const std::vector<FunctionType>& fs, int r, int c)
{
    auto n = std::make_shared<Node>(op, Matrix<FunctionType>(r, c, fs));
    for (auto& x : ins) n->FillInput(x);
    return n;
}
// 读伴随的 value 分量 (约定: adjoint 只用 value)
static double AV(const P& n, int i = 0) { return n->adjoint.Data()[i].value; }

// ===========================================================================
// [A] 最小闭环 y = sum(w ⊙ x) + b —— 三方对拍
// ===========================================================================
static double MakeY(const std::vector<double>& xv)
{
    auto x = Leaf(MV(xv)); auto w = Leaf(MV({5, 7})); auto b = Leaf(MV({0.5}));
    auto m = Op(GradOp::CWiseMult, {x, w}, 2, 1);
    auto s = Op(GradOp::SumAll, {m}, 1, 1);
    auto y = Op(GradOp::Add, {s, b}, 1, 1);
    ComputationGraph<D> g(x, y);
    g.Forward();
    return y->matrix(0, 0).value;
}

static void TestA_MinLoop()
{
    std::printf("\n[A] 最小闭环 y = sum(w ⊙ x) + b   x=[2,3] w=[5,7] b=0.5\n");
    auto x = Leaf(MV({2, 3}));
    auto w = Leaf(MV({5, 7}));
    auto b = Leaf(MV({0.5}));
    auto m = Op(GradOp::CWiseMult, {x, w}, 2, 1);
    auto s = Op(GradOp::SumAll, {m}, 1, 1);
    auto y = Op(GradOp::Add, {s, b}, 1, 1);

    ComputationGraph<D> g(x, y);
    g.Forward();
    ck.check(std::fabs(y->matrix(0, 0).value - 31.5) < 1e-14, std::string("[A1] 前向 y = 31.5"));

    g.Backward();      // 原版此处直接 throw
    ck.check(std::fabs(AV(x, 0) - 5.0) < 1e-13, std::string("[A2] dy/dx0 = w0 = 5"));
    ck.check(std::fabs(AV(x, 1) - 7.0) < 1e-13, std::string("[A3] dy/dx1 = w1 = 7"));
    ck.check(std::fabs(AV(w, 0) - 2.0) < 1e-13, std::string("[A4] dy/dw0 = x0 = 2"));
    ck.check(std::fabs(AV(w, 1) - 3.0) < 1e-13, std::string("[A5] dy/dw1 = x1 = 3"));
    ck.check(std::fabs(AV(b, 0) - 1.0) < 1e-13, std::string("[A6] dy/db  = 1"));

    // 独立判据 1: 中心差分
    {
        const double h = 1e-6;
        ck.check(std::fabs((MakeY({2+h,3}) - MakeY({2-h,3})) / (2*h) - 5.0) < 1e-8,
                 std::string("[A7] dy/dx0 vs 中心差分"));
    }
    // 独立判据 2: 前向模式 (Dual seed) 对拍 —— 两条独立实现路径
    {
        auto x2 = Leaf(MV({2,3})); auto w2 = Leaf(MV({5,7})); auto b2 = Leaf(MV({0.5}));
        auto m2 = Op(GradOp::CWiseMult, {x2, w2}, 2, 1);
        auto s2 = Op(GradOp::SumAll, {m2}, 1, 1);
        auto y2 = Op(GradOp::Add, {s2, b2}, 1, 1);
        ComputationGraph<D> g2(x2, y2);
        x2->matrix(0, 0).derive = 1.0;      // 只给 x0 播种
        g2.Forward();
        ck.check(std::fabs(y2->matrix(0, 0).derive - 5.0) < 1e-13,
                 std::string("[A8] dy/dx0 vs 前向模式"));
    }
}

// ===========================================================================
// [B] 逐算子 Backward: 12 个可用算子
// ===========================================================================
static void TestB_PerOperator()
{
    std::printf("\n[B] 逐算子 Backward (手动按形状初始化 adjoint)\n");

    {   // Add
        auto x1 = Leaf(MV({1,2})), x2 = Leaf(MV({2,3}));
        auto y = Op(GradOp::Add, {x1,x2}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x1->adjoint = MV({0,0}); x2->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x1,0)==1 && AV(x1,1)==1 && AV(x2,0)==1 && AV(x2,1)==1,
                 std::string("[B1] Add  -> [1,1] / [1,1]"));
    }
    {   // Sub
        auto x1 = Leaf(MV({1,2})), x2 = Leaf(MV({2,3}));
        auto y = Op(GradOp::Sub, {x1,x2}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x1->adjoint = MV({0,0}); x2->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x1,0)==1 && AV(x2,0)==-1 && AV(x2,1)==-1,
                 std::string("[B2] Sub  -> [1,1] / [-1,-1]"));
    }
    {   // Mult (原版缺 case, 静默 no-op)
        auto x1 = Leaf(MV({2,3})), x2 = Leaf(MV({5,7}));
        auto y = Op(GradOp::Mult, {x1,x2}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x1->adjoint = MV({0,0}); x2->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x1,0)==5 && AV(x1,1)==7 && AV(x2,0)==2 && AV(x2,1)==3,
                 std::string("[B3] Mult -> [5,7] / [2,3]  (原来静默 no-op)"));
    }
    {   // CWiseMult
        auto x1 = Leaf(MV({2,3})), x2 = Leaf(MV({5,7}));
        auto y = Op(GradOp::CWiseMult, {x1,x2}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x1->adjoint = MV({0,0}); x2->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x1,0)==5 && AV(x2,1)==3,
                 std::string("[B4] CWiseMult -> [5,7] / [2,3]"));
    }
    {   // Divide (原版 node1 用矩阵乘, 抛异常) —— 期望 1/x2 与 -x1/x2^2
        auto x1 = Leaf(MV({6,8})), x2 = Leaf(MV({2,4}));
        auto y = Op(GradOp::Divide, {x1,x2}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x1->adjoint = MV({0,0}); x2->adjoint = MV({0,0});
        y->Backward();
        ck.check(std::fabs(AV(x1,0)-0.5)<1e-15 && std::fabs(AV(x1,1)-0.25)<1e-15 &&
                 std::fabs(AV(x2,0)+1.5)<1e-15 && std::fabs(AV(x2,1)+0.5)<1e-15,
                 std::string("[B5] Divide -> 1/x2=[0.5,0.25] / -x1/x2^2=[-1.5,-0.5]"));
    }
    {   // Transpose (原版用原值 + `=`, 静默错)
        auto x = Leaf(MV({1,2}));                  // 2x1 列
        auto y = Op(GradOp::Transpose, {x}, 1, 2); // 输出 1x2
        y->Forward();
        y->adjoint = MM({{1,1}});                  // 1x2
        x->adjoint = MV({0,0});                    // 2x1
        y->Backward();
        ck.check(x->adjoint.Rows()==2 && x->adjoint.Cols()==1 && AV(x,0)==1 && AV(x,1)==1,
                 std::string("[B6] Transpose -> [1,1] 且形状 (2x1) 正确"));
    }
    {   // SumAll
        auto x = Leaf(MV({2,3}));
        auto y = Op(GradOp::SumAll, {x}, 1, 1); y->Forward();
        y->adjoint = MV({1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x,0)==1 && AV(x,1)==1, std::string("[B7] SumAll -> [1,1]"));
    }
    {   // MultAll (原版维度不匹配, 抛异常) —— 期望 [x1, x0] = [3, 2]
        auto x = Leaf(MV({2,3}));
        auto y = Op(GradOp::MultAll, {x}, 1, 1); y->Forward();
        y->adjoint = MV({1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x,0)==3 && AV(x,1)==2, std::string("[B8] MultAll -> [x1,x0]=[3,2]"));
    }
    {   // MatrixMult (列向量): y = W·x, W 2x3, x 3x1
        //   期望 dW = ā_y·xᵀ (2x3 全 1), dx = Wᵀ·ā_y = [5,7,9]
        //   旧版是 1x2 · 2x1 的点积, 只覆盖了 Row·Col 一种形状; 现改成非方阵, 两侧都查。
        auto W = Leaf(MM({{1,2,3},{4,5,6}}));
        auto x = Leaf(MV({1,1,1}));
        auto y = Op(GradOp::MatrixMult, {W,x}, 2, 1); y->Forward();
        y->adjoint = MV({1,1});
        W->adjoint = MM({{0,0,0},{0,0,0}});
        x->adjoint = MV({0,0,0});
        y->Backward();
        ck.check(W->adjoint.Rows()==2 && W->adjoint.Cols()==3
                 && AV(W,0)==1 && AV(W,5)==1
                 && std::fabs(AV(x,0)-5)<1e-14 && std::fabs(AV(x,1)-7)<1e-14
                 && std::fabs(AV(x,2)-9)<1e-14,
                 std::string("[B9] MatrixMult W(2x3)·x(3x1) -> dW 全 1 / dx=[5,7,9]"));
    }
    {   // Function (原版矩阵乘, 抛异常; 且 Derive(Apply) 依赖种子)
        auto x = Leaf(MV({2,3}));
        auto y = OpF(GradOp::Function, {x}, {FunctionType::Sin}, 1, 1); y->Forward();
        y->adjoint = MV({1,1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(std::fabs(AV(x,0)-std::cos(2.0))<1e-14 &&
                 std::fabs(AV(x,1)-std::cos(3.0))<1e-14,
                 std::string("[B10] Function(sin) -> [cos2, cos3]"));
    }
    {   // MatrixFunction
        auto x = Leaf(MV({2,3}));
        auto y = OpF(GradOp::MatrixFunction, {x}, {FunctionType::Sin, FunctionType::Exp}, 2, 1);
        y->Forward();
        y->adjoint = MV({1,1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(std::fabs(AV(x,0)-std::cos(2.0))<1e-14 &&
                 std::fabs(AV(x,1)-std::exp(3.0))<1e-14,
                 std::string("[B11] MatrixFunction([sin,exp]) -> [cos2, exp3]"));
    }
    {   // Output (原版缺 case, 静默 no-op)
        auto x = Leaf(MV({1,2}));
        auto y = Op(GradOp::Output, {x}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x,0)==1 && AV(x,1)==1,
                 std::string("[B12] Output 透传 -> [1,1]  (原来静默 no-op)"));
    }
    {   // 叶子与常量: 显式 case, 无副作用
        auto x = Leaf(MV({1,2}));
        auto y = Op(GradOp::Constant, {x}, 2, 1); y->Forward();
        y->adjoint = MV({1,1}); x->adjoint = MV({0,0});
        y->Backward();
        ck.check(AV(x,0)==0 && AV(x,1)==0, std::string("[B13] Constant -> [0,0] (不回传)"));
    }
}

// ===========================================================================
// [C] adjoint 形状初始化 (原版此处必抛 Matrix dimension mismatch)
// ===========================================================================
static void TestC_AutoInit()
{
    std::printf("\n[C] ComputationGraph::Backward 自动做形状初始化\n");
    auto x = Leaf(MV({2,3}));
    auto y = Op(GradOp::Output, {x}, 2, 1);
    ComputationGraph<D> g(x, y);
    g.Forward();

    ck.check(x->adjoint.Size() == 0, std::string("[C1] Forward 后 adjoint 仍是空的 (与原来一样)"));

    bool threw = false;
    try { g.Backward(); } catch (const std::exception&) { threw = true; }
    ck.check(!threw, std::string("[C2] Backward 不再抛 (原版必抛)"));
    ck.check(x->adjoint.Size() == 2, std::string("[C3] adjoint 形状已由 Backward 自动建立"));
    ck.check(AV(x,0)==1 && AV(x,1)==1, std::string("[C4] 默认播种全 1 -> 透传出 [1,1]"));
}

// ===========================================================================
// [D] 幂等: 连跑两次结果不变
// ===========================================================================
static void TestD_Idempotent()
{
    std::printf("\n[D] 连跑两次 Backward\n");
    auto x = Leaf(MV({2, 3}));
    auto w = Leaf(MV({5, 7}));
    auto m = Op(GradOp::CWiseMult, {x, w}, 2, 1);
    auto y = Op(GradOp::SumAll, {m}, 1, 1);
    ComputationGraph<D> g(x, y);
    g.Forward();
    g.Backward();
    const double first = AV(x, 0);
    g.Backward();
    ck.check(std::fabs(AV(x,0) - first) < 1e-15 && std::fabs(first - 5.0) < 1e-13,
             std::string("[D1] 二次 Backward 仍得 dy/dx0 = 5"));
}

// ===========================================================================
// [E] 占位算子: 如实断言"仍不可用", 不假装可用
//     (将来实现后, 本段应改为数值对拍 —— FAIL 即提醒更新用例)
// ===========================================================================
static void TestE_Placeholders()
{
    std::printf("\n[E] 占位算子 (如实记录现状)\n");

    {   // CWisePow: Forward 与 Backward 都未实现
        auto x = Leaf(MV({2}));
        auto y = Op(GradOp::CWisePow, {x}, 1, 1);
        bool fwd_threw = false, bwd_threw = false;
        try { y->Forward(); } catch (const std::logic_error&) { fwd_threw = true; }
        y->adjoint = MV({1}); x->adjoint = MV({0});
        try { y->Backward(); } catch (const std::logic_error&) { bwd_threw = true; }
        ck.check(fwd_threw && bwd_threw, std::string("[E1] CWisePow: 前后向都明确 throw (不再静默)"));
    }
    {   // Inverse: Forward 仍为占位 (注释"待 linalg::Inverse 落地"已过期)
        auto x = Leaf(MV({3}));
        auto y = Op(GradOp::Inverse, {x}, 1, 1);
        bool threw = false;
        try { y->Forward(); } catch (const std::logic_error&) { threw = true; }
        ck.check(threw, std::string("[E2] Inverse: Forward 仍 throw (占位)"));
    }
    {   // Sigmoid: Forward throw, Backward 是**空 body** ⇒ 静默 0
        auto x = Leaf(MV({0.5}));
        auto y = Op(GradOp::Sigmoid, {x}, 1, 1);
        bool threw = false;
        try { y->Forward(); } catch (const std::logic_error&) { threw = true; }
        y->adjoint = MV({1}); x->adjoint = MV({0});
        y->Backward();          // 不抛
        ck.check(threw && AV(x, 0) == 0.0,
                 std::string("[E3] Sigmoid: Forward throw + Backward 空 body ⇒ 静默 0 (待实现)"));
    }
}

// ===========================================================================
// [F] 分域: 同一份 grad_graph 在 Dual<float> / Dual<Complex<double>> 下各实例化一次
//     (graphics 3.3 第 1 条: 模板 bug 只在实例化时暴露; 单域探针连续误判过 5 次)
//     2026-10-17 这两个域只在 /tmp 探针里跑过、没进 ctest —— 本条把证据固化。
// ===========================================================================
template <typename T>
static void RunMinLoopDomain(Checker& ck, const char* tag, double tol)
{
    using DQ = Dual<T>;
    using N  = GradNode<DQ>;

    auto mk = [](double v) { return ScalarTraits<T>::fromReal(v); };

    auto leaf = [&](const std::vector<T>& vs) {
        std::vector<DQ> d;
        d.reserve(vs.size());
        for (const T& v : vs) d.emplace_back(v, ScalarTraits<T>::zero());
        return std::make_shared<N>(GradOp::Input,
                                   Matrix<DQ>(static_cast<int>(vs.size()), 1, d));
    };

    auto x = leaf({mk(2), mk(3)});
    auto w = leaf({mk(5), mk(7)});
    auto b = leaf({mk(0.5)});

    auto m = std::make_shared<N>(GradOp::CWiseMult, Matrix<DQ>(2, 1));
    m->FillInput(x); m->FillInput(w);
    auto s = std::make_shared<N>(GradOp::SumAll, Matrix<DQ>(1, 1));
    s->FillInput(m);
    auto y = std::make_shared<N>(GradOp::Add, Matrix<DQ>(1, 1));
    y->FillInput(s); y->FillInput(b);

    ComputationGraph<DQ> g(x, y);
    g.Forward();
    g.Backward();

    auto near = [&](const T& got, double want) {
        if constexpr (std::is_same_v<T, Complex<double>>) {
            return std::fabs(static_cast<double>(got.real) - want) < tol
                && std::fabs(static_cast<double>(got.imag)) < tol;
        } else {
            return std::fabs(static_cast<double>(got) - want) < tol;
        }
    };

    const bool ok = near(x->adjoint.Data()[0].value, 5.0)
                 && near(x->adjoint.Data()[1].value, 7.0)
                 && near(w->adjoint.Data()[0].value, 2.0)
                 && near(w->adjoint.Data()[1].value, 3.0)
                 && near(b->adjoint.Data()[0].value, 1.0);
    ck.check(ok, std::string("[F] ") + tag + " 最小闭环 dx=[5,7] dw=[2,3] db=1");
}

// ===========================================================================

int main()
{
    TestA_MinLoop();
    TestB_PerOperator();
    TestC_AutoInit();
    TestD_Idempotent();
    TestE_Placeholders();

    std::printf("\n[F] 分域 (同一份源码换域再实例化)\n");
    RunMinLoopDomain<float>(ck, "Dual<float>            ", 1e-4);
    RunMinLoopDomain<Complex<double>>(ck, "Dual<Complex<double>>  ", 1e-12);

    return CheckReport(ck.ok, ck.fail, "autodiff/backward");
}
