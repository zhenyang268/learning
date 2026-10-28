#pragma once

#include <map>
#include <deque>
#include <stack>
#include <memory>
#include <utility>
#include <cassert>
#include <string>
#include <stdexcept>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "autodiff/jacobian.hpp"

// 使用无作用域枚举会与 GradOp 的枚举值(Input/Data/Function)冲突, 改为作用域枚举
enum class Type
{
    Input,      // 输入层,
    Data,       // 中间数据
    Parameter,  // 参数节点, 比如常数和常量矩阵
    Operation,  // 加减乘除累加累乘等等
    Function,   // 函数, sin cos等等分为应用函数和矩阵函数两种 
};

using std::map;
using std::deque;
using std::stack;
using std::pair;
using std::string;
using std::vector;

// 二元算子表: 与 GradNode::Forward 的 switch 一一对应 (四则 + 逐元素幂 + 矩阵乘)。
// 节点 arity 由 inputs.size() 运行时携带, 缺输入时
// node1/node2 是 nullptr —— 靠这张表在 Forward 入口先拦住, 不然解引用就是段错误。
inline bool IsBinaryOp(GradOp op)
{
    switch (op)
    {
        case GradOp::Add:
        case GradOp::Sub:
        case GradOp::Mult:
        case GradOp::CWiseMult:
        case GradOp::CWisePow:
        case GradOp::Divide:
        case GradOp::MatrixMult:
            return true;
        default:
            return false;
    }
}

/* 布局: 列向量 + 分子布局雅可比 (D19 第二次修订, 2026-10-19)
//
// 假设：y = W * x，其中 x ∈ R^{n×1}, W ∈ R^{m×n}, y ∈ R^{m×1}
// 向量 y 对向量 x 的导数是一个 m×n 的矩阵（雅可比）, 这里直接就是 W。
//
// 分子布局 (Numerator Layout) —— 本库采用:
//   J_ij = ∂y_i/∂x_j, 形状 = [y 的维度, x 的维度] = m × n
//   链式法则自然序: J_{f∘g} = J_f · J_g
//
// 分母布局 (Denominator Layout) —— 旧口径, 2026-10-19 翻案前使用:
//   J_ij = ∂y_j/∂x_i, 形状 = n × m
//   链式法则镜像: J_{f∘g} = J_g · J_f
//
// ⚠ 本文件 Backward 的伴随公式**与布局无关**: y = AB 的 Ā = ā·Bᵀ、B̄ = Aᵀ·ā
// 在两种约定下形式完全相同 (只换"谁是 A 谁是 B")。翻案改的是**调用方建图的顺序**
// (以前写 X*W, 现在写 W*x) 与 jacobian.hpp 的 Jacobian() 组装, 不是这里的公式。
*/

template <typename T>
class GradNode : public std::enable_shared_from_this<GradNode<T>>
{
public:
    string name;    // "op+层数"用于debug
    GradOp op;      // 描述节点操作
    Type nodeType;  // 节点的类型
    // 不要单值, 直接存到matrix里面, 加减乘除都要对矩阵逐元素操作
    Matrix<T> matrix; // matrix存储所有的值, 函数值, 对某个xi的偏导数值
    Matrix<T> adjoint; // 反向时使用
    
    Matrix<FunctionType> matrixF; // 存储操作函数, 分为1个和n*n个函数两种
    vector<std::shared_ptr<GradNode<T>>> inputs; // 依赖节点

private:
    void op2type()
    {
        if (op == GradOp::Input || op == GradOp::Output) {
            nodeType = Type::Data;
        } else if (op == GradOp::Constant) {
            nodeType = Type::Parameter;
        } else {
            nodeType = Type::Function;
        }
    }

public:
    GradNode(GradOp op, Matrix<T> matrix):
            op(op), matrix(std::move(matrix)) { op2type(); }

    GradNode(GradOp op, Matrix<FunctionType> matrixF):
            op(op), matrixF(std::move(matrixF)) { op2type(); }

    std::shared_ptr<GradNode<T>> GetSelf()
    {
        return this->shared_from_this();
    }
    
    void FillInput(const std::shared_ptr<GradNode<T>> &input)
    {
        inputs.push_back(input);
    }

    void Forward()
    {
        // arity 校验放在最前: 缺输入时下面拿到的是 nullptr, 不拦就是段错误(无从定位)
        if (IsBinaryOp(op)) {
            if (inputs.size() != 2)
                throw std::invalid_argument(
                    "GradNode::Forward: 二元算子需要 2 个输入, 实际 " +
                    std::to_string(inputs.size()));
        } else if (op != GradOp::Input && op != GradOp::Constant) {
            if (inputs.size() != 1)
                throw std::invalid_argument(
                    "GradNode::Forward: 一元/透传算子需要 1 个输入, 实际 " +
                    std::to_string(inputs.size()));
        }

        std::shared_ptr<GradNode<T>> node1 = nullptr, node2 = nullptr;
        if (inputs.size() == 2) {
            node1 = inputs[0];
            node2 = inputs[1];
        } else if (inputs.size() == 1) {
            node1 = inputs[0];
        }
        switch (op) {
            // 输入节点value, derive以及要求的xi都由外部提前设置好
            case GradOp::Input: 
                break;

            // 常数/常量矩阵节点: 值由外部预置, 前向不动
            case GradOp::Constant:
                break;

            // 输出透传节点: 原样搬运输入矩阵(含 value 与 derive)
            case GradOp::Output:
                matrix = node1->matrix;
                break;
            
            // 交叉熵: 未实现
            case GradOp::CrossEntropy:
                throw std::logic_error("CrossEntropy: 未实现");

            case GradOp::Add:
                matrix = std::move(node1->matrix + node2->matrix);
                break;
            case GradOp::Sub:
                matrix = std::move(node1->matrix - node2->matrix);
                break;
            // 逐元素乘法 (Mult 与 CWiseMult 同义: 标量/同形矩阵; 四则内建的乘)
            case GradOp::Mult:
            case GradOp::CWiseMult:
                matrix = std::move(node1->matrix.cwiseProduct(node2->matrix));
                break;

            // 逐元素幂运算: 未实现 (应建图拆链为 exp(y*ln x))
            case GradOp::CWisePow:
                throw std::logic_error("CWisePow: 未实现");
            
            // 矩阵乘法: y = node1 · node2 —— 列向量约定下 node1 = W 在左, node2 = x 在右。
            // Dual 的乘积法则对矩阵乘同样成立 (value/derive 各自走矩阵乘), 无需特殊处理。
            case GradOp::MatrixMult:
                matrix = std::move(node1->matrix * node2->matrix);
                break;
            case GradOp::Divide:
                matrix = std::move(node1->matrix.cwiseDivision(node2->matrix));
                break;

            // transpose也是线性变换
            case GradOp::Transpose:
                matrix = std::move(node1->matrix.Transpose());
                break;

            // inverse要通过伴随内积配对来算
            // <a_x, dx> = <a_y, dy>
            // z = f(x) = g(y) (y和x是换元)
            case GradOp::Inverse:
                // 未实现。注意: `linalg::Inverse` **早已落地** (V1.2, 2026-10-14),
                // 原注释"待 linalg::Inverse 落地"已过期 —— 现在的阻塞不是依赖,
                // 而是本节点没接入 (矩阵求逆的伴随公式本身已写好, 见 Backward 同名 case)。
                throw std::logic_error("Inverse: 本节点未接入 (linalg::Inverse 已可用)");

            // 激活函数: 未实现 (Sigmoid/Tanh 应拆链,
            // Relu/Softmax 因数值保留 fused 枚举)
            case GradOp::Sigmoid:
                throw std::logic_error("Sigmoid: 未实现");
            case GradOp::Relu:
                throw std::logic_error("Relu: 未实现");
            case GradOp::Tanh:
                throw std::logic_error("Tanh: 未实现");
            case GradOp::Softmax:
                throw std::logic_error("Softmax: 未实现");
            

            case GradOp::SumAll:
                matrix = Matrix<T>(1, 1);
                matrix[0][0] = SumAll(node1->matrix);
                break;
            case GradOp::MultAll:
                matrix = Matrix<T>(1, 1);
                matrix[0][0] = MultAll(node1->matrix);
                break;
            case GradOp::LnMultAll:
                throw std::logic_error("LnMultAll: 未实现");
            
            /**
             * n维输入, n维输出, 但是雅可比矩阵是n*n
             * 所以这里每个函数的微分都在矩阵的对角线上, 可以用原n维矩阵存储
             */
            case GradOp::Function:
                if (matrixF.Size() == 0)
                    throw std::logic_error("Function: matrixF 未设置 (空 0x0), 不能取 matrixF[0][0]");
                matrix = std::move(Apply(node1->matrix, matrixF[0][0]));
                break;
            case GradOp::MatrixFunction:
                if (matrixF.Size() == 0)
                    throw std::logic_error("MatrixFunction: matrixF 未设置 (空 0x0)");
                matrix = std::move(Apply(node1->matrix, matrixF));
                break;
        }

    }

    // 反向模式尚未接入 V1 基线: 本函数只保证可编译, 语义待 tape 立项后再定
    void Backward()
    {
        std::shared_ptr<GradNode<T>> node1 = nullptr, node2 = nullptr;
        Matrix<T> adjoint1, adjoint2;
        if (inputs.size() == 2) {
            node1 = inputs[0];
            node2 = inputs[1];
            
        } else if (inputs.size() == 1) {
            node1 = inputs[0];
        }
        /**
         * adjoint只使用value, 不用derive, 加上traits分离太麻烦
         */
        switch (op) {
            // ---- 叶子与常量: 无下游可传 (显式列出, 也为消 -Wswitch) ----
            case GradOp::Input:
            case GradOp::Constant:
                break;

            case GradOp::CrossEntropy:
                break;

            // ---- 透传节点: 伴随原样交给输入 ----
            case GradOp::Output:
                node1->adjoint += adjoint;
                break;

            case GradOp::Add:
                adjoint1 = Matrix<T>::Ones(node1->matrix.Rows(),
                                           node1->matrix.Cols());
                adjoint2 = Matrix<T>::Ones(node2->matrix.Rows(),
                                           node2->matrix.Cols());
                node1->adjoint += adjoint1.cwiseProduct(adjoint);
                node2->adjoint += adjoint2.cwiseProduct(adjoint);
                break;
            case GradOp::Sub:
                adjoint1 = Matrix<T>::Ones(node1->matrix.Rows(),
                                           node1->matrix.Cols());
                adjoint2 = Matrix<T>::Ones(node2->matrix.Rows(),
                                           node2->matrix.Cols())
                               * (-ScalarTraits<T>::one());
                node1->adjoint += adjoint1.cwiseProduct(adjoint);
                node2->adjoint += adjoint2.cwiseProduct(adjoint);
                break;
            // Mult 与 CWiseMult 同义 (Forward 里就是同一个 cwiseProduct) ——
            // 原先 Backward 只写了 CWiseMult, 于是 `Mult` 节点反向**静默 no-op**。
            case GradOp::Mult:
            case GradOp::CWiseMult:
                // y = A cwiseproduct X
                // 逐元素乘法, a_x = a_y cwiseproduct A
                node1->adjoint += node2->matrix.cwiseProduct(adjoint);
                node2->adjoint += node1->matrix.cwiseProduct(adjoint);
                break;
            case GradOp::MatrixMult:
                // y = A·B  (列向量约定: A 在左, 例如 W·x; 行向量时代是 X·W)
                // dy = dA·B + A·dB
                // <ā_y, dy> = tr(ā_yᵀ dA B) + tr(ā_yᵀ A dB)
                //           = tr(B ā_yᵀ dA) + tr(ā_yᵀ A dB)
                // ⇒ Ā = ā_y·Bᵀ,   B̄ = Aᵀ·ā_y
                // **这两行与布局无关**: y = X·W (行) 与 y = W·x (列) 共用同一式, 变的只是
                // 哪个节点当 A、哪个当 B —— 所以 2026-10-19 的翻案没有动这里。
                node1->adjoint += adjoint * node2->matrix.Transpose();
                node2->adjoint += node1->matrix.Transpose() * adjoint;
                break;
            case GradOp::Divide:
                // 逐元素除法: ∂(x1/x2)/∂x1 = 1/x2 ⇒ 逐元素除, 不是矩阵乘
                node1->adjoint += adjoint.cwiseDivision(node2->matrix);
                adjoint2 = Matrix<T>::Ones(node2->matrix.Rows(),
                                           node2->matrix.Cols())
                               * (-ScalarTraits<T>::one());
                node2->adjoint += adjoint2.cwiseProduct(node1->matrix)
                                      .cwiseDivision(node2->matrix)
                                      .cwiseDivision(node2->matrix);
                break;
            case GradOp::Transpose:
                // 原式 `node1->adjoint = node1->matrix.Transpose()` 两头都错:
                // 既用原值 (而非上游伴随) 又用 `=` 覆盖 (而非累加)。
                node1->adjoint += adjoint.Transpose();
                break;

            // a_x = - Y^T * a_y * Y^T
            case GradOp::Inverse:
                adjoint1 = matrix.Transpose() * adjoint * matrix.Transpose();
                adjoint1 = adjoint1 * (-ScalarTraits<T>::one());
                node1->adjoint += adjoint1;
                break;

            case GradOp::Sigmoid:
                break;
            case GradOp::Relu:
                break;
            case GradOp::Tanh:
                break;
            case GradOp::Softmax:
                break;            

            // a_x = J * a_y, J是全1矩阵
            case GradOp::SumAll:
                adjoint1 = Matrix<T>::Ones(node1->matrix.Rows(),
                                           node1->matrix.Cols());
                node1->adjoint += adjoint1 * adjoint[0][0];
                break;

            // a_x = J * a_y, Jij是除xij的累乘
            case GradOp::MultAll:
                // MultAllBackward 返回的 Dual 把「总积」放在 value、∂总积/∂x_i 放在 derive,
                // 故必须先 Derive 把偏导提到 value 分量, 再与上游伴随逐元素乘。
                adjoint1 = Derive(MultAllBackward(node1->matrix));
                // MultAll 的输出恒为 1x1 ⇒ 上游伴随是**标量**, 这里是「标量乘向量」而非逐元素乘
                node1->adjoint += adjoint1 * adjoint[0][0];
                break;
            case GradOp::LnMultAll:
                break;

            case GradOp::CWisePow:
                throw std::logic_error("CWisePow::Backward: 未实现 (Forward 亦未实现)");
            
            // 逐元素函数的雅可比是对角阵 (注释里本来就写着) ⇒ 必须**逐元素乘**;
            // 且局部导数要用 LocalDeriv (与种子无关), 不能拿 Apply 的 derive (= f'(x)·种子)。
            case GradOp::Function:
                node1->adjoint += adjoint.cwiseProduct(LocalDeriv(node1->matrix, matrixF[0][0]));
                break;
            case GradOp::MatrixFunction:
                node1->adjoint += adjoint.cwiseProduct(LocalDeriv(node1->matrix, matrixF));
                break;
        }
    }
};

/**
 * 计算拓扑顺序
 * 1. 前向拓扑顺序, 计算函数值以及记录每个节点对前序节点的梯度矩阵
 * 2. 后向拓扑排序, 根据梯度矩阵链式法则计算对input和中间参数的梯度
 */
template <typename T>
class ComputationGraph
{
using GradNodePtr = std::shared_ptr<GradNode<T>>;
public:
    GradNodePtr input;
    GradNodePtr output;

    std::deque<GradNodePtr> forward_list;
    std::deque<GradNodePtr> backward_list;

    ComputationGraph(GradNodePtr input, GradNodePtr output):
        input(input), output(output)
    {
        // 初始化两个列表
        Topology(output);
        assert(!forward_list.empty());
        // 反向列表 = 正向列表逆序(反向模式尚未实现, 见 GradNode::Backward)
        backward_list.assign(forward_list.rbegin(), forward_list.rend());
    }

private:
    /**
     * 删除后续节点outputs, 直接reverse正向排序实现
     * 栈模拟dfs, 返回input
     */
    void Topology(GradNodePtr start)
    {
        stack<pair<GradNodePtr, bool>> st;
        st.push({start, false});
        /**
         * 核心逻辑
         * 处理st.top [parent.input[i]]
         * 当st.top 子节点都就绪时就可以出栈插入list
         * 当st.top 未就绪时, 就按访问顺序插入栈中
         */
        while (!st.empty())
        {
            auto [cur, finished] = st.top();
            st.pop();
            if (finished) {
                forward_list.push_back(cur); // input -> output
            } else {
                st.push({cur, true}); // 标记节点已完成, 按顺序访问子节点
                for (auto it = cur->inputs.rbegin(); it != cur->inputs.rend(); ++it) {
                    st.push({*it, false});
                }
            }
        }
    }

public:
    void Forward()
    {
        // 先input最后output
        for (auto node : forward_list)
        {
            node->Forward();
        }
    }

    void Backward()
    {
        Backward(Matrix<T>());      // 默认播种: 输出伴随填全 1 (= ∂Σy/∂y)
    }

    // seed: 输出节点的伴随初值 (形状须与 output->matrix 一致); 传 0x0 表示用全 1
    void Backward(const Matrix<T>& seed)
    {
        // (1) 形状初始化 —— 这正是原版"一跑就抛"的根因:
        //     adjoint 默认 0x0, 而下面所有分支都写 `node->adjoint += ...`, 形状不匹配必抛。
        //     Forward 已跑过 ⇒ 各节点 matrix 形状已定, 按它重建即可 (顺带清零)。
        for (auto& node : forward_list) {
            node->adjoint = Matrix<T>(node->matrix.Rows(), node->matrix.Cols());
        }

        // (2) 输出节点播种
        if (output->matrix.Size() > 0) {
            if (seed.Size() == output->matrix.Size()) {
                output->adjoint = seed;
            } else {
                output->adjoint = Matrix<T>::Ones(output->matrix.Rows(),
                                                  output->matrix.Cols());
            }
        }

        // (3) 拓扑逆序逐节点回传
        for (auto& node : backward_list) {
            node->Backward();
        }
    }
};
// ===========================================================================
// 四、数学公式与符号约定
// ===========================================================================
//
// **本文件所有公式集中在这里, 实现处只留"为什么这么写"的工程理由。**
//
// ---------------------------------------------------------------------------
// 符号表 (代码里出现的每个符号都在这里)
// ---------------------------------------------------------------------------
//   G              计算图 (ComputationGraph)
//   Node           图中的一个节点 (GradNode)
//   GradOp         节点运算类型枚举 (Add / Mul / MatMul / Sin / …)
//   x, y, z        节点 (或张量)
//   L              标量损失 (loss), 最终要对其求梯度
//   ∂L/∂x          L 对 x 的**偏导** (梯度)
//   x̄, ȳ           adjoint (伴随): 记作 ∂L/∂x, ∂L/∂y
//   v              上游传来的梯度 (upstream gradient)
//   W              矩阵 (MatMul 的左操作数)
//   J              雅可比
//   ε              机器 epsilon
//
// **关于 adjoint 记号**: 本文件代码里叫 `adjoint` 的那个量,
// 数学上就是 **∂L/∂(该节点)** —— 即"总损失对这个节点输出的敏感度"。
// 它不是"伴随矩阵"(adjugate matrix), 两者只是中文都译作"伴随", 别混。
//
// ---------------------------------------------------------------------------
// (G1) 反向模式的基本原理 (链式法则的递归形式)
// ---------------------------------------------------------------------------
//        设 L = f(y),  y = g(x)。链式法则:
//              ∂L/∂x = (∂L/∂y) · (∂y/∂x)
//
//        反向模式就是**从 L 往回走**, 每到一条边就把上游梯度乘上
//        这条边的局部导数, 再累加到下游节点。
//
//        形式化: 对每个节点 n, 累加
//              ∂L/∂n = Σ_{m: n 是 m 的输入}  (∂L/∂m) · (∂m/∂n)
//
//        这就是为什么一个节点被多个节点使用时, 梯度要**累加**而不是覆盖。
//
// ---------------------------------------------------------------------------
// (G2) 反向 vs 前向: 什么时候用哪个
// ---------------------------------------------------------------------------
//        前向模式 (dual.hpp): 一次求值得到 J 的**一列**   → 代价 O(n)
//        反向模式 (本文件): 一次回传得到 J 的**一行**   → 代价 O(m)
//
//        其中 n = 输入维数, m = 输出维数。
//
//        ⇒ **n >> m 用反向** (这是深度学习的情形: 参数几百万, 损失一个标量)。
//          **m >> n 用前向** (这是本仓多数几何运算的情形)。
//
//        两者得到的 J 是**同一个矩阵**, 只是填充方向不同。
//
// ---------------------------------------------------------------------------
// (G3) 各 GradOp 的局部导数 (反传公式)
// ---------------------------------------------------------------------------
//        记上游传来的梯度为 v = ∂L/∂(本节点输出)。
//
//        Add:      z = x + y        ⇒  ∂L/∂x = v,        ∂L/∂y = v
//        Sub:      z = x − y        ⇒  ∂L/∂x = v,        ∂L/∂y = −v
//        Mul:      z = x · y        ⇒  ∂L/∂x = v·y,      ∂L/∂y = v·x
//        Div:      z = x / y        ⇒  ∂L/∂x = v/y,      ∂L/∂y = −v·x/y²
//        MatMul:   z = W · x        ⇒  ∂L/∂W = v·xᵀ,     ∂L/∂x = Wᵀ·v
//        Transpose:z = xᵀ           ⇒  ∂L/∂x = vᵀ
//        Neg:      z = −x           ⇒  ∂L/∂x = −v
//        Sin:      z = sin x        ⇒  ∂L/∂x = v·cos x
//        Cos:      z = cos x        ⇒  ∂L/∂x = −v·sin x
//        Exp:      z = exp x        ⇒  ∂L/∂x = v·exp x
//        Log:      z = log x        ⇒  ∂L/∂x = v/x
//        Sqrt:     z = √x           ⇒  ∂L/∂x = v/(2√x)
//        Inverse:  z = x⁻¹         ⇒  ∂L/∂x = −x⁻ᴴ·v·x⁻ᴴ
//
//        **这组公式与布局无关** —— 见 (G5)。
//
// ---------------------------------------------------------------------------
// (G4) 为什么 MatMul 的 ∂L/∂W = v·xᵀ (形状推导)
// ---------------------------------------------------------------------------
//        设 W 是 m×n, x 是 n×1, z = W·x 是 m×1。
//        上游 v = ∂L/∂z 是 m×1。
//
//        ∂L/∂W 必须与 W **同形状** (m×n)。而 v·xᵀ 是 (m×1)·(1×n) = m×n ✓
//
//        ∂L/∂x 必须与 x 同形状 (n×1)。而 Wᵀ·v 是 (n×m)·(m×1) = n×1 ✓
//
//        ⇒ **形状唯一确定了这两个公式** —— 不需要背, 从形状就能推出来。
//          维度对不上的候选公式一定是错的。
//
// ---------------------------------------------------------------------------
// (G5) 反传公式与布局无关 (2026-10-19 的结论)
// ---------------------------------------------------------------------------
//        (G3) 里的公式写的是**单个向量**的链式法则, 不涉及"行向量还是列向量"
//        的约定。所以:
//
//        · 前向 (D19 改列向量 + 分子布局) 只影响 **J 的形状与填充方向**,
//          不影响 adjoint 的递推关系。
//        · 本文件的 Backward 代码在布局翻案时**一行没改** —— 这不是巧合,
//          是因为它算的是 adjoint 而不是 J。
//
//        ⚠ 容易混的地方: `MatrixMult` 节点的 Forward 写 `W 在左, x 在右`,
//          这是 **W·x**; 如果哪天要支持 x·W, 那是另一个 GradOp, 不是换个参数。
//
// ---------------------------------------------------------------------------
// (G6) 为什么必须有 tape (不能"一直往前传")
// ---------------------------------------------------------------------------
//        前向传播时如果只保留"当前值", 回传时就需要**重新算出**每个中间值。
//        对于 sin(x) 这类节点, 重新算 cos(x) 是可行的; 但对于
//        被多个下游节点共用的值 (fan-out), 重算会指数爆炸。
//
//        tape 的作用是**记住前向的执行顺序与中间值**, 于是:
//        · 回传只需一次线性扫描 (逆序)
//        · 每个节点的梯度只算一次, 然后累加到所有下游
//
//        ⇒ tape 不是"为了省事", 而是把 O(指数) 的重算降到 O(节点数)。
//
// ---------------------------------------------------------------------------
// (G7) 梯度检查 (独立判据)
// ---------------------------------------------------------------------------
//        ∂L/∂x_i ≈ ( L(x + h·e_i) − L(x − h·e_i) ) / (2h)      ← 中心差分
//
//        这是验证反向实现的标准方法: **它完全不用 tape, 不用 adjoint**,
//        只把 x 扰动一下重算前向, 所以是真正的独立路径。
//
//        h 的选取: double 下 h = 1e-6, 误差约 1e-10 量级。
//        ⚠ 这不是"用有限差分代替反向" —— 只是拿它当判据。
//          有限差分的代价是 O(n) 次前向, 反向是 O(1) 次, 所以生产环境只用反向。
//
// ---------------------------------------------------------------------------
// 实现口径 (跟公式直接相关)
// ---------------------------------------------------------------------------
//   · adjoint 初值: 只有最终节点 (L) 的 adjoint 置 1, 其余全 0。
//     ⚠ 如果初始化成"全 1"或别的值, 结果会错但**不一定明显** ——
//       这是本项目踩过的坑之一。判据: 对 L = x 这种恒等图, 梯度必须恰好是 1。
//   · 累加而不是覆盖 (见 G1): 一个节点被多处使用时, 每次回传都要 +=。
//   · 节点 arity (输入个数) 由 inputs.size() 运行时决定, 不写死在枚举里。
// ===========================================================================
