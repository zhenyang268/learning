#pragma once

#include <map>
#include <deque>
#include <stack>
#include <memory>
#include <utility>
#include <cassert>
#include <string>
#include <stdexcept>
#include "autodiff/dual.hpp"
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
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
// 节点 arity 由 inputs.size() 运行时携带(matrix.md 7.5 第 6 条), 缺输入时
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

/* 使用行向量+分母布局, ∂z/∂x = ∂z/∂y x ∂y/∂x = W2 * W1
// 假设：Y = X * W，其中 X ∈ R^{1×n}, W ∈ R^{n×m}, Y ∈ R^{1×m}
// 向量 y 对向量 x 的导数是一个矩阵（雅可比矩阵）, 这里直接就是 W

// 分子布局（Numerator Layout）：
// ∂y/∂x 的形状 = [y的维度, x的维度] = m × n

// 分母布局（Denominator Layout）：
// ∂y/∂x 的形状 = [x的维度, y的维度] = n × m
// 这就是转置的来源！
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

            // 逐元素幂运算: 未实现(按 matrix.md 7.5 第 6 条应建图拆链 exp(y*ln x))
            case GradOp::CWisePow:
                throw std::logic_error("CWisePow: 未实现");
            
            // 矩阵乘法也满足derive的计算规则
            // 但是需要注意node1和node2的顺序, node1入参数行主导
            // node2的derive都是0
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
                // TODO: 阶段 1 接 linalg::Inverse
                throw std::logic_error("Inverse: 待 linalg::Inverse 落地 (阶段 1)");

            // 激活函数: 未实现(按 matrix.md 7.5 第 6 条, Sigmoid/Tanh 应拆链,
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
                // dY = dX * W + x * dW
                // <a_x, dx> = <a_y, dy> = tr(dy^T * a_y), dW = 0
                // <J^*  * a_y, dx> = <a_y, dx * J>
                // 推出a_X = a_Y * W^T
                // 同理dX = 0
                // a_W = X^T * a_y
                // 这里如果是X * X, 恰好也成立
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