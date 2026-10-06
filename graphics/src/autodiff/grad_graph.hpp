#pragma once

#include <map>
#include <deque>
#include <stack>
#include <memory>
#include <utility>
#include <cassert>
#include <string>
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
        if (op >= GradOp::Input && op < GradOp::Constant) {
            nodeType = Type::Data;
        } else if (op == GradOp::Constant) {
            nodeType = Type::Parameter;
        } else {
            nodeType = Type::Function;
        }
        //     nodeType = Type::Operation; // 本质和函数没区别, 但是这里固定住不修改
    }

public:
    GradNode(GradOp op, Matrix<T> matrix):
            op(op), matrix(std::move(matrix)) { op2type(); }

    GradNode(GradOp op, Matrix<DualFunc<T>> matrixF):
            op(op), matrixF(matrixF) { op2type(); }

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
        std::shared_ptr<GradNode<T>> node1 = nullptr, node2 = nullptr;
        if (inputs.size() == 2) {
            node1 = inputs[0];
            node2 = inputs[1];
        } else if (inputs.size() == 1) {
            node1 = inputs[0];
        }
        switch (op) {
            // 输入节点value, derive以及要求的xi都由外部提前设置好
            case Input: 
                break;
            
            case CrossEntropy:
                break;

            case Add:
                matrix = std::move(node1->matrix + node2->matrix);
                break;
            case Sub:
                matrix = std::move(node1->matrix - node2->matrix);
                break;
            case CWiseMult:
            // 逐元素乘法
                matrix = std::move(node1->matrix.CWiseProduct(node2->matrix));
                break;
            
            // 矩阵乘法也满足derive的计算规则
            // 但是需要注意node1和node2的顺序, node1入参数行主导
            // node2的derive都是0
            case MatrixMult:
                matrix = std::move(node1->matrix * node2->matrix);
                break;
            case Devide:
                matrix = std::move(node1->matrix / node2->matrix);
                break;

            // transpose也是线性变换
            case Transpose:
                matrix = std::move(node1->Transpose());
                break;

            // inverse要通过伴随内积配对来算
            // <a_x, dx> = <a_y, dy>
            // z = f(x) = g(y) (y和x是换元)
            case Inverse:
                matrix = std::move(node1->Inverse());
                break;

            case Sigmoid:
                break;
            case Relu:
                break;
            case Tanh:
                break;
            case Softmax:
                break;
            

            case SumAll:
                matrix = Matrix(1, 1);
                matrix[0][0] = SumLine(node1->matrix);
                break;
            case MultAll:
                matrix = Matrix(1, 1);
                matrix[0][0] = MultLine(node1->matrix);
                break;
            case LnMultAll:
                break;
            
            /**
             * n维输入, n维输出, 但是雅可比矩阵是n*n
             * 所以这里每个函数的微分都在矩阵的对角线上, 可以用原n维矩阵存储
             */
            case Function:
                matrix = std::move(Apply(node1->matrix, matrixF[0][0]));
                break;
            case MatrixFunction:
                matrix = std::move(Apply(node1->matrix, matrixF));
                break;
        }

    }

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
            case CrossEntropy:
                break;
            case Add:
                adjoint1 = Matrix(node1->matrix.Rows(),
                                        node1->matrix.Cols(),
                                        ScalarTraits<T>::one());
                adjoint2 = Matrix(node2->matrix.Rows(),
                                        node2->matrix.Cols(),
                                        ScalarTraits<T>::one());
                node1->adjoint += adjoint1.CWiseProduct(adjoint);
                node2->adjoint += adjoint2.CWiseProduct(adjoint);
                break;
            case Sub:
                adjoint1 = Matrix(node1->matrix.Rows(), 
                                        node1->matrix.Cols(),
                                        ScalarTraits<T>::one());
                adjoint2 = Matrix(node2->matrix.Rows(),
                                        node2->matrix.Cols(),
                                        ScalarTraits<T>::zero() - ScalarTraits<T>::one());
                node1->adjoint += adjoint1.CWiseProduct(adjoint);
                node2->adjoint += adjoint2.CWiseProduct(adjoint);
                break;
            case CWiseMult:
                // y = A cwiseproduct X
                // 逐元素乘法, a_x = a_y cwiseproduct A
                node1->adjoint += node2->matrix.CWiseProduct(adjoint);
                node2->adjoint += node1->matrix.CWiseProduct(adjoint);
                break;
            case MatrixMult:
                // dY = dX * W + x * dW
                // <a_x, dx> = <a_y, dy> = tr(dy^T * a_y), dW = 0
                // <J^*  * a_y, dx> = <a_y, dx * J>
                // 推出a_X = a_Y * W^T
                // 同理dX = 0
                // a_W = X^T * a_y
                // 这里如果是X * X, 恰好也成立
                node1->adjoint += adjoint * node2->matrix.transpose();
                node2->adjoint += node1->matrix.transpose() * adjoint;
                break;
            case Devide:
                node1->adjoint += adjoint * node2->matrix.transpose();
                adjoint2 = Matrix(node2->matrix.Rows(),
                                        node2->matrix.Cols(),
                                        ScalarTraits<T>::zero() - ScalarTraits<T>::one());
                node2->adjoint += adjoint2.CWiseProduct(node1->matrix) / node2->matrix / node2->matrix;
                break;
            case Transpose:
                node1->adjoint = std::move(node1->Transpose());
                break;
            
            // a_x = - Y^T * a_y * Y^T
            case Inverse:
                adjoint1 = matrix.Transpose() * adjoint * matrix.Transpose(); 
                adjoint1.CWisePoint(ScalarTraits<T>::zero() - ScalarTraits<T>::one());
                node1->adjoint += adjoint1;
                break;

            case Sigmoid:
                break;
            case Relu:
                break;
            case Tanh:
                break;
            case Softmax:
                break;            

            // a_x = J * a_y, J是全1矩阵
            case SumAll:
                adjoint1 = Matrix(node1->matrix.Rows(),
                                        node1->matrix.Cols(),
                                        ScalarTraits<T>::one());
                node1->adjoint += adjoint1.CWiseProduct(adjoint[0][0]);
                break;
            
            // a_x = J * a_y, Jij是除xij的累乘
            case MultAll:
                adjoint1 = Derive(MultAll(node1->matrix));
                node1->adjoint += adjoint1.CWiseProduct(adjoint);
                break;
            case LnMultAll:
                break;
            
            case Function:
                node1->adjoint += adjoint * Derive(Apply(node1->matrix, matrixF[0][0]));
                break;
            case MatrixFunction:
                node1->adjoint += adjoint * Derive(Apply(node1->matrix, matrixF));
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
        // 先output最后input
        GradNodePtr last = nullptr;
        for (auto node : backward_list)
        {
            node->Backward(last);
            last = node;
        }
    }
};