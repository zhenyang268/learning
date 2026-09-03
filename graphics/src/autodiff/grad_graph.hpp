#pragma once

#include <map>
#include <deque>
#include <stack>
#include <memory>
#include <utility>
#include <cassert>
#include <string>
#include "dual.hpp"
#include "core/common.hpp"
#include "core/matrix.hpp"

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
// 向量 y 对向量 x 的导数是一个矩阵（雅可比矩阵）

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
    bool isDesent = false; // 是否是梯度下降, 只对常量矩阵求偏微分
    string name;    // "op+层数"用于debug
    GradOp op;      // 描述节点操作
    Type nodeType;  // 节点的类型
    // 不要单值, 直接存到matrix里面, 加减乘除都要对矩阵逐元素操作
    Matrix<T> matrix; // matrix存储所有的值, 函数值, 偏导数值
    Matrix<T> desent; // 这个是进行gradient desent时使用
    /**
     * 雅可比矩阵可以通过vjp方法简化, 比如Y = X * W, 雅可比矩阵就是W, 没必要再存一遍
     * Y = sin(X), 雅可比矩阵是对角矩阵Aii = cos(xi) 也可以当场计算
     * 
     * 当场计算而不是额外存储, 当然我们只使用一点vjp方法的思想
     */
    
    /**
        * 
        * 再从损失函数自顶向下计算到当前节点的偏微分值
        * matrix是链式法则计算用的关键矩阵, 保存节点操作对应的(node->value, node->derive)
        * desent只有在记录自变量梯度时才有用, x_input节点里面, 也应该用desent存一份
        * desent就是(v, 1)矩阵
        * 参数 x_input都是输入的第一层入口, 参数只有一条路径, x可能多条
    */
    /**
     * 从input正向计算节点的函数值
     * 计算每个节点对input的偏微分值, 也就是雅可比矩阵
     * Y = X * W
     * 在计算对W的偏微分用于梯度下降时就相当于W是input, 固定x
     * desent矩阵正向时用来保存雅可比矩阵, 反向时用来保存链式法则后的结果
     * 损失函数只有一个y, 所以雅可比矩阵我们不严谨, y对W矩阵的偏微分并不是一列
     * 但是其实正好矩阵形状就是我们想要的
     */
    
    
    Matrix<DualFunc<T>> matrixF; // 存储操作函数
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
            case Input:
                // 叶子节点: value 由外部设置, derive不能让外部设置, 固定为1
                // (forward 模式逐输入求梯度, 这里不做覆盖)
                // value预设, derive设为1
                Apply(matrix, [](Dual<T>& d) { return setDeriveOne(d); });
                break;
            case Output:
                matrix = Matrix<T>(1, 1);
                matrix[0][0] = node1->matrix[0][0];
                // 损失函数y拓展成n维, 用于梯度计算
                // 输出节点只保留一个值y
                break;
            case CrossEntropy:
                break;
            case Constant:
                // 常量: matrix(v, 0), desent(v, 1)
                // setDeriveZero/One 是 Dual 的友元函数, 只能靠 ADL 在调用点解析,
                // 因此用 lambda 包装成 DualFunc 再传入 Apply
                Apply(matrix, [](Dual<T>& d) { return setDeriveZero(d); });
                Apply(desent, [](Dual<T>& d) { return setDeriveOne(d); });
                break;

            case Add:
                matrix = std::move(node1->matrix + node2->matrix);
                break;
            case Sub:
                matrix = std::move(node1->matrix - node2->matrix);
                break;
            case Mult:
            // 逐元素乘法
            // z1 = x1 * y1
            // dz1 = x1.d * y1.v + x1.v * y1.d
                matrix = std::move(node1->matrix.multEach(node2->matrix));
                break;
            case MatrixMult:
                matrix = std::move(node->matrix * node2->matrix);
                break;
            case Devide:
                matrix = std::move(node1->matrix / node2->matrix);
                break;

            // case MatrixMult:
            case Transpose: // 理论上不应该有转置
                matrix = std::move(node1->Transpose());
                break;
            case Inverse:   // 理论上不应该有逆矩阵
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
                matrix = std::move(SumLine(node1->matrix));
                break;
            case MultAll:
                matrix = std::move(MultLine(node1->matrix));
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

    /**
     * 重要, 当前值和雅克比矩阵维数不匹配时, 直接拓展value位
     * 比如 y = multall(X), y 1维拓展成n维
     * ([y, dy/dx1], [y, dy/dx2], ..., [y, dy/dxn]) 变成 n维矩阵
     */
    // 反向模式(vjp)尚未实现; 当前梯度计算走 forward 模式(Forward 内的 Dual 传播)
    void Backward()
    {
        /*
        1. 函数矩阵, 一维当二维对角矩阵使用
        2. 线性运算乘法MatrixMult, 对X是W, 对W是X, 加法对B是1
        总结就是雅可比矩阵可以通过matrix里面的偏导以及inputs得到
        */
        
        std::shared_ptr<GradNode<T>> node1 = nullptr, node2 = nullptr;
        if (inputs.size() == 2) {
            node1 = inputs[0];
            node2 = inputs[1];
            
        } else if (inputs.size() == 1) {
            node1 = inputs[0];
        }
        // 每一个节点计算自己的梯度, 然后设置到变量
        // 前向是算d(data)/d(input), 后向是算d(loss)/d(data)
        // 因为loss单值, 所以雅可比矩阵一直只需要记录一列, 但是我们不做假设
        // 一直当正常matrix计算, 当然以后可以做一味的专属优化, 就是不走矩阵乘法
        switch (op) {
            // 数据节点不计算, 由上层帮自己设置
            case Input:
                break;
            case Output:
                // 最后一层不需要保存梯度, 就是1
                break;
            case Constant:
                break;
            case CrossEntropy:
                break;

            // node1, node2两个n维input
            // 所以要分别对node1和node2求雅可比矩阵
            // 矩阵只要用T.value, T.derive无关变量
            case Add:
                // node1 + node2
                // d(node1 + node2)/dnode1 = 1, d(node1 + node2)/dnode2 = 1
                node1->desent = desent;
                node2->desent = desent;
                break;
            case Sub:
                node1->desent = desent;
                node2->desent = Matrix(desent.Rows(), desent.Cols()) - desent;
                break;
            case Mult:
                // node1 * node2 (n1.v1 * n2.v1, n1.v2 * n2.v2, ...)
                // d(node1 * node2)/dnode1 = node2, d(node1 * node2)/dnode2 = node1
                // 雅可比矩阵是对角矩阵(n2.v1, n2.v2, n2.v3, ...)
                // 很明显其他位置都为0
                Matrix<T> jacobian1(node2->matrix.Rows(),
                                    node2->matrix.Cols());
                Matrix<T> jacobian2(node1->matrix.Rows(),
                                    node1->matrix.Cols());

                for (int i = 0; i < node1->matrix.Rows(); i++) {
                    jacobian1[i][i] = node2->matrix[0][i];
                    jacobian2[i][i] = node1->matrix[0][i];
                }
                node1->desent = jacobian1 * desent;
                node2->desent = jacobian2 * desent;
                
                break;
            case MatrixMult:
                if (node1->nodeType == Type::Parameter) {
                    auto tmp = node1;
                    node1 = node2;
                    node2 = tmp;
                }
                // X * W, 对X求偏微分是W, 对W求偏微分是X拓展n列
                node1->desent = node2->matrix * desent;
                node2->desent = Matrix<T>(node1->matrix.Rows(),
                                         node1->matrix.Cols());
                for (int i = 0; i < node1->matrix.Rows(); i++) {
                    for (int j = 0; j < node1->matrix.Cols(); j++) {
                        node2->desent[i][j] = node1->matrix[0][i] * desent[j][0];
                    }
                }
                break;
            case Devide:
                // node1 [m1, m2] / node2 [n1, n2]
                // d(node1/node2)/dnode1 = 1/n2, d(node1/node2)/dnode2 = -m1/n1^2
                Matrix<T> jacobian1(node2->matrix.Rows(),
                                    node2->matrix.Cols());
                Matrix<T> jacobian2(node1->matrix.Rows(),
                                    node1->matrix.Cols());
                for (int i = 0; i < node1->matrix.Rows(); i++) {
                    jacobian1[i][i] = 1.0f / node2->matrix[0][i];
                    jacobian2[i][i] = -node1->matrix[0][i] / (node2->matrix[0][i] * node2->matrix[0][i]);
                }
                node1->desent = jacobian1 * desent;
                node2->desent = jacobian2 * desent;
                break;

            // case MatrixMult:
            case Transpose: // 理论上不应该有转置
                matrix = std::move(node1->Transpose());
                break;
            case Inverse:   // 理论上不应该有逆矩阵
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
                matrix = std::move(SumLine(node1->matrix));
                break;
            case MultAll:
                matrix = std::move(MultLine(node1->matrix));
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
};
/* 使用行向量模式, 最优化内存排列和性能

1. x 8维 拆分成3维和5维, 也可以拓展更多维数
    | 1 0 0 |                               | x1  x1  x1 |
x   | 0 1 0 |  MatrixMultiply    求微分易得   | x2  x2  x2 |
    | 0 0 1 |                               | x3  x3  x3 |
这里我们将矩阵视为线性参数矩阵, yT = WT * xT

2. MatrixFunction 只保留一种, 只要是变量都应该是1维

行(x x^2 x^3) 应用在 3维 (x1, x2, x3) 上得到 1x3 新的 (x1, x2^2, x3^3)

如果需要更多的子函数, 就可以先拓展 (x1, x2, x3)
    | 1 0 0 1 0 0|
    | 0 1 0 0 1 1|   得到 (x1, x2, x3, x1, x2, x2)
    | 0 0 1 0 0 0|  
*/


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