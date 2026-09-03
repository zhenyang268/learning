#pragma once


#include <list>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>
#include "common.hpp"
#include "autodiff/dual.hpp"

/*
存储matrix, float**, m * n
    返回transpose转置
    求解特征值和特征向量(矩阵)
    SVD分解
    返回inverse求逆
        矩阵类型	求逆方法	时间复杂度	特点
        旋转矩阵	转置	O(n²)	正交矩阵特性
        缩放矩阵	对角元倒数	O(n)	对角矩阵
        平移矩阵	向量取反	O(1)	特殊结构
        刚体变换	分块公式	O(3³)	利用结构
        投影矩阵	解析公式	O(1)	已知结构
        一般矩阵	LU分解	O(n³)	最常用
        病态矩阵	SVD	O(n³)	最稳定
        对称正定	Cholesky	O(n³/3)	最快通用法

旋转矩阵：转置求逆
旋转矩阵是标准正交矩阵，满足 \(R^T R = I\)，因此它的逆矩阵就等于自身的转置 \(R^{-1}=R^T\)。
这是图形学中最高频的求逆操作，比如坐标系转换、相机视角反转，性能远高于通用求逆。
万向角问题:XYZ
就是说绕x旋转, 绕y 90度, 再绕z旋转等价于, 一开始绕x旋转或者最后绕z旋转
x轴旋转后和原先的z轴重合, 旋转是依赖于原先的坐标系的
理解就是每个旋转后的坐标都是原坐标系下的

缩放矩阵：对角元倒数
纯缩放矩阵是对角矩阵，对角线上为各轴缩放系数。对角矩阵的逆只需将每个对角元素取倒数，无需做任何矩阵运算。
注意：若某轴缩放系数为 0，则矩阵不可逆。

平移矩阵：向量取反
齐次坐标下的平移矩阵有固定结构（左上角为单位矩阵，最后一列是平移向量）。它的逆矩阵只需把平移向量取反，其余元素完全不变，计算量为常数级。
常用于模型局部坐标与世界坐标的互相转换。


刚体变换：分块公式
刚体变换（仅包含旋转 + 平移，无缩放 / 切变）的 4×4 齐次矩阵，可以拆分为 3×3 旋转分块和 3 维平移分块。
利用旋转矩阵转置求逆的特性，可直接推导出分块逆公式：
若变换矩阵 \(T=\begin{bmatrix}R & t \\ 0 & 1\end{bmatrix}\)，则 \(T^{-1}=\begin{bmatrix}R^T & -R^T t \\ 0 & 1\end{bmatrix}\)。
性能远优于完整的 4 阶矩阵通用求逆，视图矩阵、模型矩阵的求逆都适用。


投影矩阵：解析公式
正交投影、透视投影矩阵都有固定的数学形式，可以直接推导出逆矩阵的解析表达式，代入数值即可得到结果，无需数值分解。
典型应用：屏幕坐标反推世界空间射线、鼠标拾取算法。


一般矩阵：LU 分解
将矩阵分解为下三角矩阵 L 和上三角矩阵 U，再通过三角矩阵快速求解逆矩阵，是通用方阵求逆的工业标准。
无特殊结构的矩阵默认使用该方法，也是 Eigen 等线性代数库的默认实现。


病态矩阵：SVD 分解
病态矩阵指条件数很大、接近奇异（不可逆）的矩阵，LU 分解会出现数值不稳定、除以零等问题。
SVD 奇异值分解对病态矩阵鲁棒性最强，可以稳定求解伪逆，是稳定性优先场景的首选，但计算量比 LU 分解更大。
常用于点云配准、最小二乘拟合、数值优化。


对称正定矩阵：Cholesky 分解
对称正定矩阵可分解为下三角矩阵与其转置的乘积 \(LL^T\)，分解计算量约为 LU 分解的 1/3，求逆速度显著更快。
常见于物理仿真的质量矩阵、协方差矩阵、法方程求解场景。

拓展四元数转化
拓展slerp插值
    
*/

template <typename T>
class Matrix
{
private:
    int rows;
    int cols;
    
    std::vector<T> data;
    std::vector<T*> matrix; // 行标识

    int _errno;

public:
    Matrix(): rows(0), cols(0), data(0), matrix(0) {}

    // 初始化对角矩阵
    Matrix(int m, T initVal) : rows(m), cols(m), data(m * m), matrix(m)
    {
        for (int i = 0; i < rows; i++) {
            matrix[i] = &data[cols * i];
        }
        for (int i = 0; i < rows; i++) {
            matrix[i][i] = initVal;
        }
    }

    Matrix(int m, int n): rows(m), cols(n), data(m * n), matrix(m)
    {
        // Quote: 二维数组[][]必须通过赋值方式额外存储
        for (int i = 0; i < rows; i++) {
            matrix[i] = &data[cols * i];
        }
    }

    Matrix(int m, int n, std::vector<T> values): rows(m), cols(n), data(values), matrix(m)
    {
        // Quote: 二维数组[][]必须通过赋值方式额外存储
        
        for (int i = 0; i < rows; i++) {
            matrix[i] = &data[cols * i]; 
        }
        // std::copy(values.first, values.first + m * n, data.first);
    }

    // 从initializer_list构造（方便初始化）
    Matrix(std::initializer_list<std::initializer_list<T>> list) {
        rows = list.size();
        cols = (list.begin())->size();
        
        data = std::vector<T>(rows * cols);
        matrix = std::vector<T*>(rows);
        
        int i = 0;
        for (const auto& row : list) {
            matrix[i] = &data[cols * i];
            int j = 0;
            for (T val : row) {
                matrix[i][j++] = val;
            }
            i++;
        }
    }

    // ---------- Rule of Three：深拷贝 ----------
    // 注意: matrix 是 data 的行指针数组, 拷贝时必须重建指针, 否则会指向原对象 data
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols),
        data(other.data), matrix(rows)
    {
        for (int i = 0; i < rows; i++) {
            matrix[i] = &data[cols * i];
        }
    }

    Matrix& operator=(const Matrix& other)
    {
        if (this == &other) return *this;

        rows = other.rows;
        cols = other.cols;
        data = other.data;
        matrix.resize(rows);
        for (int i = 0; i < rows; i++) {
            matrix[i] = &data[cols * i];
        }
        return *this;
    }

    // ---------- 移动构造函数获取所有权 ---------- // 防止重复delete
    // 注意: noexcept 必须位于 mem-initializer 列表之前
    Matrix(Matrix&& other) noexcept : rows(other.rows), cols(other.cols),
        matrix(std::move(other.matrix)), data(std::move(other.data))
    {
        other.rows = 0;
        other.cols = 0;
    }

    // ---------- 移动赋值函数获取所有权,并且清空自己资源 ---------- // 防止重复delete
    Matrix& operator=(Matrix&& other) noexcept
    {
        if (this == &other) return *this;

        rows = other.rows;
        cols = other.cols;
        matrix = std::move(other.matrix);
        data = std::move(other.data);
        
        other.rows = 0;
        other.cols = 0;
        return *this;
    }

    // 外部访问接口A[][]
    T* operator[](int row) { return matrix[row]; }
    const T* operator[](int row) const { return matrix[row]; }

    ~Matrix();
    
private:
    // 加速用缓存, 关键点是矩阵如果是变化的, 就不能用缓存
    bool couldCached = false;
    mutable std::unique_ptr<Matrix> cachedTranspose; // 转置
    mutable std::unique_ptr<Matrix> cachedInverse;  // 逆矩阵

    mutable std::list<std::unique_ptr<Matrix>> decomposition; // lu和矩阵乘法使用

public:
    // 如果矩阵stable初始化后不变, 就可以使用缓存
    void EnableCache() { couldCached = true; }
    int Rows() { return rows; }
    int Cols() { return cols; }
    int Size() { return rows * cols; }

    // **** 重点, operator都是逐元素运算符号, 乘法用multiply()
    // operate: transpose, svd, inverse, lu, etc.
    // operate: + += *
    Matrix operator+(const Matrix& other)
    {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] + other.data[i];
        }
        return result;
    }

    Matrix operator-(const Matrix& other)
    {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] - other.data[i];
        }
        return result;
    }

    Matrix operator/(const Matrix& other)
    {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] / other.data[i];
        }
        return result;
    }

    Matrix& operator+=(const Matrix& other)
    {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        for (int i = 0; i < rows * cols; i++) {
            this->data[i] = this->data[i] + other.data[i];
        }
        return *this;
    }

    Matrix operator*(const Matrix& other)
    {
        return multiply_common(*this, other);
    }
    
    // 逐元素乘法
    Matrix MultEach(const Matrix& other)
    {
        if (this->rows != other.rows || this->cols != other.cols) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] * other.data[i];
        }
        return result;
    }

    // 计算雅可比矩阵, 但是用n维存储对角线版
    Matrix MultDiagnal(const Matrix& other)
    {
        // 特殊运算规则, 函数雅可比矩阵运算必须n*n维
        if (this->rows != other.rows || this->rows != other.rows) {
            throw std::invalid_argument("Matrix dimension mismatch!");
        }

        Matrix result(this->rows, this->cols);
        for (int i = 0; i < this->rows; i++) {
            result[i][i] = this->matrix[0][i] * other.matrix[i][i];
        }
        return result;
    }

    Matrix Transpose() const;
    Matrix Inverse() const;

    friend Matrix<T> SumLine(Matrix<T> &A); // 累加
    // 仅限Dual类用的微分矩阵    
    friend Matrix<T> MultLine(Matrix<T> &A);  // 累乘,  通过e^ln(x1*....*xn)计算
    friend Matrix<T> Apply(Matrix<T> &A, DualFunc<T> f); // 对矩阵所有元素应用f
    friend Matrix<T> Apply(Matrix<T> &A, Matrix<DualFunc<T>> &matrixf); // 逐元素应用f

private:
    Matrix multiply_common(const Matrix &A, const Matrix &B);
    Matrix multiply_tiling(const Matrix &A, const Matrix &B);

    bool is_orthogonal();           // 正交
    bool is_symmetric();            // 对称
    bool is_positive_definite();    // 正定

    Matrix lu_inverse();            // 求解梯度矩阵lu
    Matrix diagnal_inverse();       // 对角矩阵
    Matrix cholesky_inverse();      // 对称正定 LLT, 下三角
    Matrix qr_inverse();            // QR分解, 正交*上三角

};

// ---------- 模板定义必须放在头文件, 否则跨 TU 无法实例化 ----------

template <typename T>
Matrix<T>::~Matrix()
{
    // std::cout << "Matrix destructor called" << std::endl;
    // matrix和data是vector, 会自动释放内存
}

// 标准矩阵乘法: A(m×k) * B(k×n) = C(m×n)
template <typename T>
Matrix<T> Matrix<T>::multiply_common(const Matrix &A, const Matrix &B)
{
    if (A.cols != B.rows) {
        throw std::invalid_argument("Matrix dimension mismatch!");
    }

    Matrix<T> result(A.rows, B.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int k = 0; k < A.cols; k++) {
            for (int j = 0; j < B.cols; j++) {
                result[i][j] = result[i][j] + A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// 对矩阵所有元素应用同一函数
template <typename T>
Matrix<T> Apply(Matrix<T> &A, DualFunc<T> f)
{
    Matrix<T> B(A.rows, A.cols);
    for (int i = 0; i < A.rows * A.cols; i++) {
        B.data[i] = f(A.data[i]);
    }

    return B;
}

// 逐元素应用函数矩阵
template <typename T>
Matrix<T> Apply(Matrix<T> &A, Matrix<DualFunc<T>> &matrixf)
{
    if (matrixf.rows == 1 && matrixf.cols == 1) {
        return Apply(A, matrixf[0][0]); // 单函数退化为统一应用
    }

    if (A.rows != 1 || matrixf.rows != 1 || A.cols != matrixf.cols) {
        throw std::invalid_argument("Matrix apply dimension mismatch!");
    }

    Matrix<T> B(A.rows, A.cols);
    for (int i = 0; i < A.rows * A.cols; i++) {
        B.data[i] = matrixf.data[i](A.data[i]);
    }

    return B;
}

// 累加: 对第一行做 sumLine, 输出与 A 同形
template <typename T>
Matrix<T> SumLine(Matrix<T> &A)
{
    Matrix<T> B = A;
    sumLine(B[0], B.Size());
    return B;
}

// 累乘
template <typename T>
Matrix<T> MultLine(Matrix<T> &A)
{
    Matrix<T> B = A;
    multLine(B[0], B.Size());
    return B;
}

template <typename T>
Matrix<T> Matrix<T>::Transpose() const
{
    Matrix<T> result(this->cols, this->rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = this->matrix[i][j];
        }
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::Inverse() const
{

}

using Matrixf = Matrix<float>;
using Matrixd = Matrix<double>;
using MatrixDual = Matrix<Dual<float>>;