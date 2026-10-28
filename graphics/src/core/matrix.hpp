#pragma once

// ===========================================================================
// core/matrix.hpp —— 矩阵（哑容器）
// ===========================================================================
//
// 【状态】已实现。
//
// ---------------------------------------------------------------------------
// 一、API 契约（调用方看这段就够）
// ---------------------------------------------------------------------------
//
// 形状: m×n, 扁平行主序存储 (元素 (i,j) 在 data[i*n + j])。
//
// 构造:
//   Matrix()                空 (0×0)
//   Matrix(m, n)            零初始化
//   Matrix(m, n, vector)    逐元素填充; **要求 vector.size() == m*n**, 否则 throw
//   Matrix{{...}, {...}}    initializer_list; 要求每行等长且非空, 否则 throw
//   静态工厂: Zero(m,n) / Ones(m,n) / Identity(n) / Diagonal(diag)
//             (取代有歧义的 Matrix(m,n,initVal) —— 分不清第二个 n 是列数还是填充值)
//
// 访问:
//   operator[](i)    返回**即时派生的行指针** (data + i*cols), 不缓存, 永不过期
//   operator()(i,j)  逐元素读写, 不查越界
//   at(i,j)          查越界, 抛 std::out_of_range
//   Data()           裸指针出口
//
// 运算 (白名单, 其余一律走 linalg/ 自由函数):
//   +  -  +=         逐元素; 形状不一致抛 std::invalid_argument
//   * (T)            标量乘
//   * (Matrix)       **矩阵乘** —— operator* 唯一的矩阵语义
//   cwiseProduct     逐元素乘 (Hadamard)
//   cwiseDivision    逐元素除
//   Transpose()      A^T
//   Adjoint()        A^H (共轭转置; 实域退化为 Transpose)
//
// 不变量:
//   · rows/cols 非负; data.size() == rows*cols (空矩阵两者皆 0)。
//   · 本类**不是**运算层: 求逆 / 分解 / 范数 / 元素级量度都不在这里。
// 域:
//   T 由 linalg 层的 `CommutativeScalar` concept 约束; 本类自身不设限,
//   因此 Matrix<Anything> 都能编过 —— 域约束在调用点, 不在容器。
//
// ---------------------------------------------------------------------------
// 二、调用顺序与内部机制（改实现时看）
// ---------------------------------------------------------------------------
//   · 存储: 单个 vector<T>, 无缓存层、无行指针数组 (见下文「为什么」)。
//   · 拷贝: 深拷贝 (vector 语义); 移动: 移动 vector 后把对方尺寸清零。
//   · multiply_common: 朴素 i-k-j 三重循环。i-k-j 而非 i-j-k 是因为
//     内层沿 j 走时访问的是 result[i][j] 与 B[k][j] —— 都是行内连续,
//     对 cache 友好。i-j-k 会每步跨行跳。
//
// ---------------------------------------------------------------------------
// 三、数学推导 / 踩坑记录（学习与考古）
// ---------------------------------------------------------------------------
// 【为什么是哑容器】
//   Matrix 只负责「存 + 取 + 逐元素四则」。求逆、QR、SVD、expm 全部是
//   linalg/ 下的自由函数。理由: 分解算法要能自由组合 (qr_inverse 用 QR 的结果,
//   但不该让 Matrix 知道 QR 存在), 而成员函数只能靠虚函数或 if 分派。
//   判据: 加一个新分解, 是否需要改 Matrix? 需要 = 设计错了。
//
// 【为什么不存行指针 vector<T*>】
//   行指针是 data 的函数。任何一次 resize / 移动 / 赋值都会让它过期,
//   而“过期”不会报错 —— 只会静默读错内存。改成 operator[](i) 即时派生,
//   代价是每次一次乘法, 收益是**不存在过期的可能**。
//
// 【为什么不存运行期结构标签】
//   operator()(i,j) 返回 T& 是公开出口 —— 所有分解都靠它写元素。只要存在
//   这条通路, 任何缓存的结构标签 (is_symmetric / is_diagonal) 都能被外部
//   改动悄悄作废。结构判断必须每次现算 (见 linalg/inverse.hpp 的 is_* 函数)。
//
// 【容差为什么不进这个类】
//   容差是**算法的参数**, 不是容器的属性。同一个矩阵在不同问题里可能是
//   “良态”或“病态”, 取决于右端项与尺度。写死在 Matrix 里等于替调用方
//   做了一个它不知道的决定。见 numeric-contract.md §3–§4。
//
// 【A^H 而不是 A^T —— 踩过的坑】
//   Complex 域下 Transpose 只转置不共轭, 用它替代 Adjoint 会让所有酉相关的
//   判据 (Q^H·Q = I、Hermitian 判定) 静默失效。这个错误在实域测试里
//   完全看不见 (double 下 conj 恒等)。已两次踩: givens 旋转与 qr_inverse。
//   ⇒ 写任何 QR / 酉相关代码, **必查此处用的是 Transpose 还是 Adjoint**。
// ===========================================================================

#include <vector>
#include <utility>
#include <stdexcept>
#include "autodiff/dual.hpp"

template <typename T>
class Matrix
{
private:
    int rows;
    int cols;
    
    std::vector<T> data;

public:
    Matrix(): rows(0), cols(0) {}

    Matrix(int m, int n) : rows(m), cols(n), data(m * n)
    {
    }

    // 带初始数据的构造 (校验: values.size() == m*n; 见文件头 API 契约)
    Matrix(int m, int n, std::vector<T> vals)
        : rows(m), cols(n), data(std::move(vals))
    {
        if ((int)data.size() != m * n)
            throw std::invalid_argument("Matrix: values size != m*n");
    }

    // ---- 静态工厂 (取代有歧义的 Matrix(m,n,initVal); 见文件头 API 契约) ----
    static Matrix Zero(int m, int n) { return Matrix(m, n); }

    static Matrix Ones(int m, int n)
    {
        Matrix r(m, n);
        T one = ScalarTraits<T>::one();
        for (int i = 0; i < m * n; i++) r.data[i] = one;
        return r;
    }

    static Matrix Identity(int n)
    {
        Matrix r = Zero(n, n);
        T one = ScalarTraits<T>::one();
        for (int i = 0; i < n; i++) r.data[i * n + i] = one;
        return r;
    }

    static Matrix Diagonal(std::vector<T> diag)
    {
        int n = (int)diag.size();
        Matrix r = Zero(n, n);
        for (int i = 0; i < n; i++) r.data[i * n + i] = std::move(diag[i]);
        return r;
    }

    // 从 initializer_list 构造 (带校验: 空表 / 行长一致; 见文件头 API 契约)
    Matrix(std::initializer_list<std::initializer_list<T>> list) {
        rows = (int)list.size();
        if (rows == 0)
            throw std::invalid_argument("Matrix: empty initializer_list");
        cols = (int)(list.begin())->size();

        data = std::vector<T>(rows * cols);

        int i = 0;
        for (const auto& row : list) {
            if ((int)row.size() != cols)
                throw std::invalid_argument("Matrix: inconsistent row length");
            int j = 0;
            for (T val : row) {
                data[cols * i + j++] = val;
            }
            i++;
        }
    }

    // ---------- Rule of Three：深拷贝 ----------
    // 注意: matrix 是 data 的行指针数组, 拷贝时必须重建指针, 否则会指向原对象 data
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols),
        data(other.data)
    {
    }

    Matrix& operator=(const Matrix& other)
    {
        if (this == &other) return *this;

        rows = other.rows;
        cols = other.cols;
        data = other.data;
        return *this;
    }

    // ---------- 移动构造函数获取所有权 ---------- // 防止重复delete
    // 注意: noexcept 必须位于 mem-initializer 列表之前
    Matrix(Matrix&& other) noexcept : rows(other.rows), cols(other.cols),
        data(std::move(other.data))
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
        data = std::move(other.data);
        
        other.rows = 0;
        other.cols = 0;
        return *this;
    }

    // 外部访问接口 A[][]: 即时派生行指针, 永不过期 (理由见文件末尾「为什么不存行指针」)
    T* operator[](int row) { return data.data() + row * cols; }
    const T* operator[](int row) const { return data.data() + row * cols; }

    // 二维访问: A(i, j)
    T&       operator()(int i, int j)       { return data[i * cols + j]; }
    const T& operator()(int i, int j) const { return data[i * cols + j]; }

    // 安全版: at(i, j) —— 检查后直接算
    T& at(int i, int j) {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Matrix::at");
        return data[i * cols + j];
    }
    const T& at(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Matrix::at");
        return data[i * cols + j];
    }

    // 裸数据指针出口 (成员名 data 被占用, 故大写)
    T*       Data()       { return data.data(); }
    const T* Data() const { return data.data(); }

    ~Matrix();

public:
    int Rows() const { return rows; }
    int Cols() const { return cols; }
    int Size() const { return rows * cols; }

    // **** operator 口径 (见文件头 API 契约): + - += 逐元素; 标量乘 operator*(T);
    // 逐元素乘/除显式命名 cwiseProduct / cwiseDivision; operator* 只给矩阵乘 ----
    Matrix operator+(const Matrix& other) const
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

    Matrix operator-(const Matrix& other) const
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

    // 标量乘
    Matrix operator*(const T& s) const
    {
        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] * s;
        }
        return result;
    }

    // 矩阵乘 (operator* 唯一的矩阵语义)
    Matrix operator*(const Matrix& other) const
    {
        return multiply_common(*this, other);
    }

    // 逐元素乘法, Hadamard product
    Matrix cwiseProduct(const Matrix& other) const
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

    // 逐元素除法 (autodiff Divide 反传的现役需求, 2026-10-07 补)
    Matrix cwiseDivision(const Matrix& other) const
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

    Matrix Transpose() const;
    Matrix Adjoint() const;


private:
    static Matrix multiply_common(const Matrix &A, const Matrix &B);

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

template <typename T>
Matrix<T> Matrix<T>::Transpose() const
{
    Matrix<T> result(this->cols, this->rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = data[i * cols + j];
        }
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::Adjoint() const
{
    Matrix<T> result(this->cols, this->rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = ScalarTraits<T>::conj(data[i * cols + j]);
        }
    }
    return result;
}

using Matrixf = Matrix<float>;
using Matrixd = Matrix<double>;
using MatrixDual = Matrix<Dual<float>>;