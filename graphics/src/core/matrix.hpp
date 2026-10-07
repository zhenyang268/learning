#pragma once

#include <list>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>
#include "common.hpp"
#include "core/scalar.hpp"      // ScalarTraits: zero/one/conj/abs2
#include "autodiff/dual.hpp"
#include "core/scalar.hpp"

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

    // 带初始数据的构造 (校验: values.size() == m*n, matrix.md §2.2)
    Matrix(int m, int n, std::vector<T> vals)
        : rows(m), cols(n), data(std::move(vals))
    {
        if ((int)data.size() != m * n)
            throw std::invalid_argument("Matrix: values size != m*n");
    }

    // ---- 静态工厂 (matrix.md §2.2: 取代有歧义的 Matrix(m,n,initVal)) ----
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

    // 从initializer_list构造（带校验: 空表/行长一致, matrix.md §2.2）
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

    // 外部访问接口A[][], 语法糖 —— 即时派生行指针, 永不过期 (matrix.md §2.3)
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

    // 裸数据指针出口 (matrix.md §2.2; 成员名 data 被占用, 故大写)
    T*       Data()       { return data.data(); }
    const T* Data() const { return data.data(); }

    ~Matrix();

public:
    int Rows() const { return rows; }
    int Cols() const { return cols; }
    int Size() const { return rows * cols; }

    // **** operator 口径 (matrix.md §2.2): + - += 逐元素; 标量乘 operator*(T);
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

    // 矩阵乘 (operator* 唯一的矩阵语义, matrix.md §2.2)
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

    // 逐元素除法 (2026-10-07 补: autodiff Divide 反传的现役需求, matrix.md §2.2)
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