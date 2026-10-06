#pragma once

#include <list>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>
#include "common.hpp"
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

    Matrix(int m, int n, T initVal): rows(m), cols(n), data(m * n, initVal)
    {
    }

    Matrix(int m, int n) : rows(m), cols(n), data(m * n)
    {
    }

    // 从initializer_list构造（方便初始化）
    Matrix(std::initializer_list<std::initializer_list<T>> list) {
        rows = list.size();
        cols = (list.begin())->size();
        
        data = std::vector<T>(rows * cols);
        std::vector<T*> matrix(rows);
        
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

    // 外部访问接口A[][], 语法糖
    T* operator[](int row) { return &data[row * cols]; }
    const T* operator[](int row) const { return &data[row * cols]; }

    // 二维访问: A(i, j)
    T&       operator()(int i, int j)       { return data[i * cols_ + j]; }
    const T& operator()(int i, int j) const { return data[i * cols_ + j]; }

    // 安全版: at(i, j) —— 检查后直接算, 不经过任何语法糖
    // 不靠谱的调用只能通过这个接口
    T& at(int i, int j) {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
            throw std::out_of_range("Matrix::at");
        return data[i * cols_ + j];        // ← 内部直算
    }

    vector<T> GetData() { return data; }

    ~Matrix();
    
public:
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
    
    // 逐元素乘法, Hadamard product
    Matrix CWiseProduct(const Matrix& other)
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

    // 逐元素乘法, 乘上统一的标量
    Matrix CWiseProduct(const T& other)
    {
        Matrix result(*this);
        for (int i = 0; i < rows * cols; i++) {
            result.data[i] = result.data[i] * other;
        }
        return result;
    }


    Matrix Transpose() const;
    Matrix Adjoint() const;
    

private:
    Matrix multiply_common(const Matrix &A, const Matrix &B);
    Matrix multiply_tiling(const Matrix &A, const Matrix &B);  

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