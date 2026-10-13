#pragma once

#include <cmath>
#include <vector>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"
#include "core/complex.hpp"

template <CommutativeScalar T>
T Dot(const Matrix<T>& other)
{
    if (this->rows != other.rows || this->cols != other.cols) {
        throw std::invalid_argument("Matrix dimension mismatch!");
    }
    T result = 0;
    for (int i = 0; i < rows * cols; i++) {
        result += this->data[i] * other.data[i];
    }
    return result;
}

template <CommutativeScalar T>
T Abs2(Matrix<T>& matrix)
{
    T result = 0;
    for (int i = 0; i < matrix.Rows() * matrix.Cols(); i++) {
        result += ScalarTraits<T>::abs2(matrix.Data[i]);
    }
    return result;
}