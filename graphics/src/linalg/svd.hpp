#pragma once

#include <stdexcept>
#include "core/common.hpp"
#include "core/matrix.hpp"
#include "core/scalar.hpp"

/**
 * SVD奇异值分解
 */

template <CommutativeScalar T>
struct SVDcomposition
{
    Matrix<T> U;
    Matrix<T> V;
    Matrix<T> S;
};

