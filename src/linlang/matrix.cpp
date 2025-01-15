#include "matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <initializer_list>
#include <limits>
#include <new>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace math {
namespace linlang {

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, T(0)) {
  if (rows == 0 || cols == 0)
    throw std::invalid_argument("Matrix dimensions must be positive");
}

template <typename T>
Matrix<T>::Matrix(const Matrix &other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init) {
  rows_ = init.size();
  if (rows_ == 0)
    throw std::invalid_argument("Empty initializer list");

  cols_ = init.begin()->size();
  data_.reserve(rows_ * cols_);
  for (const auto &row : init) {
    if (row.size() != cols_)
      throw std::invalid_argument("Irregular matrix initialization");
    data_.insert(data_.end(), row.begin(), row.end());
  }
}

template <typename T> Matrix<T> Matrix<T>::Indentify(size_t size) {
  Matrix<T> result(size, size);
  for (size_t i = 0; i < size; i++)
    result(i, i) = T(1);
  return result;
}

template <typename T> Matrix<T> Matrix<T>::Rotation2D(T angle) {
  Matrix<T> R(2, 2);
  T cos_theta = std::cos(angle);
  T sin_theta = std::sin(angle);
  R(0, 0) = cos_theta;
  R(0, 1) = -sin_theta;
  R(1, 0) = sin_theta;
  R(1, 1) = -cos_theta;

  return R;
}

template <typename T>
Matrix<T> Matrix<T>::Rotation3D(T angle, const std::vector<T> &axis) {
  Matrix<T> R(3, 3);
  T cos_theta = std::cos(angle);
  T sin_theta = std::sin(angle);
  T one_minus_cos = T(1) - cos_theta;

  std::vector<T> a = axis.normalized();
  R(0, 0) = cos_theta + a.x() * a.x() * one_minus_cos;
  R(0, 1) = a.x() * a.y() * one_minus_cos - a.z() * sin_theta;
  R(0, 2) = a.x() * a.z() * one_minus_cos + a.y() * sin_theta;

  R(1, 0) = a.y() * a.x() * one_minus_cos + a.z() * sin_theta;
  R(1, 1) = cos_theta + a.y() * a.y() * one_minus_cos;
  R(1, 2) = a.y() * a.z() * one_minus_cos - a.x() * sin_theta;

  R(2, 0) = a.z() * a.x() * one_minus_cos - a.y() * sin_theta;
  R(2, 1) = a.z() * a.y() * one_minus_cos + a.x() * sin_theta;
  R(2, 2) = cos_theta + a.z() * a.z() * one_minus_cos;

  return R;
}

// Basic operations
template <typename T> Matrix<T> &Matrix<T>::operator=(const Matrix &other) {
  if (this != &other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    data_ = other.data_;
  }
  return *this;
}

template <typename T> Matrix<T> &Matrix<T>::operator=(Matrix &other) noexcept {
  if (this != &other) {
    rows_ = other.rows_;
    data_ = std::move(other.data_);
    other.rows_ = 0;
    other.cols_ = 0;
  }
  return *this;
}

template <typename T> T &Matrix<T>::operator()(size_t row, size_t col) {
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("Matrix indices out of bounds");
  return data_[index(row, col)];
}

template <typename T>
const T &Matrix<T>::operator()(size_t row, size_t col) const {
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("Matrix indices out of range");
  return data_[index(row, col)];
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) const {
  checkDimension(other);
  Matrix<T> result(rows_, cols_);

#ifdef __AVX__
  for (size_t i = 0; i < data_.size(); i += 4) {
    __m256d a = _mm256_load_pd(&data_[i]);
    __m256d b = _mm256_load_pd(&other.data_[i]);
    __m256d sum = _mm256_add_pd(a, b);
    _mm256_store_pd(&result.data_[i], sum);
  }
#else
  std::transform(data_.begin(), data_.end(), other.data_.begin(),
                 result.data_.begin(), std::plus<T>());
#endif
  return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) const {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("Invalid matrix dimensions for multiplication");
  }

  Matrix<T> result(rows_, other.cols_);

  constexpr size_t BLOCK_SIZE = 32;
  for (size_t i = 0; i < rows_; i += BLOCK_SIZE) {
    for (size_t j = 0; j < other.cols_; j += BLOCK_SIZE) {
      for (size_t k = 0; k < cols_; k += BLOCK_SIZE) {
        for (size_t bi = i; bi < std::min(i + BLOCK_SIZE, rows_); ++bi) {
          for (size_t bj = j; bj < std::min(j + BLOCK_SIZE, other.cols_);
               ++bj) {
            T sum = T(0);
            for (size_t bk = k; bk < std::min(k + BLOCK_SIZE, cols_); ++bk) {
              sum += (*this)(bi, bk) * other(bk, bj);
            }
            result(bi, bj) += sum;
          }
        }
      }
    }
  }
  return result;
}

template <typename T> Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> result(cols_, rows_);

  constexpr size_t BLOCK_SIZE = 32;
  for (size_t i = 0; i < rows_; i += BLOCK_SIZE) {
    for (size_t j = 0; j < cols_; j += BLOCK_SIZE) {
      for (size_t bi = i; bi < std::min(i + BLOCK_SIZE, rows_); ++bi) {
        for (size_t bj = j; bj < std::min(j + BLOCK_SIZE, cols_); ++bj)
          result(bj, bi) = (*this)(bi, bj);
      }
    }
  }
  return result;
}

template <typename T> T Matrix<T>::determinant() const {
  if (!isSquare())
    throw std::invalid_argument("Determinant only defined for square matrices");

  if (rows_ == 1)
    return data_[0];
  if (rows_ == 2) {
    return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
  }
  auto [I, U] = lu();
  T det = T(1);
  for (size_t i = 0; i < rows_; ++i)
    det *= U(i, i);

  return det;
}

template <typename T>
void Matrix<T>::checkDimension(const Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Matrix dimensions must match");
}

template <typename T> size_t Matrix<T>::index(size_t row, size_t col) const {
  return row * cols_ + col;
}

// Continuing from previous implementation...

template <typename T> Matrix<T> Matrix<T>::inverse() const {
  if (!isSquare()) {
    throw std::invalid_argument("Matrix must be square to compute inverse");
  }

  // For 2x2 matrices, use analytical formula for better performance
  if (rows_ == 2) {
    T det = determinant();
    if (std::abs(det) < std::numeric_limits<T>::epsilon()) {
      throw std::runtime_error("Matrix is singular, cannot compute inverse");
    }

    Matrix<T> result(2, 2);
    result(0, 0) = (*this)(1, 1) / det;
    result(0, 1) = -(*this)(0, 1) / det;
    result(1, 0) = -(*this)(1, 0) / det;
    result(1, 1) = (*this)(0, 0) / det;
    return result;
  }

  // For larger matrices, use LU decomposition with partial pivoting
  size_t n = rows_;
  Matrix<T> result = Identity(n);
  auto [L, U, P] = luDecompositionWithPivoting();

  // Solve LY = P for Y
  for (size_t j = 0; j < n; ++j) {
    for (size_t i = 0; i < n; ++i) {
      T sum = result(P[i], j);
      for (size_t k = 0; k < i; ++k) {
        sum -= L(i, k) * result(P[k], j);
      }
      result(P[i], j) = sum;
    }
  }

  // Solve UX = Y for X
  for (size_t j = 0; j < n; ++j) {
    for (int i = n - 1; i >= 0; --i) {
      T sum = result(P[i], j);
      for (size_t k = i + 1; k < n; ++k) {
        sum -= U(i, k) * result(P[k], j);
      }
      if (std::abs(U(i, i)) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("Matrix is singular, cannot compute inverse");
      }
      result(P[i], j) = sum / U(i, i);
    }
  }

  return result;
}

template <typename T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>>
Matrix<T>::luDecompositionWithPivoting() const {
  if (!isSquare()) {
    throw std::invalid_argument("LU decomposition requires square matrix");
  }

  size_t n = rows_;
  Matrix<T> L = Identity(n);
  Matrix<T> U = *this;
  std::vector<size_t> P(n);
  for (size_t i = 0; i < n; ++i)
    P[i] = i;

  for (size_t k = 0; k < n - 1; ++k) {
    // Find pivot
    size_t pivot_row = k;
    T max_val = std::abs(U(k, k));

    for (size_t i = k + 1; i < n; ++i) {
      if (std::abs(U(i, k)) > max_val) {
        max_val = std::abs(U(i, k));
        pivot_row = i;
      }
    }

    // Swap rows if necessary
    if (pivot_row != k) {
      U.swapRows(k, pivot_row);
      L.swapRows(k, pivot_row);
      std::swap(P[k], P[pivot_row]);
    }

    // Perform elimination
    for (size_t i = k + 1; i < n; ++i) {
      T factor = U(i, k) / U(k, k);
      L(i, k) = factor;

      for (size_t j = k; j < n; ++j) {
        U(i, j) -= factor * U(k, j);
      }
    }
  }

  return {L, U, P};
}

template <typename T> Matrix<T> Matrix<T>::pseudoInverse() const {
  // Compute SVD decomposition
  auto [U, S, V] = svd();

  // Create reciprocal of singular values with threshold
  Matrix<T> S_inv = Matrix<T>::Zeros(cols_, rows_);
  T threshold =
      std::numeric_limits<T>::epsilon() * std::max(rows_, cols_) * S(0, 0);

  for (size_t i = 0; i < std::min(rows_, cols_); ++i) {
    if (S(i, i) > threshold) {
      S_inv(i, i) = T(1) / S(i, i);
    }
  }

  // Compute pseudoinverse as V * S^+ * U^T
  return V * S_inv * U.transpose();
}

template <typename T> std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::qr() const {
  size_t m = rows_;
  size_t n = cols_;

  Matrix<T> Q = Identity(m);
  Matrix<T> R = *this;

  // Householder QR decomposition
  for (size_t k = 0; k < std::min(m - 1, n); ++k) {
    // Construct Householder vector
    std::vector<T> x = R.getColumn(k).subvector(k, m - 1);
    T norm_x = x.norm();

    if (norm_x > std::numeric_limits<T>::epsilon()) {
      if (x(0) > 0)
        norm_x = -norm_x;

      x(0) -= norm_x;
      T norm_u = x.norm();

      if (norm_u > std::numeric_limits<T>::epsilon()) {
        x *= (T(1) / norm_u);

        // Apply Householder reflection to R and Q
        for (size_t j = k; j < n; ++j) {
          T sum = T(0);
          for (size_t i = k; i < m; ++i) {
            sum += R(i, j) * x(i - k);
          }
          sum *= T(2);

          for (size_t i = k; i < m; ++i) {
            R(i, j) -= sum * x(i - k);
          }
        }

        for (size_t j = 0; j < m; ++j) {
          T sum = T(0);
          for (size_t i = k; i < m; ++i) {
            sum += Q(j, i) * x(i - k);
          }
          sum *= T(2);

          for (size_t i = k; i < m; ++i) {
            Q(j, i) -= sum * x(i - k);
          }
        }
      }
    }
  }

  return {Q, R};
}

template <typename T> Matrix<T> Matrix<T>::cholesky() const {
  if (!isSquare() || !isSymmetric() || !isPositiveDefinite()) {
    throw std::invalid_argument("Matrix must be symmetric positive definite "
                                "for Cholesky decomposition");
  }

  size_t n = rows_;
  Matrix<T> L = Zeros(n, n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      T sum = 0;

      if (j == i) {
        for (size_t k = 0; k < j; ++k) {
          sum += L(j, k) * L(j, k);
        }
        T val = (*this)(j, j) - sum;
        if (val <= 0) {
          throw std::runtime_error("Matrix is not positive definite");
        }
        L(j, j) = std::sqrt(val);
      } else {
        for (size_t k = 0; k < j; ++k) {
          sum += L(i, k) * L(j, k);
        }
        L(i, j) = ((*this)(i, j) - sum) / L(j, j);
      }
    }
  }

  return L;
}

template <typename T> std::vector<T> Matrix<T>::solve(const std::vector<T> &b) const {
  if (rows_ != b.size()) {
    throw std::invalid_argument("Invalid dimensions for solving linear system");
  }

  if (isSquare()) {
    auto [L, U, P] = luDecompositionWithPivoting();

    // Solve Ly = Pb
    std::vector<T> y(rows_);
    for (size_t i = 0; i < rows_; ++i) {
      T sum = b[P[i]];
      for (size_t j = 0; j < i; ++j) {
        sum -= L(i, j) * y[j];
      }
      y[i] = sum;
    }

    // Solve Ux = y
    std::vector<T> x(rows_);
    for (int i = rows_ - 1; i >= 0; --i) {
      T sum = y[i];
      for (size_t j = i + 1; j < rows_; ++j) {
        sum -= U(i, j) * x[j];
      }
      x[i] = sum / U(i, i);
    }

    return x;
  } else {
    auto [Q, R] = qr();
    return R.solve(Q.transpose() * b);
  }
}

template <typename T> T Matrix<T>::norm(const std::string &type) const {
  if (type == "frobenius") {
    T sum = T(0);
    for (const auto &val : data_) {
      sum += val * val;
    }
    return std::sqrt(sum);
  } else if (type == "1") {
    // Column sum norm
    T max_sum = T(0);
    for (size_t j = 0; j < cols_; ++j) {
      T sum = T(0);
      for (size_t i = 0; i < rows_; ++i) {
        sum += std::abs((*this)(i, j));
      }
      max_sum = std::max(max_sum, sum);
    }
    return max_sum;
  } else if (type == "infinity") {
    // Row sum norm
    T max_sum = T(0);
    for (size_t i = 0; i < rows_; ++i) {
      T sum = T(0);
      for (size_t j = 0; j < cols_; ++j) {
        sum += std::abs((*this)(i, j));
      }
      max_sum = std::max(max_sum, sum);
    }
    return max_sum;
  } else {
    throw std::invalid_argument("Unsupported norm type");
  }
}

template <typename T> T Matrix<T>::condition() const {
  auto [U, S, V] = svd();
  T max_singular = S(0, 0);
  T min_singular = S(0, 0);

  for (size_t i = 1; i < std::min(rows_, cols_); ++i) {
    max_singular = std::max(max_singular, S(i, i));
    min_singular = std::min(min_singular, S(i, i));
  }

  if (min_singular < std::numeric_limits<T>::epsilon()) {
    return std::numeric_limits<T>::infinity();
  }

  return max_singular / min_singular;
}

template <typename T> bool Matrix<T>::isSymmetric() const {
  if (!isSquare())
    return false;

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = i + 1; j < cols_; ++j) {
      if (std::abs((*this)(i, j) - (*this)(j, i)) >
          std::numeric_limits<T>::epsilon()) {
        return false;
      }
    }
  }
  return true;
}

template <typename T> bool Matrix<T>::isOrthogonal() const {
  if (!isSquare())
    return false;

  Matrix<T> product = transpose() * (*this);
  Matrix<T> identity = Identity(rows_);

  T tolerance = std::numeric_limits<T>::epsilon() * rows_;

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      if (std::abs(product(i, j) - identity(i, j)) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

template <typename T> bool Matrix<T>::isPositiveDefinite() const {
  if (!isSquare() || !isSymmetric())
    return false;

  try {
    cholesky();
    return true;
  } catch (const std::runtime_error &) {
    return false;
  }
}

template class Matrix<float>;
template class Matrix<double>;

} // namespace linlang

} // namespace math
