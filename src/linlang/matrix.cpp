#include "matrix.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <new>
#include <stdexcept>
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







template class Matrix<float>;
template class Matrix<double>;

} // namespace linlang

} // namespace math
